#include "voxels.h"
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <math.h>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <unordered_set>

// IO Struct for Input/Output Geometry Parameters
	struct geometry_bin_io
{
    real voxel_size_io;
	int32_t size_x_io;
	int32_t size_y_io;
	int32_t size_z_io;
	real sim_depth_io;
	int32_t sample_height_io;
	std::vector<int16_t> ini_geom_io; 
};

//General IO Helper Functions

//Slicing function
template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}

//Flatten function
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

namespace nbl { namespace geometry {

template<bool gpu_flag>
inline voxels<gpu_flag>::voxels(real voxel_size, vec3 shape, std::vector<int> initial_geometry, int max_save_height, real sim_depth)
{
	_voxel_size = voxel_size;
	_AABB_min = vec3{ 0, 0, 0 };
	//Change __AABB_max to allow scattering under under sample_height --allow bulk material to scatter
	_AABB_max = vec3{ shape.x * voxel_size, shape.y * voxel_size, (shape.z+sim_depth) * voxel_size };
	
	_sample_height = max_save_height-1;
	sample_height_ext = _sample_height*voxel_size; //external sample height in nm 
	_size_x = (int)shape.x;
	_size_y = (int)shape.y;
	_size_z = (int)shape.z;
	_sim_depth = (real)sim_depth; 
	_total_voxels = _size_x*_size_y*_size_z;

	_max_save_height = max_save_height;
	_min_save_height = max_save_height - 5; // save at least 5 layers of voxels
	
	const vec3 m = _AABB_max - _AABB_min;
	_max_extent = magnitude(m);

	//Nebula Dynamic
	_mat_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	_tag_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	_e_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	_species_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0); // old classificiation system 

	//NBD21 Additions --
	_new_species_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0); //New Species Grids
	_residence_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	_surface_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0); //Surface grid encapsulates full geometry for now(can be made more performant)


	if (initial_geometry.size() != _total_voxels)
	{
		//throw std::invalid_argument("initial geometry of wrong shape");
		std::clog << "Warning: initial geometry of wrong shape\n";
	}
	
	_mat_grid = initial_geometry;

	//Adsorpotion and Surface Diffusion Variables (Should be ported to options and condensed)
	//real _diffusion_constant = 6.0667e-8; //surface diffusion consant [cm]2 [s]-1
	real _diffusion_constant =  1e-9; //surface diffusion consant [cm]2 [s]-1
	real beam_current = 100e-12; //beam current  [A]
	real particle_charge = 1.602e-19; //charge of electron  [C]
	_particles_between_gas_routines = 1; //amount of electrons between gas handling routines
	real time_between_particles = particle_charge/beam_current; // time between electrons based on beam current
	_time_between_gas_routines = _particles_between_gas_routines*time_between_particles; //gas handling time step
	_sticking_c = 1; // sticking coefficient
	real voxel_size_cm = voxel_size/1e7; //voxel_size in cm from 0.3nm --> cm
	_jumps_routine = sqrt(4*_diffusion_constant*_time_between_gas_routines)/voxel_size_cm; //jumps per diffusion routine
	_jumps = 0; //
	_jumps += _jumps_routine; //fractional jumps available per routine, initial increment
	real precursor_flux = 8.223e19; //Approximate local gas pressure [cm]2 [s]-1 8.223e17(MTL), 8,.223e20(RRL) Smith 2007
	//real precursor_flux = 5.7575e20; //Approximate local gas pressure [cm]2 [s]-1 table 1 Smith 2008(Seems wrong tbh--This makes it RRL)
	_mean_residence_time = 1;
	_flux_constant = precursor_flux*voxel_size_cm*voxel_size_cm*_time_between_gas_routines;//product of all other constants --> multiply by number of open sites to get number of adsorbed molecules for time_step
	_adsorbs = 0; //fractional adsorbtions pet routine
	int init_site_id = 2; //configures initial site id
	int site_index;
	int initial_surface_layer_size = _size_x*_size_y;
	
	// Gas Handling 
	_adsorption = true; //toggles adsorption on/off
	_desorption = false; //remove this for performance when comparing to Smith et al. --> mean residence time =1s super high (1 million electrons ~ 0.0053s at 30pA)
	_surface_diffusion = true; //toggles surface diffusion on/off
	bool monolayer_toggle = false; //toggle between initial monolayer and no intitial monolayer
	_gas_material_toggle = _adsorption||_surface_diffusion; //for now only activate with adsorption or surface diffusion
	bool new_simul = true; //false - reconstructs surface characterstics with scans
	_smith_source = false; // Allows only 5000 adsorbates to orignate from sources --> one can change this number by chaning _start_amount
	_smith_seq = false; //sequential diffusion per adsorbate for full time --> only for debugging 1 electron inputs(Depreciated - Functionality scrapped)
	_surface_update_override = false; //forces updates of surface on deposition despite other options(Use this to test surface reconstruction of reconstructed geometries)
	_diffusion_optimization = false; //allows movement of empty sites instead of full sites (WIP- Alot of this code is removed)
	_diffusion_concentric = false; // bool to toggle diffusion update from center of interaction region outwards in order to improve overall diffusion efficacy (WIP-Alot of this code is removed)

	//Source Initialization and Adsorbate Tracking
	_central_source = false; //in the center
	_side_source = true; //this is on the surface --> create square around original surface layer
	_perm_source = true; //keep sources filled
	_surface_track = false; //track certain original surface values

	//Generate nearest_neighbours
	generate_neighbour_incs();
	
	//Bulk Propagation
	_bulk_mirrors = true; //enables bulk mirrors (turn off for depositions to prevent mirrored electrons interfereing in growth
	_bulk_transport = true; //enables transport of electrons from simulation depth bulk to voxel substrate layer. \

	//Surface IDs, (2/3) useful for optimization(let empty spots or full spots diffuse )
	// 0 - No Trait
	// 1 - Material Voxel
	// 2 - Open Site
	// 3 - Occupied Site
	// 4 - Source -- Will always stay occupied with gas (no deposition allowed)
	// 5 - Tracked Adsorbate

	//Initial surface id for adsrobates
	if (monolayer_toggle)
	{
		init_site_id = 3;
	}
	
	//reserve total amount of space for 1 layer
	_residence_vec.reserve(initial_surface_layer_size);
	_open_vec.reserve(initial_surface_layer_size);

	//Surface Initialization 
	if (new_simul)
	{
		//2D Surface Flat Surface
		for (int i = 0; i < initial_surface_layer_size; i++) 
		{
			site_index = i+(_sample_height-1)*(_size_x*_size_y); //index of vacuum voxels directly above surface
			_surface_grid.at(i+_sample_height*_size_x*_size_y) = 1; // Initialize surface voxels
			_surface_grid.at(site_index) = init_site_id; // Initialize surface grid with open or occupied sites
			if (monolayer_toggle)
			{
				_residence_vec.push_back(site_index);
				if ((_gas_material_toggle)&& (_mat_grid.at(site_index) == -123)) //check vaccuum(prevents mirror/detector material from being changed)
				{
					_mat_grid.at(site_index) = 1; //changes material values for adsorbate locations
				}
			}
			else
			{
				_open_vec.push_back(site_index);
			}
		}
	}
	else
	{
		std::vector<int> start_search{160,160,50};//index from top/+-y+-x of simulation domain to start looking for material -- Note: useful if geometry of depositon geometry known beforehand)
		int scan_depth = 3; //depth of scanning
		reconstruct_surface_grid(start_search,scan_depth); //reconstructs surface grid from geometry
	}

	// Source Initialization
	_start_amount = 5000; //caps amount of adsorbates that appear from  point source

	//Central Source
	int source_index;
	if (_central_source)
	{
		source_index = (int)(_size_x/2 + (_size_y/2 * _size_y) + (_sample_height-1)*( _size_x * _size_y)); //Middle of Domain
		_surface_grid.at(source_index) = 4; //source designation
		_mat_grid.at(source_index) = 1; //adsorbate designation
		if (_smith_source)
		{
			
			_track_map.insert(std::pair<int, int>(source_index, 0)); //position, tracking id
			_track_vec.push_back(0);
			for (int i = 0; i < _start_amount; i++)
			{
				_residence_vec.push_back(source_index); //allow 5000 "molecules" to occupy center at program start
			}
		}
		else
		{
			_residence_vec.push_back(source_index); //unlimited source
		}	
	}
	//Side Source
	if (_side_source)
	{
		//create sources along x-axis,y-axis 
		for (int i = 0; i < _size_x; i++)
		{
			int index_x_m = (i) + (0) * _size_y + (_sample_height-1) * _size_x * _size_y; //at y = 0
			int index_x_p = (i) + (_size_y-1) * _size_y + (_sample_height-1) * _size_x * _size_y; //at y = size_y
			int index_y_m = (0) + (i) * _size_y + (_sample_height-1) * _size_x * _size_y; //at x = 0 
			int index_y_p = (_size_x-1) + (i) * _size_y + (_sample_height-1) * _size_x * _size_y; //at x = size_x
			_surface_grid.at(index_x_m) = 4; //source
			_surface_grid.at(index_x_p) = 4;
			_surface_grid.at(index_y_m) = 4;
			_surface_grid.at(index_y_p) = 4;
			if (_surface_track) //track original adsorbates
			{
				//Tracking Map
				_track_map.insert(std::pair<int, int>(index_x_m, i));
				_track_map.insert(std::pair<int, int>(index_x_p, i+_size_x));
				_track_map.insert(std::pair<int, int>(index_y_m, i+2*_size_x));
				_track_map.insert(std::pair<int, int>(index_y_p, i+3*_size_x));

				//Start Map
				_start_map.insert(std::pair<int, int>(i,index_x_m));
				_start_map.insert(std::pair<int, int>(i+_size_x,index_x_p));
				_start_map.insert(std::pair<int, int>(i+2*_size_x,index_y_m));
				_start_map.insert(std::pair<int, int>(i+3*_size_x,index_y_p ));

			}
			// Add adsorbates to proper vectors
			_residence_vec.push_back(index_x_m);
			_residence_vec.push_back(index_x_p);
			_residence_vec.push_back(index_y_m);
			_residence_vec.push_back(index_y_p);
			_track_vec.push_back(0);
			_track_vec.push_back(0);
			_track_vec.push_back(0);
			_track_vec.push_back(0);
		}

		// Remove duplicates in corners (expensive but only happens once)
		sort( _residence_vec.begin(), _residence_vec.end() );
		_residence_vec.erase( unique( _residence_vec.begin(), _residence_vec.end() ), _residence_vec.end() );
	}
}


template<bool gpu_flag>
CPU void voxels<gpu_flag>::destroy(voxels<gpu_flag> & geometry)
{
	detail::voxels_factory<gpu_flag>::destroy(geometry);
}

template <bool gpu_flag>
void voxels<gpu_flag>::reconstruct_surface_grid(std::vector<int> start_search,int scan_depth)
{
	// Simple Deposition Geometry - Scans Surface of Deposition and finds appropriate surface values.
	//Build _residence and _open vectors based on input geometry -- Note: Does not work for deposition surfaces with cavities at the moment --> can be improved through clever mapping
	// manual assignment of start_search_index improves performance for known geometries
	//Also still quite expensive(especially deep scans) at the moment for large geometries(only needed once though!) (alternative is demand surface information as program argument or improve current method)
	int search_index;
	int material_index;
	bool surface_found;
	std::vector<int> pos_index (3);

	//Higher scan depth values can be used for more complex geometries(-z and +-y scans can be added as well)
	int scan_counter = 0;

	// Top Down Scan
	for (int i = 0; i < _size_x; i++) 
	{
		for(int j = 0; j < _size_y; j++)
		{
			surface_found = false;
			for (int k = start_search.at(2);k < _sample_height+1; k++)
			{
				search_index = i+j*(_size_y)+ k*(_size_x*_size_y);
				material_index = _mat_grid.at(search_index);
				//searches for the first surface voxel with material = 0 from the top of the simulation domain
				if (material_index == 0)
				{
					_surface_grid.at(search_index) = 1; //give surface voxel designation
					update_surface_grid(search_index,false); //update all voxels around this area
					//Save Height to Geometry Height
					if (k<_min_save_height)
					{
						_min_save_height = k;
					}
					scan_counter+=1;
					if (scan_counter >= scan_depth)
					{
						break;
					}
					
				}
			}
			scan_counter = 0; 
		}
	} 

	int x_lim = _size_x-start_search.at(0); //unused atm (the total scanning time only takes a second or so anyway)
	
	// Positive X Scan
	scan_counter = 0;
	for (int k = start_search.at(2); k < _sample_height; k++) 
	{
		for(int j = 0; j < _size_y; j++)
		{
			surface_found = false;
			for (int i = 0;i < _size_x; i++)
			{
				search_index = i+j*(_size_y)+ k*(_size_x*_size_y);
				material_index = _mat_grid.at(search_index);
				//searches for the first surface voxel with material = 1 from the +x of the simulation domain
				if (material_index == 0)
				{
					_surface_grid.at(search_index) = 1; //give surface voxel designation
					update_surface_grid(search_index,false); //update all voxels around this area
					scan_counter+=1;
					if (scan_counter >= scan_depth)
					{
						break;
					}
					
				}
			}
			scan_counter = 0;  
		}
	} 

	// Negative X Scan
	scan_counter = 0;
	for (int k = start_search.at(2); k < _sample_height; k++) 
	{
		for(int j = 0; j < _size_y; j++)
		{
			surface_found = false;
			for (int i = _size_x;i < 1; i--)
			{
				search_index = i+j*(_size_y)+ k*(_size_x*_size_y);
				material_index = _mat_grid.at(search_index);
				//searches for the first surface voxel with material = 1 from the -x side of the simulation domain
				if (material_index == 0)
				{
					_surface_grid.at(search_index) = 1; //give surface voxel designation
					update_surface_grid(search_index,false); //update all voxels around this area
					scan_counter+=1;
					if (scan_counter >= scan_depth)
					{
						break;
					}
					
				}
			}
			scan_counter = 0;  
		}
	} 

	// Positive Y Scan
	scan_counter = 0;
	for (int k = start_search.at(2); k < _sample_height; k++) 
	{
		for(int i = 0; i < _size_x; i++)
		{
			surface_found = false;
			for (int j = 0;j < _size_y; j++)
			{
				search_index = i+j*(_size_y)+ k*(_size_x*_size_y);
				material_index = _mat_grid.at(search_index);
				//searches for the first surface voxel with material = 1 from the -x side of the simulation domain
				if (material_index == 0)
				{
					_surface_grid.at(search_index) = 1; //give surface voxel designation
					update_surface_grid(search_index,false); //update all voxels around this area
					scan_counter+=1;
					if (scan_counter >= scan_depth)
					{
						break;
					}
					
				}
			}
			scan_counter = 0;  
		}
	} 

	// Negative Y  Scan
	scan_counter = 0;
	for (int k = start_search.at(2); k < _sample_height; k++) 
	{
		for(int i = 0; i < _size_x; i++)
		{
			surface_found = false;
			for (int j = _size_y;j < 1; j--)
			{
				search_index = i+j*(_size_y)+ k*(_size_x*_size_y);
				material_index = _mat_grid.at(search_index);
				//searches for the first surface voxel with material = 1 from the -x side of the simulation domain
				if (material_index == 0)
				{
					_surface_grid.at(search_index) = 1; //give surface voxel designation
					update_surface_grid(search_index,false); //update all voxels around this area
					scan_counter+=1;
					if (scan_counter >= scan_depth)
					{
						break;
					}
					
				}
			}
			scan_counter = 0;  
		}
	}

	// Remove duplicates in open_vec (expensive but only happens once)
	sort( _open_vec.begin(), _open_vec.end() );
	_open_vec.erase( unique( _open_vec.begin(), _open_vec.end() ), _open_vec.end() ); 
}

template<bool gpu_flag>
PHYSICS bool voxels<gpu_flag>::in_domain(vec3 pos)
{
	return ((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
		&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
		&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z));
}

template<bool gpu_flag>
PHYSICS intersect_event voxels<gpu_flag>::propagate(vec3 start, vec3 direction, real distance,
	triangle const* ignore_triangle, int ignore_material) const
{
	intersect_event evt{ distance, nullptr };

	intersect_event evt_debug{ distance, nullptr }; //debug intesection

	const real x = start.x / _voxel_size; // create vars for the location elements
	const real y = start.y / _voxel_size;
	const real z = start.z / _voxel_size;
	
	const real dx = direction.x; // create vars for the direction elements
	const real dy = direction.y;
	const real dz = direction.z;
	
	const vec3 dr = { dx, dy, dz };
	vec3 delta_S = { 0, 0, 0 };

	real delta_s_min = distance / _voxel_size; // holds the shortest pathlength to an intersection

	const int start_mat = ignore_material;

	// calculate the distances to the first planes in the 3 dimensions
	real delta_x;	// delta_x is the distance (perpendicular to the plane) to the next plane in the x direction
	if (dx > 0) {
		delta_x = std::ceil(x) - x;
		delta_S.x = delta_x / dx; // delta_S.x now is the path length to the next plane in the x direction
	}
	else if(dx < 0){
		delta_x = x - std::floor(x);
		delta_S.x = -delta_x / dx;
	}
	else // dx == 0
	{
		delta_S.x = distance / _voxel_size;
	}

	real delta_y;
	if (dy > 0) {
		delta_y = std::ceil(y) - y;
		delta_S.y = delta_y / dy;
	}
	else if(dy < 0){
		delta_y = y - std::floor(y);
		delta_S.y = -delta_y / dy;
	}
	else // dy == 0
	{
		delta_S.y = distance / _voxel_size;
	}

	real delta_z;
	if (dz > 0) {
		delta_z = std::ceil(z) - z;
		delta_S.z = delta_z / dz;
	}
	else if(dz < 0) {
		delta_z = z - std::floor(z);
		delta_S.z = -delta_z / dz;
	}
	else // dz == 0
	{
		delta_S.z = distance / _voxel_size;
	}

	if (delta_S.x < 0.000001)
	{
		delta_S.x += std::abs(1 / dx); // add one path length between two planes to delta_S.x
	}
	if (delta_S.y < 0.000001)
	{
		delta_S.y += std::abs(1 / dy);
	}
	if (delta_S.z < 0.000001)
	{
		delta_S.z += std::abs(1 / dz);
	}
	
	bool bulk_exception_flag = true; // this disables bulk exception code if new conncetion established
	uint64_t isect_id;
	bool outside_bool;


	while (distance / _voxel_size >= delta_s_min) {
		
		//std::clog << "\n" << delta_s_min ;

		// Determine minimum path length from delta_S
		delta_s_min = std::min(std::min(delta_S.x, delta_S.y), delta_S.z);

		int min_index = 0; // min_index is 0 for een intersection with the x-plane, 1 for an intersection with the y-plane
		//and 2 for an intersection with the z-plane
		if(delta_s_min == delta_S.y)
		{
			min_index = 1;
		}
		else if (delta_s_min == delta_S.z)
		{
			min_index = 2;
		}
		
		int min_i = min_index; // store min_index 

		vec3 new_pos = start / _voxel_size + (delta_s_min /*+ 0.001*/) * dr; // new position if there is an intersection in voxels

		vec3 pos = new_pos * _voxel_size; // new position in nm, for check

		if(!((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
			&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
			&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z))) 
			// check whether the electron is still in the simulation domain, if not, return evt and let it scatter outside of geometry and terminate
		{
			return evt;
		}

		int k; // indices of the material grid
		int l;
		int m;

		real dx_sgn; // signs of dx, dy and dz
		real dy_sgn;
		real dz_sgn;

		// determine the indices of the next cell that the electron will enter
		switch (min_i)
		{
		case 0: // intersection with x-plane
			dx_sgn = dx / std::abs(dx); // determine sign of dx
			k = static_cast<int>(new_pos.x + 0.1 * dx_sgn); // calculate the indices by flooring the new position plus a
			// little inward 
			l = static_cast<int>(new_pos.y);
			m = static_cast<int>(new_pos.z);
			
			/* Moved this to bulk propagation
			if (k >= _size_x || k < 0)
			{
				return evt;
			}
			*/
			break;

		case 1: // intersection with y-plane
			dy_sgn =  dy / std::abs(dy);
			k = static_cast<int>(new_pos.x);
			l = static_cast<int>(new_pos.y + 0.1 * dy_sgn);
			m = static_cast<int>(new_pos.z);
			/*
			if (l >= _size_y || l < 0)
			{
				return evt;
			}
			*/
			break;

		default: // intersection with z-plane
			dz_sgn = dz / std::abs(dz);
			k = static_cast<int>(new_pos.x);
			l = static_cast<int>(new_pos.y);
			m = static_cast<int>(new_pos.z + 0.1 *dz_sgn);
			/*
			if(m >= _size_z || m < 0)
			{
				return evt;
			}
			*/
			break;
		}

		//Bulk Propagation
		outside_bool = (k >= _size_x || k < 0)  || (l>= _size_y || l<0) || (m >= _size_z || m < 0);
		if ((outside_bool || z > _size_z) && (bulk_exception_flag)) // activate if within bulk substrate,outside of voxel boundaries or propagating into bulk substrate
		{
			bulk_exception_flag = false; //prevents bulk exception code from being run again when transfer to voxel geometry complete
			if ((_bulk_mirrors) || (_bulk_transport))
			{
				real x_intersect;
				real y_intersect;
				real z_intersect;
				real x_t;
				real y_t;
				real z_t;
				int voxel_side;
				int z_transport = _sample_height; //transport layer, vacuum layer directly above voxel substrate surface
				if (dx>0)
				{
					x_intersect = std::abs(_size_x-x)*_voxel_size; //distance to x boundary in nm
				}
				else
				{
					x_intersect = std::abs(x)*_voxel_size; 
				}
				if (dy>0)
				{
					y_intersect = std::abs(_size_y-y)*_voxel_size; //distance to y boundary in nm
				}
				else
				{
					y_intersect = std::abs(y)*_voxel_size;
				}
				//Bulk Transport -- Attempts to smoothly bridge gap between voxel geometry and larger bulk subtrate
				if (dz<0) //z<0 goes up towards substrate
				{
					z_intersect = std::abs(z-z_transport)*_voxel_size; // //distance to sample_height in nm
				}
				else //z>0 goes down
				{
					z_intersect = std::abs(_sim_depth+_size_z-z)*_voxel_size; // distance to bottom of domain in nm
				}
				// how fast each boundary is reached 
				x_t = std::abs(x_intersect/dx);
				y_t = std::abs(y_intersect/dy);
				z_t = std::abs(z_intersect/dz);
				real bulk_boundary_t;
				if (_bulk_transport)
				{
					bulk_boundary_t = std::min(std::min(x_t,y_t),z_t); // closest boundary in nm
				}
				else
				{
					bulk_boundary_t = std::min(x_t,y_t);
					
				}
				bool bulk_mirror = (y_t<z_t) || (x_t<z_t);
				real intersect_distance;
				
				if (bulk_boundary_t == x_t) //reach x boundary fastest --> distance is x_intersect
				{
					intersect_distance = x_intersect;
				}
				else if (bulk_boundary_t == y_t) //reach y boundary fastest --> distance is y_intersect
				{
					intersect_distance = y_intersect;
				}
				else
				{
					intersect_distance = z_intersect; //reach z boundary fastest --> distance is y_intersect
				}
				//intersect_distance = std::sqrt(std::pow(bulk_boundary_t*dx,2)+std::pow(bulk_boundary_t*dy,2)+std::pow(bulk_boundary_t*dz,2));
				if ((intersect_distance<=distance) && (bulk_mirror)) //Force Mirror intersection
				{
					if (x_t<y_t)
					{
						if (dx>0)
						{
							voxel_side = 1;
						}
						else
						{
							voxel_side = 2;
						}
					}
					else
					{
						if (dy>0)
						{
							voxel_side = 3;
						}
						else
						{
							voxel_side = 4;
						}
					}
					
					evt.isect_distance = intersect_distance;
					reinterpret_cast<int32_t*>(&isect_id)[0] = -122; //force mirror (-122) --> should import from materials manager in future
					reinterpret_cast<int32_t*>(&isect_id)[1] = voxel_side;
					evt.isect_triangle = reinterpret_cast<triangle*>(isect_id);
					return evt;
				}
				else if  ((intersect_distance<=distance) && (_bulk_transport) && (!bulk_mirror))
				{
					
					if (dz>0) //going down
					{
						return evt; //return a scattering event 
					}
					else //dz <0 going up into 
					{
						//Rebuild k,l,m values and delta_s_min to inject back into voxel geometry propagation at z_transport
						delta_S.x = std::abs(x_intersect/dx/_voxel_size);
						delta_S.y = std::abs(y_intersect/dy/_voxel_size);
						delta_S.z = std::abs(z_intersect/dz/_voxel_size);
						delta_s_min = std::min(std::min(delta_S.x, delta_S.y), delta_S.z);
						new_pos = start / _voxel_size + (delta_s_min) * dr;
						dz_sgn = dz / std::abs(dz);
						k = static_cast<int>(new_pos.x);
						l = static_cast<int>(new_pos.y);
						m = static_cast<int>(new_pos.z);	
						min_i = 2;													
					}
				}
				else
				{
					return evt; //scatter
				}
			}
		}
		
		int new_mat = _mat_grid[k + l * _size_x + m * _size_x * _size_y]; // determine material using the material indices
	
		//std::clog << "   " << new_mat << "   " << start_mat;
		
		if (new_mat != start_mat) { // if thcere is een intersection, return the intersection event
			
			evt.isect_distance = (delta_s_min /*+ 0.005*/) * _voxel_size; // set the distance to the intersection

			// Determine voxel side
			int	voxel_side;
			switch (min_i)
			{
			case 0:
				if (dx > 0) {
					voxel_side = 1;
				}
				else {
					voxel_side = 2;
				}
				break;
			case 1:
				if (dy > 0) {
					voxel_side = 3;
				}
				else {
					voxel_side = 4;
				}
				break;
			default:
				if (dz > 0) {
					voxel_side = 5;
				}
				else {
					voxel_side = 6;
				}
				break;
			}

			// put an identifier for the new material and a voxel side identifier in the triangle pointer
			reinterpret_cast<int32_t*>(&isect_id)[0] = new_mat;
			reinterpret_cast<int32_t*>(&isect_id)[1] = voxel_side;
			evt.isect_triangle = reinterpret_cast<triangle*>(isect_id);

			//Debugging to prevent longer intesect_distance than scattering distance 
			if (distance <= evt.isect_distance)
			{
				return evt_debug;
			}


			return evt; // return the event
		}

		// if there is no intersection:
		// Calculate new value of delta_S.i for next iteration
		switch (min_i)
		{
		case 0:
			delta_S.x += std::abs(1 / dx);
			break;

		case 1:
			delta_S.y += std::abs(1 / dy);
			break;

		default:
			delta_S.z += std::abs(1 / dz);
			break;
		}
	}
	//std::clog << "\n turn";
	
	// if no intersection was found, return the standard event. 
	return evt;
	
	
}

template <bool gpu_flag>
int voxels<gpu_flag>::get_material(int position) const
{
	return _mat_grid[position];
}

template <bool gpu_flag>
int voxels<gpu_flag>::get_particles_between_gas_routines() const
{
	return _particles_between_gas_routines;
}

template<bool gpu_flag>
PHYSICS real voxels<gpu_flag>::get_sim_depth() const
{
	return _sim_depth;
}

//Surface Functions
template <bool gpu_flag>
void voxels<gpu_flag>::update_surface_grid(int index,bool deposition_update)
{
	_surface_grid.at(index) = 1;

	//Update Surface Grid (27 voxels in 3d neighbourhood,19 in corner resricted neighbourhood incl deposition voxel)
	int surface_id;
	int material_id;

	if (deposition_update) //on if one or other model is on
	{
		_residence_grid.at(index) = 0;
		_residence_vec.erase(std::remove(_residence_vec.begin(), _residence_vec.end(), index), _residence_vec.end());
	}

	int index_i;

	//No corners code(dont allow cube corners to be given surface characteristics)
	for (int i = 0; i< 18; i++)
	{
		index_i = index + _nearest_neighbour_incs[i];
		material_id = _mat_grid[index_i];
		if (material_id == -123) //vacuum or mirror material
		{
			_surface_grid.at(index_i) = 2; //turn into open site any vacuum location next to new deposition(neighbours)
			_open_vec.push_back(index_i);
		}
		if (material_id == 0) //bulk material
		{

			_surface_grid.at(index_i) = 1; //turn into material voxel
		}
	}
}


template <bool gpu_flag>
void voxels<gpu_flag>::generate_neighbour_incs()
{
	int x_c;
	int y_c;
	int z_c;
	int t_c;

	for (int i = 0; i< 18; i++)
	{
		int x_c = _corner_restricted_values_array[i][0];
		int y_c =_corner_restricted_values_array[i][1];
		int z_c = _corner_restricted_values_array[i][2];

		t_c = x_c + y_c*_size_y + z_c*(_size_x*_size_y);

		_nearest_neighbour_incs.push_back(t_c);
	}
}



template <bool gpu_flag>
std::vector<int> voxels<gpu_flag>::surface_check(int index, bool empty_toggle,bool full_toggle)
{
	std::vector<int> grid_indices;
	grid_indices.reserve(18);
	int material_id;
	int surface_id;
	int index_i;

	for (int i = 0; i< 18; i++)
	{
		index_i = index + _nearest_neighbour_incs[i];
		material_id = _mat_grid[index_i];
		surface_id = _surface_grid[index_i];

		if (empty_toggle) //allows only checking for empty sites --> for use when doing surface diffusion(may change return value to integer of voxel index)
		{
			if (surface_id == 2)
			{
				grid_indices.emplace_back(index_i);
			}
		}
		else if (full_toggle) //allows only checking for full(occupied) sites -->  for use when doing surface diffusion(may change return value to integer of voxel index)
		{
			if (surface_id == 3)
			{
				grid_indices.emplace_back(index_i);
			}
		}
		else
		{
			if ((material_id == -123) || (material_id = 1)) // returns vacuum sites/adsorbate sites
			{
				grid_indices.emplace_back(index_i);
			}
		}
	}
	return grid_indices; // no adsorbates or vacuum voxels found
}

// Adsorption Routine
template <bool gpu_flag>
void voxels<gpu_flag>::adsorption_routine()
{

	//Adsorptions Checks
	int total_sites = _residence_vec.size() + _open_vec.size(); //does not change for subroutine
	real gas_deposits_it = _flux_constant*total_sites; //changes as deposition geometry changes
	_adsorbs += gas_deposits_it;
	int current_adsorbs = floor(_adsorbs);

	//Comment out if desorption on (Here to speed up MTL simulations)
	if (current_adsorbs == 0)
	{
		return;
	}

	// Temporary Sets --Faster than multiple O(n) vector deletions/additions
	std::unordered_set<int> residence_set(_residence_vec.begin(),_residence_vec.end()); 
	std::unordered_set<int> open_set(_open_vec.begin(),_open_vec.end()); 

	float site_ratio; //takes ratio 
	int site_index;
	real residence;

	// Desorption Subroutine (When mean residence time >> simulation time this can be turned off for computational performance)
	if (_desorption)
	{
		for (auto& it :	_residence_vec)
		{
			site_index = it;
			if (_surface_grid.at(site_index) == 4) //ignore sources
			{
				continue;
			}
			residence = _residence_grid.at(site_index) + _time_between_gas_routines;
			if (residence < _mean_residence_time) //does not exceed maxium residence time
			{
				_residence_grid.at(site_index) = residence;
			}
			else
			{
				residence_set.erase(site_index);
				_residence_grid.at(site_index) = 0;
				open_set.insert(site_index);
				_mat_grid.at(site_index) = -123; //turns to vacuum
				_surface_grid.at(site_index) = 2; //open site
			}
		}
	}

	/*
	if (current_adsorbs == 0)
	{
		return;
	}
	*/

	int random_index;
	int random_site;

	// Adsorption Subroutine
	site_ratio = open_set.size()/(float)total_sites; //update ratio --> optimization #of gas_deposits << number of open sites per moment
	for (int i = 0; i < current_adsorbs; i++) 
	{
		if ((_rand_state.unit() <= site_ratio) && (_rand_state.unit() <= _sticking_c))//adsorbption falls in open set and passes sticking check
		{
			random_index = int(_rand_state.unit()*open_set.size()); //random index
			random_site = _open_vec[random_index];
			open_set.erase(random_site);
			residence_set.insert(random_site);
			_surface_grid[random_site] = 3;
			_mat_grid[random_site] = 1; //put adsorbate material here
		} 
	}

	//Fractional Sum Decreased
	_adsorbs -= current_adsorbs;

	//Reassign vectors
	_residence_vec.assign(residence_set.begin(), residence_set.end());
	_open_vec.assign(open_set.begin(), open_set.end());

}

//Surface Diffusion
template <bool gpu_flag>
void voxels<gpu_flag>::surface_diffusion_routine() //Iteration over maps for performance
{
	//Stop diffusion if jumps insufficient
	if (_jumps<1)
	{
		_jumps += _jumps_routine;
		return;
	}

	int random_index; 
	int random_site;
	int site;
	std::vector<int> site_indices;
	int track_id;

	//Temporary Sets --Faster than multiple O(n) vector deletions
	//std::unordered_set<int> residence_set(_residence_vec.begin(),_residence_vec.end()); 
	std::set<int> open_set(_open_vec.begin(),_open_vec.end()); 
	int res_sz = _residence_vec.size();
	int open_sz = _open_vec.size();

	bool mv_bool = true; //move residents -true, move open sites -false (WIP)
	if ((open_sz<res_sz) && (_diffusion_optimization))
	{
		mv_bool = false;
	}

	if (_diffusion_concentric)
	{
		//WIP (Very Slow version discarded -- To reader--> 
		//Make custom surface data structure that keeps track of distance from center ect. to prevent std::set/std::vector conversions and allow performant concentric diffusion)

	}

	std::random_shuffle (_residence_vec.begin(), _residence_vec.end()); //better than no randomization (optimization)

	if (mv_bool)
	{
		//while _jumps > 1 --> keeps track of fractional jumps.
		while (_jumps >1)
		{
			//std::random_shuffle (_residence_vec.begin(), _residence_vec.end()); //Shuffle residence vector
			res_sz = _residence_vec.size();

			// 1 Diffusion Step per Adsorbate 
			for(int i=0; i < res_sz; i++)
			{
				site = _residence_vec.at(i); //location on grid
				site_indices = surface_check(site,true,false);
				if (site_indices.size()<1)
				{
					continue; //skip adsorbate
				}
				random_index = rand() % site_indices.size();		
				random_site = site_indices.at(random_index);

				//tracking adsorbates
				if (_surface_track)
				{
					if (_track_map.count(site)!=0)
					{
						track_id = _track_map.at(site);
						// spontaneus
						if ((track_id == 0) && (_central_source)) 
						{
							track_id = _track_vec.size(); //create new id(starts from 1)
							_track_vec.push_back(0); //push back new distance value for id
						}
						_track_map.insert(std::pair<int, int>(random_site, track_id)); //move tracker to new site
						_track_map.erase(site); //remove tracker old site
						_track_vec.at(track_id)+=1; //increase tracker id distance travelled by 1

						//creates 0 id tracker at until _start_amount reached --> first id not used
						if  ((_track_vec.size()<_start_amount+1) && (_central_source))
						{
							_track_map.insert(std::pair<int, int>(random_site, 0));
						}
					} 
				}
				//Move Values
				_residence_vec.at(i) = random_site ;
				open_set.erase(random_site);
				
				//Update grids
				_surface_grid[random_site] = 3;
				_mat_grid[random_site] = 1;

				if (_surface_grid[site] != 4) //If not source
				{
					open_set.insert(site); //using set as we dont know position of int
					_surface_grid[site] = 2;
					_mat_grid[site] = -123;
				}
				else
				{
					//allow infinite source (put a randomness check here for varied replacement at source locations(removed for performance reasons))
					if (_perm_source)
					{
						_residence_vec.emplace_back(site);
					}
				}
			}
			//jump_cycle completion
			_jumps-=1;
		}
	}
	else
	{
		//Under Construction
		/*
		//Shuffle open vector
		std::random_shuffle (_open_vec.begin(), _open_vec.end());
		
		// 1 Diffusion Step per Empty Site
		for(int i=0; i < res_sz; i++) //uses 
		{
			site = _open_vec.at(i%);
			vec_site =  grid_index_to_vec3(site);
			site_indeces = surface_check(vec_site,false,true);
			if (site_indeces.size()<1)
			{
				continue;
			}
			random_index = (int)(site_indeces.size()*_rand_state.unit());		
			random_site = site_indeces.at(random_index);

			//tracking adsorbates
			if (_surface_track)
			{
				if (_track_map.count(random_site)!=0)
				{
					track_id = _track_map.at(random_site);
					// spontaneus
					if ((track_id == 0) && (_central_source)) 
					{
						track_id = _track_vec.size(); //create new id(starts from 1)
						_track_vec.push_back(0); //push back new distance value for id
					}
					_track_map.insert(std::pair<int, int>(site, track_id)); //move tracker to new site
					_track_map.erase(random_site); //remove tracker old site
					_track_vec.at(track_id)+=1; //increase tracker id distance travelled by 1

					//creates 0 id tracker at until _start_amount reached --> first id not used
					if  ((_track_vec.size()<_start_amount+1) && (_central_source))
					{
						_track_map.insert(std::pair<int, int>(site, 0));
					}
				} 
			}
			//Switch Values
			_residence_vec.at(random_index) = site ;
			_open_vec.at(i) = random_site;
			//Update grids
			_surface_grid.at(site) = 3;
			_mat_grid.at(site) = 1;

			if (_surface_grid.at(random_site) != 4)
			{
				_surface_grid.at(random_site) = 2;
				_mat_grid.at(random_site) = -123;
			}
			else
			{
				//allow infinite source
				if (_perm_source)
				{
					_residence_vec.push_back(random_site);
				}
			}
		}
		*/
	}

	//Fractional jumps increase
	_jumps += _jumps_routine;

	//Reassign vectors
	_open_vec.assign(open_set.begin(), open_set.end());

}

//Surface Diffusion with Sets
template <bool gpu_flag>
void voxels<gpu_flag>::surface_diffusion_routine_set() 
{
	int random_index; 
	int random_site;
	int site;
	real residence;
	std::vector<int> site_indices;
	int track_id;
	bool mv_bool = true;

	std::set<int> open_set(_open_vec.begin(),_open_vec.end()); 
	std::set<int> residence_set(_residence_vec.begin(),_residence_vec.end()); 

	std::random_shuffle (_residence_vec.begin(), _residence_vec.end()); //better than no randomization (optimization)

	if (mv_bool)
	{
		//while _jumps > 1 --> keeps track of fractional jumps.
		while (_jumps >1)
		{
			/*
			// 1 Diffusion Step per Adsorbate 
			for(auto elem)
			{
				site = _residence_vec.at(i); //location on grid
				vec_site =  grid_index_to_vec3(site); // location on grid to indices k,l,m
				site_indices = surface_check(vec_site,true,false);
				if (site_indices.size()<1)
				{
					continue; //skip adsorbate
				}
				random_index = (int)((site_indices.size()-1)*_rand_state.unit());		
				random_site = site_indices.at(random_index);

				//tracking adsorbates
				if (_surface_track)
				{
					if (_track_map.count(site)!=0)
					{
						track_id = _track_map.at(site);
						// spontaneus
						if ((track_id == 0) && (_central_source)) 
						{
							track_id = _track_vec.size(); //create new id(starts from 1)
							_track_vec.push_back(0); //push back new distance value for id
						}
						_track_map.insert(std::pair<int, int>(random_site, track_id)); //move tracker to new site
						_track_map.erase(site); //remove tracker old site
						_track_vec.at(track_id)+=1; //increase tracker id distance travelled by 1

						//creates 0 id tracker at until _start_amount reached --> first id not used
						if  ((_track_vec.size()<_start_amount+1) && (_central_source))
						{
							_track_map.insert(std::pair<int, int>(random_site, 0));
						}
					} 
				}
				//Move Values
				_residence_vec[i] = random_site ;
				open_set.erase(random_site);
				
				//Update grids
				_surface_grid[random_site] = 3;
				_mat_grid[random_site] = 1;

				if (_surface_grid[site] != 4) //If not source
				{
					open_set.insert(site); //using set as we dont know position of int
					_surface_grid[site] = 2;
					_mat_grid[site] = -123;
				}
				else
				{
					//allow infinite source (put a randomness check here for varied replacement at source locations(removed for performance reasons))
					if (_perm_source)
					{
						_residence_vec.emplace_back(site);
					}
				}
			}
			*/

			//jump_cycle completion
			_jumps-=1;

		}
	}
	//Fractional jumps increase
	_jumps += _jumps_routine;
	//Reassign vectors
	_open_vec.assign(open_set.begin(), open_set.end());
	_residence_vec.assign(residence_set.begin(), residence_set.end());

}

template <bool gpu_flag>
std::vector<int> voxels<gpu_flag>::grid_index_to_vec3(int index)
{
	int count = index;
	int z = floor(count/(_size_x*_size_y));
	count -= z*(_size_x*_size_y);
	int y = floor(count/(_size_y));
	count -= y*(_size_y);
	int x = floor(count);
	std::vector<int> vector_pos{x,y,z};
	return vector_pos;

}


//Deposition
template <bool gpu_flag>
void voxels<gpu_flag>::deposit(vec3 position, vec3 normal, int material, int PE_tag, real energy, uint8_t species, std::deque<char> new_species, vec3 initial_vel)
{
	vec3 pos = position / _voxel_size + 0.1 * normal;
	
	int k = static_cast<int>(std::floor(pos.x));
	int l = static_cast<int>(std::floor(pos.y));
	int m = static_cast<int>(std::floor(pos.z));

	int index = k+ l*_size_y + m*(_size_x*_size_y);

	//Make sure that deposition location is proper when gas handling enabled
	if ((_mat_grid[index] != 1) && (_gas_material_toggle)) 
	{
		pos = position / _voxel_size - 0.1 * normal;
		k = static_cast<int>(std::floor(pos.x));
		l = static_cast<int>(std::floor(pos.y));
		m = static_cast<int>(std::floor(pos.z));
		index = k+ l*_size_y + m*(_size_x*_size_y);
	}
	//Make sure that deposition location is proper when gas handling disabled
	if ((_mat_grid[index] != -123) && (!_gas_material_toggle)) 
	{
		pos = position / _voxel_size - 0.1 * normal;
		k = static_cast<int>(std::floor(pos.x));
		l = static_cast<int>(std::floor(pos.y));
		m = static_cast<int>(std::floor(pos.z));
		index = k+ l*_size_y + m*(_size_x*_size_y);
	}
	if (m<0) //small fix if height limit exceeded --> make geometry larger to accomadate larger growths
	{
		std::clog << "Geometry Size Insufficient..." << std::endl;
	}

	_mat_grid[index] = material;
	_tag_grid[index] = PE_tag + 1; // PE_tag begins at 0, so we add 1 to distinguish the deposition from the first electron from the non-deposits
	_e_grid[index] = energy;
	_species_grid[index] = species + 1;

	// New Species Internal Classification -Need to increase performance later
	int new_species_current;
	int event_length = new_species.size();
	std::vector<char> new_species_vec = {new_species.begin(), new_species.end()};
	std::vector<char> new_species_vec_begin = slice(new_species_vec,0,3);
	std::vector<char> new_species_vec_end = slice(new_species_vec,0,4);
	std::vector<char> new_species_vec_end_3 = slice(new_species_vec,0,4);
	if (event_length > 4 )
	{
		new_species_vec_end = slice(new_species_vec,4,7);
		new_species_vec_end_3 = slice(new_species_vec,5,7);
	}
	

	//New Species Check Values
	std::vector<char> new_species_check_1 = {0,'V','V','F'};
	std::vector<char> new_species_check_2 = {0,'D','V','F'};
	std::vector<char> new_species_check_2_b = {0,'D','D','F'};
	std::vector<char> new_species_check_3 = {0,'D','V','B'};
	std::vector<char> new_species_check_4 = {0,'S','V','B'};
	std::vector<char> new_species_check_5 = {0,'V','D','F'};
	std::vector<char> new_species_check_5_b = {0,'V','S','F'};
	std::vector<char> new_species_check_6 = {0,'D','V','B',0,'V','D','B'};
	std::vector<char> new_species_check_7 = {0,'S','V','B',0,'V','D','B'};

	std::vector<char> new_species_check_8 = {'S','V','B'};
	std::vector<char> new_species_check_9 = {'D','V','B'};
	std::vector<char> new_species_check_10 = {'D','V','F'};

	//New Species Check values Gas Handling
	std::vector<char> new_species_check_3_b = {0,'D','D','B'};

	// Primary Electron Classification --Need to update for Gas Handling --(Back scatters dont track well)
	//-(Change max event length and new_species checks to create custom classifications based on recent boundary-intersection event history)
	if (new_species[0] == 0)
	{
		if ((new_species_vec_begin == new_species_check_1))  //[0,V,V,F]
		{
			new_species_current = 1; //PE
		}
		else if ((new_species_vec_end == new_species_check_5) || (new_species_vec_end == new_species_check_5_b)) //[0,V,D,F]
		{
			new_species_current = 5; //PE_FSE_DVD
		}
		else if ((new_species_vec_end == new_species_check_2) || (new_species_vec_end == new_species_check_2_b)) //[0,D,V,F]
		{
			new_species_current = 2; //PE_FSE_D
		} 
		else if ((new_species_vec_end == new_species_check_3) || (new_species_vec_end == new_species_check_3_b))  //[0,D,V,B]
		{
			new_species_current = 3; //PE_BSE_D
		}
		else if (new_species_vec_end == new_species_check_4) //[0,S,V,B]
		{
			new_species_current  = 4; //PE_BSE_S
		}
		else if (new_species_vec == new_species_check_6)  //[0,D,V,B,0,V,D,B]
		{
			new_species_current  = 6; //PE_BSE_DVD
		}
		else if (new_species_vec == new_species_check_7) //[0,S,V,B,0,V,D,B]
		{
			new_species_current = 7; //PE_BSE_SVD
		}
		else
		{
			new_species_current = 12; //uncategorized
		}
	}
	//Secondary Electron Classification
	else if (new_species[0] == 1) 
	{
		new_species_current = 8; //SE_PE
	}
	else if (new_species[0] == 2) 
	{
		new_species_current = 9; //SE_FSE
	}
	else if (new_species[0] == 3) 
	{
		new_species_current = 10; //SE_BSE_D
	}
	else if (new_species[0] == 4) 
	{
		new_species_current = 11; //SE_BSE_S
	}
    else
	{
        new_species_current = 12; //Uncategorized
	}

	if (new_species_current>12)
	{
		std::clog << "Error in Classification..." << std::endl;
	}
	_new_species_grid[index] = new_species_current;
	
	/* External Python Version

	#Current Classification (WIP)
	['0 Empty', white
    ' 1 PE', red
    '2 PE_FSE_D', green
    '3  PE_BSE_D', blue
    '4  PE_BSE_S', blue
    '5  PE_FSE_DVD', #dc5858
    '6  PE_BSE_DVD', blue
    '7  PE_BSE_SVD',blue
    '8  SE_PE', cyan
    '9  SE_FSE', '#80ff80'
    '10  SE_BSE_D', cyan
    '11  SE_BSE_S', #ffa020
    '12  Uncategorized']
	*/

	if (m > _max_save_height)
	{
		_max_save_height = m;
	}
	else if (m < _min_save_height)
	{
		_min_save_height = m;
	}

	std::vector<int> posindices{k,l,m};
	//Update Surface Voxels
	if (_surface_diffusion || _adsorption || _surface_update_override)
	{
		update_surface_grid(index,true); // updates local designations
	}
	

}

template<bool gpu_flag>
PHYSICS real voxels<gpu_flag>::get_max_extent() const
{
	return _max_extent;
}

template<bool gpu_flag>
PHYSICS real voxels<gpu_flag>::get_total_voxels() const
{
	return _total_voxels;
}

template<bool gpu_flag>
inline PHYSICS vec3 voxels<gpu_flag>::AABB_min() const
{
	return _AABB_min;
}

template<bool gpu_flag>
inline PHYSICS vec3 voxels<gpu_flag>::AABB_max() const
{
	return _AABB_max;
}

template<bool gpu_flag>
inline PHYSICS bool voxels<gpu_flag>::get_adsorption() 
{
	return _adsorption;
}

template<bool gpu_flag>
inline PHYSICS bool voxels<gpu_flag>::get_surface_diffusion() 
{
	return _surface_diffusion;
}

//Output Functions

template <bool gpu_flag>
void voxels<gpu_flag>::save_bin_output(const std::string file_name)
{
	bool enable_e_grid = false;
	//Save Height Restrictions
	int lower_save_index = _min_save_height * _size_x * _size_y;
	int upper_save_index = ((_max_save_height-1) * _size_x * _size_y)-1; //alternative to get rid of all bulk material in output
	int32_t _save_height = (_max_save_height-1)-_min_save_height; 

	//Save Slices
	std::vector<real> _e_grid_slice;
	std::vector<int> _mat_grid_slice = slice(_mat_grid,lower_save_index,upper_save_index);
	std::vector<int> _tag_grid_slice = slice(_tag_grid,lower_save_index,upper_save_index);
	std::vector<int> _species_grid_slice = slice(_species_grid,lower_save_index,upper_save_index);
	std::vector<int> _new_species_grid_slice = slice(_new_species_grid,lower_save_index,upper_save_index);

	//Save Slice Conversion
	std::vector<int16_t> _mat_grid_slice_int16(_mat_grid_slice.begin(),_mat_grid_slice.end());
	std::vector<int32_t> _tag_grid_slice_int32(_tag_grid_slice.begin(),_tag_grid_slice.end());
	std::vector<int16_t> _species_grid_slice_int16(_species_grid_slice.begin(),_species_grid_slice.end());
	std::vector<int16_t> _new_species_grid_slice_int16(_new_species_grid_slice.begin(),_new_species_grid_slice.end());

	if (enable_e_grid)
	{
		_e_grid_slice = slice(_e_grid,lower_save_index,upper_save_index);
	}
	
	//Open and determine length of integer stream (for now we save everything --> might need to change to save space/speed-->change int32_t)
	std::ofstream output_bin_file(file_name, std::ios::binary);
	int64_t len = _mat_grid_slice.size();
	output_bin_file.write( (char*)&len, sizeof(len) );

	output_bin_file.write( (char*)&_voxel_size, sizeof(real) );
	output_bin_file.write( (char*)&_size_x, sizeof(int32_t) );
	output_bin_file.write( (char*)&_size_y, sizeof(int32_t) );
	output_bin_file.write( (char*)&_save_height, sizeof(int32_t) );

	//Vectors
	output_bin_file.write( (char*)&_mat_grid_slice_int16[0], len * sizeof(int16_t) );
	output_bin_file.write( (char*)&_tag_grid_slice_int32[0], len * sizeof(int32_t) );
	if (enable_e_grid)
	{
		output_bin_file.write( (char*)&_e_grid_slice[0], len * sizeof(real) );
	}
	output_bin_file.write( (char*)&_species_grid_slice_int16[0], len * sizeof(int16_t) );
	output_bin_file.write( (char*)&_new_species_grid_slice_int16[0], len * sizeof(int16_t) );

    output_bin_file.close();
}

//Save Surface Information Seperately --> best for debugging new models

template <bool gpu_flag>
void voxels<gpu_flag>::save_bin_surface(const std::string file_name)
{
	//Save Height Restrictions
	int lower_save_index = _min_save_height * _size_x * _size_y;
	int upper_save_index = ((_max_save_height-1) * _size_x * _size_y)-1; 
	int32_t _save_height = (_max_save_height-1)-_min_save_height; 

	// Tracked Vector IO
	std::vector<int> track_ids;
	std::vector<int> start_pos;
	std::vector<int> end_pos;
	int id;
	int grid_index;
	//int index_loss = (_size_z-_save_height)* _size_x * _size_y; unused

	//Tracking loop
	for (auto& it : _track_map)
	{
		//Identify tracked adsorbates on surface grid
		grid_index = it.first;
		_surface_grid.at(grid_index) = 5;

		//Id vectors
		id = it.second;
		track_ids.push_back(id);
		start_pos.push_back(_start_map[id]);
		end_pos.push_back(grid_index);
	}

	//Save Slices
	std::vector<int> _mat_grid_slice = slice(_mat_grid,lower_save_index,upper_save_index);
	std::vector<int> _surface_grid_slice = slice(_surface_grid,lower_save_index,upper_save_index);

	//Save Slice Conversion
	std::vector<int16_t> _mat_grid_slice_int16(_mat_grid_slice.begin(),_mat_grid_slice.end());
	std::vector<int16_t> _surface_grid_slice_int16(_surface_grid_slice.begin(),_surface_grid_slice.end());
	std::vector<int16_t> _track_vec_slice_int16(_track_vec.begin(),_track_vec.end());

	//Open and determine length of integer stream 
	std::ofstream output_bin_file(file_name, std::ios::binary);
	int64_t len = _mat_grid_slice.size();
	int64_t track_sz = _track_vec.size();
	output_bin_file.write( (char*)&len, sizeof(len) );

	output_bin_file.write( (char*)&_voxel_size, sizeof(real) );
	output_bin_file.write( (char*)&_size_x, sizeof(int32_t) );
	output_bin_file.write( (char*)&_size_y, sizeof(int32_t) );
	output_bin_file.write( (char*)&_save_height, sizeof(int32_t) );

	//Grid Vectors
	output_bin_file.write( (char*)&_mat_grid_slice_int16[0], len * sizeof(int16_t) );
	output_bin_file.write( (char*)&_surface_grid_slice_int16[0], len * sizeof(int16_t) );

	//Tracking
	output_bin_file.write( (char*)&track_sz, sizeof(track_sz) );
	output_bin_file.write( (char*)&track_ids[0], track_sz * sizeof(int) );
	output_bin_file.write( (char*)&start_pos[0], track_sz * sizeof(int) );
	output_bin_file.write( (char*)&end_pos[0], track_sz * sizeof(int) );

	//New for Surface Tracking certain adsorbates
	//output_bin_file.write( (char*)&len_distance, sizeof(int32_t) );
	//output_bin_file.write( (char*)&_track_vec_slice_int16[0], len * sizeof(int16_t) );
	

    output_bin_file.close();
}

template<bool gpu_flag>
CPU void voxels<gpu_flag>::set_AABB(vec3 min, vec3 max) 
{
	_AABB_min = min;
	_AABB_max = max;

	const vec3 m = max - min;
	_max_extent = magnitude(m);
}

template<bool gpu_flag>
CPU void voxels<gpu_flag>::set_rand_state(util::random_generator<false> rand_state_in) 
{
	_rand_state = rand_state_in;
}

}} // namespace nbl::geometry


