#ifndef __GEOMETRY_VOXELS_H_
#define __GEOMETRY_VOXELS_H_

#include "../core/triangle.h"
#include "../core/events.h"
#include "../common/util/random.h"
#include <deque>
#include <map>
#include <set>
#include <tuple>

namespace nbl { namespace geometry {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 */
	template<bool gpu_flag>
	struct voxels_factory;
}


/**
 * \brief Stores geometry as a list of voxel cells.
 *
 * This class is responsible for the collision detection system. It holds the
 * simulation domain (a finite axis-aligned box) and a vector which holds the materials inside the voxels.
 *
 * It is allocated by the static {@link create} and {@link destroy} functions.
 * There is a simple straightforward constructor, but no destructor.
 */

template<bool gpu_flag = true>
class voxels
{
public:
	using triangle_index_t = uint32_t; ///< Type for indexing the triangles

	/**
	 * \brief Constructor of the voxels class
	 *
	 * \param vox_size Size of a side of the voxels in nanometer
	 * \param shape Shape of the simulation domain in voxels, as vec3 {size_x, size_y, size_z}
	 * \param initial_geometry The initial geometry of the sample, as a std::vector<int> of length size_x*size_y*size_z
	 */
	voxels(real voxel_size, vec3 shape, std::vector<int> initial_geometry, int max_save_height, real sim_depth);
	/**
	 * \brief Allocate memory for the triangles on the correct device.
	 *
	 * \param triangles List of triangles to be used in the simulation.
	 */
	
	static CPU void destroy(voxels & geometry);

	/**
	 * \brief Check whether a certain position is part of the simulation domain.
	 *
	 * \param position Position to check.
	 */
	inline PHYSICS bool in_domain(vec3 position);

	/**
	 * \brief Try to move a particle, checking for collisions with triangles.
	 *
	 * Try to propagate from \p start to \p start + \p distance * \p direction.
	 * If there is no collision with a triangle, the returned intersect_event
	 * contains a `nullptr` triangle pointer.
	 *
	 * \param start           The particle's starting position
	 * \param direction       Direction the particle is going in
	 * \param distance        Distance the particle travels, in units of direction
	 * \param ignore_triangle Triangle to ignore. This is not used in this voxel version.
	 * \param ignore_material Destination material to ignore (most likely the
	 *                        particle's current material).
	 *
	 * TODO ignore_material datatype
	 */
	inline PHYSICS intersect_event propagate(vec3 start, vec3 direction, real distance,
		triangle const * ignore_triangle, int ignore_material) const;

	/**
	 * \brief Set the material of a certain voxel - Disabled
	 *
	 * Set the material in the voxel that contains position to value material. 
	 *
	 * \param position           Where to set the material
	 * \param material			 The material to which the voxel will be set
	 */

	//inline PHYSICS void set_material(vec3 position, int material, int PE_tag, real energy, uint8_t species, std::deque<char> new_species);

	/**
	 * \brief Get the material of a voxel
	 */
	inline PHYSICS int get_material(int position) const;

	/**
	 * \brief Deposit at position position
	 */
	inline PHYSICS void deposit(vec3 position, vec3 normal, int material, int PE_tag, real energy, uint8_t species, std::deque<char> new_species,vec3 intial_vel);

	/**
	 * \brief Update surface voxels based on new deposition position
	 * 
	 * \param pos indices positions of voxel which should update nieghbours
	 * \param deposition_update indicates whether update is due to deposition or not
	 */
	inline PHYSICS void update_surface_grid(int index, bool deposition_update);

	/**
	 * \brief checks if there are any sites(empty or not) within corner-restricted 3d neighbourhood and returns their grid_indeces
	 * 
	 * \param empty_toggle send back only empty sites(vaccum) or also sites with adsorbates(vaccum/precursor)
	 * \param pos indices positions of voxel which should check nieghbours
	 * 
	 */
	inline PHYSICS std::vector<int> surface_check(int index,bool empty_toggle, bool full_toggle);


	/**
	 * \brief performs a desorption and adsorption step for surface molecules.
	 */
	inline PHYSICS std::vector<int> grid_index_to_vec3(int index);

	/**
	 * \brief performs a desorption and adsorption step for surface molecules.
	 */
	inline PHYSICS void adsorption_routine();

	/**
	 * \brief performs a surface diffision step
	 */
	inline PHYSICS void surface_diffusion_routine();

	/**
	 * \brief performs a surface diffision step
	 */
	inline PHYSICS void surface_diffusion_routine_set();

	/**
	 * \brief reconstructs surface grid from input geometry through scans from +z,+-x
	 * 
	 * \param start_search gives values from of minimum checking range given knowledge about input geometry
	 * \param scan_depth how many material voxels in depth to continue calling surface updates
	 */
	inline PHYSICS void reconstruct_surface_grid(std::vector<int> start_search,int scan_depth);

	/**
	 * \brief Get the maximum distance that can be travelled inside the
	 *        simulation domain.
	 */
	inline PHYSICS real get_max_extent() const;

	/**
	 * \brief Get the maximum distance that can be travelled inside the
	 *        simulation domain.
	 */
	inline PHYSICS void generate_neighbour_incs();

	/**
	 * \brief Get the maximum distance that can be travelled inside the
	 *        simulation domain.
	 */
	inline PHYSICS real get_total_voxels() const;

		/**
	 * \brief Get the amount of electrons between gas routines
	 */
	inline PHYSICS int get_particles_between_gas_routines() const;

	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline PHYSICS vec3 AABB_min() const;

	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline PHYSICS vec3 AABB_max() const;

	/**
	 * \brief Get the adsorption toggle
	 */
	inline PHYSICS bool get_adsorption();

	/**
	 * \brief Get the surface diffusion toggle
	 */
	inline PHYSICS bool get_surface_diffusion();

	/**
	 * \brief Save grids into binary for plotting and analysis:
	 * 
	 */
	
	inline PHYSICS real get_sim_depth() const;

	/**
	 * \brief Save geometry grids and other parameters in binary format
	 * 
	 */
	inline CPU void save_bin_output(std::string file_name);

	/**
	 * \brief Set random generator state for geometry
	 */
	inline CPU void set_rand_state(util::random_generator<false> rand_state_in);

	/**
	 * \brief Get save surface grid in binary format
	 * 
	 */
	inline CPU void save_bin_surface(std::string file_name);

public:
	float sample_height_ext;
private:
	CPU void set_AABB(vec3 min, vec3 max);

	// These vectors are implicit 3D vectors, that represent a voxel grid
	std::vector<int> _mat_grid; // a voxel grid with nebula material codes. 
	std::vector<int> _tag_grid; // .. with PE tags
	std::vector<real> _e_grid; // .. with dissociation energies
	std::vector<int> _species_grid;  // .. with electron species 
	std::vector<int> _new_species_grid; // .. with new electron species --New Internal Classification


	// Variables Gas Handling - (Port Outside of voxels.h at some point to own module)
	std::vector<int> _surface_grid; //tracks deposition surface and adsorbate sites 0 --> none 1 --> surface voxels 2--> open sites 3 -- closed sites 4--> sources
	std::vector<int> _residence_grid; //keeps track of residence times >0 means adsorbate has remaining residence time.
	std::vector<int> _open_vec; // unordered vector of open site indices
	std::vector<int> _residence_vec; //unorderd vector of full site indices
	std::map<int,int> _track_map;// position, tracking_id
	std::vector<int> _track_vec; // tracks distance travelled, index is tracking_id
	std::map<int,int> _start_map; //id,start_pos

	//Gas Handling booleans and variables
	int _particles_between_gas_routines;
	real _time_between_gas_routines; //time between gas handling steps
	real _flux_constant; //flux constant(product with number of surface sites results in amount of adsorbed molecules per electron)
	real _jumps_routine; //total jumps per timestep for surface diffusion
	real _jumps; //jumps per timestep
	bool _jump_restriction;
	real _sticking_c; //stick coefficient determining sucess of adsorption events
	real _mean_residence_time; // mean residence time for adsorption
	util::random_generator<false> _rand_state;
	bool _desorption; 
	bool _adsorption;
	bool _gas_material_toggle;
	bool _surface_diffusion;
	bool _bulk_mirrors;
	bool _bulk_transport;
	bool _side_source;
	bool _central_source;
	bool _perm_source;
	bool _smith_source;
	bool _smith_seq;
	bool _surface_track;
	int _start_amount;
	bool _surface_update_override;
	bool _diffusion_optimization;
	bool _diffusion_concentric;
	real _adsorbs;
	int _corner_restricted_values_array[18][3] =   {{1,0,1},{0,1,1},{-1,0,1},{0,-1,1},{0,0,1},
													{0,1,0},{1,0,0},{1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{0,-1,0},{-1,0,0},
													{1,0,-1},{0,1,-1},{-1,0,-1},{0,-1,-1},{0,0,-1}}; //corner restricted array 
	std::vector<int> _nearest_neighbour_incs;


	//General Variables
	real _voxel_size; // voxel size in nm
	int _size_x; // simulation domain size in x direction in voxels
	int _size_y; // .. y direction
	int _size_z; // .. z direction
	int _sample_height; // .. used by nebula dynamic construction
	real _sim_depth; // .. used by nebula dynamic construction
	int _total_voxels; // integer giving total voxel amount
	int _min_save_height = 0; // voxels with z-index < _min_save_height will not be saved
	int _max_save_height; //  voxels with z-index > _max_save_height will not be saved
	vec3 _AABB_min       = { 0, 0, 0 };
	vec3 _AABB_max       = { 0, 0, 0 };
	real _max_extent     = 0;
	

	friend struct detail::voxels_factory<gpu_flag>;
};

}} // namespace nbl::geometry

#include "voxels.inl"

#endif // __GEOMETRY_VOXELS_H_
