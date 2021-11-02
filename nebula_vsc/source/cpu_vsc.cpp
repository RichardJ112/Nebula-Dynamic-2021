#include "config/config.h"
#include "physics_config.h"
#include "core/material.h"
#include "common/cli_params.h"
#include "common/work_pool.h"
#include "common/time_log.h"

#include "io/output_stream.h"
#include "io/load_tri_file.h"
#include "io/load_pri_file.h"
#include "io/load_vox_file.h"

#include "drivers/cpu/simple_cpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"
#include "geometry/voxels.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>
#include <ctime>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/stdc++.h>
#include <cstdlib>



// Main typedefs
using geometry_t = nbl::geometry::voxels<false>;
using material_t = material<scatter_physics<false>>;

using driver = nbl::drivers::simple_cpu_driver<
	scatter_physics<false>,
	intersect_t,
	geometry_t
>;

// Get maximal energy accepted by all material files
real get_max_energy(std::vector<material_t> const & materials)
{
	if (materials.size() == 0)
		return 0;

	real max_energy = std::numeric_limits<real>::infinity();
	for (auto const & mat : materials)
	{
		const real e = mat.get_max_energy();
		if (e < max_energy)
			max_energy = e;
	}
	return max_energy;
}

int main(int argc, char** argv)
{
	// Print version information
	std::clog << "This is Nebula 2021 version "
		<< VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\n\n"
		"Physics models:\n";
	scatter_physics<false>::print_info(std::clog);
	intersect_t::print_info(std::clog);

	// Settings
	cli_params p("[options] <geometry.bin> <primaries.pri> [material0.mat] .. [materialN.mat]");
	p.add_option("energy-threshold", "Lowest energy to simulate", 0);
	p.add_option("seed", "Random seed", 0x78c7e39b);
	// Additional Settings NB2021 Deposition Type, Toggle Deposition -Add here in future to circumvent physics config file ect.

	p.parse(argc, argv);
	const real energy_threshold = p.get_flag<real>("energy-threshold");
	const unsigned int seed = p.get_flag<unsigned int>("seed");


	// Setup time logging
	time_log timer;

	// Interpret command-line options
	std::vector<std::string> pos_flags = p.get_positional();
	if (pos_flags.size() < 3)
	{
		p.print_usage(std::clog);
		return 1;
	}

	std::mt19937 random_generator(seed);

	// Load geometry
	//std::clog << "Loading geometry..." << std::endl;
	timer.start();

// Sanity check with number of materials
	{
		int max_material = 1; // in the first version of the dynamic version there is only one material 

		if (max_material > int(pos_flags.size())-3)
		{
			std::clog << "Error: not enough materials provided for this geometry!\n"
				"  Expected " << (max_material+1) << " materials, " << (pos_flags.size()-2) << " provided.\n";
			return 1;
		}
		if (max_material < int(pos_flags.size())-3)
		{
			std::clog << "Warning: too many materials provided for this geometry!\n"
				"  Expected " << (max_material+1) << " materials, " << (pos_flags.size()-2) << " provided.\n";
		}
	}
	
	timer.start();
	
	geometry_t geometry = nbl::load_vox_file(pos_flags[0]); 
	timer.stop("Loading geometry");


	auto surface_diffusion = geometry.get_surface_diffusion(); //gas handling bools currently set in voxels.inl --> move to seperate encapsulated code
	auto adsorption = geometry.get_adsorption();

	std::clog << std::boolalpha <<
	" * Gas Handling Routines\n" //This is temporary positioning --> this needs to be brought elsewhere and integrated
	"	Options:\n"
	"     - Adsorption: " << adsorption << "\n"
	"     - Surface Diffusion: " << surface_diffusion << "\n";
	std::clog << "\n" << std::string(80, '-') << "\n\n";

	//Automatic Naming for Geometry
	std::string geometry_path = pos_flags[0];
	std::size_t botDirPos_geom = geometry_path.find_last_of("/");
	std::string geometry_file = geometry_path.substr(botDirPos_geom, geometry_path.length());
	size_t last_index_geom = geometry_file.find_first_of("."); 
	std::string geometry_str = geometry_file.substr(1, last_index_geom-1);

	//Additional Parameter Information --Use Systematic file names to get certain parameters (should be commented/changed to accommadate custom naming)
	size_t first_index_sd = geometry_str.find("sd");
	size_t first_index_sh = geometry_str.find("sh");
	size_t start_geo_index = geometry_str.find("_",first_index_sh+3);
	int sim_depth_str_length = first_index_sh-first_index_sd-4;
	std::string sim_depth_str = geometry_str.substr(first_index_sd+3, sim_depth_str_length); //assume scale in micrometers thus length of 5 --> >10,000 voxels
	std::string sim_height_str = geometry_str.substr(first_index_sh+3, start_geo_index-first_index_sh-3); //atm we assume length of 3 --> between 100 and 999 voxels(30nm - 300nm)
	std::string geometry_type_str = geometry_str.substr(start_geo_index+1, 50); //atm we assume length of 3 --> between 100 and 999 voxels(30nm - 300nm)

	//Automatic Naming for Electrons
	std::string electron_path = pos_flags[1];
	std::size_t botDirPos_electron = electron_path.find_last_of("/");
	std::string electron_file = electron_path.substr(botDirPos_electron, electron_path.length());
	size_t last_index_electron = electron_file.find_last_of("."); 
	std::string electron_str = electron_file.substr(1, last_index_electron-1);

	//Automatic Naming for Materials
	std::string material_path = pos_flags[2];
	std::size_t botDirPos_mat = material_path.find_last_of("/");
	std::string material_file = material_path.substr(botDirPos_mat, material_path.length());
	size_t last_index_mat = material_file.find_first_of("."); 
	std::string material_str = material_file.substr(1, last_index_mat-1);

	//Output for Plotting and Analysis
	//std::string sim_depth_str = std::to_string(static_cast<int>(geometry.get_sim_depth())); //ad hoc solution

	//Debug Information
	std::clog << "Loaded geometry" << std::endl;
	std::clog << "Geometry Filename: " << geometry_str << ":\n";

	//DateTime Naming for automatic output to relevant folders
	time_t now = time(0);
	tm *ltm = localtime(&now);
	int month_int = 1+ ltm->tm_mon;

	std::string date_string = std::to_string(ltm->tm_mday)+"_"+std::to_string(month_int)+"_"+std::to_string(1900 + ltm->tm_year);
	
	std::string month_name;
	switch(month_int) {
		case 1  :
			month_name = "January";
			break;
		case 2 :
			month_name = "February";
			break;
		case 3 :
			month_name = "March";
			break;
		case 4 :
			month_name = "April";
			break;
		case 5  :
			month_name = "May";
			break;
		case 6 :
			month_name = "June";
			break;
		case 7 :
			month_name = "July";
			break;
		case 8 :
			month_name = "August";
			break;
		case 9  :
			month_name = "September";
			break;
		case 10 :
			month_name = "October";
			break;
		case 11 :
			month_name = "November";
			break;
		case 12 :
			month_name = "December";
			break;
		
	}
	// This document string is for WSL laptop --> change to your own output directory (document_str --> empty str will disable this function)
	//std::string document_str = "/mnt/c/Users/richa/Documents/repos/nebula_test_files/output/"+month_name+"/"+date_string+"/";
	std::string document_str = "/mnt/c/Users/Richard/source/repos/Nebula/nebula_dynamic_2021/nebula_test_files/output/"+month_name+"/"+date_string+"/";	//This is for personal WSL Desktop
	//std::string document_str = "/home/richarddejong/nebula_test_files/output/"+month_name+"/"+date_string+"/"";	//This is for HPC Server

	std::string CS_choice = "_SM_"; //cross section choice --> should be further integrated into run options/boundary intersect funtion, Still manual only for saving--> change physics_config too!

	std::string general_output_str = document_str + electron_str+"_sd_"+sim_depth_str+"_sh_"+sim_height_str+"_"+geometry_type_str+"_"+material_str+CS_choice;
	std::string output_geometry_bin_str = general_output_str+"output.bin";
	std::string output_surface_bin_str = general_output_str+"surface.bin";
	std::string output_cascade_bin_str = general_output_str+"cascade.bin";
	std::string output_energy_dist_bin_str = general_output_str+"energy_dist.bin";
	std::string output_detector_str = general_output_str+"detector.det";


	const int dir_err = mkdir(document_str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1 == dir_err)
	{
		printf("Error creating directory!n");
	}

// Load materials - Adjusted to Achieve additional
	std::clog << "Loading materials..." << std::endl;
	timer.start();
	nbl::cpu_material_manager<material_t> materials;
	for (size_t parameter_idx = 2; parameter_idx < pos_flags.size(); ++parameter_idx)
	{
		nbl::hdf5_file material(pos_flags[parameter_idx]);
		materials.add(material);

		std::clog << "  Material " << (parameter_idx - 2) << ":\n"
			"    Name: " << material.get_property_string("name") << "\n"
			"    cstool version: " << material.get_property_string("cstool_version") << "\n";
	}
	timer.stop("Loading materials");


	// Load primaries
	std::clog << "\nLoading primary electrons..." << std::endl;
	timer.start();
	std::vector<particle> primaries; std::vector<int2> pixels;
	std::tie(primaries, pixels) = nbl::load_pri_file(pos_flags[1], geometry.AABB_min(), geometry.AABB_max(), materials.get_max_energy());
	timer.stop("Loading primary electrons");

	if (primaries.empty())
	{
		std::clog << "Error: could not load primary electrons!\n";
		p.print_usage(std::clog);
		return 1;
	}

	// The driver only accepts uint32 tags. So we make a map: simulation tag is
	// the index of the primary particle in the "primaries" / "pixels" array.
	std::vector<uint32_t> gpu_tags(primaries.size());
	std::iota(gpu_tags.begin(), gpu_tags.end(), 0); // Fill with 0, 1, ... tags.size()-1

	// This manages the work to be done (thread-safe).
	work_pool pool(primaries.data(), gpu_tags.data(), primaries.size());

	intersect_t inter{ &geometry};

	// Print debug data
	std::clog << "\n"
		<< "Primary Filename: " << electron_str << "\n"
		<< "Loaded " << geometry.get_total_voxels() << " voxels.\n"
		<< "  min = {" << geometry.AABB_min().x << ", " << geometry.AABB_min().y << ", " << geometry.AABB_min().z << "}\n"
		<< "  max = {" << geometry.AABB_max().x << ", " << geometry.AABB_max().y << ", " << geometry.AABB_max().z << "}\n"
		<< "Loaded " << primaries.size() << " primaries.\n"
		<< "Loaded " << materials.size() << " materials.\n\n" << std::flush;

	// Prepare output file
	output_stream out_file(output_detector_str); //changed to detector string

	// Simulation loop
	auto sim_loop = [&pool, &out_file, &pixels,
		&geometry, &inter, &materials,&output_cascade_bin_str,&output_energy_dist_bin_str, energy_threshold](uint32_t seed)
	{
		driver d(
			inter, materials, geometry,
			energy_threshold, materials.get_max_energy(), seed);
		output_buffer buff(out_file, 1024*(7*sizeof(float) + 2*sizeof(int)));

		for (;;)
		{
			// Push new particles
			auto work_data = pool.get_work(1);

			if (std::get<2>(work_data) == 0)
				break;

			d.push(
				std::get<0>(work_data),  // particle*
				std::get<1>(work_data),  // tag*
				std::get<2>(work_data)); // number

			// Simulate a little
			d.simulate_to_end();

			// Flush output data
			d.flush_detected([&buff,&pixels](particle p, uint32_t t)
			{
			buff.add(std::array<float, 7>{
					static_cast<float>(p.pos.x), static_cast<float>(p.pos.y), static_cast<float>(p.pos.z),
						static_cast<float>(p.dir.x), static_cast<float>(p.dir.y), static_cast<float>(p.dir.z), static_cast<float>(p.kin_energy)});
				buff.add(std::array<int, 2>{
					pixels[t].x, pixels[t].y});
			});
		}
		d.output_cascade(output_cascade_bin_str); //work here
		d.output_boundary_energy_dist(output_energy_dist_bin_str);

		buff.flush();
	};

	// Simulation
	bool multithreading = false; //only toggle on when no deposition is taking place otherwise things will break 
	auto n_threads = std::thread::hardware_concurrency();
	if (multithreading)
	{
		 n_threads = std::thread::hardware_concurrency();
	}
	else 
	{
		 n_threads = 1;
	}
	
	std::clog << "Creating " << n_threads << " CPU drivers" << std::endl;
	std::vector<std::thread> threads;

	timer.start();
	for (unsigned int i = 0; i < n_threads; ++i)
		threads.push_back(std::thread(sim_loop, random_generator()));

	// Progress indicator
	unsigned long long pgpto = 0; // for testing purposes only
	for (;;)
	{
		std::this_thread::sleep_for(std::chrono::seconds(1));
		auto primaries_to_go = pool.get_primaries_to_go();
		/*
		if(pgpto == primaries_to_go && pgpto != 0)
		{
			std::clog << " This is slow :( ";
		}
		*/
		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size())) << "%";
		if (primaries_to_go == 0)
			break;
		pgpto = primaries_to_go;
	}
	timer.stop("Simulation");
	std::clog << "\nSaving Files..." << std::endl;
	timer.start();

	geometry.save_bin_output(output_geometry_bin_str); // save the final geometry and grids for plotting (bin format)\
	//new output for surface information (turn off when not using diffusion/)
	geometry.save_bin_surface(output_surface_bin_str); //saves the output surface information to verify surface functionality
	
	
	for (auto& t : threads)
		t.join();

	timer.stop("Save file(s)");

	std::clog << "\n\n";
	timer.print(std::clog);
	return 0;
}
