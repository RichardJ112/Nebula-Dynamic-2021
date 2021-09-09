#include "../config/config.h"
#include "load_vox_file.h"
#include "../geometry/voxels.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>


namespace nbl {

	geometry::voxels<false> load_vox_file(std::string const & filename)
	{
		// IO vars
			real voxel_size_io; //Changing to double for IO
			int size_x_io;
			int size_y_io;
			int size_z_io;
			real sim_depth_io;
			int sample_height_io;
			
		// Binary Implementation

		geometry_bin_io input_bin_io;

		std::ifstream input_file(filename, std::ios::binary);
		int64_t len;
		//Reading Binary
		input_file.read( (char*)&len, sizeof(len) ); //vector length
		input_bin_io.ini_geom_io.resize(len);

		//Parameters
		input_file.read( (char*)&input_bin_io.voxel_size_io, sizeof(real) );
		input_file.read( (char*)&input_bin_io.size_x_io, sizeof(int32_t) );
		input_file.read( (char*)&input_bin_io.size_y_io, sizeof(int32_t) );
		input_file.read( (char*)&input_bin_io.size_z_io, sizeof(int32_t) );
		input_file.read( (char*)&input_bin_io.sim_depth_io, sizeof(real) );

		input_file.read( (char*)&input_bin_io.sample_height_io, sizeof(int32_t) );
		//Vectors
		input_file.read( (char*)&input_bin_io.ini_geom_io[0], len * sizeof(int16_t) );
		input_file.close();

		//Read
		voxel_size_io = input_bin_io.voxel_size_io;
		size_x_io = input_bin_io.size_x_io;
		size_y_io = input_bin_io.size_y_io;
		size_z_io = input_bin_io.size_z_io;
		sim_depth_io = input_bin_io.sim_depth_io;
		sample_height_io = input_bin_io.sample_height_io;
		std::vector<int> ini_geom_io(input_bin_io.ini_geom_io.begin(), input_bin_io.ini_geom_io.end());

		// Value transfter to allow for const de_io }; clarations
		vec3 shape_io = { (float)size_x_io, (float)size_y_io, (float)size_z_io};

		//Creation Geometry Structure
		geometry::voxels<false> geometry(voxel_size_io, shape_io, ini_geom_io, sample_height_io + 1, sim_depth_io);

		//Debug saving of initial geometry, binary
		//geometry.save_bin_io("input_io.bin");
		

		return geometry;
	}

} // namespace nbl
