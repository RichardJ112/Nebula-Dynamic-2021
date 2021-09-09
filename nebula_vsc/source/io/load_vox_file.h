#ifndef __LOAD_VOX_FILE_H_
#define __LOAD_VOX_FILE_H_

#include <vector>
#include <string>
#include "../core/triangle.h"
#include "../geometry/voxels.h"

namespace nbl {

	/**
	 * \brief Load a voxel file
	 *
	 * Loading a vox file from a binary file that has all data required to satisfy voxel constructor.
	 *
	 * \return An initialized voxel geometry
	 */
	//nbl::geometry::voxels<false> load_vox_file(std::string const & filename);
	nbl::geometry::voxels<false> load_vox_file(std::string const & filename);
	 

} // namespace nbl

#endif // __LOAD_TRI_FILE_H_
