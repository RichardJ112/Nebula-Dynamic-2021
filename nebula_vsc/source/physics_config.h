#ifndef __PHYSICS_CONFIG_H_
#define __PHYSICS_CONFIG_H_

#include <iostream>
#include "core/particle.h"
#include "common/util/table_1D.h"
#include "common/util/table_2D.h"
#include "common/util/table_3D.h"
#include "common/util/random.h"
#include "core/material.h"
#include "physics/kieft/inelastic.h"
#include "physics/kieft/elastic.h"
#include "physics/full_penn.h"
#include "physics/boundary_intersect.h"

/*
 * Physics definitions below
 */

// Penn inelastic model
template<bool gpu_flag>
using inelastic_scatter = nbl::scatter::full_penn<gpu_flag,
	true  // Generate secondaries
>;

// Kieft & Bosch inelastic model
//template<bool gpu_flag>
//using inelastic_scatter = nbl::scatter::kieft_inelastic<gpu_flag,
//	true, // Optical phonon loss
//	true, // Generate secondaries
//	true, // Random instantaneous momentum for secondary
//	true  // Momentum conservation
//>;

// Kieft & Bosch elastic model
template<bool gpu_flag>
using elastic_scatter = nbl::scatter::kieft_elastic<gpu_flag,
	true, // Acoustic phonon loss
	true  // Atomic recoil loss
>;

// Material boundary intersection
// Cross-Section ID's:
// 0 - Alman, 1 - Winters, 2- Alman-Winters, 3 - Smith
using intersect_t = boundary_intersect<
	true, // Quantum-mechanical (probabilistic) transmission
	true, // Interface refraction
	false, // Kieft & Bosch empirical interface absorption
	true, // Deposition
	true, //Deposition Model(No Gas [ 2-Step True 1-Step False])/( Gas - [Bulk to Adsorbate True/ Only vacuum intersect False])
	3, // Dissociation Cross-Section (0 -Alman 1- Winters 2-Alman-Winters(AW) 3-Smith(Dingemanse) 4-Smith(Comparitive Fit) 5- AW-Me)
	false, //electron termination on deposition
	false // Saving Material outbound particle energies
	//true // using material 1 for adsorbates
>;

/*
// Diffusion and Adsorption Model -- Port functionality away from voxels.inl in future versions
using gas_handling = struct gas_parameters<
	false, //Adsorption Routine
	false //Surface Diffusion
>;
*/

// Putting it all together
template<bool gpu_flag>
using scatter_physics = scatter_list<
	inelastic_scatter<gpu_flag>,
	elastic_scatter<gpu_flag>
>;

#endif
