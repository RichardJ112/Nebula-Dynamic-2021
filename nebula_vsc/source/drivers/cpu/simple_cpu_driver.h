#ifndef __SIMPLE_CPU_DRIVER_H_
#define __SIMPLE_CPU_DRIVER_H_

#include "cpu_driver.h"
#include <iostream> 
#include <thread> 

namespace nbl { namespace drivers {

/**
 * \brief CPU driver that only tries to perform the simulation and nothing else.
 * \see cpu_driver
 */

template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
class simple_cpu_driver
	: public cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>
{
public:
	using base_t = cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>;

	using typename base_t::particle_index_t;
	using typename base_t::material_t;
	using typename base_t::material_manager_t;
	using typename base_t::seed_t;

	/// Constructor - No const geometry
	simple_cpu_driver(
		intersect_t intersect,
		material_manager_t const & materials,
		geometry_manager_t& geometry,
		real min_energy, real max_energy,
		seed_t seed = util::random_generator<false>::default_seed)
	: base_t(intersect, materials, geometry, min_energy, max_energy, seed)
	{}

	/// Perform a single iteration of the simulation for all particles.
	void do_iteration()
	{
		const auto particle_count = this->_particles.get_total_count();
		for (particle_index_t particle_idx = 0; particle_idx < particle_count; ++particle_idx)
		{
			if (!this->_particles.active(particle_idx))
				continue;

			this->init(particle_idx);
			this->intersect(particle_idx);
			this->scatter(particle_idx);
		}
	}

	/// Keep simulating until there are no particles left.
	void simulate_to_end()
	{
		this->_particles.set_sample_height(this-> _geometry.sample_height_ext);
		for (particle_index_t particle_idx = 0;
			particle_idx < this->_particles.get_total_count(); ++particle_idx)
		{
			while (this->_particles.active(particle_idx))
			{
				this->init(particle_idx);
				this->intersect(particle_idx);
				this->scatter(particle_idx);
			}
		}
	}

	//New Function for exporting cascades
	void output_cascade(const std::string file_name)
	{
	this ->_particles.output_cascade_diagram(file_name); //offloads functionality to cpu_particle_manager
	}

	//New Function for exporting Boundary Energy Distribution
	void output_boundary_energy_dist(const std::string file_name)
	{
	this ->_particles.output_boundary_energy_dist(file_name); //offloads functionality to cpu_particle_manager
	}

};

}} // namespace nbl::drivers

#endif // __SIMPLE_CPU_DRIVER_H_
