// Strucure for storing electron cascades

	struct cascade_list
{
	std::vector<uint32_t> primaries; //primary tags
	std::vector<uint32_t> secondaries; //secondary tags
	std::vector<real> x_pos;
	std::vector<real> y_pos;
	std::vector<real> z_pos;

};

cascade_list cascade_diagram;

// Toggle to enable cascade tracking and output --> false results in empty cascade output file
bool cascade_toggle = false; // seems to cause stops/slowdowns(too many scattering events causes memory issues) --> dont use with large depositions (>10kpp)

namespace nbl { namespace drivers {

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
//remove const requirement for dynamic geometry
cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::cpu_driver(
	intersect_t intersect,
	material_manager_t const & materials,
	geometry_manager_t  & geometry,
	real min_energy, real max_energy,
	seed_t seed
) :
	_min_energy(min_energy), _max_energy(max_energy),
	_particles(particle_manager_t::create()),
	_materials(materials),
	_geometry(geometry),
	_intersect(intersect),
	rand_state(seed)
{
	//distribute rand_state to geometry and import gas handling booleans from voxels.inl --> more streamlined/central options should be created here through a config file
	_geometry.set_rand_state(rand_state);
	_adsorption = _geometry.get_adsorption();
	_surface_diffusion = _geometry.get_surface_diffusion();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::~cpu_driver()
{
	particle_manager_t::destroy(_particles);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::push(
	particle* particles,
	primary_tag_t* tags,
	particle_index_t N
) -> particle_index_t
{
	return _particles.push(particles, tags, N);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::get_running_count() const
-> particle_index_t
{
	return _particles.get_running_count();
}
template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::get_detected_count() const
-> particle_index_t
{
	return _particles.get_detected_count();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
template<typename detect_function>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::flush_detected(
	detect_function func)
{
	_particles.flush_detected(func);
	_particles.flush_terminated();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::init(particle_index_t particle_idx)
{
	// Get data from memory
	auto this_particle = _particles[particle_idx];

	// Gas Handling Routines
	bool gas_handling_toggle = true; //enables/disables gas handling routines here
	int p_bet = _geometry.get_particles_between_gas_routines();
	_current_primary = (int)_particles.get_primary_tag(particle_idx)+1;
	int prim_dif = _current_primary - _previous_primary;
	if ((prim_dif >= p_bet) && (gas_handling_toggle))
	{
		if (_adsorption)
		{
			_geometry.adsorption_routine();
		}
		if (_surface_diffusion)
		{
			_geometry.surface_diffusion_routine(); //vector based
			//_geometry.surface_diffusion_routine_set(); //set based (Depreciated)
		}
		_previous_primary = _current_primary; 
		
	} 
	//Saving electron positions for cascade diagram

	if (cascade_toggle)
	{
		_particles.update_cascade_diagram(particle_idx);
		uint32_t primary = _particles.get_primary_tag(particle_idx);
		uint32_t secondary = _particles.get_secondary_tag(particle_idx);
		cascade_diagram.primaries.push_back(primary);
		cascade_diagram.secondaries.push_back(secondary);
		cascade_diagram.x_pos.push_back(this_particle.pos.x);
		cascade_diagram.y_pos.push_back(this_particle.pos.y);
		cascade_diagram.z_pos.push_back(this_particle.pos.z);
	}

	// If not in domain, terminate
	if (!_geometry.in_domain(this_particle.pos))
	{
		_particles.terminate(particle_idx);
		return;
	}

	// Next scattering event
	scatter_event next_scatter{
		0,
		_geometry.get_max_extent()
	};

	// If not in a vacuum, get next scatter event
	auto this_material_idx = _particles.get_material_index(particle_idx);

	if (_materials.is_physical(this_material_idx))
	{
		const auto this_material = _materials[this_material_idx];

		// Terminate if we are below the energy threshold
		// (which is with respect to the vacuum energy)
		if (this_particle.kin_energy < this_material.barrier + _min_energy ||
			this_particle.kin_energy > this_material.barrier + _max_energy)
		{
			_particles.terminate(particle_idx);
			return;
		}

		// Sample next scattering event
		// TODO: case of no scattering events!
		next_scatter = this_material.sample_path(this_particle, rand_state);
	}

	// Move particle to next event, unless there is a triangle in the way
	normalise(this_particle.dir);
	intersect_event next_intersect = _geometry.propagate(
		this_particle.pos, this_particle.dir, next_scatter.distance,
		_particles.get_last_triangle(particle_idx),
		_particles.get_material_index(particle_idx)
	);

	if (next_intersect.isect_triangle == nullptr)
	{
		// No triangle intersection: move to scattering position.
		// Scatter there later (after sorting)
		_particles.set_scatter_event(particle_idx, next_scatter);
	}
	else
	{
		// Triangle intersection: move to triangle position.
		_particles.set_intersect_event(particle_idx, next_intersect);
	} 
}


template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::intersect(particle_index_t particle_idx)
{
	// ignore all particles except those with an intersect event.
	if (!_particles.next_intersect(particle_idx))
		return;

	_intersect.execute(_materials, _particles, particle_idx, rand_state);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::scatter(particle_index_t particle_idx)
{
	// ignore all particles except those with an inelastic event.
	if (!_particles.next_scatter(particle_idx))
		return;

	// forget last intersected triangle. This event might cause us to scatter
	// back into that triangle and we don't want to ignore that triangle if so.
	_particles.forget_last_triangle(particle_idx);

	_materials[_particles.get_material_index(particle_idx)].execute(
		_particles.get_next_scatter(particle_idx), _particles, particle_idx, rand_state);
}

}} // namespace nbl::drivers
