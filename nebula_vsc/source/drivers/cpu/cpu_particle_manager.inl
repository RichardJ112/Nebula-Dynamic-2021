#include <algorithm>
#include <deque>

namespace nbl { namespace drivers {

template<typename material_manager_t>
cpu_particle_manager<material_manager_t>
	cpu_particle_manager<material_manager_t>::create()
{
	cpu_particle_manager<material_manager_t> manager{};
	return manager;
}

template<typename material_manager_t>
void cpu_particle_manager<material_manager_t>::destroy(
	cpu_particle_manager<material_manager_t> & manager)
{
	manager.particles.clear();
	manager.particles.shrink_to_fit();
	manager.cascades.clear();
}

template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::push(
	particle* primary_particles, primary_tag_t* tags, particle_index_t N)
-> particle_index_t
{
	particles.reserve(particles.size() + N);
	for (particle_index_t i = 0; i < N; ++i)
	{
		particles.push_back({
			NO_EVENT,
			0,
			-123,  // TODO: vacuum
			primary_particles[i],
			tags[i],
			0,
			nullptr,
			0,
			{0,'V','V','F'} //new_species
		});

		cascades.insert(std::make_pair(tags[i], cascade_struct{1, 1}));
	}
	return N;
}

template<typename material_manager_t>
template<typename detect_function>
void cpu_particle_manager<material_manager_t>::flush_detected(detect_function func)
{
	for (auto& this_particle : particles)
	{
		if (this_particle.status == DETECTED)
		{
			func(this_particle.particle_data, this_particle.primary_tag);
			this_particle.status = TERMINATED;
		}
	}
}

template<typename material_manager_t>
void cpu_particle_manager<material_manager_t>::flush_terminated()
{
	particles.erase
	(
		std::remove_if(particles.begin(), particles.end(),
			[](particle_struct const & x) -> bool
			{ return x.status == TERMINATED; }),
		particles.end()
	);

	// C++20: can use std::erase_if here.
	for (auto i = cascades.begin(); i != cascades.end(); )
	{
		if (i->second.running_count == 0)
			i = cascades.erase(i);
		else
			++i;
	}
}

template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::get_total_count() const -> particle_index_t
{
	return particles.size();
}
template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::get_running_count() const -> particle_index_t
{
	return static_cast<particle_index_t>(
	std::count_if(particles.begin(), particles.end(),
		[](particle_struct const & x) -> bool
		{ return x.status != TERMINATED && x.status != DETECTED; }));
}
template<typename material_manager_t>
auto cpu_particle_manager<material_manager_t>::get_detected_count() const -> particle_index_t
{
	return std::count_if(particles.begin(), particles.end(),
		[](particle_struct const & x) -> bool
		{ return x.status == DETECTED; });
}

template<typename material_manager_t>
PHYSICS particle & cpu_particle_manager<material_manager_t>::operator[](particle_index_t i)
{
	return particles[i].particle_data;
}
template<typename material_manager_t>
PHYSICS particle const & cpu_particle_manager<material_manager_t>::operator[](particle_index_t i) const
{
	return particles[i].particle_data;
}

template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::exists(particle_index_t i) const
{
	return i < particles.size();
}

template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::active(
	particle_index_t i) const
{
	switch (particles[i].status)
	{
		case TERMINATED:
		case DETECTED:
			return false;
		default:
			return true;
	}
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::get_material_index(particle_index_t i) const
-> material_index_t
{
	return particles[i].current_material;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_material_index(
	particle_index_t particle_idx, material_index_t new_material_idx)
{
	particles[particle_idx].current_material = new_material_idx;
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::get_primary_tag(particle_index_t i) const
-> primary_tag_t
{
	return particles[i].primary_tag;
}

//Added secondary electron return value
template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::get_secondary_tag(particle_index_t i) const
-> primary_tag_t
{
	return particles[i].secondary_tag;
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::get_species(particle_index_t i)
-> uint8_t
{
	if(particles[i].species != 0)
	{
		return particles[i].species;
	}

	if(particles[i].particle_data.dir.z > 0.999)
	{
		return 0; // PE
	}
	if(particles[i].particle_data.dir.z > 0)
	{
		return 1; // FSE
	}
	return 2; // BSE
	
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::set_species(particle_index_t i, uint8_t species)
-> void
{
	particles[i].species = species;
}

template<typename material_manager_t>
PHYSICS triangle const * cpu_particle_manager<material_manager_t>::get_last_triangle(
	particle_index_t i) const
{
	return particles[i].last_triangle;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::forget_last_triangle(
	particle_index_t i)
{
	particles[i].last_triangle = nullptr;
}

template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::next_scatter(particle_index_t i) const
{
	return particles[i].status == SCATTER_EVENT;
}
template<typename material_manager_t>
PHYSICS uint8_t cpu_particle_manager<material_manager_t>::get_next_scatter(particle_index_t i) const
{
	return particles[i].next_scatter;
}
template<typename material_manager_t>
PHYSICS bool cpu_particle_manager<material_manager_t>::next_intersect(particle_index_t i) const
{
	return particles[i].status == INTERSECT_EVENT;
}

template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::create_secondary(
	particle_index_t primary_idx, particle secondary_particle)
{
	uint8_t species;

	particle PE = particles[primary_idx].particle_data;

	if (particles[primary_idx].species == 0)
	{
		if (PE.dir.z > 0.999)
		{
			species = 4; // SE_PE
		}
		else if (PE.dir.z > 0)
		{
			species = 5; // SE_FSE
		}
		else
		{
			species = 6; // SE_BSE
		}
	}
	else if (particles[primary_idx].species == 3)
	{
		species = 7; //SE_VE
	}
	else
	{
		species = particles[primary_idx].species;
	}


	//New Species Section
	
	//Character Conversion
	char dir_char;
	char mat_char;
	if (secondary_particle.dir.z>0)
	{
		dir_char = 'F';
	}
	else
	{
		dir_char = 'B';
	}

	if (secondary_particle.pos.z > sample_height)
	{
		mat_char = 'S';
	}
	else
	{
		mat_char = 'D';
	}
	const auto primary_tag = particles[primary_idx].primary_tag;

	// Origination Tags
	std::deque<char> primary_chars = particles[primary_idx].new_species;
	char origin_tag;
	std::vector<char> primary_chars_check_1 = {0,'V','V','F'};
	// Start with simple PE,FSE,BSE Distinction (1,2,3)
	if (primary_chars[1] == primary_chars_check_1[1] && primary_chars[2] == primary_chars_check_1[2])
	{
		origin_tag = 1; //PE
	}
	else if (primary_chars.back() == 'F')
	{
		origin_tag = 2; //FSE
	}
	else if (secondary_particle.pos.z < sample_height)
	{
		origin_tag = 3; //BSE_D
	}
	else
	{
		origin_tag = 4; //BSE_S
	}

	char new_species[] = {origin_tag,mat_char,mat_char,dir_char};

	particles.push_back({
		NO_EVENT,
		0,
		get_material_index(primary_idx),
		secondary_particle,
		primary_tag,
		cascades[primary_tag].next_secondary_tag++,
		nullptr,
		species,
		{new_species[0],new_species[1],new_species[2],new_species[3]}
	});

	//New Species

	++cascades[primary_tag].running_count;
}

template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::terminate(particle_index_t i)
{
	particles[i].status = TERMINATED;
	--cascades[particles[i].primary_tag].running_count;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::detect(particle_index_t i)
{
	particles[i].status = DETECTED;
	--cascades[particles[i].primary_tag].running_count;
}

// TODO: the two functions below recalculate the normalization of "dir"
// which has already been done by the driver...
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_scatter_event(
	particle_index_t i, scatter_event event)
{
	if (event.type != 0)
	{
		particles[i].status = SCATTER_EVENT;
		particles[i].next_scatter = event.type;
	}
	else
	{
		particles[i].status = NO_EVENT;
	}
	particles[i].particle_data.pos += normalised(particles[i].particle_data.dir) * event.distance;
}
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_intersect_event(
	particle_index_t i, intersect_event event)
{
	particles[i].status = INTERSECT_EVENT;
	particles[i].last_triangle = event.isect_triangle;
	particles[i].particle_data.pos += normalised(particles[i].particle_data.dir) * event.isect_distance;
}
//Cascade_Diagram Functions

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::update_cascade_diagram(particle_index_t i)
-> void
{
	cascade_diagram.primaries.push_back(particles[i].primary_tag);
	cascade_diagram.secondaries.push_back(particles[i].secondary_tag);
	cascade_diagram.x_pos.push_back(particles[i].particle_data.pos.x);
	cascade_diagram.y_pos.push_back(particles[i].particle_data.pos.y);
	cascade_diagram.z_pos.push_back(particles[i].particle_data.pos.z);
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::output_cascade_diagram(const std::string file_name)
-> void
{
	//Length 
	std::ofstream output_bin_file(file_name, std::ios::binary);
	int64_t len = cascade_diagram.primaries.size();
	output_bin_file.write( (char*)&len, sizeof(len) );

	//Vectors
	output_bin_file.write( (char*)&cascade_diagram.primaries[0], len * sizeof(uint32_t) );
	output_bin_file.write( (char*)&cascade_diagram.secondaries[0], len * sizeof(uint32_t) );
	output_bin_file.write( (char*)&cascade_diagram.x_pos[0], len * sizeof(float_t) );
	output_bin_file.write( (char*)&cascade_diagram.y_pos[0], len * sizeof(float_t) );
	output_bin_file.write( (char*)&cascade_diagram.z_pos[0], len * sizeof(float_t) );

    output_bin_file.close();
}
// Boundary Energy Distribution

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::update_boundary_energy_dist(particle_index_t i,real energy)
-> void
{
	boundary_energy_dist.boundary_energies.push_back(energy);
	boundary_energy_dist.x_pos.push_back(particles[i].particle_data.pos.x);
	boundary_energy_dist.y_pos.push_back(particles[i].particle_data.pos.y);
	boundary_energy_dist.z_pos.push_back(particles[i].particle_data.pos.z);
}

template<typename material_manager_t>
PHYSICS auto cpu_particle_manager<material_manager_t>::output_boundary_energy_dist(const std::string file_name)
-> void
{
	//Length 
	std::ofstream output_bin_file(file_name, std::ios::binary);
	int64_t len = boundary_energy_dist.boundary_energies.size();
	output_bin_file.write( (char*)&len, sizeof(len) );

	//Vectors
	output_bin_file.write( (char*)&boundary_energy_dist.boundary_energies[0], len * sizeof(float_t) );
	output_bin_file.write( (char*)&boundary_energy_dist.x_pos[0], len * sizeof(float_t) );
	output_bin_file.write( (char*)&boundary_energy_dist.y_pos[0], len * sizeof(float_t) );
	output_bin_file.write( (char*)&boundary_energy_dist.z_pos[0], len * sizeof(float_t) );

    output_bin_file.close();
}


//New Classification System Functions

//New Species Formatting --> using seperate functions and adding functionality
template<typename material_manager_t>
PHYSICS std::deque<char> cpu_particle_manager<material_manager_t>::get_new_species(particle_index_t i) 
{
		return particles[i].new_species;
}

//Update existing event_chars for species classification
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::update_new_species(particle_index_t i,int secondary_tag,int material_in,int material_out,float z_dir,float z_pos)
{
	int event_length = 2; //current history
	char material_char_in;
	char material_char_out;
	char dir_char;

	// Origination Value
	char origin_tag = particles[i].new_species[0];
	//Keep event_length within memory bounds
	while (particles[i].new_species.size() > (event_length-1)*4 ) // 8 bytes: 1 byte per char. 
	{ 
		particles[i].new_species.pop_front();
	}
	//Character Conversion
	switch (material_in) 
		{
		case -123: //vacuum,
			material_char_in = 'V';
			break;
		case 1:
			material_char_in = 'D';
			break; //full sites always part of subtrate level
		case -122: //Mirror
				if (z_pos > sample_height)
			{
				material_char_in = 'S';
				break;
			}
			else
			{
				material_char_in = 'D';
				break;
			}
		case -126: //detector
			material_char_in = 'S';
			break;
		case 0:
			if (z_pos > sample_height)
			{
				material_char_in = 'S';
				break;
			}
			else
			{
				material_char_in = 'D';
				break;
			}
		default:
			throw std::runtime_error("Invalid Material");
		}
	
	switch (material_out) 
		{
		case -123:
			material_char_out = 'V';
			break;
		case 1:
			material_char_out = 'D'; 
			break;
		case -122: //Mirror
				if (z_pos > sample_height)
			{
				material_char_in = 'S';
				break;
			}
			else
			{
				material_char_in = 'D';
				break;
			}
		case -126: //detector
			material_char_in = 'S';
			break;
		case 0:
			if (z_pos > sample_height) // +0.3(size of 1 voxel) Ad-hoc solution for now 
			{
				material_char_out = 'S'; 
				break;
			}
			else
			{
				material_char_out = 'D';
				break;
			}
		default:
			throw std::runtime_error("Invalid Material");
		}

		if (z_dir > 0)
		{
			dir_char = 'F';
		}
		else
		{
			dir_char = 'B';
		}
	

	char event_chars[] = {origin_tag,material_char_in,material_char_out,dir_char}; 
	for (auto& event_char : event_chars)
	{
		particles[i].new_species.push_back(event_char);
	}
	
}

//Allows sample height to  be set externally (retains encapsulation)
template<typename material_manager_t>
PHYSICS void cpu_particle_manager<material_manager_t>::set_sample_height(float sample_height_ext)
{
	sample_height = sample_height_ext;
}

}} // namespace nbl::drivers
