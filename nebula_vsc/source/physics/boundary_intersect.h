#ifndef __BOUNDARY_INTERSECT_H_
#define __BOUNDARY_INTERSECT_H_

#include "../core/triangle.h"
#include "../core/cpu_material_manager.h"
#include "../geometry/voxels.h"
#include <cmath>

/**
 * \brief Material boundary intersection.
 *
 * See  (Kieft paper)
 *
 * \tparam quantum_transmission Consider probabilistic transmission
 * \tparam interface_refraction Consider refraction at interface
 * \tparam interface_absorption Empirical interface absorption (Kieft et al.
 *                                  doi:10.1088/0022-3727/41/21/215310)
 * \tparam deposition Consider deposition
 * \tparam deposition_model Consider 2-Step(True) or 1-Step(False) Deposition
 * \tparam cross_id Consider certain dissociation cross-sections
 * \tparam electron_termination_on_deposition Consider electron termination after causing deposition
 * \tparam save_energy_dist Save Energy distribution of outbound electrons
 * \tparam gas_material Consider Adsorbate Gas of another material
 */
template<
	bool quantum_transmission = true,
	bool interface_refraction = true,
	bool interface_absorption = false,
	bool deposition = true,
	bool deposition_model = true,
	int cross_id = 3,
	bool electron_termination_on_deposition = false,
	bool save_energy_dist = true>
struct boundary_intersect
{
	/**
	 * \brief Print diagnostic info
	 */
	nbl::geometry::voxels<false>* geometry;
	
	static void print_info(std::ostream& stream)
	{
		//Cross-Section Model Parsing
		std::string cross_section_str;
		switch(cross_id) {
		case 0  :
			cross_section_str = "Alman";
			break;
		case 1  :
			cross_section_str = "Winters";
			break;
		case 2 :
			cross_section_str = "Alman-Winters";
			break;
		case 3 :
			cross_section_str = "Smith-Dingemanse";
			break;
		case 4 :
			cross_section_str = "Smith-Comparitive";
			break;
		case 5 :
			cross_section_str = "Alman-Winters-MeCpPtMe3";
		}
		//Deposition Type Parsing
		std::string dep_type_str;
		switch(deposition_model) {
		case false  :
			dep_type_str = "1 - Step";
			break;
		case true :
			dep_type_str = "2 - Step";
			break;
		}

		stream << std::boolalpha <<
			" * Material boundary crossing model\n"
			"   Options:\n"
			"	- Quantum mechanical transmission: " << quantum_transmission << "\n"
			"	- Interface refraction: " << interface_refraction << "\n"
			"	- Empirical interface absorption: " << interface_absorption << "\n"
			"	- Deposition: " << deposition << "\n"
			"	- Deposition Type: " << dep_type_str << "\n"
			"	- Dissociation Cross-Section: " << cross_section_str << "\n"
			"	- Electron Termination on Deposition: " << electron_termination_on_deposition << "\n"
			"	- Save Outbound Energy Distribution: " << save_energy_dist << "\n";
	}

	/**
	 * \brief Perform intersection event.
	 */
	template<typename particle_manager, typename material_manager, bool gpu_flag>
	PHYSICS void execute(material_manager& material_mgr,
	                     particle_manager& particle_mgr, typename particle_manager::particle_index_t particle_idx,
	                     nbl::util::random_generator<gpu_flag>& rng) const
	{
		using material_index_t = typename material_manager::material_index_t;
		
		// Get particle data from particle manager
		auto this_particle = particle_mgr[particle_idx];
		//Debugging copy
		auto copy_particle = particle_mgr[particle_idx];

		//const triangle this_triangle = *(particle_mgr.get_last_triangle(particle_idx));
		material_index_t material_idx_in = particle_mgr.get_material_index(particle_idx);

		// Extract data from triangle pointer (code from luc)
		uint64_t isect_id = reinterpret_cast<uint64_t>(particle_mgr.get_last_triangle(particle_idx)); 
		material_index_t material_idx_out = reinterpret_cast<int32_t*>(&isect_id)[0];
		int voxel_side = reinterpret_cast<int32_t*>(&isect_id)[1];

		//Bulk Exception i.e. same material(happens during bulk transport to voxel substrate level and when using asorbate material equivalent to substrate material)
		if (material_idx_in == material_idx_out)
		{
			return;
		}

		// Get particle direction, normal of the intersected triangle
		auto normalised_dir = normalised(this_particle.dir);
		auto normalised_reflection = normalised(this_particle.dir); //extra reflection debug
		vec3 last_triangle_normal;
		
		// New Classification System
		particle_mgr.update_new_species(particle_idx,particle_mgr.get_secondary_tag(particle_idx),material_idx_in,material_idx_out,this_particle.dir.z,this_particle.pos.z);

		// determine the normal of the voxel side using its number in range (1..6)
		switch (voxel_side) 
		{
		case 1:
			last_triangle_normal = { 1, 0, 0 };
			break;

		case 2:
			last_triangle_normal = { -1, 0, 0 };
			break;

		case 3:
			last_triangle_normal = { 0, 1, 0 };
			break;

		case 4:
			last_triangle_normal = { 0, -1, 0 };
			break;

		case 5:
			last_triangle_normal = { 0, 0, 1 };
			break;

		case 6:
			last_triangle_normal = { 0, 0, -1 };
			break;
			
		default:	
			throw std::runtime_error("Invalid voxel side number!!!");
		}

		// Get angle between direction of motion and triangle
		const real cos_theta = dot_product(last_triangle_normal, normalised_dir);

		// manage special cases for electron detection, electron mirrors and terminators.
		//   DETECTOR always detects.
		//   DETECTOR_LT/GE50 detect under certain circumstances.
		//     If not detected, they pass through as if no intersection event has taken place.
		switch (material_idx_out) {
		case nbl::special_materials::DETECTOR:
			particle_mgr.detect(particle_idx);
			return;
		case nbl::special_materials::DETECTOR_LT50:
			if (this_particle.kin_energy < 50)
				particle_mgr.detect(particle_idx);
			return;
		case nbl::special_materials::DETECTOR_GE50:
			if (this_particle.kin_energy >= 50)
				particle_mgr.detect(particle_idx);
			return;
		case nbl::special_materials::TERMINATOR:
			particle_mgr.terminate(particle_idx);
			return;
		case nbl::special_materials::MIRROR:
			//this_particle.dir = normalised_dir - 2*last_triangle_normal*cos_theta;
			normalised_reflection =  normalised_dir - 2*last_triangle_normal*cos_theta;
			this_particle.dir = normalised_reflection;
			//small "nudge" launches electrons along their proposed path for 1nm to prevent infinite reflections in between mirror voxels
			this_particle.pos += this_particle.dir; //1nm nudge back into domain
			particle_mgr[particle_idx] = this_particle;
			return;
		case nbl::special_materials::NOP:
			return;
		//Special case for adosrbate gas(here still given same material properties)
		default:
			break;
		}
		
		// determine the change in energy `dU` (in eV) when passing through the interface
		// see thesis T.V. Eq. 3.136
		real dU = 0;
		if (material_mgr.is_physical(material_idx_in)) {
			dU -= material_mgr[material_idx_in].barrier;
		}
		if (material_mgr.is_physical(material_idx_out)) {
			dU += material_mgr[material_idx_out].barrier;
		}

		// Material Boundary Shift( when using seperate material for gas adsorbates)
		auto adsorption = geometry->get_adsorption();
		auto surface_diffusion = geometry->get_surface_diffusion();
		bool gas_material_toggle = adsorption||surface_diffusion;//for now only activate with adsorption or surface diffusion
		bool deposition_condition;
		bool deposition_condition_1;
		bool deposition_condition_2;

		if (gas_material_toggle)
		{
			deposition_condition_1 = (material_idx_in == 0) && (material_idx_out == 1); //From bulk to adsorbate
			deposition_condition_2 = (material_idx_in == nbl::special_materials::VACUUM) && (material_idx_out == 1); //From vacuum to adsorbate
			deposition_condition = deposition_condition_1 || deposition_condition_2;
		}
		else
		{
			//Without Gas Handling
			deposition_condition = ((material_idx_in == 0) && (material_idx_out == nbl::special_materials::VACUUM)) || ((material_idx_in == nbl::special_materials::VACUUM) && (material_idx_out == 0));
		}

		
		// 1 - Step Ad Hoc Energy Dist
		/*
		if (save_energy_dist && material_idx_out == nbl::special_materials::VACUUM)
		{
			particle_mgr.update_boundary_energy_dist(particle_idx,this_particle.kin_energy);
		}
		*/

		// Alman cross section values
		real E_TH = 3.5; // dissosiation threshold energy in eV(waarden uit(2008))
		real E_MAX = 18; // maximum dissosiation cross - section energy in eV
		real LAMBDA_0 = 77; // lambda_0 in eV
		real SIGMA_MAX = 1; // sigma_max
		real dissociation_energy = 3.5; // The energy that an electron looses when dissociating an precursor molecule.
		real deposition_prob; // the probability of deposition according to a cross-section	

		if (deposition_model) //True - 2-Step, False - 1-Step
		{
			//2 - Step Deposition - First Transmission Check then Depoosition

			// determine transmission probability (only if energy suffices)
			// see thesis T.V. Eq. 3.145
			if (this_particle.kin_energy*cos_theta*cos_theta + dU > 0)
			{
				const real s = sqrtr(1 + dU / (this_particle.kin_energy*cos_theta*cos_theta));
				const real T = (quantum_transmission
					? 4 * s / ((1 + s)*(1 + s))
					: 1);
				if (rng.unit() < T)
				{
					//Kinetic Energy Adjustment due to transmission (Here we take the energy or not)
					
					this_particle.kin_energy += dU;

					if (interface_refraction)
					{
						// determine the angle of refraction
						// see thesis T.V. Eq. 3.139
						this_particle.dir = (normalised_dir - last_triangle_normal * cos_theta)
							+ s * last_triangle_normal * cos_theta;
					}
					// update the current material index and EXIT.
									
					particle_mgr.set_material_index(particle_idx, material_idx_out);

					// Boundary Intersection Energy Tracker only for outgoing electrons (Advise: use detection output instead)
					
					if (save_energy_dist && material_idx_out == nbl::special_materials::VACUUM)
					{
						particle_mgr.update_boundary_energy_dist(particle_idx,this_particle.kin_energy);
					}	
					
					if (deposition) //Deposition toggle is in physics.config
					{
						
						if (deposition_condition)
						{
							
							const real E = this_particle.kin_energy;						

							//Alman
							if (cross_id == 0)
							{
								if (E <= E_TH) 
								{
									deposition_prob = 0;
								}
								else if (E < E_MAX)
								{
									deposition_prob = SIGMA_MAX * (1 - (std::pow(E_MAX - E,2) / std::pow(E_MAX - E_TH, 2)) );
								}
								else
								{
									deposition_prob = SIGMA_MAX * std::exp(-(E - E_MAX) / LAMBDA_0);
								}
							}
							
							// Winters
							if (cross_id == 1)
							{
								const real E_MAX_w = 100;
								if(E < 36.8)
								{
									deposition_prob = 0;
								}
								else
								{
									deposition_prob = 100*std::log(E / 36.7879) / E;
								}
							}

							// Alman-Winters -- Dingemanse
							if (cross_id == 2) 
							{
								if (E <= E_TH)
								{
									deposition_prob = 0;
								}
								else if (E < E_MAX)
								{
									deposition_prob = SIGMA_MAX * (1 - (std::pow(E_MAX - E,2) / std::pow(E_MAX - E_TH, 2)) );
								}
								else
								{
									deposition_prob = E_MAX*std::log(E / 6.621829941) / E;
								}
							}
							
							// Smith this is for WF6 precursor only, new fit must be made for other precursor gasses based on ionization cs
							//Dingemanse Fit
							if (cross_id == 3)
							{
								if (E < 7)
								{
									deposition_prob = 0;
								}
								else if (E < 100)
								{
									deposition_prob = (1208 * (1 - 1 / E) + -1064 * std::pow((1 - 1 / E), 2) + -15.68 * std::log(E) + -797.4 * std::log(E) / E) / E; //strange order of terms logE/E before log E
								}
								else
								{
									deposition_prob = (-16540 * (1 - 1 / E) + 15970 * std::pow((1 - 1 / E), 2) + 108.8 * std::log(E) + 5885 * std::log(E) / E) / E;
								}

								//Smith deposition probability scaling to a maximum of 1
								bool smith_comp_toggle = true;
								if (smith_comp_toggle)
								{
									//Fix Jonathan
									real E_MAX_Smith = 1.0508; //scaling issue with fucntional fit jonathan
									deposition_prob = deposition_prob/E_MAX_Smith;
								}
							}
							//Smith Comparitive Fit (max amplitude ~0.25)
							if (cross_id == 4)
							{
								if (E < 7)
								{
									deposition_prob = 0;
								}
								else if (E < 90)
								{
									deposition_prob = (308.96 * (1 - 1 / E) + -271.84 * std::pow((1 - 1 / E), 2) + -202.75 * std::log(E) / E + -4.386 * std::log(E)) / E; //changed order of variables to reflect literature equation
								}
								else
								{
									deposition_prob = (-676.83 * (1 - 1 / E) + 577.017 * std::pow((1 - 1 / E), 2) + 410.11 * std::log(E) / E + 20.77 * std::log(E)) / E;
								}
							}
							
							// Alman-Winters -- MeCpPtMe3
						
							if (cross_id == 5) 
							{
								//Change Values
								E_TH = 1.7;
								E_MAX = 18.6;
								SIGMA_MAX = 1;
								LAMBDA_0 = 77;
								dissociation_energy = 1.7;
								real constant_e = 2.718282/E_MAX;
								if (E <= E_TH)
								{
									deposition_prob = 0;
								}
								else if (E < E_MAX)
								{
									deposition_prob = SIGMA_MAX * (1 - (std::pow(E_MAX - E,2) / std::pow(E_MAX - E_TH, 2)) );
								}
								else
								{
									deposition_prob = E_MAX*std::log(E*constant_e) / E;
								}
							}

							if (rng.unit() < deposition_prob)
							{
								geometry->deposit(this_particle.pos, last_triangle_normal, 0, particle_mgr.get_primary_tag(particle_idx), this_particle.kin_energy, particle_mgr.get_species(particle_idx),particle_mgr.get_new_species(particle_idx),copy_particle.dir);
								this_particle.kin_energy -= dissociation_energy;

								//particle_mgr.update_cascade_diagram(particle_idx); //Cascade Diagram Update

								//Termination Toggle to prevent termination
								if (electron_termination_on_deposition)
								{
									particle_mgr.terminate(particle_idx); // After a deposition, the electron is terminated 
								}

							}
						}
					}
					particle_mgr[particle_idx] = this_particle;
					/*
					//Not necessary (fixed issue) but you can use a breakpoint here to check for negative energies
					if (this_particle.kin_energy<0)
					{
						particle_mgr.terminate(particle_idx);
					}
					*/
					return;
				}
			}
		}
		else
		{
			//1-Step model --> deposition first then worry about transmission and internal reflection
			bool e_termination = false;
			/*
			// Boundary Intersection Energy Tracker only for outgoing electrons --No transmission needed
			if (save_energy_dist && material_idx_out == nbl::special_materials::VACUUM)
			{
				particle_mgr.update_boundary_energy_dist(particle_idx,this_particle.kin_energy);
			}
			*/

			if (deposition)
			{
				// deposit a voxel of material 0
				if (deposition_condition)
				{
					const real E = this_particle.kin_energy;

					// Take Dingemanse Fit for now -- I reccomend to refactor code if more research towards 1/2 step deposition models is wanted.
					if (E < 7)
					{
						deposition_prob = 0;
					}
					else if (E < 100)
					{
						deposition_prob = (1208 * (1 - 1 / E) + -1064 * std::pow((1 - 1 / E), 2) + -15.68 * std::log(E) + -797.4 * std::log(E) / E) / E; //strange order of terms logE/E before log E
					}
					else
					{
						deposition_prob = (-16540 * (1 - 1 / E) + 15970 * std::pow((1 - 1 / E), 2) + 108.8 * std::log(E) + 5885 * std::log(E) / E) / E;
					}

					//Smith deposition probability scaling to a maximum of 1
					bool smith_comp_toggle = true;
					if (smith_comp_toggle)
					{
						//Fix Jonathan
						real E_MAX_Smith = 1.0508; //scaling issue with functional fit jonathan --> found maximum value in external script
						deposition_prob = deposition_prob/E_MAX_Smith;
					}
					
					if (rng.unit() < deposition_prob)
					{
						geometry->deposit(this_particle.pos, last_triangle_normal, 0, particle_mgr.get_primary_tag(particle_idx), this_particle.kin_energy, particle_mgr.get_species(particle_idx), particle_mgr.get_new_species(particle_idx),copy_particle.dir);

						this_particle.kin_energy -= dissociation_energy;

						particle_mgr.update_cascade_diagram(particle_idx); //Cascade Diagram Update

						//Termination Toggle to prevent termination
						if (electron_termination_on_deposition)
						{
							particle_mgr.terminate(particle_idx); // After a deposition, the electron is terminated 
							e_termination = true; 
						}
					}
				}
			}
			//Now we have had deposition check so if particle has not been terminated by deposition we move onto transmission through boundary
			if (this_particle.kin_energy*cos_theta*cos_theta + dU > 0 && e_termination == false)
			{
				const real s = sqrtr(1 + dU / (this_particle.kin_energy*cos_theta*cos_theta));
				const real T = (quantum_transmission
					? 4 * s / ((1 + s)*(1 + s))
					: 1);
				if (rng.unit() < T)
				{
					//Kinetic Energy Adjustment due to transmission
					this_particle.kin_energy += dU; //work function
		
					if (interface_refraction)
					{
						// determine the angle of refraction
						// see thesis T.V. Eq. 3.139
						this_particle.dir = (normalised_dir - last_triangle_normal * cos_theta)
							+ s * last_triangle_normal * cos_theta;
					}
					
					// update the current material index and EXIT.
									
					particle_mgr.set_material_index(particle_idx, material_idx_out);

					if (material_idx_out == nbl::special_materials::VACUUM)
					{
						particle_mgr.set_species(particle_idx, 3); // VE
						
					}
					
					particle_mgr[particle_idx] = this_particle;
					return;
				}
			}
		}
		//No transmission --> move on to surface absorption or internal reflection.

		// surface absorption? (this is in accordance with Kieft & Bosch code)
		// note that the default behaviour has this feature disabled
		if (interface_absorption)
		{
			if (dU < 0 && rng.unit() < expr(1 + 0.5_r*this_particle.kin_energy / dU))
			{
				particle_mgr.terminate(particle_idx);
				return;
			}
		}

		// the only remaining case is internal reflection
		this_particle.dir = normalised_dir - 2 * last_triangle_normal*cos_theta;
		particle_mgr[particle_idx] = this_particle;
	}
};

#endif // __BOUNDARY_INTERSECT_H_
