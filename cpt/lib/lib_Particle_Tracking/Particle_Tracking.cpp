#include "Particle_Tracking.h"
#include "Utils.h"
#include <cmath>
#include <iostream>
#include <tuple>
#include <string>
#include <fstream>
#include <map>

#ifdef _WIN32
#include <ciso646> //for `and` and `or` keyword
#endif

constexpr double PI = 3.141592653589793238463;

Particle_Tracking::Particle_Tracking(
	bool opt_original_ASG,
	bool opt_bin_cell_check_center_only,
	int half_width_search_window,
	int encompassing_cell_spiral_level,
	double nondim_criterion_area_ratio,
	int props_eval_spiral_level,
	bool opt_output,
	bool opt_verbose,
	int num_threads
) :
	opt_original_ASG(opt_original_ASG),
	opt_bin_cell_check_center_only(opt_bin_cell_check_center_only),
	half_width_search_window(half_width_search_window),
	encompassing_cell_spiral_level(encompassing_cell_spiral_level),
	nondim_criterion_area_ratio(nondim_criterion_area_ratio),
	props_eval_spiral_level(props_eval_spiral_level),
	opt_output(opt_output),
	opt_verbose(opt_verbose),
	num_threads(num_threads)
{
	grid = new Grid(
		opt_original_ASG, 
		opt_bin_cell_check_center_only,
		half_width_search_window, 
		encompassing_cell_spiral_level, 
		nondim_criterion_area_ratio,
		opt_output,
		opt_verbose
	);

	printf("Particle_Tracking Constructed!\n");
	num_particles = 0;
}

Particle_Tracking::~Particle_Tracking()
{
	delete grid;
	std::vector<Particle>().swap(particle);
	printf("Particle_Tracking Destructed.\n");
}


void Particle_Tracking::read_particle_file(
	const char* path_particle_file,
	const char* delimiter
	)
{
	using std::cout;
	using std::endl;
	using std::ifstream;
	using std::string;
	using std::size_t;

	std::string delimiter_str(delimiter);
	string line, token;
	size_t pos;
	std::vector<string> stack_token;
	std::map<string, int> map_column;
	std::map<string, int>::iterator it_pos_x, it_pos_y,
		it_diameter, it_density, it_vel_x, it_vel_y,
		it_opt_release_at_flow_speed;
	double* this_pos, * this_vel, 
		this_density, this_diameter;
	this_pos = new double[ndim];
	this_vel = new double[ndim];
	bool this_opt_release_at_flow_speed;
	//------------------------------------------------

	ifstream infile;
	infile.open(path_particle_file);

	if (infile.is_open())
	{
		// verbose
		if (opt_verbose)
			cout << "Reading the particle file starts ..." << endl;

		for (int sn_line = 1; getline(infile, line); sn_line++)
		{
			// Skipping empty lines
			if (line.length() == 0)
			{
				continue;
			}

			stack_token.clear();
			pos = 0;
			while ((pos = line.find(delimiter)) != std::string::npos)
			{
				token = line.substr(0, pos);
				token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
				line.erase(0, pos + delimiter_str.length());
				stack_token.push_back(token);
			}
			// last token
			line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
			if (line.length() > 0)
				stack_token.push_back(line);

			// first line: header
			if (sn_line == 1)
			{
				for (int sn_item = 0; sn_item < (int)stack_token.size(); sn_item++)
					map_column[stack_token.at(sn_item)] = sn_item;

				// iterators to keys -- can be extended
				it_pos_x = map_column.find("pos.x");
				it_pos_y = map_column.find("pos.y");
				it_diameter = map_column.find("diameter");
				it_density = map_column.find("density");
				it_vel_x = map_column.find("vel.x");
				it_vel_y = map_column.find("vel.y");
				it_opt_release_at_flow_speed = map_column.find("opt_release_at_flow_speed");

				// coordinates of particle
				if (it_pos_x == map_column.end() or it_pos_y == map_column.end())
				{
					if (opt_verbose)
						cout << "\nError:\nCoordinates must be included in the particle file.\n\n";
					exit(1);
				}
			}
			// body
			else
			{

				// #pragma region Coordinates
				if (it_pos_x != map_column.end())
					this_pos[0] = std::stod(stack_token.at(it_pos_x->second));

				if (it_pos_y != map_column.end())
					this_pos[1] = std::stod(stack_token.at(it_pos_y->second));

				if (it_diameter != map_column.end())
					this_diameter = std::stod(stack_token.at(it_diameter->second));

				if (it_density != map_column.end())
					this_density = std::stod(stack_token.at(it_density->second));

				if (it_vel_x != map_column.end())
					this_vel[0] = std::stod(stack_token.at(it_vel_x->second));
				
				if (it_vel_y != map_column.end())
					this_vel[1] = std::stod(stack_token.at(it_vel_y->second));

				if (it_opt_release_at_flow_speed != map_column.end())
				{
					if (stack_token.at(it_opt_release_at_flow_speed->second) == "1" or
						stack_token.at(it_opt_release_at_flow_speed->second) == "True" or
						stack_token.at(it_opt_release_at_flow_speed->second) == "true")
						this_opt_release_at_flow_speed = true;
					else
						this_opt_release_at_flow_speed = false;
				}

				if (it_opt_release_at_flow_speed == map_column.end() and
					(it_vel_x == map_column.end() or it_vel_y == map_column.end()))
					set_particles(//set a single stationary particle
						this_pos,
						this_diameter,
						this_density,
						false);
				else if (it_opt_release_at_flow_speed != map_column.end())
					if (this_opt_release_at_flow_speed or 
						it_vel_x == map_column.end() or it_vel_y == map_column.end())
						// The case that `opt_release_at_flow_speed` if true, 
						// or neither of `vel.x` and `vel.y` is provided
						// so even when `vel.x` and `vel.y` are present, `opt_release_at_flow_speed` can
						// supersede them.
						set_particles(//set a single particle with an option to release at flow speed or stationary.
							this_pos,
							this_diameter,
							this_density,
							this_opt_release_at_flow_speed);
					else
						//the case of full header: `opt_release_at_flow_speed`, `vel.x`, and `vel.y` are available
						set_particles(//set a single particle with a given velocity
							this_pos,
							this_diameter,
							this_density,
							this_vel);
				else if (it_vel_x != map_column.end() and it_vel_y != map_column.end())
					set_particles(//set a single particle with a given velocity
						this_pos,
						this_diameter,
						this_density,
						this_vel);
				else {
					if (opt_verbose)
						cout << "\nError:\nCould not configure the velocity for a particle.\n\n";
					exit(1);
				}
			}
		}
	
		// end of parsing
		if (opt_verbose)
			cout << "Reading the particles file ends.\n\n";
	}
	else
	{
		if (opt_verbose)
			cout << "\nCould not open the particles file: " << path_particle_file << endl;
	}
}

void Particle_Tracking::set_grid_and_sol(
	const char *path_mesh,
	const char *path_sol_nodes,
	const char *path_sol_cells,
	const char *path_report_mesh,
	const char *path_report_mesh_ordered_connectivity,
	const char *path_report_sol_nodes,
	const char *path_report_sol_cells,
	const char *path_report_bins,
	bool opt_read_sol_at_cells_center,
	int *n_bins,
	double allowed_tolerance_node_position_snap,
	double nondim_allowed_tolerance_cell_center_position_snap,
	bool strict_mode,
	bool opt_report_mesh,
	bool opt_report_mesh_ordered,
	bool opt_report_sol_nodes,
	bool opt_report_sol_cells,
	bool opt_report_bins,
	bool opt_report_nodes_of_zones_one_file_per_zone,
	bool opt_report_nodes_of_zones_single_file,
	const char* dir_nodes_of_zones,
	const char* fname_base_nodes_of_zones,
	const char* file_extension_nodes_of_zones,
	const char* delimiter_nodes_of_zones,
	int report_every_num_line,
	double bin_len_dimensionless
	)
{
	set_grid(
		path_mesh,
		path_report_mesh,
		path_report_mesh_ordered_connectivity,
		path_report_bins,
		n_bins,
		opt_report_mesh,
		opt_report_mesh_ordered,
		opt_report_bins,
		opt_report_nodes_of_zones_one_file_per_zone,
		opt_report_nodes_of_zones_single_file,
		dir_nodes_of_zones,
		fname_base_nodes_of_zones,
		file_extension_nodes_of_zones,
		delimiter_nodes_of_zones,
		report_every_num_line,
		bin_len_dimensionless
	);

	set_sol(
		path_sol_nodes,
		path_sol_cells,
		path_report_sol_nodes,
		path_report_sol_cells,
		opt_read_sol_at_cells_center,
		allowed_tolerance_node_position_snap,
		nondim_allowed_tolerance_cell_center_position_snap,
		strict_mode,
		opt_report_sol_nodes,
		opt_report_sol_cells,
		report_every_num_line);
}

void Particle_Tracking::set_grid(
	const char *path_mesh,
	const char *path_report_mesh,
	const char *path_report_mesh_ordered_connectivity,
	const char *path_report_bins,
	int *n_bins,
	bool opt_report_mesh,
	bool opt_report_mesh_ordered,
	bool opt_report_bins,
	bool opt_report_nodes_of_zones_one_file_per_zone,
	bool opt_report_nodes_of_zones_single_file,
	const char* dir_nodes_of_zones,
	const char* fname_base_nodes_of_zones,
	const char* file_extension_nodes_of_zones,
	const char* delimiter_nodes_of_zones,
	int report_every_num_line,
	double bin_len_dimensionless
	)
{
	// read mesh and store positions and connectivity info
	grid->read_ansys_fluent_mesh_line_parsing(path_mesh, report_every_num_line);
	// report of extracted mesh
	if (opt_report_mesh)
		grid->report_mesh(path_report_mesh);
	// order faces and nodes on each cell based on right-hand rule
	grid->order_cell_connectivity();
	if (opt_report_mesh_ordered)
		grid->report_mesh(path_report_mesh_ordered_connectivity);

	// storing cell area as needed in particle tracking to determine relaxation time of fluid flow in a cell
	grid->calc_cells_area();

	//stats of grid -- needed for automatic formation of bins based on the length of grid
	grid->eval_stats_faces();

	//calc pos of cell center from nodes -- needed for bin.cell population, which itself is needed for original ASG variant
	//passing no arg, by default leads to evaluation of position of cell. The same method with appropriate args will later
	//be used in `set_sol` method to evaluate properties of cell from nodes.
	grid->calc_cell_props_from_nodes();

	// set up bins
	grid->set_up_bins(n_bins, path_report_bins, opt_report_bins, bin_len_dimensionless);

	// set ndim
	ndim = grid->ndim;

	//report nodes of type zones
	if (opt_report_nodes_of_zones_one_file_per_zone)
		grid->report_nodes_of_zones_one_file_per_zone(
			dir_nodes_of_zones,
			fname_base_nodes_of_zones,
			file_extension_nodes_of_zones,
			delimiter_nodes_of_zones);

	if (opt_report_nodes_of_zones_single_file)
		grid->report_nodes_of_zones_single_file(
			dir_nodes_of_zones,
			fname_base_nodes_of_zones,
			file_extension_nodes_of_zones,
			delimiter_nodes_of_zones);
}

void Particle_Tracking::set_sol(
	const char *path_sol_nodes,
	const char *path_sol_cells,
	const char *path_report_sol_nodes,
	const char *path_report_sol_cells,
	bool opt_read_sol_at_cells_center,
	double allowed_tolerance_node_position_snap,
	double nondim_allowed_tolerance_cell_center_position_snap,
	bool strict_mode,
	bool opt_report_sol_nodes,
	bool opt_report_sol_cells,
	int report_every_num_line)
{
	// Note: sol at nodes MUST be provided.
	grid->read_ansys_fluent_sol_at_nodes(
		path_sol_nodes,
		allowed_tolerance_node_position_snap,
		report_every_num_line);
	
	// sol at cells center
	// Note: sol at nodes MUST be provided. 
	// sol at cells can also be provided (optional).
	// in the case sol at cells is not provided, it 
	// will be interpolated from the nodes. 
	// This option was written solely to check the interpolations accuracy
	// and whether they match the values provided by Ansys solution file at cell centers.
	// In summary, `opt_read_sol_at_cells_center` can be always set to false safely.
	if (opt_read_sol_at_cells_center)
	{
		// Important to find cell center position first to determine the correct cell from coordinates provided in sol_cells file
		// The following is commented as it is already included in `set_grid` method
		//grid->calc_cell_props_from_nodes(true, false, false, false, false);

		grid->read_ansys_fluent_sol_at_cells(
			path_sol_cells,
			nondim_allowed_tolerance_cell_center_position_snap,
			report_every_num_line,
			strict_mode);

		// aux props of cells, e.g. vel_mag
		grid->calc_cells_aux_props();
	}
	else
	{
		grid->calc_cell_props_from_nodes(true, true, true, true, true);

		// aux props of cells, e.g. vel_mag
		grid->calc_cells_aux_props();
	}

	// reports
	if (opt_report_sol_nodes)
		grid->report_sol_nodes(path_report_sol_nodes);
	if (opt_report_sol_cells)
		grid->report_sol_cells(path_report_sol_cells);
}

void Particle_Tracking::set_particle_source(
	int n_particles,
	double *source_pos,
	double range_diameter[2],
	double range_density[2],
	bool opt_release_at_flow_speed)
{
	double **range_vel;
	range_vel = new double *[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		range_vel[ind_dir] = new double[2];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		for (int ind_level = 0; ind_level < 2; ind_level++)
			range_vel[ind_dir][ind_level] = 0.;

	set_particle_source(
		n_particles,
		source_pos,
		range_vel,
		range_diameter,
		range_density,
		opt_release_at_flow_speed);
}

void Particle_Tracking::set_particle_source(
	int n_particles,
	double* source_pos,
	double** range_vel,
	double range_diameter[2], //[min, max]
	double range_density[2],  //[min, max]
	bool opt_release_at_flow_speed
)
{
	double* dia, * rho;
	dia = new double[n_particles];
	rho = new double[n_particles];
	double** pos, ** vel, ** acc, ** force;
	pos = new double* [n_particles];
	vel = new double* [n_particles];
	acc = new double* [n_particles];
	force = new double* [n_particles];
	for (int ind_particle = 0; ind_particle < n_particles; ind_particle++)
	{
		pos[ind_particle] = new double[ndim];
		vel[ind_particle] = new double[ndim];
		acc[ind_particle] = new double[ndim];
		force[ind_particle] = new double[ndim];
		//------------------------------------------
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			pos[ind_particle][ind_dir] = source_pos[ind_dir];
			acc[ind_particle][ind_dir] = 0;
			force[ind_particle][ind_dir] = 0;
		}

		if (n_particles == 1)
		{
			dia[ind_particle] = range_diameter[0];
			rho[ind_particle] = range_density[0];
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				vel[ind_particle][ind_dir] = range_vel[ind_dir][0];
		}
		else
		{
			// linearly varying props
			dia[ind_particle] = range_diameter[0] + ind_particle * (range_diameter[1] - range_diameter[0]) / (n_particles - 1);
			rho[ind_particle] = range_density[0] + ind_particle * (range_density[1] - range_density[0]) / (n_particles - 1);
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				vel[ind_particle][ind_dir] = range_vel[ind_dir][0] +
				ind_particle * (range_vel[ind_dir][1] - range_vel[ind_dir][0]) / (n_particles - 1);
		}
	}

	set_particles(n_particles, dia, rho, pos, vel, acc, force, opt_release_at_flow_speed);
}


void Particle_Tracking::set_streakline(
	int n_particles,
	double* pos_end_1,
	double* pos_end_2,
	double * vel_end_1,
	double* vel_end_2,
	double range_diameter[2], //[min, max]
	double range_density[2],  //[min, max]
	bool opt_release_at_flow_speed
)
{
	double *dia, *rho;
	dia = new double[n_particles];
	rho = new double[n_particles];
	double **pos, **vel, **acc, **force;
	pos = new double *[n_particles];
	vel = new double *[n_particles];
	acc = new double *[n_particles];
	force = new double *[n_particles];
	for (int ind_particle = 0; ind_particle < n_particles; ind_particle++)
	{
		pos[ind_particle] = new double[ndim];
		vel[ind_particle] = new double[ndim];
		acc[ind_particle] = new double[ndim];
		force[ind_particle] = new double[ndim];
		//------------------------------------------
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			acc[ind_particle][ind_dir] = 0;
			force[ind_particle][ind_dir] = 0;
		}
		
		if (n_particles == 1)
		{
			dia[ind_particle] = range_diameter[0];
			rho[ind_particle] = range_density[0];
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			{
				pos[ind_particle][ind_dir] = pos_end_1[ind_dir];
				vel[ind_particle][ind_dir] = vel_end_1[ind_dir];
			}
		}
		else
		{
			// linearly varying props
			dia[ind_particle] = range_diameter[0] + ind_particle * (range_diameter[1] - range_diameter[0]) / (n_particles - 1);
			rho[ind_particle] = range_density[0] + ind_particle * (range_density[1] - range_density[0]) / (n_particles - 1);
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			{
				pos[ind_particle][ind_dir] = pos_end_1[ind_dir] +
					ind_particle * (pos_end_2[ind_dir] - pos_end_1[ind_dir]) / (n_particles - 1);
				vel[ind_particle][ind_dir] = vel_end_1[ind_dir] +
					ind_particle * (vel_end_2[ind_dir] - vel_end_1[ind_dir]) / (n_particles - 1);
			}
		}
	}

	set_particles(n_particles, dia, rho, pos, vel, acc, force, opt_release_at_flow_speed);
}

void Particle_Tracking::set_streakline(
	int n_particles,
	double* pos_end_1,
	double* pos_end_2,
	double range_diameter[2], //[min, max]
	double range_density[2],  //[min, max]
	bool opt_release_at_flow_speed
)
{
	double* vel_end_1, * vel_end_2;
	vel_end_1 = new double [ndim];
	vel_end_2 = new double[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		vel_end_1[ind_dir] = 0.;
		vel_end_2[ind_dir] = 0.;
	}

	set_streakline(
		n_particles,
		pos_end_1,
		pos_end_2,
		vel_end_1,
		vel_end_2,
		range_diameter,
		range_density,
		opt_release_at_flow_speed
	);
}

void Particle_Tracking::set_particles(
	int n_particles,
	double* dia,
	double* rho,
	double** pos,
	double** vel,
	double** acc,
	double** force,
	bool opt_release_at_flow_speed
	)
{
	/*
	 * Most generic function to set particle(s)
	 * `set_particles` can be called multiple times such that
	 * one or multiple particles can be added to the system in
	 * each call.
	 */

	 // declarations
	 //--------------------------------
	int sn_cell, this_id_particle = 0;
	bool flag;
	Local_Flow_Props local_props(ndim);
	//--------------------------------

	if (opt_verbose)
		std::cout << "Setting particles...\n";

	num_particles += n_particles;
	for (int ind_particle = 0; ind_particle < n_particles; ind_particle++)
	{
		this_id_particle = (int)particle.size(); // 0-based id
		particle.push_back(Particle());

		particle[this_id_particle].set_id(this_id_particle);
		particle[this_id_particle].set_diameter(dia[ind_particle]);
		particle[this_id_particle].set_density(rho[ind_particle]);

		// set position
		particle[this_id_particle].set_pos(pos[ind_particle]);

		// get flow velocity if particle is to be released at flow speed
		if (opt_release_at_flow_speed)
		{
			std::tie(sn_cell, flag) = eval_local_props(
				pos[ind_particle],
				local_props,
				props_eval_spiral_level
				);
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				vel[ind_particle][ind_dir] = local_props.vel[ind_dir];
		}

		// set vel
		particle[this_id_particle].set_vel(vel[ind_particle]);

		// Set acc and force if needed for an application.
		// Currently, they are not being used within the time integration.
		particle[this_id_particle].set_acc(acc[ind_particle]);
		particle[this_id_particle].set_force(force[ind_particle]);

		if (opt_verbose)
			std::cout << "Setting particle " << this_id_particle << " is done.\n";
	}
	if (opt_verbose)
		std::cout << "Setting particles is done.\n\n";
}

void Particle_Tracking::set_particles(double *single_pos, bool opt_release_at_flow_speed)
{
	int n_particles = 1;
	double dia[1] = {0};
	double rho[1] = {0};
	double **pos, **vel, **acc, **force;
	pos = new double *[1];
	vel = new double *[1];
	acc = new double *[1];
	force = new double *[1];
	pos[0] = new double[ndim];
	vel[0] = new double[ndim];
	acc[0] = new double[ndim];
	force[0] = new double[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		pos[0][ind_dir] = single_pos[ind_dir];
		vel[0][ind_dir] = 0;
		acc[0][ind_dir] = 0;
		force[0][ind_dir] = 0;
	}

	set_particles(n_particles, dia, rho, pos, vel, acc, force, opt_release_at_flow_speed);
}

void Particle_Tracking::set_particles(
	double *single_pos,
	double diameter,
	double density,
	bool opt_release_at_flow_speed)
{
	/*
	* set a single particle with an option to release at flow speed or stationary.
	*/
	int n_particles = 1;

	double *dia, *rho;
	dia = new double[n_particles];
	rho = new double[n_particles];

	dia[0] = diameter;
	rho[0] = density;

	double **pos, **vel, **acc, **force;
	pos = new double *[n_particles];
	vel = new double *[n_particles];
	acc = new double *[n_particles];
	force = new double *[n_particles];
	pos[0] = new double[ndim];
	vel[0] = new double[ndim];
	acc[0] = new double[ndim];
	force[0] = new double[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		pos[0][ind_dir] = single_pos[ind_dir];
		vel[0][ind_dir] = 0;
		acc[0][ind_dir] = 0;
		force[0][ind_dir] = 0;
	}

	set_particles(n_particles, dia, rho, pos, vel, acc, force, opt_release_at_flow_speed);
}

void Particle_Tracking::set_particles(
	double* single_pos,
	double diameter,
	double density,
	double* single_vel)
{
	/*
	* set a single particle with a given velocity
	*/
	int n_particles = 1;

	double* dia, * rho;
	dia = new double[n_particles];
	rho = new double[n_particles];

	dia[0] = diameter;
	rho[0] = density;

	double** pos, ** vel, ** acc, ** force;
	pos = new double* [n_particles];
	vel = new double* [n_particles];
	acc = new double* [n_particles];
	force = new double* [n_particles];
	pos[0] = new double[ndim];
	vel[0] = new double[ndim];
	acc[0] = new double[ndim];
	force[0] = new double[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		pos[0][ind_dir] = single_pos[ind_dir];
		vel[0][ind_dir] = single_vel[ind_dir];
		acc[0][ind_dir] = 0;
		force[0][ind_dir] = 0;
	}

	set_particles(n_particles, dia, rho, pos, vel, acc, force);
}

void Particle_Tracking::track_particle(
	const char *dir_base_particle_hist,
	double time_duration,
	double time_step_over_taw_relax,
	// std::vector<int> wall_zone_type,
	int* wall_zone_type,
	// std::vector<int> particle_id_vec,
	int *particle_id_ptr,
	double L_ref,
	double Vel_ref,
	int num_wall_zone_type,
	int n_particles_to_track,
	int max_iter,
	int verbose_every_num_step,
	double elastic_reflection_coeff,
	bool opt_time_integration_for_small_stokes,
	double stokes_number_threshold,
	double fixed_time_step,
	int num_points_on_perimeter,
	bool opt_one_step_euler,
	bool opt_rk4,
	const char *fname_base,
	const char *file_extension,
	const char *delimiter)
{
	
	//---------------------------------------------
	// declaretions
	//---------------------------------------------
	std::vector<int> particle_id_vec;
	//---------------------------------------------

	// config the default case
	if (n_particles_to_track <= 0)
	{
		particle_id_vec.clear();
		for (int ind_particle = 0; ind_particle < num_particles; ind_particle++)
		{
			particle_id_vec.push_back(particle[ind_particle].id);
		}
	}
	else
	{
		particle_id_vec.clear();
		for (int ind_particle = 0; ind_particle < n_particles_to_track; ind_particle++)
		{
			particle_id_vec.push_back(particle_id_ptr[ind_particle]);
		}
	}
	
	// dev. performance measurement
	Event_Execution_Time eet(
		0,
		"default_tag",
		"./event_execution_time_performance",
		true
	);

	eet.start_event("outermost loop");

	if (opt_verbose)
		std::cout << "Start of particle tracking." << std::endl;

	
//#pragma omp parallel for
#pragma omp parallel for num_threads(num_threads)
	for (int ind_particle = 0; ind_particle < (int)particle_id_vec.size(); ind_particle++)
	{
		if (omp_get_num_threads() > 1)
			eet.active = false;

		track_single_particle(
			particle[particle_id_vec[ind_particle]],
			dir_base_particle_hist,
			time_duration,
			time_step_over_taw_relax,
			wall_zone_type,
			L_ref,
			Vel_ref,
			num_wall_zone_type,
			max_iter,
			verbose_every_num_step,
			elastic_reflection_coeff,
			opt_time_integration_for_small_stokes,
			stokes_number_threshold,
			fixed_time_step,
			num_points_on_perimeter,
			opt_one_step_euler,
			opt_rk4,
			fname_base,
			file_extension,
			delimiter,
			eet
		);
	} // end of loop on particles

	// dev. performance measurement
	eet.active = true;
	eet.stop_event("outermost loop");
	eet.write();
	//------------------------------------------------
}

void Particle_Tracking::track_single_particle(
	Particle& this_particle,
	const char* dir_base_particle_hist,
	double time_duration,
	double time_step_over_taw_relax,
	int* wall_zone_type,
	double L_ref,
	double Vel_ref,
	int num_wall_zone_type,
	int max_iter,
	int verbose_every_num_step,
	double elastic_reflection_coeff,
	bool opt_time_integration_for_small_stokes,
	double stokes_number_threshold,
	double fixed_time_step,
	int num_points_on_perimeter,
	bool opt_one_step_euler,
	bool opt_rk4,
	const char* fname_base,
	const char* file_extension,
	const char* delimiter,
	const Event_Execution_Time& eet
	)
{
	//-------------------------------------------------------------------
	// config of sanity check at the beginning of run
	//-------------------------------------------------------------------
	bool opt_check_circumferential_points_are_inside_domain = true;
	double theta_p, * pos_circumferential_point;
	pos_circumferential_point = new double[ndim];
	int sn_cell_circumferential_point, flag_circumferential_point;
	std::vector<int> cell_search_stack_circumferential_point;
	std::vector<double> dist_from_nodes_circumferential_point;
	int flag_finding_encompassing_cell_normal_return_circumferential_point = 0;
	int flag_finding_encompassing_cell_not_finding_encomp_cell_circumferential_point = -1;
	int flag_finding_encompassing_cell_not_finding_close_node_circumferential_point = -2;
	//-------------------------------------------------------------------

	// To concatenate paths later,
	// trimming any chars of '\' (for windows),
	//' ', i.e spaces,
	// and '/' (for linux) from the end of path
	std::string dir_base_particle_hist_str(dir_base_particle_hist);
	const char* endl_chars_to_be_trimmed = "\\ /";
	trim_string_endl(dir_base_particle_hist_str, endl_chars_to_be_trimmed);
#if defined(_WIN32)
	dir_base_particle_hist_str += "\\";
#else
	dir_base_particle_hist_str += "/";
#endif

	//---------------------------------------------
	// declarations
	//---------------------------------------------
	Local_Flow_Props local_props(ndim);

	double* pos_p, * pos_p_new,
		* vel_p, * vel_p_mid, * vel_p_new,
		* acc_p, * acc_p_new,
		* pos_after_bc, * vel_after_bc,
		taw_relax_p, taw_relax_f;
	pos_p = new double[ndim];
	pos_p_new = new double[ndim];
	vel_p = new double[ndim];
	vel_p_mid = new double[ndim];
	vel_p_new = new double[ndim];
	acc_p = new double[ndim];
	acc_p_new = new double[ndim];
	pos_after_bc = new double[ndim];
	vel_after_bc = new double[ndim];

	double** rk_pos, ** rk_vel;
	rk_pos = new double* [4];
	rk_vel = new double* [4];
	for (int rk_step = 0; rk_step < 4; rk_step++)
	{
		rk_pos[rk_step] = new double[ndim];
		rk_vel[rk_step] = new double[ndim];
	}

	bool stokes_number_is_small = false;
	double time_step = 0, elapsed_time = 0, start_time = 0, EPSILON = 1e-12, stokes_number;
	int sn_iter, sn_cell, num_line_print = 0;
	int flag, flag_normal_return = 0, flag_out_of_domain = -1, flag_return_bc;
	std::string delimiter_verbose = "\t";

	std::vector<int> wall_zone_type_vec;
	for (int ind_wall_zone_type = 0; ind_wall_zone_type < num_wall_zone_type; ind_wall_zone_type++)
		wall_zone_type_vec.push_back(wall_zone_type[ind_wall_zone_type]);
	//-------------------------------------------------------------------------------------------

	// start time for the current particle
	if (this_particle.hist.time.size() > 0)
		start_time = this_particle.hist.time.back();
	else
	{
		start_time = 0;
		// add the initial state to history
		this_particle.add_to_hist(0.);
	}

	elapsed_time = start_time; // set to the elapsed time of the current particle
	sn_iter = 0;			   // reset to 0 for each round of tracking of particle

	// #pragma region tracking each particle
	//while (elapsed_time < time_duration and sn_iter < max_iter)
	while (elapsed_time < time_duration)
	{
		if (sn_iter >= max_iter)
		{
			if (opt_verbose)
				std::cout << "\nMax number of iterations, i.e. " << max_iter << ", is reached.\n\n";
			break;
		}

		// using the current particle pos to find the new pos
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			// `pos_p` and `vel_p` are pos and vel at time `n`. 
			// `pos_p_new` and `vel_p_new` are pos and vel at time `n+1`. 
			pos_p[ind_dir] = this_particle.props.pos[ind_dir];
			vel_p[ind_dir] = this_particle.props.vel[ind_dir];
		}

		//locating the hosting cell and evaluation of flow props only for the first iteration
		//this process is performed at the last step of time integration so the status of 
		//fluid hosting cell and flow props will be updated at the beginning of next time 
		//step (here).
		if (sn_iter == 0)
		{
			// event: start
			const_cast<Event_Execution_Time&>(eet).start_event(std::string("first iter: eval_local_props"));

			// #pragma region interpolating fluid flow properties
			std::tie(sn_cell, flag) = eval_local_props(
				pos_p,
				local_props,
				props_eval_spiral_level,
				flag_normal_return,
				flag_out_of_domain);

			const_cast<Event_Execution_Time&>(eet).stop_event("first iter: eval_local_props");
			// event: stop

			if (flag == flag_out_of_domain)
			{
				if (opt_verbose)
					std::cout << "Main time iteration loop | Particle out of domain.\n\n";
				break;
			}

			// event: start				
			const_cast<Event_Execution_Time&>(eet).start_event("first iter: check_circumferential_points");

			//check if all circumferential points are inside the domain
			if (opt_check_circumferential_points_are_inside_domain)
			{
				for (int ind_p = 0; ind_p < num_points_on_perimeter; ind_p++)
				{
					// #pragma region finding current and old positions of this circumferential point
					theta_p = ind_p * 2. * PI / num_points_on_perimeter;
					for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
						pos_circumferential_point[ind_dir] = pos_p[ind_dir];
					pos_circumferential_point[0] += this_particle.props.diameter / 2. * cos(theta_p);
					pos_circumferential_point[1] += this_particle.props.diameter / 2. * sin(theta_p);

					// #pragma region Attempting to finding the encompassing cell
					std::tie(sn_cell_circumferential_point,
						dist_from_nodes_circumferential_point,
						flag_circumferential_point,
						cell_search_stack_circumferential_point) = grid->find_cell_encompassing_pos(
							pos_circumferential_point,
							flag_finding_encompassing_cell_normal_return_circumferential_point,
							flag_finding_encompassing_cell_not_finding_encomp_cell_circumferential_point,
							flag_finding_encompassing_cell_not_finding_close_node_circumferential_point,
							eet
						);
					// #pragma endregion

					// #pragma region Not finding close node
					if (flag_circumferential_point == flag_finding_encompassing_cell_not_finding_close_node_circumferential_point)
					{
						if (opt_verbose)
						{
							std::cout << "\nError\nNot finding a close node for a circumferential point on \n"
								<< "\nparticle position: ("
								<< pos_p[0] << ", "
								<< pos_p[1] << ")"
								<< "\ncircumferential point position: ("
								<< pos_circumferential_point[0] << ", "
								<< pos_circumferential_point[1] << ")\n";
							std::cout << "A potential remedy could be to increase the `half_width_search_window` which currently is\n";
							std::cout << half_width_search_window << "\n";
						}
						break;
					}
					else if (flag_circumferential_point == flag_finding_encompassing_cell_not_finding_encomp_cell_circumferential_point)
					{
						if (opt_verbose)
							std::cout << "\nError\nUnacceptable release location for Particle " << this_particle.id
							<< "\nparticle position: ("
							<< pos_p[0] << ", "
							<< pos_p[1] << ")"
							<< "\nA hosting cell could not be found for the circumferential point: ("
							<< pos_circumferential_point[0] << ", "
							<< pos_circumferential_point[1] << ")\n\n";
						break;
					}
				}
				if (flag_circumferential_point == flag_finding_encompassing_cell_not_finding_close_node_circumferential_point or
					flag_circumferential_point == flag_finding_encompassing_cell_not_finding_encomp_cell_circumferential_point)
					break;
			}
			const_cast<Event_Execution_Time&>(eet).stop_event("first iter: check_circumferential_points");
			// event: stop
		}//end of locating the hosting cell and evaluation of flow props only for the first iteration

		// #pragma endregion

		// #pragma region integrator
		if (local_props.viscosity <= 0)
		{
			if (opt_verbose)
			{
				std::cout << "\n\nError:\n";
				std::cout << "Invalide viscosity: " << local_props.viscosity << "\n";
				std::cout
					<< "sn_iter: " << sn_iter << " | "
					<< "x, y: "
					<< pos_p[0] << ", "
					<< pos_p[1] << " | "
					<< "fluid v.x, v.y: " << local_props.vel[0] << ", " << local_props.vel[1]
					<< "\n\n";
			}
			exit(1);
		}

		//-------------------------------------------------------------
		// Stokes-based time integration scheme
		//-------------------------------------------------------------

		//local relaxation time of particle
		taw_relax_p = this_particle.props.density
			* std::pow(this_particle.props.diameter, 2.) / 18.0 / local_props.viscosity;

		if (L_ref > 0 and Vel_ref > 0)//for the case that L_ref and Vel_ref are both provided
			taw_relax_f = L_ref / Vel_ref;
		else if (L_ref > 0)//for the case that only L_ref is provided; Vel_ref is the local cell-based velocity magnitude
			taw_relax_f = L_ref / std::sqrt(std::pow(local_props.vel[0], 2.) + std::pow(local_props.vel[1], 2.));
		else
			taw_relax_f = fixed_time_step / time_step_over_taw_relax;

		stokes_number = taw_relax_p / (taw_relax_f + EPSILON);
		if (stokes_number < stokes_number_threshold or opt_time_integration_for_small_stokes)
			stokes_number_is_small = true;
		else
			stokes_number_is_small = false;

		time_step = time_step_over_taw_relax * taw_relax_f;
		//Time integration
		if (stokes_number_is_small)
		{
			if (opt_one_step_euler)
			{
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					pos_p_new[ind_dir] = pos_p[ind_dir] + time_step * vel_p[ind_dir];
				this_particle.set_pos(pos_p_new);
			}
			else if (opt_rk4)
			{
				//Step 1
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				{
					rk_pos[0][ind_dir] = pos_p[ind_dir];
					rk_vel[0][ind_dir] = local_props.vel[ind_dir];
				}

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 1 -- apply_bc");

				//bc
				std::tie(pos_after_bc, vel_after_bc, flag_return_bc) = apply_bc(
					this_particle.hist.pos.back(),
					rk_pos[0],
					rk_vel[0],
					this_particle.props.diameter,
					wall_zone_type_vec,
					elastic_reflection_coeff,
					num_points_on_perimeter,
					eet
				);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 1 -- apply_bc");
				// event: stop

				if (flag_return_bc < 0)
				{
					if (opt_verbose)
						std::cout << "A circumferential point seems to be out of domain (?)\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				{
					rk_pos[0][ind_dir] = pos_after_bc[ind_dir];
					rk_vel[0][ind_dir] = vel_after_bc[ind_dir];
				}

				//step 2
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					rk_pos[1][ind_dir] = pos_p[ind_dir] + time_step / 2 * rk_vel[0][ind_dir];

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 2 -- eval_local_props");

				std::tie(sn_cell, flag) = eval_local_props(
					rk_pos[1],
					local_props,
					props_eval_spiral_level,
					flag_normal_return,
					flag_out_of_domain);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 2 -- eval_local_props");
				// event: stop


				if (flag == flag_out_of_domain)
				{
					if (opt_verbose)
						std::cout << "Main time iteration loop | Particle out of domain.\n\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					rk_vel[1][ind_dir] = local_props.vel[ind_dir];

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 2 -- apply_bc");

				//bc
				std::tie(pos_after_bc, vel_after_bc, flag_return_bc) = apply_bc(
					this_particle.hist.pos.back(),
					rk_pos[1],
					rk_vel[1],
					this_particle.props.diameter,
					wall_zone_type_vec,
					elastic_reflection_coeff,
					num_points_on_perimeter,
					eet
				);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 2 -- apply_bc");
				// event: stop

				if (flag_return_bc < 0)
				{
					if (opt_verbose)
						std::cout << "A circumferential point seems to be out of domain (?)\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				{
					rk_pos[1][ind_dir] = pos_after_bc[ind_dir];
					rk_vel[1][ind_dir] = vel_after_bc[ind_dir];
				}

				//step 3
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					rk_pos[2][ind_dir] = pos_p[ind_dir] + time_step / 2 * rk_vel[1][ind_dir];

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 3 -- eval_local_props");

				std::tie(sn_cell, flag) = eval_local_props(
					rk_pos[2],
					local_props,
					props_eval_spiral_level,
					flag_normal_return,
					flag_out_of_domain);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 3 -- eval_local_props");
				// event: stop

				if (flag == flag_out_of_domain)
				{
					if (opt_verbose)
						std::cout << "Main time iteration loop | Particle out of domain.\n\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					rk_vel[2][ind_dir] = local_props.vel[ind_dir];

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 3 -- apply_bc");

				//bc
				std::tie(pos_after_bc, vel_after_bc, flag_return_bc) = apply_bc(
					this_particle.hist.pos.back(),
					rk_pos[2],
					rk_vel[2],
					this_particle.props.diameter,
					wall_zone_type_vec,
					elastic_reflection_coeff,
					num_points_on_perimeter,
					eet
				);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 3 -- apply_bc");
				// event: stop

				if (flag_return_bc < 0)
				{
					if (opt_verbose)
						std::cout << "A circumferential point seems to be out of domain (?)\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				{
					rk_pos[2][ind_dir] = pos_after_bc[ind_dir];
					rk_vel[2][ind_dir] = vel_after_bc[ind_dir];
				}

				//step 4
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					rk_pos[3][ind_dir] = pos_p[ind_dir] + time_step * rk_vel[2][ind_dir];

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 4 -- eval_local_props");

				std::tie(sn_cell, flag) = eval_local_props(
					rk_pos[3],
					local_props,
					props_eval_spiral_level,
					flag_normal_return,
					flag_out_of_domain);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 4 -- eval_local_props");
				// event: stop

				if (flag == flag_out_of_domain)
				{
					if (opt_verbose)
						std::cout << "Main time iteration loop | Particle out of domain.\n\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					rk_vel[3][ind_dir] = local_props.vel[ind_dir];

				// event: start
				const_cast<Event_Execution_Time&>(eet).start_event("rk4: step 4 -- apply_bc");

				//bc
				std::tie(pos_after_bc, vel_after_bc, flag_return_bc) = apply_bc(
					this_particle.hist.pos.back(),
					rk_pos[3],
					rk_vel[3],
					this_particle.props.diameter,
					wall_zone_type_vec,
					elastic_reflection_coeff,
					num_points_on_perimeter,
					eet
				);

				const_cast<Event_Execution_Time&>(eet).stop_event("rk4: step 4 -- apply_bc");
				// event: stop

				if (flag_return_bc < 0)
				{
					if (opt_verbose)
						std::cout << "A circumferential point seems to be out of domain (?)\n";
					break;
				}

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				{
					rk_pos[3][ind_dir] = pos_after_bc[ind_dir];
					rk_vel[3][ind_dir] = vel_after_bc[ind_dir];
				}

				//-------------------------
				//new pos
				//-------------------------
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					pos_p_new[ind_dir] = pos_p[ind_dir] + time_step / 6 *
					(rk_vel[0][ind_dir] + 2 * rk_vel[1][ind_dir] + 2 * rk_vel[2][ind_dir] + rk_vel[3][ind_dir]);

				this_particle.set_pos(pos_p_new);
			}

		}//end of first part of small stokes mode integration
		else
		{
			acc_p = calc_forces(
				local_props,
				vel_p,
				this_particle.props.diameter,
				this_particle.props.density
			);

			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			{
				vel_p_mid[ind_dir] = vel_p[ind_dir] + acc_p[ind_dir] * time_step / 2.;
				pos_p_new[ind_dir] = pos_p[ind_dir] + vel_p_mid[ind_dir] * time_step;
			}

			this_particle.set_pos(pos_p_new);
			this_particle.set_vel(vel_p_mid);
		}//end of first part of velociy verlet for large stokes mode

		// event: start
		const_cast<Event_Execution_Time&>(eet).start_event("after time integratin part 1 -- apply_bc");

		// apply boundary condition
		std::tie(pos_after_bc, vel_after_bc, flag_return_bc) = apply_bc(
			this_particle.hist.pos.back(),
			this_particle.props.pos,
			this_particle.props.vel,
			this_particle.props.diameter,
			wall_zone_type_vec,
			elastic_reflection_coeff,
			num_points_on_perimeter,
			eet
		);

		const_cast<Event_Execution_Time&>(eet).stop_event("after time integratin part 1 -- apply_bc");
		// event: stop

		if (flag_return_bc < 0)
		{
			if (opt_verbose)
				std::cout << "A circumferential point seems to be out of domain (?)\n";
			break;
		}
		this_particle.set_pos(pos_after_bc);
		this_particle.set_vel(vel_after_bc);

		// event: start
		const_cast<Event_Execution_Time&>(eet).start_event("after time integration part 1 -- eval_local_props");

		//find local fluid flow props: velocity, viscosity
		std::tie(sn_cell, flag) = eval_local_props(
			this_particle.props.pos,
			local_props,
			props_eval_spiral_level,
			flag_normal_return,
			flag_out_of_domain);

		const_cast<Event_Execution_Time&>(eet).stop_event("after time integration part 1 -- eval_local_props");
		// event: stop

		if (flag == flag_out_of_domain)
		{
			if (opt_verbose)
			{
				std::cout << "Main time iteration loop\n"
					<< "between Step 1, determination of pos ^ (n + 1), and Step 2, vel ^ (n + 1)\n"
					<< "after `eval_local_props` | Particle out of domain.\n\n";
				if (this_particle.hist.pos.size() > 1)
					std::cout << "position before boundary treatment: (" << this_particle.hist.pos.end()[-2][0] << ", " << this_particle.hist.pos.end()[-2][1] << ")\n";
				std::cout << "position after  boundary treatment: (" << this_particle.hist.pos.end()[-1][0] << ", " << this_particle.hist.pos.end()[-1][1] << ")\n\n";
			}
			break;
		}

		if (stokes_number_is_small)
		{
			//update particle vel based on new position
			this_particle.set_vel(local_props.vel);
		}
		else
		{
			// event: start
			const_cast<Event_Execution_Time&>(eet).start_event("time integration part 2 -- large Stokes -- calc_forces");

			//second step of velocity verlet
			// Update acc based on midstep particle vel and fluid vel.
			// As of now, forces are position-independent.
			acc_p_new = calc_forces(
				local_props,
				this_particle.props.vel, //vel_p_mid,
				this_particle.props.diameter,
				this_particle.props.density
			);

			const_cast<Event_Execution_Time&>(eet).stop_event("time integration part 2 -- large Stokes -- calc_forces");
			// event: stop

			//new vel
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				//vel_p_new[ind_dir] = vel_p_mid[ind_dir] + acc_p_new[ind_dir] * time_step / 2.;
				vel_p_new[ind_dir] = this_particle.props.vel[ind_dir] + acc_p_new[ind_dir] * time_step / 2.;

			//update particle vel
			this_particle.set_vel(vel_p_new);
		}

		// #pragma endregion

		// update time
		sn_iter++;
		elapsed_time += time_step;

		// add the current state to history
		this_particle.add_to_hist(elapsed_time);

		// verbose
		if (opt_verbose)
			if (sn_iter == 1 or sn_iter % verbose_every_num_step == 0)
			{
				if (sn_iter == 1 or num_line_print % 50 == 0)
					std::cout << "\nparticle id" << delimiter_verbose
					<< "iter" << delimiter_verbose
					<< "time_step" << delimiter_verbose
					<< "time" << delimiter_verbose
					<< "x" << delimiter_verbose
					<< "y" << delimiter_verbose
					<< "vel_x" << delimiter_verbose
					<< "vel_y" << delimiter_verbose
					<< "Stokes_number" << delimiter_verbose
					<< std::endl;

				std::cout << this_particle.id << delimiter_verbose
					<< sn_iter << delimiter_verbose
					<< time_step << delimiter_verbose
					<< elapsed_time << delimiter_verbose
					<< this_particle.props.pos[0] << delimiter_verbose
					<< this_particle.props.pos[1] << delimiter_verbose
					<< this_particle.props.vel[0] << delimiter_verbose
					<< this_particle.props.vel[1] << delimiter_verbose
					<< stokes_number << delimiter_verbose
					<< std::endl;
				num_line_print++;
			}
	} // end of sn_iter loop

	if (opt_verbose)
	{
		std::cout << std::endl;
		if (elapsed_time >= time_duration)
			std::cout << "\nMax allowed time duration is " << time_duration << " seconds and is already reached.\n\n";
	}

	// event: start
	const_cast<Event_Execution_Time&>(eet).start_event("report_hist");

	if (opt_output)
		this_particle.report_hist(
			dir_base_particle_hist_str.c_str(),
			file_extension,
			delimiter,
			true,
			fname_base);
	// #pragma endregion

	const_cast<Event_Execution_Time&>(eet).stop_event("report_hist");
	// event: stop

	//-----------------------------------------------------------------
	// Release memory
	//-----------------------------------------------------------------
	std::vector<int>().swap(cell_search_stack_circumferential_point);
	std::vector<double>().swap(dist_from_nodes_circumferential_point);
	std::vector<int>().swap(wall_zone_type_vec);
	
}


std::tuple<int, int> Particle_Tracking::eval_local_props(
	double *pos,
	Local_Flow_Props &lp,
	int props_eval_spiral_level,
	int flag_normal_return,
	int flag_out_of_domain)
{
	//-----------------------------------
	// Declarations
	//-----------------------------------
	double EPSILON = 1e-12;

	int flag_finding_encompassing_cell_normal_return = 0;
	int flag_finding_encompassing_cell_not_finding_encomp_cell = -1;
	int flag_finding_encompassing_cell_not_finding_close_node = -2;

	std::vector<int> cell_search_stack, cell_stack_for_props_eval, node_stack_for_props_eval;
	std::vector<double> dist_from_nodes;
	double sum_dist_inv;
	int sn_cell, sn_node, flag_finding_encompassing_cell;
	//-----------------------------------

	// #pragma region Finding the encompassing cell
	std::tie(sn_cell, dist_from_nodes, flag_finding_encompassing_cell, cell_search_stack) = grid->find_cell_encompassing_pos(
		pos,
		flag_finding_encompassing_cell_normal_return,
		flag_finding_encompassing_cell_not_finding_encomp_cell,
		flag_finding_encompassing_cell_not_finding_close_node);
	if (flag_finding_encompassing_cell == flag_finding_encompassing_cell_not_finding_close_node)
	{
		if (opt_verbose)
			std::cout << "`eval_local_props` | Particle out of domain (wrong release coordinates introduction?).\n\n";
		return std::make_tuple(-1, flag_out_of_domain);
	}

	if (flag_finding_encompassing_cell == flag_finding_encompassing_cell_not_finding_encomp_cell)
	{
		/*
		* Potential causes:
		* 1. Particle has crossed a non-wall boundary and is out of domain. 
		* Note: as of now, only wall boundary condition is implemented for particle motion.
		* 2. The applied boundary condition may not work properly.
		* 3. The code for finding encomp cell may have a bug.
		*/
		if (opt_verbose)
			std::cout << "`eval_local_props` | Particle out of domain.\n\n";
		return std::make_tuple(-1, flag_out_of_domain);
	}
	// #pragma endregion

	cell_stack_for_props_eval = grid->get_cells_near_cell(sn_cell, props_eval_spiral_level);
	node_stack_for_props_eval = grid->get_nodes_from_cells(cell_stack_for_props_eval);
	dist_from_nodes = grid->get_dist_from_nodes(pos, node_stack_for_props_eval);

	// weights based on distance from nodes of cell encompassing particle position
	sum_dist_inv = 0.;
	for (int ind_n = 0; ind_n < (int)node_stack_for_props_eval.size(); ind_n++)
		sum_dist_inv += 1 / (EPSILON + dist_from_nodes[ind_n]);

	// init
	lp.viscosity = 0;
	lp.density = 0;
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		lp.vel[ind_dir] = 0.;

	// accumulation
	for (int ind_n = 0; ind_n < (int)node_stack_for_props_eval.size(); ind_n++)
	{
		sn_node = node_stack_for_props_eval[ind_n];
		lp.viscosity += grid->node[sn_node].props.viscosity / (EPSILON + dist_from_nodes[ind_n]);
		lp.density += grid->node[sn_node].props.density / (EPSILON + dist_from_nodes[ind_n]);

		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			lp.vel[ind_dir] += grid->node[sn_node].props.vel[ind_dir] / (EPSILON + dist_from_nodes[ind_n]);
	}

	// avg.
	lp.viscosity /= sum_dist_inv;
	lp.density /= sum_dist_inv;
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		lp.vel[ind_dir] /= sum_dist_inv;

	return std::make_tuple(sn_cell, flag_normal_return);
}


std::tuple<double*, double*, int> Particle_Tracking::apply_bc(
	double* particle_center_pos_old,
	double* particle_center_pos_current,
	double* particle_vel_current,
	double diameter_p,
	std::vector<int> wall_zone_type,
	double elastic_reflection_coeff,
	int num_points_on_perimeter,
	const Event_Execution_Time &eet
)
{
	if (diameter_p == 0)
		num_points_on_perimeter = 1;

	// #pragma region declarations
	struct Deepest_Penetration
	{
		Deepest_Penetration(int val) : ndim(val)
		{
			pos_point = new double[ndim];
			pos_point_old = new double[ndim];
		}
		int ind_p = -1, sn_face = -1;
		double depth = -1;
		double *pos_point, *pos_point_old;

	private:
		int ndim;
	};

	//-------------------------------------------------------
	// Declarations
	//-------------------------------------------------------
	Deepest_Penetration deepest_penetration(ndim);

	std::vector<int> cell_search_stack;
	std::vector<double> dist_from_nodes;
	std::vector<int> face_stack_wall_type;
	double theta_p, this_penetration_depth, vel_normal_comp;
	int sn_cell, sn_face, flag_finding_encompassing_cell;
	//
	double *vec_displacement, *vec_perpendicular_to_face, 
		*pos_p, //position of circumferential points on particle
		*pos_p_old, //old position of circumferential points on particle
		*intersection_coords, *mirrored_point, *face_normal_vec,
		*pos_new, *vel_new;
	vec_displacement = new double[ndim];
	vec_perpendicular_to_face = new double[ndim];
	face_normal_vec = new double[ndim];
	mirrored_point = new double[ndim];
	pos_p = new double[ndim];
	pos_p_old = new double[ndim];
	pos_new = new double[ndim];
	vel_new = new double[ndim];
	intersection_coords = new double[ndim];
	//
	bool has_intersection;
	//

	// init pos and vel so they evaluate to the current vals 
	// in case no reflection is applied, the current pos and vel should be returned
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		pos_new[ind_dir] = particle_center_pos_current[ind_dir];
		vel_new[ind_dir] = particle_vel_current[ind_dir];
	}

	int flag_finding_encompassing_cell_normal_return = 0;
	int flag_finding_encompassing_cell_not_finding_encomp_cell = -1;
	int flag_finding_encompassing_cell_not_finding_close_node = -2;

	// dependent on mesh can be different from one case to another.
	// should be an attr of Particle_Tracking class, or a parameter of function
	//int wall_zone_type = 7;
	//-------------------------------------------------------
	// #pragma endregion

	// #pragma region Checking circumferential points and potentially find deepest penetration
	//pos_old = prt.hist.pos.back();

	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		vec_displacement[ind_dir] = particle_center_pos_current[ind_dir] - particle_center_pos_old[ind_dir];

	for (int ind_p = 0; ind_p < num_points_on_perimeter; ind_p++)
	{
		// #pragma region finding current and old positions of this circumferential point
		theta_p = ind_p * 2. * PI / num_points_on_perimeter;
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			pos_p[ind_dir] = particle_center_pos_current[ind_dir];
		pos_p[0] += diameter_p / 2. * cos(theta_p);
		pos_p[1] += diameter_p / 2. * sin(theta_p);

		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			pos_p_old[ind_dir] = pos_p[ind_dir] - vec_displacement[ind_dir];

		// event: start
		const_cast<Event_Execution_Time&>(eet).start_event("apply_bc | find_cell_encompassing_pos");

		// #pragma region Attempting to finding the encompassing cell
		std::tie(sn_cell, dist_from_nodes, flag_finding_encompassing_cell, cell_search_stack) = grid->find_cell_encompassing_pos(
			pos_p,
			flag_finding_encompassing_cell_normal_return,
			flag_finding_encompassing_cell_not_finding_encomp_cell,
			flag_finding_encompassing_cell_not_finding_close_node,
			eet
			);
		// #pragma endregion
		
		const_cast<Event_Execution_Time&>(eet).stop_event("apply_bc | find_cell_encompassing_pos");
		// event: stop

		// #pragma region Not finding close node
		if (flag_finding_encompassing_cell == flag_finding_encompassing_cell_not_finding_close_node)
		{
			if (opt_verbose)
			{
				std::cout << "\nError\nNot finding a close node for a circumferential point on \n"
					<< "\nparticle position: ("
					<< particle_center_pos_current[0] << ", "
					<< particle_center_pos_current[1] << ")"
					<< "\ncircumferential point position: ("
					<< pos_p[0] << ", "
					<< pos_p[1] << ")\n";
				std::cout << "A potential remedy could be to increase the `half_width_search_window` which currently is\n";
				std::cout << half_width_search_window << "\n";
			}
			return std::make_tuple(pos_new, vel_new, -1);//flag_return_bc=-1
		}
		// #pragma endregion

		// #pragma region Not being able to find encomp cell, i.e. pos out of domain
		if (flag_finding_encompassing_cell == flag_finding_encompassing_cell_not_finding_encomp_cell)
		{
			// event: start
			const_cast<Event_Execution_Time&>(eet).start_event("apply_bc | face intersection");

			// getting a stack of all faces of the wall type zone
			face_stack_wall_type = grid->get_faces_from_cells(cell_search_stack, wall_zone_type);

			// iterating over all faces of the wall type zone to check if pos
			for (int ind_f = 0; ind_f < (int)face_stack_wall_type.size(); ind_f++)
			{
				sn_face = face_stack_wall_type[ind_f];

				std::tie(has_intersection, intersection_coords) = grid->face_intersection_with_line_segment(
					sn_face,
					pos_p,
					pos_p_old);

				if (has_intersection)
				{
					// Currently, criterion is based on normal distance of point from face. The distance between the point and intersection may be another good option.
					this_penetration_depth = grid->distance_point_from_face(sn_face, pos_p);

					if (this_penetration_depth > deepest_penetration.depth)
					{
						deepest_penetration.ind_p = ind_p;
						deepest_penetration.sn_face = sn_face;
						deepest_penetration.depth = this_penetration_depth;
						for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
						{
							deepest_penetration.pos_point[ind_dir] = pos_p[ind_dir];
							deepest_penetration.pos_point_old[ind_dir] = pos_p_old[ind_dir];
						}
					}
				}
			}

			// event: stop
			const_cast<Event_Execution_Time&>(eet).stop_event("apply_bc | face intersection");
		}
		// #pragma endregion
	}
	// #pragma endregion

	// #pragma region Update particle pos and vel
	if (deepest_penetration.depth > 0)
	{
		// event: start
		const_cast<Event_Execution_Time&>(eet).start_event("apply_bc | reflection");

		mirrored_point = grid->mirror_point_wrt_face(deepest_penetration.sn_face, deepest_penetration.pos_point);

		//vec_perpendicular_to_face pointing towards the fluid domain
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			vec_perpendicular_to_face[ind_dir] = mirrored_point[ind_dir] - deepest_penetration.pos_point[ind_dir];

		// normal vec of face
		face_normal_vec = normal_vec(vec_perpendicular_to_face, ndim);

		// normal component of vel
		vel_normal_comp = dot(particle_vel_current, face_normal_vec);

		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			// update particle pos
			pos_new[ind_dir] += vec_perpendicular_to_face[ind_dir];

			// update particle vel for elastic reflection
			//if (opt_bc_update_vel)
			//	vel_new[ind_dir] -= 2. * vel_normal_comp * face_normal_vec[ind_dir];

			// generalized reflection scheme
			// The acceptable range of `elastic_reflection_coeff` is 0 to 1 (inclusive).
			// 0: inelastic reflection
			// 1: fully elastic reflection
			vel_new[ind_dir] -= 2. * elastic_reflection_coeff * vel_normal_comp * face_normal_vec[ind_dir];
		}
		
		// event: stop
		const_cast<Event_Execution_Time&>(eet).stop_event("apply_bc | reflection");
	}
	// #pragma endregion

	return std::make_tuple(pos_new, vel_new, 0);//flag_return_bc=0
}

void Particle_Tracking::report_particles_hist(const char *fname, const char *delimiter)
{
	if (opt_verbose)
		std::cout << "Perparing particle history file.\n";
	
	for (int ind_particle = 0; ind_particle < num_particles; ind_particle++)
		particle[ind_particle].report_hist(fname, delimiter);
	
	if (opt_verbose)
		std::cout << "End of prep.\n\n";
}


double* Particle_Tracking::calc_forces(
	Local_Flow_Props lp,
	double* vel_p,
	double diameter_p,
	double density_p
	)
{
	//declarations
	double* force_per_unit_mass_drag, * slip_vel, * slip_vel_unit_vec, * acc_p;
	slip_vel = new double[ndim];
	slip_vel_unit_vec = new double[ndim];
	force_per_unit_mass_drag = new double[ndim];
	acc_p = new double[ndim];
	
	double reynolds_p, slip_vel_mag, drag_coeff, drag_force_mag_per_unit_mass;
	//---------------------------------------

	//slip velocity: vel_fluid - vel_particle
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		slip_vel[ind_dir] = lp.vel[ind_dir] - vel_p[ind_dir];
	slip_vel_mag = std::sqrt(std::pow(slip_vel[0], 2.) + std::pow(slip_vel[1], 2.));
	slip_vel_unit_vec = normal_vec(slip_vel, ndim);

	
	//nondim numbers
	reynolds_p = lp.density * slip_vel_mag * diameter_p / lp.viscosity;

	//forces
	drag_coeff = 24. / reynolds_p * (1. + 0.15 * std::pow(reynolds_p, 0.687));
	drag_force_mag_per_unit_mass = 0.75 * lp.viscosity * drag_coeff * reynolds_p * slip_vel_mag
		/ density_p / std::pow(diameter_p, 2.);

	/*
	* The drag force magnitude is not needed
	drag_force_mag = 0.5 * local_props.density
		* std::pow(slip_vel_mag, 2.)
		* PI / 4. * std::pow(this_particle.props.diameter, 2.)
		* drag_coeff;
	*/

	//apply different forces; as of now, only drag included
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		//drag
		//force_per_unit_mass_drag[ind_dir] = taw_relax_p * slip_vel[ind_dir];
		force_per_unit_mass_drag[ind_dir] = drag_force_mag_per_unit_mass * slip_vel_unit_vec[ind_dir];
	}

	//calc acceleration of particle
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		acc_p[ind_dir] = 0;
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		//drag
		acc_p[ind_dir] += force_per_unit_mass_drag[ind_dir];
	}

	return acc_p;
}

//-------------------------------------------------------
// Getters
//-------------------------------------------------------

int Particle_Tracking::get_num_particles()
{
	return num_particles;
}

int Particle_Tracking::get_ndim()
{
	return ndim;
}
