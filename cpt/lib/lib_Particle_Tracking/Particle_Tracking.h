#pragma once
#include "Particle.h"
#include "Grid.h"

// event execution time
#include "Event_Execution_Time.h"

#ifndef struct_local_flow_props
#define struct_local_flow_props
struct Local_Flow_Props
{
	Local_Flow_Props(int val) : ndim(val)
	{
		vel = new double[ndim];
	}
	double density = -1, viscosity = -1;
	double *vel;

private:
	int ndim;
};
#endif // struct_local_flow_props


/**
Particle tracking class with methods to read and set up configurations of problem. 

//---------------------------------------------------------------------------------------------
// Stokes regimes-- different integration schemes & Stokes-dependent adaptive time step configs
//---------------------------------------------------------------------------------------------
Small Stokes time integration scheme is applied if `(stokes_number < stokes_number_threshold or opt_time_integration_for_small_stokes)`
By default, `opt_time_integration_for_small_stokes=false`, wherein small Stokes time integration scheme is applied only if `stokes_number < stokes_number_threshold`
If `opt_time_integration_for_small_stokes` is set to true,  small Stokes time integration scheme is enforced throughout the simulation.

`fixed_time_step` is used if adaptive time step is not applicable.
`L_ref` is needed for Stokes-dependent adaptive time step. In case `L_ref` is provided:
	1. If `Vel_ref` is provided, it is used to evaluate a fixed Stokes number
	2. Otherwise, local cell-based velocity magnitude is used to evaluate a local Stokes number
*/
class Particle_Tracking
{

public:

	/**
	
	@param opt_original_ASG If true, the original variant of ASG search is used to locate hosting cell.Otherwise, the modified variant
	(two - step spiral) is used.
	@param opt_bin_cell_check_center_only Two implementation approaches to add cells to bins to be used in the original ASG variant:
	1. (`opt_bin_cell_check_center_only`=true) center of cell to be examined: each cell resides only in one bin.
	2. (`opt_bin_cell_check_center_only`=false) each node of cell to be examined: each cell can reside in one or multiple bins.
	@param half_width_search_window Search window bins when searching using ASG or its modified variant.
	@param encompassing_cell_spiral_level Adding extra layer of cells to search stack when looking for the cell encompassing a position.
	@param nondim_criterion_area_ratio Nondim criterion on area ratio when checking if a position is inside a constituent triangle of a cell.
	@param props_eval_spiral_level Number of stacks of cells around an encompassing cell to be considered for local properties interpolation. 
	`props_eval_spiral_level=0` means that only the encompassing cell is considered for properties interpolation.
	@param opt_output Whether to write files, which can be turned off when evaluating the computational performance of algorithm.
	@param opt_verbose Whether to log using std::cout, which can be turned off when evaluating the computational performance of algorithm.
	@param int num_threads Number of threads to be used for parallel processing.
	*/
	Particle_Tracking(
		bool opt_original_ASG = false,
		bool opt_bin_cell_check_center_only = false,
		int half_width_search_window = 1,
		int encompassing_cell_spiral_level = 1,
		double nondim_criterion_area_ratio = 1e-11,
		int props_eval_spiral_level=0,
		bool opt_output = true,
		bool opt_verbose = true,
		int num_threads=1
	);

	~Particle_Tracking();

	/**
	Reading a configuration file for particles to be introduced into the system.

	The following shows an example of particles configuration file with single space as delimiter:

		pos.x pos.y diameter density vel.x vel.y
		1e-6 5e-6 1e-6 1e3 0 0
		4e-6 7.5e-6 2e-6 1e3 0 0

	The following shows an example of particles configuration file with three spaces as delimiter:

		pos.x   pos.y   diameter   density   vel.x   vel.y 
		1e-6   5e-6   1e-6   1e3   0   0
		4e-6   7.5e-6   2e-6   1e3   0   0
		4e-6   2.5e-6   0.5e-6   1e3   -1e-3   -1e-3 

	The following shows an example of particles configuration file with a tab ("\t") as delimiter:

		pos.x	pos.y	diameter	density	vel.x	vel.y	opt_release_at_flow_speed
		1e-6	4.9e-06	1e-06	1.0	0	0	1
		1e-6	4.95e-06	1e-06	1.0	0	0	1
		1e-6	5e-06	1e-06	1.0	0	0	1
		1e-6	5.05e-06	1e-06	1.0	0	0	1
		1e-6	5.1e-06	1e-06	1.0	0	0	1

	@param path_particle_file filename.
	@param delimiter delimiter used in the file separating columns.
	*/
	void read_particle_file(
		const char* path_particle_file, 
		const char* delimiter=","
		);
	
	/**
	Reading and setting up grid and solution files.

	@param path_mesh (path including) filename for grid/mesh, e.g., "mesh", "./mesh", "./inp/mesh", etc.
	@param path_sol_nodes (path including) filename for solution of continuous phase solution fields on nodes of grid, e.g., "sol_nodes", "./sol_nodes", "./inp/sol_nodes", etc.
	@param path_sol_cells (path including) filename for solution of continuous phase solution fields on cells of grid, e.g., "sol_cells", "./sol_cells", "./inp/sol_cells", etc. 
	Typically, one of path_sol_nodes and path_sol_cells is sufficient.
	@param path_report_mesh (path including) filename for reporting the parsed mesh that can be used to make sure the mesh is read correctly. Default: "./report/mesh".
	@param path_report_mesh_ordered_connectivity (path including) filename for reporting the parsed mesh after ordering the connectivity of faces. 
	Default: "./report/mesh_ordered_connectivity".
	@param path_report_sol_nodes (path including) filename for reporting the continuous phase solution fields on nodes to make sure the solution is read correctly. 
	Default: "./report/sol_nodes".
	@param path_report_sol_cells (path including) filename for reporting the continuous phase solution fields on cells to make sure the solution is read correctly. 
	Default: "./report/sol_cells".
	@param path_report_bins (path including) filename for reporting the bins and nodes of main mesh residing within each bin. Default: "./report/bins".
	@param opt_read_sol_at_cells_center If true, continuous phase solution fields on cells need to be provided. Default: false.
	@param n_bins Dimensions of Auxiliary Structured Grid (ASG). Default: {20, 20}.
	@param allowed_tolerance_node_position_snap Tolerance used for detecting a node. Any location within this distance from a node is considered the same node. Default: 1e-12 (meter).
	@param nondim_allowed_tolerance_cell_center_position_snap Maximum non-dimensional distance of a cell center provided within the file related to continuous phase solution fields 
	for cells from that calculated by considering the nodes of the cell. Only relevant when `opt_read_sol_at_cells_center` is true. 
	If the distance is larger than this threshold, the grid and solution files are inconsistent. Default: 0.1.
	@param strict_mode Whether to run a quality check to ensure the consistency of grid and soluion files. 
	If true, it inspects cell center coordinates provided within solution file for cells to make sure the right cell exists amongst the cells provided in the mesh file. 
	Only relevant when `opt_read_sol_at_cells_center` is true. Default: true.
	@param opt_report_mesh If yes, a report file is produced for the grid, which can be used to ensure the correctness of reading the provided mesh file. Default: false.
	@param opt_report_mesh_ordered If yes, a report file is produced for the grid after ordering the connectivity of faces, which can be used for troubleshooting as needed. 
	Default: false.
	@param opt_report_sol_nodes If yes, a report file is produced for the continuous phase solution fields at grid nodes, which can be used to ensure the correctness of 
	reading the solution file. Default: false.
	@param opt_report_sol_cells If yes, a report file is produced for the continuous phase solution fields at grid cells, which can be used to ensure the correctness of 
	reading the solution file. Default: false.
	@param opt_report_bins If yes, a report file is produced for ASG, which can be used for troubleshooting and pos-simulation analysis. Default: false.
	@param opt_report_nodes_of_zones_one_file_per_zone If yes, a report file is produced for nodes belonging to each zone (zones defined in the mesh file), which can be used for 
	troubleshooting and pos-simulation analysis. Here, each file is related to a specific zone. Default: false.
	@param opt_report_nodes_of_zones_single_file If yes, a report file is produced for nodes belonging to each zone (zones defined in the mesh file), which can be used for 
	troubleshooting and pos-simulation analysis. Here, each line of the file is related to a specific zone. Default: false.
	@param dir_nodes_of_zones Path to directory to save the report file(s) for nodes of znoes. Default: <dir_base>/report.
	@param fname_base_nodes_of_zones Base filename for the report file(s) for nodes of zones. Default: "nodes_of_zone".
	@param file_extension_nodes_of_zones Extension for the report file(s) for nodes of zones. Default: "", i.e., files would have no extension.
	@param delimiter_nodes_of_zones Delimiter separating node indices in nodes of zones files. Default: ",".
	@param report_every_num_line Every this number of lines, the code will verbose when reading and setting up the files. Default: 100,000.
	@param bin_len_dimensionless Dimensionless length of sides of bins of ASG w.r.t. grid dimensions. To be applied if any item of `n_bins` is zero.
	*/
	void set_grid_and_sol(
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
		double bin_len_dimensionless = 1.5
	);

	/**
	Reading and setting up grid file.

	@param path_mesh (path including) filename for grid/mesh, e.g., "mesh", "./mesh", "./inp/mesh", etc.
	@param path_report_mesh (path including) filename for reporting the parsed mesh that can be used to make sure the mesh is read correctly. Default: "./report/mesh".
	@param path_report_mesh_ordered_connectivity (path including) filename for reporting the parsed mesh after ordering the connectivity of faces. 
	Default: "./report/mesh_ordered_connectivity".
	@param path_report_bins (path including) filename for reporting the bins and nodes of main mesh residing within each bin. Default: "./report/bins".
	@param n_bins Dimensions of Auxiliary Structured Grid (ASG). Default: {20, 20}.
	@param opt_report_mesh If yes, a report file is produced for the grid, which can be used to ensure the correctness of reading the provided mesh file. Default: false.
	@param opt_report_mesh_ordered If yes, a report file is produced for the grid after ordering the connectivity of faces, which can be used for troubleshooting as needed. 
	Default: false.
	@param opt_report_bins If yes, a report file is produced for ASG, which can be used for troubleshooting and pos-simulation analysis. Default: false.
	@param opt_report_nodes_of_zones_one_file_per_zone If yes, a report file is produced for nodes belonging to each zone (zones defined in the mesh file), 
	which can be used for troubleshooting and pos-simulation analysis. Here, each file is related to a specific zone. Default: false.
	@param opt_report_nodes_of_zones_single_file If yes, a report file is produced for nodes belonging to each zone (zones defined in the mesh file), 
	which can be used for troubleshooting and pos-simulation analysis. Here, each line of the file is related to a specific zone. Default: false.
	@param dir_nodes_of_zones Path to directory to save the report file(s) for nodes of znoes. Default: <dir_base>/report.
	@param fname_base_nodes_of_zones Base filename for the report file(s) for nodes of zones. Default: "nodes_of_zone".
	@param file_extension_nodes_of_zones Extension for the report file(s) for nodes of zones. Default: "", i.e., files would have no extension.
	@param delimiter_nodes_of_zones Delimiter separating node indices in nodes of zones files. Default: ",".
	@param report_every_num_line Every this number of lines, the code will verbose when reading and setting up the files. Default: 100,000.
	@param bin_len_dimensionless Dimensionless length of sides of bins of ASG w.r.t. grid dimensions. To be applied if any item of `n_bins` is zero.
	*/
	void set_grid(
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
		double bin_len_dimensionless=1.5
		);

	/**
	Reading and setting up solution files.

	@param path_sol_nodes (path including) filename for solution of continuous phase solution fields on nodes of grid, e.g., "sol_nodes", "./sol_nodes", "./inp/sol_nodes", etc.
	@param path_sol_cells (path including) filename for solution of continuous phase solution fields on cells of grid, e.g., "sol_cells", "./sol_cells", "./inp/sol_cells", etc. 
	Typically, one of path_sol_nodes and path_sol_cells is sufficient.
	@param path_report_sol_nodes (path including) filename for reporting the continuous phase solution fields on nodes to make sure the solution is read correctly. 
	Default: "./report/sol_nodes".
	@param path_report_sol_cells (path including) filename for reporting the continuous phase solution fields on cells to make sure the solution is read correctly. 
	Default: "./report/sol_cells".
	@param opt_read_sol_at_cells_center If true, continuous phase solution fields on cells need to be provided. Default: false.
	@param allowed_tolerance_node_position_snap Tolerance used for detecting a node. Any location within this distance from a node is considered the same node. Default: 1e-12 (meter).
	@param nondim_allowed_tolerance_cell_center_position_snap Maximum non-dimensional distance of a cell center provided within the file related to continuous phase solution fields 
	for cells from that calculated by considering the nodes of the cell. Only relevant when `opt_read_sol_at_cells_center` is true. 
	If the distance is larger than this threshold, the grid and solution files are inconsistent. Default: 0.1.
	@param strict_mode Whether to run a quality check to ensure the consistency of grid and soluion files. 
	If true, it inspects cell center coordinates provided within solution file for cells to make sure the right cell exists amongst the cells provided in the mesh file. 
	Only relevant when `opt_read_sol_at_cells_center` is true. Default: true.
	@param opt_report_sol_nodes If yes, a report file is produced for the continuous phase solution fields at grid nodes, which can be used to ensure the correctness of 
	reading the solution file. Default: false.
	@param opt_report_sol_cells If yes, a report file is produced for the continuous phase solution fields at grid cells, which can be used to ensure the correctness of 
	reading the solution file. Default: false.
	@param report_every_num_line Every this number of lines, the code will verbose when reading and setting up the files. Default: 100,000.
	*/
	void set_sol(
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
		int report_every_num_line);

	/**
	Setting up a streakline for injection of particles.

	@param n_particles Number of particles on streakline to be injected.
	@param pos_end_1 Position of first end of streakline.
	@param pos_end_2 Position of second end of streakline.
	@param vel_end_1 Velocity of particle on first end of streakline.
	@param vel_end_2 Velocity of particle on second end of streakline. Particles velocity varies lineary between those associated with first and second ends of streakline.
	@param range_diameter Range of particles diameter linearly varying between first and second ends of streakline. For example, {1e-6, 40e-6} for diameters 
	between 1-micron and 40-micron.
	@param range_density Range of particles density linearly varying between first and second ends of streakline. For example, {1, 1e3} for densities between 1 and 1e3 kg/m^3.
	@param opt_release_at_flow_speed Whether to release particles at fluid flow speed. It can override other specified velocities. Default: false.
	*/ 
	void set_streakline(
		int n_particles,
		double* pos_end_1,
		double* pos_end_2,
		double* vel_end_1,
		double* vel_end_2,
		double range_diameter[2],
		double range_density[2], 
		bool opt_release_at_flow_speed = false
	);

	/**
	\overload void Particle_Tracking::set_streakline(int n_particles,double* pos_end_1,double* pos_end_2,double range_diameter[2],double range_density[2],bool opt_release_at_flow_speed)

	In particular, vel_end_1 and vel_end_2 are not passed here, which will be set to zero.
	 */
	void set_streakline(
		int n_particles,
		double* pos_end_1,
		double* pos_end_2,
		double range_diameter[2],
		double range_density[2],
		bool opt_release_at_flow_speed = false
	);

	/**
	Setting up a point source for injection of particles.

	@param n_particles Number of particles to be injected.
	@param source_pos Position of point source of injection.
	@param range_vel Range of particles velocity varies lineary between those associated with first and second specified values. For example, {{0,0}, {0.001, 0}}.
	@param range_diameter Range of particles diameter linearly varying between first and second specified values. For example, {1e-6, 40e-6} for diameters between 
	1-micron and 40-micron.
	@param range_density Range of particles density linearly varying between first and second specified values. For example, {1, 1e3} for densities between 1 and 1e3 kg/m^3.
	@param opt_release_at_flow_speed Whether to release particles at fluid flow speed. It can override other specified velocities. Default: false.
	*/
	void set_particle_source(
		int n_particles,
		double *source_pos,
		double **range_vel,
		double range_diameter[2], 
		double range_density[2],  
		bool opt_release_at_flow_speed = false
	);

	/**
	\overload void Particle_Tracking::set_particle_source(int n_particles,double *source_pos,double range_diameter[2],double range_density[2],bool opt_release_at_flow_speed = false)

	In particular, range_vel is not passed here, for which zero velocities will be assumed.
	 */
	void set_particle_source(
		int n_particles,
		double *source_pos,
		double range_diameter[2],
		double range_density[2],
		bool opt_release_at_flow_speed = false);

	/**
	Setting up particles in the system from given arrays of attributes, e.g., diameter, density, position, velocity, etc.

	@param n_particles Number of particles.
	@param dia Array of diameters to be considered for particles to be registered into the system.
	@param rho Array of densities to be considered for particles to be registered into the system.
	@param pos Array of initial positions to be considered for particles to be registered into the system.
	@param vel Array of initial velocities to be considered for particles to be registered into the system.
	@param acc Array of initial acceleratoins to be considered for particles to be registered into the system.
	@param force Array of initial forces to be considered for particles to be registered into the system. 
	The force attribute is currently not fully connected with the time integration, as the current formulation is mainly based on force per unit mass, i.e., acceleration.
	@param opt_release_at_flow_speed Whether to release particles at fluid flow speed. It can override other specified velocities. Default: false.
	*/
	void set_particles(//most generic function to set particle(s)
		int n_particles,
		double *dia,
		double *rho,
		double **pos,
		double **vel,
		double **acc,
		double **force,
		bool opt_release_at_flow_speed = false
	);

	/**
	\overload void Particle_Tracking::void set_particles(double *single_pos,bool opt_release_at_flow_speed)

	It sets up a single pointwise particle (density and diameter equal to 0) with an option to release at flow speed or stationary.
	*/
	void set_particles(
		double *single_pos,
		bool opt_release_at_flow_speed = false);

	/**
	\overload void Particle_Tracking::void set_particles(double *single_pos,double diameter,double density,bool opt_release_at_flow_speed)

	It sets up a single set a single particle of given diameter and density with an option to release at flow speed or stationary.
	*/
	void set_particles(
		double *single_pos,
		double diameter,
		double density,
		bool opt_release_at_flow_speed = false);

	/**
	\overload void Particle_Tracking::void set_particles(double* single_pos,double diameter,double density,double* single_vel)

	It sets up a single set a single particle of given diameter and density with a given velocity.
	*/
	void set_particles(
		double* single_pos,
		double diameter,
		double density,
		double* single_vel);

	/**
	Main method for tracking particles.

	@param dir_base_particle_hist Path to directory for reporting particles history file. Default: <dir_base>/"particle".
	@param time_duration Maximum time duration for tracking particles.
	@param time_step_over_taw_relax Non-dimensional time step w.r.t. local instantaneous fluid flow relaxation time to enable time integration with adaptive time step.
	@param wall_zone_type Array of zones that are of type wall.
	@param particle_id_vec Vector of particle indices to be tracked.
	@param L_ref Charateristic length of system.
	@param Vel_ref Charateristic velocity magnitude of system.
	@param num_wall_zone_type Number of zones of type of wall, i.e., length of wall_zone_type.
	@param n_particles_to_track Number of particles from the beginning of particle_id_vec to be tracked. Default: -1. A value of 0 or negative implies tracking all particles.
	@param max_iter Maximum number of time steps.
	@param verbose_every_num_step Every this number of steps, the code will verbose.
	@param elastic_reflection_coeff Elastic reflection coefficient. 0 corresponds to a pure inelastic reflection and 1 corresponds to a pure elastic reflection. Default: 1.
	@param opt_time_integration_for_small_stokes If true, time integration will be acceleration-independent throughout the simulation, because of which inertia effects are not captured. 
	Only suitable for low-Stokes regime, wherein particle precisely follows streamlines. Default: false.
	@param stokes_number_threshold A threshold defining the transition from low-Stokes regime to high-Stokes regime. Default: 0.1.
	@param fixed_time_step A fixed time step in case adaptive time step is not configurable. Default: 1e-3 (second).
	@param num_points_on_perimeter Number of circumferential points on perimeter of particle. Default: 8.
	@param opt_one_step_euler If true, one-step explicit Euler scheme will be used for time integration. Default: true.
	@param opt_rk4 `opt_one_step_euler` set to false, fourth-order Runge-Kutta technique will be used for time integration. Default: false.
	@param fname_base Base filename to be used when saving trajectory of particles. Trajectory filenames would be of the format of `<fname_base>_<particle_id>(.<file_extension>)`. 
	Default: "particle". 
	@param file_extension Extension for particle trajectory files. Default: "", i.e., no extension.
	@param delimiter Delimiter for separating attributes to be written in each row, e.g., iter,time,x,y,vel_x,vel_y,.... Default: ",".
	*/
	void track_particle(
		const char* dir_base_particle_hist,
		double time_duration,
		double time_step_over_taw_relax,
		int* wall_zone_type,	
		int* particle_id_vec,
		double L_ref = 0,	
		double Vel_ref = 0, 
		int num_wall_zone_type = 1, 
		int n_particles_to_track = -1, 
		int max_iter = 1e3,
		int verbose_every_num_step = 1,
		double elastic_reflection_coeff=1, 
		bool opt_time_integration_for_small_stokes=false,
		double stokes_number_threshold = 1e-1,
		double fixed_time_step=1e-3,
		int num_points_on_perimeter = 8,
		bool opt_one_step_euler=true,
		bool opt_rk4=false,
		const char *fname_base = "particle",
		const char *file_extension = "",
		const char *delimiter = ","
	);

	void track_single_particle(
		Particle& this_particle,
		const char* dir_base_particle_hist,
		double time_duration,
		double time_step_over_taw_relax,
		int* wall_zone_type, 
		double L_ref = 0,	
		double Vel_ref = 0, 
		int num_wall_zone_type = 1, 
		int max_iter = 1e3,
		int verbose_every_num_step = 1,
		double elastic_reflection_coeff = 1, 
		bool opt_time_integration_for_small_stokes = false,
		double stokes_number_threshold = 1e-1,
		double fixed_time_step = 1e-3,
		int num_points_on_perimeter = 8,
		bool opt_one_step_euler = true,
		bool opt_rk4 = false,
		const char* fname_base = "particle",
		const char* file_extension = "",
		const char* delimiter = ",",
		const Event_Execution_Time& eet = Event_Execution_Time(false)	// Event Execution Time object to measure performance
	);

	/**
	Calculate forces exerted on a particle with given diameter, density, and velocity subject to a given set of local fluid flow properties.

	@param lp A structure datatype providing the local fluid flow properties.
	@param vel_p Velocity of particle.
	@param diameter_p Diameter of particle.
	@param density_p Density of particle.
	*/
	double* calc_forces(
		Local_Flow_Props lp,
		double* vel_p,
		double diameter_p,
		double density_p
		);

	/**
	Applying particle-wall collision model.

	@param particle_center_pos_old Old position of paticle.
	@param particle_center_pos_current Current position of paticle.
	@param particle_vel_current Current velocity of paticle.
	@param diameter_p Diameter of paticle.
	@param wall_zone_type Array of zones that are of type wall.
	@param elastic_reflection_coeff Elastic reflection coefficient. 0 corresponds to a pure inelastic reflection and 1 corresponds to a pure elastic reflection. Default: 1.
	@param num_points_on_perimeter Number of circumferential points on perimeter of particle. Default: 8.
	@param Event Execution Time object to measure performance
	*/
	std::tuple<double*, double*, int> apply_bc(
		double* particle_center_pos_old,
		double* particle_center_pos_current,
		double* particle_vel_current,
		double diameter_p,
		std::vector<int> wall_zone_type,
		double elastic_reflection_coeff,
		int num_points_on_perimeter,
		const Event_Execution_Time &eet=Event_Execution_Time(false)
	);


	/**
	Evaluate local properties of fluid flow by interpolating the corresponding values from those related to the nearby nodes.

	@param pos Target position.
	@param local_props A structure datatype to store local fluid flow properties.
	@param props_eval_spiral_level Number of stacks of cells around an encompassing cell to be considered for local properties interpolation. 
	`props_eval_spiral_level=0` means that only the encompassing cell is considered for properties interpolation.
 	*/
	std::tuple<int, int> eval_local_props(
		double *pos,
		Local_Flow_Props &local_props,
		int props_eval_spiral_level,
		int flag_normal_return = 0,
		int flag_out_of_domain = -1);

	/**
	Report history of all particles.

	@param fname Filename to be used when saving trajectory of particle. 
	@param delimiter Delimiter for separating attributes to be written in each row, e.g., iter,time,x,y,vel_x,vel_y,.... Default: ",".
	*/ 
	void report_particles_hist(
		const char *fname, 
		const char *delimiter = ","
	);

	
	//-----------------------------------------------------------
	// getters for data attributes
	//-----------------------------------------------------------

	/**
	Returns the number of particles in the system.
	*/
	int get_num_particles();

	/**
	Returns the number of dimensions in the system.
	*/
	int get_ndim();

	//-----------------------------------------------------------
	// data attributes
	//-----------------------------------------------------------
	Grid *grid; /**<Grid for continuous phase solution fields.*/

	std::vector<Particle> particle; /**<Particles available in the system.*/

	//----------------------------------
	// config algorithm
	//----------------------------------
	bool opt_original_ASG; 
	bool opt_bin_cell_check_center_only;
	int half_width_search_window; 
	int encompassing_cell_spiral_level;
	double nondim_criterion_area_ratio;
	int props_eval_spiral_level; 

	//----------------------------------
	// config verbose
	//----------------------------------
	bool opt_verbose;
	bool opt_output;

	//----------------------------------
	// config OpenMP
	//----------------------------------
	int num_threads;

private:
	int ndim;
	int num_particles;
};
