#pragma once

#include "Cell.h"
#include "Node.h"
#include "Face.h"
#include "Zone.h"
#include "Bin.h"
#include "Props.h"

#include <tuple>
#include <vector>
#include <algorithm>

// event execution time
#include "Event_Execution_Time.h"

//parallel loop
#include <omp.h>


/**
Grid class with methods to read and set up configurations of mesh and continuous phase solution fields.

Important:

In Ansys Fluent, the index of cell, node, face, and zone starts from 1. To be consistent with that, we add a dummy item at 
index of 0 for the corresponding arrays.

At the same time, the code is devised so that `cell[0]` holds the information related to the "boundary cells".

For example, cell[0].face and cell[0].node holds the index of faces and nodes that are related to "boundary cells", respectively.
These boundary cells show up with index of 0 in other attributes such as face and node.

For example, face[sn_face].cell may hold the following indices for cells related to the given face: 882529 0, wherein 
the index of 0 indicates that the face is on a boundary.

Similarly, node[sn_node].cell may hold the following indices for cells related to the given node: 882529 882530 0, wherein 
the index of 0 indicates that the node is on a boundary.
*/
class Grid
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
	@param opt_output Whether to write files, which can be turned off when evaluating the computational performance of algorithm.
	@param opt_verbose Whether to log using std::cout, which can be turned off when evaluating the computational performance of algorithm.
	*/
	Grid(bool opt_original_ASG=false,
		bool opt_bin_cell_check_center_only=false,
		int half_width_search_window=1,
		int encompassing_cell_spiral_level = 1, 
		double nondim_criterion_area_ratio = 1e-11,
		bool opt_output = true,
		bool opt_verbose = true
		);
	~Grid();

	/**
	Read Ansys Fluent mesh file.

	@param fname Filename.
	@param report_every_num_line Verbose every this number of lines. Default: 1.
	*/
	void read_ansys_fluent_mesh_line_parsing(const char *fname, int report_every_num_line = 1); // fast as it parses line-by-line

	/**
	Read Ansys Fluent solution file at nodes.

	@param fname Filename.
	@param allowed_tolerance_node_position_snap Tolerance used for detecting a node. Any location within this distance from a node is considered the same node. Default: 1e-11 (meter).
	@param report_every_num_line Verbose every this number of lines. Default: 1.
	@param delimiter Delimiter used in Ansys file.
	*/
	void read_ansys_fluent_sol_at_nodes(
		const char *fname,
		double allowed_tolerance_node_position_snap = 1e-11,
		int report_every_num_line = 1,
		const char *delimiter = ",");

	/**
	Read Ansys Fluent solution file at cells.

	@param fname Filename.
	@param nondim_allowed_tolerance_cell_center_position_snap Maximum non-dimensional distance of a cell center provided within the file related to continuous phase solution fields
	for cells from that calculated by considering the nodes of the cell. If the distance is larger than this threshold, the grid and solution files are inconsistent. Default: 0.1.
	@param report_every_num_line Verbose every this number of lines. Default: 1.
	@param strict_mode Whether to run a quality check to ensure the consistency of grid and soluion files.
	If true, it inspects cell center coordinates provided within solution file for cells to make sure the right cell exists amongst the cells provided in the mesh file.
	Default: true.
	@param delimiter Delimiter used in Ansys file.
	*/
	void read_ansys_fluent_sol_at_cells(
		const char *fname,
		double nondim_allowed_tolerance_cell_center_position_snap = 1e-1, // distance nondim by cell size
		int report_every_num_line = 1,
		bool strict_mode = true,
		const char *delimiter = ",");


	/**
	Sets up bins of ASG.

	@param n_bins Number of bins of ASG in each direction.
	@param fname Filename for reporting the ASG.
	@param opt_report_bins If true, a report of ASG will be produced.
	@param bin_len_dimensionless Dimensionless length of sides of bins of ASG w.r.t. grid dimensions. To be applied if any item of `n_bins` is zero.
	@param delimiter Delimiter to be used for separating items in each row of report file.
	*/
	void set_up_bins(
		int *n_bins,
		const char *fname,
		bool opt_report_bins,
		double bin_len_dimensionless = 1.5,
		const char *delimiter = ","
	);

	/**
	Determines the bounding box of grid so that ASG can be constructed accordingly.
	*/
	void set_up_bounding_box();

	/**
	A setter function to set the data field `bin_num`.

	@param n_bins Number of bins of ASG in each direction.
	*/
	void set_bin_num(int *n_bins);

	/**
	A setter function to set the data field `bin_num`.

	@param bin_len Desired length of sides of bins of ASG in each direction.
	*/
	void set_bin_num(double bin_len);

	/**
	Produce a report for ASG.

	@param fname Filename for reporting the ASG.
	@param delimiter Delimiter to be used for separating items in each row of report file.
	*/
	void report_bins(const char *fname, const char *delimiter = ",");
	
	/**
	Initialize the ASG.
	*/
	void init_bins();

	/**
	Determine the nodes residing within each bin of ASG.
	*/
	void populate_bins();

	/**
	Produce a report for grid.

	@param fname Filename of the report.
	*/
	void report_mesh(const char *fname);
	
	/**
	Produce a report for continuous phase solution fields at nodes.

	@param fname Filename of the report.
	*/
	void report_sol_nodes(const char *fname, const char *delimiter = ",");

	/**
	Produce a report for continuous phase solution fields at cells.

	@param fname Filename of the report.
	*/
	void report_sol_cells(const char *fname, const char *delimiter = ",");

	/**
	Produce a report for nodes belonging to each zone (zones defined in the mesh file).
	Each line of the file will be related to a specific zone.

	@param dir Path to directory where report needs to be created.
	@param fname_base Base filename for the report file.
	@param file_extension Extension for the report file.
	@param delimiter Delimiter separating node indices within the report. Default: ",".
	*/
	void report_nodes_of_zones_single_file(//generates one file, wherein each line corresponds to a zone
		const char* dir,
		const char* fname_base,
		const char* file_extension,
		const char* delimiter = ",");

	/**
	Produce report(s) for nodes belonging to each zone (zones defined in the mesh file). A report file is produced for nodes belonging to each zone.

	@param dir Path to directory where report(s) need to be created.
	@param fname_base Base filename for the report file(s).
	@param file_extension Extension for the report file(s).
	@param delimiter Delimiter separating node indices within the report(s). Default: ",".
	*/
	void report_nodes_of_zones_one_file_per_zone(//generates one file per zone
		const char* dir,
		const char* fname_base,
		const char* file_extension,
		const char* delimiter = ",");


	/**
	Determines whether a position is within a convex polygon by splitting it into triangles and check for each triangle individually.
	
	It returns true if `pos` is inside cell #`sn_cell`. Otherwise, it returns false.

	@param pos Target position.
	@param sn_cell Index of cell of interest.
	@param nondim_criterion_area_ratio A non-dimensional area criterion for detecting whether a position falls inside a polygon. Default: 1e-11 (square meter).
	*/
	bool pos_inside_cell(double *pos,
						 int sn_cell,
						 double nondim_criterion_area_ratio = 1e-11); // true if pos is within cell, otherwise false


	/**
	To calculate and set fluid flow properties at the center of cell by applying interpolation to the corresponding values on its nodes.

	@param pos If yes, position at cell center will be set to the corresponsing value obtained from interpolation.
	@param vel If yes, velocity at cell center will be set to the corresponsing value obtained from interpolation.
	@param density If yes, density at cell center will be set to the corresponsing value obtained from interpolation.
	@param pressure If yes, pressure at cell center will be set to the corresponsing value obtained from interpolation.
	@param viscosity If yes, viscosity at cell center will be set to the corresponsing value obtained from interpolation.
	*/
	void calc_cell_props_from_nodes(
		bool pos = true,
		bool vel = false,
		bool density = false,
		bool pressure = false,
		bool viscosity = false);

	/**
	Returns area of cell.

	@param sn_cell Index of cell of interest.
	*/
	double calc_cell_area(int sn_cell);

	/**
	Calculate and set auxiliary properties of cell, e.g., vel_mag, and store it in a structure datatype `props_aux`.
	*/
	void calc_cells_aux_props();

	/**
	Calculate and set area of all cells.
	*/
	void calc_cells_area();

	/**
	Order faces and nodes on each cell based on right-hand rule.
	*/
	void order_cell_connectivity();

	/**
	Find the closest node to a given position (`pos`). The nodes to be inspected include those residing within the central bin of ASG (wherein the given position `pos` resides) 
	and bins around the central bin according to a given 
	half width search window (`half_width_search_window`).

	@param pos Position of interest.
	@param Event Execution Time object to measure performance
	*/
	std::tuple<int, double, int> find_node_close_to_pos(
		double *pos,
		double allowed_tolerance_node_position_snap,
		int flag_normal_return = 0,						 // flag normal return: found node withint the specified range
		int flag_far_to_snap = -1,						 // found node outside the specified range
		int flag_out_of_domain_or_small_search_win = -2, // flag out of domain or small search window
		int flag_out_of_bbox = -3,						 // flag out of bounding box
		const Event_Execution_Time & eet = Event_Execution_Time(false)	// Event Execution Time object to measure performance
	);
	
	//to be used in read solution at nodes
	std::tuple<int, double, int> find_node_close_to_pos_spiral_algorithm(
		double* pos,
		double allowed_tolerance_node_position_snap,
		int flag_normal_return = 0,						 // flag normal return: found node withint the specified range
		int flag_far_to_snap = -1,						 // found node outside the specified range
		int flag_out_of_domain_or_small_search_win = -2, // flag out of domain or small search window
		int flag_out_of_bbox = -3						 // flag out of bounding box
	);

	/**
	Applying the two-step search algorithm to locate hosting cell for a given position.

	@param pos Target position for which a hosting cell is sought.
	@param encompassing_cell_spiral_level Maximum number of stacks of cells around the central stack of cells to be considered when searching to locate hosting cell.
		A value of 0 means only cells that the closest node to `pos` is one of their vertices will be considered for search (central stack of cells).
		By increasing `spiral_level`, more outermost stack of cells will be added to the considered cells.
		Default: 4. The default value is conservative and allows dealing with low quality grids with highly skewed cells. For high quality grids, a value as low as 1 may be applied.
	@param nondim_criterion_area_ratio A non-dimensional area criterion for detecting whether a position falls inside a polygon. Default: 1e-11 (square meter).
	@param Event Execution Time object to measure performance
	*/
	std::tuple<int, std::vector<double>, int, std::vector<int>> find_cell_encompassing_pos(
		double *pos,
		int flag_normal_return = 0,					// normal return flag
		int flag_not_finding_encomp_cell = -1,		// flag of failing to find encompassing cell around the found close node
		int flag_not_finding_close_node = -2,		// flag not finding any close node
		const Event_Execution_Time& eet= Event_Execution_Time(false)	// Event Execution Time object to measure performance
	);

	/**
	Returns a vector of distances from a point to a given vector of nodes

	@param pos Position of interest.
	@param node_stack Vector consisting of index of nodes of interest.
	*/
	std::vector<double> get_dist_from_nodes(double *pos, std::vector<int> node_stack);

	/**
	Returning index of bin encompassing the position `pos`.

	@param pos Desired position.
	@param opt_snap_edge If true, position on the right/top-most edge of bounding box is considered in the last bin, and
	slightly before the left/bottom-most edge (may deal with due to lack of precision of double datatype) is snapped into the first bin.
	*/
	int* get_bin_from_pos(double* pos, bool opt_snap_edge=false);

	/**
	Returning a vector of cells that are sufficiently close to a bin.

	@param central_bin Index of bin.
	*/
	std::vector<int> get_cells_near_bin(
		int* central_bin,
		int half_width_search_window
	);

	/**
	Returning a vector of nodes that are sufficiently close to a bin.

	@param central_bin Index of bin.
	*/
	std::vector<int> get_nodes_near_bin(
		int* central_bin,
		int half_width_search_window
	);

	/**
	Returning a vector of cells that are sufficiently close to a given node.

	@param sn_node Index of node.
	@param spiral_level How many stacks of cells around the node to be considered. A spiral_level of 0 means that only cells with the given node as one of their vertices will be considered. By adding to `spiral_level`, more stacks of cells will be added to the vector of cells.
	*/
	std::vector<int> get_cells_near_node(int sn_node, int spiral_level);

	/**
	Returning a vector of cells that are sufficiently close to a given cell.

	@param sn_cell_central Index of target cell.
	@param spiral_level How many stacks of cells around the central cell to be considered. A spiral_level of 0 means that only the central cell will be considered. By adding to `spiral_level`, more stacks of cells will be added to the vector of cells.
	*/
	std::vector<int> get_cells_near_cell(int sn_cell_central, int spiral_level);

	/**
	Returns a vector of nodes from a vector of cells.

	@param cell_stack Vector of indices of cells of interest.
	*/
	std::vector<int> get_nodes_from_cells(std::vector<int> cell_stack);



	/**
	Returns a vector of nodes from a vector of cells with a constraint on zone of faces.

	@param cell_stack Vector of indices of cells of interest.
	@param zone_of_face A vector specifying one or multiple desired zones of face.
	*/
	std::vector<int> get_faces_from_cells(std::vector<int> cell_stack, std::vector<int> zone_of_face);

	/**
	Returns a vector of nodes from a vector of cells with a constraint on zone of faces.

	@param cell_stack Vector of indices of cells of interest.
	@param zone_of_face An integer specifying the desired zone of face.
	*/
	std::vector<int> get_faces_from_cells(std::vector<int> cell_stack, int zone_of_face);

	/**
	Returns a vector of faces from a vector of cells.

	@param cell_stack Vector of indices of cells of interest.
	*/
	std::vector<int> get_faces_from_cells(std::vector<int> cell_stack);



	/**
	Returns a vector of nodes related to faces and/or cells with a zone amongst a given vector of zones.
	Note: Nodes themselves are not directly assigned a zone, but faces and cells are.

	@param zone_of_node Vector of zones of interest.
	*/
	std::vector<int> get_nodes_of_zone_type(std::vector<int> zone_of_node);
	
	/**
	Returns a vector of faces with a zone amongst a given vector of zones.

	@param zone_of_face Vector of zones of interest.
	*/
	std::vector<int> get_faces_of_zone_type(std::vector<int> zone_of_face);

	/**
	Returns a vector of cells with a zone amongst a given vector of zones.

	@param zone_of_cell Vector of zones of interest.
	*/
	std::vector<int> get_cells_of_zone_type(std::vector<int> zone_of_cell);

	/**
	Find intersection of a face and a segment line connecting two points.

	@param sn_face Index of face of interest.
	@param seg_v1 Coordinates of first end of line segment.
	@param seg_v2 Coordinates of second end of line segment.
	*/
	std::tuple<bool, double *> face_intersection_with_line_segment(int sn_face, double *seg_v1, double *seg_v2);

	/**
	Find distance of a point from a face.

	@param sn_face Index of face of interest.
	@param point Coordinates of point of interest.
	*/
	double distance_point_from_face(int sn_face, double *point);

	/**
	Find mirror of a point w.r.t. a face.

	@param sn_face Index of face of interest.
	@param point Coordinates of point of interest.
	*/
	double *mirror_point_wrt_face(int sn_face, double *point);

	/**
	Evaluates the statistics of faces and stores the results in the `stat_face` struct.
	*/
	void eval_stats_faces();

	//-----------------------------------------------------------
	// getters for data attributes
	//-----------------------------------------------------------
	//get the solution field
	std::tuple<
		std::vector<double>, std::vector<double>,
		std::vector<double>, std::vector<double>,
		std::vector<double>, std::vector<double>,
		std::vector<double>> get_sol_full();

	//get only positions from the solution field
	std::tuple<std::vector<double>, std::vector<double>> get_sol_pos();

	//get only velocities from the solution field
	std::tuple<std::vector<double>, std::vector<double>> get_sol_vel();

	//get only pressure from the solution field
	std::vector<double> get_sol_pressure();
	void get_sol_pressure(double* output_array);

	//get only density from the solution field
	std::vector<double> get_sol_density();
	void get_sol_density(double* output_array);

	//get only viscosity from the solution field
	std::vector<double> get_sol_viscosity();
	void get_sol_viscosity(double* output_array);

	//get only x-component of position from the solution field
	std::vector<double> get_sol_pos_x();
	void get_sol_pos_x(double* output_array);

	//get only y-component of position from the solution field
	std::vector<double> get_sol_pos_y();
	void get_sol_pos_y(double* output_array);

	//get only x-component of velocity from the solution field
	std::vector<double> get_sol_vel_x();
	void get_sol_vel_x(double* output_array);

	//get only y-component of velocity from the solution field
	std::vector<double> get_sol_vel_y();
	void get_sol_vel_y(double* output_array);

	//getters
	/**
	Returns number of nodes in the grid.
	*/
	int get_num_nodes();

	//-----------------------------
	// data attributes
	//-----------------------------
	Cell *cell; /**<Cells available in the grid.*/
	Node *node; /**<Nodes available in the grid.*/
	Face *face; /**<Faces available in the grid.*/
	std::vector<Zone> zone;  /**<Zones extracted from the mesh file. Note: Ansys Fluent's zone index starts from 1, so a dummy item is created at index of 0.*/
	Bin **bin;				/**<Array of bins that forms the auxiliary structured grid (ASG) used when searching to locate hosting cell.*/

	int num_cells, /**<Number of cells of the grid.*/
		num_nodes, /**<Number of nodes of the grid.*/
		num_faces, /**<Number of faces of the grid.*/
		num_zones, /**<Number of zones of the grid.*/
		ndim, /**<Number of dimensions of the grid.*/
		*bin_num; /**<Number of bins of ASG in each direction.*/
	double **bounding_box, /**<Bounding box of grid used to form ASG.*/
		*bin_len; /**<Length of bins of ASG in each direction.*/
	Stat stat_face; /**<Statistics of faces of grid, e.g., min, max, and average.*/

	//------------------------------------------------
	// configs of search algorithm
	//------------------------------------------------
	bool opt_original_ASG;
	bool opt_bin_cell_check_center_only; 
	int half_width_search_window; 
	int encompassing_cell_spiral_level; 
	double nondim_criterion_area_ratio; 

	//----------------------------------
	// verbose
	//----------------------------------
	bool opt_verbose;
	bool opt_output;
};
