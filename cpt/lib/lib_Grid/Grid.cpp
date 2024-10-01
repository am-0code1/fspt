#include "Grid.h"
#include "Utils.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <map>
#include <tuple>
#include <algorithm>
#include <cstring>

#ifdef _WIN32
#include <ciso646> //for `and` and `or` keywords
#endif

#if defined(_WIN32)
#define Dir_Separator_Char "\\" //"/" also works on this system
#else
#define Dir_Separator_Char "/"
#endif

Grid::Grid(
	bool opt_original_ASG,
	bool opt_bin_cell_check_center_only,
	int half_width_search_window,
	int encompassing_cell_spiral_level,
	double nondim_criterion_area_ratio,
	bool opt_output,
	bool opt_verbose
	): 
	opt_original_ASG(opt_original_ASG),
	opt_bin_cell_check_center_only(opt_bin_cell_check_center_only),
	half_width_search_window(half_width_search_window),
	encompassing_cell_spiral_level(encompassing_cell_spiral_level),
	nondim_criterion_area_ratio(nondim_criterion_area_ratio),
	opt_output(opt_output),
	opt_verbose(opt_verbose)
{
	num_cells = 0;
	num_nodes = 0;
	num_faces = 0;
	num_zones = 0;
}

Grid::~Grid()
{

	delete[] cell;
	delete[] node;
	delete[] face;
	for (int ind_bin_x = 0; ind_bin_x < bin_num[0]; ind_bin_x++)
		delete[] bin[ind_bin_x];
	std::vector<Zone>().swap(zone);
	
	printf("Grid Destructed\n");
}

void Grid::calc_cell_props_from_nodes(
	bool pos,
	bool vel,
	bool density,
	bool pressure,
	bool viscosity)
{
	if (opt_verbose)
		std::cout << "Beginning of calc of cell props from nodes.\n";

	int sn_node;

	// index of sn_cell=0 refers to a dummy item
	for (int sn_cell = 1; sn_cell <= num_cells; sn_cell++)
	{
		// init
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			if (pos)
				cell[sn_cell].props.pos[ind_dir] = 0;
			if (vel)
				cell[sn_cell].props.vel[ind_dir] = 0;
		}
		if (density)
			cell[sn_cell].props.density = 0;
		if (pressure)
			cell[sn_cell].props.pressure = 0;
		if (viscosity)
			cell[sn_cell].props.viscosity = 0;

		// accumulation
		for (int ind_n = 0; ind_n < (int)cell[sn_cell].node.size(); ind_n++)
		{
			sn_node = cell[sn_cell].node[ind_n];

			if (density)
				cell[sn_cell].props.density += node[sn_node].props.density;
			if (pressure)
				cell[sn_cell].props.pressure += node[sn_node].props.pressure;
			if (viscosity)
				cell[sn_cell].props.viscosity += node[sn_node].props.viscosity;

			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			{
				if (pos)
					cell[sn_cell].props.pos[ind_dir] += node[sn_node].props.pos[ind_dir];

				// avg. props can also be calculated in the future based on the inverse of distance from nodes
				if (vel)
					cell[sn_cell].props.vel[ind_dir] += node[sn_node].props.vel[ind_dir];
			}
		}

		// averaging
		if (cell[sn_cell].node.size() > 0)
		{
			if (density)
				cell[sn_cell].props.density /= cell[sn_cell].node.size();

			if (pressure)
				cell[sn_cell].props.pressure /= cell[sn_cell].node.size();

			if (viscosity)
				cell[sn_cell].props.viscosity /= cell[sn_cell].node.size();

			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			{
				if (pos)
					cell[sn_cell].props.pos[ind_dir] /= cell[sn_cell].node.size();

				if (vel)
					cell[sn_cell].props.vel[ind_dir] /= cell[sn_cell].node.size();
			}
		}
	}

	if (opt_verbose)
		std::cout << "End of calc.\n\n";
}

void Grid::calc_cells_aux_props()
{
	if (opt_verbose)
		std::cout << "\nCalc. cell aux props...";
	
	double this_vel_mag;

	// index of sn_cell=0 refers to a dummy item
	for (int sn_cell = 1; sn_cell <= num_cells; sn_cell++)
	{
		this_vel_mag = 0;
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			this_vel_mag += std::pow(cell[sn_cell].props.vel[ind_dir], 2);
		cell[sn_cell].props_aux.vel_mag = std::sqrt(this_vel_mag);
	}
	if (opt_verbose)
		std::cout << "\nEnd of calc.\n\n";
}

double Grid::calc_cell_area(int sn_cell)
{
	double area = 0;
	int sn_node_1, sn_node_2, sn_node_3;

	for (int ind_tri = 0; ind_tri < (int)cell[sn_cell].node.size() - 2; ind_tri++)
	{
		// first vertex of triangle is fixed, while the other two shift until the whole polygon is covered.
		sn_node_1 = cell[sn_cell].node.at(0);
		sn_node_2 = cell[sn_cell].node.at(1 + ind_tri);
		sn_node_3 = cell[sn_cell].node.at(2 + ind_tri);

		area += area_triangle(
			node[sn_node_1].props.pos,
			node[sn_node_2].props.pos,
			node[sn_node_3].props.pos);
	}
	return area;
}

void Grid::calc_cells_area()
{
	if (opt_verbose)
		std::cout << "Beginning of calc of cells area.\n";
	
	// skipping the first dummy cell
	for (int sn_cell = 1; sn_cell <= num_cells; sn_cell++)
		cell[sn_cell].area = calc_cell_area(sn_cell);
	
	if (opt_verbose)
		std::cout << "End of calc.\n\n";
}

void Grid::set_up_bins(
	int *n_bins,
	const char *fname,
	bool opt_report_bins,
	double bin_len_dimensionless,
	const char *delimiter
)
{
	if (opt_verbose)
		std::cout << "Setting up bins ..." << std::endl;

	set_up_bounding_box();
	
	if (n_bins[0]<1 or n_bins[1] < 1)
		set_bin_num(bin_len_dimensionless * stat_face.max);
	else
		set_bin_num(n_bins);

	init_bins();
	populate_bins();
	if (opt_report_bins)
		report_bins(fname, delimiter);

	if (opt_verbose)
		std::cout << "End of setting up bins.\n\n";
}

void Grid::report_bins(const char *fname, const char *delimiter)
{
	using std::cout;
	using std::endl;

	if (opt_verbose)
		cout << "Preparing report bins ... " << endl;
	std::ofstream fout(fname);
	
	// #pragma region report bin.node
	//  reporting bin.node
	fout << "i" << delimiter
		 << "j" << delimiter
		 << "vertex.x" << delimiter
		 << "vertex.y" << delimiter
		 << "Nodes"
		 << endl;
	for (int i = 0; i < bin_num[0]; i++)
	{
		for (int j = 0; j < bin_num[1]; j++)
		{
			fout << i << delimiter
				 << j << delimiter
				 << bin[i][j].vertex[0] << delimiter
				 << bin[i][j].vertex[1] << delimiter;
			for (int ind_n = 0; ind_n < (int)bin[i][j].node.size(); ind_n++)
				fout << bin[i][j].node[ind_n] << delimiter;
			fout << endl;
		}
		fout << endl;
	}
	fout.close();
	// #pragma endregion

	// #pragma region report bin.cell
	fout.clear();
	fout.open(std::string(fname) + "_cells");

	fout << "i" << delimiter
		<< "j" << delimiter
		<< "vertex.x" << delimiter
		<< "vertex.y" << delimiter
		<< "Cells"
		<< endl;
	for (int i = 0; i < bin_num[0]; i++)
	{
		for (int j = 0; j < bin_num[1]; j++)
		{
			fout << i << delimiter
				<< j << delimiter
				<< bin[i][j].vertex[0] << delimiter
				<< bin[i][j].vertex[1] << delimiter;
			for (int ind_n = 0; ind_n < (int)bin[i][j].cell.size(); ind_n++)
				fout << bin[i][j].cell[ind_n] << delimiter;
			fout << endl;
		}
		fout << endl;
	}
	// #pragma endregion
}

void Grid::set_bin_num(double bin_len)
{
	bin_num = new int[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		bin_num[ind_dir] = (int)std::floor((bounding_box[1][ind_dir] - bounding_box[0][ind_dir]) / bin_len);
}

void Grid::set_bin_num(int *n_bins)
{
	bin_num = new int[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		bin_num[ind_dir] = n_bins[ind_dir]; // default value for number of bins in each direction
}

void Grid::set_up_bounding_box()
{
	// bbox: [lower_left[x_1, x_2,..., x_ndim], uper_right[x_1, x_2,..., x_ndim]]
	bounding_box = new double *[2];
	for (int i = 0; i < 2; i++)
		bounding_box[i] = new double[ndim];
	for (int i = 0; i < ndim; i++)
	{
		bounding_box[0][i] = node[1].props.pos[i];
		bounding_box[1][i] = node[1].props.pos[i];
	}

	// node index start from 1; index of 0 is dummy to match the 1-based style of ansys fluent mesh
	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
	{
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			if (node[sn_node].props.pos[ind_dir] > bounding_box[1][ind_dir])
				bounding_box[1][ind_dir] = node[sn_node].props.pos[ind_dir];

			if (node[sn_node].props.pos[ind_dir] < bounding_box[0][ind_dir])
				bounding_box[0][ind_dir] = node[sn_node].props.pos[ind_dir];
		}
	}
}

void Grid::init_bins()
{
	// bins
	bin_len = new double[ndim];
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		bin_len[ind_dir] = (bounding_box[1][ind_dir] - bounding_box[0][ind_dir]) / bin_num[ind_dir];
	}

	bin = new Bin *[bin_num[0]];
	for (int i = 0; i < bin_num[0]; i++)
		bin[i] = new Bin[bin_num[1]];

	// can be extended to 3D
	for (int i = 0; i < bin_num[0]; i++)
		for (int j = 0; j < bin_num[1]; j++)
			bin[i][j].vertex = new double[ndim];

	for (int i = 0; i < bin_num[0]; i++)
		for (int j = 0; j < bin_num[1]; j++)
		{
			bin[i][j].vertex[0] = bounding_box[0][0] + i * bin_len[0];
			bin[i][j].vertex[1] = bounding_box[0][1] + j * bin_len[1];
		}
}

void Grid::populate_bins()
{
	int this_bin[2];

	// important to not start from sn_node=0 as it is a dummy item.
	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
	{
		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		{
			this_bin[ind_dir] = (int)std::floor((node[sn_node].props.pos[ind_dir] - bounding_box[0][ind_dir]) / bin_len[ind_dir]);
			// Nodes residing right on the rightmost edge to be included in the last bin
			if (this_bin[ind_dir] == bin_num[ind_dir])
				this_bin[ind_dir] = bin_num[ind_dir] - 1;
		}

		// A sanity check to ensure a node is included in a bin only once at max. Seems unnecessary for this implementation.
		if (std::find(
				bin[this_bin[0]][this_bin[1]].node.begin(),
				bin[this_bin[0]][this_bin[1]].node.end(), sn_node) ==
			bin[this_bin[0]][this_bin[1]].node.end())
		{
			bin[this_bin[0]][this_bin[1]].node.push_back(sn_node);
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				node[sn_node].bin[ind_dir] = this_bin[ind_dir];
		}
	}

	//declarations
	int sn_node;
	
	if (opt_bin_cell_check_center_only)
	{
		//important to not start from sn_cell=0 as it is a dummy item.
		for (int sn_cell = 1; sn_cell <= num_cells; sn_cell++)
		{
			//adding cells based on their center position
			for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			{
				this_bin[ind_dir] = (int)std::floor((cell[sn_cell].props.pos[ind_dir] - bounding_box[0][ind_dir]) / bin_len[ind_dir]);
				// Cells with a center residing right on the rightmost edge to be included in the last bin
				if (this_bin[ind_dir] == bin_num[ind_dir])
					this_bin[ind_dir] = bin_num[ind_dir] - 1;
			}

			// A sanity check to ensure a cell is included in a bin only once at max. Seems unnecessary for this implementation.
			if (std::find(
				bin[this_bin[0]][this_bin[1]].cell.begin(),
				bin[this_bin[0]][this_bin[1]].cell.end(), sn_cell) ==
				bin[this_bin[0]][this_bin[1]].cell.end())
			{
				bin[this_bin[0]][this_bin[1]].cell.push_back(sn_cell);
				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
					cell[sn_cell].bin[ind_dir] = this_bin[ind_dir];
			}
		}
	}
	else
	{
		//important to not start from sn_cell=0 as it is a dummy item.
		for (int sn_cell = 1; sn_cell <= num_cells; sn_cell++)
		{
			for (int ind_node = 0; ind_node < cell[sn_cell].node.size(); ind_node++)
			{
				sn_node = cell[sn_cell].node[ind_node];

				for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
				{
					this_bin[ind_dir] = (int)std::floor((node[sn_node].props.pos[ind_dir] - bounding_box[0][ind_dir]) / bin_len[ind_dir]);
					// Nodes residing right on the rightmost edge to be included in the last bin
					if (this_bin[ind_dir] == bin_num[ind_dir])
						this_bin[ind_dir] = bin_num[ind_dir] - 1;
				}

				// A sanity check to ensure a cell is included in a bin only once at max. It is necessary for this implementation as a cell may have
				// already been assigned to a bin due to its other nodes.
				if (std::find(
					bin[this_bin[0]][this_bin[1]].cell.begin(),
					bin[this_bin[0]][this_bin[1]].cell.end(), sn_cell) ==
					bin[this_bin[0]][this_bin[1]].cell.end())
				{
					bin[this_bin[0]][this_bin[1]].cell.push_back(sn_cell);

					//What about the `bin` attr of `cell`?
					//Do we need to keep track of all bins that all nodes of cell reside in?
					//currently, only one bin is stored and that is determined by the last node.
					for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
						cell[sn_cell].bin[ind_dir] = this_bin[ind_dir];
				}
			}
		}
	}
}

std::tuple<int, double, int> Grid::find_node_close_to_pos_spiral_algorithm(
	double *pos,
	double allowed_tolerance_node_position_snap,
	int flag_normal_return,
	int flag_far_to_snap,
	int flag_out_of_domain_or_small_search_win,
	int flag_out_of_bbox)
{
	int central_bin[2], sn_node, sn_node_closest, count_examined_nodes = 0;
	double dist_min, dist;
	double EPSILON = 1e-12;

	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		central_bin[ind_dir] = (int)std::floor((pos[ind_dir] - bounding_box[0][ind_dir]) / bin_len[ind_dir]);

		// positions and/or nodes that are right on the rightmost edge or slightly beyond that due to loss of precision of double datatype
		if (central_bin[ind_dir] == bin_num[ind_dir])
			central_bin[ind_dir] = bin_num[ind_dir] - 1;

		// positions and/or nodes that are slightly befor the leftmost edge due to loss of precision of double datatype
		if (central_bin[ind_dir] == -1)
		{
			if (opt_verbose)
				std::cout << "Warning: A position in bin[0]=-1 is assigned bin[0]=0\n";
			central_bin[ind_dir] = 0;
		}
	}

	// #pragma region Out of bounding box
	//  position out of bounding box
	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
		if (central_bin[ind_dir] < 0 or central_bin[ind_dir] >= bin_num[ind_dir])
		{
			/*
			 * Position out of bounding box.
			 * Returning invalid values as flag to be handled by the caller.
			 * This flag cannot be tolerated when reading the solution file where
			 * the props need to be attributed to nodes at exact positions.
			 * On the other hand, the same flag can be an indication of particles crossing
			 * a boundary when tracking particles.
			 */
			sn_node_closest = flag_out_of_bbox;
			dist_min = -2;
			return std::make_tuple(sn_node_closest, dist_min, flag_out_of_bbox);
		}
	// #pragma endregion

	dist_min = 1 / EPSILON; // initialize with a large number; larger than minimum spacing between nodes

	// #pragma region Central bin
	//  central bin
	for (int sn_item = 0; sn_item < (int)bin[central_bin[0]][central_bin[1]].node.size(); sn_item++)
	{
		count_examined_nodes++;
		sn_node = bin[central_bin[0]][central_bin[1]].node.at(sn_item);

		dist = distance_between_two_points(pos, node[sn_node].props.pos, ndim);

		if (dist < dist_min)
		{
			dist_min = dist;
			sn_node_closest = sn_node;
		}

		// assuming the minimum spacing between nodes is greater than epsilon
		if (dist_min < EPSILON)
		{
			return std::make_tuple(sn_node_closest, dist_min, flag_normal_return);
		}
	}
	// #pragma endregion

	// spiral out
	for (int dist_bin_from_central_bin = 1; dist_bin_from_central_bin <= half_width_search_window; dist_bin_from_central_bin++)
	{
		// #pragma region left and right edges
		//  left and right edges
		for (int this_bin_x = central_bin[0] - dist_bin_from_central_bin;
			 this_bin_x <= central_bin[0] + dist_bin_from_central_bin;
			 this_bin_x += 2 * dist_bin_from_central_bin)
		{
			if (this_bin_x < 0 or this_bin_x >= bin_num[0])
				continue;

			for (int this_bin_y = central_bin[1] - dist_bin_from_central_bin;
				 this_bin_y <= central_bin[1] + dist_bin_from_central_bin;
				 this_bin_y++)
			{
				if (this_bin_y < 0 or this_bin_y >= bin_num[1])
					continue;

				for (int sn_item = 0; sn_item < (int)bin[this_bin_x][this_bin_y].node.size(); sn_item++)
				{
					count_examined_nodes++;
					sn_node = bin[this_bin_x][this_bin_y].node.at(sn_item);
					dist = distance_between_two_points(pos, node[sn_node].props.pos, ndim);

					if (dist < dist_min)
					{
						dist_min = dist;
						sn_node_closest = sn_node;
					}

					// assuming the minimum spacing between nodes is greater than epsilon
					if (dist_min < EPSILON)
					{
						return std::make_tuple(sn_node_closest, dist_min, flag_normal_return);
					}
				}
			}
		}
		// #pragma endregion

		// #pragma region bottom and top edges
		//  bottom and top edges
		for (int this_bin_y = central_bin[1] - dist_bin_from_central_bin;
			 this_bin_y <= central_bin[1] + dist_bin_from_central_bin;
			 this_bin_y += 2 * dist_bin_from_central_bin)
		{
			if (this_bin_y < 0 or this_bin_y >= bin_num[1])
				continue;

			for (int this_bin_x = central_bin[0] - dist_bin_from_central_bin + 1;
				 this_bin_x <= central_bin[0] + dist_bin_from_central_bin - 1;
				 this_bin_x++)
			{
				if (this_bin_x < 0 or this_bin_x >= bin_num[0])
					continue;

				for (int sn_item = 0; sn_item < (int)bin[this_bin_x][this_bin_y].node.size(); sn_item++)
				{
					count_examined_nodes++;
					sn_node = bin[this_bin_x][this_bin_y].node.at(sn_item);
					dist = distance_between_two_points(pos, node[sn_node].props.pos, ndim);

					if (dist < dist_min)
					{
						dist_min = dist;
						sn_node_closest = sn_node;
					}

					// assuming the minimum spacing between nodes is greater than epsilon
					if (dist_min < EPSILON)
					{
						return std::make_tuple(sn_node_closest, dist_min, flag_normal_return);
					}
				}
			}
		}
	}
	// #pragma endregion

	if (count_examined_nodes == 0)
	{
		if (opt_verbose)
			std::cout << "No nodes close to the position: (" << pos[0] << ", " << pos[1] << ").\n";
		/*
		 * Position within the bounding box, but at least one of the following may be the case:
		 *
		 * 1. `half_width_search_window` is not sufficiently large to cover an area with a node
		 *
		 * 2. poisition may be out of mesh, e.g. a hole within the domain.
		 * Returning invalid values as flag to be handled based on application,
		 * e.g. this flag cannot be tolerated when reading the solution file where
		 * the props need to be attributed to nodes at exact positions.
		 * On the other hand, the same flag can be an indication of particles crossing
		 * a boundary when tracking particles.
		 */
		sn_node_closest = flag_out_of_domain_or_small_search_win;
		dist_min = -1;
		return std::make_tuple(sn_node_closest, dist_min, flag_out_of_domain_or_small_search_win);
	}
	else if (dist_min > allowed_tolerance_node_position_snap)
	{
		if (opt_verbose)
		{
			std::cout << "\nError" << std::endl;
			std::cout
				<< "Provided coordinates: ("
				<< pos[0] << ", "
				<< pos[1] << ")\n"
				<< "Closest coordinates found within the mesh: ("
				<< node[sn_node_closest].props.pos[0] << ", "
				<< node[sn_node_closest].props.pos[1] << ")\n"
				<< "Difference between points: " << dist_min << "\n"
				<< "Make sure the provided `allowed_tolerance_node_position_snap` is not too restrictive.\n"
				<< "`allowed_tolerance_node_position_snap` currently is: "
				<< allowed_tolerance_node_position_snap << "\n";
		}
		return std::make_tuple(sn_node_closest, dist_min, flag_far_to_snap);
	}
	else
		return std::make_tuple(sn_node_closest, dist_min, flag_normal_return);
}


std::tuple<int, double, int> Grid::find_node_close_to_pos(
	double* pos,
	double allowed_tolerance_node_position_snap,
	int flag_normal_return,
	int flag_far_to_snap,
	int flag_out_of_domain_or_small_search_win,
	int flag_out_of_bbox,
	const Event_Execution_Time& eet
	)
{
	//-----------------------------------------
	// declarations
	//-----------------------------------------
	int *central_bin, sn_node_closest, index_min_elem;
	double dist_min;
	std::vector<int> node_stack;
	std::vector<double> dist_stack;
	//-----------------------------------------

	// event: start
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).start_event("find_node_close_to_pos | prep");

	central_bin = get_bin_from_pos(pos, true);
	node_stack=get_nodes_near_bin(central_bin, half_width_search_window);
	
	if (node_stack.size() == 0)
	{
		if (opt_verbose)
			std::cout << "No nodes close to the position: (" << pos[0] << ", " << pos[1] << ").\n";
		/*
		 * Position within the bounding box, but at least one of the following may be the case:
		 *
		 * 1. `half_width_search_window` is not sufficiently large to cover an area with a node
		 *
		 * 2. poisition may be out of mesh, e.g. a hole within the domain.
		 * Returning invalid values as flag to be handled based on application,
		 * e.g. this flag cannot be tolerated when reading the solution file where
		 * the props need to be attributed to nodes at exact positions.
		 * On the other hand, the same flag can be an indication of particles crossing
		 * a boundary when tracking particles.
		 */
		sn_node_closest = flag_out_of_domain_or_small_search_win;
		dist_min = -1;
		return std::make_tuple(sn_node_closest, dist_min, flag_out_of_domain_or_small_search_win);
	}


	//dynamic init of vector
	dist_stack.reserve(node_stack.size());
	dist_stack.insert(dist_stack.begin(), node_stack.size(), -1);

	// event: stop
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).stop_event("find_node_close_to_pos | prep");
	
	// event: start
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).start_event("find_node_close_to_pos | vectorized");

#pragma omp parallel for 
	for (int ind_node = 0; ind_node < node_stack.size(); ind_node++)
	{		
		dist_stack[ind_node]=distance_between_two_points(pos, node[node_stack[ind_node]].props.pos, ndim);
	}

	// event: stop
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).stop_event("find_node_close_to_pos | vectorized");

	// event: start
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).start_event("find_node_close_to_pos | find_min");

	auto it = std::min_element(dist_stack.begin(), dist_stack.end());
	index_min_elem = std::distance(dist_stack.begin(), it);
	sn_node_closest = node_stack.at(index_min_elem);
	dist_min = *it;

	// event: stop
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).stop_event("find_node_close_to_pos | find_min");
	

	if (dist_min > allowed_tolerance_node_position_snap)
	{
		if (opt_verbose)
		{
			std::cout << "\nError" << std::endl;
			std::cout
				<< "Provided coordinates: ("
				<< pos[0] << ", "
				<< pos[1] << ")\n"
				<< "Closest coordinates found within the mesh: ("
				<< node[sn_node_closest].props.pos[0] << ", "
				<< node[sn_node_closest].props.pos[1] << ")\n"
				<< "Difference between points: " << dist_min << "\n"
				<< "Make sure the provided `allowed_tolerance_node_position_snap` is not too restrictive.\n"
				<< "`allowed_tolerance_node_position_snap` currently is: "
				<< allowed_tolerance_node_position_snap << "\n";
		}
		return std::make_tuple(sn_node_closest, dist_min, flag_far_to_snap);
	}
	else
		return std::make_tuple(sn_node_closest, dist_min, flag_normal_return);
}



std::tuple<int, std::vector<double>, int, std::vector<int>> Grid::find_cell_encompassing_pos(
	double *pos,
	int flag_normal_return,
	int flag_not_finding_encomp_cell,
	int flag_not_finding_close_node,
	const Event_Execution_Time& eet
)
{
	//---------------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------------
	double dist, EPSILON = 1e-12;
	int sn_node_closest=-1, sn_cell = -1, flag_finding_close_node;
	std::vector<double> dist_from_nodes;
	std::vector<int> cell_search_stack;
	int* central_bin;
	//---------------------------------------------------------------------

	// event: start
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).start_event("find_cell_encompassing_pos | find_node_close_to_pos_spiral_algorithm");

	int flag_far_to_snap = flag_not_finding_close_node - 1;
	int finding_close_node_normal_return = flag_not_finding_close_node + 1;

	if (!opt_original_ASG)
	{
		//std::tie(sn_node_closest, dist, flag_finding_close_node) = find_node_close_to_pos_spiral_algorithm(
		std::tie(sn_node_closest, dist, flag_finding_close_node) = find_node_close_to_pos(
			pos,
			1 / EPSILON,
			finding_close_node_normal_return,
			flag_far_to_snap,
			flag_not_finding_close_node,
			flag_not_finding_close_node
			//eet
			);

		// #pragma region Handling exceptions of not finding a valid close node
		if (flag_finding_close_node == flag_not_finding_close_node)
		{
			sn_cell = flag_not_finding_close_node;
			dist_from_nodes.push_back(-1);
			return std::make_tuple(sn_cell, dist_from_nodes, flag_not_finding_close_node, cell_search_stack);
		}
		else if (flag_finding_close_node == flag_far_to_snap)
		{
			if (opt_verbose)
			{
				std::cout << "\nError in `find_cell_encompassing_pos` after calling `find_node_close_to_pos_spiral_algorithm`\n";
				std::cout << "Close node found but is too far. The distance is " << dist << std::endl;
			}
			exit(1);
		}
		// #pragma endregion
	}

	// event: stop
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).stop_event("find_cell_encompassing_pos | find_node_close_to_pos_spiral_algorithm");


	// event: start
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).start_event("find_cell_encompassing_pos | get_cells_near_node");

	if (!opt_original_ASG)
	{
		// prep search stack of cells
		cell_search_stack = get_cells_near_node(sn_node_closest, encompassing_cell_spiral_level);
	}
	else
	{
		central_bin = get_bin_from_pos(pos, true);

		for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
			if (central_bin[ind_dir] < 0 or central_bin[ind_dir] >= bin_num[ind_dir])
				return std::make_tuple(sn_cell, dist_from_nodes, flag_not_finding_encomp_cell, cell_search_stack);
				
		cell_search_stack = get_cells_near_bin(central_bin, half_width_search_window);
	}

	// event: stop
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).stop_event("find_cell_encompassing_pos | get_cells_near_node");



	// event: start
	if (eet.active)
		const_cast<Event_Execution_Time&>(eet).start_event("find_cell_encompassing_pos | loop over near cells");

	// #pragma region iterating over the search stack looking for the encompassing cell
	//  when an encompassing cell is found, the loop is broken out of, so sn_cell
	//  holds the sn of encompassing cell.
	for (int ind_c = 0; ind_c < (int)cell_search_stack.size(); ind_c++)
	{
		sn_cell = cell_search_stack[ind_c];

		// sn_cell=0 is a dummy cell related to boundary faces.
		if (sn_cell == 0)
			continue;
		else if (pos_inside_cell(pos, sn_cell, nondim_criterion_area_ratio))
		{
			// event: stop
			if (eet.active)
				const_cast<Event_Execution_Time&>(eet).stop_event("find_cell_encompassing_pos | loop over near cells");
			break;
		}
		else if (ind_c == (int)cell_search_stack.size() - 1)
		{
			/*
			 * Not finding an encompassing cell while a vlid close node has been found.
			 * Potential scenarios:
			 *
			 * 1. spiral level may be small. There are cases with skewed cells that spiral
			 * level of 0 fail as the encompassing cell may not be amongst the cells of the
			 * found closest node. For a normal mesh, a spiral level of 1 should be sufficiently
			 * large to deal with those cases. However, for grids with stack of high aspect ratio
			 * cells and/or skewed cells, a larger level may be needed. For example, for a grid with
			 * `n` inflation layers, the spiral level should be at least `n`.
			 *
			 * 2. Position may be out of domain, e.g. a hole within the mesh that the position falls
			 * within.
			 *
			 * 3. criterion for a position to be inside a cell is based on comparing the cell area and
			 * the sum area created by position and every pair of nodes of cell. The criterion might
			 * have been too tight. A threshold of 0.01 should be a good starting point.
			 */

			sn_cell = flag_not_finding_encomp_cell;

			dist_from_nodes.push_back(sn_node_closest); // sn_node_closest is returned to be used by caller as needed
			
			// event: stop
			if (eet.active)
				const_cast<Event_Execution_Time&>(eet).stop_event("find_cell_encompassing_pos | loop over near cells");
			
			return std::make_tuple(sn_cell, dist_from_nodes, flag_not_finding_encomp_cell, cell_search_stack);
		}
	}
	// #pragma endregion

	// prep distance vector that can be used as a weight factor when interpolating props of a loc within cell
	dist_from_nodes = get_dist_from_nodes(pos, cell[sn_cell].node);

	return std::make_tuple(sn_cell, dist_from_nodes, flag_normal_return, cell_search_stack);
}



std::vector<double> Grid::get_dist_from_nodes(double *pos, std::vector<int> node_stack)
{
	// Declarations
	//---------------------------------
	int sn_node;
	std::vector<double> dist_from_nodes;
	//---------------------------------

	for (int ind_n = 0; ind_n < (int)node_stack.size(); ind_n++)
	{
		sn_node = node_stack.at(ind_n);
		dist_from_nodes.push_back(distance_between_two_points(pos, node[sn_node].props.pos, ndim));
	}

	return dist_from_nodes;
}

std::vector<int> Grid::get_cells_near_cell(int sn_cell_central, int spiral_level)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------
	std::vector<int> cell_stack, cell_stack_to_add;
	int sn_cell, cell_spiral, sn_node;
	//---------------------------------------------------------------

	// include central cell -- skip sn_cell=0 as it is a dummy item
	if (sn_cell_central != 0)
		cell_stack.push_back(sn_cell_central);

	int ind_c_start = 0;
	for (int ind_spiral_level = 0; ind_spiral_level < spiral_level; ind_spiral_level++)
	{
		cell_stack_to_add.clear();
		for (int ind_c = ind_c_start; ind_c < (int)cell_stack.size(); ind_c++)
		{
			sn_cell = cell_stack[ind_c];
			for (int ind_n = 0; ind_n < (int)cell[sn_cell].node.size(); ind_n++)
			{
				sn_node = cell[sn_cell].node[ind_n];
				for (int ind_c_out = 0; ind_c_out < (int)node[sn_node].cell.size(); ind_c_out++)
				{
					cell_spiral = node[sn_node].cell[ind_c_out];

					// skip sn_cell=0 as it is a dummy item
					if (cell_spiral == 0)
						continue;

					// add the cell to the search stack if not already in either stack
					if (std::find(
							cell_stack.begin(),
							cell_stack.end(), cell_spiral) == cell_stack.end() and
						std::find(
							cell_stack_to_add.begin(),
							cell_stack_to_add.end(), cell_spiral) == cell_stack_to_add.end())
					{
						cell_stack_to_add.push_back(cell_spiral);
					}
				}
			}
		}

		// setting start index of search stack for the next spiral level
		// Note: it needs to be done before adding found cells to stack.
		ind_c_start = (int)cell_stack.size();

		// adding found cells to stack
		for (int ind_c = 0; ind_c < (int)cell_stack_to_add.size(); ind_c++)
			cell_stack.push_back(cell_stack_to_add[ind_c]);
	}

	return cell_stack;
}

int* Grid::get_bin_from_pos(double* pos, bool opt_snap_edge)
{
	int* central_bin;
	central_bin = new int[ndim];

	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		central_bin[ind_dir] = (int)std::floor((pos[ind_dir] - bounding_box[0][ind_dir]) / bin_len[ind_dir]);

		if (opt_snap_edge)
		{
			// positions and/or nodes that are right on the rightmost edge or slightly beyond that due to loss of precision of double datatype
			if (central_bin[ind_dir] == bin_num[ind_dir])
				central_bin[ind_dir] = bin_num[ind_dir] - 1;

			// positions and/or nodes that are slightly befor the leftmost edge due to loss of precision of double datatype
			if (central_bin[ind_dir] == -1)
			{
				if (opt_verbose)
					std::cout << "Warning: A position in bin[0]=-1 is assigned bin[0]=0\n";
				central_bin[ind_dir] = 0;
			}
		}
	}

	return central_bin;
}

std::vector<int> Grid::get_cells_near_bin(
	int* central_bin,
	int half_width_search_window
)
{
	//declarations
	std::vector<int> cell_stack;
	int sn_cell;
	//--------------------------------------

	for (int this_bin_x = central_bin[0] - half_width_search_window; this_bin_x <= central_bin[0] + half_width_search_window; this_bin_x++)
	{
		if (this_bin_x < 0 or this_bin_x >= bin_num[0])
			continue;
		for (int this_bin_y = central_bin[1] - half_width_search_window; this_bin_y <= central_bin[1] + half_width_search_window; this_bin_y++)
		{
			if (this_bin_y < 0 or this_bin_y >= bin_num[1])
				continue;
			for (int ind_cell = 0; ind_cell < (int)bin[this_bin_x][this_bin_y].cell.size(); ind_cell++)
			{
				sn_cell = bin[this_bin_x][this_bin_y].cell[ind_cell];

				//a check to avoid inclusion of items more than once -- not necessary.
				/*
				if (std::find(
					cell_stack.begin(),
					cell_stack.end(), sn_cell) ==
					cell_stack.end())
				{
					cell_stack.push_back(sn_cell);
				}
				*/

				cell_stack.push_back(sn_cell);
			}
		}
	}
	
	return cell_stack;
}

std::vector<int> Grid::get_nodes_near_bin(
	int* central_bin,
	int half_width_search_window
)
{
	//declarations
	std::vector<int> node_stack;
	int sn_node;
	//--------------------------------------

	for (int this_bin_x = central_bin[0] - half_width_search_window; this_bin_x <= central_bin[0] + half_width_search_window; this_bin_x++)
	{
		if (this_bin_x < 0 or this_bin_x >= bin_num[0])
			continue;
		for (int this_bin_y = central_bin[1] - half_width_search_window; this_bin_y <= central_bin[1] + half_width_search_window; this_bin_y++)
		{
			if (this_bin_y < 0 or this_bin_y >= bin_num[1])
				continue;
			for (int ind_node = 0; ind_node < (int)bin[this_bin_x][this_bin_y].node.size(); ind_node++)
			{
				sn_node = bin[this_bin_x][this_bin_y].node[ind_node];
				node_stack.push_back(sn_node);
			}
		}
	}

	return node_stack;
}


std::vector<int> Grid::get_cells_near_node(
	int sn_node, 
	int spiral_level 
)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------
	std::vector<int> cell_search_stack, cell_search_stack_to_add;
	int sn_cell, cell_spiral;
	//---------------------------------------------------------------

	for (int ind_c = 0; ind_c < (int)node[sn_node].cell.size(); ind_c++)
	{
		sn_cell = node[sn_node].cell[ind_c];
		// skip sn_cell=0 as it is a dummy item
		if (sn_cell == 0)
			continue;
		cell_search_stack.push_back(sn_cell);
	}

	int ind_c_start = 0;
	for (int ind_spiral_level = 0; ind_spiral_level < spiral_level; ind_spiral_level++)
	{
		cell_search_stack_to_add.clear();
		for (int ind_c = ind_c_start; ind_c < (int)cell_search_stack.size(); ind_c++)
		{
			sn_cell = cell_search_stack[ind_c];
			for (int ind_n = 0; ind_n < (int)cell[sn_cell].node.size(); ind_n++)
			{
				sn_node = cell[sn_cell].node[ind_n];
				for (int ind_c_out = 0; ind_c_out < (int)node[sn_node].cell.size(); ind_c_out++)
				{
					cell_spiral = node[sn_node].cell[ind_c_out];

					// skip sn_cell=0 as it is a dummy item
					if (cell_spiral == 0)
						continue;

					// add the cell to the search stack if not already in either stack
					if (std::find(
							cell_search_stack.begin(),
							cell_search_stack.end(), cell_spiral) == cell_search_stack.end() and
						std::find(
							cell_search_stack_to_add.begin(),
							cell_search_stack_to_add.end(), cell_spiral) == cell_search_stack_to_add.end())
					{
						cell_search_stack_to_add.push_back(cell_spiral);
					}
				}
			}
		}

		// setting start index of search stack for the next spiral level
		// Note: it needs to be done before adding found cells to stack.
		ind_c_start = (int)cell_search_stack.size();

		// adding found cells to search stack
		for (int ind_c = 0; ind_c < (int)cell_search_stack_to_add.size(); ind_c++)
			// add the cell to the search stack if not already in the stack
			// the following check is not necessary as is already implemented at the time of preparing
			// cell_search_stack_to_add
			// The redundant check is removed on Sept. 15, 2024
			/*
			if (std::find(
					cell_search_stack.begin(),
					cell_search_stack.end(), cell_search_stack_to_add[ind_c]) == cell_search_stack.end())
			{
				cell_search_stack.push_back(cell_search_stack_to_add[ind_c]);
			}
			*/
			cell_search_stack.push_back(cell_search_stack_to_add[ind_c]);
	}

	return cell_search_stack;
}

std::vector<int> Grid::get_nodes_from_cells(std::vector<int> cell_stack)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------
	int sn_cell, sn_node;
	std::vector<int> node_stack;
	//---------------------------------------------------------------

	for (int ind_c = 0; ind_c < (int)cell_stack.size(); ind_c++)
	{
		sn_cell = cell_stack[ind_c];
		// skip sn_cell=0 as it is a dummy item
		if (sn_cell == 0)
			continue;

		for (int ind_n = 0; ind_n < (int)cell[sn_cell].node.size(); ind_n++)
		{
			sn_node = cell[sn_cell].node[ind_n];
			if (std::find(
					node_stack.begin(),
					node_stack.end(), sn_node) == node_stack.end())
				node_stack.push_back(sn_node);
		}
	}

	return node_stack;
}

std::vector<int> Grid::get_faces_from_cells(std::vector<int> cell_stack)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------
	int sn_cell, sn_face;
	std::vector<int> face_stack;
	//---------------------------------------------------------------

	for (int ind_c = 0; ind_c < (int)cell_stack.size(); ind_c++)
	{
		sn_cell = cell_stack[ind_c];
		// skip sn_cell=0 as it is a dummy item
		if (sn_cell == 0)
			continue;

		for (int ind_f = 0; ind_f < (int)cell[sn_cell].face.size(); ind_f++)
		{
			sn_face = cell[sn_cell].face[ind_f];
			if (std::find(
					face_stack.begin(),
					face_stack.end(), sn_face) == face_stack.end())
				face_stack.push_back(sn_face);
		}
	}

	return face_stack;
}

std::vector<int> Grid::get_faces_from_cells(std::vector<int> cell_stack, int zone_of_face)
{
	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------
	int sn_face;
	std::vector<int> face_stack_full, face_stack;
	//---------------------------------------------------------------

	face_stack_full = get_faces_from_cells(cell_stack);

	for (int ind_f = 0; ind_f < (int)face_stack_full.size(); ind_f++)
	{
		sn_face = face_stack_full[ind_f];

		// if (face_bc_type == face[sn_face].bc_type)
		if (std::find(
				face[sn_face].zone.begin(),
				face[sn_face].zone.end(), zone_of_face) != face[sn_face].zone.end())
			face_stack.push_back(sn_face);
	}

	return face_stack;
}

std::vector<int> Grid::get_faces_from_cells(std::vector<int> cell_stack, std::vector<int> zone_of_face)
{
	/*
	* returns a stck of faces with a zone type existing in the passed array of `zones_of_face`
	* 
	* getting a vector of zone_of_face to be able to search for an element inside the vector
	* pointer array seems to not be a vlid type of argument for `begin`. 
	* simple array is valid, but its length would need to be specified which is restricting.
	*/

	//---------------------------------------------------------------
	// Declarations
	//---------------------------------------------------------------
	int sn_face, sn_zone;
	std::vector<int> face_stack_full, face_stack;

	//---------------------------------------------------------------

	face_stack_full = get_faces_from_cells(cell_stack);

	for (int ind_f = 0; ind_f < (int)face_stack_full.size(); ind_f++)
	{
		sn_face = face_stack_full[ind_f];

		for (int ind_zone = 0; ind_zone < (int)face[sn_face].zone.size(); ind_zone++)
		{
			sn_zone = face[sn_face].zone[ind_zone];
			if (std::find(
				std::begin(zone_of_face),
				std::end(zone_of_face), sn_zone) != std::end(zone_of_face))
			{
				face_stack.push_back(sn_face);
				break; //to avoid potential duplications in face_stack if a face has multiple zonetypes existing in `zone_of_face`
			}
		}
	}

	return face_stack;
}

std::vector<int> Grid::get_nodes_of_zone_type(std::vector<int> zone_of_node)
{
	int sn_zone, sn_face, sn_cell, this_sn_node;
	std::vector<int> node_stack, face_stack, cell_stack;

	for (int sn_node = 0; sn_node <= num_nodes; sn_node++)
	{
		for (int ind_zone = 0; ind_zone < (int)node[sn_node].zone.size(); ind_zone++)
		{
			sn_zone = node[sn_node].zone[ind_zone];

			if (std::find(
				std::begin(zone_of_node),
				std::end(zone_of_node), sn_zone) != std::end(zone_of_node))
			{
				node_stack.push_back(sn_node);
				break; //to avoid potential duplications in node_stack if a node has multiple zonetypes existing in `zone_of_node`
			}
		}
	}

	// The nodes of faces and cells of the given zone type need to be considered
	// There are not overlap. For example, zone=2 may be assigned to a group of cells, but no face or node is assigned the same group.
	// Similarly, a group of faces may be assigned zone of 5, but no cells or nodes are assigned the same group.
	// Or, a group of nodes is assigned a zone id of 3, but no cells or faces have the same zone id.
	// As a result, the search for existance of nodes belonging to `face_stack` and `cell_stack` in `node_stack`
	// is not necessary, but is kept in the following. It can be removed later after more examinations. 
	face_stack = get_faces_of_zone_type(zone_of_node);

	for (int ind_face = 0; ind_face < (int)face_stack.size(); ind_face++)
	{
		sn_face = face_stack[ind_face];
		for (int ind_node = 0; ind_node < (int)face[sn_face].node.size(); ind_node++)
		{
			this_sn_node = face[sn_face].node[ind_node];
			if (std::find(
				std::begin(node_stack),
				std::end(node_stack), this_sn_node) == std::end(node_stack))
				node_stack.push_back(this_sn_node);
		}
	}
	
	
	cell_stack = get_cells_of_zone_type(zone_of_node);
	for (int ind_cell = 0; ind_cell < (int)cell_stack.size(); ind_cell++)
	{
		sn_cell = cell_stack[ind_cell];
		for (int ind_node = 0; ind_node < (int)cell[sn_cell].node.size(); ind_node++)
		{
			this_sn_node = cell[sn_cell].node[ind_node];
			if (std::find(
				std::begin(node_stack),
				std::end(node_stack), this_sn_node) == std::end(node_stack))
				node_stack.push_back(this_sn_node);
		}
	}

	return node_stack;
}

std::vector<int> Grid::get_faces_of_zone_type(std::vector<int> zone_of_face)
{
	int sn_zone;
	std::vector<int> face_stack;

	for (int sn_face = 0; sn_face <= num_faces; sn_face++)
	{
		for (int ind_zone = 0; ind_zone < (int)face[sn_face].zone.size(); ind_zone++)
		{
			sn_zone = face[sn_face].zone[ind_zone];

			if (std::find(
				std::begin(zone_of_face),
				std::end(zone_of_face), sn_zone) != std::end(zone_of_face))
			{
				face_stack.push_back(sn_face);
				break; //to avoid potential duplications in node_stack if a node has multiple zonetypes existing in `zone_of_node`
			}
		}
	}

	return face_stack;
}

std::vector<int> Grid::get_cells_of_zone_type(std::vector<int> zone_of_cell)
{
	int sn_zone;
	std::vector<int> cell_stack;

	for (int sn_cell = 0; sn_cell <=num_cells; sn_cell++)
	{
		for (int ind_zone = 0; ind_zone < (int)cell[sn_cell].zone.size(); ind_zone++)
		{
			sn_zone = cell[sn_cell].zone[ind_zone];

			if (std::find(
				std::begin(zone_of_cell),
				std::end(zone_of_cell), sn_zone) != std::end(zone_of_cell))
			{
				zone_of_cell.push_back(sn_cell);
				break; //to avoid potential duplications in node_stack if a node has multiple zonetypes existing in `zone_of_node`
			}
		}
	}

	return zone_of_cell;
}



bool Grid::pos_inside_cell(
	double *pos,
	int sn_cell,
	double nondim_criterion_area_ratio)
{
	/*
	 * Determines whether a position is within a convex polygon
	 * by splitting it into triangles and check for each triangle
	 * individually.
	 */

	int sn_node_1, sn_node_2, sn_node_3;
	bool condition = false;

	for (int ind_tri = 0; ind_tri < (int)cell[sn_cell].node.size() - 2; ind_tri++)
	{
		// first vertex of triangle is fixed, while the other two shift until the whole polygon is covered.
		sn_node_1 = cell[sn_cell].node.at(0);
		sn_node_2 = cell[sn_cell].node.at(1 + ind_tri);
		sn_node_3 = cell[sn_cell].node.at(2 + ind_tri);

		if (PointInTriangle(pos,
							node[sn_node_1].props.pos,
							node[sn_node_2].props.pos,
							node[sn_node_3].props.pos,
							nondim_criterion_area_ratio))
			return true; // if pos resides in one tri it resides in polygon as well.
	}

	return condition;
}

std::tuple<bool, double *> Grid::face_intersection_with_line_segment(
	int sn_face,
	double *p1,
	double *p2)
{
	//----------------------------------
	// declarations
	//----------------------------------
	int sn_node_1, sn_node_2;
	bool has_intersection = false;
	double *intersection_coords;
	intersection_coords = new double[ndim];
	//----------------------------------

	sn_node_1 = face[sn_face].node[0];
	sn_node_2 = face[sn_face].node[1];
	std::tie(has_intersection, intersection_coords) = intersection_of_two_line_segments(
		node[sn_node_1].props.pos,
		node[sn_node_2].props.pos,
		p1, p2);

	return std::make_tuple(has_intersection, intersection_coords);
}

double Grid::distance_point_from_face(int sn_face, double *point)
{
	double dist;
	int sn_node_1 = face[sn_face].node[0];
	int sn_node_2 = face[sn_face].node[1];

	dist = distance_point_from_line_segment(
		node[sn_node_1].props.pos,
		node[sn_node_2].props.pos,
		point);

	return dist;
}

double *Grid::mirror_point_wrt_face(int sn_face, double *point)
{
	double *mirrored_point;
	int sn_node_1 = face[sn_face].node[0];
	int sn_node_2 = face[sn_face].node[1];

	mirrored_point = mirror_point_wrt_line_segment(
		node[sn_node_1].props.pos,
		node[sn_node_2].props.pos,
		point);

	return mirrored_point;
}

void Grid::order_cell_connectivity()
{
	if (opt_verbose)
		std::cout << "Ordering cells connectivity begins ...\n";

	/*
	 * Ordering faces and nodes on each cell based on right-hand rule.
	 * This function is not tested for 3d.
	 */
	using std::cout;
	using std::endl;

	std::vector<int> this_cell_node_ordered, this_cell_face_ordered;
	std::vector<int> this_node_vec;

	int sn_face, last_node_of_last_face_added;

	// excluding sn_cell=0 as it refers to a dummy item
	for (int sn_cell = 1; sn_cell <= num_cells; sn_cell++)
	{
		// init
		this_cell_face_ordered.clear();
		this_cell_node_ordered.clear();

		// processing each cell
		// in each iteration an appropriate face to be found followed by appending the face and its ordered nodes
		for (int sn_iter = 0; sn_iter < (int)cell[sn_cell].face.size(); sn_iter++)
		{
			/*
			 * check for the correct orientation (right - hand rule)
			 * Thumb will point to the "correct" side of a face.
			 * In a 2d case, first cell is on the left (thumb direction),
			 * and second cell is on the right of face.
			 */

			// iterating on faces of cell to find the correct face to append
			for (int ind_f = 0; ind_f < (int)cell[sn_cell].face.size(); ind_f++)
			{
				// init the node vector to store ordered node ids for this face
				this_node_vec.clear();

				sn_face = cell[sn_cell].face[ind_f];

				// skip if this face has already been added to stack
				if (std::find(
						this_cell_face_ordered.begin(),
						this_cell_face_ordered.end(), sn_face) != this_cell_face_ordered.end())
					continue;

				if (face[sn_face].cell[0] == sn_cell)
				{
					for (int ind_n = 0; ind_n < (int)face[sn_face].node.size(); ind_n++)
						this_node_vec.push_back(face[sn_face].node[ind_n]);
				}
				else
				{
					for (int ind_n = 0; ind_n < (int)face[sn_face].node.size(); ind_n++)
						this_node_vec.push_back(face[sn_face].node[face[sn_face].node.size() - ind_n - 1]);
				}

				// add the face if it is the first item
				if (this_cell_face_ordered.size() == 0)
					break;

				// otherwise, check whether this face has an edge (node for 2d) in common with the last added face
				last_node_of_last_face_added = this_cell_node_ordered.back();

				// Important step: if the correct face is found break out of loop
				if (last_node_of_last_face_added == this_node_vec[0])
					break;
			}

			// add the face if not already in the stack
			if (std::find(
					this_cell_face_ordered.begin(),
					this_cell_face_ordered.end(), sn_face) == this_cell_face_ordered.end())
			{
				this_cell_face_ordered.push_back(sn_face);

				// and nodes associated with the face
				for (int ind_n = 0; ind_n < (int)this_node_vec.size(); ind_n++)
					if (std::find(
							this_cell_node_ordered.begin(),
							this_cell_node_ordered.end(), this_node_vec[ind_n]) == this_cell_node_ordered.end())
						this_cell_node_ordered.push_back(this_node_vec[ind_n]);
			}
		}

		if (cell[sn_cell].face.size() != this_cell_face_ordered.size() or
			cell[sn_cell].node.size() != this_cell_node_ordered.size())
		{
			if (opt_verbose)
				cout << "Error in ordering faces/nodes of cell: " << sn_cell << endl;
			exit(1);
		}
		else
		{
			// reassigning nodes and faces of cell from the ordered vectors
			cell[sn_cell].face.clear();
			cell[sn_cell].node.clear();

			for (int ind_f = 0; ind_f < (int)this_cell_face_ordered.size(); ind_f++)
				cell[sn_cell].face.push_back(this_cell_face_ordered[ind_f]);

			for (int ind_n = 0; ind_n < (int)this_cell_node_ordered.size(); ind_n++)
				cell[sn_cell].node.push_back(this_cell_node_ordered[ind_n]);
		}
	}
	if (opt_verbose)
		std::cout << "End of process.\n\n";
}

void Grid::read_ansys_fluent_mesh_line_parsing(const char *fn, int report_every_num_line)
{
	/*
	 * This function is parsing line-by-line, so performs faster than
	 * `read_ansys_fluent_mesh_char_parsing`, which parses char-by-char.
	 * However, `read_ansys_fluent_mesh_char_parsing` can directly operate on
	 * an ascii Fluent mesh, while this function, `read_ansys_fluent_mesh_line_parsing`,
	 * would need the ascii Fluent mesh file to be saved as .txt file to avoid locale issues
	 * when parsing.
	 */

	// #pragma region Declarations
	using std::cout;
	using std::endl;
	using std::ifstream;
	using std::istringstream;
	using std::size_t;
	using std::string;

	ifstream infile;
	infile.open(fn);

	//*******************************************************************
	// declaration of vars
	//*******************************************************************
	string line = "", word[3] = {"", "", ""}, segment_level_1 = "", str[5] = {"", "", "", "", ""}, this_seg = "";
	size_t space_pos, this_start, this_len;
	istringstream ss[2];
	int zone_id, first_index, last_index, type, nd, elem_type,
		bc_type, face_type, num_zones_to_add;
	string zone_name, zone_type;
	int bodyFirstLine = -1, this_space_pos;
	int sn_node, sn_face, sn_cell;
	bool flag_nodes_body = false, flag_faces_body = false, flag_cells_body = false,
		 flag_comment = false;
	int tmp_ind[3], count_space, this_start_pos, this_elem_type;
	std::vector<std::string> tmp_stack;
	//*******************************************************************
	// #pragma endregion

	if (opt_verbose)
		cout << "Reading the mesh file starts ..." << endl;

	if (infile.is_open())
	{
		for (int sn_line = 1; getline(infile, line); sn_line++)
		{
			// trim endline of string
			trim_string_endl(line);

			// #pragma region Empty line
			//  Skipping empty lines
			if (line.length() == 0)
			{
				if (sn_line % report_every_num_line == 0)
					if (opt_verbose)
						cout << "line: " << sn_line << " | EMPTY" << endl;
				continue;
			}
			// #pragma endregion

			// #pragma region End of comments
			if (line == ")")
			{
				flag_comment = false;
				if (sn_line % report_every_num_line == 0)
					if (opt_verbose)
						cout << "line: " << sn_line << " | Comments section " << endl;
				continue;
			}
			// #pragma endregion

			// #pragma region Prep
			//  count number of spaces to determine number of items
			//  to be modified so a potential space at the end of line is not counted
			count_space = 0;
			for (int this_pos_sn = 0; this_pos_sn < (int)line.length(); this_pos_sn++)
				if (isspace(line.at(this_pos_sn)))
					count_space++;

			// index of first space in the line
			space_pos = line.find(' ');

			// first word of line
			word[0] = line.substr(0, space_pos);

			// last character of line
			word[1] = line[line.length() - 1];

			// segment within parentheses for index: 10, 12, 13
			this_start = space_pos + 2;
			this_len = line.length() - space_pos - 4;
			segment_level_1 = line.substr(this_start, this_len);

			// #pragma endregion

			// #pragma region End of a body section
			if (word[0] == "))")
			{
				if (flag_nodes_body)
				{
					if (sn_node != last_index)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error" << endl;
							cout << "line: " << sn_line << " | Index of last node extracted: " << sn_node << "\nis not consistent with the provided first and last indices: " << endl;
							cout << first_index << "\t" << last_index << endl;
						}
						exit(1);
					}
					flag_nodes_body = false;
					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | End of this nodes body section" << endl;
				}
				else if (flag_faces_body)
				{
					if (sn_face != last_index)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error" << endl;
							cout << "Index of last face extracted: " << sn_face << "\nis not consistent with the provided first and last indices: " << endl;
							cout << first_index << "\t" << last_index << endl;
						}
						exit(1);
					}
					flag_faces_body = false;
					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | End of this faces body section" << endl;
				}
				else if (flag_cells_body)
				{

					if (sn_cell != last_index)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error" << endl;
							cout << "Index of last cell extracted: " << sn_cell << "\nis not consistent with the provided first and last indices: " << endl;
							cout << first_index << "\t" << last_index << endl;
						}
						exit(1);
					}
					flag_cells_body = false;
					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | End of this cells body section" << endl;
				}
				else if (flag_comment)
				{
					if (opt_verbose)
						cout << "line: " << sn_line << " | Considered as comment " << endl;
				}
				else
				{
					if (opt_verbose)
						cout << "line: " << sn_line << " | Error in parsing." << endl;
					exit(1);
				}
				continue;
			}
			// #pragma endregion

			// #pragma region Index: 0 | Comments
			// currently, considers any line after '(0' as comments until a line with ')' reaches.
			// modifications are needed to support single line comments like: (0 Variables:)
			if (word[0] == "(0" or flag_comment)
			{
				flag_comment = true;
				if (sn_line % report_every_num_line == 0)
					if (opt_verbose)
						cout << "line: " << sn_line << " | Comments section " << endl;
			}
			// #pragma endregion

			// #pragma region Index: 2 | Dimensions
			else if (word[0] == "(2")
			{
				ndim = std::stoi(line.substr(space_pos + 1, line.length() - space_pos - 2));
				if (sn_line % report_every_num_line == 0)
					if (opt_verbose)
						cout << "line: " << sn_line << " | Dimensions | " << ndim << endl;
			}
			// #pragma endregion

			// #pragma region Index: 10 | Nodes
			else if (word[0] == "(10" or flag_nodes_body)
			{
				// #pragma region Index:10 | Declaration
				if (word[1] == ")")
				{
					ss[1].str("");
					ss[1].clear();
					ss[1].str(segment_level_1);
					ss[1] >> str[0] >> str[1] >> str[2] >> str[3];
					ss[1].str("");
					ss[1].clear();

					zone_id = std::stoi(str[0], 0, 16);
					first_index = std::stoi(str[1], 0, 16);
					last_index = std::stoi(str[2], 0, 16);
					type = std::stoi(str[3], 0, 16);

					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Nodes declaration | "
							 << zone_id << " " << first_index << " " << last_index << " " << type << endl;

					if (zone_id != 0 or first_index != 1 or type != 0)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error" << endl;
							cout << "Acceptable syntax for node declaration should be:" << endl;
							cout << "(10 (0 1 total_number_of_nodes 0))" << endl;
						}
						exit(1);
					}
					//------------------------------
					// process
					//------------------------------
					num_nodes = last_index;
					node = new Node[num_nodes + 1]; // Note: ansys fluent node index starts from 1, so we add a dummy node at index of 0.
				}
				// #pragma endregion
				//  #pragma region Index: 10 | Body | firstLine
				else if (word[1] == "(")
				{
					bodyFirstLine = sn_line;
					flag_nodes_body = true;

					ss[1].str("");
					ss[1].clear();
					ss[1].str(segment_level_1);
					ss[1] >> str[0] >> str[1] >> str[2] >> str[3] >> str[4];
					ss[1].str("");
					ss[1].clear();

					zone_id = std::stoi(str[0], 0, 16);
					first_index = std::stoi(str[1], 0, 16);
					last_index = std::stoi(str[2], 0, 16);
					type = std::stoi(str[3], 0, 16);
					nd = std::stoi(str[4]);

					if (zone_id == 0 or (last_index - first_index + 1) > num_nodes or type != 1 or nd != ndim)
					{
						if (opt_verbose)
							cout << "line: " << sn_line << " | Error in parsing nodes ...";
						exit(1);
					}
					//------------------------------
					// process
					//------------------------------
					if (zone_id >= (int)zone.size())
					{
						// Note: ansys fluent zone index starts from 1, so a dummy item is added at index of 0.
						// That is the reason for adding 1 to the right side.
						num_zones_to_add = zone_id - (int)zone.size() + 1;
						for (int this_zone_to_add_id = 0; this_zone_to_add_id < num_zones_to_add; this_zone_to_add_id++)
							zone.push_back(Zone());
					}
					for (int this_node_index = first_index; this_node_index <= last_index; this_node_index++)
					{
						// assuming a correct format file provided, the expensive check below is unnecessary
						// if (std::find(zone[zone_id].node.begin(), zone[zone_id].node.end(), this_node_index) == zone[zone_id].node.end())
						zone[zone_id].node.push_back(this_node_index);

						// add zone_id to node `zone` if not exists already
						if (std::find(node[this_node_index].zone.begin(), node[this_node_index].zone.end(), zone_id) == node[this_node_index].zone.end())
							node[this_node_index].zone.push_back(zone_id);

						if (node[this_node_index].type == node[this_node_index].type_init)
							node[this_node_index].type = type;
						else if (node[this_node_index].type != type)
						{
							if (opt_verbose)
							{
								cout << "line: " << sn_line << " | Error:" << endl;
								cout << "Different types given for the same node ...." << endl;
							}
							exit(1);
						}
					}

					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Nodes body first line | " << zone_id << " " << first_index << " " << last_index << " " << type << endl;
				}
				// #pragma endregion
				//  #pragma region Index: 10 | Body | Coordinates
				else
				{
					//---------------------------------------------------
					// parse coordinates of nodes
					//---------------------------------------------------
					double this_x, this_y;

					this_space_pos = (int)line.find(' ');
					this_x = std::stod(line.substr(0, this_space_pos));
					this_y = std::stod(line.substr(this_space_pos + 1));

					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Nodes body | " << this_x << " " << this_y << std::endl;

					sn_node = first_index + sn_line - bodyFirstLine - 1;
					node[sn_node].props.pos[0] = this_x;
					node[sn_node].props.pos[1] = this_y;
				}
				// #pragma endregion
			}
			// #pragma endregion

			// #pragma region Index: 12 | Cells
			else if (word[0] == "(12" or flag_cells_body)
			{
				if (word[1] == ")")
				{
					// #pragma region Index: 12 | Declaration
					if (count_space == 4)
					{
						ss[1].str("");
						ss[1].clear();
						ss[1].str(segment_level_1);
						ss[1] >> str[0] >> str[1] >> str[2] >> str[3];
						ss[1].str("");
						ss[1].clear();

						zone_id = std::stoi(str[0], 0, 16);
						first_index = std::stoi(str[1], 0, 16);
						last_index = std::stoi(str[2], 0, 16);
						type = std::stoi(str[3], 0, 16);

						if (sn_line % report_every_num_line == 0)
							if (opt_verbose)
								cout << "line: " << sn_line << " | Cells declaration | "
								 << zone_id << " " << first_index << " " << last_index << " " << type << endl;

						if (zone_id != 0 or first_index != 1 or type != 0)
						{
							if (opt_verbose)
							{
								cout << "line: " << sn_line << " | Error" << endl;
								cout << "Acceptable syntax for cell declaration should be:" << endl;
								cout << "(12 (0 1 total_number_of_cells 0))" << endl;
							}
							exit(1);
						}
						//------------------------------
						// process
						//------------------------------
						num_cells = last_index;
						cell = new Cell[num_cells + 1]; // Note: ansys fluent cell index starts from 1, so we add a dummy cell at index of 0.
					}
					// #pragma endregion
					//  #pragma region Index: 12 | Body (non-mixed mesh)
					else if (count_space == 5)
					{
						ss[1].str("");
						ss[1].clear();
						ss[1].str(segment_level_1);
						ss[1] >> str[0] >> str[1] >> str[2] >> str[3] >> str[4];
						ss[1].str("");
						ss[1].clear();

						zone_id = std::stoi(str[0], 0, 16);
						first_index = std::stoi(str[1], 0, 16);
						last_index = std::stoi(str[2], 0, 16);
						type = std::stoi(str[3], 0, 16);
						elem_type = std::stoi(str[4], 0, 16);

						if (sn_line % report_every_num_line == 0)
							if (opt_verbose)
								cout << "line: " << sn_line << " | Cells body (non-mixed mesh) | "
								 << zone_id << " " << first_index << " " << last_index << " " << type << " " << elem_type << endl;

						if (zone_id == 0 or (last_index - first_index + 1) > num_cells or type != 1 or elem_type == 0)
						{
							if (opt_verbose)
								cout << "line: " << sn_line << " | Error in parsing non-mixed cell element types.";
							exit(1);
						}
						//------------------------------
						// process
						//------------------------------
						if (zone_id >= (int)zone.size())
						{
							num_zones_to_add = zone_id - (int)zone.size() + 1;
							for (int this_zone_to_add_id = 0; this_zone_to_add_id < num_zones_to_add; this_zone_to_add_id++)
								zone.push_back(Zone());
						}
						for (int this_cell_index = first_index; this_cell_index <= last_index; this_cell_index++)
						{
							// assuming a correct format file provided, the expensive check below is unnecessary
							// if (std::find(zone[zone_id].cell.begin(), zone[zone_id].cell.end(), this_cell_index) == zone[zone_id].cell.end())
							zone[zone_id].cell.push_back(this_cell_index);

							// add zone_id to cell `zone` if not exists already
							if (std::find(cell[this_cell_index].zone.begin(), cell[this_cell_index].zone.end(), zone_id) == cell[this_cell_index].zone.end())
								cell[this_cell_index].zone.push_back(zone_id);

							// add element type to cell
							if (cell[this_cell_index].type == cell[this_cell_index].type_init)
								cell[this_cell_index].type = elem_type;
							else if (cell[this_cell_index].type != elem_type)
							{
								if (opt_verbose)
								{
									cout << "line: " << sn_line << " | Error:" << endl;
									cout << "Different element types given for the same cell." << endl;
								}
								exit(1);
							}
						}
					}
					// #pragma endregion
					else
					{
						if (opt_verbose)
							cout << "line: " << sn_line << " | Error in parsing." << endl;
						exit(1);
					}
				}
				// #pragma region Index: 12 | Body (mixed mesh) | firstLine
				else if (word[1] == "(")
				{
					bodyFirstLine = sn_line;
					flag_cells_body = true;

					ss[1].str("");
					ss[1].clear();
					ss[1].str(segment_level_1);
					ss[1] >> str[0] >> str[1] >> str[2] >> str[3] >> str[4];
					ss[1].str("");
					ss[1].clear();

					zone_id = std::stoi(str[0], 0, 16);
					first_index = std::stoi(str[1], 0, 16);
					last_index = std::stoi(str[2], 0, 16);
					type = std::stoi(str[3], 0, 16);
					elem_type = std::stoi(str[4], 0, 16);

					if (zone_id == 0 or (last_index - first_index + 1) > num_cells or type != 1 or elem_type != 0)
					{
						if (opt_verbose)
							cout << "line: " << sn_line << " | Error in parsing mixed cells ..." << endl;
						exit(1);
					}
					//------------------------------
					// process
					//------------------------------
					if (zone_id >= (int)zone.size())
					{
						// Note: ansys fluent zone index starts from 1, so a dummy item is added at index of 0.
						// That is the reason for adding 1 to the right side.
						num_zones_to_add = zone_id - (int)zone.size() + 1;
						for (int this_zone_to_add_id = 0; this_zone_to_add_id < num_zones_to_add; this_zone_to_add_id++)
							zone.push_back(Zone());
					}
					for (int this_cell_index = first_index; this_cell_index <= last_index; this_cell_index++)
					{
						// assuming a correct format file provided, the expensive check below is unnecessary
						// if (std::find(zone[zone_id].cell.begin(), zone[zone_id].cell.end(), this_cell_index) == zone[zone_id].cell.end())
						zone[zone_id].cell.push_back(this_cell_index);
					}

					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Cells body (mixed mesh) firstLine | "
							 << zone_id << " " << first_index << " " << last_index
							 << " " << type << " " << elem_type << std::endl;
				}
				// #pragma endregion
				//  #pragma region Index: 12 | Body (mixed mesh) | elem types
				else
				{
					//---------------------------------------------------
					// parse mixed element-types
					//---------------------------------------------------

					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Cells body (mixed mesh)" << endl;

					this_start_pos = 0; // start and end pos of elem_type
					sn_cell = first_index - 1;
					for (int this_pos = 0; this_pos < (int)line.length(); this_pos++)
					{
						if (isspace(line.at(this_pos)))
						{
							// this_seg: string related to an elem_type
							this_seg = line.substr(this_start_pos, this_pos - this_start_pos);
							if (this_seg.length() > 0)
							{
								sn_cell++;
								this_elem_type = std::stoi(this_seg, 0, 16);

								// add zone_id to cell `zone` if not exists already
								if (std::find(cell[sn_cell].zone.begin(), cell[sn_cell].zone.end(), zone_id) == cell[sn_cell].zone.end())
									cell[sn_cell].zone.push_back(zone_id);

								// add element type to cell
								if (cell[sn_cell].type == cell[sn_cell].type_init)
									cell[sn_cell].type = this_elem_type;
								else if (cell[sn_cell].type != this_elem_type)
								{
									if (opt_verbose)
									{
										cout << "line: " << sn_line << " | Error:" << endl;
										cout << "Different element types given for the same cell." << endl;
									}
									exit(1);
								}
							}
							this_start_pos = this_pos + 1; // start pos for next pair of coords is right after space char.
						}
					}
				}
			}
			// #pragma endregion
			// #pragma endregion

			// #pragma region Index: 13 | Faces
			else if (word[0] == "(13" or flag_faces_body)
			{
				// #pragma region Index: 13 | Declaration
				if (word[1] == ")")
				{
					ss[1].str("");
					ss[1].clear();
					ss[1].str(segment_level_1);
					ss[1] >> str[0] >> str[1] >> str[2] >> str[3];
					ss[1].str("");
					ss[1].clear();

					zone_id = std::stoi(str[0], 0, 16);
					first_index = std::stoi(str[1], 0, 16);
					last_index = std::stoi(str[2], 0, 16);
					bc_type = std::stoi(str[3], 0, 16);

					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Face declaration section | "
							 << zone_id << " " << first_index << " " << last_index << " " << bc_type << endl;

					if (zone_id != 0 or first_index != 1 or type != 0)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error ...";
							cout << "Acceptable syntax for face declaration should be:" << endl;
							cout << "(13 (0 1 total_number_of_faces 0))" << endl;
						}
						exit(1);
					}
					//------------------------------
					// process
					//------------------------------
					num_faces = last_index;
					face = new Face[num_faces + 1]; // Note: ansys fluent face index starts from 1, so we add a dummy face at index of 0.
				}
				// #pragma endregion
				//  #pragma region Index: 13 | Body | firstline
				else if (word[1] == "(")
				{
					bodyFirstLine = sn_line;
					flag_faces_body = true;

					ss[1].str("");
					ss[1].clear();
					ss[1].str(segment_level_1);
					ss[1] >> str[0] >> str[1] >> str[2] >> str[3] >> str[4];
					ss[1].str("");
					ss[1].clear();

					zone_id = std::stoi(str[0], 0, 16);
					first_index = std::stoi(str[1], 0, 16);
					last_index = std::stoi(str[2], 0, 16);
					bc_type = std::stoi(str[3], 0, 16);
					face_type = std::stoi(str[4], 0, 16);

					if (zone_id == 0 or (last_index - first_index + 1) > num_faces) // or type != 1 or elem_type != 0)
					{
						if (opt_verbose)
							cout << "line: " << sn_line << " | Error in parsing faces ..." << endl;
						exit(1);
					}
					//------------------------------
					// process
					//------------------------------
					if (zone_id >= (int)zone.size())
					{
						num_zones_to_add = zone_id - (int)zone.size() + 1;
						for (int this_zone_to_add_id = 0; this_zone_to_add_id < num_zones_to_add; this_zone_to_add_id++)
							zone.push_back(Zone());
					}
					for (int this_face_index = first_index; this_face_index <= last_index; this_face_index++)
					{
						// assuming a correct format file provided, the expensive check below is unnecessary
						// if (std::find(zone[zone_id].face.begin(), zone[zone_id].face.end(), this_face_index) == zone[zone_id].face.end())
						zone[zone_id].face.push_back(this_face_index);

						// add zone_id to face `zone` if not exists already
						if (std::find(face[this_face_index].zone.begin(), face[this_face_index].zone.end(), zone_id) == face[this_face_index].zone.end())
							face[this_face_index].zone.push_back(zone_id);
					}
					if (sn_line % report_every_num_line == 0)
						if (opt_verbose)
							cout << "line: " << sn_line << " | Face body first line | "
							 << zone_id << " " << first_index << " " << last_index
							 << " " << bc_type << " " << face_type << endl;
				}
				// #pragma endregion
				//  #pragma region Index: 13 | Body | facesConnectivity
				else
				{
					//---------------------------------------------------
					// parse faces section
					//---------------------------------------------------

					// stack of items in a line
					tmp_ind[0] = 0;
					tmp_stack.clear();
					for (int this_pos_sn = 0; this_pos_sn < (int)line.length(); this_pos_sn++)
						if (isspace(line.at(this_pos_sn)))
						{
							if (this_pos_sn - tmp_ind[0] > 0)
								tmp_stack.push_back(line.substr(tmp_ind[0], this_pos_sn - tmp_ind[0]));
							tmp_ind[0] = this_pos_sn + 1;
						}

					if (sn_line % report_every_num_line == 0)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Faces connectivity | ";
							for (int sn_item = 0; sn_item < (int)tmp_stack.size(); sn_item++)
								cout << std::stoi(tmp_stack.at(sn_item), 0, 16) << " ";
							cout << endl;
						}
					}

					// add connectivity information
					// Note: the last two items in stack of items are cell ids across the face
					sn_face = first_index + sn_line - bodyFirstLine - 1;
					for (int sn_item = 0; sn_item < (int)tmp_stack.size() - 2; sn_item++)
					{
						sn_node = std::stoi(tmp_stack.at(sn_item), 0, 16);

						// adding node to face if not already exists
						if (std::find(face[sn_face].node.begin(), face[sn_face].node.end(), sn_node) == face[sn_face].node.end())
							face[sn_face].node.push_back(sn_node);

						// adding face to node if not already exists
						if (std::find(node[sn_node].face.begin(), node[sn_node].face.end(), sn_face) == node[sn_node].face.end())
							node[sn_node].face.push_back(sn_face);

						for (int this_sn = 0; this_sn < 2; this_sn++)
						{
							// Note: for boundary faces, one of cells is zero
							sn_cell = std::stoi(tmp_stack.at(tmp_stack.size() - 2 + this_sn), 0, 16);

							// adding cell to face if not already exists
							if (std::find(face[sn_face].cell.begin(), face[sn_face].cell.end(), sn_cell) == face[sn_face].cell.end())
								face[sn_face].cell.push_back(sn_cell);

							// adding face to cell if not already exists
							if (std::find(cell[sn_cell].face.begin(), cell[sn_cell].face.end(), sn_face) == cell[sn_cell].face.end())
								cell[sn_cell].face.push_back(sn_face);

							// adding cell to node if not already exists
							if (std::find(node[sn_node].cell.begin(), node[sn_node].cell.end(), sn_cell) == node[sn_node].cell.end())
								node[sn_node].cell.push_back(sn_cell);

							// adding node to cell if not already exists
							if (std::find(cell[sn_cell].node.begin(), cell[sn_cell].node.end(), sn_node) == cell[sn_cell].node.end())
								cell[sn_cell].node.push_back(sn_node);
						}
					}

					// add boundary condition type to face
					if (face[sn_face].bc_type == face[sn_face].bc_type_init)
						face[sn_face].bc_type = bc_type;
					else if (face[sn_face].bc_type != bc_type)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error:" << endl;
							cout << "Different boundary condition types given for the same face ...." << endl;
						}
						exit(1);
					}

					// add face type to face
					if (face[sn_face].face_type == face[sn_face].face_type_init)
						face[sn_face].face_type = face_type;
					else if (face[sn_face].face_type != face_type)
					{
						if (opt_verbose)
						{
							cout << "line: " << sn_line << " | Error:" << endl;
							cout << "Different face types given for the same face ...." << endl;
						}
						exit(1);
					}
				}
				// #pragma endregion
			}
			// #pragma endregion

			// #pragma region Index:39 or 45 | Zones
			else if (word[0] == "(45" or word[0] == "(39")
			{
				// segment within parentheses for index: 10, 12, 13
				this_start = space_pos + 2;
				this_len = line.length() - space_pos - 6;
				segment_level_1 = line.substr(this_start, this_len);

				// stack of items in a line
				tmp_ind[0] = 0;
				tmp_stack.clear();
				for (int this_pos_sn = 0; this_pos_sn < (int)segment_level_1.length(); this_pos_sn++)
					if (isspace(segment_level_1.at(this_pos_sn)))
					{
						if (this_pos_sn - tmp_ind[0] > 0)
							tmp_stack.push_back(segment_level_1.substr(tmp_ind[0], this_pos_sn - tmp_ind[0]));
						tmp_ind[0] = this_pos_sn + 1;
					}
				if (!isspace(segment_level_1.at(segment_level_1.length() - 1)))
					tmp_stack.push_back(segment_level_1.substr(tmp_ind[0], segment_level_1.length() - tmp_ind[0]));

				zone_id = std::stoi(tmp_stack.at(0));
				zone_type = tmp_stack.at(1);
				zone_name = tmp_stack.at(2);

				if (zone_id >= (int)zone.size())
				{
					num_zones_to_add = zone_id - (int)zone.size() + 1;
					for (int this_zone_to_add_id = 0; this_zone_to_add_id < num_zones_to_add; this_zone_to_add_id++)
						zone.push_back(Zone());
				}

				// adding zone type
				if (zone.at(zone_id).type == zone.at(zone_id).type_init)
					zone.at(zone_id).type = zone_type;
				else if (zone.at(zone_id).type != zone_type)
				{
					if (opt_verbose)
					{
						cout << "line: " << sn_line << " | Error" << endl;
						cout << "Multiple different types have been provided for the same zone with id "
							<< zone_id << endl;
					}
					exit(1);
				}

				// adding zone name
				if (zone.at(zone_id).name == zone.at(zone_id).name_init)
					zone.at(zone_id).name = zone_name;
				else if (zone.at(zone_id).name != zone_name)
				{
					if (opt_verbose)
					{
						cout << "line: " << sn_line << " | Error" << endl;
						cout << "Multiple different names have been provided for the same zone with id "
							<< zone_id << endl;
					}
					exit(1);
				}

				// verbose
				if (sn_line % report_every_num_line == 0)
				{
					if (opt_verbose)
					{
						cout << "line: " << sn_line << " | Zone | ";
						for (int sn_item = 0; sn_item < (int)tmp_stack.size(); sn_item++)
							cout << tmp_stack.at(sn_item) << " ";
						cout << endl;
					}
				}
			}
			// #pragma endregion

			// #pragma region Invalid index
			else
			{
				if (opt_verbose)
					cout << "line: " << sn_line << " | Error in parsing." << endl;
				exit(1);
			}
			// #pragma endregion
		}
	}
	else
	{
		if (opt_verbose)
			std::cout << "\nCannot open mesh file: " << fn << "\n";
		exit(1);
	}

	num_zones = (int)zone.size() - 1;
	if (opt_verbose)
		cout << "Reading the mesh file ends." << endl
		 << endl;

}

void Grid::read_ansys_fluent_sol_at_nodes(
	const char *fn,
	double allowed_tolerance_node_position_snap,
	int report_every_num_line,
	const char *delimiter)
{
	/*
	 * This function reads the ansys fluent solution file at nodes
	 * and store the results in props of nodes.
	 * The props of cells can then calculated from the corresponding
	 * values of constituent nodes and/or be read from another cell
	 * center-based solution file.
	 */

	// #pragma region Declarations
	using std::cout;
	using std::endl;
	using std::ifstream;
	using std::string;
	// using std::istringstream;
	// using std::ofstream;
	using std::size_t;

	ifstream infile;
	infile.open(fn);

	//*******************************************************************
	// declaration of vars
	//*******************************************************************
	std::string delimiter_str(delimiter);
	string line, token;
	//string str[10];
	size_t pos;
	std::vector<string> stack_token;
	std::map<string, int> map_column;
	std::map<string, int>::iterator it_node_number, it_pos_x, it_pos_y,
		it_pressure, it_density, it_vel_x, it_vel_y, it_viscosity;

	int sn_node;
	double dist;
	double *this_pos;
	this_pos = new double[ndim];
	int flag_finding_close_node;
	//*******************************************************************
	// #pragma endregion

	// verbose
	if (opt_verbose)
	{
		cout << "Reading the solution file at nodes starts ..." << endl;
		cout << "sn_line" << delimiter
			<< "nodenumber" << delimiter
			<< "x" << delimiter
			<< "y" << delimiter
			<< "pressure" << delimiter
			<< "density" << delimiter
			<< "vel_x" << delimiter
			<< "vel_y" << delimiter
			<< "viscosity" << endl;
	}

	if (infile.is_open())
	{
		for (int sn_line = 1; getline(infile, line); sn_line++)
		{
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
				it_node_number = map_column.find("nodenumber");
				it_pos_x = map_column.find("x-coordinate");
				it_pos_y = map_column.find("y-coordinate");
				it_pressure = map_column.find("pressure");
				it_density = map_column.find("density");
				it_vel_x = map_column.find("x-velocity");
				it_vel_y = map_column.find("y-velocity");
				it_viscosity = map_column.find("viscosity-lam");

				// coordinates must be included in the file. Note: nodenumber is not reliable as it can be different from the mesh file.
				if (it_pos_x == map_column.end() or it_pos_y == map_column.end())
				{
					if (opt_verbose)
						cout << "Coordinates must be included in the file." << endl;
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
				// #pragma endregion

				// #pragma region Determining the corresponding node
				std::tie(sn_node, dist, flag_finding_close_node) = find_node_close_to_pos_spiral_algorithm(
				//std::tie(sn_node, dist, flag_finding_close_node) = find_node_close_to_pos(
					this_pos,
					allowed_tolerance_node_position_snap);

				if (flag_finding_close_node < 0)
				{
					if (opt_verbose)
					{
						std::cout << "\nError in `read_ansys_fluent_sol_at_nodes` after calling `find_node_close_to_pos_spiral_algorithm`:\n";
						std::cout << "Not finding a node close to position:\n";
						std::cout << "(" << this_pos[0] << ", " << this_pos[1] << ")\n";
						std::cout << "A potential remedy could be to increase the `half_width_search_window` which currently is " << endl;
						std::cout << half_width_search_window << "\n";
						std::cout << "Also, make sure the provided `allowed_tolerance_node_position_snap` is not too restrictive.\n"
							<< "`allowed_tolerance_node_position_snap` currently is: "
							<< allowed_tolerance_node_position_snap << "\n.";
					}
					exit(1);
				}

				// #pragma endregion

				// #pragma region pressure
				if (it_pressure != map_column.end())
					node[sn_node].props.pressure = std::stod(stack_token.at(it_pressure->second));
				// #pragma endregion

				// #pragma region density
				if (it_density != map_column.end())
					node[sn_node].props.density = std::stod(stack_token.at(it_density->second));
				// #pragma endregion

				// #pragma region vel_x
				if (it_vel_x != map_column.end())
					node[sn_node].props.vel[0] = std::stod(stack_token.at(it_vel_x->second));
				// #pragma endregion

				// #pragma region vel_y
				if (it_vel_y != map_column.end())
					node[sn_node].props.vel[1] = std::stod(stack_token.at(it_vel_y->second));
				// #pragma endregion

				// #pragma region viscosity
				if (it_viscosity != map_column.end())
					node[sn_node].props.viscosity = std::stod(stack_token.at(it_viscosity->second));
				// #pragma endregion

				// verbose
				if (sn_line % report_every_num_line == 0)
					if (opt_verbose)
						cout << sn_line << delimiter
						 << sn_node << delimiter
						 << node[sn_node].props.pos[0] << delimiter
						 << node[sn_node].props.pos[1] << delimiter
						 << node[sn_node].props.pressure << delimiter
						 << node[sn_node].props.density << delimiter
						 << node[sn_node].props.vel[0] << delimiter
						 << node[sn_node].props.vel[1] << delimiter
						 << node[sn_node].props.viscosity << endl;
			}
		}
	}
	else
	{
		if (opt_verbose)
			cout << "\nCould not open the nodes solution file: " << fn << endl;
		exit(1);
	}
	// end of parsing
	if (opt_verbose)
		cout << "Reading the solution file at nodes ends.\n\n";
}

void Grid::read_ansys_fluent_sol_at_cells(
	const char *fn,
	double nondim_allowed_tolerance_cell_center_position_snap,
	int report_every_num_line,
	bool strict_mode,
	const char *delimiter)
{
	/*
	 * This function reads the ansys fluent solution file at cell centers
	 * and store the results in props of cells.
	 */

	// #pragma region Declarations
	using std::cout;
	using std::endl;
	using std::ifstream;
	using std::size_t;
	using std::string;

	ifstream infile;
	infile.open(fn);

	//*******************************************************************
	// declaration of vars
	//*******************************************************************
	std::string delimiter_str(delimiter);
	string line, str[10], token;
	size_t pos;
	std::vector<string> stack_token;
	std::map<string, int> map_column;
	std::map<string, int>::iterator it_cell_number, it_pos_x, it_pos_y,
		it_pressure, it_density, it_vel_x, it_vel_y, it_viscosity, it_vol;

	int sn_cell;
	double dist, diff_area_nondim, diff_area, this_area;
	std::vector<double> dist_from_nodes;
	double *this_pos;
	this_pos = new double[ndim];
	std::vector<int> cell_search_stack;
	double EPSILON = 1e-9; // criterion for "NON-DIM" diff of calc. area of a cell and that read from a file.
						   //*******************************************************************
						   // #pragma endregion

	// verbose
	if (opt_verbose)
	{
		cout << "Reading the solution file at cells starts ..." << endl;
		cout << "sn_line" << delimiter
			<< "cellnumber" << delimiter
			<< "x" << delimiter
			<< "y" << delimiter
			<< "pressure" << delimiter
			<< "density" << delimiter
			<< "vel_x" << delimiter
			<< "vel_y" << delimiter
			<< "viscosity" << delimiter
			<< "volume" << endl;
	}

	int flag_normal_return = 0;
	int flag_not_finding_encomp_cell = -1;
	int flag_not_finding_close_node = -2;
	int flag_finding_encompassing_cell;

	if (infile.is_open())
	{
		for (int sn_line = 1; getline(infile, line); sn_line++)
		{
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
				it_cell_number = map_column.find("cellnumber");
				it_pos_x = map_column.find("x-coordinate");
				it_pos_y = map_column.find("y-coordinate");
				it_pressure = map_column.find("pressure");
				it_density = map_column.find("density");
				it_vel_x = map_column.find("x-velocity");
				it_vel_y = map_column.find("y-velocity");
				it_viscosity = map_column.find("viscosity-lam");
				it_vol = map_column.find("cell-volume");

				// coordinates must be included in the file. Note: nodenumber is not reliable as it can be different from the mesh file.
				if (it_pos_x == map_column.end() or it_pos_y == map_column.end())
				{
					if (opt_verbose)
						cout << "Coordinates must be included in the file." << endl;
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
				// #pragma endregion

				// #pragma region Determining the corresponding cell
				std::tie(sn_cell, dist_from_nodes, flag_finding_encompassing_cell, cell_search_stack) = find_cell_encompassing_pos(
					this_pos,
					flag_normal_return,
					flag_not_finding_encomp_cell,
					flag_not_finding_close_node
				);

				if (flag_finding_encompassing_cell == flag_not_finding_close_node)
				{
					if (opt_verbose)
						std::cout << "Error in `read_ansys_fluent_sol_at_cells` after calling `find_cell_encompassing_pos`:\n"
						<< "line " << sn_line
						<< "\nNot finding a close node."
						<< "\nPosition: ("
						<< this_pos[0] << ", "
						<< this_pos[1] << ")\n"
						<< "A remedy could be to increase the search window half width `half_width_search_window` which currently is:\n"
						<< half_width_search_window << std::endl;
					exit(1);
				}

				if (flag_finding_encompassing_cell == flag_not_finding_encomp_cell)
				{
					if (opt_verbose)
						std::cout << "Error in `read_ansys_fluent_sol_at_cells` after calling `find_cell_encompassing_pos`:\n"
						<< "line " << sn_line
						<< "\nPosition: ("
						<< this_pos[0] << ", "
						<< this_pos[1] << ")\n"
						<< "closest node found: " << dist_from_nodes[0]
						<< "\nNot finding an encompassing cell.\n";
					exit(1);
				}

				// strict mode: 
				// checking the cell center coordinates provided within the solution file
				// to make sure the right cell exists amongst the cells provided in the mesh file.
				if (strict_mode and it_vol != map_column.end())
				{
					this_area = std::stod(stack_token.at(it_vol->second));
					dist = distance_between_two_points(cell[sn_cell].props.pos, this_pos, ndim);

					// #pragma region checking the calc area of the cell and that provided from the sol file
					if (cell[sn_cell].area > 0)
					{
						diff_area = std::abs(cell[sn_cell].area - this_area);
						diff_area_nondim = diff_area / cell[sn_cell].area;

						if (diff_area_nondim > EPSILON)
						{
							if (opt_verbose)
							{
								cout << "line " << sn_line << " | Error" << endl;
								cout << "Cell " << sn_cell
									<< "\ncalc. area: " << cell[sn_cell].area
									<< "\nread area: " << this_area
									<< "\ndiff area: " << diff_area
									<< "\ndiff area nondim: " << diff_area_nondim
									<< endl;

								cout << "Cell " << sn_cell
									<< " with a center at:\n("
									<< cell[sn_cell].props.pos[0] << ", " << cell[sn_cell].props.pos[1] << "),\n"
									<< "found for position:\n("
									<< this_pos[0] << ", " << this_pos[1] << "),\n"
									<< "distance is\n"
									<< dist << endl
									<< "nondim distance is\n"
									<< dist / std::sqrt(this_area) << endl;
							}
							exit(1);
						}
					}
					// #pragma endregion

					// #pragma region Checking the discrepancy between the calc center of cell and that given from the sol file
					if (dist / std::sqrt(this_area) > nondim_allowed_tolerance_cell_center_position_snap)
					{
						if (opt_verbose)
							cout << "\nCell " << sn_cell
							<< " with a center at:\n("
							<< cell[sn_cell].props.pos[0] << ", " << cell[sn_cell].props.pos[1] << "),\n"
							<< "was found for position:\n("
							<< this_pos[0] << ", " << this_pos[1] << "),\n"
							<< "but the distance is\n"
							<< dist << endl
							<< "and nondim distance is\n"
							<< dist / std::sqrt(this_area) << endl
							<< "which is larger than the specified `nondim_allowed_tolerance_cell_center_position_snap`:"
							<< nondim_allowed_tolerance_cell_center_position_snap << endl;

						// exit(1);
					}
					// #pragma endregion
				}
				// #pragma endregion

				// #pragma region pressure
				if (it_pressure != map_column.end())
					cell[sn_cell].props.pressure = std::stod(stack_token.at(it_pressure->second));
				// #pragma endregion

				// #pragma region density
				if (it_density != map_column.end())
					cell[sn_cell].props.density = std::stod(stack_token.at(it_density->second));
				// #pragma endregion

				// #pragma region vel_x
				if (it_vel_x != map_column.end())
					cell[sn_cell].props.vel[0] = std::stod(stack_token.at(it_vel_x->second));
				// #pragma endregion

				// #pragma region vel_y
				if (it_vel_y != map_column.end())
					cell[sn_cell].props.vel[1] = std::stod(stack_token.at(it_vel_y->second));
				// #pragma endregion

				// #pragma region viscosity
				if (it_viscosity != map_column.end())
					cell[sn_cell].props.viscosity = std::stod(stack_token.at(it_viscosity->second));
				// #pragma endregion

				// #pragma region volume
				if (it_vol != map_column.end())
				{
					cell[sn_cell].area = std::stod(stack_token.at(it_vol->second));
				}
				// #pragma endregion

				// verbose
				if (sn_line % report_every_num_line == 0)
					if (opt_verbose)
						cout << sn_line << delimiter
						<< sn_cell << delimiter
						<< cell[sn_cell].props.pos[0] << delimiter
						<< cell[sn_cell].props.pos[1] << delimiter
						<< cell[sn_cell].props.pressure << delimiter
						<< cell[sn_cell].props.density << delimiter
						<< cell[sn_cell].props.vel[0] << delimiter
						<< cell[sn_cell].props.vel[1] << delimiter
						<< cell[sn_cell].props.viscosity << delimiter
						<< cell[sn_cell].area << endl;
			}
		}
	}
	else
	{
		if (opt_verbose)
			cout << "\nCould not open the solution file for cells center: " << fn << endl;
		exit(1);
	}
	if (opt_verbose)
		cout << "Reading the solution file at cells end.\n\n";
}

void Grid::report_sol_nodes(const char *fname, const char *delimiter)
{
	using std::cout;
	using std::endl;

	if (opt_verbose)
		cout << "Preparing report file for sol at nodes ... " << endl;

	std::ofstream fout(fname);

	fout << "nodenumber" << delimiter
		 << "x" << delimiter
		 << "y" << delimiter
		 << "pressure" << delimiter
		 << "density" << delimiter
		 << "vel_x" << delimiter
		 << "vel_y" << delimiter
		 << "viscosity"
		 << endl;

	for (int i = 1; i <= num_nodes; i++)
	{
		fout << i << delimiter
			 << node[i].props.pos[0] << delimiter
			 << node[i].props.pos[1] << delimiter
			 << node[i].props.pressure << delimiter
			 << node[i].props.density << delimiter
			 << node[i].props.vel[0] << delimiter
			 << node[i].props.vel[1] << delimiter
			 << node[i].props.viscosity
			 << endl;
	}
	if (opt_verbose)
		std::cout << "End of prep.\n\n";
}

void Grid::report_sol_cells(const char *fname, const char *delimiter)
{
	using std::cout;
	using std::endl;

	if (opt_verbose)
		cout << "Preparing report file for sol at cells ... " << endl;

	std::ofstream fout(fname);

	fout << "cellnumber" << delimiter
		 << "x" << delimiter
		 << "y" << delimiter
		 << "pressure" << delimiter
		 << "density" << delimiter
		 << "vel_x" << delimiter
		 << "vel_y" << delimiter
		 << "viscosity" << delimiter
		 << "cell-volume" << delimiter
		 << "vel_mag"
		 << endl;

	for (int i = 1; i <= num_cells; i++)
	{
		fout << i << delimiter
			 << cell[i].props.pos[0] << delimiter
			 << cell[i].props.pos[1] << delimiter
			 << cell[i].props.pressure << delimiter
			 << cell[i].props.density << delimiter
			 << cell[i].props.vel[0] << delimiter
			 << cell[i].props.vel[1] << delimiter
			 << cell[i].props.viscosity << delimiter
			 << cell[i].area << delimiter
			 << cell[i].props_aux.vel_mag
			 << endl;
	}
	if (opt_verbose)
		std::cout << "End of prep.\n\n";
}

void Grid::report_nodes_of_zones_single_file(
	const char* dir,
	const char* fname_base,
	const char* file_extension,
	const char* delimiter)
{
	/*
	* Note:
	* Zone numbers are 1-based, i.e., the first zone is 1.
	* Zone number of 0 is a dummy zone.
	*/

	std::ofstream fout;
	std::string fname;
	std::vector<int> zone_type, node_stack;
	
	//set the end char of dir
	std::string dir_str(dir);
	const char* endl_chars_to_be_trimmed = "\\ /";
	trim_string_endl(dir_str, endl_chars_to_be_trimmed);
	dir_str += Dir_Separator_Char;
	fname = dir_str + std::string(fname_base) + file_extension;

	if (opt_verbose)
		std::cout << "Preparing report file for nodes of zones (single file) ...\n";
	
	fout.open(fname);
	for (int ind_zone = 0; ind_zone < num_zones + 1; ind_zone++)//zone ids are 1-based, e.g. for a sys with 7 zones, the last zone's id is 7
	{
		zone_type.push_back(ind_zone);
		node_stack = get_nodes_of_zone_type(zone_type);
				
		for (int ind_node = 0; ind_node < (int)node_stack.size(); ind_node++)
			fout << node_stack[ind_node] << delimiter;
		
		zone_type.clear();
		node_stack.clear();
		fout << "\n";
	}
	fout.close();

	if (opt_verbose)
		std::cout << "End of prep.\n\n";
}


void Grid::report_nodes_of_zones_one_file_per_zone(
	const char* dir,
	const char* fname_base, 
	const char* file_extension, 
	const char* delimiter)
{
	/*
	* Note: 
	* Zone numbers are 1-based, i.e., the first zone is 1. 
	* Zone number of 0 is a dummy zone.
	*/

	std::ofstream fout;
	std::string fname;

	std::vector<int> zone_type, node_stack;


	//set the end char of dir
	std::string dir_str(dir);
	const char* endl_chars_to_be_trimmed = "\\ /";
	trim_string_endl(dir_str, endl_chars_to_be_trimmed);
	dir_str += Dir_Separator_Char;

	if (opt_verbose)
		std::cout << "Preparing report file for nodes of zones (one file per zone) ...\n";
	for (int ind_zone = 0; ind_zone < num_zones+1; ind_zone++)//zone ids are 1-based, e.g. for a sys with 7 zones, the last zone's id is 7
	{
		fname = dir_str + std::string(fname_base) + "_" + std::to_string(ind_zone) + file_extension;
		zone_type.push_back(ind_zone);
		node_stack = get_nodes_of_zone_type(zone_type);
		
		fout.open(fname);
		for (int ind_node = 0; ind_node < (int)node_stack.size(); ind_node++)
			fout << node_stack[ind_node] << delimiter;
		fout.close();
		fout.clear();
		zone_type.clear();
		node_stack.clear();
	}

	if (opt_verbose)
		std::cout << "End of prep.\n\n";
}


void Grid::report_mesh(const char *fname)
{
	using std::cout;
	using std::endl;

	if (opt_verbose)
		cout << "Preparing report mesh ... " << endl;

	std::ofstream fout(fname);

	fout << "number of nodes: " << num_nodes << endl;
	fout << "number of cells: " << num_cells << endl;
	fout << "number of face: " << num_faces << endl;
	fout << "number of zones: " << num_zones << endl;

	fout << endl;
	if (num_nodes <= 0)
	{
		if (opt_verbose)
			std::cout << "\nError:\nNumber of nodes: " << num_nodes << "\n";
		exit(1);
	}
	for (int i = 0; i <= num_nodes; i++)
	{
		fout << "node " << i << " | (x, y): ("
			 << node[i].props.pos[0] << ", "
			 << node[i].props.pos[1] << ")";

		fout << " | type: " << node[i].type;

		fout << " | Zones: ";
		for (int sn = 0; sn < (int)node[i].zone.size(); sn++)
			fout << node[i].zone.at(sn) << " ";

		fout << " | Faces: ";
		for (int sn = 0; sn < (int)node[i].face.size(); sn++)
			fout << node[i].face.at(sn) << " ";

		fout << " | Cells: ";
		for (int sn = 0; sn < (int)node[i].cell.size(); sn++)
			fout << node[i].cell.at(sn) << " ";

		fout << endl;
	}

	fout << endl;
	for (int i = 0; i <= num_cells; i++)
	{
		fout << "cell " << i << " | type: " << cell[i].type;

		fout << " | Zones: ";
		for (int sn = 0; sn < (int)cell[i].zone.size(); sn++)
			fout << cell[i].zone.at(sn) << " ";

		fout << " | Faces: ";
		for (int sn = 0; sn < (int)cell[i].face.size(); sn++)
			fout << cell[i].face.at(sn) << " ";

		fout << " | Nodes: ";
		for (int sn = 0; sn < (int)cell[i].node.size(); sn++)
			fout << cell[i].node.at(sn) << " ";

		fout << endl;
	}

	fout << endl;
	for (int i = 0; i <= num_faces; i++)
	{
		fout << "face " << i << " | bc_type: " << face[i].bc_type;

		fout << " | Zones: ";
		for (int sn = 0; sn < (int)face[i].zone.size(); sn++)
			fout << face[i].zone.at(sn) << " ";

		fout << " | Cells: ";
		for (int sn = 0; sn < (int)face[i].cell.size(); sn++)
			fout << face[i].cell.at(sn) << " ";

		fout << " | Nodes: ";
		for (int sn = 0; sn < (int)face[i].node.size(); sn++)
			fout << face[i].node.at(sn) << " ";

		fout << endl;
	}

	fout << endl;
	for (int i = 0; i < (int)zone.size(); i++)
	{
		fout << "zone " << i
			 << " | type: " << zone[i].type
			 << " | name: " << zone[i].name
			 << endl;
	}

	if (opt_verbose)
		cout << "Writing the mesh report file ends." << endl << endl;
}

//-------------------------------------------------------
// Getters
//-------------------------------------------------------

std::tuple<
	std::vector<double>, std::vector<double>,
	std::vector<double>, std::vector<double>,
	std::vector<double>, std::vector<double>,
	std::vector<double>>
Grid::get_sol_full()
{
	std::vector<double> x, y, p, rho, u, v, mu;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
	{
		x.push_back(node[sn_node].props.pos[0]);
		y.push_back(node[sn_node].props.pos[1]);
		p.push_back(node[sn_node].props.pressure);
		rho.push_back(node[sn_node].props.density);
		u.push_back(node[sn_node].props.vel[0]);
		v.push_back(node[sn_node].props.vel[1]);
		mu.push_back(node[sn_node].props.viscosity);
	}
	return std::make_tuple(x, y, p, rho, u, v, mu);
}

std::tuple<std::vector<double>, std::vector<double>> Grid::get_sol_pos()
{
	std::vector<double> x, y;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
	{
		x.push_back(node[sn_node].props.pos[0]);
		y.push_back(node[sn_node].props.pos[1]);
	}
	return std::make_tuple(x, y);
}

std::tuple<std::vector<double>, std::vector<double>> Grid::get_sol_vel()
{
	std::vector<double> u, v;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
	{
		u.push_back(node[sn_node].props.vel[0]);
		v.push_back(node[sn_node].props.vel[1]);
	}
	return std::make_tuple(u, v);
}

std::vector<double> Grid::get_sol_pressure()
{
	std::vector<double> p;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		p.push_back(node[sn_node].props.pressure);

	return p;
}

std::vector<double> Grid::get_sol_density()
{
	std::vector<double> rho;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		rho.push_back(node[sn_node].props.density);

	return rho;
}

std::vector<double> Grid::get_sol_viscosity()
{
	std::vector<double> mu;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		mu.push_back(node[sn_node].props.viscosity);

	return mu;
}

std::vector<double> Grid::get_sol_pos_x()
{
	std::vector<double> x;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		x.push_back(node[sn_node].props.pos[0]);

	return x;
}

std::vector<double> Grid::get_sol_pos_y()
{
	std::vector<double> y;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		y.push_back(node[sn_node].props.pos[1]);

	return y;
}

std::vector<double> Grid::get_sol_vel_x()
{
	std::vector<double> u;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		u.push_back(node[sn_node].props.vel[0]);

	return u;
}

std::vector<double> Grid::get_sol_vel_y()
{
	std::vector<double> v;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		v.push_back(node[sn_node].props.vel[1]);

	return v;
}

void Grid::get_sol_pos_x(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.pos[0]);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

void Grid::get_sol_pos_y(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.pos[1]);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

void Grid::get_sol_vel_x(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.vel[0]);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

void Grid::get_sol_vel_y(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.vel[1]);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

void Grid::get_sol_pressure(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.pressure);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

void Grid::get_sol_density(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.density);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

void Grid::get_sol_viscosity(double *output_array)
{
	std::vector<double> vec;

	for (int sn_node = 1; sn_node <= num_nodes; sn_node++)
		vec.push_back(node[sn_node].props.viscosity);

	std::memcpy(
		output_array,
		vec.data(),
		num_nodes * sizeof(double));
}

int Grid::get_num_nodes()
{
	printf("from inside get_num_nodes!\nNumber of nodes: %d\n", num_nodes);
	return num_nodes;
}

void Grid::eval_stats_faces()
{
	std::vector<double> face_len;
	int sn_node_1, sn_node_2;
	double dist;

	for (int sn_face = 1; sn_face <= num_faces; sn_face++)
	{
		sn_node_1 = face[sn_face].node[0];
		sn_node_2 = face[sn_face].node[1];
		dist = distance_between_two_points(node[sn_node_1].props.pos, node[sn_node_2].props.pos);
		face_len.push_back(dist);
	}

	auto min_max = std::minmax_element(begin(face_len), end(face_len));
	stat_face.min = *min_max.first;
	stat_face.max = *min_max.second;
	stat_face.avg = average(face_len);
}
