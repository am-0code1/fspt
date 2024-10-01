// This file contains the 'main' function.

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>

#ifdef _WIN32
#include <ciso646> //for `and` and `or` keyword
#endif

#if defined(_WIN32)
#define Dir_Separator_Char "\\" //"/" also works on this system
#else
#define Dir_Separator_Char "/"
#endif

#if not defined(VERSION)
#define VERSION "1.0"
#endif


#include "Particle_Tracking.h"
#include "Utils.h"

//--- dev performance evaluation
#include <chrono>

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    using std::ifstream;

    //----------------------------------------------------
    //  default variables
    //----------------------------------------------------
    std::string input_filename = "config.txt";

    //OpenMP
    int num_threads = 1;
    
    //verbose 
    bool opt_output = true;
    bool opt_verbose = true;

    // search bins
    int n_bins[2] = {0, 0};
    double bin_len_dimensionless = 1.5;
    //congig search algorithm
    bool opt_original_ASG = false; //If true, the original variant of ASG search is used to locate hosting cell.Otherwise, the modified variant (two - step spiral) is used.
    bool opt_bin_cell_check_center_only = false; //Two implementation approaches to add cells to bins to be used in the original ASG variant: 1. (`opt_bin_cell_check_center_only`=true) center of cell to be examined : each cell resides only in one bin. 2. (`opt_bin_cell_check_center_only`=false) each node of cell to be examined : each cell can reside in one or multiple bins.
    // config for finding closest node to a position
    int half_width_search_window = 1;
    double allowed_tolerance_node_position_snap = 1e-11;
    // finding cell encompassing a position
    int encompassing_cell_spiral_level = 1;
    double nondim_criterion_area_ratio = 1e-11;
    // reading sol files
    bool opt_read_sol_at_cells_center = false;
    bool strict_mode = true;                                            // needed only when opt_read_sol_at_cells_center==true
    double nondim_allowed_tolerance_cell_center_position_snap = 1.0e-1; // needed only when opt_read_sol_at_cells_center==true
    // report files
    bool opt_report_mesh = false;
    bool opt_report_mesh_ordered = false;
    bool opt_report_sol_nodes = false;
    bool opt_report_sol_cells = false;
    bool opt_report_bins = false;
    bool opt_report_nodes_of_zones_one_file_per_zone = false;
    bool opt_report_nodes_of_zones_single_file = false;
    std::string dir_nodes_of_zones="";
    std::string fname_base_nodes_of_zones = "nodes_of_zone";
    std::string file_extension_nodes_of_zones = "";
    std::string delimiter_nodes_of_zones=",";
    // config verbose
    int report_every_num_line = int(1e5);
    // base directory
    std::string dir_base;
    // filenames
    std::string filename_mesh = "mesh.msh";
    std::string filename_sol_nodes = "sol_nodes";
    std::string filename_sol_cells = "sol_cells";
    //particle tracking related variables
    double time_duration = 1;
    double time_step_over_taw_relax = 1e0;
    double L_ref = 0;	//characteristic lenght of system (minimum feature size)
    double Vel_ref = 0; //characteristic velocity magnitude of system
    int max_iter = int(1e4);
    int verbose_every_num_step = 100;
    int props_eval_spiral_level = 0;
    double elastic_reflection_coeff = 1; //elastic reflection coeff. 0: inelastic, 1: pure elastic, between them: partially elastic reflection
    bool opt_time_integration_for_small_stokes = false;
    double stokes_number_threshold = 1e-1;
    double fixed_time_step = 1e-3;
    int num_points_on_perimeter = 8;
    bool opt_one_step_euler = true;
    bool opt_rk4 = false;
    //-----------------------------------
    //source of particle injection
    //for reading the particles file
    //-----------------------------------
    std::string particles_config_filename ="config_particles.txt";
    std::string particles_config_delimiter = ",";
    bool opt_release_at_flow_speed = false;
    int n_particles = 0;
    double range_diameter[2];
    double range_density[2];
    double source_pos[2], pos_end_1[2], pos_end_2[2];
    bool def_source_pos=false, def_pos_end_1 = false, def_pos_end_2 = false;
    //wall zone type(s)
    std::vector<int>wall_zone_type_vec = { 7 };



    //------------------------------------------------------
    // Parsing potential arguments from terminal
    //------------------------------------------------------
    // #pragma region Parsing arguments from terminal
    if (argc > 1)
    {
        cout << "Parsing the terminal arguments begins.\n";
        for (int ind_arg = 0; ind_arg < argc; ind_arg++)
            cout << argv[ind_arg] << " ";
        cout << endl;

        for (int ind_arg = 1; ind_arg < argc; ind_arg++)
        {
            if (strcmp(argv[ind_arg], "--filename") == 0 or
                strcmp(argv[ind_arg], "--fname") == 0 or
                strcmp(argv[ind_arg], "-f") == 0 or
                strcmp(argv[ind_arg], "-F") == 0)
            {
                input_filename = argv[ind_arg + 1];
                ind_arg++;
            }
            else if (strcmp(argv[ind_arg], "--version") == 0 or
                strcmp(argv[ind_arg], "-v") == 0 or
                strcmp(argv[ind_arg], "-V") == 0)
            {
                cout << "pt version: " << VERSION << endl;
                exit(0);
            }
        }
        cout << "Parsing the terminal arguments ends.\n\n";
    }
    // #pragma endregion




    //------------------------------------------------------
    // Parsing the input file
    //------------------------------------------------------
    // declaration for parsing input file
    std::string line = "", word[5] = { "", "", "", "", "" }, key = "";
    int tmp_ind = 0;
    std::vector<std::string> stack_items;
    //-------------------------------------

    // #pragma region Parsing input file
    ifstream infile;
    infile.open(input_filename);
    if (infile.is_open())
    {
        // parsing
        cout << "Parsing the input file begins... " << endl;

        for (int sn_line = 1; getline(infile, line); sn_line++)
        {
            // trim endline of string
            trim_string_endl(line);

            // Skipping empty lines
            if (line.length() == 0)
            {
                continue;
            }

            // Skipping comments
            if (line[0] == '#' or
                line[0] == '/')
            {
                continue;
            }

            // stack of args in a line
            tmp_ind = 0;
            stack_items.clear();
            for (int this_pos_sn = 0; this_pos_sn < (int)line.length(); this_pos_sn++)
                if (isspace(line.at(this_pos_sn)))
                {
                    if (this_pos_sn - tmp_ind > 0)
                        stack_items.push_back(line.substr(tmp_ind, this_pos_sn - tmp_ind));
                    tmp_ind = this_pos_sn + 1;
                }
            if (tmp_ind < (int)line.length())
                stack_items.push_back(line.substr(tmp_ind, (int)line.length() - tmp_ind));

            if (stack_items.size() < 2)
                continue;

            key = stack_items[0];

            if (key == "bin" or key == "bins")
            {
                assert((stack_items.size() == 3) && "`bins` requires 2 arguments for x and y directions.");
                n_bins[0] = (int)std::stod(stack_items[1]);
                n_bins[1] = (int)std::stod(stack_items[2]);
            }

            else if (key == "bin_len_dimensionless")
            {
                assert((stack_items.size() == 2) && "`bin_len_dimensionless` requires 1 argument.");
                bin_len_dimensionless = std::stod(stack_items[1]);
            }

            else if (key == "opt_verbose")
            {
                assert((stack_items.size() == 2) && "`opt_verbose` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_verbose = true;
                else
                    opt_verbose = false;
            }

            else if (key == "opt_output")
            {
                assert((stack_items.size() == 2) && "`opt_output` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_output = true;
                else
                    opt_output = false;
            }

            else if (key == "opt_original_ASG")
            {
                assert((stack_items.size() == 2) && "`opt_original_ASG` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_original_ASG = true;
                else
                    opt_original_ASG = false;
            }

            else if (key == "opt_bin_cell_check_center_only")
            {
                assert((stack_items.size() == 2) && "`opt_bin_cell_check_center_only` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_bin_cell_check_center_only = true;
                else
                    opt_bin_cell_check_center_only = false;
            }

            else if (key == "num_threads")
            {
                assert((stack_items.size() == 2) && "`num_threads` requires 1 argument.");
                num_threads = (int)std::stod(stack_items[1]);
            }

            else if (key == "half_width_search_window")
            {
                assert((stack_items.size() == 2) && "`half_width_search_window` requires 1 argument.");
                half_width_search_window = (int)std::stod(stack_items[1]);
            }

            else if (key == "allowed_tolerance_node_position_snap")
            {
                assert((stack_items.size() == 2) && "`allowed_tolerance_node_position_snap` requires 1 argument.");
                allowed_tolerance_node_position_snap = std::stod(stack_items[1]);
            }

            else if (key == "encompassing_cell_spiral_level")
            {
                assert((stack_items.size() == 2) && "`encompassing_cell_spiral_level` requires 1 argument.");
                encompassing_cell_spiral_level = (int)std::stod(stack_items[1]);
            }

            else if (key == "nondim_criterion_area_ratio")
            {
                assert((stack_items.size() == 2) && "`nondim_criterion_area_ratio` requires 1 argument.");
                nondim_criterion_area_ratio = std::stod(stack_items[1]);
            }

            else if (key == "opt_read_sol_at_cells_center")
            {
                assert((stack_items.size() == 2) && "`opt_read_sol_at_cells_center` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_read_sol_at_cells_center = true;
                else
                    opt_read_sol_at_cells_center = false;
            }

            else if (key == "strict_mode")
            {
                assert((stack_items.size() == 2) && "`strict_mode` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    strict_mode = true;
                else
                    strict_mode = false;
            }

            else if (key == "nondim_allowed_tolerance_cell_center_position_snap")
            {
                assert((stack_items.size() == 2) && "`nondim_allowed_tolerance_cell_center_position_snap` requires 1 argument.");
                nondim_allowed_tolerance_cell_center_position_snap = std::stod(stack_items[1]);
            }

            else if (key == "opt_report_mesh")
            {
                assert((stack_items.size() == 2) && "`opt_report_mesh` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_mesh = true;
                else
                    opt_report_mesh = false;
            }

            else if (key == "opt_report_mesh_ordered")
            {
                assert((stack_items.size() == 2) && "`opt_report_mesh_ordered` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_mesh_ordered = true;
                else
                    opt_report_mesh_ordered = false;
            }

            else if (key == "opt_report_sol_nodes")
            {
                assert((stack_items.size() == 2) && "`opt_report_sol_nodes` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_sol_nodes = true;
                else
                    opt_report_sol_nodes = false;
            }

            else if (key == "opt_report_sol_cells")
            {
                assert((stack_items.size() == 2) && "`opt_report_sol_cells` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_sol_cells = true;
                else
                    opt_report_sol_cells = false;
            }

            else if (key == "opt_report_bins")
            {
                assert((stack_items.size() == 2) && "`opt_report_bins` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_bins = true;
                else
                    opt_report_bins = false;
            }

            else if (key == "opt_report_nodes_of_zones_one_file_per_zone")
            {
                assert((stack_items.size() == 2) && "`opt_report_nodes_of_zones_one_file_per_zone` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_nodes_of_zones_one_file_per_zone = true;
                else
                    opt_report_nodes_of_zones_one_file_per_zone = false;
            }

            else if (key == "opt_report_nodes_of_zones_single_file")
            {
                assert((stack_items.size() == 2) && "`opt_report_nodes_of_zones_single_file` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_report_nodes_of_zones_single_file = true;
                else
                    opt_report_nodes_of_zones_single_file = false;
            }
            
            else if (key == "report_every_num_line")
            {
                assert((stack_items.size() == 2) && "`report_every_num_line` requires 1 argument.");
                report_every_num_line = (int)std::stod(stack_items[1]);
            }

            else if (key == "dir_base")
            {
                assert((stack_items.size() == 2) && "`dir_base` requires 1 argument.");
                dir_base = stack_items[1];
            }

            else if (key == "filename_mesh")
            {
                assert((stack_items.size() == 2) && "`filename_mesh` requires 1 argument.");
                filename_mesh = stack_items[1];
            }

            else if (key == "filename_sol_nodes")
            {
                assert((stack_items.size() == 2) && "`filename_sol_nodes` requires 1 argument.");
                filename_sol_nodes = stack_items[1];
            }

            else if (key == "filename_sol_cells")
            {
                assert((stack_items.size() == 2) && "`filename_sol_cells` requires 1 argument.");
                filename_sol_cells = stack_items[1];
            }

            else if (key == "particles_config_filename")
            {
                assert((stack_items.size() == 2) && "`particles_config_filename` requires 1 argument.");
                particles_config_filename = stack_items[1];
            }

            else if (key == "dir_nodes_of_zones")
            {
                assert((stack_items.size() == 2) && "`dir_nodes_of_zones` requires 1 argument.");
                dir_nodes_of_zones = stack_items[1];
            }

            else if (key == "fname_base_nodes_of_zones")
            {
                assert((stack_items.size() == 2) && "`fname_base_nodes_of_zones` requires 1 argument.");
                fname_base_nodes_of_zones = stack_items[1];
            }

            else if (key == "file_extension_nodes_of_zones")
            {
                assert((stack_items.size() == 2) && "`file_extension_nodes_of_zones` requires 1 argument.");
                file_extension_nodes_of_zones = stack_items[1];
            }

            else if (key == "delimiter_nodes_of_zones")
            {
                assert((stack_items.size() == 2) && "`delimiter_nodes_of_zones` requires 1 argument.");
                delimiter_nodes_of_zones = stack_items[1];
            }

            else if (key == "particles_config_delimiter")
            {
                assert((stack_items.size() >= 2) && "`particles_config_delimiter` requires 1 argument.");
                particles_config_delimiter = line;
                particles_config_delimiter.erase(particles_config_delimiter.find(key), key.size());
                ltrim(particles_config_delimiter);
                rtrim(particles_config_delimiter);
                particles_config_delimiter.erase(std::remove(particles_config_delimiter.begin(), particles_config_delimiter.end(), '\"'), particles_config_delimiter.end());
                particles_config_delimiter.erase(std::remove(particles_config_delimiter.begin(), particles_config_delimiter.end(), '\''), particles_config_delimiter.end());
            }

            //particle_tracking related keys           
            else if (key == "verbose_every_num_step")
            {
                assert((stack_items.size() == 2) && "`verbose_every_num_step` requires 1 argument.");
                verbose_every_num_step = (int)std::stod(stack_items[1]);
            }
            else if (key == "max_iter")
            {
                assert((stack_items.size() == 2) && "`max_iter` requires 1 argument.");
                max_iter = (int)std::stod(stack_items[1]);
            }
            else if (key == "time_duration")
            {
                assert((stack_items.size() == 2) && "`time_duration` requires 1 argument.");
                time_duration = std::stod(stack_items[1]);
            }
            else if (key == "time_step_over_taw_relax")
            {
                assert((stack_items.size() == 2) && "`time_step_over_taw_relax` requires 1 argument.");
                time_step_over_taw_relax = std::stod(stack_items[1]);
            }
            else if (key == "L_ref")
            {
                assert((stack_items.size() == 2) && "`L_ref` requires 1 argument.");
                L_ref = std::stod(stack_items[1]);
            }
            else if (key == "Vel_ref")
            {
                assert((stack_items.size() == 2) && "`Vel_ref` requires 1 argument.");
                Vel_ref = std::stod(stack_items[1]);
            }            
            else if (key == "props_eval_spiral_level")
            {
                assert((stack_items.size() == 2) && "`props_eval_spiral_level` requires 1 argument.");
                props_eval_spiral_level = (int)std::stod(stack_items[1]);
            }
            else if (key == "num_points_on_perimeter")
            {
                assert((stack_items.size() == 2) && "`num_points_on_perimeter` requires 1 argument.");
                num_points_on_perimeter = (int)std::stod(stack_items[1]);
            }
            else if (key == "elastic_reflection_coeff")
            {
                assert((stack_items.size() == 2) && "`elastic_reflection_coeff` requires 1 argument.");
                elastic_reflection_coeff = std::stod(stack_items[1]);
            }
            else if (key == "opt_one_step_euler")
            {
                assert((stack_items.size() == 2) && "`opt_one_step_euler` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_one_step_euler = true;
                else
                    opt_one_step_euler = false;
            }
            else if (key == "opt_rk4")
            {
                assert((stack_items.size() == 2) && "`opt_rk4` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_rk4 = true;
                else
                    opt_rk4 = false;
            }
            else if (key == "opt_time_integration_for_small_stokes")
            {
                assert((stack_items.size() == 2) && "`opt_time_integration_for_small_stokes` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_time_integration_for_small_stokes = true;
                else
                    opt_time_integration_for_small_stokes = false;
            }
            else if (key == "stokes_number_threshold")
            {
                assert((stack_items.size() == 2) && "`stokes_number_threshold` requires 1 argument.");
                stokes_number_threshold = std::stod(stack_items[1]);
            }
            else if (key == "fixed_time_step")
            {
                assert((stack_items.size() == 2) && "`fixed_time_step` requires 1 argument.");
                fixed_time_step = std::stod(stack_items[1]);
            }
            else if (key == "opt_release_at_flow_speed")
            {
                assert((stack_items.size() == 2) && "`opt_release_at_flow_speed` requires 1 argument.");
                if (stack_items[1] == "1" or stack_items[1] == "True" or stack_items[1] == "true")
                    opt_release_at_flow_speed = true;
                else
                    opt_release_at_flow_speed = false;
            }
            else if (key == "n_particles")
            {
                assert((stack_items.size() == 2) && "`n_particles` requires 1 argument.");
                n_particles = (int)std::stod(stack_items[1]);
            }
            else if (key == "range_diameter")
            {
                assert((stack_items.size() == 3) && "`range_diameter` requires 2 argument.");
                range_diameter[0] = std::stod(stack_items[1]);
                range_diameter[1] = std::stod(stack_items[2]);
            }
            else if (key == "range_density")
            {
                assert((stack_items.size() == 3) && "`range_density` requires 2 argument.");
                range_density[0] = std::stod(stack_items[1]);
                range_density[1] = std::stod(stack_items[2]);
            }
            else if (key == "source_pos")
            {
                assert((stack_items.size() == 3) && "`source_pos` requires 2 argument.");
                source_pos[0] = std::stod(stack_items[1]);
                source_pos[1] = std::stod(stack_items[2]);
                def_source_pos = true;
            }
            else if (key == "pos_end_1")
            {
                assert((stack_items.size() == 3) && "`pos_end_1` requires 2 argument.");
                pos_end_1[0] = std::stod(stack_items[1]);
                pos_end_1[1] = std::stod(stack_items[2]);
                def_pos_end_1 = true;
            }
            else if (key == "pos_end_2")
            {
                assert((stack_items.size() == 3) && "`pos_end_2` requires 2 argument.");
                pos_end_2[0] = std::stod(stack_items[1]);
                pos_end_2[1] = std::stod(stack_items[2]);
                def_pos_end_2 = true;
                }
            else if (key == "wall_zone_type")
            {
                assert((stack_items.size() >= 2) && "`wall_zone_type` requires 1 or more argument.");
                wall_zone_type_vec.clear();
                for (int sn_item=1; sn_item <stack_items.size(); sn_item++)
                    wall_zone_type_vec.push_back((int)std::stod(stack_items[sn_item]));
            }
            
            else
            {
                cout << "Invalid keyword: " << key << endl;
                return (1);
            }
        }
    }
    else
    {
        cout << "\nCould not open the input file: " << input_filename << endl;
        return (1);
    }
    // end of parsing
    cout << "Parsing the input file ends.\n\n";

    //set wall_zone_type array to be passed to track_particle method later
    int* wall_zone_type, num_wall_zone_type;
    num_wall_zone_type = (int)wall_zone_type_vec.size();
    wall_zone_type = new int(num_wall_zone_type);
    for (int ind_wzt=0; ind_wzt< num_wall_zone_type; ind_wzt++)
        wall_zone_type[ind_wzt]=wall_zone_type_vec[ind_wzt];
    // #pragma endregion




    //------------------------------------------------------
    // config the paths
    //------------------------------------------------------
    //set the end char of dir
    const char* endl_chars_to_be_trimmed = "\\ /";
    trim_string_endl(dir_base, endl_chars_to_be_trimmed);
    dir_base += Dir_Separator_Char;

    // input
    std::string path_mesh = dir_base + "inp" + Dir_Separator_Char + filename_mesh;
    std::string path_sol_nodes = dir_base + "inp" + Dir_Separator_Char + filename_sol_nodes;
    std::string path_sol_cells = dir_base + "inp" + Dir_Separator_Char + filename_sol_cells;
    std::string path_particle_config_file = dir_base + "inp" + Dir_Separator_Char + particles_config_filename;
    // report
    std::string path_report_mesh = dir_base + "report" + Dir_Separator_Char + "mesh";
    std::string path_report_mesh_ordered_connectivity = dir_base + "report" + Dir_Separator_Char + "mesh_ordered_connectivity";
    std::string path_report_bins = dir_base + "report" + Dir_Separator_Char + "bins";
    std::string path_report_sol_nodes = dir_base + "report" + Dir_Separator_Char + "nodes";
    std::string path_report_sol_cells = dir_base + "report" + Dir_Separator_Char + "cells";
    if (dir_nodes_of_zones.length() == 0)
        dir_nodes_of_zones = dir_base + "report";
    // particles
    std::string dir_base_particle_hist = dir_base + "particle" + Dir_Separator_Char;





    //------------------------------------------------------
    // Beginning of particle tracking part
    //------------------------------------------------------
    Particle_Tracking pt(
        opt_original_ASG,
        opt_bin_cell_check_center_only,
        half_width_search_window,
        encompassing_cell_spiral_level,
        nondim_criterion_area_ratio,
        props_eval_spiral_level,
        opt_output,
        opt_verbose,
        num_threads
    );
    pt.set_grid_and_sol(
        path_mesh.c_str(),
        path_sol_nodes.c_str(),
        path_sol_cells.c_str(),
        path_report_mesh.c_str(),
        path_report_mesh_ordered_connectivity.c_str(),
        path_report_sol_nodes.c_str(),
        path_report_sol_cells.c_str(),
        path_report_bins.c_str(),
        opt_read_sol_at_cells_center,
        n_bins,
        allowed_tolerance_node_position_snap,
        nondim_allowed_tolerance_cell_center_position_snap,
        strict_mode,
        opt_report_mesh,
        opt_report_mesh_ordered,
        opt_report_sol_nodes,
        opt_report_sol_cells,
        opt_report_bins,
        opt_report_nodes_of_zones_one_file_per_zone,
        opt_report_nodes_of_zones_single_file,
        dir_nodes_of_zones.c_str(),
        fname_base_nodes_of_zones.c_str(),
        file_extension_nodes_of_zones.c_str(),
        delimiter_nodes_of_zones.c_str(),
        report_every_num_line,
        bin_len_dimensionless
        );

    /*
     * Particle tracking
     */

    //Introducing particles using two mechanisms: 1. streakline, and 2. single point injection source.
    //For both cases, `n_particles` needs to be explicitly defined and to be larger than 0.
    if (n_particles > 0)
    {
        //streakline injection
        if (def_pos_end_1 and def_pos_end_2)
            pt.set_streakline(
                n_particles,
                pos_end_1,
                pos_end_2,
                range_diameter,
                range_density,
                opt_release_at_flow_speed);
        
        //single point source of injection
        else if (def_source_pos)
            pt.set_particle_source(
                n_particles,
                source_pos,
                range_diameter,
                range_density,
                opt_release_at_flow_speed);
    }
    
    //Introducing particles using a file.
    //By default, a particles file is read, but no error is thrown if the file does not exist or cannot be read.
    pt.read_particle_file(path_particle_config_file.c_str(), particles_config_delimiter.c_str());

    int particle_id_arr[] = {0};   // id of particles to be tracked
    int n_particles_to_track = -1; // 0 or a negative val --> all particles to be tracked, in which case `particle_id_arr` will be superseded.


    //--- dev performance measurement
    using picoseconds = std::chrono::duration<long long, std::pico>;
    auto start = std::chrono::high_resolution_clock::now();

    pt.track_particle(
        dir_base_particle_hist.c_str(),
        time_duration,
        time_step_over_taw_relax,
        wall_zone_type,
        particle_id_arr,
        L_ref,
        Vel_ref,
        num_wall_zone_type,
        n_particles_to_track,
        max_iter,
        verbose_every_num_step,
        elastic_reflection_coeff,
        opt_time_integration_for_small_stokes,
        stokes_number_threshold,
        fixed_time_step,
        num_points_on_perimeter,
        opt_one_step_euler,
        opt_rk4
    );

    //--- dev performance measurement
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_ps = picoseconds{ stop - start };
    cout << "\nElapsed time total: " << duration_ps.count() << " picoseconds" << endl;
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << "Elapsed time total: " << duration_us.count() << " microseconds" << endl;
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout << "Elapsed time total: " << duration_ms.count() << " milliseconds" << endl;

    return 0;
}
