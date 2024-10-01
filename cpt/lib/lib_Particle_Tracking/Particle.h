#pragma once
#include "Props.h"
#include <vector>

class Particle
{

public:
	Particle();

	void set_id(int);
	void set_diameter(double);
	void set_density(double);
	void set_pos(double *);
	void set_vel(double *);
	void set_acc(double *);
	void set_force(double *);
	void add_to_hist(double);
	double get_mass();
	double get_volume();
	void report_hist(
		const char *fname_path,
		const char *file_extension = "",
		const char *delimiter = ",",
		bool opt_append_id_to_file_name = false,
		const char *fname_base_ONLYIF_opt_append_id_to_file_name = "particle");

	// vars
	Props_particle props;
	Hist hist;

	int id;
};
