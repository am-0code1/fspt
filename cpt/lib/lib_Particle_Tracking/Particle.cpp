#include "Particle.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iostream>
// include <string>

constexpr double PI = 3.141592653589793238463;

Particle::Particle()
{
	int my_ID = 0;
	double my_diameter = 0;
	double my_density = 1000;

	set_id(my_ID);
	set_diameter(my_diameter);
	set_density(my_density);
	int ndim = props.NDIM;

	// props
	props.pos = new double[ndim];
	props.vel = new double[ndim];
	props.acc = new double[ndim];
	props.force = new double[ndim];
}

void Particle::set_id(int myID)
{
	id = myID;
}

double Particle::get_mass()
{
	return props.density * get_volume();
}

double Particle::get_volume()
{
	return PI / 6. * std::pow(props.diameter, 3.);
}

void Particle::add_to_hist(double this_time)
{
	int ndim = props.NDIM;

	double *this_pos, *this_vel, *this_acc, *this_force;
	this_pos = new double[ndim];
	this_vel = new double[ndim];
	this_acc = new double[ndim];
	this_force = new double[ndim];

	for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
	{
		this_pos[ind_dir] = props.pos[ind_dir];
		this_vel[ind_dir] = props.vel[ind_dir];
		this_acc[ind_dir] = props.acc[ind_dir];
		this_force[ind_dir] = props.force[ind_dir];
	}

	hist.time.push_back(this_time);

	hist.pos.push_back(this_pos);
	hist.vel.push_back(this_vel);
	hist.acc.push_back(this_acc);
	hist.force.push_back(this_force);
}

void Particle::set_force(double *this_force)
{
	for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
		props.force[ind_dir] = this_force[ind_dir];
}

void Particle::set_acc(double *this_acc)
{
	for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
		props.acc[ind_dir] = this_acc[ind_dir];
}

void Particle::set_vel(double *this_vel)
{
	for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
		props.vel[ind_dir] = this_vel[ind_dir];
}

void Particle::set_pos(double *this_pos)
{
	for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
		props.pos[ind_dir] = this_pos[ind_dir];
}

void Particle::set_density(double rho)
{
	props.density = rho;
}

void Particle::set_diameter(double d)
{
	props.diameter = d;
}

void Particle::report_hist(
	const char *fname_path,
	const char *file_extension,
	const char *delimiter,
	bool opt_append_id_to_file_name,
	const char *fname_base_ONLYIF_opt_append_id_to_file_name)
{
	using std::cout;
	using std::endl;

	cout << "Preparing particles " << id << "  history ... " << endl;

	std::string this_fname = std::string(fname_path);
	if (opt_append_id_to_file_name)
		this_fname += std::string(fname_base_ONLYIF_opt_append_id_to_file_name) + "_" + std::to_string(id);
	this_fname += file_extension;

	std::ofstream fout(this_fname);

	fout << "id" << delimiter << id << delimiter
		 << "diameter" << delimiter << props.diameter << delimiter
		 << "density" << delimiter << props.density
		 << endl;
	fout << "iter" << delimiter
		 << "time" << delimiter
		 << "x" << delimiter
		 << "y" << delimiter
		 << "vel_x" << delimiter
		 << "vel_y" << delimiter
		 << "acc_x" << delimiter
		 << "acc_y" << delimiter
		 << "force_x" << delimiter
		 << "force_y" << delimiter
		 << endl;

	for (int ind_iter = 0; ind_iter < (int)hist.time.size(); ind_iter++)
	{
		fout << ind_iter << delimiter
			 << hist.time[ind_iter] << delimiter;

		for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
			fout << hist.pos[ind_iter][ind_dir] << delimiter;

		for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
			fout << hist.vel[ind_iter][ind_dir] << delimiter;

		for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
			fout << hist.acc[ind_iter][ind_dir] << delimiter;

		for (int ind_dir = 0; ind_dir < props.NDIM; ind_dir++)
			fout << hist.force[ind_iter][ind_dir] << delimiter;

		fout << endl;
	}
}