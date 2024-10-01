#pragma once
#include <vector>

#ifndef struct_props
#define struct_props
struct Props
{
	int NDIM = 2;
	double* pos, * vel, density, pressure, viscosity;
};
#endif //struct_props

#ifndef struct_props_aux
#define struct_props_aux
struct Props_aux
{
	double vel_mag;
};
#endif //struct_props_aux

#ifndef struct_props_particle
#define struct_props_particle
struct Props_particle
{
	int NDIM = 2;
	double* pos, * vel, * acc, * force, density, diameter;
};
#endif //struct_props_particle

#ifndef struct_hist
#define struct_hist
struct Hist
{
	int NDIM = 2;
	std::vector<double> time;
	std::vector<double*> pos; 
	std::vector<double*> vel;
	std::vector<double*> acc;
	std::vector<double*> force;
};
#endif //struct_hist

#ifndef struct_stat
#define struct_stat
struct Stat
{
	double min, max, avg;
};
#endif //struct_stat

