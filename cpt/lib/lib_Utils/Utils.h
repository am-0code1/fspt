#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <tuple>
#include <vector>
#include <numeric>

// Most math utility functions are written
// for a 2-D space.
//-----------------------------------------
const int NDIM = 2;
//------------------------------------------

void ltrim(std::string& s); // trim from start (in place)
void rtrim(std::string& s); // trim from end (in place)

void trim_string_endl(std::string &str, const char endl_char[] = "\n\r");

double distance_between_two_points(double *p1, double *p2, int dim = NDIM);
double area_triangle(double *v1, double *v2, double *v3);
bool PointInTriangle(double *p, double *v1, double *v2, double *v3, double EPSILON = 1e-12);

std::tuple<double, double, double> get_line_param(double *p1, double *p2); // returns the parameters of line in standard form of Ax + By = C

std::tuple<bool, double *> intersection_of_two_lines(double A1, double B1, double C1, double A2, double B2, double C2);           // intersection of two standard fomr lines
std::tuple<bool, double *> intersection_of_two_line_segments(double *seg1_v1, double *seg1_v2, double *seg2_v1, double *seg2_v2); // intersection of two line segments in 2D

double distance_point_from_line(double A, double B, double C, double *point);           // distance of a point from a line
double distance_point_from_line_segment(double *seg_v1, double *seg_v2, double *point); // distance of a point from a line segment

double *projection_of_point_on_line(double A, double B, double C, double *point); // projection of a point onto a line

double *mirror_point_wrt_line(double A, double B, double C, double *point);           // mirror a point wrt a line
double *mirror_point_wrt_line_segment(double *seg_v1, double *seg_v2, double *point); // mirror a point wrt a line segment

double len_vec(double *vec, int ndim = NDIM);
double *normal_vec(double *vec, int ndim = NDIM);

double dot(double *vec_1, double *vec_2, int ndim = NDIM);

double average(std::vector<double> const& v);

#endif
