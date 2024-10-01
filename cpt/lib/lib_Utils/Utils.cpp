#include "Utils.h"
#include <cmath>
#include <iomanip>
#include <iostream>

#ifdef _WIN32
#include <ciso646> //for `and` and `or` keyword
#endif

// trim from start (in place)
void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
        }));
}

// trim from end (in place)
void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
        }).base(), s.end());
}

void trim_string_endl(std::string &str, const char endl_char[])
{
    // A function to trim endline chars of a string to make code compatible with
    // both Windows OS, wherein lines are terminated with "\r\n"
    // and Linux OS, wherein lines are terminated with "\n".
    // static const char endl_char[] = "\n\r"

    str.erase(str.find_last_not_of(endl_char) + 1);
}

double distance_between_two_points(double *p1, double *p2, int dim)
{
    using std::pow;
    using std::sqrt;
    double val = 0;

    for (int ind_comp = 0; ind_comp < dim; ind_comp++)
        val += pow(p2[ind_comp] - p1[ind_comp], 2.);

    return sqrt(val);
}

// area of triangle
double area_triangle(double *v1, double *v2, double *v3)
{
    return std::abs((v1[0] * (v2[1] - v3[1]) + v2[0] * (v3[1] - v1[1]) + v3[0] * (v1[1] - v2[1])) / 2.0);
}

// determine whether point p resides inside triange with vertices v1, v2, and v3
bool PointInTriangle(double *p, double *v1, double *v2, double *v3, double EPSILON)
{
    double area, area_1, area_2, area_3, sum_area, diff_area, diff_area_nondim;

    area = area_triangle(v1, v2, v3);

    area_1 = area_triangle(p, v1, v2);
    area_2 = area_triangle(p, v1, v3);
    area_3 = area_triangle(p, v2, v3);
    sum_area = area_1 + area_2 + area_3;

    diff_area = std::abs(sum_area - area);
    diff_area_nondim = diff_area / area;

    if (diff_area_nondim < EPSILON)
        return true;
    else
        return false;
}

// intersection of two line segments in 2D
std::tuple<bool, double *> intersection_of_two_line_segments(double *p1, double *p2, double *p3, double *p4)
{
    //---------------------------------------
    // Declarations
    //---------------------------------------
    double A1, B1, C1, A2, B2, C2;
    bool has_intersection = false;
    double *intersection_coords;
    double EPSILON = 1e-12; // to deal with the double precision issues
    //---------------------------------------

    // line 1
    std::tie(A1, B1, C1) = get_line_param(p1, p2);

    // line 2
    std::tie(A2, B2, C2) = get_line_param(p3, p4);

    std::tie(has_intersection, intersection_coords) = intersection_of_two_lines(
        A1, B1, C1,
        A2, B2, C2);

    if (has_intersection)
    {
        // segment 1: first 4 conditions
        // segment 2: second 4 conditions
        if (
            !(
                intersection_coords[0] >= std::min(p1[0], p2[0]) - EPSILON and
                intersection_coords[0] <= std::max(p1[0], p2[0]) + EPSILON and
                intersection_coords[1] >= std::min(p1[1], p2[1]) - EPSILON and
                intersection_coords[1] <= std::max(p1[1], p2[1]) + EPSILON and
                intersection_coords[0] >= std::min(p3[0], p4[0]) - EPSILON and
                intersection_coords[0] <= std::max(p3[0], p4[0]) + EPSILON and
                intersection_coords[1] >= std::min(p3[1], p4[1]) - EPSILON and
                intersection_coords[1] <= std::max(p3[1], p4[1]) + EPSILON))
            has_intersection = false;
    }

    return std::make_tuple(has_intersection, intersection_coords);
}

std::tuple<bool, double *> intersection_of_two_lines(
    double A1, double B1, double C1,
    double A2, double B2, double C2)
{
    //---------------------------------------
    // Declarations
    double det;
    bool has_intersection = false;
    double *intersection_coords;
    intersection_coords = new double[2];
    //---------------------------------------

    det = A1 * B2 - A2 * B1;
    if (det == 0)
    {
        // Lines are parallel. discarding the case of infinite number of solutions.
        has_intersection = false;
    }
    else
    {
        intersection_coords[0] = (B2 * C1 - B1 * C2) / det;
        intersection_coords[1] = (A1 * C2 - A2 * C1) / det;
        has_intersection = true;
    }

    return std::make_tuple(has_intersection, intersection_coords);
}

std::tuple<double, double, double> get_line_param(double *p1, double *p2)
{
    // returns the parameters of line in standard form of
    // Ax + By = C

    double A = p2[1] - p1[1];
    double B = p1[0] - p2[0];
    double C = A * p1[0] + B * p1[1];

    return std::make_tuple(A, B, C);
}

double distance_point_from_line_segment(double *p1, double *p2, double *p3) // finding the distance of a point from a line segment
{
    // p1 and p2 are endpoints of a line segment

    double A, B, C, dist;
    std::tie(A, B, C) = get_line_param(p1, p2);

    dist = distance_point_from_line(A, B, C, p3);

    return dist;
}

double distance_point_from_line(double A, double B, double C, double *point) // finding the distance of a point from a line
{
    double dist = std::abs(A * point[0] + B * point[1] - C) / std::sqrt(std::pow(A, 2) + std::pow(B, 2));

    return dist;
}

double *projection_of_point_on_line(double A, double B, double C, double *point)
{
    //--------------------------
    // declarations
    //--------------------------
    bool has_intersection;
    double *proj;
    //--------------------------

    // a line perpendicular to the line Ax+By=C has a standard form of
    // Bx-Ay=D, where D can be calculated as
    double D = B * point[0] - A * point[1];

    // projection point is the intersection of two lines
    std::tie(has_intersection, proj) = intersection_of_two_lines(
        A, B, C,
        B, -A, D);

    return proj;
}

double *mirror_point_wrt_line_segment(double *p1, double *p2, double *p3)
{
    // p1 and p2 are endpoints of the mirror line segment

    double A, B, C, *mirrored_point;
    std::tie(A, B, C) = get_line_param(p1, p2);

    mirrored_point = mirror_point_wrt_line(A, B, C, p3);

    return mirrored_point;
}

double *mirror_point_wrt_line(double A, double B, double C, double *point)
{
    //-----------------------------
    // declarations
    //-----------------------------
    double *proj, *mirrored_point;
    mirrored_point = new double[2];
    //-----------------------------

    proj = projection_of_point_on_line(A, B, C, point);

    for (int ind_dir = 0; ind_dir < 2; ind_dir++)
        mirrored_point[ind_dir] = point[ind_dir] + 2 * (proj[ind_dir] - point[ind_dir]);

    return mirrored_point;
}

double len_vec(double *vec, int ndim)
{
    double len = 0;

    for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
        len += pow(vec[ind_dir], 2);
    len = sqrt(len);

    return len;
}

double *normal_vec(double *vec, int ndim)
{
    // declaration
    double len, *normaliezed_vec;
    normaliezed_vec = new double[ndim];
    //--------------------------------

    len = len_vec(vec, ndim);
    for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
        normaliezed_vec[ind_dir] = vec[ind_dir] / len;

    return normaliezed_vec;
}

double dot(double *vec_1, double *vec_2, int ndim)
{
    double val = 0;

    for (int ind_dir = 0; ind_dir < ndim; ind_dir++)
        val += vec_1[ind_dir] * vec_2[ind_dir];

    return val;
}

double average(std::vector<double> const& v) {
    if (v.empty()) {
        return 0;
    }

    auto const count = static_cast<float>(v.size());

    //accumulate can also be used, but is not as stable and fast as reduce.
    //reduce requires at least C++17
    //
    // std::accumulate(v.begin(), v.end(), 0.0) / count;

    return  std::reduce(v.begin(), v.end()) / count;
}
