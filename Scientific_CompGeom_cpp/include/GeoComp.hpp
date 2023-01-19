#include <linalg.hpp>
#include <fstream>
#include <algorithm>

struct mesh {
    vector<vector<double>> coords;
    vector<vector<int>> elems;
    vector<vector<int>> sibhfs;
};

struct Triangulation {
    vector<vector<double>> coords;
    vector<double> param;
    vector<vector<int>> elems;
    vector<vector<int>> sibhfs;
    int nelems;
    vector<vector<int>> facets;
    vector<bool> delete_elem;
    vector<int> badtris;
    vector<int> vedge;
    vector<bool> bwork;
    vector<bool> on_boundary;
};

struct Spline {
    vector<vector<double>> coords;
    int nv;
    int degree;
    vector<vector<double>> xweights;
    vector<vector<double>> yweights;
    vector<double> params;
};

// spline functions
Spline spline_init(const vector<vector<double>> &xs, const vector<bool> &corners, int degree);
vector<double> spline_var(Spline* spl, double t, int order=0);
vector<double> spline_point_segment(Spline* spl, double a, double b, double ratio);

// mesh smoothing functions
void mesh_smoothing_2d(Triangulation* mesh, vector<bool> no_move, function<double(vector<double>)> r_ref, double mu);

// mesh refinement functions
void GeoComp_refine(Triangulation* DT, function<double(vector<double>)> r_ref, Spline* spl);
void GeoComp_refine(Triangulation* DT, double r_ref, Spline* spl);

// delaunay triangulation functions
Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs, vector<double> &params);
Triangulation GeoComp_Delaunay_Triangulation(const vector<vector<int>> &segs, vector<vector<double>> &xs, vector<double> &params);
Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs);
void Flip_Insertion(Triangulation* DT, int* vid, int tri_s);
void Flip_Insertion_segment(Triangulation* DT, int vid, int hfid, Spline* spl);
void Bowyer_watson2d(Triangulation* DT, int vid, int tri_s,bool refine);
void delete_tris(Triangulation* DT);
void delete_tris(Triangulation* DT, int* tri);
bool check_sibhfs(Triangulation* DT);
bool check_jacobians(Triangulation* DT);
double check_minangle(Triangulation* DT);
vector<bool> find_boundary_nodes(Triangulation* DT);

// writeing to file functions
void WriteObj_mesh(const mesh &msh);
void WrtieVtk_tri(const mesh &msh);

void WrtieVtk_tri(const Triangulation &msh);
void WrtieVtk_tri(const Triangulation &msh, const vector<double> &data);
