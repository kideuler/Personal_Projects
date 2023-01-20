#include <linalg.hpp>
#include <fstream>
#include <algorithm>

struct Mesh {
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
Spline spline_init(const vector<vector<double>> &xs, const vector<bool> &corners);
vector<double> spline_var(Spline* spl, double t, int order=0);
vector<double> spline_point_segment(Spline* spl, double a, double b, double ratio);

// mesh smoothing functions
void mesh_smoothing_2d(Mesh* mesh, vector<bool> no_move, function<double(vector<double>)> r_ref, double mu);

// mesh refinement functions
void GeoComp_refine(Mesh* DT, function<double(vector<double>)> r_ref, Spline* spl);
void GeoComp_refine(Mesh* DT, double r_ref, Spline* spl);

// delaunay Mesh functions 2d
Mesh GeoComp_Delaunay_Mesh(vector<vector<double>> &xs, vector<double> &params);
Mesh GeoComp_Delaunay_Mesh(const vector<vector<int>> &segs, vector<vector<double>> &xs, vector<double> &params);
Mesh GeoComp_Delaunay_Mesh(vector<vector<double>> &xs);
void Flip_Insertion(Mesh* DT, int* vid, int tri_s);
void Flip_Insertion_segment(Mesh* DT, int vid, int hfid, Spline* spl);
void Bowyer_watson2d(Mesh* DT, int vid, int tri_s,bool refine);
void delete_tris(Mesh* DT);
void delete_tris(Mesh* DT, int* tri);
bool check_sibhfs(Mesh* DT);
bool check_jacobians(Mesh* DT);
double check_minangle(Mesh* DT);
vector<bool> find_boundary_nodes(Mesh* DT);

// delaunay Mesh functions 3d
Mesh GeoComp_Delaunay_Mesh3d(vector<vector<double>> &xs);

// writeing to file functions
void WrtieVtk_tri(const Mesh &msh);
void WrtieVtk_tet(const Mesh &msh);
void WrtieVtk_tri(const Mesh &msh, const vector<double> &data);

// small functions to be used in multiple files
vector<double> min_array(const vector<vector<double>> &xs);
vector<double> max_array(const vector<vector<double>> &xs);
int hfid2eid(int hfid);
int hfid2lid(int hfid);
int elids2hfid(int eid, int lid);

struct stack
{
    int hfid;
    struct stack *next;
};
void push_stack(stack** head, int hfid);
void pop_stack(stack** head);


