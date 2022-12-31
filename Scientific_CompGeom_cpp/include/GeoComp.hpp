#include <linalg.hpp>
#include <fstream>
#include <algorithm>

struct Elems
{
    vector<int> elem;
    struct Elems *next;
    struct Elems *prev;
    double alpha;
};

struct mesh {
    vector<vector<double>> coords;
    vector<vector<int>> elems;
    vector<vector<int>> sibhfs;
};

struct Triangulation {
    vector<vector<double>> coords;
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

void push_elem(Elems** head, vector<int> elem);
Elems* pop_elem(Elems* head);
void insert_elem(Elems* prev_elem, vector<int> elem);
void delete_elem(Elems** head, Elems** elem);
void print_LLElems(Elems* head);

void GeoComp_refine(Triangulation* DT, function<double(vector<double>)> r_ref);
void GeoComp_refine(Triangulation* DT, double r_ref);
Triangulation GeoComp_Delaunay_Triangulation(const vector<vector<int>> &segs, vector<vector<double>> &xs);
Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs);
void Bowyer_watson2d(Triangulation* DT, int vid, int tri_s,bool refine);
void delete_tris(Triangulation* DT);
void delete_tris(Triangulation* DT, int* tri);
bool check_sibhfs(Triangulation* DT);
bool check_jacobians(Triangulation* DT);

void WriteObj_mesh(const mesh &msh);
void WrtieVtk_tri(const mesh &msh);

void WrtieVtk_tri(const Triangulation &msh);
void WrtieVtk_tri(const Triangulation &msh, const vector<double> &data);
