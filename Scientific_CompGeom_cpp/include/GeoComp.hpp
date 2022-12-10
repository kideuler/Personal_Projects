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
    vector<vector<int>> elems;
    vector<vector<int>> sibhfs;
    int nelems;
    vector<int> badtris;
    vector<int> vedge;
    vector<bool> deletes;
    vector<vector<int>> facets;
    vector<bool> bwork;
    vector<bool> on_boundary;
};

// encoding and decoding for sibhfs array
int elids2hfid(int eid, int lid);
int hfid2eid(int hfid);
int hfid2lid(int hlid);

Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs);
void Bowyer_watson_insert_point2d(Triangulation &DT, int vid, int starting_tri);
void find_circumtri(Triangulation &DT, int tri, int &nbad, int vid);
void reorder(vector<vector<double>> &xs);
bool inside_tri(const Triangulation &DT, int tri, int vertex);
double detv(const vector<double> &u, const vector<double> &v);


void WriteStl(const Triangulation &msh);
void WrtieVtk_tri(const Triangulation &msh);
void WrtieVtk_tri(const Triangulation &msh, const vector<double> &data);