#include <linalg.hpp>
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
    vector<bool> delete_elems;
    vector<vector<int>> facets;
    vector<bool> bwork;
    vector<bool> on_boundary;
};

// encoding and decoding for sibhfs array
int elids2hfid(int eid, int lid);
int hfid2eid(int hfid);
int hfid2lid(int hlid);

Triangulation GeoComp_Delaunay_Triangulation(const vector<vector<double>> &xs);