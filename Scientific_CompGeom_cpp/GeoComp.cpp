#include "GeoComp.hpp"
using namespace std;

int elids2hfid(int eid, int lid){
    return (eid << 8) + lid;
}

int hfid2eid(int hfid){
    return hfid >> 8;
}

int hfid2lid(int hfid){
    return (hfid & 255);
}

vector<double> min_array(const vector<vector<double>> &xs){
    int ndims = xs[0].size();
    vector<double> min_(ndims);

    for (int j = 0; j<ndims; j++){
        min_[j] = xs[0][j];
    }

    for (int i = 1; i<xs.size(); i++){
        for (int j = 0; j<ndims; j++){
            if (xs[i][j] < min_[j]){
                min_[j] = xs[i][j];
            }
        }
    }

    return min_;
}

vector<double> max_array(const vector<vector<double>> &xs){
    int ndims = xs[0].size();
    vector<double> min_(ndims);

    for (int j = 0; j<ndims; j++){
        min_[j] = xs[0][j];
    }

    for (int i = 1; i<xs.size(); i++){
        for (int j = 0; j<ndims; j++){
            if (xs[i][j] > min_[j]){
                min_[j] = xs[i][j];
            }
        }
    }

    return min_;
}


Triangulation GeoComp_Delaunay_Triangulation(const vector<vector<double>> &xs){
    // size checking
    int nv = xs.size();
    int ndims = xs[0].size();
    assert(ndims == 2);

    // establishing upper bound
    int upper_bound = 2*nv;

    Triangulation DT;
    DT.coords = Zeros(nv+3,2);
    DT.elems.reserve(upper_bound);
    DT.sibhfs.reserve(upper_bound);

    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);

    DT.coords[nv][0] = a[0];
    DT.coords[nv][1] = a[1];
    DT.coords[nv+1][0] = a[0] + 2*(b[0]-a[0]);
    DT.coords[nv+1][1] = a[1];
    DT.coords[nv+2][0] = a[0];
    DT.coords[nv+2][1] = a[1] + 2*(b[1]-a[1]);

    DT.elems.push_back({nv,nv+1,nv+2});
    DT.nelems = 1;

    printMat(DT.coords);
    printMat(DT.elems);
    

    return DT;
}