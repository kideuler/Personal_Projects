#include <GeoComp.hpp>

using namespace std;

vec uexpr(Mat xs){
    int nv = xs.size();
    vec u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = xs[i][0]*xs[i][1];
    }
    return u;
}

int main( ){
    Mat xs = rMat(10,2);
    Triangulation DT =  GeoComp_Delaunay_Triangulation(xs);
    //vec u = uexpr(DT.coords);
    printMat(DT.coords);
    //WrtieVtk_tri(DT,u);
}

