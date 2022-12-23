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

int main(){
    Mat xs = rMat(100,2);
    Triangulation DT = GeoComp_Delaunay_Triangulation(xs);
    mesh msh = GeoComp_DT2mesh(DT);
    printMat(msh.coords);
    WrtieVtk_tri(msh);
}

