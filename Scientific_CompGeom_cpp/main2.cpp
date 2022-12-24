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
    Mat xs = rMat(500,2);
    cout << "created random points" << endl;
    Triangulation DT = GeoComp_Delaunay_Triangulation(xs);
    cout << "finished delaunay triangulation" << endl;
    mesh msh = GeoComp_DT2mesh(DT);
    WrtieVtk_tri(msh);
    cout << "finished writing to file" << endl;
}

