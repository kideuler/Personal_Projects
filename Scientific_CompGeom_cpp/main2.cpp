#include <GeoComp.hpp>

using namespace std;


int main( ){
    Mat xs = rMat(10,2);
    printMat(xs);
    cout << endl;
    Triangulation DT =  GeoComp_Delaunay_Triangulation(xs);
}

