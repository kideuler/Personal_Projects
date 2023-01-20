#include <GeoComp.hpp>
#include <cmath>
#include <chrono>

#define FAIL "\033[0;31mFAIL\033[0m"
#define PASS "\033[0;32mPASS\033[0m"

using namespace std;

int main(){
    Mat xs = {{0,0},{0,-1}};
    Mesh DT = GeoComp_Delaunay_Mesh3d(xs);
    WrtieVtk_tet(DT);
}