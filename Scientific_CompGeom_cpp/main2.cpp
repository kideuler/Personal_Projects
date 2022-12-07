#include <SciComp.hpp>

using namespace std;



// redo derivative stuff in multiple dimensions as function<vec(vec)> f, vec x and
// function<double(vec), vec x>



int main(){
    Mat us = Tri_natcoords(2);
    printMat(us);
    Mat V = gen_vander_2d(us,2);
    printMat(V);
}

