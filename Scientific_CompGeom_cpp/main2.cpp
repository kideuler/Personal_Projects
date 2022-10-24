#include <SciComp.hpp>

using namespace std;


vec a = {0.0, 2.0};

// redo derivative stuff in multiple dimensions as function<vec(vec)> f, vec x and
// function<double(vec), vec x>

double func(vec x){return sin(x[0])+pow(x[1],4);}

vec func2(vec x){return {sin(2*x[0])+cos(x[1]), x[0]*x[1]-sin(x[0])};}

int main(){
    printMat(Hessian(func,rvec(2)));
    printMat(Jacobian(func2,rvec(2)));

    return 0;
}

