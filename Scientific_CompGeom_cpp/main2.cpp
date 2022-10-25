#include <SciComp.hpp>

using namespace std;


vec a = {0.0, 2.0};

// redo derivative stuff in multiple dimensions as function<vec(vec)> f, vec x and
// function<double(vec), vec x>

vec func(vec x){vec y(2);
y[0] = pow(x[1]-1,2);
y[1] = pow(x[0]-2,2);
return y;}

vec func2(vec &x){return {sin(2*x[0])+cos(x[1]), x[0]*x[1]-sin(x[0])};}

double f(vec x){return sin(2*(x[0])) + cos(3*x[1]);};

int main(){
    Mat init;
    init.push_back({0,0});
    vec Jacobian(&func2,{0,0});
    //Mat Solve(func2,init);

    return 0;
}

