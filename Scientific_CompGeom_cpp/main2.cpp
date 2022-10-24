#include <SciComp.hpp>

using namespace std;


vec a = {0.0, 2.0};

// redo derivative stuff in multiple dimensions as function<vec(vec)> f, vec x

double func(double x){return sin(2*x)+cos(x);}

int main(){
    vec init(2); 
    init[0] = -5;
    init[1] = 10;
    vec sol = Solve(func,5,init);
    cout << " " << endl;
    printvec(sol);
}

