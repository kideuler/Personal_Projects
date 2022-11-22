#include <SciComp.hpp>

using namespace std;



// redo derivative stuff in multiple dimensions as function<vec(vec)> f, vec x and
// function<double(vec), vec x>


double f(vec x) {return pow(x[0]-1,2) + pow(x[1]-2,2);}

double a = 1.0;
double b = 100.0;
double Rosenbrock(vec x) {return (a - x[0]*x[0]) + b*pow(x[1]-x[0]*x[0],2);}
double Matyas(vec x) {return 0.26*(x[0]*x[0] + x[1]*x[1]) - 0.48*x[0]*x[1];}
double Beale(vec x) {return pow((1.5 - x[0] + x[0]*x[1]),2) + pow((2.25 - x[0] + x[0]*x[1]*x[1]),2) + pow(2.625 - x[0] + x[0]*x[1]*x[1]*x[1],2);}

int main(){
    vec x = Minimize(Rosenbrock, {2.50,1}, 1, 1e-6);
    cout << "minimum occurs at: " << endl;
    printvec(x);
    return 0;
}

