#include <linalg.hpp>
#include <functional>
#include <cmath>

typedef vector<function<double(double)>> function1;
typedef vector<function<double(double,double)>> function2;
typedef vector<function<double(double,double,double)>> function3;

// derivative operations
double Dx(function<double(double)> f, double x);
double Dx(function<double(double, double)> f, double x, double y);
double Dx(function<double(double,double,double)> f, double x, double y, double z);

double Dy(function<double(double,double)> f, double x, double y);
double Dy(function<double(double,double,double)> f, double x, double y, double z);

double Dz(function<double(double,double,double)> f, double x, double y, double z);


// Gradients
vec Grad(function<double(double,double)> f, double x, double y);
vec Grad(function<double(double,double,double)> f, double x, double y, double z);

// Divergence
double Div(function<double(double,double)> f, double x, double y);
double Div(function<double(double,double,double)> f, double x, double y, double z);

// Curl
vec Curl(function<double(double,double,double)> f, double x, double y, double z);

// Second Derivatives
double DDx(function<double(double)> f, double x);
double DDx(function<double(double,double)> f, double x, double y);
double DDx(function<double(double,double,double)> f, double x, double y, double z);
double DDy(function<double(double,double)> f, double x, double y);
double DDy(function<double(double,double,double)> f, double x, double y, double z);
double DDz(function<double(double,double,double)> f, double x, double y, double z);
double DxDy(function<double(double,double)> f, double x, double y);
double DxDy(function<double(double,double,double)> f, double x, double y, double z);
double DxDz(function<double(double,double,double)> f, double x, double y, double z);
double DyDz(function<double(double,double,double)> f, double x, double y, double z);

// Laplacian
double Laplacian(function<double(double,double)> f, double x, double y);
double Laplacian(function<double(double,double,double)> f, double x, double y, double z);

// Hessian
Mat Hessian(function<double(double,double)> f, double x, double y);
Mat Hessian(function<double(double,double,double)> f, double x, double y, double z);

// Convection Diffusion Reaction Operator (cdr)
double cdr(function<double(double)> f, function<double(double)> v, function<double(double)> rho, double x);
double cdr(function<double(double,double)> f, function<vec(double,double)> v, function<double(double,double)> rho, double x, double y);
double cdr(function<double(double,double,double)> f, function<vec(double,double,double)> v, function<double(double,double,double)> rho, double x, double y, double z);


// Integration
double Integrate(function<double(double)> f, double a, double b, int nsteps=1);


// Equation solver
vec Solve(function<double(double)> f, int nsolutions, vec initial);
