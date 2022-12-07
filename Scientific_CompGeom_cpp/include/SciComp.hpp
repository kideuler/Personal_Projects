#include <linalg.hpp>
#include <functional>
#include <cmath>
#include <cstring>

typedef vector<function<double(double)>> function1;
typedef vector<function<double(double,double)>> function2;
typedef vector<function<double(double,double,double)>> function3;

// derivative operations
double Dx(function<double(double)> f, double x);
double Dxi(function<double(vec)> f, vec x, int var);
vec Dxi(function<vec(vec)> f, vec x, int var);

// Gradients
vec Grad(function<double(vec)> f, vec x);
Mat Jacobian(function<vec(vec)> f, vec x);

// Divergence
double Div(function<double(vec)> f, vec x);

// Curl
vec Curl(function<double(vec)> f, vec x);

// Second Derivatives
double DDxi(function<double(vec)> f, vec x, int var);
double DxiDxj(function<double(vec)> f, vec x, int var1, int var2);

// Laplacian
double Laplacian(function<double(vec)> f, vec x);

// Hessian
Mat Hessian(function<double(vec)> f, vec x);

// Convection Diffusion Reaction Operator (cdr)
double cdr(function<double(vec)> f, function<vec(vec)> v, function<double(vec)> rho, vec x);

// Integration
double Integrate(function<double(double)> f, double a, double b, int nsteps=1);


// Equation solver
vec Solve(function<double(double)> f, int nsolutions, vec initial);
Mat Solve(function<vec(vec)> f, Mat initial);

// 1d optimization
double Minimize(function<double(double)> f, vec bounds, double Tol = 1e-10);

// Multidimensional Optimization
vec Minimize(function<double(vec)> f, vec x_initial, int method, double Tol = 1e-6);

// Interpolation functions
Mat gen_vander_1d(Mat us, int degree);
Mat gen_vander_2d(Mat us, int degree);
Mat Tri_natcoords(int degree);
Mat Line_natcoords(int degree);
