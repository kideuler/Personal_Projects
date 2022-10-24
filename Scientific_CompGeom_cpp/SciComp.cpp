#include<SciComp.hpp>

#define eps 2.2e-16
using namespace std;


template<typename T>
struct args;

template<typename R, typename ...Args> 
struct args<function<R(Args...)>>
{
    static const size_t nargs = sizeof...(Args);
};


// derivitives of a function pointers
double Dx(function<double(double)> f, double x){
    double h = sqrt(eps)*x;
    if (abs(f(x)) > 1e-15){
        return (0.5*f(x+h) - 0.5*f(x-h))/h;
    } else {
        return (-1.5*f(x) + 2*f(x+h) - 0.5*f(x+2*h))/h;
    }
}
double Dx(function<double(double,double)> f, double x, double y){
    double h = sqrt(eps)*x;
    if (abs(f(x,y)) > 1e-15){
        return (0.5*f(x+h,y) - 0.5*f(x-h,y))/h;
    } else {
        return (-1.5*f(x,y) + 2*f(x+h,y) - 0.5*f(x+2*h,y))/h;
    }
}
double Dx(function<double(double,double,double)> f, double x, double y, double z){
    double h = sqrt(eps)*x;
    if (abs(f(x,y,z)) > 1e-15){
        return (0.5*f(x+h,y,z) - 0.5*f(x-h,y,z))/h;
    } else {
        return (-1.5*f(x,y,z) + 2*f(x+h,y,z) - 0.5*f(x+2*h,y,z))/h;
    }
}
double Dy(function<double(double,double)> f, double x, double y){
    double h = sqrt(eps)*y;
    if (abs(f(x,y)) > 1e-15){
        return (0.5*f(x,y+h) - 0.5*f(x,y-h))/h;
    } else {
        return (-1.5*f(x,y) + 2*f(x,y+h) - 0.5*f(x,y+2*h))/h;
    }
}
double Dy(function<double(double,double,double)> f, double x, double y, double z){
    double h = sqrt(eps)*y;
    if (abs(f(x,y,z)) > 1e-15){
        return (0.5*f(x,y+h,z) - 0.5*f(x,y-h,z))/h;
    } else {
        return (-1.5*f(x,y,z) + 2*f(x,y+h,z) - 0.5*f(x,y+2*h,z))/h;
    }
}
double Dz(function<double(double,double,double)> f, double x, double y, double z){
    double h = sqrt(eps)*z;
    if (abs(f(x,y,z)) > 1e-15){
        return (0.5*f(x,y,z+z) - 0.5*f(x,y,z-h))/h;
    } else {
        return (-1.5*f(x,y,z) + 2*f(x,y,z+h) - 0.5*f(x,y,z+2*h))/h;
    }
}

// Gradients
vec Grad(function<double(double,double)> f, double x, double y){
    vec G(2);
    G[0] = Dx(f,x,y);
    G[1] = Dy(f,x,y);
    return G;
}
vec Grad(function<double(double,double,double)> f, double x, double y, double z){
    vec G(3);
    G[0] = Dx(f,x,y,z);
    G[1] = Dy(f,x,y,z);
    G[2] = Dz(f,x,y,z);
    return G;
}

// divergence
double Div(function<double(double,double)> f, double x, double y){
    return Dx(f,x,y) + Dy(f,x,y);
}
double Div(function<double(double,double,double)> f, double x, double y, double z){
    return Dx(f,x,y,z) + Dy(f,x,y,z) + Dz(f,x,y,z);
}

// curl
vec Curl(function<double(double,double,double)> f, double x, double y, double z){
    vec C(3);
    C[0] = Dy(f,x,y,z)-Dz(f,x,y,z);
    C[1] = Dz(f,x,y,z)-Dx(f,x,y,z);
    C[2] = Dx(f,x,y,z)-Dy(f,x,y,z);
    return C;
}

// second derivatives
double DDx(function<double(double)> f, double x){
    double h = sqrt(eps)*x;
    return (f(x+h)-2*f(x)+f(x-h))/(h*h);
}
double DDx(function<double(double,double)> f, double x, double y){
    double h = sqrt(eps)*x;
    return (f(x+h,y)-2*f(x,y)+f(x-h,y))/(h*h);
}
double DDx(function<double(double,double,double)> f, double x, double y, double z){
    double h = sqrt(eps)*x;
    return (f(x+h,y,z)-2*f(x,y,z)+f(x-h,y,z))/(h*h);
}
double DDy(function<double(double,double)> f, double x, double y){
    double h = sqrt(eps)*y;
    return (f(x,y-h)-2*f(x,y)+f(x,y+h))/(h*h);
}
double DDy(function<double(double,double,double)> f, double x, double y, double z){
    double h = sqrt(eps)*y;
    return (f(x,y-h,z)-2*f(x,y,z)+f(x,y+h,z))/(h*h);
}
double DDz(function<double(double,double,double)> f, double x, double y, double z){
    double h = sqrt(eps)*z;
    return (f(x,y,z-h)-2*f(x,y,z)+f(x,y,z+h))/(h*h);
}

double DxDy(function<double(double,double)> f, double x, double y){
    double hx = sqrt(eps)*x;
    double hy = sqrt(eps)*y;
    return (f(x+hx,y+hy)-f(x+hx,y-hy)-f(x-hx,y+hy)+f(x-hx,y-hy))/(4*hx*hy);
}
double DxDy(function<double(double,double,double)> f, double x, double y, double z){
    double hx = sqrt(eps)*x;
    double hy = sqrt(eps)*y;
    return (f(x+hx,y+hy,z)-f(x+hx,y-hy,z)-f(x-hx,y+hy,z)+f(x-hx,y-hy,z))/(4*hx*hy);
}
double DxDz(function<double(double,double,double)> f, double x, double y, double z){
    double hx = sqrt(eps)*x;
    double hz = sqrt(eps)*z;
    return (f(x+hx,y,z+hz)-f(x+hx,y,z-hz)-f(x-hx,y,z+hz)+f(x-hx,y,z-hz))/(4*hx*hz);
}
double DyDz(function<double(double,double,double)> f, double x, double y, double z){
    double hy = sqrt(eps)*y;
    double hz = sqrt(eps)*z;
    return (f(x,y+hy,z+hz)-f(x,y+hy,z-hz)-f(x,y-hy,z+hz)+f(x,y-hy,z-hz))/(4*hy*hz);
}

// Laplacian
double Laplacian(function<double(double,double)> f, double x, double y){
    return DDx(f,x,y) + DDy(f,x,y);
}
double Laplacian(function<double(double,double,double)> f, double x, double y, double z){
    return DDx(f,x,y,z) + DDy(f,x,y,z) + DDz(f,x,y,z);
}

// Hessian Matrix
Mat Hessian(function<double(double,double)> f, double x, double y){
    Mat H = Zeros(2,2);
    H[0][0] = DDx(f,x,y);
    H[1][0] = DxDy(f,x,y);
    H[0][1] = H[1][0];
    H[1][1] = DDy(f,x,y);
    return H;
}
Mat Hessian(function<double(double,double,double)> f, double x, double y, double z){
    Mat H = Zeros(3,3);
    H[0][0] = DDx(f,x,y,z);
    H[1][0] = DxDy(f,x,y,z);
    H[0][1] = H[1][0];
    H[0][2] = DxDz(f,x,y,z);
    H[2][0] = H[0][2];
    H[1][1] = DDy(f,x,y,z);
    H[1][2] = DyDz(f,x,y,z);
    H[2][1] = H[1][2];
    H[2][2] = DDz(f,x,y,z);
    return H;
}

// CDR
double cdr(function<double(double)> f, function<double(double)> v, function<double(double)> rho, double x){
    return -DDx(f,x) + v(x)*Dx(f,x) + rho(x)*f(x);
}
double cdr(function<double(double,double)> f, function<vec(double,double)> v, function<double(double,double)> rho, double x, double y){
    return -Laplacian(f,x,y) + inner(v(x,y), Grad(f,x,y)) + rho(x,y)*f(x,y);
}
double cdr(function<double(double,double,double)> f, function<vec(double,double,double)> v, function<double(double,double,double)> rho, double x, double y, double z){
    return -Laplacian(f,x,y,z) + inner(v(x,y,z), Grad(f,x,y,z)) + rho(x,y,z)*f(x,y,z);
}


// Integration
double Integrate(function<double(double)> f, double a, double b, int nsteps){
    // using Simpsons formula
    double I = 0.0;
    double x = a;
    double h = (b-a)/((double) nsteps);
    for (int i = 0; i<nsteps; i++){
        I+= (f(x) + 4*f(x+h/2) + f(x+h));
        x += h;
    }
    I = I*(h/6);
    cout << x << " ";
    return I;
}




// Nonlinear equation solving
vec Solve(function<double(double)> f, int nsolutions, vec initial){
    int max_iters = 1000;
    int nsegs, i, n, iter;
    double x,x0,x_1,x_new,h,Tol,D;

    nsegs = 100;
    vec solutions(0);
    for (int n = 0; n<nsolutions; n++){

        // creating new function by eliminating old solutions
        auto F = [solutions,f](double x){
            double denom = 1.0;
            for(int i = 0; i<solutions.size();i++){
                denom = denom*(x - solutions[i]);
            }
            return f(x)/denom;
        };

        // use search to find rough estimate within bounds        
        h = (initial[1]-initial[0])/((double) nsegs);
        x0 = initial[0];
        x = initial[0];
        for (int i = 0; i<nsegs; i++){
            if (abs(F(x)) < abs(F(x0))){
                x0 = x;
            }
            x = x + h;
        }


        // Finding solution using Newtons Method (Secant method in rare cases)
        iter = 0;
        Tol = sqrt(eps)*abs(x0);
        x = x0;
        x_1 = x + h;
        while (abs(F(x)) >= Tol && iter < max_iters){
            D = Dx(F,x);
            if (abs(D) > Tol){
                x_1 = x;
                x = x - F(x)/Dx(F,x);
            } else {
                x_new = (x_1*F(x) - x*F(x_1))/(F(x) - F(x_1));
                x_1 = x;
                x = x_new;
            }           
            iter++;
            Tol = sqrt(eps)*abs(x);
        }

        if (iter == max_iters){
            cout << "SciComp.cpp::Solve, No solution found for sol " << n+1 << endl;
        } else {
            solutions.push_back(x);
        }
    }
    return solutions;
}


// Rootfinding for multivariate problems Newtons / Broydens

// 1d opt Golden Search

// Multivariate Optimization (Steepest Descent / Newtons method / Conjugate Gradient / BFGS)

// Interpolation

