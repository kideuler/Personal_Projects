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

double Dxi(function<double(vec)> f, vec x, int var){
    double h = sqrt(eps)*x[var];
    return (0.5*f(add_i(x,h,var)) - 0.5*f(add_i(x,-h,var)))/h;
}
vec Dxi(function<vec(vec)> f, vec x, int var){
    double h = sqrt(eps)*x[var];
    return (0.5*f(add_i(x,h,var)) - 0.5*f(add_i(x,-h,var)))/h;
}

// Gradients
vec Grad(function<double(vec)> f, vec x){
    int n = x.size();
    vec G(n);
    for (int i = 0; i<n; i++){
        G[i] = Dxi(f,x,i);
    }
    return G;
}
Mat Jacobian(function<vec(vec)> f, vec x){
    int m = f(x).size();
    int n = x.size();
    vec G(n);
    Mat J = Zeros(m,n);
    for (int i = 0; i<n; i++){
        G = Dxi(f,x,i);
        for (int j = 0; j<m; j++){
            J[j][i] = G[j];
        }
    }
    return J;
}

// divergence
double Div(function<double(vec)> f, vec x){
    double D = 0.0;
    for (int i = 0; i<x.size(); i++){
        D += Dxi(f,x,i);
    }
    return D;
}

// curl
vec Curl(function<double(vec)> f, vec x){
    int n = x.size();
    assert(n == 3);
    vec G = Grad(f,x);
    return {G[2]-G[3],G[3]-G[1],G[1]-G[2]};
}


// second derivatives
double DDxi(function<double(vec)> f, vec x, int var){
    double h = sqrt(eps)*x[var];
    return (f(add_i(x,h,var)) - 2*f(x) + f(add_i(x,-h,var)))/(h*h);
}
double DxiDxj(function<double(vec)> f, vec x, int var1, int var2){
    double h1 = sqrt(eps)*x[var1];
    double h2 = sqrt(eps)*x[var2];
    return (f(add_i(add_i(x,h1,var1),h2,var2)) - \
    f(add_i(add_i(x,h1,var1),-h2,var2)) - \ 
    f(add_i(add_i(x,-h1,var1),h2,var2)) + \
    f(add_i(add_i(x,-h1,var1),-h2,var2)))/(4*h1*h2);
}

// Laplacian
double Laplacian(function<double(vec)> f, vec x){
    double L = 0.0;
    for (int i = 0; i<x.size(); i++){
        L += DDxi(f,x,i);
    }
    return L;
}

// Hessian Matrix
Mat Hessian(function<double(vec)> f, vec x){
    int n = x.size();
    Mat H = Zeros(n,n);
    int i,j;
    for (i = 0; i<n; i++){
        for (j = 0; j<n; j++){
            if (i==j){
                H[i][i] = DDxi(f,x,i);
            } else{
                H[i][j] = DxiDxj(f,x,i,j);
            }
        }
    }
    return H;
}

// CDR
double cdr(function<double(vec)> f, function<vec(vec)> v, function<double(vec)> rho, vec x){
    return -Laplacian(f,x) + inner(v(x),Grad(f,x)) + rho(x)*f(x);
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

    nsegs = (int) ceil(pow(1000.0,1/( (double) initial.size())));
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
Mat Solve(function<vec(vec)> f, Mat initial){
    // f: vector valued function
    // initial: initial starting values which give n solutions (distinct)
    int max_iters  = 1000;
    int nsegs, i,n, iter;
    double Tol;
    bool solves;
    vec x,x_1,x_new,s,b;
    Mat J;

    Mat solutions;
    function<vec(vec)> F = f;
    for (n = 0; n<initial.size(); n++){
        // creating new function by eliminating old solutions
        if (n > 0) {
            function<vec(vec)> F = [solutions,f](vec x){
                double denom = 1.0;
                for(int i = 0; i<solutions.size();i++){
                    for (int j = 0; j<solutions[i].size();j++){
                        denom = denom*(x[j] - solutions[i][j]);
                    }
                }
                return f(x)/denom;
            };
        }

        // main algorithm
        iter = 0;
        Tol = sqrt(eps)*norm(row(initial,n));
        x = row(initial,n);
        x_1 = x;
        while (norm(F(x)) >= Tol && iter < max_iters){
            J = Jacobian(F,x);
            b = F(x);
            s = QR_solve(J,-b,solves);
            if (solves){
                // use Newtons Method
                x_1 = x;
                x = x + s;
            } else {
                cout << "bad Jacobian" << endl;
            }           
            iter++;
            Tol = sqrt(eps)*norm(x);
        }

        // add solutions to solution matrix
        if (iter == max_iters){
            cout << "SciComp.cpp::Solve, No solution found for sol " << n+1 << endl;
        } else {
            solutions.push_back(x);
        }
    }
    return solutions;
}

// 1d opt Golden Search

// Multivariate Optimization (Steepest Descent / Newtons method / Conjugate Gradient / BFGS)

// Interpolation

