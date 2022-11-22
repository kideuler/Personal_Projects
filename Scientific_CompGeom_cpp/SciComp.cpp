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
    double h = sqrt(eps)*max(abs(x[var]),1e-4);
    if (abs(f(x)) > 1e-15){
        return (0.5*f(add_i(x,h,var)) - 0.5*f(add_i(x,-h,var)))/h;
    } else {
        return (-1.5*f(x) + 2*f(add_i(x,h,var)) - 0.5*f(add_i(x,2*h,var)))/h;
    }
}
vec Dxi(function<vec(vec)> f, vec x, int var){
    double h = sqrt(eps)*max(abs(x[var]),1e-4);
    vec D = (0.5*f(add_i(x,h,var)) - 0.5*f(add_i(x,-h,var)))/h;

    if ( isnan(D[0]) ){
        D = (-1.5*f(x) + 2*f(add_i(x,h,var)) - 0.5*f(add_i(x,2*h,var)))/h;
    } 

    return D;
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


// Rootfinding for multivariate problems Newtons / Broydens (kinda works)
Mat Solve(function<vec(vec)> f, Mat initial){
    // f: vector valued function
    // initial: initial starting values which give n solutions (distinct)
    int max_iters  = 50;
    int nsegs, i,n, iter;
    double Tol;
    bool solves;
    vec x,x_1,x_new,s,b,y;
    Mat J,B;

    Mat solutions;
    function<vec(vec)> F = f;
    for (n = 0; n<initial.size(); n++){
        // creating new function by eliminating old solutions
        if (n > 0) {
            function<vec(vec)> F = [solutions,f](vec x){
                vec fx = f(x);
                for(int i = 0; i<solutions.size();i++){
                    for (int j = 0; j<solutions[i].size();j++){
                        fx[j] = fx[j]/(x[j] - solutions[i][j]);
                    }
                }
                return fx;
            };
        }

        // main algorithm
        iter = 0;
        Tol = sqrt(eps)*norm(row(initial,n));
        x = row(initial,n);
        printvec(x);
        x_1 = x;
        B = Jacobian(F,x+Tol);
        while (norm(F(x)) >= Tol && iter < max_iters){
            J = Jacobian(F,x);
            b = F(x);
            s = QR_solve(J,-b,solves);
            if (solves){
                // use Newtons Method
                x_1 = x;
                x = x + s;
                B = J;
                printvec(s);
            } else {
                // Broydens method
                s = QR_solve(B,-b,solves);
                x = x + s;
                y = F(x) - b;
                B = B + (outer(y-B*s,s))/(inner(s,s));
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

// 1d optimization Golden Search
double Minimize(function<double(double)> f, vec bounds, double Tol){
    double phi = (sqrt(5)-1)/2;
    double a = bounds[0];
    double b = bounds[1];
    double x1 = a + (1-phi)*(b-a);
    double f1 = f(x1);
    double x2 = a + phi*(b-a);
    double f2 = f(x2);
    while ((b-a) > Tol){
        if (f1 > f2){
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi*(b-a);
            f2 = f(x2);
        } else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1-phi)*(b-a);
            f1 = f(x1);
        }
    }
    return ((x1+x2)/2);
}

// Multivariate Optimization (Steepest Descent / Newtons method / Conjugate Gradient / BFGS)
// use Rosenbrock / Ackley function for testing
vec Minimize(function<double(vec)> f, vec x_initial, int method, double Tol){
    
    if (method  == 1){ // Use Gradient Descent
        double beta = 0.8;
        vec s = Grad(f,x_initial);
        double alpha_max;
        vec x = x_initial;
        vec x_1 = s;
        double alpha;
        while (norm(x-x_1) > Tol){
            s = Grad(f,x);
            s = s/(norm(s));
            alpha = 1.0;
            while (f(x-alpha*s) >= f(x) + alpha* 1e-4 *inner(s,s) && inner(s,Grad(f,x-alpha*s)) < f(x) + 0.325*inner(s,s)){
                alpha = beta*alpha;
            }
            x_1 = x;
            x = x - alpha*s;
            printvec(x);
        }
        return x;
    }

    return vec(0);
}

// Interpolation

