#include <GeoComp.hpp>
#include <cmath>

using namespace std;

vec uexpr(Mat xs){
    int nv = xs.size();
    vec u(nv);
    for (int i = 0; i<nv; i++){
        u[i] = xs[i][0]*xs[i][1];
    }
    return u;
}

Mat Circle(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.5*cos(t)+0.5;
        xs[i][1] = 0.5*sin(t)+0.5;
    }
    return xs;
}

Mat Ellipse(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.5*cos(t)+0.5;
        xs[i][1] = 0.25*sin(t)+0.5;
    }
    return xs;
}

Mat Flower(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = (0.25 + 0.1*sin(5*t))*cos(t)+0.5;
        xs[i][1] = (0.25 + 0.1*sin(5*t))*sin(t)+0.5;
    }
    return xs;
}

Mat Box(int npoints){
    Mat xs = Zeros(npoints*npoints,2);
    int k = 0;
    for (int i =0;i<npoints;i++){
        for (int j=0; j<npoints;j++){
            xs[k][0] = ((double) j)/((double) npoints-1);
            xs[k][1] = ((double) i)/((double) npoints-1);
            k++;
        }
    }
    return xs;
}

double Gradiate(vec xs){
    if ((pow(xs[0]-0.5,2)+pow(xs[1]-0.5,2))< .1*.1 || abs(xs[0]-xs[1]) < 0.05){
        return 0.6*(sqrt(3)/3)*M_PI/(99);
    } else{
        return (sqrt(3)/3)*M_PI/(99);
    }
}

int main(){
    int n = 7;
    Mat xs = Box(n);
    cout << "created points" << endl;
    vector<vector<int>> segs = Zerosi(n,2);
    double h = 0.0;
    for (int i = 0; i<n; i++){
        segs[i][0] = i;
        segs[i][1] = (i+1)%n;
        h = h + norm(xs[(i+1)%n]-xs[i]);
    }
    h = (sqrt(3)/3)*(h/(double(n)));
    Triangulation DT = GeoComp_Delaunay_Triangulation(xs);
    cout << h << endl;
    //GeoComp_refine(&DT, h);
    cout << "finished delaunay triangulation" << endl;
    WrtieVtk_tri(DT);
    cout << "finished writing to file" << endl;
}