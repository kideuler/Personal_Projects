#include <GeoComp.hpp>
#include <cmath>

#define FAIL "\033[0;31mFAIL\033[0m"
#define PASS "\033[0;32mPASS\033[0m"

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

Mat Star(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.5*pow(cos(t),3)+0.5;
        xs[i][1] = 0.5*pow(sin(t),3)+0.5;
    }
    return xs;
}

Mat Cardiod(int npoints){
    Mat xs = Zeros(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs[i][0] = 0.5*(2*cos(t) - cos(2*t))+0.5;
        xs[i][1] = 0.5*(2*sin(t) - sin(2*t))+0.5;
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
        return 0.4*(sqrt(3)/3)*M_PI/(99);
    } else{
        return (sqrt(3)/3)*M_PI/(99);
    }
}
/**
 * STILL TO DO:
 *  - Add Test cases from CAD from scratch and get them working
 *  - Bucket sort point data for optimal insertion in intial DT (SLOAN PAPER)
 *  - Add constriained delaunay triangulation support
 *  - Add energy based nodal smoothing
 */

int main(){
    Mat xs = {{0,7},{-5,5},{5,5},{-2,3},{3,1},{-4,-1},{1,-2},{-6,-4},{5,-4}};
    Triangulation DT = GeoComp_Delaunay_Triangulation(xs);
    check_jacobians(&DT);

    xs = {{1,1},{3,4},{-2,3},{-2,2},{-1,-1},{-2,-3},{4,-2}};
    DT = GeoComp_Delaunay_Triangulation(xs);
    check_jacobians(&DT);

    xs = {{0,7},{-5,5},{5,5},{-2,3},{3,1},{-4,-1},{1,-2},{-6,-4},{5,-4}};
    vector<vector<int>> segs = {{4,6}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);
    check_jacobians(&DT);
    
    segs = {{4,5}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);
    check_jacobians(&DT);
    

    segs = {{6,0}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);
    

    segs = {{2,7}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);
    WrtieVtk_tri(DT);

    segs = {{1,8}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);

    segs = {{4,5},{2,3}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);

    segs = {{1,2},{2,5}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);

    segs = {{7,4},{4,8}};
    DT = GeoComp_Delaunay_Triangulation(segs, xs);
    

    /*
    int n = 50;
    xs  = Flower(n);
    cout << "created points" << endl;
    segs = Zerosi(n,2);
    double h = 0.0;
    for (int i = 0; i<n; i++){
        segs[i][0] = i;
        segs[i][1] = (i+1)%n;
        h = h + norm(xs[(i+1)%n]-xs[i]);
    }
    h = 1.5*(sqrt(3)/3)*(h/(double(n)));

    DT = GeoComp_Delaunay_Triangulation(segs,xs);
    
    cout << h << endl;
    GeoComp_refine(&DT, h);
    cout << "finished delaunay triangulation" << endl;
    WrtieVtk_tri(DT);
    cout << "finished writing to file" << endl;
    */
}