#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdint>
#include <cmath>
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
using namespace std;

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

int Mandel(double x0, double y0);
int Julia(double x0, double y0, double R, const double cx, const double cy);
void Julia_draw(const string imid, bool show);

// image resolution
const int dimx = 500;
const int dimy = 500;


double cx = -0.75;
double cy = 0.01;
double R = 0.5 + sqrt(1 + 4*sqrt(cx*cx +cy*cy))/2;

int main() {
    cy = 0.1;
    int iter = 1;
    while (cy<0.5){
        cy+=0.005;
        R = 0.5 + sqrt(1 + 4*sqrt(cx*cx +cy*cy))/2;
        Julia_draw(to_string(iter), false);
        iter++;
    }
    system(("convert -delay 5 -loop 0 Julia*.png Julia.gif"));
    system("rm -rf *.png");
    system("xdg-open Julia.gif");

    return EXIT_SUCCESS;
}

void Julia_draw(const string imid, bool show){
    double hx = (2*R)/ ( (double) dimx-1);
    double hy = (2*R)/ ( (double) dimy-1);

    string filename = "Julia"+ imid;

    ofstream ofs(filename+".ppm", ios_base::out | ios_base::binary);
    ofs << "P6" << endl << dimx << ' ' << dimy << endl << "255" << endl;

    double x = -R;
    double y = -R;
    double progress = 0.0;
    
    bool tf = false;
    int val;
    for (auto j = 0u; j < dimy; ++j){
        x = -R;
        for (auto i = 0u; i < dimx; ++i){
            val = Julia(x,y,R,cx,cy);
            ofs << (char) (val) << (char) (val) << (char) (val);       // red, green, blue
            x += hx;
        }    
        progress = (double) (j+1)/dimy;
        printProgress(progress);
        y += hy;
    }
    ofs.close();
    system(("pnmtopng " + filename+".ppm > "+ filename + ".png").c_str());
    system(("rm -rf " + filename+".ppm").c_str());
    if (show) {system(("xdg-open "+filename+".png").c_str());}
}


bool circle(double x, double y){
    return pow(x+0.765,2) + pow(y,2) <= 0.25;
}

int Mandel(double x0, double y0){
    int iter = 0;
    int max_iter = 255;
    double x = 0.0;
    double y = 0.0;
    double xtemp;
    while (x*x + y*y <= 4 && iter < max_iter){
        xtemp = x*x - y*y + x0;
        y = 2*x*y + y0;
        x = xtemp;
        iter++;
    }

    return iter;
}

int Julia(double x, double y, const double R, const double cx, const double cy){
    int iter = 0;
    int max_iter = 51;
    double xtemp;
    while (x*x + y*y <= R*R && iter <= max_iter){
        xtemp = x*x - y*y + cx;
        y = 2*x*y + cy;
        x = xtemp;
        iter++;
    }

    return iter*5;
}