#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdint>
#include <cmath>
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

int Mandel(double x0, double y0);

const double xmin = -2.0;
const double xmax = 0.47;
const double ymin = -1.12;
const double ymax = 1.12;

int main() {
    constexpr auto dimx = 1000u, dimy = 1000u;

    double hx = (xmax - xmin)/ ( (double) dimx-1);
    double hy = (ymax - ymin)/ ( (double) dimy-1);

    using namespace std;
    ofstream ofs("Mandelbrot.ppm", ios_base::out | ios_base::binary);
    ofs << "P6" << endl << dimx << ' ' << dimy << endl << "255" << endl;

    double x = xmin;
    double y = ymin;
    double progress = 0.0;
    bool tf = false;
    int val;
    for (auto j = 0u; j < dimy; ++j){
        x = xmin;
        for (auto i = 0u; i < dimx; ++i){
            val = Mandel(x,y);
            ofs << (char) (val) << (char) (val) << (char) (val);       // red, green, blue
            x += hx;
        }    
        progress = (double) (j+1)/dimy;
        printProgress(progress);
        y += hy;
    }
    ofs.close();

    return EXIT_SUCCESS;
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
        iter = iter + 1;
    }

    return iter;
}