#include "Mandelbrot.hpp"


int main() {
    constexpr auto dimx = 800u, dimy = 800u;

    using namespace std;
    ofstream ofs("first.ppm", ios_base::out | ios_base::binary);
    ofs << "P6" << endl << dimx << ' ' << dimy << endl << "255" << endl;

    for (auto j = 0u; j < dimy; ++j)
        for (auto i = 0u; i < dimx; ++i)
            ofs << (char) (i % 256) << (char) (j % 256) << (char) ((i * j) % 256);       // red, green, blue

    ofs.close();

    return EXIT_SUCCESS;
}