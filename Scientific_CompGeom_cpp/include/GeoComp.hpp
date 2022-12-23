#include <linalg.hpp>
#include <fstream>
#include <algorithm>

struct Elems
{
    vector<int> elem;
    struct Elems *next;
    struct Elems *prev;
};

struct mesh {
    vector<vector<double>> coords;
    vector<vector<int>> elems;
    vector<vector<int>> sibhfs;
};

struct Triangulation {
    vector<vector<double>> coords;
    Elems* elem_head;
    Elems* sibhfs_head;
    int nelems;
    vector<vector<int>> facets;
    vector<bool> bwork;
    vector<bool> on_boundary;
};

void push_elem(Elems** head, vector<int> elem);
Elems* pop_elem(Elems* head);
void insert_elem(Elems* prev_elem, vector<int> elem);
void delete_elem(Elems** head, Elems** elem);
void print_LLElems(Elems* head);


Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs);
void Bowyer_watson_insert_point2d(Triangulation &DT, int vid);
mesh GeoComp_DT2mesh(Triangulation &DT);
void WrtieVtk_tri(const mesh &msh);
/*
static void WriteStl(const Triangulation &msh);
void WrtieVtk_tri(const Triangulation &msh);
void WrtieVtk_tri(const Triangulation &msh, const vector<double> &data);
*/