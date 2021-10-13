.SUFFIXES:
FC = gfortran
FLAGS = -Wall -fcheck=all
main: 
	$(FC) -g -L/usr/local/lib/ init_mesh.f90 Mesh_smoothing.f90 marching_cubes.f90 surface.f90 utilities.f90 main.f90 -llapack -lblas
.PHONY: clean
clean:
	-rm -f *.o *.mod *.out
.PHONY: run
run:
	./a.out