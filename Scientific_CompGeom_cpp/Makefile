CC = g++
CFLAGS = -Iinclude -I ~/eigen/
DEPS = linalg.cpp Meshgen2d.cpp Meshgen3d.cpp MeshSmooth.cpp Splines.cpp Blossom.cpp Meshutils.cpp SurfaceMesh.cpp
DEPS_O = linalg.o Meshgen2d.o Meshgen3d.o MeshSmooth.o Splines.o Blossom.o Meshutils.o SurfaceMesh.o
OBJ = main

make:
	$(CC) $(CFLAGS) -c -g  $(DEPS)

.PHONY: clean
clean:
	rm -rf *.o

.PHONY: run
run:
	$(CC) $(CFLAGS)  -g main2.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ)
	sudo cp test.vtk /mnt/c/Users/Jacob/Documents/test.vtk
	sudo cp bunny.obj /mnt/c/Users/Jacob/Documents/bunny.obj
