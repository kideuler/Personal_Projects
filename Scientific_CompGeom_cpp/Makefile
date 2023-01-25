CC = g++
CFLAGS = -Iinclude -I ~/eigen/
DEPS = linalg.cpp SciComp.cpp Meshgen2d.cpp Meshgen3d.cpp MeshSmooth.cpp Splines.cpp
DEPS_O = linalg.o SciComp.o Meshgen2d.o Meshgen3d.o MeshSmooth.o Splines.o
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
