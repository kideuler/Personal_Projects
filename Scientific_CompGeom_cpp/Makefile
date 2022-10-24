CC = g++
CFLAGS = -Iinclude
DEPS = linalg.cpp SciComp.cpp
DEPS_O = linalg.o SciComp.o
OBJ = main

make:
	$(CC) $(CFLAGS) -c -g  $(DEPS)
	$(CC) $(CFLAGS)  -g main2.cpp $(DEPS_O) -o $(OBJ)
	./$(OBJ)

.PHONY: clean
clean:
	rm -rf *.o