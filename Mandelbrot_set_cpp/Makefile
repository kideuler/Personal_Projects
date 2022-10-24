CC = g++
CFLAGS = -Iinclude
DEPS = Mandelbrot.cpp
OBJ = main

make:
	rm -rf *.png *.ppm
	$(CC) $(CFLAGS)  -g Mandelbrot.cpp -o $(OBJ)
	./$(OBJ)
.PHONY: clean
clean:
	rm -rf *.o