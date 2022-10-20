CC = g++
CFLAGS = -Iinclude
DEPS = Mandelbrot.cpp
OBJ = main

make:
	$(CC) $(CFLAGS)  -g Mandelbrot.cpp -o $(OBJ)
	./$(OBJ)
	pnmtopng Mandelbrot.ppm > Mandelbrot.png

.PHONY: clean
clean:
	rm -rf *.o