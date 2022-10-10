CC = g++
CFLAGS = -Iinclude
DEPS = linalg.cpp
OBJ = main

make:
	$(CC) $(CFLAGS) -c -Wall -g  $(DEPS)
	$(CC) $(CFLAGS)  -g main.cpp linalg.o -o $(OBJ)
	./$(OBJ)

.PHONY: clean
clean:
	rm -rf *.o