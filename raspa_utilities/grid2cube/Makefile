default: grid2cube.x

grid2cube.o: grid2cube.c 
	gcc -c grid2cube.c -o grid2cube.o

grid2cube.x: grid2cube.o
	gcc grid2cube.o -o grid2cube.x -lm

clean:
	-rm -f grid2cube.o
	-rm -f grid2cube.x
