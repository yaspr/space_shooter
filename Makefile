all:
	gcc -std=c99 -g3 -Ofast flame.c main.c -o main -lm -lX11

clean:
	rm -Rf *~ main
