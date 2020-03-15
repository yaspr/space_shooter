gcc:
	gcc -std=c99 -g3 -Ofast flame.c main.c -o main -lm -lX11

icc:
	icc -std=c99 -Ofast flame.c main.c -o main -mkl -lX11

clean:
	rm -Rf *~ main
