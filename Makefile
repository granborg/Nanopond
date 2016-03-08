npx: nanopond2.c
	gcc --verbose 									\
		-g -Wall									\
		-O3 -msse2 nanopond2.c -o npx					\
		-lSDL -fopenmp

clean:
	rm -f ./npx

test:  npx
	/usr/bin/time ./npx
		

