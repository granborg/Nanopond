npx: nanopond.c
	gcc --verbose 									\
		-g -Wall									\
		-O3 -msse2 nanopond.c -o npx					\
		-fopenmp

clean:
	rm -f ./npx

test:  npx
	/usr/bin/time ./npx
		

