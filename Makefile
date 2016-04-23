npo: nanopond.c
	gcc --verbose 					\
		-g -Wall				\
		-O3 -msse2 nanopond.c -o npo		\
		-lSDL -fopenmp

npsoa: nanopondSOA.c
	gcc --verbose					\
		-g -Wall				\
		-O3 -msse2 nanopondSOA.c -o npsoa	\
		-lSDL2 -fopenmp


clean:
	rm -f ./npo ./npn

testold:  npx
	/usr/bin/time ./npo

testnew: npn
	/usr/bin/time ./npn

