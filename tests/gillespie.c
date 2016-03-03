#include <assert.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdio.h>
#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include "SFMT.h"	// Fast random number generation.  See http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
#include <immintrin.h>  // Intel intrinsics, see https://software.intel.com/sites/landingpage/IntrinsicsGuide/


/**************************************************************************************************************************/
// Useful #defines
// 
/**************************************************************************************************************************/
#define START_TIMER gettimeofday(&start, NULL)
#define STOP_TIMER gettimeofday(&end, NULL)
#define DELTA_TIMER ((end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/1000000.0)
#define POND_SIZE_X (8192)
#define POND_SIZE_Y (8192)

/**************************************************************************************************************************/
// Data types and global variables.
//
/**************************************************************************************************************************/
/*
static struct cell {
	uint64_t genome;		// Doesn't do anything yet.
	uint64_t ts_last_executed;	// When was the last time it was chosen at random.
	uint64_t ts_last_modified;	// When was the last time the genome was rewritten.
}[POND_SIZE_X][POND_SIZE_Y] __attribute__ ((aligned (128)));
*/
// List of random numbers.
static uint64_t rnd[ POND_SIZE_X * POND_SIZE_Y ] __attribute__ ((aligned (4096)));
static struct timeval start, end;

/**************************************************************************************************************************/
// Function declarations
//
/**************************************************************************************************************************/
static        void     init_genrand(unsigned long s);	// Original nanopond code
static inline uint32_t genrand_int32();			// 
static        void     sequential_write_test();

/**************************************************************************************************************************/
// main
//
/**************************************************************************************************************************/
int main(){
	sequential_write_test();
}

static void sequential_write_test(){
	uint64_t i, j;
	sfmt_t sfmt;	
	sfmt_init_gen_rand(&sfmt, 13);	// 13 is random number seed.
	init_genrand(13);

	fprintf(stdout, "rnd = %p\n", rnd);

	for(j=0; j<5; j++){

		START_TIMER;
		for(i=0; i<(POND_SIZE_X * POND_SIZE_Y); i++){
			rnd[i] = i;
		}
		STOP_TIMER;
		fprintf(stdout, "%lfs,  Sequential Non-random value write-only fill.\n", DELTA_TIMER );

		START_TIMER;
		for(i=0; i<(POND_SIZE_X * POND_SIZE_Y); i++){
			rnd[i] = genrand_int32();
		}
		STOP_TIMER;
		fprintf(stdout, "%lfs,  Sequential Random write-only fill (nanopond 32bit).\n", DELTA_TIMER );

		START_TIMER;
		for(i=0; i<(POND_SIZE_X * POND_SIZE_Y); i++){
			rnd[i] = ((((uint64_t)genrand_int32()) << 32) ^ ((uint64_t)genrand_int32()));
		}
		STOP_TIMER;
		fprintf(stdout, "%lfs,  Sequential Random write-only fill (nanopond 64bit).\n", DELTA_TIMER );

		START_TIMER;
		for(i=0; i<(POND_SIZE_X * POND_SIZE_Y); i++){
			rnd[i] = sfmt_genrand_uint64(&sfmt);
		}
		STOP_TIMER;
		fprintf(stdout, "%lfs,  Sequential Random write-only fill (sfmt single).\n", DELTA_TIMER );

		START_TIMER;
		sfmt_init_gen_rand(&sfmt, 13);	// 13 is random number seed.  Must reinitialize.
		sfmt_fill_array64(&sfmt, rnd, POND_SIZE_X * POND_SIZE_Y);
		STOP_TIMER;
		fprintf(stdout, "%lfs,  Sequential Random write-only fill (sfmt array).\n", DELTA_TIMER );

		
#if 0
		assert( _may_i_use_cpu_feature(_FEATURE_RDRND) );
		assert( _rdrand64_step( rnd ) );
		START_TIMER;
		for(i=0; i<(POND_SIZE_X * POND_SIZE_Y); i++){
			_rdrand64_step( &rnd[i] );	
		}
		STOP_TIMER;
		fprintf(stdout, "%lfs,  Sequential Random write-only fill (intel intrinsics).\n", DELTA_TIMER );
#endif 
/*
 * icc results
 *
 * 0.118801s,  Sequential Non-random value write-only fill.
 * 0.216810s,  Sequential Random write-only fill (nanopond 32bit).
 * 0.435843s,  Sequential Random write-only fill (nanopond 64bit).
 * 0.204480s,  Sequential Random write-only fill (sfmt single).
 * 0.107825s,  Sequential Random write-only fill (sfmt array).
 * 5.891248s,  Sequential Random write-only fill (intel intrinsics).
 *
 * gcc results
 *
 * 0.060865s,  Sequential Non-random value write-only fill.
 * 0.228625s,  Sequential Random write-only fill (nanopond 32bit).
 * 0.449827s,  Sequential Random write-only fill (nanopond 64bit).
 * 0.184499s,  Sequential Random write-only fill (sfmt single).
 * 0.067045s,  Sequential Random write-only fill (sfmt array).
 */
		fprintf(stdout, "\n");

	}

	return 0;
}


/**************************************************************************************************************************/
//
// Original nanopond random number code.
//	Licence not appended:  do not redistribute.
//	See original for comments.
//
/**************************************************************************************************************************/

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   
#define UPPER_MASK 0x80000000UL 
#define LOWER_MASK 0x7fffffffUL 

static unsigned long mt[N]; 
static int mti=N+1; 

static void init_genrand(unsigned long s) {
	mt[0]= s & 0xffffffffUL;
	for (mti=1; mti<N; mti++) {
		mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
		mt[mti] &= 0xffffffffUL;
	}
}

static inline uint32_t genrand_int32() {
	uint32_t y;
	static uint32_t mag01[2]={0x0UL, MATRIX_A};

	if (mti >= N) { // generate N words at one time 
		int kk;

		for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

