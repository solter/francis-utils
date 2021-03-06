
/* Francis algorithm in C; see README */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "bulge.h"
#include "util.h"

/* NOTE: BLAS routine conventions: target is last, size comes before operand,
 * incX comes after */



const size_t Anorm_N = 7, Ahess_N = 7;

/** A is a random test matrix with eigenvalues 2,3,5,7,11,13,17 */
double Anorm[] = {
1.2205,-0.9778,-3.5875,6.3239,-18.6134,-11.5920,-27.5207,
15.4760,14.0774,9.7466,-9.1991,32.6386,10.2165,44.3215,
-5.3925,-6.9141,-0.3162,-5.2980,-8.8870,-4.6269,-7.1133,
-3.2413,1.2535,-0.0849,11.0866,-6.0826,-1.4200,-12.6913,
-6.0873,-0.6122,-5.6896,9.1950,-7.6241,-4.1090,-21.6177,
-7.5542,-8.1024,-4.4474,-2.2788,-5.1283,6.1540,-2.6075,
12.7377,9.4829,10.2203,-1.8157,21.6712,10.9930,33.4019,
};



/** Ahess is the hessenberg form of A */
double Ahess[] = {
1.2205,7.1720,-4.6229,-10.5124,21.2274,-11.3153,-23.0235,
-23.1400,27.3588,-5.6390,-17.8708,38.2222,-38.6213,-42.2665,
0,5.4755,-2.2132,-7.4969,5.9252,0.8973,-20.5604,
0,0,1.3863,5.5655,-2.7729,1.4370,-2.3967,
0,0,0,-1.7075,13.9845,-4.1728,-4.6832,
0,0,0,0,2.1182,4.8584,-4.7575,
0,0,0,0,0,-3.3022,7.2255,
};

double Simp[] = {
	2,2,2,2,
	1,2,2,2,
	0,1,2,2,
	0,0,1,2,
};



void test_two_way(void) {
	double *M = Ahess;
	size_t N = Ahess_N;

	struct bulge_info forward, backward;

	double forward_shifts[] = {2};
	double backward_shifts[] = {3};

	ssmd("original matrix", N, M);

	form_bulge(&forward, N, M, 1, forward_shifts, CHASE_FORWARD);
	ssmd("with forward shift", N, M);

	form_bulge(&backward, N, M, 1, backward_shifts, CHASE_BACKWARD);
	ssmd("with backward shift as well", N, M);

	chase_bulge_step(&forward);
	ssmd("chased forward", N, M);

	chase_bulge_step(&backward);
	ssmd("chased backward", N, M);

	chase_bulge_step(&forward);
	ssmd("chased forward", N, M);

	chase_bulge_step(&backward);
	ssmd("chased backward (nonsensical)", N, M);
}

void test_qr_algorithm(void) {
	double *M = Ahess;
	size_t N = 7;
	double error, tol = 0.000001;

	size_t i;
	struct bulge_info bi;

	double shifts[2];

	do {

		/* use last two sub-diagonal entries as shifts */
		shifts[0] = M[N-2 + (N-1)*N];
		shifts[1] = M[N-3 + (N-2)*N];

		/* use last two diagonal entries as shifts */
		//shifts[0] = M[N-1 + (N-1)*N];
		//shifts[1] = M[N-2 + (N-2)*N];

		build_and_chase_bulge(N, M, 2, shifts, CHASE_FORWARD);

		error = cblas_dnrm2(N-1, &M[N], N+1);

		printf("iteration complete with shifts %f %f, error = %6.18f\n", shifts[0], shifts[1], error);
		ssm(N, M);
		printf("press <enter> to do another iteration, or q<enter> to quit\n");
		if (getchar() == 'q') return;

	} while (error > tol);

	ssmd("output of QR algorithm", N, M);
}

void test_development(void) {
	double shifts[] = {0, 0, 0, 1,1,2,24,24,2,42,4};
	double *M = Ahess;
	size_t N = 7;
	int i;

	struct bulge_info b;

	ssmd("%original matrix", N, M);
#if 0
	puts("\n\n\n\n\n\nRUNNING FORWARD...");
	i = form_bulge(&b, N, M, 2, shifts, CHASE_FORWARD);
#else
	puts("\n\n\n\n\n\nRUNNING BACKWARD...");
	i = form_bulge(&b, N, M, 2, shifts, CHASE_BACKWARD);
#endif
	ssmd("%new shiny bulge", b.order, b.M);
	printf("%%need to chase it for %u steps; start the chase!\n\n", i);
	do {
		i = chase_bulge_step(&b);
		printf("%%have %u steps to go, just completed number %lu, here's the result:\n\n", i, b.steps_chased);
		ssm(b.order, b.M);
		printf("eig(M)\n");
	} while (i > 0);
}



int main(int argc, char *argv[]) {
	//test_two_way();
	test_qr_algorithm();
	return EXIT_SUCCESS;
}

