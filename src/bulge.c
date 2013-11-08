
/* Francis algorithm in C; see README */

// stdio for debugging
#include <stdio.h>

#include <cblas.h>
#include "bulge.h"
#include "util.h"


/* for this file only */
#define MYSIGN(x) ((x > 0) - (x < 0))



/** Create the Householder symmetric matrix v*v^T
 *
 * This matrix is used to process the rows and columns of the input matrix.
 */
static void create_house_matrix_packed(size_t order, double shift, double *source, size_t incs, double *hhp) {
	double h[order];
	int i;

	/* zero out the destination hhp (it's packed aka triangular, so the size is
	 * non-square) */
	for (i = 0; i < order * (order + 1) / 2; i++) {
		hhp[i] = 0;
	}

	/* create and normalize householder vector h */
	cblas_dcopy(order, source, incs, h, 1);
	h[0] += MYSIGN(h[0]) * cblas_dnrm2(order, h, 1);
	h[0] -= shift;
	cblas_dscal(order, 1.0 / cblas_dnrm2(order, h, 1), h, 1);

	/* hhp = h h^T */
	cblas_dspr(CblasRowMajor, CblasUpper, order, 1.0, h, 1, hhp);
}



/** Create bulge
 *
 * Creates the bulge structure for matrix M. Create the bulge with as many shifts as
 * possible.
 *
 * Returns the number of chasing steps needed to chase the bulge or -1 on error.
 */
int form_bulge(struct bulge_info *bi, size_t order, double *M, size_t nshifts,
               double *shifts, enum chase_direction direction) {
	size_t bulge_size = nshifts + 2;
	size_t bulge_position;
	int read_direction;
	double *tempM;

	double vv[bulge_size * (bulge_size + 1)/2];

	size_t shiftidx, r, c;


	/* populate bulge_info structure now so it can serve as a useful "return
	 * value" */
	bi->order = order;
	bi->M = M;
	bi->nshifts = nshifts;
	bi->shifts = shifts;
	bi->nshifts_applied = 0;
	bi->steps_chased = 0;
	bi->direction = direction;

	/* apply as many shifts as we can to "build the bulge" */
	for (shiftidx = 0; (shiftidx < nshifts) && (shiftidx <= order - 1); shiftidx++) {
		bulge_size = shiftidx + 2;

		switch (direction) {
		case CHASE_FORWARD:
			bulge_position = 0;
			read_direction = 1;

			/* build up vv by pulling out v from top to bottom */
			r = 0;
			c = 0;
			tempM = &M[c + r*order];
			create_house_matrix_packed(bulge_size, shifts[shiftidx], tempM, order, vv);
			break;

		case CHASE_BACKWARD:
			bulge_position = order - bulge_size;
			read_direction = -1;

			/* build up vv by pulling out v from right to left (the vector gets
			 * read backwards--i.e. ending at the indicated point--when the
			 * stride is negative) */
			r = order - 1;
			c = order - 2 - shiftidx;
			tempM = &M[c + r*order];
			create_house_matrix_packed(bulge_size, shifts[shiftidx], tempM, -1, vv);
			break;

		default:
			return FORM_BULGE_ERROR;
		}

		/* use vv to process each small col and row which intersects with the bulge zone */
		for (c = 0; c < order; c++) {
			r = bulge_position;
			tempM = &M[c + r*order];
			cblas_dspmv(CblasRowMajor, CblasUpper, bulge_size, -2.0, vv,
				tempM, read_direction*order,
				1.0,
				tempM, read_direction*order);
		}
		for (r = 0; r < order; r++) {
			c = bulge_position;
			tempM = &M[c + r*order];
			cblas_dspmv(CblasRowMajor, CblasUpper, bulge_size, -2.0, vv,
				tempM, read_direction*1,
				1.0,
				tempM, read_direction*1);
		}

		/* OPTIMIZATION: we can unroll the first few hits to the above loops
		 * into a level 3 BLAS operation */

		/* in case we don't complete all shifts for some reason... */
		bi->nshifts_applied = shiftidx + 1;
	}

	return (bi->order - 2); /* number of shifts needed to eradicate the bulge */
}



/** Incrementally chase bulge
 *
 * returns number of steps left to do
 */
int chase_bulge(struct bulge_info *bi) {
	size_t bulge_size = bi->nshifts_applied + 2;
	size_t bulge_position;
	double vv[bulge_size * (bulge_size + 1)/2];
	int r, c;
	
	int steps_remaining = bi->order - 2 - bi->steps_chased;
	
	if (steps_remaining <= 0) {
		return 0;
	}

	switch (bi->direction) {
	case CHASE_FORWARD:
		bulge_position = bi->steps_chased;

		r = bulge_position + 1;
		c = bulge_position;
		printf("bp=%u bs=%u r=%u c=%u\n", bulge_position, bulge_size, r, c);
		create_house_matrix_packed(bulge_size - 1, 0.0, &bi->M[c + r*bi->order], bi->order, vv);
		break;
	case CHASE_BACKWARD:
		bulge_position = bi->order - bulge_size - bi->steps_chased;

		r = bulge_position + bulge_size - 1;
		c = bulge_position + bulge_size - 2;
		printf("bp=%u bs=%u r=%u c=%u\n", bulge_position, bulge_size, r, c);
		create_house_matrix_packed(bulge_size - 1, 0.0, &bi->M[c + r*bi->order], -1, vv);
		break;
	}

	/* use vv to process each small col and row which intersects with the bulge zone */
	for (c = 0; c < bi->order; c++) {
		r = bulge_position + 1;
		/* this is still running "backwards" */
		cblas_dspmv(CblasRowMajor, CblasUpper, bulge_size - 1, -2.0, vv, &bi->M[c + r*bi->order], bi->order, 1.0, &bi->M[c + r*bi->order], bi->order);
	}
	for (r = 0; r < bi->order; r++) {
		c = bulge_position + 1;
		/* this is still running the wrong way when we go backwards... */
		cblas_dspmv(CblasRowMajor, CblasUpper, bulge_size - 1, -2.0, vv, &bi->M[c + r*bi->order], 1,         1.0, &bi->M[c + r*bi->order], 1);
	}


	/* keep an accurate count */
	bi->steps_chased++;
	steps_remaining--;

	return steps_remaining;
}

