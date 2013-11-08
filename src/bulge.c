
/* Francis algorithm in C; see README */

// stdio for debugging
#include <stdio.h>

#include <cblas.h>
#include "bulge.h"
#include "util.h"


/* for this file only */
#define MYSIGN(x) ((x > 0) - (x < 0))



/** Create bulge
 *
 * Creates the bulge structure for matrix M. Create the bulge with as many shifts as
 * possible.
 *
 * Returns the number of steps needed to chase the bulge.
 */
int form_bulge(struct bulge_info *bi, size_t order, double *M, size_t nshifts, double *shifts) {
	size_t bulge_size = nshifts + 2;

	double v[bulge_size]; /* work vector for householder vectors and chasing the bulge */
	const size_t vv_sz = bulge_size * (bulge_size + 1)/2;
	double vv[vv_sz]; /* for holding v * v^T */

	size_t i, j, r, c;

	/* populate bulge_info structure now so it can serve as a useful "return
	 * value" */
	bi->order = order;
	bi->M = M;
	bi->nshifts = nshifts;
	bi->shifts = shifts;
	bi->nshifts_applied = 0;
	bi->steps_chased = 0;

	/* apply as many shifts as we can to "build the bulge" */
	for (i = 0; (i < nshifts) && (i <= order - 1); i++) {
		bulge_size = i + 2;

		/* build normalized Householder vector */
		cblas_dcopy(bulge_size, M, order, v, 1);
		v[0] += MYSIGN(v[0]) * cblas_dnrm2(bulge_size, M, order);
		v[0] -= shifts[i];
		cblas_dscal(bulge_size, 1.0 / cblas_dnrm2(bulge_size, v, 1), v, 1);

		/* vv = v v^T */
		for (j = 0; j < vv_sz; j++) vv[j] = 0; /* zero out vv */
		cblas_dspr(CblasRowMajor, CblasUpper, bulge_size, 1.0, v, 1, vv);

		/* use vv to process each small col and row which intersects with the bulge zone */
		for (c = 0; c < order; c++) {
			r = 0;
			cblas_dspmv(CblasRowMajor, CblasUpper, bulge_size, -2.0, vv, &M[c + r*order], order, 1.0, &M[c + r*order], order);
		}
		for (r = 0; r < order; r++) {
			c = 0;
			cblas_dspmv(CblasRowMajor, CblasUpper, bulge_size, -2.0, vv, &M[c + r*order], 1,     1.0, &M[c + r*order], 1);
		}

		/* OPTIMIZATION: we can unroll the first few hits to the above loops */

		/* in case we don't complete all shifts for some reason... */
		bi->nshifts_applied = i + 1;
	}

	return (bi->order - 2); /* number of shifts needed to eradicate the bulge */
}



/** Incrementally chase bulge
 *
 * returns number of steps left to do
 */
int chase_bulge(struct bulge_info *bi) {
	size_t bulge_size = bi->nshifts_applied + 2;

	double v[bulge_size]; /* work vector for householder vectors and chasing the bulge */
	const size_t vv_sz = bulge_size * (bulge_size + 1)/2;
	double vv[vv_sz]; /* for holding v * v^T */
	size_t b_N = bi->nshifts_applied + 1;

	int bulge_offset, j, r, c;
	
	if (bi->steps_chased < bi->order - 2) {
		bulge_offset = bi->steps_chased;
		b_N = bi->nshifts_applied + 1;

		/* build normalized vector */
		cblas_dcopy(b_N, &(bi->M[bulge_offset + (bulge_offset+1)*bi->order]), bi->order, v, 1);
		v[0] += MYSIGN(v[0]) * cblas_dnrm2(b_N, v, 1);
		cblas_dscal(b_N, 1.0 / cblas_dnrm2(b_N, v, 1), v, 1);

		/* vv = v v^T */
		for (j = 0; j < vv_sz; j++) vv[j] = 0; /* zero out vv */
		cblas_dspr(CblasRowMajor, CblasUpper, b_N, 1.0, v, 1, vv);

		/* use vv to process each small col and row which intersects with the bulge zone */
		for (c = 0; c < bi->order; c++) {
			r = bulge_offset + 1;
			cblas_dspmv(CblasRowMajor, CblasUpper, b_N, -2.0, vv, &bi->M[c + r*bi->order], bi->order, 1.0, &bi->M[c + r*bi->order], bi->order);
		}
		for (r = 0; r < bi->order; r++) {
			c = bulge_offset + 1;
			cblas_dspmv(CblasRowMajor, CblasUpper, b_N, -2.0, vv, &bi->M[c + r*bi->order], 1,     1.0, &bi->M[c + r*bi->order], 1);
		}

		bi->steps_chased++;
	}
	return (bi->order - 2 - bi->steps_chased);
}

