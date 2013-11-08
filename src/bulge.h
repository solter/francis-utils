
/* Francis algorithm in C; see README */

#ifndef BULGE_H
#define BULGE_H
#include <stddef.h>

struct bulge_info {
	size_t nshifts;
	size_t nshifts_applied;
	double *shifts;
	size_t order;
	double *M;
	size_t steps_chased;
};

int form_bulge(struct bulge_info *bi,
               size_t order,
               double *M,
               size_t nshifts,
               double *shifts);

int chase_bulge(struct bulge_info *bi);

#endif

