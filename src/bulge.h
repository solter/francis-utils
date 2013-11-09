
/* Francis algorithm in C; see README */

#ifndef BULGE_H
#define BULGE_H
#include <stddef.h>



#define FORM_BULGE_ERROR (-1)



/** Parameter to form_bulge(...) which determines whether the bulge steps
 * top-left to bottom-right or the opposite direction each time form_bulge(...)
 * is called. */
enum chase_direction {
	CHASE_FORWARD = (int) 'F',
	CHASE_BACKWARD = (int) 'B'
};



/** Housekeeping structure for form_bulge and chase_bulge_step. */
struct bulge_info {
	size_t nshifts;
	size_t nshifts_applied;
	double *shifts;
	size_t order;
	double *M;
	enum chase_direction direction;
	size_t steps_chased;
};



int form_bulge(struct bulge_info *bi,
               size_t order,
               double *M,
               size_t nshifts,
               double *shifts,
               enum chase_direction direction);

int chase_bulge_step(struct bulge_info *bi);

int build_and_chase_bulge(size_t order, double *M, size_t nshifts, double *shifts,
                           enum chase_direction direction);



#endif

