
/* Francis algorithm in C; see README */

#include "util.h"
#include <stdio.h>



/** Display a matrix as follows:
 *
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     ...            ...            ...
 *
 */
void showmat(size_t nrows, size_t ncols, double *M) {
	int r, c;

	printf("M=[");
	for (r = 0; r < nrows; r++) {
		for (c = 0; c < ncols; c++) {
			printf("\t%+3.6f", M[c + r*ncols]);

			/* formatting and frill */
			if (c == ncols - 1) {
				if (r == nrows - 1) {
					printf("];\n");
				} else {
					printf(";\n");
				}
			} else {
				putchar(',');
			}
		}
	}
}



/** Display a square matrix as follows:
 *
 * <description>:
 *
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     +nnn.nnnnnn    +nnn.nnnnnn    +nnn.nnnnnn    ...
 *     ...            ...            ...
 *
 */
void ssm(size_t order, double *M) {
	showmat(order, order, M);
	printf("\n");
}
void ssmd(const char *description, size_t order, double *M) {
	printf("%s:\n\n", description);
	showmat(order, order, M);
	printf("\n");
}

void scvd(const char *description, size_t length, double *v) {
	printf("%s:\n\n", description);
	showmat(length, 1, v);
}

