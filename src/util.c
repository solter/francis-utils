
/* Francis algorithm in C; see README */

#include "util.h"
#include <stdio.h>



/** Display a matrix as follows:
 * M=[   +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *            ...           ...                   ...
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn];
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



/** Display a matrix and a description as follows:
 *
 * <description>:
 *
 * M=[   +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *            ...           ...                   ...
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn];
 */
void smd(const char *description, size_t nrows, size_t ncols, double *M) {
	printf("\n%s:\n\n", description);
	showmat(nrows, ncols, M);
}



/** Display a square matrix as follows:
 * M=[   +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *            ...           ...                   ...
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn];
 */
void ssm(size_t order, double *M) {
	showmat(order, order, M);
}



/** Display a square matrix and a description as follows:
 *
 * <description>:
 *
 * M=[   +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn,
 *            ...           ...                   ...
 *       +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn];
 */
void ssmd(const char *description, size_t order, double *M) {
	smd(description, order, order, M);
}



/** Display a row vector as follows:
 * M=[   +nnn.nnnnnn,   +nnn.nnnnnn,   ...,   +nnn.nnnnnn];
 */
void srv(size_t length, double *v) {
	showmat(1, length, v);
}



/** Display a column vector and a description as follows:
 *
 * <description>:
 *
 * M=[   +nnn.nnnnnn,
 *       +nnn.nnnnnn,
 *            ...
 *       +nnn.nnnnnn];
 */
void scvd(const char *description, size_t length, double *v) {
	smd(description, length, 1, v);
}

