
/* Francis algorithm in C; see README */

#ifndef UTIL_H
#define UTIL_H
#include <stddef.h>

void showmat(size_t nrows, size_t ncols, double *M);
void ssm(size_t order, double *M);
void ssmd(const char *description, size_t order, double *M);
void scvd(const char *description, size_t length, double *v);

#endif

