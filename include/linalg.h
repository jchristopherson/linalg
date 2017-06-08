// linalg.h
#ifndef LINALG_H_DEFINED
#define LINALG_H_DEFINED

#include <stdbool.h>
#include <complex.h>
#include "ferror/include/ferror.h"

#ifdef __cplusplus
extern "C" {
#endif

void mtx_mult(int m, int n, int k, double alpha, const double *a, 
              const double *b, double beta, double *c);
void diag_mtx_mult(bool trans, int m, int n, double alpha, int na, 
                   const double *a, int mb, int nb, const double *b, 
                   double beta, double *c);
void diag_mtx_mult_cmplx(bool trans, int m, int n, double alpha, int na,
                         const double complex *a, int mb, int nb, 
                         const double *b, double beta, double complex *c);

#ifdef __cplusplus
}
#endif
#endif // LINALG_H_DEFINED