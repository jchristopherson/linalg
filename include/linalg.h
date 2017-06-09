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
void rank1_update(int m, int n, double alpha, const double *x, const double *y,
                  double *a);
double trace(int m, int n, const double *x);
int mtx_rank(int m, int n, double *a, errorhandler err);
double det(int n, double *a, errorhandler err);
void swap(int n, double *x, double *y);
void lu_factor(int m, int n, double *a, int ni, int *ipvt, errorhandler err);
void form_lu(int n, double *lu, const int *ipvt, double *u, double *p);
void qr_factor(int m, int n, double *a, int nt, double *tau, errorhandler err);
void qr_factor_pivot(int m, int n, double *a, int nt, double *tau, int *jpvt,
                     errorhandler err);
void form_qr(int m, int n, double *r, int nt, const double *tau, double *q,
             errorhandler err);
void form_qr_pivot(int m, int n, double *r, int nt, const double *tau,
                   const int *pvt, double *q, double *p, errorhandler err);
void mult_qr(bool trans, int m, int n, double *q, int nt, const double *tau,
             double *c, errorhandler err);
void qr_rank1_update(int m, int n, double *q, double *r, double *u, double *v, 
                     errorhandler err);
void cholesky_factor(int n, double *a, bool upper, errorhandler err);
void cholesky_rank1_update(int n, double *r, double *u, errorhandler err);
void rz_factor(int m, int n, double *a, double *tau, errorhandler err);
void mult_rz(bool trans, int m, int n, double *a, const double *tau, double *c,
             errorhandler err);
void svd(int m, int n, double *a, int ns, double *s, double *u, double *vt,
         errorhandler err);

#ifdef __cplusplus
}
#endif
#endif // LINALG_H_DEFINED