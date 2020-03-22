// c_linalg_eigen_example.c

/*
 * This is an example illustrating the use of the eigenvalue and eigenvector 
 * routines to solve a free vibration problem of 3 masses connected by springs.
 *
 *     k1           k2           k3           k4
 * |-\/\/\-| m1 |-\/\/\-| m2 |-\/\/\-| m3 |-\/\/\-|
 *
 * As illustrated above, the system consists of 3 masses connected by springs.
 * Spring k1 and spring k4 connect the end masses to ground.  The equations of 
 * motion for this system are as follows.
 *
 * | m1  0   0 | |x1"|   | k1+k2  -k2      0  | |x1|   |0|
 * | 0   m2  0 | |x2"| + |  -k2  k2+k3    -k3 | |x2| = |0|
 * | 0   0   m3| |x3"|   |   0    -k3    k3+k4| |x3|   |0|
 *
 * Notice: x1" = the second time derivative of x1.
 */

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "linalg.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define INDEX(i, j, m) ((j) * (m) + (i))
void normalize_array(int n, double *x);

int main() {
    // Constants
    const int ndof = 3;
    const double pi = 3.14159265359;
    const double m1 = 0.5;
    const double m2 = 2.5;
    const double m3 = 0.75;
    const double k1 = 5.0e6;
    const double k2 = 10.0e6;
    const double k3 = 10.0e6;
    const double k4 = 5.0e6;

    // Local Variables
    int i, j, flag;
    double m[9], k[9], beta[3], natFreq[3], modeShapes[9];
    double complex alpha[3], vals[3], vecs[9];

    // Build the system matrices
    m[0] = m1;          m[3] = 0.0;         m[6] = 0.0;
    m[1] = 0.0;         m[4] = m2;          m[7] = 0.0;
    m[2] = 0.0;         m[5] = 0.0;         m[8] = m3;

    k[0] = k1 + k2;     k[3] = -k2;         k[6] = 0.0;
    k[1] = -k2;         k[4] = k2 + k3;     k[7] = -k3;
    k[2] = 0.0;         k[5] = -k3;         k[8] = k3 + k4;

    // Compute the eigenvalues and eigenvectors
    flag = la_eigen_gen(true, ndof, k, ndof, m, ndof, alpha, beta, vecs, ndof);

    // Compute the eigenvalues from their components
    for (i = 0; i < ndof; ++i) vals[i] = alpha[i] / beta[i];

    // Sort the eigenvalues and eigenvectors
    flag = la_sort_eigen_cmplx(true, ndof, vals, vecs, ndof);

    // Compute the natural frequencies and extract the mode shape info
    for (i = 0; i < ndof; ++i) {
        natFreq[i] = sqrt(creal(vals[i])) / (2.0 * pi);
        
        // Extract the real components - the imaginary component is zero
        for (j = 0; j < ndof; ++j) {
            modeShapes[INDEX(j,i,ndof)] = creal(vecs[INDEX(j,i,ndof)]);
        }

        // Normalize the mode shape
        normalize_array(ndof, &modeShapes[INDEX(0,i,ndof)]);
    }

    // Print out each mode shape
    printf("Modal Information:\n");
    for (i = 0; i < ndof; ++i) {
        printf("Mode %i: (%f Hz)\n", i + 1, natFreq[i]);
        for (j = 0; j < ndof; ++j) {
            printf("\t%f\n", modeShapes[INDEX(j, i, ndof)]);
        }
    }

    // End
    return 0;
}

void normalize_array(int n, double *x) {
    // Local Variables
    int i;
    double val, maxval;

    // Find the largest magnitude value
    maxval = x[0];
    for (i = 1; i < n; ++i) {
        val = x[i];
        if (fabs(val) > fabs(maxval)) maxval = val;
    }

    // Normalize the array
    for (i = 0; i < n; ++i) x[i] /= maxval;
}
