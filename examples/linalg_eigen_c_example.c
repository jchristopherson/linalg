// linalg_eigen_c_example.c

#include <stdio.h>
#include <math.h>
#include "linalg.h"

// A simple macro for addressing matrices stored in column-major format.
#define INDEX(i, j, m) ((i) + (j) * (m))

// This is an example illustrating the use of the eigenvalue and eigenvector 
// routines to solve a free vibration problem of 3 masses connected by springs.
//
//     k1           k2           k3           k4
// |-\/\/\-| m1 |-\/\/\-| m2 |-\/\/\-| m3 |-\/\/\-|
//
// As illustrated above, the system consists of 3 masses connected by springs.
// Spring k1 and spring k4 connect the end masses to ground.  The equations of 
// motion for this system are as follows.
//
// | m1  0   0 | |x1"|   | k1+k2  -k2      0  | |x1|   |0|
// | 0   m2  0 | |x2"| + |  -k2  k2+k3    -k3 | |x2| = |0|
// | 0   0   m3| |x3"|   |   0    -k3    k3+k4| |x3|   |0|
//
// Notice: x1" = the second time derivative of x1.
int main() {
    // Model Parameters
    const double pi = 3.14159265359,
        m1 = 0.5,
        m2 = 2.5,
        m3 = 0.75,
        k1 = 5.0e6,
        k2 = 10.0e6,
        k3 = 10.0e6,
        k4 = 5.0e6;
    
    // Matrices (column-major format)
    double m[9] = {m1, 0.0, 0.0, 0.0, m2, 0.0, 0.0, 0.0, m3},
        k[9] = {k1 + k2, -k2, 0.0, -k2, k2 + k3, -k3, 0.0, -k3, k3 + k4};
    
    // Additional Variables
    int i, j;
    double complex vals[3], modeShapes[9];
    double beta[3], natFreq[3];

    // Compute the eigenvalues and eigenvectors
    eigen_gen_(3, k, m, vals, beta, modeShapes, NULL);

    // Sort the eigenvalues and eigenvectors
    for (i = 0; i < 3; ++i) vals[i] /= beta[i]; // Compute the full eigenvalue
    sort_eigen_cmplx_(true, 3, vals, modeShapes);

    // Compute the natural frequency values (units = Hz)
    for (i = 0; i < 3; ++i) natFreq[i] = sqrt(creal(vals[i])) / (2.0 * pi);

    // Print the results
    printf("Modal Information:\n");
    for (i = 0; i < 3; ++i) {
        printf("Mode %i: (%f Hz)\n", i + 1, natFreq[i]);
        for (j = 0; j < 3; ++j) {
            printf("\t%f\n", creal(modeShapes[INDEX(j, i, 3)]));
        }
    }

    // End
    return 0;
}