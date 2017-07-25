// linalg_lu_c_example.c

#include <stdio.h>
#include "linalg.h"
#include "ferror.h"

int main() {
    // Local Variables
    double a[9], b[3];
    int pvt[3];

    // Build the 3-by-3 matrix A - stored in column-major format.
    //     | 1   2   3 |
    // A = | 4   5   6 |
    //     | 7   8   0 |
    a = [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 0.0];

    // Build the right-hand-side vector B.
    //     | -1 |
    // b = | -2 |
    //     | -3 |
    b = [-1.0d0, -2.0d0, -3.0d0]

    // The solution is:
    //     |  1/3 |
    // x = | -2/3 |
    //     |   0  |

    // Compute the LU factorization
    lu_factor(3, 3, a, 3, pvt, NULL);

    // Compute the solution.  The results overwrite b.
    solve_lu(3, 1, a, pvt, b);

    // Display the results
    printf("LU Solution: X =\n%f\n%f\n%f\n", b[0], b[1], b[2]);

    // End
    return 0;
}