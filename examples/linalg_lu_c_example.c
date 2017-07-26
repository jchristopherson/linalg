// linalg_lu_c_example.c

#include <stdio.h>
#include "linalg.h"

int main() {
    // Define local variables and populate the matrices (column-major format):
    //
    //     | 1   2   3 |
    // A = | 4   5   6 |
    //     | 7   8   0 |
    //
    //     | -1 |
    // b = | -2 |
    //     | -3 |
    double a[9] = {1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 0.0},
        b[3] = {-1.0, -2.0, -3.0};
    int pvt[3];

    // The solution is:
    //     |  1/3 |
    // x = | -2/3 |
    //     |   0  |

    // Compute the LU factorization
    lu_factor_(3, 3, a, 3, pvt, NULL);

    // Compute the solution.  The results overwrite b.
    solve_lu_(3, 1, a, pvt, b);

    // Display the results
    printf("LU Solution: X =\n%f\n%f\n%f\n", b[0], b[1], b[2]);

    // End
    return 0;
}