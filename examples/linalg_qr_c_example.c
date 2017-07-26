// linalg_qr_c_example.c

#include <stdio.h>
#include "linalg.h"
#include "ferror.h"

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
    double tau[3],
        a[9] = {1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 0.0},
        b[3] = {-1.0, -2.0, -3.0};
    int pvt[3];

    // The solution is:
    //     |  1/3 |
    // x = | -2/3 |
    //     |   0  |

    // Compute the QR factorization
    qr_factor_pivot_(3, 3, a, 3, tau, pvt, NULL);

    // Compute the solution.  The results overwrite b.
    solve_qr_pivot_(3, 3, 1, a, 3, tau, pvt, 3, b);

    // Display the results
    printf("QR Solution: X =\n%f\n%f\n%f\n", b[0], b[1], b[2]);

    // End
    return 0;
}
