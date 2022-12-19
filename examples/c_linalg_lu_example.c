#include <stdio.h>
#include <math.h>
#include "linalg.h"

#define INDEX(i, j, m) ((j) * (m) + (i))

int main() {
    // Local Variables
    int i, flag, pvt[3];

    // Build the 3-by-3 matrix A - Use column-major formating!
    //     | 1  2  3 |
    // A = | 4  5  6 |
    //     | 7  8  0 |
    double a[] = {1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 0.0};

    // Build the right-hand-side vector B.
    //     | -1 |
    // b = | -2 |
    //     | -3 |
    double b[] = {-1.0, -2.0, -3.0};

    // The solution is:
    //     |  1/3 |
    // x = | -2/3 |
    //     |   0  |

    // Compute the LU factorization
    flag = la_lu_factor(3, 3, a, 3, pvt);
    if (flag != LA_NO_ERROR) return flag;

    // Solve.  The results overwrite b.
    flag = la_solve_lu(3, 1, a, 3, pvt, b, 3);
    if (flag != LA_NO_ERROR) return flag;

    // Display the results
    printf("LU Solution: X = \n");
    for (i = 0; i < 3; i++) {
        printf("%8.4f\n", b[i]);
    }
}