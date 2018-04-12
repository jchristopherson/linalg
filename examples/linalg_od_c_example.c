// linalg_od_c_example.c

#include <stdio.h>
#include "linalg.h"

int main() {
    // Define local variables and populate the matrices (column-major format):
    //
    //     | 2   1 |
    // A = |-3   1 |
    //     |-1   1 |
    //
    //     |-1 |
    // b = |-2 |
    //     | 1 |
    double a[6] = {2.0, -3.0, -1.0, 1.0, 1.0, 1.0},
        b[3] = {-1.0, -2.0, 1.0};

    // The solution is:
    // x = [0.13158, -0.57895]**T

    // Compute the solution.  The results overwrite the first 2 elements in b.
    solve_least_squares(3, 2, 1, a, b, NULL);

    // Display the results
    printf("Least Squares Solution: X =\n%f\n%f\n", b[0], b[1]);

    // End
    return 0;
}