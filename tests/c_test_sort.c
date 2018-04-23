// c_test_lu.c

#include "c_test_core.h"

bool test_ascending_sort() {
    // Local Variables
    const int n = 1000;
    double x[n];
    bool rst;
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(n, 1, x);

    // Sort the array into ascending order
    sort_dbl(true, n, x, NULL);

    // Ensure X is now ascending
    for (i = 1; i < n; ++i) {
        if (x[i] < x[i-1]) {
            rst = false;
            printf("Test Failed: Sorted array is not in ascending order.\n");
            break;
        }
    }

    // End
    return rst;
}

bool test_descending_sort() {
    // Local Variables
    const int n = 1000;
    double x[n];
    bool rst;
    int i;

    // Initialization
    rst = true;
    make_rand_mtx(n, 1, x);

    // Sort the array into descending order
    sort_dbl(false, n, x, NULL);

    // Ensure X is now ascending
    for (i = 1; i < n; ++i) {
        if (x[i] > x[i-1]) {
            rst = false;
            printf("Test Failed: Sorted array is not in descending order.\n");
            break;
        }
    }

    // End
    return rst;
}
