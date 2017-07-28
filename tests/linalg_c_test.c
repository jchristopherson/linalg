// linalg_c_test.c

#include <stdio.h>
#include "c_test_core.h"

int main() {
    // Local Variables
    bool rst, overall = true;

    // Misc Tests
    rst = test_diagonal_mtx_mult();
    if (!rst) overall = false;

    // End
    if (overall) printf("C API LINALG TEST STATUS: PASS\n");
    else printf("C API LINALG TEST STATUS: FAIL\n");
    return 0;
}