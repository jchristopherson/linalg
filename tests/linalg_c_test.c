// linalg_c_test.c

#include <stdio.h>
#include "c_test_core.h"

int main() {
    // Local Variables
    bool rst, overall = true;

    // Misc Tests
    rst = test_diagonal_mtx_mult();
    if (!rst) overall = false;

    rst = test_rank1_update();
    if (!rst) overall = false;

    rst = test_rank();
    if (!rst) overall = false;

    rst = test_tri_mtx_mult();
    if (!rst) overall = false;

    rst = test_lu_factor();
    if (!rst) overall = false;

    rst = test_lu_solve();
    if (!rst) overall = false;

    rst = test_qr_factor();
    if (!rst) overall = false;

    rst = test_qr_factor_od();
    if (!rst) overall = false;

    rst = test_qr_factor_ud();
    if (!rst) overall = false;

    // End
    if (overall) printf("C API LINALG TEST STATUS: PASS\n");
    else printf("C API LINALG TEST STATUS: FAIL\n");
    return 0;
}