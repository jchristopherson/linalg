#include <stdbool.h>
#include "c_linalg_test.h"

int main()
{
    // Variables
    bool check;
    int flag = 0;

    // Tests
    check = test_diagonal_mtx_mult();
    if (!check) flag = 1;
    check = test_cmplx_diagonal_mtx_mult();
    if (!check) flag = -1;

    check = test_rank1_update();
    if (!check) flag = 2;
    check = test_cmplx_rank1_update();
    if (!check) flag = -2;

    check = test_trace();
    if (!check) flag = 3;
    check = test_cmplx_trace();
    if (!check) flag = -3;

    // End
    return flag;
}
