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

    // End
    return flag;
}
