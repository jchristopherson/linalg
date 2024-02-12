! linalg_test.f90

! A testing application.
program main
    use iso_fortran_env, only : int32
    use test_qr
    use test_svd_ops
    use test_mtx_inverse
    use test_cholesky
    use test_eigen
    use test_misc
    use test_lu
    use test_sort
    use test_lq
    use test_sparse

    ! Local Variables
    logical :: rst
    integer(int32) :: flag

    ! Initialization
    flag = 0

    ! Misc. Item Tests
    rst = test_diagonal_mtx_mult()
    if (.not.rst) flag = 1

    rst = test_diagonal_mtx_mult_cmplx()
    if (.not.rst) flag = 2

    rst = test_rank1_update()
    if (.not.rst) flag = 3

    rst = test_rank1_update_cmplx()
    if (.not.rst) flag = 4
    
    rst = test_rank()
    if (.not.rst) flag = 5
    
    rst = test_tri_mtx_mult_1()
    if (.not.rst) flag = 6

    rst = test_tri_mtx_mult_1_cmplx()
    if (.not.rst) flag = 7
    
    rst = test_tri_mtx_mult_2()
    if (.not.rst) flag = 8

    rst = test_tri_mtx_mult_2_cmplx()
    if (.not.rst) flag = 9

    rst = test_mtx_mult_1()
    if (.not.rst) flag = 10

    rst = test_mtx_mult_1_cmplx()
    if (.not.rst) flag = 11

    rst = test_tri_mtx_solve_1()
    if (.not.rst) flag = 12

    rst = test_tri_mtx_solve_1_cmplx()
    if (.not.rst) flag = 13
    
    ! LU Factorization Tests
    rst = test_lu_factor()
    if (.not.rst) flag = 14
    
    rst = test_lu_solve()
    if (.not.rst) flag = 15

    rst = test_lu_factor_cmplx()
    if (.not.rst) flag = 16

    rst = test_lu_solve_cmplx()
    if (.not.rst) flag = 17
    

    ! QR Factorization Tests
    rst = test_qr_factor()
    if (.not.rst) flag = 18

    rst = test_qr_factor_cmplx()
    if (.not.rst) flag = 19
    
    rst = test_qr_factor_od()
    if (.not.rst) flag = 20

    rst = test_qr_factor_od_cmplx()
    if (.not.rst) flag = 21
    
    rst = test_qr_factor_ud()
    if (.not.rst) flag = 22

    rst = test_qr_factor_ud_cmplx()
    if (.not.rst) flag = 23
    
    rst = test_qr_mult()
    if (.not.rst) flag = 24

    rst = test_qr_mult_cmplx()
    if (.not.rst) flag = 25
    
    rst = test_qr_mult_od()
    if (.not.rst) flag = 26

    rst = test_qr_mult_od_cmplx()
    if (.not.rst) flag = 27
    
    rst = test_qr_mult_ud()
    if (.not.rst) flag = 28

    rst = test_qr_mult_ud_cmplx()
    if (.not.rst) flag = 29
    
    rst = test_qr_mult_right()
    if (.not.rst) flag = 30

    rst = test_qr_mult_right_cmplx()
    if (.not.rst) flag = 31
    
    rst = test_qr_mult_right_od()
    if (.not.rst) flag = 32

    rst = test_qr_mult_right_od_cmplx()
    if (.not.rst) flag = 33
    
    rst = test_qr_mult_right_ud()
    if (.not.rst) flag = 34

    rst = test_qr_mult_right_ud_cmplx()
    if (.not.rst) flag = 35
    
    rst = test_qr_mult_vector()
    if (.not.rst) flag = 36

    rst = test_qr_mult_vector_cmplx()
    if (.not.rst) flag = 37
    
    rst = test_qr_solve_no_pivot()
    if (.not.rst) flag = 38

    rst = test_qr_solve_no_pivot_cmplx()
    if (.not.rst) flag = 39
    
    rst = test_qr_solve_pivot()
    if (.not.rst) flag = 40

    rst = test_qr_solve_no_pivot_cmplx()
    if (.not.rst) flag = 41

    rst = test_qr_solve_pivot_od()
    if (.not.rst) flag = 42

    rst = test_qr_solve_pivot_od_cmplx()
    if (.not.rst) flag = 43
    
    rst = test_qr_solve_pivot_ud()
    if (.not.rst) flag = 44

    rst = test_qr_solve_pivot_ud_cmplx()
    if (.not.rst) flag = 45
    
    rst = test_qr_update_1()
    if (.not.rst) flag = 46

    rst = test_qr_update_1_cmplx()
    if (.not.rst) flag = 47
    
    ! SVD Tests
    rst = test_svd()
    if (.not.rst) flag = 48

    rst = test_svd_cmplx()
    if (.not.rst) flag = 49
    
    rst = test_svd_od()
    if (.not.rst) flag = 49

    rst = test_svd_od_cmplx()
    if (.not.rst) flag = 50
    
    rst = test_svd_ud()
    if (.not.rst) flag = 51

    rst = test_svd_ud_cmplx()
    if (.not.rst) flag = 52
    
    ! Matrix Inverse Tests
    rst = test_pinv()
    if (.not.rst) flag = 53
    
    rst = test_pinv_od()
    if (.not.rst) flag = 54
    
    rst = test_inv()
    if (.not.rst) flag = 55

    rst = test_pinv_cmplx()
    if (.not.rst) flag = 56

    rst = test_pinv_od_cmplx()
    if (.not.rst) flag = 57

    rst = test_pinv_ud_cmplx()
    if (.not.rst) flag = 58

    rst = test_inv_cmplx()
    if (.not.rst) flag = 59
    
    ! Cholesky Factorization Tests
    rst =  test_cholesky_factor()
    if (.not.rst) flag = 60

    rst = test_cholesky_factor_cmplx()
    if (.not.rst) flag = 61
    
    rst = test_cholesky_rank1_update()
    if (.not.rst) flag = 62

    rst = test_cholesky_rank1_update_cmplx()
    if (.not.rst) flag = 63
    
    rst = test_cholesky_rank1_downdate()
    if (.not.rst) flag = 64

    rst = test_cholesky_rank1_downdate_cmplx()
    if (.not.rst) flag = 65
    
    ! Eigenvalue/Eigenvector Tests
    rst = test_eigen_symm()
    if (.not.rst) flag = 66
    
    rst = test_eigen_asymm()
    if (.not.rst) flag = 67

    rst = test_eigen_asymm_cmplx()
    if (.not.rst) flag = 68
    
    rst = test_eigen_gen()
    if (.not.rst) flag = 69

    ! Sorting Tests
    rst = test_dbl_ascend_sort()
    if (.not.rst) flag = 70

    rst = test_ascend_sort_cmplx()
    if (.not.rst) flag = 71

    rst = test_dbl_descend_sort()
    if (.not.rst) flag = 72

    rst = test_descend_sort_cmplx()
    if (.not.rst) flag = 73

    ! LQ Factorization Tests
    rst = test_lq_factor()
    if (.not.rst) flag = 74

    rst = test_lq_factor_ud()
    if (.not.rst) flag = 75

    rst = test_lq_factor_cmplx()
    if (.not.rst) flag = 76

    rst = test_lq_factor_ud_cmplx()
    if (.not.rst) flag = 77

    rst = test_lq_mult()
    if (.not.rst) flag = 78

    rst = test_lq_mult_ud()
    if (.not.rst) flag = 79

    rst = test_lq_mult_cmplx()
    if (.not.rst) flag = 80

    rst = test_lq_mult_cmplx_ud()
    if (.not.rst) flag = 81

    rst = test_lq_mult_right()
    if (.not.rst) flag = 82

    rst = test_lq_mult_right_cmplx()
    if (.not.rst) flag = 83

    rst = test_lq_mult_right_ud()
    if (.not.rst) flag = 84

    rst = test_lq_mult_right_cmplx_ud()
    if (.not.rst) flag = 85

    ! Banded Matrices
    rst = test_banded_mtx_mult_dbl()
    if (.not.rst) flag = 86

    rst = test_banded_mtx_mult_cmplx()
    if (.not.rst) flag = 87

    rst = test_band_diag_mtx_mult_dbl()
    if (.not.rst) flag = 88

    rst = test_band_diag_mtx_mult_cmplx()
    if (.not.rst) flag = 89

    ! Sparse Matrices
    rst = test_csr_1()
    if (.not.rst) flag = 90

    rst = test_csr_mult_1()
    if (.not.rst) flag = 91

    rst = test_csr_mult_2()
    if (.not.rst) flag = 92

    rst = test_csr_add_1()
    if (.not.rst) flag = 93

    rst = test_csr_subtract_1()
    if (.not.rst) flag = 94

    rst = test_csr_scalar_mult_1()
    if (.not.rst) flag = 95

    rst = test_csr_scalar_divide_1()
    if (.not.rst) flag = 96

    rst = test_csr_transpose_1()
    if (.not.rst) flag = 97

    rst = test_csr_diag_mult_1()
    if (.not.rst) flag = 98

    rst = test_csr_diag_mult_2()
    if (.not.rst) flag = 99

    rst = test_csr_sparse_direct_solve_1()
    if (.not.rst) flag = 100

    rst = test_diag_to_csr_1()
    if (.not.rst) flag = 101

    rst = test_banded_to_csr_1()
    if (.not.rst) flag = 102

    rst = test_extract_diagonal_csr_1()
    if (.not.rst) flag = 103

    rst = test_msr_1()
    if (.not.rst) flag = 104

    rst = test_csr_lu_factor_1()
    if (.not.rst) flag = 105

    ! End
    if (flag /= 0) stop flag
end program
