! linalg.f90

module linalg
    use linalg_sparse
    use linalg_basic
    use linalg_sorting
    use linalg_eigen
    use linalg_lu
    use linalg_rz
    use linalg_qr
    use linalg_tri
    use linalg_cholesky
    use linalg_lq
    use linalg_svd
    use linalg_inverse
    use linalg_least_squares
    implicit none
    private

    ! LINALG_BASIC.F90
    public :: LA_NO_OPERATION
    public :: LA_TRANSPOSE
    public :: LA_HERMITIAN_TRANSPOSE
    public :: mtx_mult
    public :: rank1_update
    public :: diag_mtx_mult
    public :: trace
    public :: mtx_rank
    public :: det
    public :: swap
    public :: recip_mult_array
    public :: tri_mtx_mult
    public :: band_mtx_mult
    public :: band_mtx_to_full_mtx
    public :: band_diag_mtx_mult
    public :: banded_to_dense
    public :: dense_to_banded
    public :: extract_diagonal

    ! LINALG_SPARSE.F90
    public :: csr_matrix
    public :: msr_matrix
    public :: size
    public :: create_empty_csr_matrix
    public :: create_empty_msr_matrix
    public :: nonzero_count
    public :: dense_to_csr
    public :: diag_to_csr
    public :: banded_to_csr
    public :: csr_to_dense
    public :: csr_to_msr
    public :: msr_to_csr
    public :: dense_to_msr
    public :: msr_to_dense
    public :: create_csr_matrix
    public :: matmul
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)
    public :: operator(/)
    public :: assignment(=)
    public :: transpose
    public :: sparse_direct_solve
    public :: pgmres_solver

    ! LINALG_SORTING.F90
    public :: sort

    ! LINALG_EIGEN.F90
    public :: eigen

    ! LINALG_LU.F90
    public :: lu_factor
    public :: form_lu
    public :: solve_lu

    ! LINALG_RZ.F90
    public :: rz_factor
    public :: mult_rz

    ! LINALG_QR.F90
    public :: qr_factor
    public :: form_qr
    public :: mult_qr
    public :: qr_rank1_update
    public :: solve_qr

    ! LINALG_TRI.F90
    public :: solve_triangular_system

    ! LINALG_CHOLESKY.F90
    public :: cholesky_factor
    public :: cholesky_rank1_update
    public :: cholesky_rank1_downdate
    public :: solve_cholesky

    ! LINALG_LQ.F90
    public :: lq_factor
    public :: form_lq
    public :: mult_lq
    public :: solve_lq

    ! LINALG_SVD.F90
    public :: svd

    ! LINALG_INVERSE.F90
    public :: mtx_inverse
    public :: mtx_pinverse

    ! LINALG_LEAST_SQUARES.F90
    public :: solve_least_squares
    public :: solve_least_squares_full
    public :: solve_least_squares_svd

end module
