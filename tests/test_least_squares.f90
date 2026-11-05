! SPDX-License-Identifier: MIT
module test_least_squares
    use, intrinsic :: iso_fortran_env, only : real64
    use linalg
    use test_core
    use fortran_test_helper
    implicit none
contains
    function test_least_squares_shapes() result(rst)
        real(real64) :: tall(4,2), wide(2,3), b_tall(4), b_wide(2)
        real(real64) :: b_tall_many(4,2), b_wide_many(2,2)
        real(real64) :: x_tall(2), x_tall_many(2,2), x_wide(3), x_wide_many(3,2)
        real(real64) :: expected_real(2), expected_real_many(2,2)
        complex(real64) :: ctall(4,2), cwide(2,3), cb_tall(4), cb_wide(2)
        complex(real64) :: cb_tall_many(4,2), cb_wide_many(2,2)
        complex(real64) :: cx_tall(2), cx_tall_many(2,2), cx_wide(3), cx_wide_many(3,2)
        complex(real64) :: expected_tall(2), expected_tall_many(2,2)
        logical :: rst
        integer :: solver

        rst = .true.
        tall = reshape([1.0_real64, 0.0_real64, 1.0_real64, 2.0_real64, &
            0.0_real64, 1.0_real64, 1.0_real64, -1.0_real64], [4,2])
        expected_real = [2.0_real64, -1.0_real64]
        expected_real_many = reshape([2.0_real64, -1.0_real64, &
            -2.0_real64, 3.0_real64], [2,2])
        b_tall = matmul(tall, expected_real)
        b_tall_many = matmul(tall, expected_real_many)
        ! Real overdetermined systems: vector/matrix RHS across all interfaces.
        x_tall = solve_least_squares(tall, b_tall)
        x_tall_many = solve_least_squares(tall, b_tall_many)
        if (.not.assert(x_tall, expected_real, tol = REAL64_TOL) .or. &
            .not.assert(x_tall_many, expected_real_many, tol = REAL64_TOL)) rst = .false.
        x_tall = solve_least_squares_full(tall, b_tall)
        x_tall_many = solve_least_squares_full(tall, b_tall_many)
        if (.not.assert(x_tall, expected_real, tol = REAL64_TOL) .or. &
            .not.assert(x_tall_many, expected_real_many, tol = REAL64_TOL)) rst = .false.
        x_tall = solve_least_squares_svd(tall, b_tall)
        x_tall_many = solve_least_squares_svd(tall, b_tall_many)
        if (.not.assert(x_tall, expected_real, tol = REAL64_TOL) .or. &
            .not.assert(x_tall_many, expected_real_many, tol = REAL64_TOL)) rst = .false.

        ! Complex overdetermined systems exercise the corresponding overloads.
        ctall = cmplx(tall, 0.25_real64 * tall, kind = real64)
        expected_tall = [(1.0_real64, 1.0_real64), (2.0_real64, -1.0_real64)]
        expected_tall_many(:,1) = expected_tall
        expected_tall_many(:,2) = [(-2.0_real64, 1.0_real64), &
            (3.0_real64, 0.5_real64)]
        cb_tall = matmul(ctall, expected_tall)
        cb_tall_many = matmul(ctall, expected_tall_many)
        cx_tall = solve_least_squares(ctall, cb_tall)
        cx_tall_many = solve_least_squares(ctall, cb_tall_many)
        if (.not.assert(cx_tall, expected_tall, tol = REAL64_TOL) .or. &
            .not.assert(cx_tall_many, expected_tall_many, tol = REAL64_TOL)) then
            print '(A)', "Least-squares failed: complex tall DGELS"
            rst = .false.
        end if
        cx_tall = solve_least_squares_full(ctall, cb_tall)
        cx_tall_many = solve_least_squares_full(ctall, cb_tall_many)
        if (.not.assert(cx_tall, expected_tall, tol = REAL64_TOL) .or. &
            .not.assert(cx_tall_many, expected_tall_many, tol = REAL64_TOL)) then
            print '(A)', "Least-squares failed: complex tall ZGELSY"
            rst = .false.
        end if
        cx_tall = solve_least_squares_svd(ctall, cb_tall)
        cx_tall_many = solve_least_squares_svd(ctall, cb_tall_many)
        if (.not.assert(cx_tall, expected_tall, tol = REAL64_TOL) .or. &
            .not.assert(cx_tall_many, expected_tall_many, tol = REAL64_TOL)) then
            print '(A)', "Least-squares failed: complex tall ZGELSS"
            rst = .false.
        end if

        ! Underdetermined full-row-rank systems take a separate LAPACK branch.
        wide = reshape([1.0_real64, 0.0_real64, 0.0_real64, &
            1.0_real64, 1.0_real64, 1.0_real64], [2,3])
        b_wide = [2.0_real64, -1.0_real64]
        b_wide_many(:,1) = b_wide
        b_wide_many(:,2) = [-1.0_real64, 3.0_real64]
        do solver = 1, 3
            select case (solver)
            case (1)
                x_wide = solve_least_squares(wide, b_wide)
                x_wide_many = solve_least_squares(wide, b_wide_many)
            case (2)
                x_wide = solve_least_squares_full(wide, b_wide)
                x_wide_many = solve_least_squares_full(wide, b_wide_many)
            case (3)
                x_wide = solve_least_squares_svd(wide, b_wide)
                x_wide_many = solve_least_squares_svd(wide, b_wide_many)
            end select
            if (.not.assert(matmul(wide, x_wide), b_wide, tol = REAL64_TOL) .or. &
                .not.assert(matmul(wide, x_wide_many), b_wide_many, tol = REAL64_TOL)) then
                print '(A,I0)', "Least-squares failed: real wide solver ", solver
                rst = .false.
            end if
        end do

        cwide = cmplx(wide, 0.2_real64 * wide, kind = real64)
        cb_wide = [(2.0_real64, 1.0_real64), (-1.0_real64, 0.5_real64)]
        cb_wide_many(:,1) = cb_wide
        cb_wide_many(:,2) = [(-1.0_real64, 2.0_real64), (3.0_real64, -1.0_real64)]
        do solver = 1, 3
            select case (solver)
            case (1)
                cx_wide = solve_least_squares(cwide, cb_wide)
                cx_wide_many = solve_least_squares(cwide, cb_wide_many)
            case (2)
                cx_wide = solve_least_squares_full(cwide, cb_wide)
                cx_wide_many = solve_least_squares_full(cwide, cb_wide_many)
            case (3)
                cx_wide = solve_least_squares_svd(cwide, cb_wide)
                cx_wide_many = solve_least_squares_svd(cwide, cb_wide_many)
            end select
            if (.not.assert(matmul(cwide, cx_wide), cb_wide, tol = REAL64_TOL) .or. &
                .not.assert(matmul(cwide, cx_wide_many), cb_wide_many, tol = REAL64_TOL)) then
                print '(A,I0)', "Least-squares failed: complex wide solver ", solver
                rst = .false.
            end if
        end do

        if (.not.rst) print '(A)', "Test Failed: Least-Squares API Shapes"
    end function

    function test_least_squares_rank_deficient() result(rst)
        real(real64) :: a(3,2), b(3), expected(2), x(2)
        complex(real64) :: ac(3,2), cb(3), cx(2), cexpected(2)
        logical :: rst

        rst = .true.
        a = reshape([1.0_real64, 2.0_real64, 3.0_real64, &
            2.0_real64, 4.0_real64, 6.0_real64], [3,2])
        b = [1.0_real64, 2.0_real64, 3.0_real64]
        expected = [0.2_real64, 0.4_real64]
        x = solve_least_squares_svd(a, b)
        if (.not.assert(x, expected, tol = REAL64_TOL)) rst = .false.

        ac = cmplx(a, 0.3_real64 * a, kind = real64)
        cb = ac(:,1)
        cexpected = cmplx(expected, 0.0_real64, kind = real64)
        cx = solve_least_squares_svd(ac, cb)
        if (.not.assert(cx, cexpected, tol = REAL64_TOL)) rst = .false.

        if (.not.rst) print '(A)', "Test Failed: Rank-Deficient Least-Squares"
    end function
end module
