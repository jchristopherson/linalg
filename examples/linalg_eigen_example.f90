! linalg_eigen_example.f90

! This is an example illustrating the use of the eigenvalue and eigenvector 
! routines to solve the problem a free vibration problem of 3 masses connected
! by springs.
!
!     k1           k2           k3           k4
! |-\/\/\-| m1 |-\/\/\-| m2 |-\/\/\-| m3 |-\/\/\-|
!
! As illustrated above, the system consists of 3 masses connected by springs.
! Spring k1 and spring k4 connect the end masses to ground.  The equations of 
! motion for this system are as follows.
!
! | m1  0   0 | |x1"|   | k1+k2  -k2      0  | |x1|   |0|
! | 0   m2  0 | |x2"| + |  -k2  k2+k3    -k3 | |x2| = |0|
! | 0   0   m3| |x3"|   |   0    -k3    k3+k4| |x3|   |0|
!
! Notice: x1" = the second time derivative of x1.
program example
    use linalg_constants, only : dp, i32
    use linalg_eigen
    implicit none

    ! Define the model parameters
    real(dp), parameter :: pi = 3.14159265359d0
    real(dp), parameter :: m1 = 0.5d0
    real(dp), parameter :: m2 = 2.5d0
    real(dp), parameter :: m3 = 0.75d0
    real(dp), parameter :: k1 = 5.0d6
    real(dp), parameter :: k2 = 10.0d6
    real(dp), parameter :: k3 = 10.0d6
    real(dp), parameter :: k4 = 5.0d6

    ! Local Variables
    integer(i32) :: i, j
    real(dp) :: m(3,3), k(3,3), natFreq(3)
    complex(dp) :: vals(3), modeShapes(3,3)

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])
    
    ! Compute the eigenvalues and eigenvectors.
    call eigen(k, m, vals, vecs = modeShapes)

    ! Compute the natural frequency values, and return them with units of Hz.  
    ! Notice, all eigenvalues and eigenvectors are real for this example.
    natFreq = sqrt(real(vals)) / (2.0d0 * pi)

    ! Display the natural frequency and mode shape values.  Notice, the eigen
    ! routine does not necessarily sort the values.
    print '(A)', "Modal Information (Not Sorted):"
    do i = 1, size(natFreq)
        print '(AI0AF8.4A)', "Mode ", i, ": (", natFreq(i), " Hz)"
        print '(F10.3)', (real(modeShapes(j,i)), j = 1, size(natFreq))
    end do
end program