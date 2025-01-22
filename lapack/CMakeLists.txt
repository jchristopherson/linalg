include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    lapack
    GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
)
FetchContent_MakeAvailable(lapack)

set(BLAS_FOUND ${BLAS_FOUND} PARENT_SCOPE)
set(BLAS_LIBRARIES ${BLAS_LIBRARIES} PARENT_SCOPE)
set(LAPACK_FOUND ${LAPACK_FOUND} PARENT_SCOPE)
set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} PARENT_SCOPE)