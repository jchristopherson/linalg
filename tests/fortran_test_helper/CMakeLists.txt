include(FetchContent)
FetchContent_Declare(
    fortran_test_helper
    GIT_TAG "origin/main"
    GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(fortran_test_helper)
set(fortran_test_helper_INCLUDE_DIR ${fortran_test_helper_BINARY_DIR}/include)
set(fortran_test_helper_LIBRARY fortran_test_helper)

set(fortran_test_helper_INCLUDE_DIR ${fortran_test_helper_INCLUDE_DIR} PARENT_SCOPE)
set(fortran_test_helper_LIBRARY ${fortran_test_helper_LIBRARY} PARENT_SCOPE)