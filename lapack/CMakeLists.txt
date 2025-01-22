include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    lapack
    GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
)
FetchContent_MakeAvailable(lapack)