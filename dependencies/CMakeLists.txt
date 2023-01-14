# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Get BLAS
add_subdirectory(BLAS)

# Get LAPACK
add_subdirectory(LAPACK)

# Get FERROR
add_subdirectory(ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)

# Get QRUpdate
add_subdirectory(QRUpdate)
set(qrupdate_LIBRARY ${qrupdate_LIBRARY} PARENT_SCOPE)
