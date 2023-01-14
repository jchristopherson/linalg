# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the LINALG source files
set(LINALG_SOURCES
    ${dir}/linalg.f90
    ${dir}/linalg_eigen.f90
    ${dir}/linalg_factor.f90
    ${dir}/linalg_solve.f90
    ${dir}/linalg_sorting.f90
    ${dir}/linalg_basic.f90
)
if (${BUILD_C_API})
    set(LINALG_C_SOURCES
        ${dir}/linalg_c_api.f90
    )
    list(APPEND LINALG_SOURCES ${LINALG_C_SOURCES})
endif()

set(LINALG_SOURCES ${LINALG_SOURCES} PARENT_SCOPE)
