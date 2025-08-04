# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the QRUPDATE source files
set(QRUPDATE_SOURCES
    ${dir}/qrupdate/caxcpy.f
    ${dir}/qrupdate/cch1dn.f
    ${dir}/qrupdate/cch1up.f
    ${dir}/qrupdate/cchdex.f
    ${dir}/qrupdate/cchinx.f
    ${dir}/qrupdate/cchshx.f
    ${dir}/qrupdate/cgqvec.f
    ${dir}/qrupdate/clu1up.f
    ${dir}/qrupdate/clup1up.f
    ${dir}/qrupdate/cqhqr.f
    ${dir}/qrupdate/cqr1up.f
    ${dir}/qrupdate/cqrdec.f
    ${dir}/qrupdate/cqrder.f
    ${dir}/qrupdate/cqrinc.f
    ${dir}/qrupdate/cqrinr.f
    ${dir}/qrupdate/cqrot.f
    ${dir}/qrupdate/cqrqh.f
    ${dir}/qrupdate/cqrshc.f
    ${dir}/qrupdate/cqrtv1.f
    ${dir}/qrupdate/dch1dn.f
    ${dir}/qrupdate/dch1up.f
    ${dir}/qrupdate/dchdex.f
    ${dir}/qrupdate/dchinx.f
    ${dir}/qrupdate/dchshx.f
    ${dir}/qrupdate/dgqvec.f
    ${dir}/qrupdate/dlu1up.f
    ${dir}/qrupdate/dlup1up.f
    ${dir}/qrupdate/dqhqr.f
    ${dir}/qrupdate/dqr1up.f
    ${dir}/qrupdate/dqrdec.f
    ${dir}/qrupdate/dqrder.f
    ${dir}/qrupdate/dqrinc.f
    ${dir}/qrupdate/dqrinr.f
    ${dir}/qrupdate/dqrot.f
    ${dir}/qrupdate/dqrqh.f
    ${dir}/qrupdate/dqrshc.f
    ${dir}/qrupdate/dqrtv1.f
    ${dir}/qrupdate/sch1dn.f
    ${dir}/qrupdate/sch1up.f
    ${dir}/qrupdate/schdex.f
    ${dir}/qrupdate/schinx.f
    ${dir}/qrupdate/schshx.f
    ${dir}/qrupdate/sgqvec.f
    ${dir}/qrupdate/slu1up.f
    ${dir}/qrupdate/slup1up.f
    ${dir}/qrupdate/sqhqr.f
    ${dir}/qrupdate/sqr1up.f
    ${dir}/qrupdate/sqrdec.f
    ${dir}/qrupdate/sqrder.f
    ${dir}/qrupdate/sqrinc.f
    ${dir}/qrupdate/sqrinr.f
    ${dir}/qrupdate/sqrot.f
    ${dir}/qrupdate/sqrqh.f
    ${dir}/qrupdate/sqrshc.f
    ${dir}/qrupdate/sqrtv1.f
    ${dir}/qrupdate/zaxcpy.f
    ${dir}/qrupdate/zch1dn.f
    ${dir}/qrupdate/zch1up.f
    ${dir}/qrupdate/zchdex.f
    ${dir}/qrupdate/zchinx.f
    ${dir}/qrupdate/zchshx.f
    ${dir}/qrupdate/zgqvec.f
    ${dir}/qrupdate/zlu1up.f
    ${dir}/qrupdate/zlup1up.f
    ${dir}/qrupdate/zqhqr.f
    ${dir}/qrupdate/zqr1up.f
    ${dir}/qrupdate/zqrdec.f
    ${dir}/qrupdate/zqrder.f
    ${dir}/qrupdate/zqrinc.f
    ${dir}/qrupdate/zqrinr.f
    ${dir}/qrupdate/zqrot.f
    ${dir}/qrupdate/zqrqh.f
    ${dir}/qrupdate/zqrshc.f
    ${dir}/qrupdate/zqrtv1.f
)

# Define the SPARSKIT source files
set(SPARSEKIT_SOURCES
    ${dir}/sparskit2/blassm.f
    ${dir}/sparskit2/matvec.f
    ${dir}/sparskit2/formats.f
    ${dir}/sparskit2/unary.f
    ${dir}/sparskit2/iters.f
    ${dir}/sparskit2/distdot.f
    ${dir}/sparskit2/ilut.f
)

# Define the LINALG source files
set(LINALG_SOURCES
    ${dir}/linalg.f90
    ${dir}/linalg_eigen.f90
    ${dir}/linalg_sorting.f90
    ${dir}/linalg_basic.f90
    ${dir}/blas.f90
    ${dir}/lapack.f90
    ${dir}/qrupdate.f90
    ${dir}/sparskit.f90
    ${dir}/linalg_sparse.f90
    ${dir}/linalg_errors.f90
    ${dir}/linalg_lu.f90
    ${dir}/linalg_rz.f90
    ${dir}/linalg_qr.f90
    ${dir}/linalg_tri.f90
    ${dir}/linalg_cholesky.f90
    ${dir}/linalg_lq.f90
    ${dir}/linalg_svd.f90
    ${dir}/linalg_inverse.f90
    ${dir}/linalg_least_squares.f90
)
if (${BUILD_LINALG_C_API})
    set(LINALG_C_SOURCES
        ${dir}/linalg_c_api.f90
    )
    list(APPEND LINALG_SOURCES ${LINALG_C_SOURCES})
endif()
list(APPEND LINALG_SOURCES ${QRUPDATE_SOURCES})
list(APPEND LINALG_SOURCES ${SPARSEKIT_SOURCES})
set(LINALG_SOURCES ${LINALG_SOURCES} PARENT_SCOPE)

include(FetchContent)

find_package(BLAS)
find_package(LAPACK)
find_package(ferror QUIET)

if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(STATUS "BLAS/LAPACK could not be found.  A reference version will be employed.")
  FetchContent_Declare(
      lapack
      GIT_REPOSITORY "https://github.com/Reference-LAPACK/lapack"
  )
  FetchContent_MakeAvailable(lapack)
  set(BLAS_LIBRARIES blas)
  set(LAPACK_LIBRARIES lapack)
endif()

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
else()
    set(FERROR_LIBRARIES ferror::ferror)
endif()

add_library(${PROJECT_NAME} ${LINALG_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${FERROR_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(LINALG_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${LINALG_MOD_DIR}")
    file(MAKE_DIRECTORY "${LINALG_MOD_DIR}")
endif()


set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${LINALG_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${LINALG_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

install(
    FILES "${PROJECT_SOURCE_DIR}/include/ferror.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(DIRECTORY ${LINALG_MOD_DIR} DESTINATION "${CMAKE_INSTALL_MODULEDIR}")
