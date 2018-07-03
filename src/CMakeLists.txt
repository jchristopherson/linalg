# Locate BLAS
option(USE_EXISTING_BLAS "Use an existing BLAS instead of the default NETLIB BLAS?" OFF)
if (USE_EXISTING_BLAS)
    message(STATUS "Using an existing BLAS.")
    find_package(blas)
    set(BLAS_LIBRARY blas)
else()
    message(STATUS "Using the default NETLIB BLAS.")
    set(blas_sources
        external/lapack/BLAS/SRC/dasum.f
        external/lapack/BLAS/SRC/daxpy.f
        external/lapack/BLAS/SRC/dcopy.f
        external/lapack/BLAS/SRC/ddot.f
        external/lapack/BLAS/SRC/dgbmv.f
        external/lapack/BLAS/SRC/dgemm.f
        external/lapack/BLAS/SRC/dgemv.f
        external/lapack/BLAS/SRC/dger.f
        external/lapack/BLAS/SRC/dnrm2.f
        external/lapack/BLAS/SRC/drot.f
        external/lapack/BLAS/SRC/drotg.f
        external/lapack/BLAS/SRC/drotm.f
        external/lapack/BLAS/SRC/drotmg.f
        external/lapack/BLAS/SRC/dsbmv.f
        external/lapack/BLAS/SRC/dscal.f
        external/lapack/BLAS/SRC/dspmv.f
        external/lapack/BLAS/SRC/dspr.f
        external/lapack/BLAS/SRC/dspr2.f
        external/lapack/BLAS/SRC/dswap.f
        external/lapack/BLAS/SRC/dsymm.f
        external/lapack/BLAS/SRC/dsymv.f
        external/lapack/BLAS/SRC/dsyr.f
        external/lapack/BLAS/SRC/dsyr2.f
        external/lapack/BLAS/SRC/dsyr2k.f
        external/lapack/BLAS/SRC/dsyrk.f
        external/lapack/BLAS/SRC/dtbmv.f
        external/lapack/BLAS/SRC/dtbsv.f
        external/lapack/BLAS/SRC/dtpmv.f
        external/lapack/BLAS/SRC/dtpsv.f
        external/lapack/BLAS/SRC/dtrmm.f
        external/lapack/BLAS/SRC/dtrmv.f
        external/lapack/BLAS/SRC/dtrsm.f
        external/lapack/BLAS/SRC/dtrsv.f
        external/lapack/BLAS/SRC/idamax.f
        external/lapack/BLAS/SRC/lsame.f
        external/lapack/BLAS/SRC/xerbla.f
        external/lapack/BLAS/SRC/dcabs1.f
        external/lapack/BLAS/SRC/dzasum.f
        external/lapack/BLAS/SRC/dznrm2.f
        external/lapack/BLAS/SRC/izamax.f
        external/lapack/BLAS/SRC/zaxpy.f
        external/lapack/BLAS/SRC/zcopy.f
        external/lapack/BLAS/SRC/zdotc.f
        external/lapack/BLAS/SRC/zdotu.f
        external/lapack/BLAS/SRC/zdscal.f
        external/lapack/BLAS/SRC/zrotg.f
        external/lapack/BLAS/SRC/zscal.f
        external/lapack/BLAS/SRC/zswap.f
        external/lapack/BLAS/SRC/zdrot.f
        external/lapack/BLAS/SRC/zgemv.f
        external/lapack/BLAS/SRC/zgbmv.f
        external/lapack/BLAS/SRC/zhemv.f
        external/lapack/BLAS/SRC/zhbmv.f
        external/lapack/BLAS/SRC/zhpmv.f
        external/lapack/BLAS/SRC/ztrmv.f
        external/lapack/BLAS/SRC/ztbmv.f
        external/lapack/BLAS/SRC/ztpmv.f
        external/lapack/BLAS/SRC/ztrsv.f
        external/lapack/BLAS/SRC/ztbsv.f
        external/lapack/BLAS/SRC/ztpsv.f
        external/lapack/BLAS/SRC/zgerc.f
        external/lapack/BLAS/SRC/zgeru.f
        external/lapack/BLAS/SRC/zher.f
        external/lapack/BLAS/SRC/zhpr.f
        external/lapack/BLAS/SRC/zher2.f
        external/lapack/BLAS/SRC/zhpr2.f
        external/lapack/BLAS/SRC/zgemm.f
        external/lapack/BLAS/SRC/zsymm.f
        external/lapack/BLAS/SRC/zsyrk.f
        external/lapack/BLAS/SRC/zsyr2k.f
        external/lapack/BLAS/SRC/ztrmm.f
        external/lapack/BLAS/SRC/ztrsm.f
        external/lapack/BLAS/SRC/zhemm.f
        external/lapack/BLAS/SRC/zherk.f
        external/lapack/BLAS/SRC/zher2k.f
    )
    add_library(blas ${blas_sources})
    set(BLAS_LIBRARY blas)
endif()

# Include only the necessary LAPACK files
set(lapack_sources
    external/lapack/SRC/dbdsqr.f
    external/lapack/SRC/dgbtf2.f
    external/lapack/SRC/dgbtrf.f
    external/lapack/SRC/dgbtrs.f
    external/lapack/SRC/dgebak.f
    external/lapack/SRC/dgebal.f
    external/lapack/SRC/dgebd2.f
    external/lapack/SRC/dgebrd.f
    external/lapack/SRC/dgecon.f
    external/lapack/SRC/dgeev.f
    external/lapack/SRC/dgehd2.f
    external/lapack/SRC/dgehrd.f
    external/lapack/SRC/dgelq2.f
    external/lapack/SRC/dgelqf.f
    external/lapack/SRC/dgels.f
    external/lapack/SRC/dgelss.f
    external/lapack/SRC/dgelsy.f
    external/lapack/SRC/dgeqp3.f
    external/lapack/SRC/dgeqr2.f
    external/lapack/SRC/dgeqrf.f
    external/lapack/SRC/dgesvd.f
    external/lapack/SRC/dgetf2.f
    external/lapack/SRC/dgetrf.f
    external/lapack/SRC/dgetri.f
    external/lapack/SRC/dgetrs.f
    external/lapack/SRC/dggbak.f
    external/lapack/SRC/dggbal.f
    external/lapack/SRC/dggev.f
    external/lapack/SRC/dgghrd.f
    external/lapack/SRC/dhgeqz.f
    external/lapack/SRC/dhseqr.f
    external/lapack/SRC/disnan.f
    external/lapack/SRC/dlabad.f
    external/lapack/SRC/dlabrd.f
    external/lapack/SRC/dlacn2.f
    external/lapack/SRC/dlacpy.f
    external/lapack/SRC/dladiv.f
    external/lapack/SRC/dlae2.f
    external/lapack/SRC/dlaev2.f
    external/lapack/SRC/dlaexc.f
    external/lapack/SRC/dlag2.f
    external/lapack/SRC/dlahqr.f
    external/lapack/SRC/dlahr2.f
    external/lapack/SRC/dlaic1.f
    external/lapack/SRC/dlaisnan.f
    external/lapack/SRC/dlaln2.f
    external/lapack/INSTALL/dlamch.f
    external/lapack/SRC/dlange.f
    external/lapack/SRC/dlanhs.f
    external/lapack/SRC/dlanst.f
    external/lapack/SRC/dlansy.f
    external/lapack/SRC/dlanv2.f
    external/lapack/SRC/dlapy2.f
    external/lapack/SRC/dlapy3.f
    external/lapack/SRC/dlaqp2.f
    external/lapack/SRC/dlaqps.f
    external/lapack/SRC/dlaqr0.f
    external/lapack/SRC/dlaqr1.f
    external/lapack/SRC/dlaqr2.f
    external/lapack/SRC/dlaqr3.f
    external/lapack/SRC/dlaqr4.f
    external/lapack/SRC/dlaqr5.f
    external/lapack/SRC/dlarf.f
    external/lapack/SRC/dlarfb.f
    external/lapack/SRC/dlarfg.f
    external/lapack/SRC/dlarft.f
    external/lapack/SRC/dlarfx.f
    external/lapack/SRC/dlartg.f
    external/lapack/SRC/dlarz.f
    external/lapack/SRC/dlarzb.f
    external/lapack/SRC/dlarzt.f
    external/lapack/SRC/dlas2.f
    external/lapack/SRC/dlascl.f
    external/lapack/SRC/dlaset.f
    external/lapack/SRC/dlasq1.f
    external/lapack/SRC/dlasq2.f
    external/lapack/SRC/dlasq3.f
    external/lapack/SRC/dlasq4.f
    external/lapack/SRC/dlasq5.f
    external/lapack/SRC/dlasq6.f
    external/lapack/SRC/dlasr.f
    external/lapack/SRC/dlasrt.f
    external/lapack/SRC/dlassq.f
    external/lapack/SRC/dlasv2.f
    external/lapack/SRC/dlaswp.f
    external/lapack/SRC/dlasy2.f
    external/lapack/SRC/dlatrd.f
    external/lapack/SRC/dlatrs.f
    external/lapack/SRC/dlatrz.f
    external/lapack/SRC/dorg2l.f
    external/lapack/SRC/dorg2r.f
    external/lapack/SRC/dorgbr.f
    external/lapack/SRC/dorghr.f
    external/lapack/SRC/dorgl2.f
    external/lapack/SRC/dorglq.f
    external/lapack/SRC/dorgql.f
    external/lapack/SRC/dorgqr.f
    external/lapack/SRC/dorgtr.f
    external/lapack/SRC/dorm2r.f
    external/lapack/SRC/dormbr.f
    external/lapack/SRC/dormhr.f
    external/lapack/SRC/dorml2.f
    external/lapack/SRC/dormlq.f
    external/lapack/SRC/dormqr.f
    external/lapack/SRC/dormr3.f
    external/lapack/SRC/dormrz.f
    external/lapack/SRC/dpotf2.f
    external/lapack/SRC/dpotrf.f
    external/lapack/SRC/dpotrf2.f
    external/lapack/SRC/dpotrs.f
    external/lapack/SRC/drscl.f
    external/lapack/SRC/dsteqr.f
    external/lapack/SRC/dsterf.f
    external/lapack/SRC/dsyev.f
    external/lapack/SRC/dsytd2.f
    external/lapack/SRC/dsytrd.f
    external/lapack/SRC/dtgevc.f
    external/lapack/SRC/dtrevc.f
    external/lapack/SRC/dtrexc.f
    external/lapack/SRC/dtrti2.f
    external/lapack/SRC/dtrtri.f
    external/lapack/SRC/dtrtrs.f
    external/lapack/SRC/dtzrzf.f
    external/lapack/SRC/ieeeck.f
    external/lapack/SRC/iladlc.f
    external/lapack/SRC/iladlr.f
    external/lapack/SRC/ilaenv.f
    external/lapack/SRC/iparmq.f
    external/lapack/SRC/dtrevc3.f
    external/lapack/SRC/dgetrf2.f
    external/lapack/SRC/zgetrf.f
    external/lapack/SRC/zgetf2.f
    external/lapack/SRC/zlaswp.f
    external/lapack/SRC/zgetrf2.f
    external/lapack/SRC/zgetrs.f
)

# Define the QRUPDATE source files
set(qrupdate_sources
    external/qrupdate/dch1dn.f
    external/qrupdate/dch1up.f
    external/qrupdate/dlup1up.f
    external/qrupdate/dqhqr.f
    external/qrupdate/dqr1up.f
    external/qrupdate/dqrot.f
    external/qrupdate/dqrqh.f
    external/qrupdate/dqrtv1.f
)

# Define the LINALG source files
set(linalg_sources
    linalg_core.f90
    linalg_constants.f90
    linalg_eigen.f90
    linalg_factor.f90
    linalg_solve.f90
    linalg_sorting.f90
    linalg_basic.f90
    linalg_immutable.f90
)
list(APPEND linalg_sources ${qrupdate_sources})
list(APPEND linalg_sources ${lapack_sources})

# Build the library
add_library(linalg ${linalg_sources})
target_link_libraries(linalg ${BLAS_LIBRARY} ${FERROR_LIBRARIES})

# https://stackoverflow.com/questions/30987311/exporting-an-imported-library

# ------------------------------------------------------------------------------
# INSTALLATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Define target information
set_property(TARGET linalg PROPERTY VERSION ${linalg_VERSION})
set_property(TARGET linalg PROPERTY SOVERSION ${linalg_VERSION_PATCH})
set_property(TARGET linalg PROPERTY INTERFACE_linalg_MAJOR_VERSION ${linalg_VERSION_MINOR})
set_property(TARGET linalg APPEND PROPERTY COMPATIBLE_INTERFACE_STRING linalg_VERSION_MAJOR)

# Locate the "include" directories
set(linalg_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

# Define the installation instructions
install(TARGETS linalg ${BLAS_LIBRARY}
    EXPORT linalgTargets
    RUNTIME DESTINATION linalg/bin
    LIBRARY DESTINATION linalg/lib
    ARCHIVE DESTINATION linalg/lib
    CONFIGURATIONS Release
    INCLUDES DESTINATION linalg/include
)
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}/linalg
    COMPONENT Devel
)

# Include the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION linalg/doc)
install(FILES ${PROJECT_SOURCE_DIR}/doc/refman.pdf DESTINATION linalg/doc)

# Define the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/linalgConfigVersion.cmake"
    VERSION ${linalg_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT linalgTargets
    FILE "${CMAKE_BINARY_DIR}/linalgTargets.cmake"
)

# Define the configuration file
configure_file(
    "${PROJECT_SOURCE_DIR}/linalgConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/linalgConfig.cmake"
    COPYONLY
)

set(ConfigPackageLocation linalg/lib/cmake/linalg)
install(
    EXPORT linalgTargets
    FILE linalgTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)
install(
    FILES
        "${CMAKE_BINARY_DIR}/linalgConfig.cmake"
        "${CMAKE_BINARY_DIR}/linalgConfigVersion.cmake"
    DESTINATION
        ${ConfigPackageLocation}
    COMPONENT
        Devel
)
