# Define the QRUPDATE source files
set(qrupdate_sources
    external/qrupdate/caxcpy.f
    external/qrupdate/cch1dn.f
    external/qrupdate/cch1up.f
    external/qrupdate/cchdex.f
    external/qrupdate/cchinx.f
    external/qrupdate/cchshx.f
    external/qrupdate/cgqvec.f
    external/qrupdate/clu1up.f
    external/qrupdate/clup1up.f
    external/qrupdate/cqhqr.f
    external/qrupdate/cqr1up.f
    external/qrupdate/cqrdec.f
    external/qrupdate/cqrder.f
    external/qrupdate/cqrinc.f
    external/qrupdate/cqrinr.f
    external/qrupdate/cqrot.f
    external/qrupdate/cqrqh.f
    external/qrupdate/cqrshc.f
    external/qrupdate/cqrtv1.f
    external/qrupdate/dch1dn.f
    external/qrupdate/dch1up.f
    external/qrupdate/dchdex.f
    external/qrupdate/dchinx.f
    external/qrupdate/dchshx.f
    external/qrupdate/dgqvec.f
    external/qrupdate/dlu1up.f
    external/qrupdate/dlup1up.f
    external/qrupdate/dqhqr.f
    external/qrupdate/dqr1up.f
    external/qrupdate/dqrdec.f
    external/qrupdate/dqrder.f
    external/qrupdate/dqrinc.f
    external/qrupdate/dqrinr.f
    external/qrupdate/dqrot.f
    external/qrupdate/dqrqh.f
    external/qrupdate/dqrshc.f
    external/qrupdate/dqrtv1.f
    external/qrupdate/sch1dn.f
    external/qrupdate/sch1up.f
    external/qrupdate/schdex.f
    external/qrupdate/schinx.f
    external/qrupdate/schshx.f
    external/qrupdate/sgqvec.f
    external/qrupdate/slu1up.f
    external/qrupdate/slup1up.f
    external/qrupdate/sqhqr.f
    external/qrupdate/sqr1up.f
    external/qrupdate/sqrdec.f
    external/qrupdate/sqrder.f
    external/qrupdate/sqrinc.f
    external/qrupdate/sqrinr.f
    external/qrupdate/sqrot.f
    external/qrupdate/sqrqh.f
    external/qrupdate/sqrshc.f
    external/qrupdate/sqrtv1.f
    external/qrupdate/zaxcpy.f
    external/qrupdate/zch1dn.f
    external/qrupdate/zch1up.f
    external/qrupdate/zchdex.f
    external/qrupdate/zchinx.f
    external/qrupdate/zchshx.f
    external/qrupdate/zgqvec.f
    external/qrupdate/zlu1up.f
    external/qrupdate/zlup1up.f
    external/qrupdate/zqhqr.f
    external/qrupdate/zqr1up.f
    external/qrupdate/zqrdec.f
    external/qrupdate/zqrder.f
    external/qrupdate/zqrinc.f
    external/qrupdate/zqrinr.f
    external/qrupdate/zqrot.f
    external/qrupdate/zqrqh.f
    external/qrupdate/zqrshc.f
    external/qrupdate/zqrtv1.f
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

# Build the library
add_library(linalg ${linalg_sources})
target_link_libraries(linalg ${ferror_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

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
install(TARGETS linalg
    EXPORT linalgTargets
    RUNTIME DESTINATION linalg/bin
    LIBRARY DESTINATION linalg/bin
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
