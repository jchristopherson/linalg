# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the QRUPDATE source files
set(qrupdate_sources
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
list(APPEND LINALG_SOURCES ${qrupdate_sources})
set(LINALG_SOURCES ${LINALG_SOURCES} PARENT_SCOPE)
