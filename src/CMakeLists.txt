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
if (${BUILD_C_API})
    set(LINALG_C_SOURCES
        ${dir}/linalg_c_api.f90
    )
    list(APPEND LINALG_SOURCES ${LINALG_C_SOURCES})
endif()
list(APPEND LINALG_SOURCES ${QRUPDATE_SOURCES})
list(APPEND LINALG_SOURCES ${SPARSEKIT_SOURCES})
set(LINALG_SOURCES ${LINALG_SOURCES} PARENT_SCOPE)