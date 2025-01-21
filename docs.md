---
project: LINALG
summary: LINALG is a linear algebra library that provides a user-friendly interface to several BLAS and LAPACK routines.
project_github: https://github.com/jchristopherson/linalg
author: Jason Christopherson
author_email: jchristopherson@hotmail.com
src_dir: ./src
exclude_dir: **/qrupdate
exclude_dir: **/sparskit2
exclude: **/blas.f90
exclude:  **/lapack.f90
exclude:  **/qrupdate.f90
exclude:  **/sparskit.f90
exclude: **/linalg_c_api.f90
output_dir: ./doc
display: public
source: true
proc_internals: true
sort: permission-alpha
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
---