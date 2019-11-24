# newmark_validation
This is a small validation program of newmark method to solve typical Finite Element Equations,like:
                                M(u'')+C(u')+K(u)=F
where, the matrix M,C and K are represented in Block Sparse Row (BSR) Compressed format or Compressed Sparse Row(CSR) format using 3 array
method.
  In Newmark method, a linear algbraic equations in the form of
                                  Ax = b
is solved by Intel Pardiso which is provided by Intel Math Kernel Library as a direct sparse solver based on LU decompositon and sparse matrix
storage.The Intel Pardiso solver is automatically speeded up using openMP.
