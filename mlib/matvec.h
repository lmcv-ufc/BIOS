// ------------------------------------------------------------------------
// matvec.h - prototypes of functions to handle matrices and vectors.
// -------------------------------------------------------------------------
// Copyright (c) 2021 LMCV/UFC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------------------------------------------------
//
// ------------------------------------------------------------------------
// Public functions:
// ------------------------------------------------------------------------
//
// double *VecAlloc(int n)
//
//   n - vector dimension                                             (in)
//
// This function allocates memory for a vector of the given size.
// ------------------------------------------------------------------------
//
// void VecFree(double *v)
//
//   v - given vector                                                 (in)
//
// This function releases the memory allocated to the given vector.
// ------------------------------------------------------------------------
//
// void VecZero(int n, double *v)
//
//   n - vector dimension                                             (in)
//   v - given vector                                             (in/out)
//
// This function assigns 0.0 to all elements of the given vector.
// ------------------------------------------------------------------------
//
// void VecAssign(int n, double *u, double *v)
//
//   n - vector dimension                                             (in)
//   u - input vector                                                 (in)
//   v - output vector                                               (out)
//
// This function copies a given vector to another vector: {v} = {u}.
// ------------------------------------------------------------------------
//
// void VecMult(int n, double a, double *u, double *v)
//
//   n - vector dimension                                             (in)
//   a - given scalar                                                 (in)
//   u - input vector                                                 (in)
//   v - output vector                                               (out)
//
// This function multiplies a given vector by a scalar: {v} = a*{u}.
// ------------------------------------------------------------------------
//
// void VecAdd(int n, double *u, double *v)
//
//   n - vector dimension                                             (in)
//   u - input vector                                                 (in)
//   v - output vector                                               (out)
//
// This function adds a given vector to another vector: {v} += {u}.
// ------------------------------------------------------------------------
//
// void VecAdd(int n, double *u, double *v, double *w)
//
//   n   - vector dimension                                           (in)
//   u,v - input vectors                                              (in)
//   w   - output vector                                             (out)
//
// This function performs the operation {w} = {u} + {v}.
// ------------------------------------------------------------------------
//
// void VecSclAdd(int n, double a, double *u, double *v, double *w)
//
//   n   - vector dimension                                           (in)
//   a   - given scalar                                               (in)
//   u,v - input vectors                                              (in)
//   w   - output vector                                             (out)
//
// This function performs the operation {w} = a*{u} + {v}.
// ------------------------------------------------------------------------
//
// void VecSub(int n, double *u, double *v)
//
//   n - vector dimension                                             (in)
//   u - input vector                                                 (in)
//   v - output vector                                               (out)
//
// This function performs the operation {v} -= {u}.
// ------------------------------------------------------------------------
//
// void VecSub(int n, double *u, double *v, double *w)
//
//   n   - vector dimension                                           (in)
//   u,v - input vectors                                              (in)
//   w   - output vector                                             (out)
//
// This function performs the operation {w} = {u} - {v}.
// ------------------------------------------------------------------------
//
// void VecSclSub(int n, double a, double *u, double *v, double *w)
//
//   n   - vector dimension                                           (in)
//   a   - given scalar                                               (in)
//   u,v - input vectors                                              (in)
//   w   - output vector                                             (out)
//
// This function performs the operation {w} = a*{u} - {v}.
// ------------------------------------------------------------------------
//
// void VecDiv(int n, double *u, double *v, double *w)
//
//   n   - vector dimension                                           (in)
//   u,v - input vectors                                              (in)
//   w   - output vector                                             (out)
//
// This function performs the operation {w} = {u}/{v} => w[i] = u[i]/v[i].
// ------------------------------------------------------------------------
//
// double VecDot(int n, double *u, double *v)
//
//   n   - vector dimension                                           (in)
//   u,v - given vectors                                              (in)
//
// This function returns the scalar (dot) product between two vectors.
// ------------------------------------------------------------------------
//
// void VecCrossProd(double *u, double *v, double *w)
//
//   u,v - given vectors                                              (in)
//   w   - result vector {w} = {u}x{v}                               (out)
//
// This function evaluates the cross product between two vectors.
// ------------------------------------------------------------------------
//
// double VecLen(int n, double *u)
//
//   n - vector dimension                                             (in)
//   u - given vector                                                 (in)
//
// This function returns the length of the given vector.
// ------------------------------------------------------------------------
//
// double VecNormInf(int n, double *u)
//
//   n - vector dimension                                             (in)
//   u - given vector                                                 (in)
//
// This function returns the infinite norm (max|u(i)|) of the given vector.
// ------------------------------------------------------------------------
//
// double **MatAlloc(int n, int m)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//
// This function allocates memory for a matrix of the given dimensions.
// ------------------------------------------------------------------------
//
// void MatZero(int n, int m, double **A)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   A - given matrix                                             (in/out)
//
// This function assigns 0.0 to all elements of the given matrix.
// ------------------------------------------------------------------------
//
// void MatAssign(int n, int m, double **A, double **B)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   A - given matrix                                                 (in)
//   B - output matrix                                               (out)
//
// This function performs the operation [B] = [A].
// ------------------------------------------------------------------------
//
// void MatAdd(int n, int m, double **A, double **B)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   A - given matrix                                                 (in)
//   B - output matrix                                               (out)
//
// This function performs the operation [B] += [A].
// ------------------------------------------------------------------------
//
// void MatAdd(int n, int m, double **A, double **B, double **C)
//
//   n   - number of rows                                             (in)
//   m   - number of colums                                           (in)
//   A,B - given matrices                                             (in)
//   C   - output matrix                                             (out)
//
// This function performs the operation [C] = [A] + [B].
// ------------------------------------------------------------------------
//
// void MatSub(int n, int m, double **A, double **B)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   A - given matrix                                                 (in)
//   B - output matrix                                               (out)
//
// This function performs the operation [B] -= [A].
// ------------------------------------------------------------------------
//
// void MatSub(int n, int m, double **A, double **B, double **C)
//
//   n   - number of rows                                             (in)
//   m   - number of colums                                           (in)
//   A,B - given matrices                                             (in)
//   C   - output matrix                                             (out)
//
// This function performs the operation [C] = [A] - [B].
// ------------------------------------------------------------------------
//
// void MatMult(int n, int m, double f, double **A, double **B)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   f - given scalar                                                 (in)
//   A - given matrix                                                 (in)
//   B - output matrix                                               (out)
//
// This function performs the operation [B] = f*[A].
// ------------------------------------------------------------------------
//
// void MatMult(int n, int m, int l, double **A, double **B, double **C)
//
//   n,m,l - dimension of given matrices                              (in)
//   A - left matrix (n x l)                                          (in)
//   B - right matrix (l x m)                                         (in)
//   C - output matrix (n x m)                                       (out)
//
// This function performs the operation [C] = [A] * [B].
//                                        nxm   nxl   lxm
// ------------------------------------------------------------------------
//
// void MatVecMult(int n, int m, double **A, double *u, double *v)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   A - given matrix                                                 (in)
//   u - given vector                                                 (in)
//   v - output vector                                               (out)
//
// This function performs the operation {v} = [A]*{u}.
// ------------------------------------------------------------------------
//
// void MatVecMultAcc(int n, int m, double **A, double *u, double *v)
//
//   n - number of rows                                               (in)
//   m - number of colums                                             (in)
//   A - given matrix                                                 (in)
//   u - given vector                                                 (in)
//   v - output vector                                            (in/out)
//
// This function performs the operation {v} += [A]*{u}.
// ------------------------------------------------------------------------
//
// void MatVecMultT(int m, int n, double **A, double *u, double *v)
//
//   m - number of colums                                             (in)
//   n - number of rows                                               (in)
//   A - given matrix                                                 (in)
//   u - given vector                                                 (in)
//   v - output vector                                               (out)
//
// This function performs the operation {v} = [A]t*{u}.
// ------------------------------------------------------------------------
//
// void MatVecMultTAcc(double a, int m, int n, double **A, double *u,
//                     double *v)
//
//   a - given scalar                                                 (in)
//   m - number of colums                                             (in)
//   n - number of rows                                               (in)
//   A - given matrix                                                 (in)
//   u - given vector                                                 (in)
//   v - output vector                                            (in/out)
//
// This function performs the operation {v} += a*[A]t*{u}.
// ------------------------------------------------------------------------
//
// void MatTripBtCB(int m, int n, double **B, double **C, double w,
//                  double **K)
//
//   m,n - matrix dimensions                                          (in)
//   B   - rectangular matrix => m x n                                (in)
//   C   - square matrix => m x m                                     (in)
//   w   - scalar                                                     (in)
//   K   - output matrix                                             (out)
//
// This function performs the operation [K] += w*[B]t*[C]*[B].
// ------------------------------------------------------------------------
//
// double **SkylAlloc(int n, int *p)
//
//   n - number of rows                                               (in)
//   p - matrix profile                                               (in)
//
// This function allocates a symmetric square matrix stored in skyline
// form. The profile vector stores the column of the first non-zero
// element of each row of the matrix.
// ------------------------------------------------------------------------
//
// void SkylFree(double **A, int n, int *p)
//
//   A - given matrix                                                 (in)
//   n - number of rows                                               (in)
//   p - matrix profile                                               (in)
//
// This function releases the memory allocated to the given profile matrix.
// ------------------------------------------------------------------------
//
// void SkylZero(int n, int *p, double **A)
//
//   n - number of rows                                               (in)
//   p - matrix profile                                               (in)
//   A - given matrix                                             (in/out)
//
// This function assigns 0.0 to all elements of the given matrix.
// ------------------------------------------------------------------------
//
// int CroutSolver(int op, double **A, double *b, int n, int *p)
//
//   op - (1) complete sol. / (2) factorization / (3) backsub.        (in)
//   A  - system matrix                                               (in)
//   b  - r.h.s. vector / solution vector                         (in/out)
//   n  - system dimension                                            (in)
//   p  - matrix skyline                                              (in)
//
// This function solves the system [A]{x} = {b} using the [L][D][L]t
// factorization. It returns (1) for a successful operation and (0)
// otherwise.
// ------------------------------------------------------------------------

#ifndef _MATVEC_H
#define _MATVEC_H

#include <string.h>

// ------------------------------------------------------------------------
// Vector functions:
//
double *VecAlloc(int n);
void    VecFree(double *v);
void    VecZero(int n, double *u);
void    VecAssign(int n, double *u, double *v);
void    VecSet(int n, double a, double *v);
void    VecMul(int n, double a, double *u, double *v);
void    VecMulAdd(int n, double a, double *u, double *v);
void    VecAdd(int n, double *u, double *v);
void    VecAdd(int n, double *u, double *v, double *w);
void    VecSclAdd(int n, double a, double *u, double *v, double *w);
void    VecSclAddInc(int n, double a, double *u, double *v, double *w);
void    VecSub(int n, double *u, double *v);
void    VecSub(int n, double *u, double *v, double *w);
void    VecSclSub(int n, double a, double *u, double *v, double *w);
void    VecSclSubInc(int n, double a, double *u, double *v, double *w);
void    VecDiv(int n, double *u, double *v, double *w);
double  VecDot(int n, double *u, double *v);
void    VecCrossProd(double *u, double *v, double *w);
double  VecLen(int n, double *u);
double  VecNormInf(int n, double *u);

// ================================ VecZero ================================

inline void VecZero(int n, double *u)
{
  memset(u, '\0', n*sizeof(double));
}

// =============================== VecAssign ===============================

inline void VecAssign(int n, double *u, double *v)
{
  memcpy(v, u, n*sizeof(double));
}

// ------------------------------------------------------------------------
// Matrix functions:
//
double **MatAlloc(int n, int m);
void     MatFree(double **A, int n);
void     MatZero(int n, int m, double **A);
void     MatAssign(int n, int m, double **A, double **B);
void     MatAdd(int n, int m, double **A, double **B);
void     MatAdd(int n, int m, double **A, double **B, double **C);
void     MatSub(int n, int m, double **A, double **B);
void     MatSub(int n, int m, double **A, double **B, double **C);
void     MatMult(int n, int m, double f, double **A, double **B);
void     MatMult(int n, int m, int l, double **A, double **B, double **C);
void     MatVecMult(int n, int m, double **A, double *u, double *v);
void     MatVecMultAcc(int n, int m, double **A, double *u, double *v);
void     MatVecMultT(int m, int n, double **A, double *u, double *v);
void     MatVecMultTAcc(double a, int m, int n, double **A, double *u,
                        double *v);
void     MatTripBtCB(int m, int n, double **B, double **C, double w,
                     double **K);
int      MatDecompLU(int n, double **A);
void     MatSolveLU (int n, double **A, double *v);

// ------------------------------------------------------------------------
// Skyline matrix functions:
//
double **SkylAlloc(int n, int *p);
void     SkylFree(double **A, int n, int *p);
void     SkylZero(int n, int *p, double **A);
void     SkylAssign(int n, int *p, double **A, double **B);
void     SkylAdd(int n, int *p, double **A, double **B);
void     SkylAdd(double a, int n, int *p, double **A, double **B);
void     SkylAdd(int n, int *p, double **A, double **B, double **C);
void     SkylSub(int n, int *p, double **A, double **B);
void     SkylSub(int n, int *p, double **A, double **B, double **C);
void     SkylMult(double a, int n, int *p, double **A);
void     SkylMult(double a, int n, int *p, double **A, double **B);
void     SkylMultAcc(double a, int n, int *p, double **A, double **B);
void     SkylVecMult(int n, int *p, double **A, double *u, double *v);
int      CroutSolver(int op, double **A, double *b, int n, int *p);

#endif
