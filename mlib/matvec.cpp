// -------------------------------------------------------------------------
// matvec.cpp - functions to handle matrices and vectors.
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
// -------------------------------------------------------------------------
// Created:      17-Nov-2000     Evandro Parente Junior
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

#include "matvec.h"

const double TOLC = 1.0e-12;

// -------------------------------------------------------------------------
// Local functions:
//

// ================================== Max ==================================

inline int Max(int a, int b)
{
  if (a > b)
    return(a);
  else
    return(b);
}

// -------------------------------------------------------------------------
// Vector functions:
//

// =============================== VecAlloc ================================

double *VecAlloc(int n)
{
  double *v = new double[n];
  if (!v)
  {
    cout << "Allocation failure in VecAlloc !!!\n";
    exit(0);
  }

  return(v);
}

// =============================== VecFree ================================

void VecFree(double *v)
{
  delete v;
}

// ================================ VecSet =================================

void VecSet(int n, double a, double *v)
{
  for (int i = 0; i < n; i++) v[i] = a;
}

// ================================ VecMul =================================

void VecMul(int n, double a, double *u, double *v)
{
  for (int i = 0; i < n; i++) v[i] = a*u[i];
}

// ============================== VecMulAdd ================================

void VecMulAdd(int n, double a, double *u, double *v)
{
  for (int i = 0; i < n; i++) v[i] += a*u[i];
}

// ================================ VecAdd =================================

void VecAdd(int n, double *u, double *v)
{
  for (int i = 0; i < n; i++) v[i] += u[i];
}

// ================================ VecAdd =================================

void VecAdd(int n, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] = u[i] + v[i];
}

// =============================== VecSclAdd ===============================

void VecSclAdd(int n, double a, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] = a*u[i] + v[i];
}

// ============================= VecSclAddInc ==============================

void VecSclAddInc(int n, double a, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] += a*u[i] + v[i];
}

// ================================ VecSub =================================

void VecSub(int n, double *u, double *v)
{
  for (int i = 0; i < n; i++) v[i] -= u[i];
}

// ================================ VecSub =================================

void VecSub(int n, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] = u[i] - v[i];
}

// =============================== VecSclSub ===============================

void VecSclSub(int n, double a, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] = a*u[i] - v[i];
}

// ============================= VecSclSubInc ==============================

void VecSclSubInc(int n, double a, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] += a*u[i] - v[i];
}

// ================================ VecDiv =================================

void VecDiv(int n, double *u, double *v, double *w)
{
  for (int i = 0; i < n; i++) w[i] = u[i]/v[i];
}

// ================================ VecDot =================================

double VecDot(int n, double *u, double *v)
{
  double dot = 0.0;

  for (int i = 0; i < n; i++) dot += u[i]*v[i];

  return(dot);
}

// ============================ VecCrossProd ===============================

void VecCrossProd(double *u, double *v, double *w)
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

// ================================ VecLen =================================

double VecLen(int n, double *u)
{
  double len = sqrt(VecDot(n, u, u));

  return(len);
}
// ================================ VecLen =================================

double VecNormInf(int n, double *u)
{
  double norm = fabs(u[0]);

  for (int i = 1; i < n; i++)
  {
    double val = fabs(u[i]);
    if (val > norm) norm = val;
  }

  return(norm);
}


// -------------------------------------------------------------------------
// Matrix functions:
//

// =============================== MatAlloc ================================

double **MatAlloc(int n, int m)
{
  double **A = new double*[n];
  if (!A)
  {
    cout << "Allocation failure in MatAlloc !!!\n";
    exit(1);
  }

  for (int i = 0; i < n; i++)
  {
    A[i] = new double[m];
    if (!A[i])
    {
      cout << "Allocation failure in MatAlloc !!!\n";
      exit(1);
    }
  }

  return(A);
}

// =============================== MatFree ================================

void MatFree(double **A, int n)
{
  for (int i = 0; i < n; i++) delete [](A[i]);

  delete []A;
}

// ================================ MatZero ================================

void MatZero(int n, int m, double **A)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) A[i][j] = 0.0;
#else
  for (int i = 0; i < n; i++) VecZero(m, A[i]);
#endif
}

// =============================== MatAssign ===============================

void MatAssign(int n, int m, double **A, double **B)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) B[i][j] = A[i][j];
#else
  for (int i = 0; i < n; i++) VecAssign(m, A[i], B[i]);
#endif
}

// ================================ MatAdd =================================

void MatAdd(int n, int m, double **A, double **B)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) B[i][j] += A[i][j];
#else
  for (int i = 0; i < n; i++) VecAdd(m, A[i], B[i]);
#endif
}

// ================================ MatAdd =================================

void MatAdd(int n, int m, double **A, double **B, double **C)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) C[i][j] = A[i][j] + B[i][j];
#else
  for (int i = 0; i < n; i++) VecAdd(m, A[i], B[i], C[i]);
#endif
}

// ================================ MatSub =================================

void MatSub(int n, int m, double **A, double **B)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) B[i][j] -= A[i][j];
#else
  for (int i = 0; i < n; i++) VecSub(m, A[i], B[i]);
#endif
}

// ================================ MatSub =================================

void MatSub(int n, int m, double **A, double **B, double **C)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) C[i][j] = A[i][j] - B[i][j];
#else
  for (int i = 0; i < n; i++) VecSub(m, A[i], B[i], C[i]);
#endif
}

// ================================ MatMult ================================

void MatMult(int n, int m, double f, double **A, double **B)
{
#if 0
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) B[i][j] = f*A[i][j];
#else
  for (int i = 0; i < n; i++) VecMul(m, f, A[i], B[i]);
#endif
}

// ================================ MatMult ================================

void MatMult(int n, int m, int l, double **A, double **B, double **C)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      C[i][j] = 0.0;
      for (int k = 0; k < l; k++) C[i][j] += A[i][k]*B[k][j];
    }
  }
}

// ============================== MatVecMult ===============================

void MatVecMult(int n, int m, double **A, double *u, double *v)
{
#if 0
  for (int i = 0; i < n; i++)
  {
    v[i] = 0.0;
    for (int j = 0; j < m; j++) v[i] += A[i][j]*u[j];
  }
#else
  for (int i = 0; i < n; i++) v[i] = VecDot(m, A[i], u);
#endif
}

// ============================= MatVecMultAcc =============================

void MatVecMultAcc(int n, int m, double **A, double *u, double *v)
{
  for (int i = 0; i < n; i++) v[i] += VecDot(m, A[i], u);
}

// ============================== MatVecMultT ==============================

void MatVecMultT(int m, int n, double **A, double *u, double *v)
{
  VecZero(m, v);
  for (int j = 0; j < n; j++)
  {
    double *Aj = &A[j][0];
    double  uj = u[j];
    for (int i = 0; i < m; i++) v[i] += Aj[i]*uj;  // v[i] += A[j][i]*u[j]
  }
}

// ============================= MatVecMultTAcc ============================

void MatVecMultTAcc(double a, int m, int n, double **A, double *u,
                    double *v)
{
  for (int j = 0; j < n; j++)
  {
    double *Aj = &A[j][0];
    double auj = a*u[j];
    for (int i = 0; i < m; i++) v[i] += Aj[i]*auj; // v[i] += a*A[j][i]*u[j]
  }
}

// ============================== MatTripBtCB ==============================

void MatTripBtCB(int m, int n, double **B, double **C, double w, double **K)
{
  double aux1,aux2;

  for (int i = 0; i < n; i++)
  {
    for (int k = 0; k < m; k++) if (B[k][i] != 0.0)
    {
      aux1 = w*B[k][i];
      for (int l = 0; l < m; l++) if (C[k][l] != 0.0)
      {
        aux2 = aux1*C[k][l];
        for (int j = 0; j < n; j++) K[i][j] += aux2*B[l][j];
      }
    }
  }
}

// ============================== MatDecompLU ==============================

int MatDecompLU(int n, double **A)
{
  int i,j,k;

  // Perform the Crout factorization [A] = [L][U], where [U] has 1?s at the
  // diagonal. The factorization is carried-out in-place, i.e. the elements
  // Lij and Uij are stored in the same place occupied previously the
  // elements Aij below and above the diagonal, respectively.

  for (i = 1; i < n; i++)
  {
    for (j = 0; j < i; j++)
    {
      for (k = 0; k < j; k++)
      {
        A[i][j] -= A[i][k]*A[k][j]; // Lower matrix
        A[j][i] -= A[j][k]*A[k][i]; // Upper matrix
      }

      if (fabs(A[j][j]) < TOLC) return(0);

      A[j][i] /= A[j][j];          // Upper matrix

      A[i][i] -= A[i][j]*A[j][i];  // Diagonal
    }
  }

  return(1);
}

#if 0
int MatDecompLU(int n, double **A)
{
  int i,j,k;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j <= i; j++)  // Lower matrix
    {
      for (k = 0; k < j; k++) A[i][j] -= A[i][k]*A[k][j];
    }

    if (fabs(A[i][i]) < TOLC) return(0);

    for (j = i+1; j < n; j++) // Upper matrix
    {
      for (k = 0; k < i; k++) A[i][j] -= A[i][k]*A[k][j];
      A[i][j] /= A[i][i];
    }
  }

  return(1);
}

#endif

// ============================== MatSolveLU ===============================

void MatSolveLU(int n, double **A, double *v)
{
  int i,j;

  for (i = 0; i < n; i++)
  {
  for (j = 0; j < i; j++) v[i] -= A[i][j]*v[j];
    v[i] /= A[i][i];
  }

  for (i = n-1; i >= 0; i--)
  {
  for (j = i+1; j < n; j++) v[i] -= A[i][j]*v[j];
  }
}

// -------------------------------------------------------------------------
// Skyline matrix functions:
//

// =============================== SkylAlloc ===============================

double **SkylAlloc(int n, int *p)
{
  double **A = new double*[n];
  if (!A)
  {
    cout << "Allocation failure in SkylAlloc !!!\n";
    exit(1);
  }

  double *v;
  for (int i = 0; i < n; i++)
  {
    v = new double[i-p[i]+1];
    if (!v)
    {
      cout << "Allocation failure in SkylAlloc !!!\n";
      exit(1);
    }
    A[i] = v - p[i];
  }

  return(A);
}

// =============================== SkylFree ================================

void SkylFree(double **A, int n, int *p)
{
  for (int i = 0; i < n; i++) delete [](A[i]+p[i]);

  delete []A;
}

// =============================== SkylZero ================================

void SkylZero(int n, int *p, double **A)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) A[i][j] = 0.0;
}

// ============================== SkylAssign ===============================

void SkylAssign(int n, int *p, double **A, double **B)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) B[i][j] = A[i][j];
}

// ================================ SkylAdd ================================

void SkylAdd(int n, int *p, double **A, double **B)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) B[i][j] += A[i][j];
}

// ================================ SkylAdd ================================

void SkylAdd(double a, int n, int *p, double **A, double **B)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) B[i][j] += a*A[i][j];
}

// ================================ SkylAdd ================================

void SkylAdd(int n, int *p, double **A, double **B, double **C)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) C[i][j] = A[i][j] + B[i][j];
}

// ================================ SkylSub ================================

void SkylSub(int n, int *p, double **A, double **B)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) B[i][j] -= A[i][j];
}

// ================================ SkylSub ================================

void SkylSub(int n, int *p, double **A, double **B, double **C)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) C[i][j] = A[i][j] - B[i][j];
}

// =============================== SkylMult ================================

void SkylMult(double a, int n, int *p, double **A)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) A[i][j] *= a;
}

// =============================== SkylMult ================================

void SkylMult(double a, int n, int *p, double **A, double **B)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) B[i][j] = a*A[i][j];
}

// ============================== SkylMultAcc ==============================

void SkylMultAcc(double a, int n, int *p, double **A, double **B)
{
  for (int i = 0; i < n; i++)
    for (int j = p[i]; j <= i; j++) B[i][j] += a*A[i][j];
}

// ============================== SkylVecMult ==============================

void SkylVecMult(int n, int *p, double **A, double *u, double *v)
{
  int i,j;

  for (i = 0; i < n; i++)
  {
    v[i] = 0.0;

    for (j = p[i]; j < i; j++)
    {
      v[i] += A[i][j]*u[j];
    }

    for (j = i; j < n; j++)
    {
      if (p[j] <= i) v[i] += A[j][i]*u[j];
    }
  }
}

// ============================== CroutSolver ==============================

int CroutSolver(int op, double **A, double *b, int n, int *p)
{
  int i,j,k;

  // Reduction of matrix [A].

  if (op != 3)
  {
    for (j = 1; j < n; j++)
    {
      for (i = p[j]+1; i < j; i++)
        for (k = Max(p[i], p[j]); k < i; k++)
          A[j][i] -= A[i][k] * A[j][k];

      for (k = p[j]; k < j; k++)
      {
        if (fabs(A[k][k]) < TOLC)
        {
          cout << "\nSingular matrix in CroutSolver.\n";
          cout << "Equation = " << k << "\n";
          return(0);
        }

        A[j][j] -= A[j][k] / A[k][k] * A[j][k];
        A[j][k] /= A[k][k];
      }
    }
  }

  if (op == 2) return(1);

  // Reduction of vector {b}.

  for (i = 1; i < n; i++)
  {
    for (k = p[i]; k < i; k++)
      b[i] -= A[i][k] * b[k];
  }

  for (i = 0; i < n; i++)
  {
    b[i] /= A[i][i];
  }

  // Back-substitution.

  for (i = n-1; i > 0; i--)
  {
    for (k = p[i]; k < i; k++)
      b[k] -= A[i][k] * b[i];
  }

  return(1);
}

// ======================================================= End of file =====
