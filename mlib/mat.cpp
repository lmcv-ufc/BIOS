// -------------------------------------------------------------------------
// mat.cpp - implementation of the matrix class.
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
// Created:  14-Jan-1998    Evandro Parente Junior
//
// Modified: 02-Aug-2000    Evandro Parente Junior
//           Use of new/delete instead of calloc/free.
//
// Modified: 21-Sep-2001    Evandro Parente Junior
//           Use of compositors to increase the effic. of some operators.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

using namespace std;

#include "mat.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Private methods:
//

// ================================ AllocMem ===============================

void cMatrix :: AllocMem(void)
{
  val = new double*[nrow];

  for (int i = 0; i < nrow; i++) val[i] = new double[ncol];
}

// =============================== ReleaseMem ==============================

void cMatrix :: ReleaseMem(void)
{
  for (int i = 0; i < nrow; i++) delete []val[i];

  delete []val;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cMatrix ==================================

cMatrix :: cMatrix(void)
{
  nrow = 0;
  ncol = 0;
  val  = 0;
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(int n, int m)
{
  assert(n > 0 && m > 0);
  nrow = n;
  ncol = m;
  AllocMem( );
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(int n, int m, double **mat)
{
  assert(n > 0 && m > 0);
  nrow = n;
  ncol = m;
  AllocMem( );
  MatAssign(nrow, ncol, mat, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const cMatrix &mat)
{
  nrow = mat.nrow;
  ncol = mat.ncol;
  AllocMem( );
  MatAssign(nrow, ncol, mat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sSclMat &op)
{
  nrow = op.mat.nrow;
  ncol = op.mat.ncol;
  AllocMem( );
  MatMult(nrow, ncol, op.scl, op.mat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sAddMat &op)
{
  assert(op.lmat.nrow == op.rmat.nrow && op.lmat.ncol == op.rmat.ncol);
  nrow = op.lmat.nrow;
  ncol = op.lmat.ncol;
  AllocMem( );
  MatAdd(nrow, ncol, op.lmat.val, op.rmat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sSubMat &op)
{
  assert(op.lmat.nrow == op.rmat.nrow && op.lmat.ncol == op.rmat.ncol);
  nrow = op.lmat.nrow;
  ncol = op.lmat.ncol;
  AllocMem( );
  MatSub(nrow, ncol, op.lmat.val, op.rmat.val, val);
}

// ============================== cMatrix ==================================

cMatrix :: cMatrix(const sMulMat &op)
{
  assert(op.lmat.ncol == op.rmat.nrow);
  nrow = op.lmat.nrow;
  ncol = op.rmat.ncol;
  AllocMem( );
  MatMult(nrow, ncol, op.lmat.ncol, op.lmat.val, op.rmat.val, val);
}

// ============================== ~cMatrix =================================

cMatrix :: ~cMatrix(void)
{
  ReleaseMem( );
}

// =============================== Resize ==================================

void cMatrix :: Resize(int n, int m)
{
  ReleaseMem( );
  nrow = n;
  ncol = m;
  AllocMem( );
}

// =============================== Print ===================================

void cMatrix :: Print(void)
{
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++) cout << val[i][j] << " ";
    cout << "\n";
  }
}

// =============================== Transp ==================================

void cMatrix :: Transp(cMatrix &T)
{
  assert(nrow == T.ncol && ncol == T.nrow);

  for (int i = 0; i < nrow; i++)
    for (int j = 0; j < ncol; j++) T.val[j][i] = val[i][j];
}

// ============================== CompInverse ==============================

int cMatrix :: CompInverse(cMatrix &B)
{
  assert(nrow == ncol && B.nrow == B.ncol && nrow == B.nrow);
  cMatrix A(*this);  // Create a temporary copy

  int dec = MatDecompLU(nrow, A.val);
  if (!dec) return(0);

  // Compute the inverse matrix solving [A][B] = [I].

  double *e = new double[nrow];
  for (int j = 0; j < ncol; j++)
  {
    VecZero(nrow, e);
    e[j] = 1.0;
    MatSolveLU(nrow, A.val, e);
    if (!dec) break;
    for (int i = 0; i < ncol; i++) B.val[i][j] = e[i];
  }
  delete []e;
  return(1);
}

// ================================= Solve =================================

int cMatrix :: Solve(cVector &b, cVector &x)
{
  int dec = MatDecompLU(nrow, val);
  if (!dec) return(0);

  x = b;
  MatSolveLU(nrow, val, x.Val( ));
  return(1);
}

// ================================= Solve =================================

int cMatrix :: DecompLU(void)
{ 
  return(MatDecompLU(nrow, val));
}

// ================================= Solve =================================

void cMatrix :: SolveLU(cVector &x) 
{ 
  MatSolveLU(nrow, val, x.Val( )); 
}

// ================================================ End of file ============
