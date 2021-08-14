// -------------------------------------------------------------------------
// mat.h - definition of the matrix class.
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
// Public methods:
// -------------------------------------------------------------------------
//
// void Print(void)
//
// Prints the matrix on the screen.
// -------------------------------------------------------------------------
//
// void Zero(void)
//
// Initializes the matrix.
// -------------------------------------------------------------------------
//
// void Resize(int n, m)
//
//   n,m  -  new dimensions                                           (in)
//
// Redefine the dimensions of the stored matrix.
// ------------------------------------------------------------------------
//
// double& operator()(int i, int j)
//
//   i  -  row of the matrix                                           (in)
//   j  -  column of the matrix                                        (in)
//
// Returns the "i,j" element of the matrix.
//
// -------------------------------------------------------------------------
// Public operators:
// -------------------------------------------------------------------------
//
// (1) Implemented in the conventional form:
//
// [B]  = [A]
// [B] += [A]
// [B] -= [A]
//
// (2) Implemented using compositors to improve the computer efficiency:
//
// [C] = a*[A]
// [C] = [A] + [B]
// [C] = [A] - [B]
// [C] = [A]*[B]
//
// Important: the matrices must have the correct dimensions required by
//            the chosen operation, otherwise errors will occurs.
//
// -------------------------------------------------------------------------
// Public friend functions:
// -------------------------------------------------------------------------
//
// friend void MatTripBtCB(const cMatrix &B, const cMatrix &C, double wgt,
//                         cMatrix &K)
//
//   B    -  given (strain-displacement) matrix                        (in)
//   C    -  given (constitutive) matrix                               (in)
//   wgt  -  given coefficient                                         (in)
//   K    -  updated (stiffness) matrix                               (out)
//
// Performs the operation [K] += wgt*[B]t*[C]*[B].
// -------------------------------------------------------------------------

#ifndef _MAT_H
#define _MAT_H

#include <assert.h>

#include "matvec.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;

// -------------------------------------------------------------------------
// Definition of compositors to improve the efficiency of some operations:
//
struct sSclMat  // Operation [C] = a*[A]
{
  const double  &scl;
  const cMatrix &mat;
  sSclMat(const double &a, const cMatrix &A) : scl(a), mat(A) { }
};

inline sSclMat operator*(const double &a, const cMatrix &A)
{
  return(sSclMat(a, A));
}

struct sAddMat  // Operation [C] = [A] + [B]
{
  const cMatrix &lmat;
  const cMatrix &rmat;
  sAddMat(const cMatrix &A, const cMatrix &B) : lmat(A), rmat(B) { }
};

inline sAddMat operator+(const cMatrix &A, const cMatrix &B)
{
  return(sAddMat(A, B));
}

struct sSubMat  // Operation [C] = [A] - [B]
{
  const cMatrix &lmat;
  const cMatrix &rmat;
  sSubMat(const cMatrix &A, const cMatrix &B) : lmat(A), rmat(B) { }
};

inline sSubMat operator-(const cMatrix &A, const cMatrix &B)
{
  return(sSubMat(A, B));
}

struct sMulMat  // Operation [C] = [A]*[B]
{
  const cMatrix &lmat;
  const cMatrix &rmat;
  sMulMat(const cMatrix &A, const cMatrix &B) : lmat(A), rmat(B) { }
};

inline sMulMat operator*(const cMatrix &A, const cMatrix &B)
{
  return(sMulMat(A, B));
}

struct sTrnMat  // Stores [A] for operations involving [A]t (transpose)
{
  const cMatrix &mat;
  sTrnMat(const cMatrix &A) : mat(A) { }
};

inline sTrnMat t(const cMatrix &A)
{
  return(sTrnMat(A));
}

// -------------------------------------------------------------------------
// Definition of cMatrix class:
//
class cMatrix
{
 friend class cVector;  // Allow cVector access the private data

 private:
  int      nrow;   // Number of rows
  int      ncol;   // Number of colums
  double **val;    // Matrix values

  void AllocMem  (void);
  void ReleaseMem(void);

 public:
                  cMatrix    (void);
                  cMatrix    (int, int);
                  cMatrix    (int, int, double **);
                  cMatrix    (const cMatrix &);
                  cMatrix    (const sSclMat &);
                  cMatrix    (const sAddMat &);
                  cMatrix    (const sSubMat &);
                  cMatrix    (const sMulMat &);
                 ~cMatrix    (void);
         int      NRow       (void) { return nrow; }
         int      NCol       (void) { return ncol; }
         double **Val        (void) { return val; }
         void     Resize     (int, int);
         void     Print      (void);
         void     Zero       (void) { MatZero(nrow, ncol, val); }
         void     Transp     (cMatrix &);
         int      CompInverse(cMatrix &);
         int      Solve      (cVector &, cVector &);
         int      DecompLU   (void);
         void     SolveLU    (cVector &);
         double*  operator[] (int i) { return val[i]; }
         double&  operator() (int i, int j) { return val[i][j]; }
         cMatrix& operator=  (const cMatrix &);
         cMatrix& operator=  (const sSclMat &);
         cMatrix& operator=  (const sAddMat &);
         cMatrix& operator=  (const sSubMat &);
         cMatrix& operator=  (const sMulMat &);
         void     operator+= (const cMatrix &);
         void     operator-= (const cMatrix &);

  friend void     MultTAcc   (double, const cMatrix &, const cVector &,
                              cVector &);
  friend void     MatTripBtCB(const cMatrix &, const cMatrix &, double,
                              cMatrix &);
};

// -------------------------------------------------------------------------
// Inline functions:
//

// ============================= operator= =================================

#if 1
inline cMatrix& cMatrix :: operator=(const cMatrix &mat)
{
  assert(nrow == mat.nrow && ncol == mat.ncol);
  MatAssign(nrow, ncol, mat.val, val);
  return(*this);
}
#endif

// ============================= operator= =================================

inline cMatrix& cMatrix :: operator=(const sSclMat &op)
{
  assert(nrow == op.mat.nrow && ncol == op.mat.ncol);
  MatMult(nrow, ncol, op.scl, op.mat.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cMatrix& cMatrix :: operator=(const sAddMat &op)
{
  assert(nrow == op.lmat.nrow && ncol == op.lmat.ncol);
  assert(nrow == op.rmat.nrow && ncol == op.rmat.ncol);
  MatAdd(nrow, ncol, op.lmat.val, op.rmat.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cMatrix& cMatrix :: operator=(const sSubMat &op)
{
  assert(nrow == op.lmat.nrow && ncol == op.lmat.ncol);
  assert(nrow == op.rmat.nrow && ncol == op.rmat.ncol);
  MatSub(nrow, ncol, op.lmat.val, op.rmat.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cMatrix& cMatrix :: operator=(const sMulMat &op)
{
  assert(op.lmat.ncol == op.rmat.nrow);
  assert(nrow == op.lmat.nrow && ncol == op.rmat.ncol);
  MatMult(nrow, ncol, op.lmat.ncol, op.lmat.val, op.rmat.val, val);
  return(*this);
}

// ============================= operator+= ================================

inline void cMatrix :: operator+=(const cMatrix &mat)
{
  assert(nrow == mat.nrow && ncol == mat.ncol);
  MatAdd(nrow, ncol, mat.val, val);
}

// ============================= operator-= ================================

inline void cMatrix :: operator-=(const cMatrix &mat)
{
  assert(nrow == mat.nrow && ncol == mat.ncol);
  MatSub(nrow, ncol, mat.val, val);
}

// ============================== MatTripBtCB ==============================

inline void MatTripBtCB(const cMatrix &B, const cMatrix &C, double wgt,
                        cMatrix &K)
{
  assert(B.nrow == C.ncol && C.nrow == C.ncol);
  MatTripBtCB(B.nrow, B.ncol, B.val, C.val, wgt, K.val);
}

#endif
