// -------------------------------------------------------------------------
// vec.h - definition of the vector class.
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
// void Zero(void)
//
// Initializes the vector.
// -------------------------------------------------------------------------
//
// double Length(void)
//
// Evaluates and returns the vector length.
// -------------------------------------------------------------------------
//
// double NormInf(void)
//
// Evaluates and returns the vector infinite norm (max. absolute value).
// -------------------------------------------------------------------------
//
// void Normalize(void)
//
// Normalize the vector, turning it in an unit vector.
// -------------------------------------------------------------------------
// void Sort(eVecSort order)
//
//   order  -  ASCENDING/DESCENDING                                    (in)
//
// Sort the vector in the chosen order.
// -------------------------------------------------------------------------
//
// double Max(void)
//
// Returns the maximum value in the vector.
// -------------------------------------------------------------------------
//
// double Min(void)
//
// Returns the minimum value in the vector.
// -------------------------------------------------------------------------
//
// void Resize(int n)
//
//   n  -  new dimension                                               (in)
//
// Redefine the dimension of the stored array.
// -------------------------------------------------------------------------
//
// void Print(const char *fmt)
//
//   fmt  -  given format                                              (in)
//
// Prints the stored data in the given format.
// -------------------------------------------------------------------------
//
// double& operator[](int i)
//
//   i  -  index of desired element                                    (in)
//
// Returns a reference to element of index "i".
// -------------------------------------------------------------------------
//
// double& operator()(int i)
//
//   i  -  index of desired element                                    (in)
//
// Returns a reference to element of index "i".
//
// -------------------------------------------------------------------------
// Public operators:
// -------------------------------------------------------------------------
//
// (1) Implemented in the conventional form:
//
// {v}  = {u}
// {v} += {u}
// {v} -= {u}
// {v} *= a
// {v} /= a
//
// (2) Implemented using compositors to improve the computer efficiency:
//
// {w}  = a*{u}
// {w} += a*{u}
// {w}  = {u} + {v}
// {w}  = {u} - {v}
// {w}  = {u}/{v}           (w[i] = u[i]/v[i])
// {w}  = a*{u} + {v}
// {w} += a*{u} + {v}
// {w}  = a*{u} - {v}
// {w} += a*{u} - {v}
// {w}  = [A]*{u}
// {w} += [A]*{u}
// {w}  = [A]t*{u}
// {w}  = [S]*{u}           (Sparse matrix)
//
// Important: the vectors and matrices must have the correct dimensions
//            required by the chosen operation, otherwise errors will
//            occurs.
//
// -------------------------------------------------------------------------
// Public friend functions:
// -------------------------------------------------------------------------
//
// friend void CrossProd(cVector &u, cVector &v, cVector &w)
//
//   u  -  given vector                                                (in)
//   v  -  given vector                                                (in)
//   w  -  result vector                                              (out)
//
// Returns the cross product {w} = {u}*{v}.
// -------------------------------------------------------------------------
//
// friend double operator*(const cVector &u, const cVector &v)
//
//   u  -  given vector                                                (in)
//   v  -  given vector                                                (in)
//
// Returns the scalar product (dot) of the vectors {u}*{v}.
// -------------------------------------------------------------------------
//
// void MultTAcc(double a, const cMatrix &B, const cVector &s, cVector &g)
//
//   a - given scalar                                                  (in)
//   B - given matrix                                                  (in)
//   s - given vector                                                  (in)
//   g - output vector                                                (out)
//
// Performs the operation: {g} += a*[B]t*{s}.
// -------------------------------------------------------------------------

#ifndef _VEC_H
#define _VEC_H

#include <assert.h>

#include "mat.h"
#include "sysmat.h"
#include "matvec.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Auxiliary types:
//
typedef enum
{
  ASCENDING,
  DESCENDING
} eVecSort;

// -------------------------------------------------------------------------
// Definition of compositors:
//
struct sSclVec  // Operation {w} = a*{u}
{
  const double  &scl;
  const cVector &vec;
  sSclVec(const double &a, const cVector &u) : scl(a), vec(u) { }
};

inline sSclVec operator*(const double &a, const cVector &v)
{
  return(sSclVec(a, v));
}

struct sAddVec  // Operation {w} = {u} + {v}
{
  const cVector &lvec;
  const cVector &rvec;
  sAddVec(const cVector &u, const cVector &v) : lvec(u), rvec(v) { }
};

inline sAddVec operator+(const cVector &u, const cVector &v)
{
  return(sAddVec(u, v));
}

struct sSubVec  // Operation {w} = {u} - {v}
{
  const cVector &lvec;
  const cVector &rvec;
  sSubVec(const cVector &u, const cVector &v) : lvec(u), rvec(v) { }
};

inline sSubVec operator-(const cVector &u, const cVector &v)
{
  return(sSubVec(u, v));
}

struct sDivVec  // Operation {w} = {u}/{v}
{
  const cVector &lvec;
  const cVector &rvec;
  sDivVec(const cVector &u, const cVector &v) : lvec(u), rvec(v) { }
};

inline sDivVec operator/(const cVector &u, const cVector &v)
{
  return(sDivVec(u, v));
}

struct sSclAddVec  // Operation {w} = a*{u} + {v}
{
  const sSclVec &prod;
  const cVector &rvec;
  sSclAddVec(const sSclVec &au, const cVector &v) : prod(au), rvec(v) { }
};

inline sSclAddVec operator+(const sSclVec &au, const cVector &v)
{
  return(sSclAddVec(au, v));
}

inline sSclAddVec operator+(const cVector &v, const sSclVec &au)
{
  return(sSclAddVec(au, v));
}

struct sSclSubVec  // Operation {w} = a*{u} - {v}
{
  const sSclVec &prod;
  const cVector &rvec;
  sSclSubVec(const sSclVec &au, const cVector &v) : prod(au), rvec(v) { }
};

inline sSclSubVec operator-(const sSclVec &au, const cVector &v)
{
  return(sSclSubVec(au, v));
}

struct sMatVec  // Operation {w} = [A]*{u}
{
  const cMatrix &mat;
  const cVector &vec;
  sMatVec(const cMatrix &A, const cVector &u) : mat(A), vec(u) { }
};

inline sMatVec operator*(const cMatrix &A, const cVector &u)
{
  return(sMatVec(A, u));
}

struct sTrnMatVec  // Operation {w} = [A]t*{u}
{
  const cMatrix &mat;
  const cVector &vec;
  sTrnMatVec(const cMatrix &A, const cVector &u) : mat(A), vec(u) { }
};

inline sTrnMatVec operator*(const sTrnMat &At, const cVector &u)
{
  return(sTrnMatVec(At.mat, u));
}

struct sSysMatVec  // Operation {w} = [S]*{u}
{
  cSysMatrix *mat;
  const cVector &vec;
  sSysMatVec(cSysMatrix *S, const cVector &u) : mat(S), vec(u) { }
};

inline sSysMatVec operator*(cSysMatrix *S, const cVector &u)
{
  return(sSysMatVec(S, u));
}

// -------------------------------------------------------------------------
// Definition of cVector class:
//
class cVector
{
 private:
  int     dim;   // Vector size
  double *val;   // Vector values

 public:
                  cVector   (void);
                  cVector   (int);
                  cVector   (int, double *);
                  cVector   (const cVector &);
                 ~cVector   (void);
         int      Dim       (void) { return dim; }  // temporario !
         double*  Val       (void) { return val; }  // temporario !
         void     Zero      (void) { VecZero(dim, val); }
         double   Length    (void) { return VecLen(dim, val); }
         double   NormInf   (void) { return VecNormInf(dim, val); }
         void     Normalize (void);
	 void     Sort      (eVecSort);
         double   Max       (void); 
         double   Min       (void); 
	 void     Resize    (int);
         void     Print     (void);
         double&  operator[](int i) { return val[i]; }
   const double&  operator[](int i) const { return val[i]; }
         double&  operator()(int i) { return val[i]; }
         cVector& operator= (const cVector &);
         cVector& operator= (const double &);
         cVector& operator= (const sSclVec &);
         cVector& operator+=(const sSclVec &);
         cVector& operator= (const sAddVec &);
         cVector& operator= (const sSubVec &);
         cVector& operator= (const sDivVec &);
         cVector& operator= (const sSclAddVec &);
         cVector& operator+=(const sSclAddVec &);
         cVector& operator= (const sSclSubVec &);
         cVector& operator+=(const sSclSubVec &);
         cVector& operator= (const sMatVec &);
         cVector& operator+=(const sMatVec &);
         cVector& operator= (const sTrnMatVec &);
         cVector& operator= (const sSysMatVec &);
         void     operator+=(const cVector &);
         void     operator-=(const cVector &);
         void     operator*=(const double &);
         void     operator/=(const double &);

  friend void     CrossProd (cVector &, cVector &, cVector &);
  friend double   operator* (const cVector &, const cVector &);
  friend void     MultTAcc  (double, const cMatrix &, const cVector &,
                             cVector &);
};

// -------------------------------------------------------------------------
// Inline functions:
//

// ============================= operator= =================================

#if 1
inline cVector& cVector :: operator=(const cVector &vec)
{
  assert(dim == vec.dim);
  if (vec.val != val)              // Avoid self copy
    VecAssign(dim, vec.val, val);
  return(*this);
}
#endif

// ============================= operator= =================================

inline cVector& cVector :: operator=(const double &scl)
{
  VecSet(dim, scl, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sSclVec &op)
{
  assert(dim == op.vec.dim);
  VecMul(dim, op.scl, op.vec.val, val);
  return(*this);
}

// ============================= operator+= ================================

inline cVector& cVector :: operator+=(const sSclVec &op)
{
  assert(dim == op.vec.dim);
  VecMulAdd(dim, op.scl, op.vec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sAddVec &op)
{
  assert(dim == op.lvec.dim && dim == op.rvec.dim);
  VecAdd(dim, op.lvec.val, op.rvec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sSubVec &op)
{
  assert(dim == op.lvec.dim && dim == op.rvec.dim);
  VecSub(dim, op.lvec.val, op.rvec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sDivVec &op)
{
  assert(dim == op.lvec.dim && dim == op.rvec.dim);
  VecDiv(dim, op.lvec.val, op.rvec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sSclAddVec &op)
{
  assert(dim == op.prod.vec.dim && dim == op.rvec.dim);
  VecSclAdd(dim, op.prod.scl, op.prod.vec.val, op.rvec.val, val);
  return(*this);
}

// ============================ operator+= =================================

inline cVector& cVector :: operator+=(const sSclAddVec &op)
{
  assert(dim == op.prod.vec.dim && dim == op.rvec.dim);
  VecSclAddInc(dim, op.prod.scl, op.prod.vec.val, op.rvec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sSclSubVec &op)
{
  assert(op.prod.vec.dim == op.rvec.dim);
  VecSclSub(dim, op.prod.scl, op.prod.vec.val, op.rvec.val, val);
  return(*this);
}

// ============================ operator+= =================================

inline cVector& cVector :: operator+=(const sSclSubVec &op)
{
  assert(op.prod.vec.dim == op.rvec.dim);
  VecSclSubInc(dim, op.prod.scl, op.prod.vec.val, op.rvec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sMatVec &op)
{
  assert(dim == op.mat.nrow && op.mat.ncol == op.vec.dim);
  MatVecMult(op.mat.nrow, op.mat.ncol, op.mat.val, op.vec.val, val);
  return(*this);
}

// ============================= operator+= ================================

inline cVector& cVector :: operator+=(const sMatVec &op)
{
  assert(dim == op.mat.nrow && op.mat.ncol == op.vec.dim);
  MatVecMultAcc(op.mat.nrow, op.mat.ncol, op.mat.val, op.vec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sTrnMatVec &op)
{
  assert(dim == op.mat.ncol && op.mat.nrow == op.vec.dim);
  MatVecMultT(op.mat.ncol, op.mat.nrow, op.mat.val, op.vec.val, val);
  return(*this);
}

// ============================= operator= =================================

inline cVector& cVector :: operator=(const sSysMatVec &op)
{
  assert(dim == op.mat->dim && op.mat->dim == op.vec.dim);
  op.mat->MultVect(op.vec.val, val);
  return(*this);
}

// ============================= operator+= ================================

inline void cVector :: operator+=(const cVector &u)
{
  assert(dim == u.dim);
  VecAdd(dim, u.val, val);
}

// ============================= operator-= ================================

inline void cVector :: operator-=(const cVector &u)
{
  assert(dim == u.dim);
  VecSub(dim, u.val, val);
}

// ============================= operator* =================================

inline double operator*(const cVector &u, const cVector &v)
{
  assert(u.dim == v.dim);
  return(VecDot(u.dim, u.val, v.val));
}

// ============================== MultTAcc =================================

inline void MultTAcc(double a, const cMatrix &B, const cVector &s,
                     cVector &g)
{
  assert(B.nrow == s.dim);
  MatVecMultTAcc(a, B.ncol, B.nrow, B.val, s.val, g.val);
}

#endif
