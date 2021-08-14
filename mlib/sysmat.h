// -------------------------------------------------------------------------
// sysmat.h - Definition of classes to handle square sparse matrices of
//            large linear systems.
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
// Remarks:
// -------------------------------------------------------------------------
// cSysMatrix => class for large scale square matrices
//   cSymSkylMatrix => symmetric skyline matrices (lower triangular)
//   cUnsymSkylMatrix => general skyline matrices
//   cSprsMatrix => general sparse matrix
//     cSymSprsMatrix => symmetric sparse matrix (lower triangular)
//
// -------------------------------------------------------------------------

#ifndef _SYSMAT_H
#define _SYSMAT_H

// -------------------------------------------------------------------------
// SysMat types:
//
typedef enum
{
  SYM_SKYLINE,
  UNSYM_SKYLINE,
  SYM_SPARSE,
  UNSYM_SPARSE
} eSMType;

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Abstract base class to handle sparse matrices:
//
class cSysMatrix
{
 friend class cVector;  // Allow cVector access the private data

 protected:
  int dim;        // Matrix dimension (rows/columns)
  int nelm;       // Number of stored elements
  int maxit;      // Max iterations for linear solver
  double tol;     // Solver tolerance

 protected:
  int IsSup(int i, int j) { return(i < j); }
  int IsInf(int i, int j) { return(i > j); }

 public:
  static cSysMatrix *CreateMatrix(eSMType, int, int *);

                   cSysMatrix(int n = 1) { dim = n; nelm = 0;
                                           maxit = n; tol = 1.0e-05; }
  virtual         ~cSysMatrix(void) { }
          int      GetNumElm (void) { return nelm; }
          int      Dim       (void) { return dim; }
          void     SetParam  (int k, double t) { maxit = k; tol = t; }
  virtual int      Symmetric (void) { return 0; }
  virtual double **Val       (void) { return 0; }
  virtual void     Zero      (void) { }
  virtual int      Solve     (cVector &, cVector &) { return 0; }
  virtual int      GetSmlPvt (void) { return 0; }
  virtual int      GetLgtPvt (void) { return 0; }
  virtual int      GetNgtPvt (void) { return 0; }
  virtual double   Get       (int, int) = 0;
  virtual void     Add       (int, int, double) = 0;
  virtual void     Print     (void) = 0;
  virtual void     MultVect  (double *, double *) = 0;
  virtual void     AddMat    (double, cSysMatrix *) { }
};

// -------------------------------------------------------------------------
// Symmetric skyline matrices:
//
class cSymSkylMatrix : public cSysMatrix
{
 protected:
  int    dec;     // Decomposition flag
  int    *skl;    // Matrix skyline
  double **val;   // Stored values

 public:
           cSymSkylMatrix(int, int *);
          ~cSymSkylMatrix(void);
  double&  operator()    (int i, int j) { return val[i][j]; }
  double **Val           (void) { return val; }
  int      Symmetric     (void) { return 1; }
  void     Zero          (void);
  int      Solve         (cVector &, cVector &);
  int      GetSmlPvt     (void);
  int      GetLgtPvt     (void);
  int      GetNgtPvt     (void);
  double   Get           (int, int);
  void     Add           (int, int, double);
  void     Print         (void);
  void     MultVect      (double *, double *);
  void     AddMat        (double, cSysMatrix *);
};

// -------------------------------------------------------------------------
// Unsymmetric matrices with symmetric skyline:
//
class cUnsymSkylMatrix : public cSysMatrix
{
 protected:
  int    dec;     // Decomposition flag
  int    *skl;    // Matrix skyline
  double **L;     // Lower matrix
  double **U;     // Upper matrix

 public:
           cUnsymSkylMatrix(int, int *);
          ~cUnsymSkylMatrix(void);
  void     Zero          (void);
  int      Solve         (cVector &, cVector &);
  int      GetSmlPvt     (void);
  int      GetLgtPvt     (void);
  int      GetNgtPvt     (void);
  double   Get           (int, int);
  void     Add           (int, int, double);
  void     Print         (void);
  void     MultVect      (double *, double *);
  void     AddMat        (double, cSysMatrix *);
};

// -------------------------------------------------------------------------
// Sparse matrix element:
//
typedef struct _spelm SpElm;

struct _spelm
{
  int    col;
  double val;
  SpElm *nxt;
};

// -------------------------------------------------------------------------
// Class to handle nonsymentric sparse matrices:
//
class cSprsMatrix : public cSysMatrix
{
 protected:
  SpElm **row;   // Array of rows (list of nonzero elements)

 public:
                  cSprsMatrix(int n = 1);
  virtual        ~cSprsMatrix(void);
  virtual double  Get        (int, int);
  virtual void    Add        (int, int, double);
  virtual void    Print      (void);
  virtual void    MultVect   (double *, double *);
};

// -------------------------------------------------------------------------
// Class to handle symmetric sparse matrices:
//
class cSymSprsMatrix : public cSprsMatrix
{
 public:
          cSymSprsMatrix(int n = 1);
         ~cSymSprsMatrix(void) { }
  int     Symmetric     (void) { return 1; }
  int     Solve         (cVector &, cVector &);
  double  Get           (int, int);
  void    Add           (int, int, double);
  void    Print         (void);
  void    MultVect      (double *, double *);
};

// -------------------------------------------------------------------------
// Generic Solvers:
//
int PCGSolver(int maxit, double tol, cSysMatrix *K, cVector &f,
              cVector &u);

int FstEigPair(int maxit, double tol, cSysMatrix *K, double *, cVector &v);

#endif
