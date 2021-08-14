// -------------------------------------------------------------------------
// sysmat.cpp - Implementation of class to handle large sparse matrices.
// -------------------------------------------------------------------------
// Copyright (c) 2015 LMCV/UFC
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

#include <math.h>
#include <iostream>
#include <iomanip>
#include <string.h>

using namespace std;

#include "sysmat.h"
#include "vec.h"

const double TOLZ = 1.0e-12;

inline int MAX(int a, int b)
{
  if (a > b) return(a);
  return(b);
}

// -------------------------------------------------------------------------
// Static methods:
//

// ============================= CreateMatrix ==============================

cSysMatrix* cSysMatrix :: CreateMatrix(eSMType type, int n, int *p)
{
  cSysMatrix *K = 0;

  switch (type)
  {
    case SYM_SKYLINE:
      K = new cSymSkylMatrix(n, p);
    break;

    case UNSYM_SKYLINE:
      K = new cUnsymSkylMatrix(n, p);
    break;

    case SYM_SPARSE:
      K = new cSymSprsMatrix(n);
    break;

    case UNSYM_SPARSE:
      K = new cSprsMatrix(n);
    break;
  }

  return(K);
}

// -------------------------------------------------------------------------
// cSymSkylMatrix class:
// -------------------------------------------------------------------------

// ============================ cSymSkylMatrix =============================

cSymSkylMatrix :: cSymSkylMatrix(int n, int *p) : cSysMatrix(n)
{
  int    r;
  double *v;

  dec  = 0;
  dim  = n;
  skl  = new int[dim];
  val  = new double*[dim];
  nelm = 0;
  for (int i = 0; i < dim; i++)
  {
    r = i - p[i] + 1;     // Evaluate the row's length (or column height)
    v = new double[r];    // Alloc the elements inside the profile
    nelm += r;            // Increment the number of allocated elements
    val[i] = v - p[i];    // Move the allocated elements
    skl[i] = p[i];        // Store the given skyline
  }
}

// ============================ ~cSymSkylMatrix ============================

cSymSkylMatrix :: ~cSymSkylMatrix(void)
{
  for (int i = 0; i < dim; i++) delete [](val[i] + skl[i]);
  delete []val;

  delete []skl;
}

// ================================== Zero =================================

void cSymSkylMatrix :: Zero(void)
{
  dec = 0;
  SkylZero(dim, skl, val);
}

// ================================== Solve ================================

int cSymSkylMatrix :: Solve(cVector &b, cVector &x)
{
  x = b;
  if (!dec)
  {
    dec = 1;
    return(CroutSolver(1, val, x.Val( ), dim, skl));
  }
  else
  {
    return(CroutSolver(3, val, x.Val( ), dim, skl));
  }
}

// =============================== GetSmlPvt ===============================

int cSymSkylMatrix :: GetSmlPvt(void)
{
  int    idx  = 0;
  double mink = fabs(val[0][0]);

  for (int i = 1; i < dim; i++)
  {
    double curk = fabs(val[i][i]);
    if (curk < mink)
    {
      idx  = i;
      mink = curk;
    }
  }

  return(idx);
}

// =============================== GetLgtPvt ===============================

int cSymSkylMatrix :: GetLgtPvt(void)
{
  int    idx  = 0;
  double maxk = fabs(val[0][0]);

  for (int i = 1; i < dim; i++)
  {
    double curk = fabs(val[i][i]);
    if (curk > maxk)
    {
      idx  = i;
      maxk = curk;
    }
  }

  return(idx);
}

// =============================== GetNgtPvt ===============================

int cSymSkylMatrix :: GetNgtPvt(void)
{
 int nnp = 0;

 for (int i = 0; i < dim; i++)
 {
   if (val[i][i] < 0.0) nnp++;
 }

 return(nnp);
}

// ================================== Get ==================================

double cSymSkylMatrix :: Get(int i, int j)
{
  if (IsSup(i,j)) return(val[j][i]);

  return(val[i][j]);
}

// ================================== Add ==================================

void cSymSkylMatrix :: Add(int i, int j, double v)
{
  if (IsSup(i,j)) return;

  val[i][j] += v;
}

// ================================= Print =================================

void cSymSkylMatrix :: Print(void)
{
  int i,j;

  cout << "Stored elements = " << nelm << "\n";
  cout << showpos << scientific << setprecision(2);
  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < skl[i]; j++) cout << 0.0 << "  "; 
    for (j = skl[i]; j <= i; j++) cout << val[i][j] << "  "; 
    cout << "\n";
  }
  cout << noshowpos << fixed;
}

// =============================== MultVect ================================

void cSymSkylMatrix :: MultVect(double *u, double *v)
{
  int i,j;

  // Initialize the given vector

  for (i = 0; i < dim; i++) v[i] = 0.0;

  // Compute the matrix/vector product

  for (i = 0; i < dim; i++)
  {
    for (j = skl[i]; j < i; j++)
    {
      v[i] += val[i][j]*u[j];  // Inferior elements
      v[j] += val[i][j]*u[i];  // Superior elements
    }
    v[i] += val[i][i]*u[i];    // Diagonal element
  }
}

// ================================ AddMat =================================

void cSymSkylMatrix :: AddMat(double a, cSysMatrix *M)
{
  for (int i = 0; i < dim; i++)
  {
    for (int j = skl[i]; j <= i; j++) val[i][j] += a*M->Get(i, j);
  }
}

// -------------------------------------------------------------------------
// cUnsymSkylMatrix class:
// -------------------------------------------------------------------------
// This class deals with large unsymmetric matrices with symmetric
// skyline (profile).
//
// The elements of below the diagonal (including the diagonal itself) are
// stored in L matrix and the elements above the diagonal are stored in Ut.
//
// Therefore:
// Lower elements (Ai,j | i >= j) are stored at Li,j => [L]
// Upper elements (Ai,j | i <  j) are stored at Uj,i => [U]t
// -------------------------------------------------------------------------

// =========================== cUnsymSkylMatrix ============================

cUnsymSkylMatrix :: cUnsymSkylMatrix(int n, int *p) : cSysMatrix(n)
{
  int    r;
  double *v;

  dec  = 0;
  dim  = n;
  skl  = new int[dim];
  L    = new double*[dim];
  U    = new double*[dim];
  nelm = 0;
  for (int i = 0; i < dim; i++)
  {
    r = i - p[i] + 1;     // Evaluate the row's length
    v = new double[r];    // Alloc the row's elements
    L[i] = v - p[i];      // Move the allocated elements
    nelm += r;            // Increment the number of allocated elements

    r = i - p[i];         // Evaluate the column´s length
    v = new double[r];    // Alloc the column's elements
    v = new double[r];    // Alloc the elements inside the profile
    U[i] = v - p[i];      // Move the allocated elements
    nelm += r;            // Increment the number of allocated elements

    skl[i] = p[i];        // Store the given skyline
  }
}

// =========================== ~cUnsymSkylMatrix ===========================

cUnsymSkylMatrix :: ~cUnsymSkylMatrix(void)
{
  for (int i = 0; i < dim; i++)
  {
    delete [](L[i] + skl[i]);
    delete [](U[i] + skl[i]);
  }
  delete []L;
  delete []U;
  delete []skl;
}

// ================================== Zero =================================

void cUnsymSkylMatrix :: Zero(void)
{
  dec = 0;
  for (int i = 0; i < dim; i++)
  {
    for (int j = skl[i]; j < i; j++)
    {
      L[i][j] = 0.0;
      U[i][j] = 0.0;
    }
    L[i][i] = 0.0;
  }
}

// ================================== Solve ================================

int cUnsymSkylMatrix :: Solve(cVector &b, cVector &x)
{
  x = b;

  // Perform the Crout's factorization [A] = [L][U], where Uii = 1.

  if (!dec)
  {
    for (int i = 1; i < dim; i++)
    {
      for (int j = skl[i]; j < i; j++)
      {
	for (int k = MAX(skl[i], skl[j]); k < j; k++)
	{
	  L[i][j] -= L[i][k]*U[j][k]; // A[i][j] -= A[i][k]*A[k][j];
	  U[i][j] -= L[j][k]*U[i][k]; // A[j][i] -= A[j][k]*A[k][i];
	}
	if (fabs(L[j][j]) < TOLZ)
        {
          cout << "Zero pivot in LU decomposition, line = " << j << ".\n";
          return(0);
        }
        U[i][j] /= L[j][j];          // A[j][i] /= A[j][j];
        L[i][i] -= L[i][j]*U[i][j];  // A[i][i] -= A[i][j]*A[j][i];
      }
    }
    dec = 1;
  }

  // Forward substitution.

  int i,j;
  for (i = 0; i < dim; i++)
  {
    for (j = skl[i]; j < i; j++) x[i] -= L[i][j]*x[j];
    x[i] /= L[i][i];
  }

  // Back-substitution.

#if 1
  for (i = dim-1; i >= 0; i--)
  {
    for (j = skl[i]; j < i; j++) x[j] -= U[i][j]*x[i];
  }
#else
  for (i = dim-1; i >= 0; i--)
  {
    for (j = i+1; j < dim; j++) x[i] -= U[j][i]*x[j];
  }
#endif
  return(1);
}

// =============================== GetSmlPvt ===============================

int cUnsymSkylMatrix :: GetSmlPvt(void)
{
  int    idx  = 0;
  double mink = fabs(L[0][0]);

  for (int i = 1; i < dim; i++)
  {
    double curk = fabs(L[i][i]);
    if (curk < mink)
    {
      idx  = i;
      mink = curk;
    }
  }

  return(idx);
}

// =============================== GetLgtPvt ===============================

int cUnsymSkylMatrix :: GetLgtPvt(void)
{
  int    idx  = 0;
  double maxk = fabs(L[0][0]);

  for (int i = 1; i < dim; i++)
  {
    double curk = fabs(L[i][i]);
    if (curk > maxk)
    {
      idx  = i;
      maxk = curk;
    }
  }

  return(idx);
}

// =============================== GetNgtPvt ===============================

int cUnsymSkylMatrix :: GetNgtPvt(void)
{
  int nnp = 0;

  for (int i = 0; i < dim; i++) if (L[i][i] < 0.0) nnp++;

  return(nnp);
}

// ================================== Get ==================================

double cUnsymSkylMatrix :: Get(int i, int j)
{
  if (IsSup(i,j)) return(U[j][i]);

  return(L[i][j]);
}

// ================================== Add ==================================

void cUnsymSkylMatrix :: Add(int i, int j, double v)
{
  if (IsSup(i,j))
    U[j][i] += v;
  else
    L[i][j] += v;
}

// ================================= Print =================================

void cUnsymSkylMatrix :: Print(void)
{
  int i,j;

  cout << "Stored elements = " << nelm << "\n";
  cout << showpos << scientific << setprecision(2);
  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < skl[i]; j++) cout << 0.0 << "  "; 
    for (j = skl[i]; j <= i; j++) cout << L[i][j] << "  ";
    for (j = i+1; j < dim; j++)
    {
      if (j >= skl[i])
	cout << U[i][j] << "  ";
      else
	cout << 0.0 << "  ";
    }
    cout << "\n";
  }
  cout << noshowpos << fixed;
}

// =============================== MultVect ================================

void cUnsymSkylMatrix :: MultVect(double *u, double *v)
{
  int i,j;

  // Initialize the given vector

  for (i = 0; i < dim; i++) v[i] = 0.0;

  // Compute the matrix/vector product

  for (i = 0; i < dim; i++)
  {
    for (j = skl[i]; j < i; j++)
    {
      v[i] += L[i][j]*u[j];
      v[j] += U[i][j]*u[i];       // v[j] += A[j][i]*u[i];
    }
    v[i] += L[i][i]*u[i];
  }
}

// ================================ AddMat =================================

void cUnsymSkylMatrix :: AddMat(double a, cSysMatrix *M)
{
  for (int i = 0; i < dim; i++)
  {
    for (int j = skl[i]; j < i; j++)
    {
      L[i][j] += a*M->Get(i, j);
      U[i][j] += a*M->Get(j, i);  // A[j][i] += a*M->Get(j, i);
    }
    L[i][i] += a*M->Get(i, i);
  }
}

// -------------------------------------------------------------------------
// cSprsMatrix class:
// -------------------------------------------------------------------------

// ============================== cSprsMatrix ==============================

cSprsMatrix :: cSprsMatrix(int n) : cSysMatrix(n)
{
  dim  = n;
  nelm = 0;

  row = new SpElm*[dim];
  for (int i = 0; i < dim; i++) row[i] = 0;
}

// ============================= ~cSprsMatrix ==============================

cSprsMatrix :: ~cSprsMatrix(void)
{
  SpElm *e,*l;

  for (int i = 0; i < dim; i++)
  {
    e = row[i];
    while (e)
    {
      l = e;
      e = e->nxt;
      delete l;
    }
  }
  delete []row;
}

// ================================== Add ==================================

void cSprsMatrix :: Add(int i, int j, double v)
{
  // Ignore a zero value

  if (v == 0.0) return;

  // Search the given column

  SpElm *e,*l;
  for (e = row[i], l = 0; e != 0 && e->col < j; e = e->nxt) l = e;

  // Add to an existing element

  if (e != 0 && e->col == j)
  {
    e->val += v;
    return;
  }

  // Create a new element

  SpElm *c = new SpElm;
  c->col = j;
  c->val = v;
  nelm++;

  // Add to the list

  if (!l)       // Add on top
  {
    c->nxt = row[i];
    row[i] = c;
  }
  else if (!e)  // Add on end
  {
    l->nxt = c;
    c->nxt = 0;
  }
  else          // Add between two elements
  {
    c->nxt = l->nxt;
    l->nxt = c;
  }
}

// ================================== Get ==================================

double cSprsMatrix :: Get(int i, int j)
{
  for (SpElm *e = row[i]; e != 0; e = e->nxt)
  {
    if (e->col == j) return(e->val);
  }

  return(0.0);
}

// ================================= Print =================================

void cSprsMatrix :: Print(void)
{
  int   j,fst,lst;
  SpElm *e;

  cout << "Nonzero elements = " << nelm << "\n";
  cout << showpos << scientific << setprecision(2);
  for (int i = 0; i < dim; i++)
  {
    fst = 0;
    for (e = row[i]; e != 0; e = e->nxt)
    {
      lst = e->col;
      for (j = fst; j < lst; j++) cout << 0.0 << "  ";
      cout << e->val << "  ";
      fst = lst+1;
    }
    for (j = fst; j < dim; j++) cout << 0.0 << "  ";
    cout << "\n";
  }
  cout << noshowpos << fixed << setprecision(2);
}

// =============================== MultVect ================================

void cSprsMatrix :: MultVect(double *u, double *v)
{
  SpElm *e;

  for (int i = 0; i < dim; i++)
  {
    v[i] = 0.0;
    for (e = row[i]; e != 0; e = e->nxt)
    {
      v[i] += e->val*u[e->col];
    }
  }
}


// -------------------------------------------------------------------------
// cSymSprsMatrix class:
// -------------------------------------------------------------------------

// ============================== cSprsMatrix ==============================

cSymSprsMatrix :: cSymSprsMatrix(int n) : cSprsMatrix(n)
{
}

// ================================== Add ==================================

void cSymSprsMatrix :: Add(int i, int j, double v)
{
  if (IsSup(i,j)) return;

  cSprsMatrix :: Add(i, j, v);
}

// ================================== Get ==================================

double cSymSprsMatrix :: Get(int i, int j)
{
  if (IsSup(i,j))
  {
    return(cSprsMatrix :: Get(j, i));
  }
  else
  {
    return(cSprsMatrix :: Get(i, j));
  }
}

// ================================= Print =================================

void cSymSprsMatrix :: Print(void)
{
  int   j,fst,lst;
  SpElm *e;

  cout << "Nonzero elements = " << nelm << "\n"; 
  cout << showpos << scientific << setprecision(2);
  for (int i = 0; i < dim; i++)
  {
    fst = 0;
    for (e = row[i]; e != 0; e = e->nxt)
    {
      lst = e->col;
      for (j = fst; j < lst; j++) cout << 0.0 << "  ";
      cout << e->val << "  ";
      fst = lst+1;
    }
    for (j = fst; j <= i; j++) cout << 0.0 << "  ";
    cout << "\n";
  }
}

// =============================== MultVect ================================

void cSymSprsMatrix :: MultVect(double *u, double *v)
{
  int   i;
  SpElm *e;

  // Initialize the given vector

//  for (i = 0; i < dim; i++) v[i] = 0.0;
  memset(v, '\0', dim*sizeof(double));

  // Compute the matrix/vector product

  for (i = 0; i < dim; i++)
  {
    for (e = row[i]; e->nxt != 0; e = e->nxt)
    {
      v[i] += e->val*u[e->col];  // Inferior elements
      v[e->col] += e->val*u[i];  // Superior elements
    }
    v[i] += e->val*u[i];         // Diagonal element
  }
}

// ================================== Solve ================================

int cSymSprsMatrix :: Solve(cVector &b, cVector &x)
{
  return(PCGSolver(maxit, tol, this, b, x));
}

// -------------------------------------------------------------------------
// Generic Solvers:
//

// =============================== PCGSolver ===============================

int PCGSolver(int maxiter, double tol, cSysMatrix *K, cVector &f,
              cVector &u)
{
  // Check the given data

  double lenf = f.Length( );
  cout << "lenf = " << showpos << scientific << lenf <<"\n";
  if (lenf < tol) return(1);
  if (lenf > 10e+10) lenf = 1.0;           // temporario
//  if (lenf > pow(10.0, log10(BIGK)/2.0)) lenf = 1.0;

  // Get memory for local vectors.

  int n = f.Dim( );
  cout << noshowpos << fixed << "n = " << n << "\n";
  cVector r(n);   // residual vector
  cVector p(n);   // search direction
  cVector d(n);   // preconditioner
  cVector w(n);   // {w} = [K]*{p}
  cVector z(n);   // {z} = inv([D]){r}

  // Evaluate the initial residual and search direction

  for (int i = 0; i < n; i++) d[i] = K->Get(i, i);
  u.Zero( );
  r = f;
  z = r/d;
  p = z;
  double dot0 = z*r;

  // Iteration loop

  int conv = 0;
  double alpha,beta;
  double dot1,err;
  for (int k = 1; k <= maxiter; k++)
  {
    w = K*p;
    alpha = dot0/(p*w);
    u = u + alpha*p;
    r = r + (-alpha)*w;     // temporario !
    z = r/d;

    err = r.Length( )/lenf;
    if (err < tol)
    {
      cout << "iter = " << setprecision(4) << k << "  err = ";
      cout << scientific << err << "\n" << fixed;
      conv = 1;
      break;
    }
    if (!(k % 10))
    {
      cout << "iter = " << setprecision(4) << k << "  err = ";
      cout << scientific << err << "\n" << fixed;
    }

    dot1 = z*r;
    beta = dot1/dot0;
    dot0 = dot1;
    p = z + beta*p;
  }

  return(conv);
}

// =============================== FstEigPair ==============================

int FstEigPair(int maxit,double tol, cSysMatrix *K, double *lbd, cVector &v)
{
  // Initialization

  cVector a(v.Dim( ));
  a = 1.0;
  double l0 = 0.0;
  double l1 = 1.0;

  // Inverse iteration loop

  double dl;
  for (int i = 0; i < maxit; i++)
  {
    // Compute new eigenpair

    K->Solve(a, v);          // Solve [K]{v} = {a}
    l1 = (a*v)/(v*v);
    v /= v.Length( );

    // Check convergence

    dl = fabs(l1 - l0)/(1.0 + fabs(l1));
    if (dl < tol && i > 0)
    {
      *lbd = l1;   // Return the correct eigenvalue
      return(1);
    }

    // Store old values

    a  = v;
    l0 = l1;
  }

  cout << "Convergence not achived in inverse iteration\n";
  cout << scientific << "dl = " << dl << "\n";
  cout << "tol = " << tol << "\n" << fixed;

  *lbd = l1;       // Return the 'best' approximation
  return(0);
}

// ======================================================= End of file =====
