// -------------------------------------------------------------------------
// vec.cpp - implementation of the vector class.
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
// Created:   13-Jan-1998    Evandro Parente Junior
//
// Modified:  02-Aug-2000    Evandro Parente Junior
//            Use of new/delete instead of calloc/free. Use of assert.
//
// Modified:  13-Sep-2001    Evandro Parente Junior
//            Use of compositors to improve the efficiency of bin. oper.
//
// Modified:  14-Jun-2002    Evandro Parente Junior
//            Added to femoop.
//
// Modified:  31-Oct-2012    Iuri Barcelos Rocha
//            Added Max() and Min()
//
// Modified:  02-Oct-2013    Evandro Parente Junior
//            Created CrossProd and Normalize.
//
// Modified:  08-Apr-2015     Elias Saraiva Barroso
//            Input/output operations using c++ streams.
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <assert.h>

using namespace std;

#include "vec.h"
#include "mat.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cVector ==================================

cVector :: cVector(void)
{
  dim = 0;
  val = 0;
}

// ============================== cVector ==================================

cVector :: cVector(int n)
{
  dim = n;
  val = new double[n];
}

// ============================== cVector ==================================

cVector :: cVector(int n, double *vec)
{
  dim = n;
  val = new double[n];
  for (int i = 0; i < dim; i++) val[i] = vec[i];
}

// ============================== cVector ==================================

cVector :: cVector(const cVector &vec)
{
  dim = vec.dim;
  val = new double[dim];
  for (int i = 0; i < dim; i++) val[i] = vec.val[i];
}

// ============================== ~cVector =================================

cVector :: ~cVector(void)
{
  delete []val;
}

// ============================== Normalize ================================

void cVector :: Normalize(void)
{
  double len = VecLen(dim, val);

  if (len != 0)
    for (int i = 0; i < dim; i++) val[i] /= len;
}

// ================================ Sort ===================================

inline int CmpAsc(const void *a, const void *b)
{
  double da = *(double *)a;
  double db = *(double *)b;
  if (da > db) return(1);
  if (da < db) return(-1);
  return(0);
}

inline int CmpDes(const void *a, const void *b)
{
  double da = *(double *)a;
  double db = *(double *)b;
  if (da > db) return(-1);
  if (da < db) return(1);
  return(0);
}

void cVector :: Sort(eVecSort order)
{
  if (order == ASCENDING)
    qsort(val, dim, sizeof(double), CmpAsc);
  else
    qsort(val, dim, sizeof(double), CmpDes);
}

// ================================ Max ====================================

double cVector :: Max(void)
{
  double max = val[0];
  
  for (int i = 1; i < dim; i++)
    if (val[i] > max) max = val[i];
  
  return max;
}

// ================================ Min ====================================

double cVector :: Min(void)
{
  double min = val[0];
  
  for (int i = 1; i < dim; i++)
    if (val[i] < min) min = val[i];
  
  return min;
}

// =============================== Resize ==================================

void cVector :: Resize(int n)
{
  delete []val;
  dim = n;
  val = new double[dim];
}

// =============================== Print ===================================

void cVector :: Print(void)
{
  for (int i = 0; i < dim; i++) cout << " " << val[i];
  cout << "\n";
}

// ============================= operator*= ================================

void cVector :: operator*=(const double &a)
{
  for (int i = 0; i < dim; i++) val[i] *= a;
}

// ============================= operator/= ================================

void cVector :: operator/=(const double &a)
{
  assert(a != 0.0);

  for (int i = 0; i < dim; i++) val[i] /= a;
}

// ============================== CrossProd ================================

void CrossProd(cVector &u, cVector &v, cVector &w)
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

// ================================================ End of file ============
