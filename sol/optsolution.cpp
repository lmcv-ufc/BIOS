// -------------------------------------------------------------------------
// optsolution.cpp - implementation of the individual class.
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
// Created:      01-Oct-2014    Elias Saraiva Barroso
//
// Modified:     10-Nov-2014    Elias Saraiva Barroso
//               Created cOptSolution, the main base class of all optimization
//               agents.
//
// Modified:     07-Abr-2016    Elias Saraiva Barroso
//               Created the pure virtual method CompVar.
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "optsolution.h"
#include "problem.h"
#include "individual.h"
#include "particle.h"
#include "food.h"
#include "utl.h"
#include "input.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Auxiliary methods:
//

// ============================ operator < =================================

istream& operator>>(istream &in,eSolType &t)
{
  // Read the individual type label.

  char label[100];
  if (!Utl::ReadString(in, label))
    Utl::Exit("Error in the input of the optimization solution type label.");

  // Set the appropriate individual type.

  if (string(label) == "IntegerVector")
    t = SOL_INT_VEC;
  else if (string(label) == "IntegerMatrix")
    t = SOL_INT_MAT;
  else if (string(label) == "DoubleVector")
    t = SOL_DBL_VEC;
  else if (string(label) == "BinaryVector")
    t = SOL_BIN_VEC;
  else
    Utl::Exit("Unknown optimization solution type: " + string(label));

  return in;
}

// -------------------------------------------------------------------------
// Local methods:
//

// ============================ operator < =================================

bool cOptSolution :: operator<(cOptSolution& a)
{
  if (Prob->GetNumConstr( ) == 0)
    return(GetObjFunc(0) < a.GetObjFunc(0) );
  else
    return(GetPenObjFunc( ) < a.GetPenObjFunc() );
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cOptSolution ====================================

cOptSolution :: cOptSolution(cProblem *prob)
{
  Prob = prob;

  Fobjs.Resize(prob->GetNumObj( ));
  Fobjs.Zero( );

  if (Prob->GetNumObj( ) < 2)                   // Single-objective problem
  {
      PenObjFunc = 0.0;
      FitFuncVal = 0.0;
  }

  int nc = Prob->GetNumConstr( );
  Constr.Resize(nc);
  Constr.Zero( );
}

// ============================ ~cOptSolution ====================================

cOptSolution :: ~cOptSolution(void)
{
}

// ============================ LamMutate ==================================

void cOptSolution :: LamMutate(double *p)
{
  Utl::Exit("cOptSolution::LamMutate is not implemented in this Derived object.");
}

// ============================ Swap =======================================

void cOptSolution :: Swap(double p)
{
  Utl::Exit("cOptSolution::Swap is not implemented in this Derived object.");
}

// ============================ Add ========================================

void cOptSolution :: Add(double p)
{
  Utl::Exit("cOptSolution::Add is not implemented in this Derived object.");
}

// ============================ Delete =====================================

void cOptSolution :: Delete(double p)
{
  Utl::Exit("cOptSolution::Delete is not implemented in this Derived object.");
}

// ============================== GetVarType ===============================

eVarType cOptSolution :: GetVarType(const eSolType &st)
{
  eVarType type = DISCRETE;

  switch(st)
  {
    case SOL_INT_VEC:
    case SOL_INT_MAT:
    case SOL_BIN_VEC:
      type = DISCRETE;
    break;

    case SOL_DBL_VEC:
      type = CONTINUOUS;
    break;
  }

  return(type);
}

// ============================== GetNormVar ===============================

void cOptSolution :: GetNormVar(cVector &x)
{
  Utl::Exit("cOptSolution::GetNormVar is not implemented in this Derived object.");
}

// ============================== Clone ====================================

cOptSolution* cOptSolution :: Clone(void)
{
  cOptSolution* clone = this->NewObj( );
  clone->Copy(this);
  
  return clone;
}

// ======================================================= End of file =====
