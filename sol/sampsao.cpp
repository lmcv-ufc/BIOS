// -------------------------------------------------------------------------
// sampsao.cpp - implementation of the particle class.
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
// Created:      08-Oct-2013    Elias Saraiva Barroso
//
// Modified:
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "sampsao.h"
#include "sao.h"
#include "problem.h"
#include "lam.h"
#include "utl.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSampSAO ==================================

cSampSAO :: cSampSAO(cProblem *prob) : cOptSolution(prob)
{
}

// ============================ ~cSampSAO ==================================

cSampSAO :: ~cSampSAO(void)
{
}

// ============================== CreateSamp ===============================
cSampSAO* cSampSAO :: CreateSamp(eSolType type, cProblem* prob, sProbAppOut &out)
{
  
  cSampSAO *smp = 0;

  switch(type)
  {

    case SOL_INT_VEC:
      smp = new cSampIntVec(prob);
    break;

    case SOL_INT_MAT:
      smp = new cSampIntMat(prob);
    break;

    case SOL_DBL_VEC:
      smp = new cSampDblVec(prob);
    break;

    case SOL_BIN_VEC:
      cout << "Error: Individual type not yet implemented!" << endl;
      exit(0);
    break;
  }

  if (smp) smp->appout = &out; 

  return(smp);
}


// -------------------------------------------------------------------------
// cSampIntVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSampIntVec ==============================
cSampIntVec :: cSampIntVec(cProblem *prob) : cSampSAO(prob)
{
  OptSolType = SOL_INT_VEC;
  int nvar = Prob->GetNumVar( );
  Var = new int[nvar];
}

// ============================ ~cSampIntVec ==============================

cSampIntVec :: ~cSampIntVec(void)
{
  delete [] Var;
}

// ================================ Init ===================================

void cSampIntVec :: Init( )
{
  // Initialize the position and velocity variables with random values.
  int low,upp;
  int nvar = Prob->GetNumVar( );

  // Initialize the variables with random values.
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    Var[i] = Utl::RandInt(low, upp);
  }
}

void cSampIntVec :: Init(const sInpSol &ivar)
{
  // Store input values.
  int nvar = Prob->GetNumVar( );
  for(int i = 0; i < nvar; ++i)
    Var[i] = static_cast<int>(ivar.CodVar[i]);
}

void cSampIntVec :: Init(const cVector &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Initialize variables with given values in [0,1] range.
  int low,upp,nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    Var[i] = static_cast<int>(low + (upp - low + 1.0)*ivar[i]);
    Var[i] = (Var[i] > upp) ? upp : Var[i];
  }
}

// =============================== GetVar ==================================

void cSampIntVec :: GetVar(int *a)
{
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
    a[i] = Var[i];
}

// ============================== Evaluate =================================

void cSampIntVec :: Evaluate(void)
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cSampIntVec :: Print( )
{
  int nvar = Prob->GetNumVar( );

  cout << "Sample data" << endl;
  cout << "NumVar = " << nvar << endl;
  cout << "Variables" << endl;
  for (int i = 0; i < nvar; i++)
    cout << Var[i] << "  ";
  cout << endl;

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "Constraints" << endl;
    Constr.Print( );
  }

  cout << "ObjectiveFunction = ";
  Fobjs.Print( );
  cout << endl;
  if (Prob->GetNumConstr( ) > 0)
    cout << "PenObjFunction    = " << PenObjFunc << endl << endl;
  else
    cout << endl;
}

// =============================== Write ===================================

void cSampIntVec :: Write(ostream &out)
{
  int nvar = Prob->GetNumVar( );

  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  for (int i = 0; i < Fobjs.Dim( ); i++)
    out << Fobjs[i] << endl;

  if (Prob->GetNumConstr() > 0)
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n" << PenObjFunc << endl;

  out << "\n%DESIGN.VARIABLES\n";
  for (int i = 0; i < nvar; i++)
    out << Var[i] << endl;

  if (Prob->GetNumConstr() > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << Constr[i] << endl;
  }
}

// ================================ Copy ===================================

void cSampIntVec :: Copy(cOptSolution *a)
{
  cSampIntVec *smp = (cSampIntVec*) a;

  Fobjs      = smp->Fobjs;
  Constr     = smp->Constr;
  PenObjFunc = smp->PenObjFunc;

  int nv =  Prob->GetNumVar();

  for (int i = 0; i < nv; i++)
    Var[i]  = smp->Var[i];
}

// ================================ CompVar ================================

bool cSampIntVec :: CompVar(cOptSolution *sol)
{
  cSampIntVec *part = dynamic_cast<cSampIntVec*> (sol);

  if (part == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
    if (Var[i] != part->Var[i])
      return false;

  return true;
}

// =============================== GetNormVar ==================================

void cSampIntVec ::  GetNormVar(cVector &x)
{
  Prob->GetNormVar(Var,x);
}

// ================================ GetSurrOutRes ==========================

void cSampIntVec :: GetSurrOutRes(cVector &x)
{
  // Store approximated objective functions.
  int id = 0;
  for(int i = 0; i < appout->GetNumAppObjFunc( ); ++i)
    x[id++] = Fobjs[appout->GetAppObjFuncID(i)];

  // Store approximated constriants.
  for(int i = 0; i < appout->GetNumAppConstr( ); ++i)
    x[id++] = Constr[appout->GetAppConstrID(i)];
}

// -------------------------------------------------------------------------
// cIndivIntMat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSampIntMat ==============================

cSampIntMat :: cSampIntMat(cProblem *prob) : cSampSAO(prob)
{
  OptSolType = SOL_INT_MAT;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  Var = new int*[n];
  for (int i = 0; i < n; i++)
  {
    Var[i] = new int[m];
  }
}

// ============================ ~cPartIntMat ==============================

cSampIntMat :: ~cSampIntMat(void)
{
  int n = Prob->VarNumRow( );
  for (int i = 0; i < n; i++)
  {
    delete []Var[i];
  }
  delete []Var;
}

// ================================ Init ===================================

void cSampIntMat :: Init( )
{
  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
      Var[i][j] = Utl::RandInt(low, upp);
    }
  }
}

void cSampIntMat :: Init(const sInpSol &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Store input values.
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Var[i][j] = static_cast<int>(ivar.CodVar[i*m+j]);
}

void cSampIntMat :: Init(const cVector &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Initialize variables with given values in [0,1] range.

  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  cVector xn = ivar;

  int c = 0;
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);

    if (low == upp)
    {
        for (int j = 0; j < m; j++)
        {
          Var[i][j] = low;
        }
    }
    else
    {
      for (int j = 0; j < m; j++)
      {
        Var[i][j] = static_cast<int>(low + (upp - low + 1)*ivar[c*m+j]);
        Var[i][j] = (Var[i][j] > upp) ? upp : Var[i][j];
      }
      c += 1;
    }
  }
}

// ================================ Copy ===================================

void cSampIntMat :: Copy(cOptSolution *a)
{
  cSampIntMat *smp = (cSampIntMat*) a;

  Fobjs      = smp->Fobjs;
  Constr     = smp->Constr;
  PenObjFunc = smp->PenObjFunc;

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
    {
      Var[i][j] = smp->Var[i][j];
    }
}

// ================================ CompVar ================================

bool cSampIntMat :: CompVar(cOptSolution *sol)
{

  cSampIntMat *part = dynamic_cast<cSampIntMat*> (sol);

  if (part == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if (Var[i][j] != part->Var[i][j])
        return false;

  return true;
}
// ============================== Evaluate =================================

void cSampIntMat :: Evaluate(void)
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cSampIntMat :: Print(void)
{
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  cout << "------------------------------------------------------------" << endl;
  cout << "Sample data" << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) cout << Var[i][j] << "  ";
    cout << endl;
  }

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "Constraints = ";
    Constr.Print( );
  }

  cout << "ObjectiveFunction = ";
  Fobjs.Print( );
  cout << endl;
  if (Prob->GetNumConstr( ) > 0)
    cout << "PenObjFunction = " << PenObjFunc << endl << endl;
  else
    cout << endl;
}

// =============================== Write ===================================

void cSampIntMat :: Write(ostream &out)
{
    out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
    for (int i = 0; i < Fobjs.Dim( ); i++)
      out << Fobjs[i] << endl;

    if (Prob->GetNumConstr() > 0)
      out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n" << PenObjFunc << endl;

    out << "\n%DESIGN.VARIABLES\n";
    int n = Prob->VarNumRow( );
    int m = Prob->VarNumCol( );
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++) out << Var[i][j] << "  ";
      out << endl;
    }

    if (Prob->GetNumConstr() > 0)
    {
      out << "\n%CONSTRAINT.VALUES\n";
      for (int i = 0; i < Constr.Dim( ); i++)
        out << Constr[i] << endl;
    }
}

// =============================== GetNormVar ==================================

void cSampIntMat ::  GetNormVar(cVector &xn)
{
  Prob->GetNormVar(Var,xn);
}

// ================================ GetSurrOutRes ==========================

void cSampIntMat :: GetSurrOutRes(cVector &x)
{
  // Store approximated objective functions.
  int id = 0;

  for(int i = 0; i < appout->GetNumAppObjFunc( ); ++i)
  {
    x[id++] = Fobjs[appout->GetAppObjFuncID(i)];
  }

  // Store approximated constriants.
  for(int i = 0; i < appout->GetNumAppConstr( ); ++i)
    x[id++] = Constr[appout->GetAppConstrID(i)];
}

/*
// ================================ Send ===================================

void cPartIntMat :: Send(void)
{
#ifdef _MPI_
  // Get necessary data.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );
  int nconstr = Prob->GetNumConstr( );

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  // Store data in contiguous memory arrays.

  int varsend[n*m];
  double objsend;
  double constrsend[nconstr];

  int bvarsend[n*m];
  double bobjsend;
  double bconstrsend[nconstr];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
    {
      varsend[m*i + j] = PosVar[i][j];
      bvarsend[m*i + j] = BestPosVar[i][j];
    }

  objsend = Fobjs[0];
  bobjsend = BestObjFuncVal;

  for (int i = 0; i < nconstr; i++)
  {
    constrsend[i] = Constr[i];
    bconstrsend[i] = BestConstr[i];
  }

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, n*m, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(bvarsend, n*m, MPI_INT, j, 2);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 3);
      MPI::COMM_WORLD.Send(&bobjsend, 1, MPI_DOUBLE, j, 4);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 5);
      MPI::COMM_WORLD.Send(bconstrsend, nconstr, MPI_DOUBLE, j, 6);
    }
#endif
}

// =============================== Receive =================================

void cPartIntMat :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nconstr = Prob->GetNumConstr( );

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  // Prepare contiguous memory arrays.

  int varrecv[n*m];
  double objrecv;
  double constrrecv[nconstr];

  int bvarrecv[n*m];
  double bobjrecv;
  double bconstrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, n*m, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(bvarrecv, n*m, MPI_INT, deme, 2);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 3);
  MPI::COMM_WORLD.Recv(&bobjrecv, 1, MPI_DOUBLE, deme, 4);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 5);
  MPI::COMM_WORLD.Recv(bconstrrecv, nconstr, MPI_DOUBLE, deme, 6);

  // Store received values.

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
    {
      PosVar[i][j] = varrecv[m*i + j];
      BestPosVar[i][j] = bvarrecv[m*i + j];
    }

  Fobjs[0] = objrecv;
  BestObjFuncVal = bobjrecv;

  for (int i = 0; i < nconstr; i++)
  {
    Constr[i] = constrrecv[i];
    BestConstr[i] = bconstrrecv[i];
  }

  PenObjFunc = BestPenObjFunc = 0.0;
  FitFuncVal = BestFitFuncVal = 0.0;
#endif
}
*/

// -------------------------------------------------------------------------
// cSampDblVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSampDblVec ===============================

cSampDblVec :: cSampDblVec(cProblem *prob) : cSampSAO(prob) 
{
  OptSolType = SOL_DBL_VEC;
  int nvar = Prob->GetNumVar( );

  Var.Resize(nvar);
  Var.Zero( );
}

// ============================ ~cSampDblVec ===============================

cSampDblVec :: ~cSampDblVec(void)
{
}

// ================================ Init ===================================

void cSampDblVec :: Init(const cVector &ivar)
{
  int nvar = Prob->GetNumVar( );
  double *low,*upp;
  low = new double[nvar];
  upp = new double[nvar];

  // Decode variables from [0,1] to problem bounds.
  Prob->GetDblBounds(low,upp);
  for (int i = 0; i < nvar; i++)
    Var[i] = low[i] + (upp[i]-low[i])*ivar[i];

  delete[] low;
  delete[] upp;
}

// ============================== Evaluate =================================

void cSampDblVec :: Evaluate( )
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cSampDblVec :: Print( )
{
  cout << "Sample data" << endl;
  cout << "NumVar = " << Var.Dim( ) << endl;
  cout << "Variables" << endl;
  Var.Print( );

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "Constraints" << endl;
    Constr.Print( );
  }

  cout << "ObjectiveFunction = ";
  Fobjs.Print( );
  cout << endl;
  if (Prob->GetNumConstr( ) > 0)
    cout << "PenObjFunction    = " << PenObjFunc << endl << endl;
  else
    cout << endl;
}

// =============================== Write ===================================

void cSampDblVec :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  for (int i = 0; i < Fobjs.Dim( ); i++)
    out << Fobjs[i] << endl;

  if (Prob->GetNumConstr() > 0)
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n" << PenObjFunc << endl;
  
  out << "\n%DESIGN.VARIABLES\n";
  for (int i = 0; i < Var.Dim( ); i++)
    out << Var[i] << endl;

  if (Prob->GetNumConstr() > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << Constr[i] << endl;
  }
}

// ================================ Copy ===================================

void cSampDblVec :: Copy(cOptSolution *a)
{
  cSampDblVec *smp = (cSampDblVec*) a;

  Fobjs      = smp->Fobjs;
  Constr     = smp->Constr;
  PenObjFunc = smp->PenObjFunc;

  int nv =  Prob->GetNumVar();

  for (int i = 0; i < nv; i++)
    Var[i]  = smp->Var[i];
}

// ================================ CompVar ================================

bool cSampDblVec :: CompVar(cOptSolution *sol)
{
  cSampDblVec *part = dynamic_cast<cSampDblVec*> (sol);

  if (part == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
    if (Var[i] != part->Var[i])
      return false;

  return true;
}

// =============================== GetNormVar ==================================

void cSampDblVec ::  GetNormVar(cVector &xn)
{
  Prob->GetNormVar(Var,xn);
}

// ================================ GetSurrOutRes ==========================

void cSampDblVec :: GetSurrOutRes(cVector &x)
{
  // Store approximated objective functions.
  int id = 0;
  for(int i = 0; i < appout->GetNumAppObjFunc( ); ++i)
    x[id++] = Fobjs[appout->GetAppObjFuncID(i)];  

  // Store approximated constriants.
  for(int i = 0; i < appout->GetNumAppConstr( ); ++i)
    x[id++] = Constr[appout->GetAppConstrID(i)];  
}

// ======================================================= End of file =====
