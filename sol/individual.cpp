// -------------------------------------------------------------------------
// individual.cpp - implementation of the individual class.
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
// Created:      12-Mar-2013    Evandro Parente Junior
//               Based on ind.cpp
//
// Modified:     10-Sep-2017    Marina Alves Maia
//               Added multiobjective initialization print and write functions
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "individual.h"
#include "problem.h"
#include "optalg.h"
#include "lam.h"
#include "utl.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cIndividual ===============================

cIndividual :: cIndividual(cProblem *prob) : cOptSolution(prob)
{
}

// ============================ ~cIndividual ===============================

cIndividual :: ~cIndividual(void)
{
}

// ============================== CreateIndiv  ===============================

cIndividual* cIndividual :: CreateIndividual(eSolType type, cProblem* prob)
{
  switch(type)
  {
    case SOL_INT_VEC:
      return new cIndivIntVec(prob);
    break;

    case SOL_INT_MAT:
      return new cIndivIntMat(prob);
    break;

    case SOL_DBL_VEC:
      return new cIndivDblVec(prob);
    break;

    case SOL_BIN_VEC:
      return new cIndivBinVec(prob);
    break;
  }

  return 0;
}

// -------------------------------------------------------------------------
// cIndivIntVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cIndivIntVec ==============================

cIndivIntVec :: cIndivIntVec(cProblem *prob) : cIndividual(prob)
{
  OptSolType = SOL_INT_VEC;
  int nvar   = Prob->GetNumVar( );
  Var        = new int[nvar];
}

// ============================ ~cIndivIntVec ==============================

cIndivIntVec :: ~cIndivIntVec(void)
{
  delete [] Var;
}

// ================================ Init ===================================

void cIndivIntVec :: Init(void)
{
  int nvar = Prob->GetNumVar( );

  // Initialize the variables with random values.

  int low,upp;
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    Var[i] = Utl::RandInt(low, upp);
  }
}

void cIndivIntVec :: Init(const sInpSol &ivar)
{
  int nvar = Prob->GetNumVar( );

  // Store input values.
  //
  for(int i = 0; i < nvar; ++i)
    Var[i] = static_cast<int>(ivar.CodVar[i]);
}

void cIndivIntVec :: Init(const cVector &ivar)
{
  int nvar = Prob->GetNumVar( );

  // Initialize variables with given value in [0,1] range.

  int low,upp;
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    Var[i] = static_cast<int>(low + (upp-low)*ivar[i]);
    Var[i] = (Var[i] > upp) ? upp : Var[i]; 
  }
}

// ================================ AddBestSample ==========================

void cIndivIntVec :: AddBestSample(cVector x, double f)
{
    Fobjs[0] = f;

    int nvar = Prob->GetNumVar();
    double step = 1.0e-3;

    double low, upp;

    cVector xt(nvar);

    for (int i = 0; i < nvar; i++)
    {
        Prob->GetVarBounds(i, low, upp);

        xt[i] = x[i]*(upp-low)+low;
        Var[i] = round((xt[i]-low)/step);
    }
}

// ================================ Init ===================================

void cIndivIntVec :: GetBestSample(cVector &x)
{
    int nvar = Prob->GetNumVar();
    cVector xsamp;

    Prob->DecodeVar(Var, xsamp);

    x.Resize(nvar);
    x = xsamp;
}

// =============================== ExploreBounds ==================================

void cIndivIntVec :: ExploreBounds(void)
{
  // Initialize the variables with lower bounds.

  int low,upp;
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    Var[i] = Utl::RandInt(low, low);
  }
}

// =============================== ExploreUpperBounds ==================================

void cIndivIntVec :: ExploreUpperBounds(void)
{
  // Initialize the variables with upper bounds.

  int low,upp;
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    Var[i] = Utl::RandInt(upp, upp);
  }
}

// =============================== GetNormVar ==================================

void cIndivIntVec ::  GetNormVar(cVector &xn)
{
  Prob->GetNormVar(Var,xn);
}

// =============================== GetVar ==================================

void cIndivIntVec :: GetVar(int *a)
{
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    a[i] = Var[i];
  }
}

// ============================== Evaluate =================================

void cIndivIntVec :: Evaluate(void)
{
  // Evaluate the objective function and constraints.

  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cIndivIntVec :: Print(void)
{
  int nv = Prob->GetNumVar( );

  cout << "------------------------------------------------------------" << endl;
  cout << "NumVar = " << nv << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;

  for (int i = 0; i < nv; i++) cout << 1+i << " " << Var[i] << endl;

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (problem)" << endl;
  cout << "------------------------------------------------------------" << endl;
  Prob->PrintVar(Var);

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "------------------------------------------------------------" << endl;
    cout << "Constraints" << endl;
    cout << "------------------------------------------------------------" << endl;
    Constr.Print( );
  }

  cout << "------------------------------------------------------------" << endl;
  if (Prob->GetNumObj( ) > 1)          // Multiobjective
  {
      cout << "Objective Function 1 = " << Fobjs[0] << endl;
      cout << "Objective Function 2 = " << Fobjs[1] << endl;
      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "------------------------------------------------------------" << endl;
          cout << "Sum of normalized constraints = " << NormConst << endl;
      }
  }
  else                                // Single-objective
  {
      cout << "ObjectiveFunction = " << Fobjs[0] << endl;
      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "PenObjFunction    = " << PenObjFunc << endl;

      }
      cout << "------------------------------------------------------------" << endl;
      cout << "FitnessFunction   = " << FitFuncVal << endl;
  }
}

// =============================== Write ===================================

void cIndivIntVec :: Write(ostream &out)
{
    if (Prob->GetNumObj( ) > 1)
    {
        out << "\n%RESULT.OBJECTIVE.FUNCTIONS\n";
        out << Fobjs[0] << endl;
        out << Fobjs[1] << endl;
    }
    else
    {
        out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
        out << Fobjs[0] << endl;
    }

    if (Prob->GetNumConstr( ) > 0)
    {
        if (Prob->GetNumObj( ) > 1)
        {
            out << "\n%RESULT.SUM.OF.NORMALIZED.CONSTRAINTS\n";
            out << NormConst << endl;
        }
        else
        {
            out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
            out << PenObjFunc << endl;
        }
    }

  out << "\n%OPTIMIZATION.VARIABLES\n";
  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++) out << Var[i] << "  ";
  out << endl;

  out << "\n%PROBLEM.VARIABLES\n";
  Prob->WriteVar(Var);

  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << setw(3) << i+1 << "  " << Constr[i] << endl;
  }
}

// ================================ Copy ===================================

void cIndivIntVec :: Copy(cOptSolution *a)
{
  cIndivIntVec *ind = dynamic_cast<cIndivIntVec*>(a);

  if (!a)
  {
      Utl::Exit("Failed dynamic_cast on cIndivIntVec::Copy");
  }

  Fobjs = ind->Fobjs;

  if (Prob->GetNumObj() < 2)
  {
      PenObjFunc = ind->PenObjFunc;
      FitFuncVal = ind->FitFuncVal;
  }

  Constr = ind->Constr;

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++) Var[i] = ind->Var[i];
}

// ================================ CompVar ================================

bool cIndivIntVec :: CompVar(cOptSolution *sol)
{
  cIndivIntVec *ind = dynamic_cast<cIndivIntVec*> (sol);

  if (ind == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
    if (Var[i] != ind->Var[i])
      return false;

  return true;
}

// =============================== Mutate ==================================

void cIndivIntVec :: Mutate(double mut_rate)
{
  int low,upp;
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    double r = Utl::RandDec( );
    if (r <= mut_rate)
    {
      Prob->GetBounds(i, &low, &upp);
      Var[i] = Utl::RandInt(low, upp);
    }
  }
}

// =============================== Hypermutation ==========================

void cIndivIntVec :: Hypermutation(double mut_rate)
{
  int low,upp;
  int nvar = Prob->GetNumVar( );
  int result;
  bool eventCheck = false;

  for (int i = 0; i < nvar; i++)
  {
    double r = Utl::RandDec( );
    if (r <= mut_rate)
    {
      Prob->GetBounds(i, &low, &upp);
      if (upp == 0)
        continue;

      result = Utl::RandInt(low, upp);

      if (result != Var[i])
        eventCheck = true;

      Var[i] = Utl::RandInt(low, upp);
    }
  }

  if (eventCheck == false)
  {
    // Get valid index

    vector<int> inds;

    for (int i = 0; i < nvar; i++)
    {
      Prob->GetBounds(i,&low,&upp);

      if (upp != 0)
        inds.push_back(i);
    }

    int valid_index = Utl::RandInt(0,inds.size()-1);

    int i = inds[valid_index];

    Prob->GetBounds(i,&low,&upp);

    // Get valid result

    do
      result = Utl::RandInt(low,upp);
    while (result == Var[i]);

    Var[i] = result;
  }
}

//=============================== Crossover ================================

void cIndivIntVec :: Crossover(eCrossType type, double r, cOptSolution *p1, cOptSolution *p2)
{
  cIndivIntVec *par1 = dynamic_cast<cIndivIntVec*>(p1);
  cIndivIntVec *par2 = dynamic_cast<cIndivIntVec*>(p2);

  if (!par1 || !par2)
    Utl::Exit("Failed dynamic_cast on cIndivIntVec::Crossover");

  int nvar = Prob->GetNumVar( );
  double s = 1.0 - r;
  double pntid = round((nvar-1)*r);

  switch (type)
  {
    case LINEAR_COMBINATION:
      for (int i = 0; i < nvar; i++)
        Var[i] = round(r*par1->Var[i] + s*par2->Var[i]);
    break;

    case CLASSICAL:
      for (int j = 0; j < pntid; j++)
       Var[j] = par1->Var[j];
      for (int j = pntid; j < nvar; j++)
       Var[j] = par2->Var[j];
    break;

    case BIN:
      for (int i = 0; i < nvar; i++)
      {
        double rr = Utl::RandDouble(0, 1);
        if (rr < r){
          Var[i] = par1->Var[i];
        }
        else {
          Var[i] = par2->Var[i];
        }
    }
    break;
  }
}

// ============================== Differetiation ================================

void cIndivIntVec :: Differentiation(eDifType type, double f, cOptSolution *lastbest, cOptSolution *p1, cOptSolution *p2, cOptSolution *p3)
{
  cIndivIntVec *par1 = dynamic_cast<cIndivIntVec*>(p1);
  cIndivIntVec *par2 = dynamic_cast<cIndivIntVec*>(p2);
  cIndivIntVec *par3 = dynamic_cast<cIndivIntVec*>(p3);

  if (!par1 || !par2 || !par3)
    Utl::Exit("Failed dynamic_cast on cIndivIntVec::Crossover");

  int nvar = Prob->GetNumVar( );

  switch (type)
  {
    case Rand1:
      for (int i = 0; i < nvar; i++)
      {
       Var[i] = (par3->Var[i] + f*(par1->Var[i] - par2->Var[i]));
      }
      break;
  }
}

// ================================ Send ===================================

void cIndivIntVec :: Send(void)
{
#ifdef _MPI_
  // Get necessary data.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );
  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Store data in contiguous memory arrays.

  int varsend[nvar];
  double objsend;
  double constrsend[nconstr];

  for (int i = 0; i < nvar; i++)
    varsend[i] = Var[i];

  objsend = ObjFuncVal;

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, nvar, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 3);
    }
#endif
}

// =============================== Receive =================================

void cIndivIntVec :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Prepare contiguous memory arrays.

  int varrecv[nvar];
  double objrecv;
  double constrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, nvar, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 3);

  // Store received values.

  for (int i = 0; i < nvar; i++)
    Var[i] = varrecv[i];

  ObjFuncVal = objrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// -------------------------------------------------------------------------
// cIndivIntMat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cIndivIntMat ==============================

cIndivIntMat :: cIndivIntMat(cProblem *prob) : cIndividual(prob)
{
  OptSolType = SOL_INT_MAT;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  Var = new int*[n];
  for (int i = 0; i < n; i++) Var[i] = new int[m];
}

// ============================ ~cIndivIntMat ==============================

cIndivIntMat :: ~cIndivIntMat(void)
{
  int n = Prob->VarNumRow( );
  for (int i = 0; i < n; i++) delete []Var[i];
  delete []Var;
}

// ================================ Init ===================================

void cIndivIntMat :: Init( )
{
  // Initialize the variables with random values.

  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
      Var[i][j] = Utl::RandInt(low, upp);
  }
}

void cIndivIntMat :: Init(const sInpSol &ivar)
{
  // Store input values.
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Var[i][j] = static_cast<int>(ivar.CodVar[i*m+j]);
}

void cIndivIntMat :: Init(const cVector &ivar)
{
  // Initialize variables with given values in [0,1] range.

  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
      Var[i][j] = static_cast<int>(low + (upp+1-low)*ivar[i*m+j]);
      Var[i][j] = (Var[i][j] > upp) ? upp : Var[i][j];
    }
  }
}

// ================================ AddBestSample ===================================

void cIndivIntMat :: AddBestSample(cVector x, double f)
{
    Fobjs[0] = f;

    int m = Prob->GetNumVar();

    for (int i = 0; i < m; i++) Var[1][i] = x[i]*2;

}

// ================================ Copy ===================================

void cIndivIntMat :: Copy(cOptSolution *a)
{
  cIndivIntMat *ind = dynamic_cast<cIndivIntMat*>(a);

  if (!ind)
    Utl::Exit("Failed dynamic_cast on cIndivIntMat::Copy");

  Fobjs = ind->Fobjs;

  if (Prob->GetNumObj( ) < 2)
  {
      PenObjFunc = ind->PenObjFunc;
      FitFuncVal = ind->FitFuncVal;
  }

  Constr = ind->Constr;

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Var[i][j] = ind->Var[i][j];
}

// ================================ CompVar ================================

bool cIndivIntMat :: CompVar(cOptSolution *sol)
{
  cIndivIntMat *ind = dynamic_cast<cIndivIntMat*> (sol);

  if (ind == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if (Var[i][j] != ind->Var[i][j])
        return false;

  return true;
}

// ============================== Evaluate =================================

void cIndivIntMat :: Evaluate(void)
{
  // Evaluate the objective function and constraints.

  Prob->Evaluate(Var, Constr, Fobjs);

}

// =============================== Print ===================================

void cIndivIntMat :: Print(void)
{
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) cout << Var[i][j] << endl;
  }

  cout << endl;

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (problem)" << endl;
  cout << "------------------------------------------------------------" << endl;
  Prob->PrintVar(Var);

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "------------------------------------------------------------" << endl;
    cout << "Constraints" << endl;
    cout << "------------------------------------------------------------" << endl;
    Constr.Print( );
  }

  cout << "------------------------------------------------------------" << endl;
  if (Prob->GetNumObj( ) > 1)                             // Multiobjective
  {
      cout << "Objective Function 1 = " << Fobjs[0] << endl;
      cout << "Objective Function 2 = " << Fobjs[1] << endl;

      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "------------------------------------------------------------" << endl;
          cout << "Sum of normalized constraints = " << NormConst << endl;
      }
  }
  else                                                  // Single-objective
  {
      cout << "ObjectiveFunction = " << Fobjs[0] << endl;

      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "PenObjFunction    = " << PenObjFunc << endl;
      }
      cout << "------------------------------------------------------------" << endl;
      cout << "FitnessFunction   = " << FitFuncVal << endl;
  }

}

// ================================ Init ===================================

void cIndivIntMat :: GetBestSample(cVector &x)
{
    cVector xsamp;
    cMatrix layup;
    int nvar = Prob->VarNumCol();
    x.Resize(nvar);

    Prob->DecodeVar(Var, xsamp, layup);

    x = xsamp;
}

// =============================== Write ===================================

void cIndivIntMat :: Write(ostream &out)
{
    if (Prob->GetNumObj() > 1)
    {
        out << "\n%RESULT.OBJECTIVE.FUNCTIONS\n";
        out << Fobjs[0] << endl;
        out << Fobjs[1] << endl;
    }
    else
    {
        out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
        out << Fobjs[0] << endl;
    }

    if (Prob->GetNumConstr( ) > 0)
    {
        if (Prob->GetNumObj( ) > 1)
        {
            out << "\n%RESULT.SUM.OF.NORMALIZED.CONSTRAINTS\n";
            out << NormConst << endl;
        }
        else
        {
            out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
            out << PenObjFunc << endl;
            out << FitFuncVal << endl;
        }
    }

  out << "\n%OPTIMIZATION.VARIABLES\n";
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) out << Var[i][j] << "  ";
    out << endl;
  }

  out << "\n%PROBLEM.VARIABLES\n";
  Prob->WriteVar(Var);

  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << setw(3) << i+1 << "  " << Constr[i] << endl;
  }
}

// =============================== Mutate ==================================

void cIndivIntMat :: Mutate(double mut_rate)
{
  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);

    for (int j = 0; j < m; j++)
    {
      double r = Utl::RandDec( );
      if (r <= mut_rate)
      {
          if (i == 2)
          {
                Var[i][j] = Utl::RandInt(low, upp);
          //    Var[i][j] = Utl::RandInt(low, 1);
          }
          else
          {
              Var[i][j] = Utl::RandInt(low, upp);
          }
      }
    }
  }
}

// =============================== ExploreBounds ==================================

void cIndivIntMat :: ExploreBounds(void)
{
  // Initialize the variables with lower bounds.

  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  for (int i = 2; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
        Var[i][j] = Utl::RandInt(low, low);
    }
  }
}

// =============================== ExploreUpperBounds ==================================

void cIndivIntMat :: ExploreUpperBounds(void)
{
   // Initialize the variables with upper bounds.

  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  for (int i = 2; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
        Var[i][j] =  Utl::RandInt(upp, upp); //1
    }
  }
}

// =============================== LamMutate ===============================

void cIndivIntMat :: LamMutate(double *mut_rate)
{
  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
      double r = Utl::RandDec( );
      if (r <= mut_rate[i]) Var[i][j] = Utl::RandInt(low, upp);
    }
  }
}

// =============================== Hypermutation ==========================

void cIndivIntMat :: Hypermutation(double mut_rate)
{
  int low, upp, result;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  bool eventCheck = false;    // Check if a valid change was peformed.

  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);

    if (upp == 0)
      continue;

    for (int j = 0; j < m; j++)
    {
      double r = Utl::RandDec( );
      if (r <= mut_rate && Var[0][j] != -1)
      {
        result = Utl::RandInt(low, upp);

	if (result != Var[i][j])
	  eventCheck = true;

        Var[i][j] = result;
      }
    }
  }

  if (eventCheck == false)
  {
    // Get valid row index

    vector<int> inds;

    for(int i = 0; i < n; i++)
    {
      Prob->GetBounds(i,&low,&upp);
      if (upp != 0)
        inds.push_back(i);
    }

    int valid_index = Utl::RandInt(0,inds.size()-1);

    int i = inds[valid_index];
    int j = Utl::RandInt(0,m-1);

    Prob->GetBounds(i,&low,&upp);

    if (Var[0][j] != -1)
    {
      // Get valid result

      do
        result = Utl::RandInt(low , upp);
      while (result == Var[i][j]);

      Var[i][j] = result;
    }

  }
}

// ================================= Swap ==================================

void cIndivIntMat :: Swap(double swp_rate)
{
  // Required data.

  cLaminated *lamprob = dynamic_cast<cLaminated*>(Prob->GetHFP( ));

  if (!lamprob)
    Utl::Exit("Failed dynamic_cast on cIndivIntMat::Swap");

  int n = lamprob->VarNumRow( );
  int m = lamprob->VarNumCol( );
  int *TempVar = new int[n];

  // Check if swap is possible.

  int nonnull = 0;
  for (int j = 0; j < m; j++)
    if (lamprob->GetThk(Var[0][j]) != 0.0)
      nonnull++;

  // Perform the swap operation.

  if (nonnull >= 2)
    for (int j = 0; j < m; j++)
      if (lamprob->GetThk(Var[0][j]) != 0.0)
      {
        double r = Utl :: RandDec( );
        if (r <= swp_rate)
        {
          for (int k = 0; k < n; k++)
            TempVar[k] = Var[k][j];

          while(1)
          {
            int k = Utl :: RandInt(0, m-1);
            if (k != j && lamprob->GetThk(Var[0][j]) != 0.0)
            {
              for (int l = 0; l < n; l++)
              {
                Var[l][j] = Var[l][k];
                Var[l][k] = TempVar[l];
              }
              break;
            }
          }
        }
      }
  delete[] TempVar;
}

// ================================= Add ===================================

void cIndivIntMat :: Add(double add_rate)
{
  // Required data.

  int m = Prob->VarNumCol( );

  // Check if addition is possible.

  int null = 0;
  for (int j = 0; j < m; j++)
    if (Var[0][j] == 0)
      null++;

  if (null == 0)
    return;

  // Perform the addition operation.

  int low, upp;
  Prob->GetBounds(0, &low, &upp);

  for (int k = 0; k < m; k++)
    if (Var[0][k] == 0)
    {
      double r = Utl :: RandDec( );
      if (r <= add_rate) Var[0][k]++;
    }
}

// ================================ Delete =================================

void cIndivIntMat :: Delete(double del_rate)
{
  // Required data.

  int m = Prob->VarNumCol( );

  // Check if deletion is possible.

  int nonnull = 0;
  for (int j = 0; j < m; j++)
    if (Var[0][j] != 0)
      nonnull++;

  if (nonnull == 0)
    return;

  // Perform the deletion operation.

  for (int k = 0; k < m; k++)
    if (Var[0][k] != 0)
    {
      double r = Utl :: RandDec( );
      if (r <= del_rate) Var[0][k]--;
    }
}

// ============================== Crossover ================================

void cIndivIntMat :: Crossover(eCrossType type, double r, cOptSolution *p1, cOptSolution *p2)
{
  cIndivIntMat *par1 = dynamic_cast<cIndivIntMat*>(p1);
  cIndivIntMat *par2 = dynamic_cast<cIndivIntMat*>(p2);

  if (!par1 || !par2)
    Utl::Exit("Failed dynamic_cast on cIndivIntMat::Crossover");

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  double s = 1.0 - r;
  double pntid = round((m-1)*r);

  switch (type)
  {
    case LINEAR_COMBINATION:
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < m; j++)
         Var[i][j] = round(r*par1->Var[i][j] + s*par2->Var[i][j]);
      }
    break;

    case CLASSICAL:
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < pntid; j++)
         Var[i][j] = par1->Var[i][j];

        for (int j = pntid; j < m; j++)
         Var[i][j] = par2->Var[i][j];
      }

    case BIN:
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < m; j++)
        {
          double rr = Utl::RandDouble(0, 1);
          if (rr < r){
            Var[i][j]= par1->Var[i][j];
          }
          else {
            Var[i][j] = par2->Var[i][j];
          }
        }
      }
    break;
  }
}

// ============================== Differetiation ================================

void cIndivIntMat :: Differentiation(eDifType type, double f, cOptSolution *lastbest, cOptSolution *p1, cOptSolution *p2, cOptSolution *p3)
{
  cIndivIntMat *par1 = dynamic_cast<cIndivIntMat*>(p1);
  cIndivIntMat *par2 = dynamic_cast<cIndivIntMat*>(p2);
  cIndivIntMat *par3 = dynamic_cast<cIndivIntMat*>(p3);

  if (!par1 || !par2 || !par3)
    Utl::Exit("Failed dynamic_cast on cIndivIntMat::Crossover");

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  switch (type)
  {
    case Rand1:
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < m; j++)
          Var[i][j] = (par3->Var[i][j] + f*(par1->Var[i][j] - par2->Var[i][j]));
      }
      break;
  }
}
// =============================== GetNormVar ==================================

void cIndivIntMat ::  GetNormVar(cVector &xn)
{
  Prob->GetHFP( )->GetNormVar(Var,xn);
}

// ================================ Send ===================================

void cIndivIntMat :: Send(void)
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

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      varsend[m*i + j] = Var[i][j];

  objsend = ObjFuncVal;

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, n*m, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 3);
    }
#endif
}

// =============================== Receive =================================

void cIndivIntMat :: Receive(int deme)
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

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, n*m, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 3);

  // Store received values.

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Var[i][j] = varrecv[m*i + j];

  ObjFuncVal = objrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// -------------------------------------------------------------------------
// cIndivDblVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cIndivDblVec ==============================

cIndivDblVec :: cIndivDblVec(cProblem *prob) : cIndividual(prob)
{
  OptSolType = SOL_DBL_VEC;
  int nv = Prob->GetNumVar( );
  Var.Resize(nv);
  Var.Zero( );
}

// ============================ ~cIndivDblVec ==============================

cIndivDblVec :: ~cIndivDblVec(void)
{
}

// ================================ Init ===================================

void cIndivDblVec :: Init( )
{
  int nvar = Prob->GetNumVar( );
  double *low,*upp;
  low = new double[nvar];
  upp = new double[nvar];

  // Initialize the variables with random values.
  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
    Var[i] = Utl::RandDouble(low[i], upp[i]);

  delete[] low;
  delete[] upp;
}

void cIndivDblVec :: Init(const sInpSol &ivar)
{
  int nvar = Prob->GetNumVar( );

  // Store input values.
  for(int i = 0; i < nvar; ++i)
    Var[i] = ivar.CodVar[i];
}

void cIndivDblVec :: Init(const cVector &ivar)
{
  int nvar = Prob->GetNumVar( );
  double *low,*upp;
  low = new double[nvar];
  upp = new double[nvar];

  // Map variables from [0,1] to problem bounds.
  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
    Var[i] = low[i] + (upp[i]-low[i])*ivar[i];

  delete[] low;
  delete[] upp;
}

// ================================ AddBestSample ===================================

void cIndivDblVec :: AddBestSample(cVector samplex, double sampley)
{

}

// =============================== GetVar ==================================

void cIndivDblVec :: GetVar(cVector &x)
{
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    x[i] = Var[i];
  }
  cout << "entrou aqui!" << endl;
}

// ================================ ExploreBounds ===================================

void cIndivDblVec :: ExploreBounds(void)
{

  // Initialize the variables with lower bound values.

  double *low,*upp;
  int nvar = Prob->GetNumVar( );

  low = new double[nvar];
  upp = new double[nvar];

  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
    Var[i] = Utl::RandDouble(low[i], low[i]);

  delete[] low;
  delete[] upp;
}

// ================================ ExploreBounds ===================================

void cIndivDblVec :: ExploreUpperBounds(void)
{

  // Initialize the variables with upper bound values.

  double *low,*upp;
  int nvar = Prob->GetNumVar( );

  low = new double[nvar];
  upp = new double[nvar];

  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
    Var[i] = Utl::RandDouble(upp[i], upp[i]);

  delete[] low;
  delete[] upp;
}

// ============================== Evaluate =================================

void cIndivDblVec :: Evaluate(void)
{
  // Evaluate the objective function and constraints.
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cIndivDblVec :: Print(void)
{
  cout << "------------------------------------------------------------" << endl;
  cout << "NumVar = " << Var.Dim( ) << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables" << endl;
  cout << "------------------------------------------------------------" << endl;
  Var.Print( );

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "------------------------------------------------------------" << endl;
    cout << "Constraints" << endl;
    cout << "------------------------------------------------------------" << endl;
    Constr.Print( );
  }

  cout << "------------------------------------------------------------" << endl;
  if (Prob->GetNumObj( ) > 1)          // Multiobjective
  {
      cout << "Objective Function 1 = " << Fobjs[0] << endl;
      cout << "Objective Function 2 = " << Fobjs[1] << endl;
      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "------------------------------------------------------------" << endl;
          cout << "Sum of normalized constraints = " << NormConst << endl;
      }
  }
  else                                // Single-objective
  {
      cout << "ObjectiveFunction = " << Fobjs[0] << endl;
      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "PenObjFunction    = " << PenObjFunc << endl;

      }
      cout << "------------------------------------------------------------" << endl;
      cout << "FitnessFunction   = " << FitFuncVal << endl;
  }
}

// ================================ Init ===================================

void cIndivDblVec :: GetBestSample(cVector &x)
{

}


// =============================== Write ===================================

void cIndivDblVec :: Write(ostream &out)
{
    if (Prob->GetNumObj( ) > 1)
    {
        out << "\n%RESULT.OBJECTIVE.FUNCTIONS\n";
        out << Fobjs[0] << endl;
        out << Fobjs[1] << endl;
    }
    else
    {
        out << "\n%RESULT.OBJECTIVE.FUNCTION\n" << Fobjs[0] << endl;
    }

    if (Prob->GetNumConstr( ) > 0)
    {
        if (Prob->GetNumObj( ) > 1)
        {
            out << "\n%RESULT.SUM.OF.NORMALIZED.CONSTRAINTS\n" << NormConst << endl;
        }
        else
        {
            out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n" << PenObjFunc << endl;
        }
     }

  out << "\n%DESIGN.VARIABLES\n";
  for (int i = 0; i < Var.Dim( ); i++)
    out << Var[i] << endl;
  out << endl;

  if (Prob->GetNumConstr() > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << Constr[i] << endl;
  }
}

// ================================ Copy ===================================

void cIndivDblVec :: Copy(cOptSolution *a)
{
  cIndivDblVec *ind = dynamic_cast<cIndivDblVec*>(a);

  if (!a)
    Utl::Exit("Failed dynamic_cast on cIndivDblVec::Copy");

  Fobjs = ind->Fobjs;

  if (Prob->GetNumObj() < 2)
  {
      PenObjFunc = ind->PenObjFunc;
      FitFuncVal = ind->FitFuncVal;
  }

  Constr = ind->Constr;

  int nv =  Prob->GetNumVar();
  for (int i = 0; i < nv; i++)
    Var[i] = ind -> Var[i];
}

// ================================ CompVar ================================

bool cIndivDblVec :: CompVar(cOptSolution *sol)
{
  cIndivDblVec *ind = dynamic_cast<cIndivDblVec*> (sol);

  if (ind == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
    if (Var[i] != ind->Var[i])
      return false;

  return true;
}

// =============================== GetNormVar ==================================

void cIndivDblVec ::  GetNormVar(cVector &xn)
{
  Prob->GetNormVar(Var,xn);
}

// =============================== Mutate ==================================

void cIndivDblVec :: Mutate(double mut_rate)
{
  double *low,*upp;
  int nvar = Prob->GetNumVar( );

  low = new double[nvar];
  upp = new double[nvar];

  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
  {
    double r = Utl::RandDec( );
    if (r <= mut_rate)
    {
      Var[i] = Utl::RandDouble(low[i], upp[i]);
    //  cout << "\nVar" << i << ": " << Var[i] << endl;
     // if (Var[i] == 0) cout << "MUT VARIAVEL ZERO ===============\n\n\nn" << endl;
    }
  }

  delete[] low;
  delete[] upp;
}

// =============================== Hypermutation ==========================

void cIndivDblVec :: Hypermutation(double mut_rate)
{
  double *low,*upp;
  int nvar = Prob->GetNumVar( );

  low = new double[nvar];
  upp = new double[nvar];

  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  bool eventCheck = false;

  for (int i = 0; i < nvar; i++)
  {
    double r = Utl::RandDec( );
    if (r <= mut_rate)
    {
      Var[i] = Utl::RandDouble(low[i], upp[i]);
      eventCheck = true;
    }
  }

  if (eventCheck == false)
  {
    int i = Utl::RandInt(0,nvar-1);
    Var[i] = Utl::RandDouble(low[i], upp[i]);
  }

  delete[] low;
  delete[] upp;
}

// ============================== Crossover ================================

void cIndivDblVec :: Crossover(eCrossType type, double r, cOptSolution *p1, cOptSolution *p2)
{
  cIndivDblVec *par1 = dynamic_cast<cIndivDblVec*>(p1);
  cIndivDblVec *par2 = dynamic_cast<cIndivDblVec*>(p2);

  if (!par1 || !par2)
    Utl::Exit("Failed dynamic_cast on cIndivDblVec::Crossover");

  int nvar = Prob->GetNumVar( );
  double s = 1.0 - r;
  double pntid = round((nvar-1)*r);

  //cout << "Taxa de Crossover: " << r << endl;

  switch (type)
  {
    case LINEAR_COMBINATION:
      for (int i = 0; i < nvar; i++)
      {
       Var[i] = (r*par1->Var[i] + s*par2->Var[i]);
       break;
      }

    case CLASSICAL:
      for (int j = 0; j < pntid; j++)
       Var[j] = par1->Var[j];
      for (int j = pntid; j < nvar; j++)
       Var[j] = par2->Var[j];
    break;
    case BIN:
      for (int i = 0; i < nvar; i++)
      {
          //cout << "Var " << i << endl;
          //cout << "Par1 Var " << par1->Var[i] << endl;
          //cout << "Par2 Var " << par2->Var[i] << endl;
          double rr = Utl::RandDouble(0, 1);
          //cout << "Random Rate: " << rr << endl;
          if (rr < r){
            //cout << "Par1 ganhou" << endl;
            Var[i] = par1->Var[i];
          }
          else {
            //cout << "Par2 ganhou" << endl;
            Var[i] = par2->Var[i];
          }
      }
      break;
  }
}

// ============================== Differetiation ================================

void cIndivDblVec :: Differentiation(eDifType type, double f, cOptSolution *lastbest, cOptSolution *p1, cOptSolution *p2, cOptSolution *p3)
{
    cIndivDblVec *par1 = dynamic_cast<cIndivDblVec*>(p1);
    cIndivDblVec *par2 = dynamic_cast<cIndivDblVec*>(p2);
    cIndivDblVec *par3 = dynamic_cast<cIndivDblVec*>(p3);
    cIndivDblVec *best = dynamic_cast<cIndivDblVec*>(lastbest);

    if (!par1 || !par2 || !par3 || !best)
      Utl::Exit("Failed dynamic_cast on cIndivDblVec::Differentiation");

    int nvar = Prob->GetNumVar( );

    double *low,*upp;

    low = new double[nvar];
    upp = new double[nvar];

    Prob->GetDblBounds(low,upp);

    switch (type)
    {
      case Rand1:
        for (int i = 0; i < nvar; i++)
        {
         //cout << "Variavel "<<i<<" par1 : " << par1->Var[i] << endl;
         //cout << "Variavel "<<i<<" par2 : " << par2->Var[i] << endl;
         //cout << "Variavel "<<i<<" par3 : " << par3->Var[i] << endl;
         Var[i] = (par3->Var[i] + f*(par1->Var[i] - par2->Var[i]));
         //cout << "Variavel "<<i<<" son  : " << Var[i] << endl;
        }
        break;
    case Loc2Best:
      for (int i = 0; i < nvar; i++)
      {
         Var[i] = (par3->Var[i] + f*(best->Var[i] - par3->Var[i]) + f*(par1->Var[i] - par2->Var[i]));
      }
      break;
    case BestJitter:
      for (int i = 0; i < nvar; i++)
      {
         double fjitter = (f + (1 - 0.9999)*Utl::RandDouble(0, 1));
         Var[i] = (best->Var[i] + fjitter*(par1->Var[i] - par2->Var[i]));
      }
      break;

      //case CLASSICAL:
      //  for (int j = 0; j < pntid; j++)
      //   Var[j] = par1->Var[j];
      //  for (int j = pntid; j < nvar; j++)
      //   Var[j] = par2->Var[j];
      //break;
    }

    // Check bounds

    double fac = 0.5;

    for (int i = 0; i < nvar; i++)
    {
      if(Var[i] > upp[i])
      {
        //cout << "Diminuição de Var " << i << endl;
        Var[i] = upp[i] - fac*(Var[i] - upp[i]);
        //cout << "Novo Var " << i << ": " << Var[i] << endl;
      }
      else if (Var[i] < low[i])
      {
        //cout << "Aumento de Var " << i << endl;
        Var[i] = low[i] - fac*(Var[i] - low[i]);
        //cout << "Novo Var " << i << ": " << Var[i] << endl;
      }
    }


}

// ================================ Send ===================================

void cIndivDblVec :: Send(void)
{
#ifdef _MPI_
  // Get necessary data.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );
  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Store data in contiguous memory arrays.

  double varsend[nvar];
  double objsend;
  double constrsend[nconstr];

  for (int i = 0; i < nvar; i++)
    varsend[i] = Var[i];

  objsend = ObjFuncVal;

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, nvar, MPI_DOUBLE, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 3);
    }
#endif
}

// =============================== Receive =================================

void cIndivDblVec :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Prepare contiguous memory arrays.

  double varrecv[nvar];
  double objrecv;
  double constrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, nvar, MPI_DOUBLE, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 3);

  // Store received values.

  for (int i = 0; i < nvar; i++)
    Var[i] = varrecv[i];

  ObjFuncVal = objrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// -------------------------------------------------------------------------
//cIndivBinVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cIndivBinVec ==============================

cIndivBinVec :: cIndivBinVec(cProblem *prob) : cIndividual(prob)
{
  OptSolType = SOL_BIN_VEC;

  // Number of design variables.

  int nvar = Prob->GetNumVar( );

  // Chromosome size.

  int low,upp;
  int srange;
  Scromo = new int[nvar];

  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    srange = (upp - low) + 1;
    Scromo[i] = round(log2(srange)+0.5);
  }

  // Total number of genes

  Sind = 0;
  for (int i = 0; i < nvar; i++)
  {
    Sind += Scromo[i];
  }

  Var = new int[nvar];
  VarBin = new int[Sind];
}

// ============================ ~cIndivBoolVec ==============================

cIndivBinVec :: ~cIndivBinVec(void)
{
  delete []VarBin;
  delete []Var;
  delete []Scromo;
}

// ================================ Init ===================================

void cIndivBinVec :: Init( )
{
  // Initialize the variables with random values.

  for (int i = 0; i < Sind; i++)
    VarBin[i] = Utl::RandInt(0, 1);
}

void cIndivBinVec :: Init(const sInpSol &ivar)
{
  // Initialize the variables with random values.

  for (int i = 0; i < Sind; i++)
    VarBin[i] = static_cast<int>(ivar.CodVar[i]);
}

// ================================ AddBestSample ===================================

void cIndivBinVec :: AddBestSample(cVector samplex, double sampley)
{

}

// ================================ ExploreBounds ===================================

void cIndivBinVec :: ExploreBounds(void)
{

}

// ================================ ExploreBounds ===================================

void cIndivBinVec :: ExploreUpperBounds(void)
{

}

// ============================== Evaluate =================================

void cIndivBinVec :: Evaluate(void)
{
  // Evaluate the objective function and constraints.

  Decoding( );
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cIndivBinVec :: Print(void)
{
  int nvar = Prob->GetNumVar( );
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;
  Decoding ( );

  for (int i = 0; i < nvar; i++) cout << Var[i] << "  ";
  cout << endl;

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (problem)" << endl;
  cout << "------------------------------------------------------------" << endl;
  Prob->PrintVar(Var);

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "------------------------------------------------------------" << endl;
    cout << "Constraints" << endl;
    cout << "------------------------------------------------------------" << endl;
    Constr.Print( );
  }

  cout << "------------------------------------------------------------" << endl;
  if (Prob->GetNumObj( ) > 1)          // Multiobjective
  {
      cout << "Objective Function 1 = " << Fobjs[0] << endl;
      cout << "Objective Function 2 = " << Fobjs[1] << endl;
      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "------------------------------------------------------------" << endl;
          cout << "Sum of normalized constraints = " << NormConst << endl;
      }
  }
  else                                // Single-objective
  {
      cout << "ObjectiveFunction = " << Fobjs[0] << endl;
      if (Prob->GetNumConstr( ) > 0)
      {
          cout << "PenObjFunction    = " << PenObjFunc << endl;

      }
      cout << "------------------------------------------------------------" << endl;
      cout << "FitnessFunction   = " << FitFuncVal << endl;
  }
}

// ================================ Init ===================================

void cIndivBinVec :: GetBestSample(cVector &x)
{

}


// =============================== Write ===================================

void cIndivBinVec :: Write(ostream &out)
{
  Decoding ( );

  if (Prob->GetNumObj( ) > 1)
  {
      out << "\n%RESULT.OBJECTIVE.FUNCTIONS\n";
      out << Fobjs[0] << endl;
      out << Fobjs[1] << endl;
  }
  else
  {
      out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
      out << Fobjs[0] << endl;
  }

  if (Prob->GetNumConstr( ) > 0)
  {
      if (Prob->GetNumObj( ) > 1)
      {
          out << "\n%RESULT.SUM.OF.NORMALIZED.CONSTRAINTS\n";
          out << NormConst << endl;
      }
      else
      {
          out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
          out << PenObjFunc << endl;
      }
  }

  out << "\n%OPTIMIZATION.VARIABLES\n";
  int nv = Prob->GetNumVar( );
  Decoding ( );
  for (int i = 0; i < nv; i++) out << Var[i] << "  ";
  out << endl;

  out << "\n%PROBLEM.VARIABLES\n";
  Prob->WriteVar(Var);

  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << setw(3) << i+1 << "  " << Constr[i] << endl;
  }
}

// ================================ Copy ===================================

void cIndivBinVec :: Copy(cOptSolution *a)
{
  cIndivBinVec *ind = dynamic_cast<cIndivBinVec*>(a);

  if (!a)
    Utl::Exit("Failed dynamic_cast on cIndivBinVec::Copy");

  Fobjs = ind->Fobjs;

  if (Prob->GetNumObj( ) < 2)
  {
      PenObjFunc = ind->PenObjFunc;
      FitFuncVal = ind->FitFuncVal;
  }

  Constr = ind->Constr;

  Sind = ind->Sind;
  for (int i = 0; i < Sind; i++)
    VarBin[i] = ind->VarBin[i];
}

// ================================ CompVar ================================

bool cIndivBinVec :: CompVar(cOptSolution *sol)
{
  cIndivBinVec *ind = dynamic_cast<cIndivBinVec*> (sol);

  if (ind == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  Sind = ind->Sind;
  for (int i = 0; i < Sind; i++)
    if (VarBin[i] != ind->VarBin[i])
      return false;

  return true;
}

// =============================== Mutate ==================================

void cIndivBinVec :: Mutate(double mut_rate)
{
  for (int i = 0; i < Sind; i++)
  {
    double r = Utl::RandDec( );
    if (r <= mut_rate)
    {
      if (VarBin[i] == 0)
        VarBin[i] = 1;
      else
        VarBin[i] = 0;
    }
  }
}

// ============================== Crossover ================================

void cIndivBinVec :: Crossover(eCrossType type, double r, cOptSolution *p1, cOptSolution *p2)
{
  cIndivBinVec *par1 = dynamic_cast<cIndivBinVec*>(p1);
  cIndivBinVec *par2 = dynamic_cast<cIndivBinVec*>(p2);

  if (!par1 || !par2)
    Utl::Exit("Failed dynamic_cast on cIndivBinVec::Crossover");

  int point = round((Sind - 1)*r);

  for (int i = 0; i < point; i++)
    VarBin[i] = par1->VarBin[i];
  for( int i = point; i < Sind; i++)
    VarBin[i]= par2->VarBin[i];
}

// ============================== Differetiation ================================

void cIndivBinVec :: Differentiation(eDifType type, double f, cOptSolution *lastbest, cOptSolution *p1, cOptSolution *p2, cOptSolution *p3)
{
  cIndivBinVec *par1 = dynamic_cast<cIndivBinVec*>(p1);
  cIndivBinVec *par2 = dynamic_cast<cIndivBinVec*>(p2);
  cIndivBinVec *par3 = dynamic_cast<cIndivBinVec*>(p3);

  if (!par1 || !par2 || !par3)
    Utl::Exit("Failed dynamic_cast on cIndivBinVec::Crossover");

  switch (type)
  {
    case Rand1:
      for (int i = 0; i < Sind; i++)
      {
       Var[i] = (par3->Var[i] + f*(par1->Var[i] - par2->Var[i]));
      }
      break;

    //case CLASSICAL:
    //  for (int j = 0; j < pntid; j++)
    //   Var[j] = par1->Var[j];
    //  for (int j = pntid; j < nvar; j++)
    //   Var[j] = par2->Var[j];
    //break;
  }
}

// ============================== Decoding ================================

void cIndivBinVec :: Decoding(void)
{
  // Transform binary to temporary integers.

  int nvar = Prob->GetNumVar( );
  int cont = Sind - 1;
  int temp[nvar];
  for(int i = nvar - 1; i >= 0; i--)
  {
    temp[i] = 0;
    for (int k = 0; k < Scromo[i]; k++)
    {
      temp[i] += VarBin[cont]*pow(2.0, k);
      cont--;
    }
  }

  // Map the variables.

  int low,upp;
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp) ;
    Var[i] = round (low + (upp-low)/(pow(2.0, Scromo[i])-1)*temp[i]);
  }
}


// ================================ Send ===================================

void cIndivBinVec :: Send(void)
{
#ifdef _MPI_
  // Get necessary data.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );
  int nconstr = Prob->GetNumConstr( );

  // Store data in contiguous memory arrays.

  int varsend[Sind];
  double objsend;
  double constrsend[nconstr];

  for (int i = 0; i < Sind; i++)
    varsend[i] = VarBin[i];

  objsend = ObjFuncVal;

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, Sind, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 3);
    }
#endif

}

// =============================== Receive =================================

void cIndivBinVec :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nconstr = Prob->GetNumConstr( );

  // Prepare contiguous memory arrays.

  int varrecv[Sind];
  double objrecv;
  double constrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, Sind, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 3);

  // Store received values.

  for (int i = 0; i < Sind; i++)
    VarBin[i] = varrecv[i];

  ObjFuncVal = objrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// ======================================================= End of file =====
