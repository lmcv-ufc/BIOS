// -------------------------------------------------------------------------
// food.cpp - implementation of the Artificial Bees Colony food class.
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
// Created:      05-Nov-2014    Elias Saraiva Barroso
//
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

#include "food.h"
#include "group.h"
#include "problem.h"
#include "lam.h"
#include "utl.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFood =====================================

cFood :: cFood(cProblem *prob) : cOptSolution(prob)
{
  Prob = prob;

  Fobjs.Resize(prob->GetNumObj( ));
  Fobjs.Zero( );

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;

  int nc = Prob->GetNumConstr( );
  Constr.Resize(nc);
  Constr.Zero( );

  FailCounter = 0;
}

// ============================ ~cFood =====================================

cFood :: ~cFood(void)
{
}

// ============================== CreateFood ===============================

cFood* cFood :: CreateFood(eSolType type, cProblem* prob)
{
  switch(type)
  {
    case SOL_INT_VEC:
      return new cFoodIntVec(prob);
      cout << "Error: Integer Food not yet implemented!" << endl;
      exit(0);
    break;

    case SOL_INT_MAT:
      return new cFoodIntMat(prob);
      cout << "Error: Integer Matrix Food not yet implemented!" << endl;
      exit(0);
    break;

    case SOL_DBL_VEC:
      return new cFoodDblVec(prob);
    break;
  }

  return(0);
}

// -------------------------------------------------------------------------
// cFoodIntVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFoodIntVec ===============================

cFoodIntVec :: cFoodIntVec(cProblem *prob) : cFood(prob)
{
  OptSolType = SOL_INT_VEC;
  int nvar = Prob->GetNumVar( );
  Var = new int[nvar];
}

// ============================ ~cFoodIntVec ==============================

cFoodIntVec :: ~cFoodIntVec(void)
{
  delete [] Var;
}

// ================================ Init ===================================

void cFoodIntVec :: Init(void)
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

void cFoodIntVec :: Init(const sInpSol &ivar)
{
  int nvar = Prob->GetNumVar( );

  // Store input values.
  for(int i = 0; i < nvar; ++i)
    Var[i] = static_cast<int>(ivar.CodVar[i]);
}

// ============================== Evaluate =================================

void cFoodIntVec :: Evaluate( )
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cFoodIntVec :: Print(void)
{
  int nv = Prob->GetNumVar( );
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;

  for (int i = 0; i < nv; i++) cout << Var[i] << "  ";
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
  cout << "ObjectiveFunction = " << Fobjs[0] << endl;

  if (Prob->GetNumConstr( ) > 0)
  { cout << "------------------------------------------------------------" << endl;
    cout << "PenObjFunction    = " << PenObjFunc << endl;
  }
 
  cout << "------------------------------------------------------------" << endl;
  cout << "FitnessFunction   = " << FitFuncVal << endl << endl;
}

// =============================== Write ===================================

void cFoodIntVec :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  out << Fobjs[0] << endl;
  
  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
    out << PenObjFunc << endl;
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

void cFoodIntVec :: Copy(cOptSolution *a)
{
  cFoodIntVec *food = dynamic_cast<cFoodIntVec*>(a);

  if (!a)
    Utl::Exit("Failed dynamic_cast on cFoodIntVec::Copy");

  Fobjs[0] = food->Fobjs[0];
  PenObjFunc = food->PenObjFunc;
  FitFuncVal = food->FitFuncVal;

  Constr = food->Constr;

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++) 
    Var[i] = food->Var[i];
}

// ================================ CompVar ================================

bool cFoodIntVec :: CompVar(cOptSolution *sol)
{
  cFoodIntVec *food = dynamic_cast<cFoodIntVec*> (sol);

  if (food == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }
  
  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++) 
    if (Var[i] != food->Var[i])
     return false;

  return true;
}

// ============================== ABCmutate ===============================

void cFoodIntVec :: ABCmutate(cFoodSource &foods, int FoodIndex)
{
  int nvar = Prob->GetNumVar( );

  int NumFood = foods.GetSize();
  int i = Utl::RandInt(0, nvar);
  int j = 0;
  
  if (NumFood > 1)
  {
    do
      j = Utl::RandInt(0, NumFood-1);
    while (j == FoodIndex);
  }
  
  double r = Utl::RandDec( );

  cFoodIntVec *food = dynamic_cast<cFoodIntVec*>(foods[j]);

  if (!food)
    Utl::Exit("Failed dynamic_cast on cFoodIntVec::ABCmutate");

  r = (Var[i] - food->Var[i] )*(r - 0.5)*2;
  
  Var[i] = Var[i] + round(r);
}

// ================================ Send ===================================

void cFoodIntVec :: Send(void)
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
  double fcounter;
  double constrsend[nconstr];

  for (int i = 0; i < nvar; i++)
    varsend[i] = Var[i];

  objsend = Fobjs[0];
  fcounter = FailCounter; 

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, nvar, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(&fcounter, 1, MPI_DOUBLE, j, 3);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 4);
    }
#endif
}

// =============================== Receive =================================

void cFoodIntVec :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Prepare contiguous memory arrays.

  int varrecv[nvar];
  double objrecv;
  double fcrecv;
  double constrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, nvar, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(&fcrecv, 1, MPI_DOUBLE, deme, 3);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 4);

  // Store received values.

  for (int i = 0; i < nvar; i++)
    Var[i] = varrecv[i];

  Fobjs[0] = objrecv;
  FailCounter = fcrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// -------------------------------------------------------------------------
// cFoodIntMat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFoodIntMat ===============================

cFoodIntMat :: cFoodIntMat(cProblem *prob) : cFood(prob)
{
  OptSolType = SOL_INT_MAT;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  Var = new int*[n];
  for (int i = 0; i < n; i++) 
    Var[i] = new int[m];
}

// ============================ ~cFoodIntMat ===============================

cFoodIntMat :: ~cFoodIntMat(void)
{
  int n = Prob->VarNumRow( );
  for (int i = 0; i < n; i++) delete []Var[i];
  delete []Var;
}

// ================================ Init ===================================

void cFoodIntMat :: Init( )
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

void cFoodIntMat :: Init(const sInpSol &ivar)
{
  // Store input values.
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Var[i][j] = static_cast<int>(ivar.CodVar[i*m+j]);
}

// ================================ Copy ===================================

void cFoodIntMat :: Copy(cOptSolution *a)
{
  cFoodIntMat *food = dynamic_cast<cFoodIntMat*>(a);

  if (!a)
    Utl::Exit("Failed dynamic_cast on cFoodIntMat::Copy");

  Fobjs[0] = food->Fobjs[0];
  PenObjFunc = food->PenObjFunc;
  FitFuncVal = food->FitFuncVal;

  Constr = food->Constr;

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) 
      Var[i][j] = food->Var[i][j];
}

// ================================ CompVar ================================

bool cFoodIntMat :: CompVar(cOptSolution *sol)
{
  cFoodIntMat *food = dynamic_cast<cFoodIntMat*> (sol);

  if (food == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) 
      if (Var[i][j] != food->Var[i][j])
        return false;

  return true;
}

// ============================== Evaluate =================================

void cFoodIntMat :: Evaluate( )
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cFoodIntMat :: Print(void)
{
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) cout << Var[i][j] << "  ";
    cout << endl;
  }

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
  cout << "ObjectiveFunction = " << Fobjs[0] << endl;
  
  if (Prob->GetNumConstr( ) > 0)
  {   
    cout << "------------------------------------------------------------" << endl;
    cout << "PenObjFunction    = " << PenObjFunc << endl;
  }

  cout << "------------------------------------------------------------" << endl;
  cout << "FitnessFunction   = " << FitFuncVal << endl << endl;
}

// =============================== Write ===================================

void cFoodIntMat :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  out << Fobjs[0] << endl;
  
  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
    out << PenObjFunc << endl;
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

// ================================= Swap ==================================

void cFoodIntMat :: Swap(double swp_rate)
{
  // Required data. 

  cLaminated *lamprob = dynamic_cast<cLaminated*>(Prob);

  if (!lamprob)
    Utl::Exit("Failed dynamic_cast on cFoodIntMat::Swap");

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

void cFoodIntMat :: Add(double add_rate)
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

void cFoodIntMat :: Delete(double del_rate)
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

// ============================== ABCmutate ===============================

void cFoodIntMat :: ABCmutate(cFoodSource &foods, int FoodIndex)
{
  int ncol = Prob->VarNumCol();
  int nrow = Prob->VarNumRow();

  // Get only valid row index

  int low,upp;
  
  vector<int> rowInds;

  for(int i=0; i<nrow; i++)
  {
    Prob->GetBounds(i,&low,&upp);

    if (upp > 0)
      rowInds.push_back(i);
  }

  // Get random position on Matrix to mutate

  int i = rowInds[Utl::RandInt(0, rowInds.size()-1)]; 
  int j = Utl::RandInt(0,ncol-1);
  int k = 0;

  int NumFood = foods.GetSize();

  if (NumFood > 1)
  {
    do
      k = Utl::RandInt(0, NumFood-1);
    while (k == FoodIndex);
  }

  cFoodIntMat *food = dynamic_cast<cFoodIntMat*>(foods[k]);

  if (!food)
    Utl::Exit("Failed dynamic_cast on cFoodIntMat::ABCmutate");

  double r = Utl::RandDec( );

  r = (Var[i][j] - food->Var[i][j])*(r - 0.5)*2;
  
  Var[i][j] = Var[i][j] + round(r);
  
  // Check if food's position leaves problem bounds
    
  Prob->GetBounds(i, &low, &upp);

  if(Var[i][j] > upp)
    Var[i][j] = upp;
  else if (Var[i][j] < low)
    Var[i][j] = low;
}

// ================================ Send ===================================

void cFoodIntMat :: Send(void)
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
  double fcounter;
  double constrsend[nconstr];

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      varsend[m*i + j] = Var[i][j];

  objsend = Fobjs[0];
  fcounter = FailCounter;

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, n*m, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(&fcounter, 1, MPI_DOUBLE, j, 3);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 4);
    }
#endif
}

// =============================== Receive =================================

void cFoodIntMat :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nconstr = Prob->GetNumConstr( );

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  // Prepare contiguous memory arrays.

  int varrecv[n*m];
  double objrecv;
  double fcrecv;
  double constrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, n*m, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(&fcrecv, 1, MPI_DOUBLE, deme, 3);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 4);

  // Store received values.

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Var[i][j] = varrecv[m*i + j];

  Fobjs[0] = objrecv;
  FailCounter = fcrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// -------------------------------------------------------------------------
// cFoodDblVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFoodDblVec ===============================

cFoodDblVec :: cFoodDblVec(cProblem *prob) : cFood(prob)
{
  OptSolType = SOL_DBL_VEC;
  int nv = Prob->GetNumVar( );
  Var.Resize(nv);
  Var.Zero( );
}

// ============================ ~cFoodDblVec ===============================

cFoodDblVec :: ~cFoodDblVec(void)
{
}

// ================================ Init ===================================

void cFoodDblVec :: Init( )
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

void cFoodDblVec :: Init(const sInpSol &ivar)
{
  int nvar = Prob->GetNumVar( );

  // Store input values.
  for(int i = 0; i < nvar; ++i)
    Var[i] = ivar.CodVar[i];
}

// ============================== Evaluate =================================

void cFoodDblVec :: Evaluate( )
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(Var, Constr, Fobjs);
}

// =============================== Print ===================================

void cFoodDblVec :: Print(void)
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
  cout << "ObjectiveFunction = " << Fobjs[0] << endl;
  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "------------------------------------------------------------" << endl;
    cout << "PenObjFunction    = " << PenObjFunc << endl;
  }

  cout << "------------------------------------------------------------" << endl;
  cout << "FitnessFunction   = " << FitFuncVal << endl << endl;
}

// =============================== Write ===================================

void cFoodDblVec :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  out << Fobjs[0] << endl;
  
  out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
  out << PenObjFunc << endl;
  
  out << "\n%DESIGN.VARIABLES\n";
  for (int i = 0; i < Constr.Dim( ); i++)
    out << Var[i] << endl;
  out << endl;

  out << "\n%CONSTRAINT.VALUES\n";
  for (int i = 0; i < Constr.Dim( ); i++)
    out << Constr[i] << endl;
}

// ================================ Copy ===================================

void cFoodDblVec :: Copy(cOptSolution *a)
{
  cFoodDblVec *food = dynamic_cast<cFoodDblVec*>(a);

  if (!a)
    Utl::Exit("Failed dynamic_cast on cFoodDblVec::Copy");

  Fobjs[0] = food->Fobjs[0];
  PenObjFunc = food->PenObjFunc;
  FitFuncVal = food->FitFuncVal;

  Constr = food->Constr;

  int nv =  Prob->GetNumVar();
  for (int i = 0; i < nv; i++) 
    Var[i] = food -> Var[i];
}

// ================================ CompVar ================================

bool cFoodDblVec :: CompVar(cOptSolution *sol)
{
  cFoodDblVec *food = dynamic_cast<cFoodDblVec*> (sol);

  if (food == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }
  
  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++) 
    if (Var[i] != food->Var[i])
     return false;

  return true;
}

// ============================== ABCmutate ===============================

void cFoodDblVec :: ABCmutate(cFoodSource &foods, int FoodIndex)
{
  int nvar = Prob->GetNumVar( );
  int NumFood = foods.GetSize();

  int vid = Utl::RandInt(0, nvar-1);     // selected variable indice
  int fid = Utl::RandInt(0, NumFood-1);  // selected food indice

  if (NumFood > 1)
  {
    do
      fid = Utl::RandInt(0, NumFood-1);
    while (fid == FoodIndex);
  }

  cFoodDblVec *food = dynamic_cast<cFoodDblVec*>(foods[fid]);

  if (!food)
    Utl::Exit("Failed dynamic_cast on cFoodDblVec::ABCmutate");

  double r = Utl::RandDec( );
 
  Var[vid] = Var[vid] + (Var[vid] - food->Var[vid])*(r - 0.5)*2;
}

// ================================ Send ===================================

void cFoodDblVec :: Send(void)
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
  double fcounter;
  double constrsend[nconstr];

  for (int i = 0; i < nvar; i++)
    varsend[i] = Var[i];

  objsend = Fobjs[0];
  fcounter = FailCounter;

  for (int i = 0; i < nconstr; i++)
    constrsend[i] = Constr[i];

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, nvar, MPI_DOUBLE, j, 1);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(&fcounter, 1, MPI_DOUBLE, j, 3);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 4);
    }
#endif
}

// =============================== Receive =================================

void cFoodDblVec :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Prepare contiguous memory arrays.

  double varrecv[nvar];
  double objrecv;
  double fcrecv;
  double constrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, nvar, MPI_DOUBLE, deme, 1);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(&fcrecv, 1, MPI_DOUBLE, deme, 3);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 4);

  // Store received values.

  for (int i = 0; i < nvar; i++)
    Var[i] = varrecv[i];

  Fobjs[0] = objrecv;

  for (int i = 0; i < nconstr; i++)
    Constr[i] = constrrecv[i];

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
#endif
}

// ======================================================= End of file =====
