// -------------------------------------------------------------------------
// particle.cpp - implementation of the particle class.
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

#ifdef _MPI_
#include "mpi.h"
#endif

#include "particle.h"
#include "problem.h"
#include "lam.h"
#include "utl.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cParticle ===============================

cParticle :: cParticle(cProblem *prob) : cOptSolution(prob)
{
  BestObjFuncVal = 0.0;
  BestPenObjFunc = 0.0;
  BestFitFuncVal = 0.0;
  
  int nc = Prob->GetNumConstr( );
  BestConstr.Resize(nc);
  BestConstr.Zero( );
}

// ============================ ~cParticle ===============================

cParticle :: ~cParticle(void)
{
}

// ============================== UpdateBestPos ===========================

void cParticle :: UpdateBestPos(bool force)
{
  // Verify if history position is worse than current position.
  // Note that force flag can block this verification.

  if (force == false)
  {
    if (Prob->GetNumConstr( ) == 0)
    {
      if(Fobjs[0] >= BestObjFuncVal)
      {
        //cout << "nao att" << endl;
        return;
      }
    }
    else
      if(PenObjFunc >= BestPenObjFunc)
        return;
  }

  // Assing object function (or penalized one) and constraints.
  BestPenObjFunc = PenObjFunc;
  BestObjFuncVal = Fobjs[0];
  BestConstr = Constr;

  // Assing project variables. 

  UpdateBestVar( );
}

// ============================== CreataParticle ===========================

cParticle* cParticle :: CreateParticle(eSolType type, cProblem* prob)
{
  switch(type)
  {
    case SOL_INT_VEC:
      return new cPartIntVec(prob);
    break;

    case SOL_INT_MAT:
      return new cPartIntMat(prob);
    break;

    case SOL_DBL_VEC:
      return new cPartDblVec(prob);
    break;

    case SOL_BIN_VEC:
      cout << "Error: Binary Particle not yet implemented!" << endl;
      exit(0);
    break;
  }

  return(0);
}

// -------------------------------------------------------------------------
// cIndivIntVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPartIntVec ==============================

cPartIntVec :: cPartIntVec(cProblem *prob) : cParticle(prob)
{
  OptSolType = SOL_INT_VEC;
  int nvar = Prob->GetNumVar( );
  PosVar = new int[nvar];
  BestPosVar = new int[nvar];
  VelVar.Resize(nvar);
  VelVar.Zero( );
}

// ============================ ~cPartIntVec ==============================

cPartIntVec :: ~cPartIntVec(void)
{
  delete [] PosVar;
  delete [] BestPosVar;
}

// ================================ Init ===================================

void cPartIntVec :: Init( )
{
  // Initialize the position and velocity variables with random values.
  int low,upp;
  int nvar = Prob->GetNumVar( );

  // Initialize the variables with random values.
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    PosVar[i] = Utl::RandInt(low, upp);
    VelVar[i] = (upp - low) *  Utl::RandDouble(-1.0, 1.0);
  }
}

void cPartIntVec :: Init(const sInpSol &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Store input values.
  int nvar = Prob->GetNumVar( );
  for(int i = 0; i < nvar; ++i)
    PosVar[i] = static_cast<int>(ivar.CodVar[i]);
}

void cPartIntVec :: Init(const cVector &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Initialize variables with given values in [0,1] range.
  int low,upp,nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    PosVar[i] = static_cast<int>(low + (upp-low)*ivar[i]);
    PosVar[i] = (PosVar[i] > upp) ? upp : PosVar[i];
  }
}

// =============================== GetVar ==================================

void cPartIntVec :: GetVar(int *a)
{
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
    a[i] = PosVar[i];
}

// =============================== GetVar ==================================

void cPartIntVec :: GetBestVar(int *a)
{
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
    a[i] = BestPosVar[i];
}
// ============================== Evaluate =================================

void cPartIntVec :: Evaluate(void)
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(PosVar, Constr, Fobjs);
}


// =============================== Print ===================================

void cPartIntVec :: Print(void)
{
  int nv = Prob->GetNumVar( );

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;
  
  for (int i = 0; i < nv; i++) cout << BestPosVar[i] << "  ";
  cout << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (problem)" << endl;
  cout << "------------------------------------------------------------" << endl;
  Prob->PrintVar(BestPosVar);

  if (Prob->GetNumConstr( ) > 0)
  { 
    cout << "------------------------------------------------------------" << endl;
    cout << "Constraints" << endl;
    cout << "------------------------------------------------------------" << endl;
    BestConstr.Print( );
  }

  cout << "------------------------------------------------------------" << endl;
  cout << "ObjectiveFunction = " << BestObjFuncVal << endl;

  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "------------------------------------------------------------" << endl;  
    cout << "PenObjFunction    = " << BestPenObjFunc << endl << endl;
  }

  else
    cout << endl;
}

// =============================== Write ===================================

void cPartIntVec :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  out << BestObjFuncVal << endl;
  
  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
    out << BestPenObjFunc << endl;
  }

  out << "\n%OPTIMIZATION.VARIABLES\n";
  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++) out << BestPosVar[i] << "  ";
  out << endl;

  out << "\n%PROBLEM.VARIABLES\n";
  Prob->WriteVar(BestPosVar);

  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << setw(3) << i+1 << "  " << BestConstr[i] << endl;
  }    
}

// ================================ Copy ===================================

void cPartIntVec :: Copy(cOptSolution *a)
{
  cPartIntVec *part = (cPartIntVec*)a;

  Fobjs[0] = part->Fobjs[0];
  PenObjFunc = part->PenObjFunc;
  FitFuncVal = part->FitFuncVal;
  BestObjFuncVal = part->BestObjFuncVal;
  BestPenObjFunc = part->BestPenObjFunc;
  BestFitFuncVal = part->BestFitFuncVal;

  Constr = part->Constr;
  BestConstr = part->BestConstr;

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
  {
    PosVar[i] = part->PosVar[i];
    BestPosVar[i] = part->BestPosVar[i];
    VelVar[i] = part->VelVar[i];
  }
}

// ================================ CompVar ================================

bool cPartIntVec :: CompVar(cOptSolution *sol)
{
  cPartIntVec *part = dynamic_cast<cPartIntVec*> (sol);

  if (part == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
    if (BestPosVar[i] != part->BestPosVar[i])
      return false;

  return true;
}

// =============================== Mutate ==================================

void cPartIntVec :: Mutate(double mut_rate)
{
  int low,upp;
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
  {
    double r = Utl::RandDec( );
    if (r <= mut_rate)
    {
      Prob->GetBounds(i, &low, &upp);
      VelVar[i] = Utl::RandDouble(low, upp) * Utl::RandInt(-1,1);
    }
  }
}


// ============================== EvalVelocity =======================================

void cPartIntVec :: EvalVelocity(double w, int num, double c[],cParticle *PVec[])
{
  cPartIntVec **part_Int = new cPartIntVec* [num];
   
  for (int i=0; i < num; i++)
    part_Int[i] = (cPartIntVec*) PVec[i];

  for (int i = 0; i < Prob->GetNumVar( ); i++)
  {
    VelVar[i] *= w;
    
    for (int j = 0; j < num; j++)
      VelVar[i] += c[j]*(part_Int[j]->BestPosVar[i] - PosVar[i]);
  }

  delete [] part_Int;
}

// =============================== GetNormVar ==================================

void cPartIntVec ::  GetNormVar(cVector &xn)
{
  Prob->GetNormVar(BestPosVar,xn);
}

// ============================== UpdatePos ===========================

void cPartIntVec :: UpdatePos(void)
{
  int low,upp;
  int nvar = Prob->GetNumVar( );
  
  for (int i = 0; i < nvar; i++)
  {
    // Update particle position

    PosVar[i] += (int)VelVar[i];
    
    // Check if particle's position leaves problem bounds
    
    Prob->GetBounds(i, &low, &upp);

    if(PosVar[i] > upp)
    {
      PosVar[i] = upp;
      VelVar[i] *= -0.5;
    }
    else if (PosVar[i] < low)
    {
      PosVar[i] = low;
      VelVar[i] *= -0.5;
    }
  }
}

// ============================== UpdateBestVar ===========================

void cPartIntVec :: UpdateBestVar(void)
{
  for(int i = 0; i < Prob->GetNumVar( ); i++)
    BestPosVar[i] = PosVar[i];
}

// ================================ Send ===================================

void cPartIntVec :: Send(void)
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

  int bvarsend[nvar];
  double bobjsend;
  double bconstrsend[nconstr];
  
  for (int i = 0; i < nvar; i++)
  {
    varsend[i] = PosVar[i];
    bvarsend[i] = BestPosVar[i];
  }

  objsend = Fobjs[0];
  bobjsend = BestObjFuncVal;

  for (int i = 0; i < nconstr; i++)
  {
    constrsend[i] = Constr[i];
    bconstrsend[i] = BestConstr[i];
  }

  // Broadcast the particle to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, nvar, MPI_INT, j, 1);
      MPI::COMM_WORLD.Send(bvarsend, nvar, MPI_INT, j, 2);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 3);
      MPI::COMM_WORLD.Send(&bobjsend, 1, MPI_DOUBLE, j, 4);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 5);
      MPI::COMM_WORLD.Send(bconstrsend, nconstr, MPI_DOUBLE, j, 6);
    }
#endif
}

// =============================== Receive =================================

void cPartIntVec :: Receive(int deme)
{
#ifdef _MPI_
  // Get necessary data.

  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Prepare memory arrays.

  int varrecv[nvar];
  double objrecv;
  double constrrecv[nconstr];

  int bvarrecv[nvar];
  double bobjrecv;
  double bconstrrecv[nconstr];

  // Receive the individual.

  MPI::COMM_WORLD.Recv(varrecv, nvar, MPI_INT, deme, 1);
  MPI::COMM_WORLD.Recv(bvarrecv, nvar, MPI_INT, deme, 2);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 3);
  MPI::COMM_WORLD.Recv(&bobjrecv, 1, MPI_DOUBLE, deme, 4);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 5);
  MPI::COMM_WORLD.Recv(bconstrrecv, nconstr, MPI_DOUBLE, deme, 6);

  // Store received values.

  for (int i = 0; i < nvar; i++)
  {
    PosVar[i] = varrecv[i];
    BestPosVar[i] = bvarrecv[i];
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

// -------------------------------------------------------------------------
// cIndivIntMat class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cIndivIntMat ==============================

cPartIntMat :: cPartIntMat(cProblem *prob) : cParticle(prob)
{
  OptSolType = SOL_INT_MAT;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  PosVar = new int*[n];
  BestPosVar = new int*[n];
  for (int i = 0; i < n; i++)
  {
    PosVar[i] = new int[m];
    BestPosVar[i] = new int[m];
  }

  VelVar.Resize(n,m);
  VelVar.Zero( );
}

// ============================ ~cPartIntMat ==============================

cPartIntMat :: ~cPartIntMat(void)
{
  int n = Prob->VarNumRow( );
  for (int i = 0; i < n; i++)
  {
    delete []PosVar[i];
    delete []BestPosVar[i];
  }
  delete []PosVar;
  delete []BestPosVar;
}

// ================================ Init ===================================

void cPartIntMat :: Init( )
{
  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
      PosVar[i][j] = Utl::RandInt(low, upp);
      VelVar[i][j] = (upp - low) *  Utl::RandDouble(-1,1); 
    }
  }
}

void cPartIntMat :: Init(const sInpSol &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Store input values.
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      PosVar[i][j] = static_cast<int>(ivar.CodVar[i*m+j]);
}

void cPartIntMat :: Init(const cVector &ivar)
{
  // Initialize velocity with random values.
  Init( );

  // Initialize variables with given values in [0,1] range.

  int low,upp;
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    Prob->GetBounds(i, &low, &upp);
    for (int j = 0; j < m; j++)
    {
      PosVar[i][j] = static_cast<int>(low + (upp-low)*ivar[i*m+j]);
      PosVar[i][j] = (PosVar[i][j] > upp) ? upp : PosVar[i][j];
    }
  }
}

// ================================ Copy ===================================

void cPartIntMat :: Copy(cOptSolution *a)
{
  cPartIntMat *part = (cPartIntMat*)a;

  Fobjs[0] = part->Fobjs[0];
  PenObjFunc = part->PenObjFunc;
  FitFuncVal = part->FitFuncVal;
 
  Constr = part->Constr;

  BestObjFuncVal = part->BestObjFuncVal;
  BestPenObjFunc = part->BestPenObjFunc;
  BestFitFuncVal = part->BestFitFuncVal;

  BestConstr = part->BestConstr;

  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
    {
      PosVar[i][j] = part->PosVar[i][j];
      BestPosVar[i][j] = part->BestPosVar[i][j];
      VelVar[i][j] = part->VelVar[i][j];
    }
}

// ================================ CompVar ================================

bool cPartIntMat :: CompVar(cOptSolution *sol)
{
  cPartIntMat *part = dynamic_cast<cPartIntMat*> (sol);

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
      if (BestPosVar[i][j] != part->BestPosVar[i][j])
        return false;

  return true;
}

// ============================== Evaluate =================================

void cPartIntMat :: Evaluate(void)
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(PosVar, Constr, Fobjs);
}

// =============================== Print ===================================

void cPartIntMat :: Print(void)
{
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (algorithm)" << endl;
  cout << "------------------------------------------------------------" << endl;

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) cout << BestPosVar[i][j] << "  ";
    cout << endl;
  }

  cout << "------------------------------------------------------------" << endl;
  cout << "Variables (problem)" << endl;
  cout << "------------------------------------------------------------" << endl;
  Prob->PrintVar(BestPosVar);

  if (Prob->GetNumConstr( ) > 0)
  { 
    cout << "------------------------------------------------------------" << endl;
    cout << "Constraints" << endl;
    cout << "------------------------------------------------------------" << endl;
    BestConstr.Print( );
  }

  cout << "------------------------------------------------------------" << endl;
  cout << "ObjectiveFunction = " << BestObjFuncVal << endl;

  if (Prob->GetNumConstr( ) > 0)
  { 
    cout << "------------------------------------------------------------" << endl;
    cout << "PenObjFunction    = " << BestPenObjFunc << endl << endl;
  }

  else
    cout << endl;
}

// =============================== Write ===================================

void cPartIntMat :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n";
  out << BestObjFuncVal << endl;
  
  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n";
    out << BestPenObjFunc << endl;
  }
  
  out << "\n%OPTIMIZATION.VARIABLES\n";
  int n = Prob->VarNumRow( );
  int m = Prob->VarNumCol( );
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) out << BestPosVar[i][j] << "  ";
    out << endl;
  }

  out << "\n%PROBLEM.VARIABLES\n";
  Prob->WriteVar(BestPosVar);

  if (Prob->GetNumConstr( ) > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < Constr.Dim( ); i++)
      out << setw(3) << i+1 << "  " << BestConstr[i] << endl;
  }    
}

// =============================== Mutate ==================================

void cPartIntMat :: Mutate(double mut_rate)
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
        PosVar[i][j] = Utl::RandInt(low, upp);
    }
  }
}

// =============================== LamMutate ===============================

void cPartIntMat :: LamMutate(double *mut_rate)
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
      if (r <= mut_rate[i]) 
        PosVar[i][j] = Utl::RandInt(low, upp);
    }
  }
}

// ============================== EvalVelocity =======================================

void cPartIntMat :: EvalVelocity(double w, int num, double c[],cParticle *PVec[])
{
  cPartIntMat **part_Mat = new cPartIntMat* [num];
   
  for (int i=0; i < num; i++)
    part_Mat[i] = (cPartIntMat*) PVec[i];

  for (int i = 0; i < Prob->VarNumRow(); i++)
  {
    for (int j = 0; j < Prob->VarNumCol(); j++)
    { 
      VelVar[i][j] *= w;

      for (int k = 0; k < num; k++)
        VelVar[i][j] += c[k]*(part_Mat[k]->BestPosVar[i][j] - PosVar[i][j]);
    }
  }

  delete [] part_Mat;
}

// ============================== UpdatePos ===========================

void cPartIntMat :: UpdatePos(void)
{
  int low,upp;
  
  for (int i = 0; i < Prob->VarNumRow(); i++)
  {
    Prob->GetBounds(i, &low, &upp);

    for (int j = 0; j < Prob->VarNumCol(); j++)
    {
      // Update particle position

      PosVar[i][j] += round(VelVar[i][j]);
          
      // Check if particle's position leaves problem bounds
      
      if(PosVar[i][j] > upp)
      {
        PosVar[i][j] = upp;
        VelVar[i][j] *= -0.5;
      }
      else if (PosVar[i][j] < low)
      {
        PosVar[i][j] = low;
        VelVar[i][j] *= -0.5;
      }
    }
  }

}

// ============================== UpdateBestVar ===========================

void cPartIntMat :: UpdateBestVar(void)
{
  for (int i = 0; i < Prob->VarNumRow(); i++)
    for (int j = 0; j < Prob->VarNumCol(); j++) 
      BestPosVar[i][j] = PosVar[i][j];   
}

// ================================= Swap ==================================

void cPartIntMat :: Swap(double swp_rate)
{
  // Required data. 

  cLaminated *lamprob = dynamic_cast<cLaminated*>(Prob->GetHFP( ));

  if (!lamprob)
    Utl::Exit("Failed dynamic_cast on cPartIntMat::Swap");

  int n = lamprob->VarNumRow( );
  int m = lamprob->VarNumCol( );
  int *TempVar = new int[n];

  // Check if swap is possible.
  
  int nonnull = 0;
  for (int j = 0; j < m; j++)
    if (lamprob->GetThk(PosVar[0][j]) != 0.0)
      nonnull++;

  // Perform the swap operation.

  if (nonnull >= 2)
    for (int j = 0; j < m; j++)
      if (lamprob->GetThk(PosVar[0][j]) != 0.0)
      {
        double r = Utl :: RandDec( );
        if (r <= swp_rate)
        {
          for (int k = 0; k < n; k++)
            TempVar[k] = PosVar[k][j];

          while(1)
          {
            int k = Utl :: RandInt(0, m-1);
            if (k != j && lamprob->GetThk(PosVar[0][j]) != 0.0)
            {
              for (int l = 0; l < n; l++)
              {
                PosVar[l][j] = PosVar[l][k];
                PosVar[l][k] = TempVar[l];
              }
              break;
            }
          }
        }
      }
  delete[] TempVar;
}

// ================================= Add ===================================

void cPartIntMat :: Add(double add_rate)
{
  // Required data.

  int m = Prob->VarNumCol( );

  // Check if addition is possible.

  int null = 0;
  for (int j = 0; j < m; j++)
    if (PosVar[0][j] == 0)
      null++;

  if (null == 0)
    return;

  // Perform the addition operation.

  int low, upp;
  Prob->GetBounds(0, &low, &upp);

  for (int k = 0; k < m; k++)
    if (PosVar[0][k] == 0)
    {
      double r = Utl :: RandDec( );
      if (r <= add_rate) PosVar[0][k]++;
    }
}

// ================================ Delete =================================

void cPartIntMat :: Delete(double del_rate)
{
  // Required data.

  int m = Prob->VarNumCol( );

  // Check if deletion is possible.

  int nonnull = 0;
  for (int j = 0; j < m; j++)
    if (PosVar[0][j] != 0)
      nonnull++;

  if (nonnull == 0)
    return;

  // Perform the deletion operation.

  for (int k = 0; k < m; k++)
    if (PosVar[0][k] != 0)
    {
      double r = Utl :: RandDec( );
      if (r <= del_rate) PosVar[0][k]--;
    }
}
// =============================== GetNormVar ==================================

void cPartIntMat ::  GetNormVar(cVector &xn)
{
  Prob->GetHFP( )->GetNormVar(PosVar,xn);
}

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

// -------------------------------------------------------------------------
// cPartDblVec class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPartDblVec ===============================

cPartDblVec :: cPartDblVec(cProblem *prob) : cParticle(prob) 
{
  OptSolType = SOL_DBL_VEC;

  int nvar = Prob->GetNumVar( );

  PosVar.Resize(nvar);
  PosVar.Zero( );
  BestPosVar.Resize(nvar);
  BestPosVar.Zero( );
  VelVar.Resize(nvar);
  VelVar.Zero( );
}

// ============================ ~cPartDblVec ===============================

cPartDblVec :: ~cPartDblVec(void)
{
}

// ================================ Init ===================================

void cPartDblVec :: Init( )
{
  // Initialize the position and velocity variables with random values.
  double *low,*upp;
  int nvar = Prob->GetNumVar( );

  low = new double[nvar];
  upp = new double[nvar];

  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

//  for (int i=0; i< nvar; i++) cout << "low: " << low[i] << " upp: " << upp[i] << endl;

  for (int i = 0; i < nvar; i++)
  {
    PosVar[i] = Utl::RandDouble(low[i], upp[i]); 
    VelVar[i] = (upp[i] - low[i]) * Utl::RandDouble(-1,1); 
  }
 
  delete[] low;
  delete[] upp;
}

void cPartDblVec :: Init(const sInpSol &ivar)
{
  // Initialize velocity variables with random values.
  Init( );

  // Store input values.
  int nvar = Prob->GetNumVar( );
  for (int i = 0; i < nvar; i++)
    PosVar[i] = ivar.CodVar[i];
}

void cPartDblVec :: Init(const cVector &ivar)
{
  int nvar = Prob->GetNumVar( );
  double *low,*upp;
  low = new double[nvar];
  upp = new double[nvar];

  // Decode variables from [0,1] to problem bounds.
  for (int i=0; i< nvar; i++)
    Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
    PosVar[i] = low[i] + (upp[i]-low[i])*ivar[i];

  delete[] low;
  delete[] upp;
}

// ============================== Evaluate =================================

void cPartDblVec :: Evaluate( )
{
  // Evaluate the objective function and constraints.
  
  Prob->Evaluate(PosVar, Constr, Fobjs);
}

// =============================== Print ===================================

void cPartDblVec :: Print( )
{
  cout << "Best particle position" << endl; 
  cout << "NumVar = " << BestPosVar.Dim( ) << endl;
  cout << "Variables" << endl;
  BestPosVar.Print( );
 
  if (Prob->GetNumConstr( ) > 0)
  {
    cout << "Constraints" << endl;
    BestConstr.Print( );
  }
 
  cout << "ObjectiveFunction = " << BestObjFuncVal << endl;
  if (Prob->GetNumConstr( ) > 0)
    cout << "PenObjFunction    = " << BestPenObjFunc << endl << endl;
  else
    cout << endl;
}

// =============================== Write ===================================

void cPartDblVec :: Write(ostream &out)
{
  out << "\n%RESULT.OBJECTIVE.FUNCTION\n" << BestObjFuncVal << endl;

  if (Prob->GetNumConstr() > 0)
    out << "\n%RESULT.PENALIZED.OBJECTIVE.FUNCTION\n" << BestPenObjFunc << endl;
  
  out << "\n%DESIGN.VARIABLES\n";
  for (int i = 0; i < BestPosVar.Dim( ); i++)
    out << BestPosVar[i] << endl;
  out << endl;

  if (Prob->GetNumConstr() > 0)
  {
    out << "\n%CONSTRAINT.VALUES\n";
    for (int i = 0; i < BestConstr.Dim( ); i++)
      out << BestConstr[i] << endl;
  }
}

// ================================ Copy ===================================

void cPartDblVec :: Copy(cOptSolution *a)
{
  cPartDblVec *part = (cPartDblVec*) a;

  Fobjs[0]     = part->Fobjs[0];
  BestObjFuncVal = part->BestObjFuncVal;
  PenObjFunc     = part->PenObjFunc;
  BestPenObjFunc = part->BestPenObjFunc;
  FitFuncVal     = part->FitFuncVal;
  BestFitFuncVal = part->BestFitFuncVal;

  Constr     = part->Constr;
  BestConstr = part->BestConstr;

  int nv =  Prob->GetNumVar();

  for (int i = 0; i < nv; i++)
  {
    PosVar[i]  = part->PosVar[i];
    BestPosVar[i] = part->BestPosVar[i];
    VelVar[i]  = part->VelVar[i];
  }
}

// ================================ CompVar ================================

bool cPartDblVec :: CompVar(cOptSolution *sol)
{
  cPartDblVec *part = dynamic_cast<cPartDblVec*> (sol);

  if (part == 0)
  {
    cout << "Warning: CompVar(cOptSolution*) is comparing objects of different";
    cout << " types!" << endl;
    return false;
  }

  int nv = Prob->GetNumVar( );
  for (int i = 0; i < nv; i++)
    if (BestPosVar[i] != part->BestPosVar[i])
      return false;

  return true;
}

// =============================== GetNormVar ==================================

void cPartDblVec ::  GetNormVar(cVector &xn)
{
  Prob->GetNormVar(BestPosVar,xn);
}

// =============================== Mutate ==================================

void cPartDblVec :: Mutate(double mut_rate)
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
      VelVar[i] = Utl::RandDouble(low[i], upp[i]) * Utl::RandInt(-1,1);
  }

  delete[] low;
  delete[] upp;
}

// ============================== EvalVelocity =======================================

void cPartDblVec :: EvalVelocity(double w, int num, double c[],cParticle *PVec[])
{
  cPartDblVec **part_Dbl = new cPartDblVec* [num];
   
  for (int i=0; i < num; i++)
    part_Dbl[i] = (cPartDblVec*) PVec[i];

  for (int i = 0; i < Prob->GetNumVar( ); i++)
  {
    VelVar[i] *= w;
    
    for (int j = 0; j < num; j++)
      VelVar[i] += c[j]*(part_Dbl[j]->BestPosVar[i] - PosVar[i]);
  }

  delete [] part_Dbl;
}

// ============================== UpdatePos ===========================

void cPartDblVec :: UpdatePos(void)
{
  double *low,*upp;
  int nvar = Prob->GetNumVar( );
 
  low = new double[nvar];
  upp = new double[nvar];

  Prob->GetDblBounds(low,upp);

  for (int i = 0; i < nvar; i++)
  {
    // Update particle position 
    
    PosVar[i] += VelVar[i];

    // Check if particle's position leaves problem bounds

    if(PosVar[i] > upp[i])
    {
      PosVar[i] = upp[i];
      VelVar[i] *= -0.5;
    }
    else if (PosVar[i] < low[i])
    {
      PosVar[i] = low[i];
      VelVar[i] *= -0.5;
    }
  }

}

// ============================== UpdateBestVar ===========================

void cPartDblVec :: UpdateBestVar(void)
{
  BestPosVar = PosVar;
}

// ================================ Send ===================================

void cPartDblVec :: Send(void)
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

  double bvarsend[nvar];
  double bobjsend;
  double bconstrsend[nconstr];

  for (int i = 0; i < nvar; i++)
  {
    varsend[i]  = PosVar[i];
    bvarsend[i] = BestPosVar[i];
  }

  objsend  = Fobjs[0];
  bobjsend = BestObjFuncVal;

  for (int i = 0; i < nconstr; i++)
  {
    constrsend[i]  = Constr[i];
    bconstrsend[i] = BestConstr[i];
  }

  // Broadcast the individual to every other deme.

  for (int j = 0; j < size; j++)
    if (rank != j)
    {
      MPI::COMM_WORLD.Send(varsend, nvar, MPI_DOUBLE, j, 1);
      MPI::COMM_WORLD.Send(bvarsend, nvar, MPI_DOUBLE, j, 2);
      MPI::COMM_WORLD.Send(&objsend, 1, MPI_DOUBLE, j, 3);
      MPI::COMM_WORLD.Send(&bobjsend, 1, MPI_DOUBLE, j, 4);
      MPI::COMM_WORLD.Send(constrsend, nconstr, MPI_DOUBLE, j, 5);
      MPI::COMM_WORLD.Send(bconstrsend, nconstr, MPI_DOUBLE, j, 6);
    }
#endif
}

// =============================== Receive =================================

void cPartDblVec :: Receive(int deme)
{
  #ifdef _MPI_
  // Get necessary data.

  int nvar = Prob->GetNumVar( );
  int nconstr = Prob->GetNumConstr( );

  // Prepare contiguous memory arrays.

  double varrecv[nvar];
  double objrecv;
  double constrrecv[nconstr];
  
  double bvarrecv[nvar];
  double bobjrecv;
  double bconstrrecv[nconstr];

  // Receive the particle data.

  MPI::COMM_WORLD.Recv(varrecv, nvar, MPI_DOUBLE, deme, 1);
  MPI::COMM_WORLD.Recv(bvarrecv, nvar, MPI_DOUBLE, deme, 2);
  MPI::COMM_WORLD.Recv(&objrecv, 1, MPI_DOUBLE, deme, 3);
  MPI::COMM_WORLD.Recv(&bobjrecv, 1, MPI_DOUBLE, deme, 4);
  MPI::COMM_WORLD.Recv(constrrecv, nconstr, MPI_DOUBLE, deme, 5);
  MPI::COMM_WORLD.Recv(bconstrrecv, nconstr, MPI_DOUBLE, deme, 6);

  // Store received values.

  for (int i = 0; i < nvar; i++)
  {
    PosVar[i]  = varrecv[i];
    BestPosVar[i] = bvarrecv[i];
  }

  Fobjs[0]     = objrecv;
  BestObjFuncVal = bobjrecv;

  for (int i = 0; i < nconstr; i++)
  {
    Constr[i] = constrrecv[i];
    BestConstr[i] = bconstrrecv[i];
  }

  PenObjFunc = 0.0;
  FitFuncVal = 0.0;
  BestPenObjFunc = 0.0;
  BestFitFuncVal = 0.0;
#endif
}

// ======================================================= End of file =====
