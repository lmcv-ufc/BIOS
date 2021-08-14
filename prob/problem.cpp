// -------------------------------------------------------------------------
// problem.cpp - implementation of the problem class.
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
// Created:      21-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     05-Mar-2013    Evandro Parente Junior
//               File name changed to problem.cpp. Created different 
//               Evaluate and GetBounds methods.
//   
// Modified:     20-Nov-2015    Elias Saraiva Barroso
//               Creation of cProblemFactory class to handle problem creating
//               process decentralized. 
//
// Modified:     03-Mar-2019    Marina Alves Maia
//               Evaluate method was changed to handle the multiobjective
//               problems by using a vector of variable size to store the
//               objective function(s) values.
// -------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <string>

using namespace std;

#include "problem.h"
#include "metaopt.h"
#include "benchmark.h"
#include "lam.h"
#include "lamplt.h"
#include "fgm.h"
#include "fgmplt.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
bool      cProblem :: MinObjFlag  = false; 
double    cProblem :: MinObjFunc  = 0;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadMinObjFunc =============================

void cProblem :: ReadMinObjFunc(void)
{
  if (!(in >> MinObjFunc))
    {
      cout << "Error in the input of the number of optimizations." << endl;
      exit(0);
    }

  MinObjFlag = true;
}

// ============================== cProblem =================================

cProblem :: cProblem(void)
{
}

// ============================== ~cProblem ================================

cProblem :: ~cProblem(void)
{
}

// =============================== GetBounds ===============================

void cProblem :: GetBounds(int i, int *low, int *upp)
{
  cout << "Error: GetBounds not defined\n";
  exit(0);
}

// =============================== GetBounds ===============================

void cProblem :: GetVarBounds(int i, double &low, double &upp)
{
  cout << "Error: GetVarBounds not defined\n";
  exit(0);
}

// ============================= GetDblBounds ==============================

void cProblem :: GetDblBounds(double *low, double *upp)
{
    cout << "Error: GetDblBounds not defined\n";
      exit(0);
}

// ============================= GetNormVar ================================

void cProblem :: GetNormVar(cVector &x, cVector &nx)
{
  double *low,*upp;
  low = new double[NumVar];
  upp = new double[NumVar];

  // Encode variables from problem bounds to [0,1].
  GetDblBounds(low,upp);
  for (int i = 0; i < NumVar; i++)
    nx[i] = (x[i] - low[i]) / (upp[i]-low[i]);

  delete[] low;
  delete[] upp; 
}

// ============================= GetNormVar ================================

void cProblem :: GetNormVar(int *x, cVector &nx)
{
  // Encode variables from problem bounds to [0,1].
  int low, upp;
  for (int i = 0; i < NumVar; i++)
  {
    GetBounds(i,&low,&upp);
    nx[i] = (double (x[i] - low)) / (upp-low);
  }
}

// ============================= GetNormVar ================================

void cProblem :: GetNormVar(int **x, cVector &nx)
{
  // Encode variables from problem bounds to [0,1].

  int low, upp;
  int n = VarNumRow( );
  int m = VarNumCol( );
  int eff = VarNumEff( );
  int c = 0;

  nx.Resize(eff);

  for (int i = 0; i < n; i++)
  {
    GetBounds(i, &low, &upp);
    if (low == upp)
      continue;
    else
    {
      for (int j = 0; j < m; j++)
      {
        nx[c*m+j] = (double (x[i][j] - low)) / (upp-low);
      }
      c += 1;
    }
  }
}

// =============================== Evaluate ================================

void cProblem :: Evaluate(int *var, cVector &constr, cVector &fobjs)
{
  cout << "Error: Evaluate not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

void cProblem :: Evaluate(int **var, cVector &constr, cVector &fobjs)
{
  cout << "Error: Evaluate not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

void cProblem :: Evaluate(cVector &var, cVector &constr, cVector &fobjs)
{
  cout << "Error: Evaluate not defined\n";
  exit(0);
}

// ========================= EvalExactConstraint ===========================

void cProblem :: EvalExactConstraint(int index, cVector& x, double &c)
{
    cout << "Error: Evaluation of single constraints not defined\n";
    exit(0);
}

// ========================= EvalExactFobj ===========================

void cProblem :: EvalExactFobj(cVector& x, double &fobj)
{
    cout << "Error: Evaluation of exact fobj not defined\n";
    exit(0);
}

// ========================= EvalExactConstraint ===========================

void cProblem :: EvalExactConstraint(int index, int *algvar, double &c)
{
    cout << "Error: Evaluation of single constraints not defined\n";
    exit(0);
}

// ========================= EvalExactFobj ===========================

void cProblem :: EvalExactFobj(int *algvar, double &fobj)
{
    cout << "Error: Evaluation of exact fobj not defined\n";
    exit(0);
}

// ========================= EvalExactConstraint ===========================

void cProblem :: EvalExactConstraint(int index, int **algvar, double &c)
{
    cout << "Error: Evaluation of single constraints not defined\n";
    exit(0);
}

// ========================= EvalExactFobj ===========================

void cProblem :: EvalExactFobj(int **algvar, double &fobj)
{
    cout << "Error: Evaluation of exact fobj not defined\n";
    exit(0);
}

// ============================= GetApproxObj ==============================

void cProblem :: GetApproxConstr(bool *cflag)
{
  for (int i = 0; i < NumConstr; i++) cflag[i] = 1;
}

// ============================= GetApproxObj ==============================

void cProblem :: GetApproxObj(bool *objflag)
{
  for (int i = 0; i < NumObj; i++) objflag[i] = 1;
}

// -------------------------------------------------------------------------
// Class cProblemFactory:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Local methods:
//

// =============================== GetFactory =============================

cProblemFactory& cProblemFactory :: GetFactory(void)
{
  static cProblemFactory factory;
  return factory;
}

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== Register ================================

bool cProblemFactory :: Register(string label, ProbCreatorFunc func, string fext){
  map<string, ProbCreatorFunc> *Map = &(GetFactory( ).CreationalMapFunc);
  map<string, string> *MapExt = &(GetFactory( ).FileExtMap);

  if (Map->find(label) == Map->end( ))
  {
    pair<string, ProbCreatorFunc> TagFunc(label,func);
    Map->insert(TagFunc);
    if (fext != "")
    {
      pair<string,string> ProbExt(label,fext);
      MapExt->insert(ProbExt);
    }
  }
  else
    Utl::Exit("Label "+label+" already added in the problem factory.");

  return true;
}

// =============================== CreateProblem ===========================

cProblem* cProblemFactory :: CreateProblem(string label, bool IsCont)
{
  map<string, ProbCreatorFunc> *Map = &(GetFactory( ).CreationalMapFunc);
  map<string, string> *MapExt = &(GetFactory( ).FileExtMap);

  if (Map->find(label) != Map->end( ))
  {
    // Create concrete problem.
    cProblem *ptr = (*Map)[label](IsCont);
    ptr->Label = label;

    return ptr;
  }
  else
    Utl::Exit("This label: "+label+" dont have any binded problem.");

  return 0;
}

// =============================== GetProbFileExt ==========================

string cProblemFactory :: GetProbFileExt(cProblem *prob)
{
  map<string, string> *MapExt = &(GetFactory( ).FileExtMap);

  if (MapExt->find(prob->Label) != MapExt->end( ))
    return (*MapExt)[prob->Label];
  else
    return string("");

  return 0;
}

// ======================================================= End of file =====
