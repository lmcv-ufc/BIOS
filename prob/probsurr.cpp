// -------------------------------------------------------------------------
// probsurr.cpp - implementation of the cProbSurr class.
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
// Created:      27-Oct-2019    Leonardo Gonçalves Ribeiro
//
//
// -------------------------------------------------------------------------

#include <cmath>
#include <math.h>
#include <iostream>

#include "probsurr.h"
#include "sao.h"
#include "gbldef.h"
#include "gblvar.h"
#include "input.h"
#include "utl.h"
#include "vec.h"
#include "mat.h"
#include "matvec.h"
#include "sysmat.h"
#include "rbf.h"
#include "krg.h"
#include "group.h"
#include "optsolution.h"

#include <vector>

using namespace std;

// -------------------------------------------------------------------------
// Class cProbSurr:
// -------------------------------------------------------------------------

cProbSurr :: cProbSurr(cSURR *smod, cProblem *hf,eProbSurrState s)
{
  HighFidelity = hf;
  SurMod       = smod;
  state        = s;

  NumVar    = HighFidelity->VarNumEff( );
  NumConstr = HighFidelity->GetNumConstr( );
  NumObj    = HighFidelity->GetNumObj( );
}

cProbSurr :: cProbSurr(cSURR *smod, cProblem *hf, sProbAppOut *appout,double cb,eProbSurrState s)
{
  HighFidelity = hf;
  SurMod       = smod;
  state        = s;

  AppOut       = appout;
  currbest     = cb;

  NumVar    = HighFidelity->GetNumVar( );
  NumConstr = HighFidelity->GetNumConstr( );
  NumObj    = HighFidelity->GetNumObj( );

  NumVarEff = HighFidelity->VarNumEff( );
  NumRow    = HighFidelity->VarNumRow( );
  NumCol    = HighFidelity->VarNumCol( );
}

// ============================= GetBounds =================================

void cProbSurr :: GetBounds(int i, int *low, int *upp)
{
  HighFidelity->GetBounds(i,low,upp);
}

// ============================= GetVarBounds ==============================

void cProbSurr :: GetVarBounds(int i, double &low, double &upp)
{
  HighFidelity->GetVarBounds(i,low,upp);
}

// ========================== GetBoundsDouble ==============================

void cProbSurr :: GetDblBounds(double *low,double *upp)
{
  HighFidelity->GetDblBounds(low,upp);
}

// =========================CONTINUOUS VARIABLE===========================

// ============================= Evaluate ================================

void cProbSurr :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  switch (state)
  {
    case EVALUATE_SURROGATE:
      EvaluateSurr(x, c, fobjs);
    break;

    case EVALUATE_EXPECTED_IMPROVEMENT:
      EvaluateEI(x, c, fobjs);
      c.Zero( );
    break;

    case EVALUATE_LOWER_CONFIDENCE_BOUND:
      EvaluateLCB(x, c, fobjs);
    break;

    case EVALUATE_PROBABILITY_IMPROVEMENT:
      EvaluatePoI(x, c, fobjs);
      c.Zero( );
    break;

  case EVALUATE_ESTIMATED_QUADRATIC_ERROR:
      EvaluateS2(x, c, fobjs);
    c.Zero( );
  break;
  }
}

// ============================= EvaluateSurr =============================

void cProbSurr :: EvaluateSurr(cVector &x, cVector &c, cVector &fobjs)
{
  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(x,xn);

  // Get surrogate responses.
  cVector y(AppOut->GetNumAppOut( ));
  SurMod->Evaluate(xn,y);


  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    fobjs[0] = y[0];
  else
    HighFidelity->EvalExactFobj(x, fobjs[0]);

  // Constraints evaluation.

  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    c[AppOut->GetAppConstrID(i)] = y[AppOut->GetAppConstrOutID(i)];

  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),x,
                                      c[AppOut->GetExactConstrID(i)]);
}

// ============================= EvaluateSurr =============================

void cProbSurr :: EvaluateLCB(cVector &x, cVector &c, cVector &fobjs)
{
  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(x,xn);

  // Get surrogate responses.
  cVector y(AppOut->GetNumAppOut( ));
  y[0] = SurMod->EvalLCB(xn);


  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    fobjs[0] = y[0];
  else
    HighFidelity->EvalExactFobj(x, fobjs[0]);

  // Constraints evaluation.

  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    c[AppOut->GetAppConstrID(i)] = y[AppOut->GetAppConstrOutID(i)];

  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),x,
                                      c[AppOut->GetExactConstrID(i)]);
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluateEI(cVector &x, cVector &c, cVector &fobjs)
{
  // Auxiliary variables.
  double ei;                  // Expected improvement.
  double hof;                 // High-fidelity exact objective function.
  double ci;                  // High-fidelity exact constraint.

  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(x,xn);

  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    ei = SurMod->EvalExpImp(xn,currbest);
  else
  {
    HighFidelity->EvalExactFobj(x, hof);
    ei = (hof < currbest) ? abs(currbest - hof) : 0.0;
  }

  // Apply problem probability of feasibility (approximated constraints).  
  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    ei *= GetConstraintFactor(xn, AppOut->GetAppConstrOutID(i)); //SurMod->EvalProbFeas(xn,AppOut->GetAppConstrOutID(i)); LEO

  // Zero EI for violated extact constraints.
  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
  {
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),x,ci);
    if (ci > 0.0) ei = 0.0;
  }

  ei += 1e-10;    // Avoid -infinity result. 
  ei = log10(ei);
  fobjs[0] = -ei; 
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluatePoI(cVector &x, cVector &c, cVector &fobjs)
{
  // Auxiliary variables.
  double poi;                 // Probability of Improvement.
  double hof;                 // High-fidelity exact objective function.
  double ci;                  // High-fidelity exact constraint.

  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(x,xn);

  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    poi = SurMod->EvalProbImp(xn,currbest);
  else
  {
    HighFidelity->EvalExactFobj(x, hof);
    poi = (hof < currbest) ? 1.0 : 0.0;
  }

  // Apply problem probability of feasibility (approximated constraints).
  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    poi *= GetConstraintFactor(xn, AppOut->GetAppConstrOutID(i)); //SurMod->EvalProbFeas(xn,AppOut->GetAppConstrOutID(i)); LEO

  // Zero EI for violated extact constraints.
  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
  {
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),x,ci);
    if (ci > 0.0) poi = 0.0;
  }

  fobjs[0] = -poi;
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluateS2(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    /*
    double tolviol;
    if (gen < 40){
        tolviol = 0.1 - 0.1/(40 - gen);
    }
    else{
        tolviol = 0;
    }*/

  // Decodification of problem variables.
  cVector xn = x;
  vector<cVector> sy;
  int no;
  int ns;
  SurMod -> GetNumOut(no);
  SurMod -> GetSampleY(sy);
  SurMod -> GetNumSample(ns);
  cVector xb(NumVar);
  cVector yb(no);
  cVector y(no);

  c[0] = 0;
  SurMod ->GetSampData( ).Normalize(xn);

  // Objective function evaluation.

  y[0] = SurMod -> SSqrSur(xn, 0);

  /*double pf = 0;
  for (int i = 0; i < NumConstr; i++){
    pf = SurMod -> EvalConstraintPF(xn, i + 1, 0);
    y[0] = pf*y[0];
  }*/

  y[0] = log10(y[0]);

  double fobj = -y[0];
  fobjs[0] = fobj;
}

// =======================INTEGER-VECTOR VARIABLE=========================

// ============================= Evaluate ================================

void cProbSurr :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  switch (state)
  {
    case EVALUATE_SURROGATE:
      EvaluateSurr(algvar, c, fobjs);
    break;

    case EVALUATE_EXPECTED_IMPROVEMENT:
      EvaluateEI(algvar, c, fobjs);
      c.Zero( );
    break;

    case EVALUATE_LOWER_CONFIDENCE_BOUND:
      EvaluateLCB(algvar, c, fobjs);
    break;

    case EVALUATE_PROBABILITY_IMPROVEMENT:
      EvaluatePoI(algvar, c, fobjs);
      c.Zero( );
    break;

  case EVALUATE_ESTIMATED_QUADRATIC_ERROR:
    EvaluateS2(algvar, c, fobjs);
    c.Zero( );
  break;
  }
}

// ============================= EvaluateSurr =============================

void cProbSurr :: EvaluateSurr(int *algvar, cVector &c, cVector &fobjs)
{
  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(algvar,xn);

  // Get surrogate responses.
  cVector y(AppOut->GetNumAppOut( ));
  SurMod->Evaluate(xn,y);


  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    fobjs[0] = y[0];
  else
    HighFidelity->EvalExactFobj(algvar, fobjs[0]);

  // Constraints evaluation.

  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    c[AppOut->GetAppConstrID(i)] = y[AppOut->GetAppConstrOutID(i)];

  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,
                                      c[AppOut->GetExactConstrID(i)]);
}

// ============================= EvaluateSurr =============================

void cProbSurr :: EvaluateLCB(int *algvar, cVector &c, cVector &fobjs)
{
  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(algvar,xn);

  // Get surrogate responses.
  cVector y(AppOut->GetNumAppOut( ));
  y[0] = SurMod->EvalLCB(xn);


  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    fobjs[0] = y[0];
  else
    HighFidelity->EvalExactFobj(algvar, fobjs[0]);

  // Constraints evaluation.

  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    c[AppOut->GetAppConstrID(i)] = y[AppOut->GetAppConstrOutID(i)];

  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,
                                      c[AppOut->GetExactConstrID(i)]);
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluateEI(int *algvar, cVector &c, cVector &fobjs)
{
  // Auxiliary variables.
  double ei;                  // Expected improvement.
  double hof;                 // High-fidelity exact objective function.
  double ci;                  // High-fidelity exact constraint.

  // Get normalized input variables.
  cVector xn(NumVar);

  HighFidelity->GetNormVar(algvar,xn);

  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    ei = SurMod->EvalExpImp(xn,currbest);
  else
  {
    HighFidelity->EvalExactFobj(algvar, hof);
    ei = (hof < currbest) ? abs(currbest - hof) : 0.0;
  }

  // Apply problem probability of feasibility (approximated constraints).
  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    ei *= GetConstraintFactor(xn, AppOut->GetAppConstrOutID(i)); //SurMod->EvalProbFeas(xn,AppOut->GetAppConstrOutID(i)); LEO

  // Zero EI for violated extact constraints.
  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
  {
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,ci);
    if (ci > 0.0) ei = 0.0;
  }

  ei += 1e-10;    // Avoid -infinity result.
  ei = log10(ei);
  fobjs[0] = -ei;
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluatePoI(int *algvar, cVector &c, cVector &fobjs)
{
  // Auxiliary variables.
  double poi;                 // Probability of Improvement.
  double hof;                 // High-fidelity exact objective function.
  double ci;                  // High-fidelity exact constraint.

  // Get normalized input variables.
  cVector xn(NumVar);
  HighFidelity->GetNormVar(algvar,xn);

  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    poi = SurMod->EvalProbImp(xn,currbest);
  else
  {
    HighFidelity->EvalExactFobj(algvar, hof);
    poi = (hof < currbest) ? 1.0 : 0.0;
  }

  // Apply problem probability of feasibility (approximated constraints).
  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    poi *= GetConstraintFactor(xn, AppOut->GetAppConstrOutID(i)); //SurMod->EvalProbFeas(xn,AppOut->GetAppConstrOutID(i)); LEO

  // Zero EI for violated extact constraints.
  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
  {
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,ci);
    if (ci > 0.0) poi = 0.0;
  }

  fobjs[0] = -poi;
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluateS2(int *algvar, cVector &c, cVector &fobjs)
{
  cVector xn(NumVar);
  HighFidelity->GetNormVar(algvar,xn);

  vector<cVector> sy;
  int no;
  int ns;
  SurMod -> GetNumOut(no);
  SurMod -> GetSampleY(sy);
  SurMod -> GetNumSample(ns);
  cVector xb(NumVar);
  cVector yb(no);
  cVector y(no);

  c[0] = 0;
  SurMod ->GetSampData( ).Normalize(xn);

  // Objective function evaluation.

  y[0] = SurMod -> SSqrSur(xn, 0);

  /*double pf = 0;
  for (int i = 0; i < NumConstr; i++){
    pf = SurMod -> EvalConstraintPF(xn, i + 1, 0);
    y[0] = pf*y[0];
  }*/

  y[0] = log10(y[0]);

  double fobj = -y[0];
  fobjs[0] = fobj;
}

// =======================INTEGER-MATRIX VARIABLE=========================

// ============================= Evaluate ================================

void cProbSurr :: Evaluate(int **algvar, cVector &c, cVector &fobjs)
{
  switch (state)
  {
    case EVALUATE_SURROGATE:
      EvaluateSurr(algvar, c, fobjs);
    break;

    case EVALUATE_EXPECTED_IMPROVEMENT:
      EvaluateEI(algvar, c, fobjs);
      c.Zero( );
    break;

    case EVALUATE_LOWER_CONFIDENCE_BOUND:
      EvaluateLCB(algvar, c, fobjs);
    break;

    case EVALUATE_PROBABILITY_IMPROVEMENT:
      EvaluatePoI(algvar, c, fobjs);
      c.Zero( );
    break;

  case EVALUATE_ESTIMATED_QUADRATIC_ERROR:
    EvaluateS2(algvar, c, fobjs);
    c.Zero( );
  break;
  }
}

// ============================= EvaluateSurr =============================

void cProbSurr :: EvaluateSurr(int **algvar, cVector &c, cVector &fobjs)
{
  // Get normalized input variables.
  cVector xn(NumVarEff);
  HighFidelity->GetNormVar(algvar,xn);

  // Get surrogate responses.
  cVector y(AppOut->GetNumAppOut( ));
  SurMod->Evaluate(xn,y);


  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    fobjs[0] = y[0];
  else
    HighFidelity->EvalExactFobj(algvar, fobjs[0]);

  // Constraints evaluation.

  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    c[AppOut->GetAppConstrID(i)] = y[AppOut->GetAppConstrOutID(i)];

  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,
                                      c[AppOut->GetExactConstrID(i)]);
}

// ============================= EvaluateSurr =============================

void cProbSurr :: EvaluateLCB(int **algvar, cVector &c, cVector &fobjs)
{
  // Get normalized input variables.
  cVector xn(NumVarEff);
  HighFidelity->GetNormVar(algvar,xn);

  // Get surrogate responses.
  cVector y(AppOut->GetNumAppOut( ));
  y[0] = SurMod->EvalLCB(xn);


  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    fobjs[0] = y[0];
  else
    HighFidelity->EvalExactFobj(algvar, fobjs[0]);

  // Constraints evaluation.

  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    c[AppOut->GetAppConstrID(i)] = y[AppOut->GetAppConstrOutID(i)];

  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,
                                      c[AppOut->GetExactConstrID(i)]);
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluateEI(int **algvar, cVector &c, cVector &fobjs)
{
  // Auxiliary variables.
  double ei;                  // Expected improvement.
  double hof;                 // High-fidelity exact objective function.
  double ci;                  // High-fidelity exact constraint.
  // Get normalized input variables.
  cVector xn(NumVarEff);

  HighFidelity->GetNormVar(algvar,xn);

  /*for (int i = 0; i < NumRow; i++)
  {
      for (int j = 0; j < NumCol; j++)
      {
          cout << algvar[i][j] << "  ";
      }
      cout << endl;
  }

  xn.Print( );*/

  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    ei = SurMod->EvalExpImp(xn,currbest);
  else
  {
    HighFidelity->EvalExactFobj(algvar, hof);
    ei = (hof < currbest) ? abs(currbest - hof) : 0.0;
  }

  // Apply problem probability of feasibility (approximated constraints).
  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    ei *= GetConstraintFactor(xn, AppOut->GetAppConstrOutID(i)); //SurMod->EvalProbFeas(xn,AppOut->GetAppConstrOutID(i)); LEO

  // Zero EI for violated extact constraints.
  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
  {
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,ci);
    if (ci > 0.0) ei = 0.0;
  }

  ei += 1e-10;    // Avoid -infinity result.
  ei = log10(ei);
  fobjs[0] = -ei;
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluatePoI(int **algvar, cVector &c, cVector &fobjs)
{
  // Auxiliary variables.
  double poi;                 // Probability of Improvement.
  double hof;                 // High-fidelity exact objective function.
  double ci;                  // High-fidelity exact constraint.

  // Get normalized input variables.
  cVector xn(NumVarEff);
  HighFidelity->GetNormVar(algvar,xn);

  // Objective function evaluation.

  if (AppOut->GetNumAppObjFunc( ))
    poi = SurMod->EvalProbImp(xn,currbest);
  else
  {
    HighFidelity->EvalExactFobj(algvar, hof);
    poi = (hof < currbest) ? 1.0 : 0.0;
  }

  // Apply problem probability of feasibility (approximated constraints).
  for (int i = 0; i < AppOut->GetNumAppConstr( ) ; i++)
    poi *= GetConstraintFactor(xn, AppOut->GetAppConstrOutID(i)); //SurMod->EvalProbFeas(xn,AppOut->GetAppConstrOutID(i)); LEO

  // Zero EI for violated extact constraints.
  for (int i = 0; i < AppOut->GetNumExactConstr( ) ; i++)
  {
    HighFidelity->EvalExactConstraint(AppOut->GetExactConstrID(i),algvar,ci);
    if (ci > 0.0) poi = 0.0;
  }

  fobjs[0] = -poi;
}

// ============================= Evaluate ================================

void cProbSurr :: EvaluateS2(int **algvar, cVector &c, cVector &fobjs)
{
  cVector xn(NumVarEff);
  HighFidelity->GetNormVar(algvar,xn);

  vector<cVector> sy;
  int no;
  int ns;
  SurMod -> GetNumOut(no);
  SurMod -> GetSampleY(sy);
  SurMod -> GetNumSample(ns);
  cVector xb(NumVarEff);
  cVector yb(no);
  cVector y(no);

  c[0] = 0;
  SurMod ->GetSampData( ).Normalize(xn);

  // Objective function evaluation.

  y[0] = SurMod -> SSqrSur(xn, 0);

  /*double pf = 0;
  for (int i = 0; i < NumConstr; i++){
    pf = SurMod -> EvalConstraintPF(xn, i + 1, 0);
    y[0] = pf*y[0];
  }*/

  y[0] = log10(y[0]);

  double fobj = -y[0];
  fobjs[0] = fobj;
}

// ============================= GetConstraintFactor ================================
// LEO
double cProbSurr :: GetConstraintFactor(cVector x, int out)
{
  double f;
  if (ConstrMethod == INFINITE_PEN)
      f = SurMod->EvalInfPen(x, out);
  else if (ConstrMethod == POF_SCHONLAU)
      f = SurMod->EvalProbFeas(x, out);
  else if (ConstrMethod == POF_TUTUM)
      f = SurMod->EvalProbFeasTutum(x, out);
  else if (ConstrMethod == POF_BAGHERI)
      f = SurMod->EvalProbFeasBagheri(x, out);
  else{
      cout << "Constraint handling function not defined." << endl;
      exit(0);
  }

  return f;
}

// =========================== End of file =================================

