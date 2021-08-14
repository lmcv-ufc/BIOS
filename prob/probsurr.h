// -------------------------------------------------------------------------
// probsurr.h - file containing the definition cProbSurr class.
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
// The aim of this class is to solve classical benchmark problems
// of numerical optimization and demonstrate how to implement a optimization
// problem using the BIOS system.
// -------------------------------------------------------------------------

#ifndef _PROBSURR_H
#define _PROBSURR_H

#include <vector>
#include "vec.h"
#include "mat.h"
#include "matvec.h"
#include "sysmat.h"
#include "rbf.h"
#include "krg.h"
#include "problem.h"
#include "group.h"

using namespace std;

// -------------------------------------------------------------------------
// Crossover Methods:
//
typedef enum
{
  INFINITE_PEN,
  POF_SCHONLAU,
  POF_TUTUM,
  POF_BAGHERI
} eConstrType;

// ------------------------------------------------------------------------
// Definition of cProbSurr class:
//
enum eProbSurrState
{
  EVALUATE_SURROGATE,
  EVALUATE_LOWER_CONFIDENCE_BOUND,
  EVALUATE_PROBABILITY_IMPROVEMENT,
  EVALUATE_EXPECTED_IMPROVEMENT,
  EVALUATE_ESTIMATED_QUADRATIC_ERROR
};

class cProbSurr : public cProblem
{
  protected:
  cProblem          *HighFidelity;
  sProbAppOut       *AppOut;
  cSURR             *SurMod;
  eProbSurrState     state;
  double             currbest;     // Current best feasible objective function.
  eConstrType        ConstrMethod; // LEO

  // For Integer Matrix variables
  int                NumVarEff;     // Number of effective variables
  int                NumRow;        // Number of variable rows
  int                NumCol;        // Number of variable columns

  public:

           cProbSurr(cSURR*,cProblem*,eProbSurrState); // LEO
           cProbSurr(cSURR*,cProblem*,sProbAppOut*,double,eProbSurrState); // LEO

  void     SetConstrMethod(eConstrType c){ConstrMethod = c;} // LEO
  double   GetConstraintFactor(cVector, int); // LEO

  int        VarNumRow(void) { return NumRow; }
  int        VarNumCol(void) { return NumCol; }
  cProblem*  GetHFP(void){ return HighFidelity; }

  // Continuous variable

  void     Evaluate(cVector &, cVector &, cVector &);

  void     EvaluateSurr(cVector &, cVector &, cVector &);
  void     EvaluateLCB(cVector &, cVector &, cVector &);
  void     EvaluateEI(cVector &, cVector &, cVector &);
  void     EvaluatePoI(cVector &, cVector &, cVector &);
  void     EvaluateS2(cVector &, cVector &, cVector &);

  void     GetDblBounds(double*, double*);

  // Integer-vector variable

  void     Evaluate(int *, cVector &, cVector &);

  void     EvaluateSurr(int *, cVector &, cVector &);
  void     EvaluateLCB(int *, cVector &, cVector &);
  void     EvaluateEI(int *, cVector &, cVector &);
  void     EvaluatePoI(int *, cVector &, cVector &);
  void     EvaluateS2(int *, cVector &, cVector &);

  void     GetVarBounds(int, double &, double&);

  // Integer-matrix variable

  void     Evaluate(int **, cVector &, cVector &);

  void     EvaluateSurr(int **, cVector &, cVector &);
  void     EvaluateLCB(int **, cVector &, cVector &);
  void     EvaluateEI(int **, cVector &, cVector &);
  void     EvaluatePoI(int **, cVector &, cVector &);
  void     EvaluateS2(int **, cVector &, cVector &);

  void     GetBounds(int, int*, int*);

  //void     PrintVar(double*); TO DO
  //void     WriteVar(double*);
  //void     EvaluateHFM(cVector &, cVector &, cVector &);
};

#endif
