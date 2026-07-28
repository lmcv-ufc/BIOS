// -------------------------------------------------------------------------
// benchmark.cpp - implementation of the benchmark class.
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
// Modified:     14-Mar-2013    Evandro Parente Junior
//               Evaluate input parameter changed to int*.
//
// Modified:     18-Mar-2013    Evandro Parente Junior
//               Created Booth and Rastrigin problems.
//
// Modified:     04-Dec-2017    Marina Alves Maia
//               Created CONSTR, TNK, ZDT1, ZDT6, SCH and KUR multiobjective
//               problems.
//
// Modified:     18-Oct-2018    Leonardo Gonçalves Ribeiro
//               Created PeaksC and PeaksSAO problem.
//
// Modified:     07-Jul-2019    Marina Alves Maia
//               Created BraninSAO, GoldSteinSAO, ColvilleSAO, TridSAO,
//               RastriginSAO, Hart3SAO and Hart6SAO problems.
//
// -------------------------------------------------------------------------

#include <cmath>
#include <math.h>
#include <iostream>

#include "benchmark.h"
#include "gbldef.h"
#include "gblvar.h"
#include "input.h"
#include "utl.h"
#include "vec.h"
#include "mat.h"
#include "matvec.h"
#include "sysmat.h"
#include "group.h"
#include "optsolution.h"

#include <vector>

using namespace std;


// -------------------------------------------------------------------------
// Static Variables:
//
int cBenchmark :: InpNumVar = 0;

// -------------------------------------------------------------------------
// Set read functions labels:
//
static const bool ReadFuncRegister[] =
{
  ProbMap( ).Insert("PROBLEM.NUMBER.VARIABLES",cBenchmark::ReadInpNumVar),
};
// -------------------------------------------------------------------------
// Register problems on the problem factory:
//
static const bool registeredProb[] =
{
  cProblemFactory :: Register("Peaks"             , MakeProb<cPeaksC,cPeaksD>),
  cProblemFactory :: Register("Branin"            , MakeProb<cBraninC,cBraninD>),
  cProblemFactory :: Register("Hartmann3"         , MakeProb<cHart3C,cHart3D>),
  cProblemFactory :: Register("Hartmann6"         , MakeProb<cHart6C,cHart6D>),
  cProblemFactory :: Register("Rastrigin"         , MakeProb<cRastriginC,cRastriginD>),
  cProblemFactory :: Register("ConstrainedBranin" , MakeProb<cConstrainedBraninC,cConstrainedBraninD>),
  cProblemFactory :: Register("Kitayama5"         , MakeProb<cKit5C,cKit5D>),
  cProblemFactory :: Register("Gano2"             , MakeProb<cGano2C>),
  cProblemFactory :: Register("G8"                , MakeProb<cG82C>),
  cProblemFactory :: Register("ThreeBarTruss"     , MakeProb<c3BarTrussC,c3BarTrussD>),
  cProblemFactory :: Register("NowackiBeam"       , MakeProb<cNowackiBeamC,cNowackiBeamD>),
  cProblemFactory :: Register("Beam"              , MakeProb<cBeamC,cBeamD>),
  cProblemFactory :: Register("TimoshenkoBeam"    , MakeProb<cBeamTimoshenkoC>),
  cProblemFactory :: Register("FGBeam"            , MakeProb<cFGBeam>),
  cProblemFactory :: Register("ColumnBuckRitz"    , MakeProb<cColumnBucklingRitzC>),
  cProblemFactory :: Register("ForresterA"        , MakeProb<cForresterAC>),
  cProblemFactory :: Register("ForresterB"        , MakeProb<cForresterBC>),
  cProblemFactory :: Register("ForresterC"        , MakeProb<cForresterCC>),
  cProblemFactory :: Register("ForresterD"        , MakeProb<cForresterDC>),
  cProblemFactory :: Register("Ackley5"           , MakeProb<cAckley5C>),
  cProblemFactory :: Register("CONSTR"            , MakeProb<cCONSTRC>),
  cProblemFactory :: Register("TNK"               , MakeProb<cTNKC>),
  cProblemFactory :: Register("ZDT6"              , MakeProb<cZDT6C>),
  cProblemFactory :: Register("ZDT1"              , MakeProb<cZDT1C>),
  cProblemFactory :: Register("SCH"               , MakeProb<cSCHC>),
  cProblemFactory :: Register("KUR"               , MakeProb<cKURC>)
};

// -------------------------------------------------------------------------
// Class cBenchmark:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cReadInpNumVar ==============================

void cBenchmark :: ReadInpNumVar(void)
{
  if (!(in >> InpNumVar))
  {
    cout << "Error in the input of the number of optimizations." << endl;
    exit(0);
  }
}

// =========================== cBenchmark ==================================

cBenchmark :: cBenchmark(void)
{
}

// ========================== cBenchmark ===================================

cBenchmark :: ~cBenchmark(void)
{
}

// -------------------------------------------------------------------------
// Class cBenchDiscrete:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cBenchDiscrete ==============================

cBenchDiscrete :: cBenchDiscrete(void)
{
}

// ========================== cBenchContinuous =============================

cBenchContinuous :: cBenchContinuous(void)
{
}

// =========================== ~cBenchDiscrete =============================

cBenchDiscrete :: ~cBenchDiscrete(void)
{
  delete []ListDim;
  delete []List;
}

// ========================== ~cBenchContinuous ============================

cBenchContinuous :: ~cBenchContinuous(void)
{
  delete[] Low;
  delete[] Upp;
}

// ============================== PrintVar =================================

void cBenchDiscrete :: PrintVar(int *algvar)
{
  for (int i = 0; i < NumVar; i++)
    cout << 1+i << " " << List[i][algvar[i]] << endl;

  cout << endl;
}

// ============================== PrintVar =================================

void cBenchDiscrete :: DecodeVar(int *algvar, cVector &xsamp)
{
  xsamp.Resize(NumVar);

  for (int i = 0; i < NumVar; i++)
      xsamp[i] = (List[i][algvar[i]]-List[i].Min())/(List[i].Max()-List[i].Min());
}

// ============================= WriteVar ==================================

void cBenchDiscrete :: WriteVar(int *algvar, ostream &out)
{
  for (int i = 0; i < NumVar; i++)
    out << List[i][algvar[i]] << "  ";

  out << endl;
}

// ============================= GetBounds =================================

void cBenchDiscrete :: GetBounds(int i, int *low, int *upp)
{
  *low = 0;
  *upp = ListDim[i]-1;
}

// ============================= GetBounds =================================

void cBenchDiscrete :: GetVarBounds(int i, double &low, double &upp)
{
  low = List[i].Min();
  upp = List[i].Max();
}

// ========================== GetBoundsDouble ==============================

void cBenchContinuous :: GetDblBounds(double *low,double *upp)
{
  for (int i = 0; i < NumVar; i++)
  {
    low[i] = Low[i];
    upp[i] = Upp[i];
  }
}

// -------------------------------------------------------------------------
// Class cPeaksC
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPeaksC ===============================

cPeaksC :: cPeaksC(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = -3;
    Upp[i] = 3;
  }
}

// ============================ Evaluate ==============================

void cPeaksC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  float a1 = 3*(1-x[0])*(1-x[0])*exp(-(x[0]*x[0])-(x[1]+1)*(x[1]+1));
  float a2 = -10*(x[0]/5 - x[0]*x[0]*x[0] - x[1]*x[1]*x[1]*x[1]*x[1])*exp(-(x[0]*x[0])-(x[1])*(x[1]));
  float a3 = -exp(-(x[0]+1)*(x[0]+1)-(x[1])*(x[1]));
  float fobj = a1+a2+a3/3;

  fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cPeaksD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPeaksD ==================================

cPeaksD :: cPeaksD(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-2;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(6.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = -3.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================= Evaluate ================================

void cPeaksD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Objective function evaluation.

  double a1 = 3*(1-x[0])*(1-x[0])*exp(-(x[0]*x[0])-(x[1]+1)*(x[1]+1));
  double a2 = -10*(x[0]/5 - x[0]*x[0]*x[0] - x[1]*x[1]*x[1]*x[1]*x[1])*exp(-(x[0]*x[0])-(x[1])*(x[1]));
  double a3 = -exp(-((x[0]+1)*(x[0]+1))-(x[1])*(x[1]));
  double fobj = a1+a2+a3/3;

  fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cBraninC
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBraninC ===============================

cBraninC :: cBraninC(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = -5;
  Upp[0] = 10;
  Low[1] = 0;
  Upp[1] = 15;
}

// ============================ Evaluate ==============================

void cBraninC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5.0/PI;
    double r = 6.0;
    double s = 10.0;
    double t = 1.0/(8.0*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cBraninD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBraninD ===============================

cBraninD :: cBraninD(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(15.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  List[0].Resize(ListDim[0]);
  List[0][0] = -5.0;
  for (int j = 1; j < ListDim[0]; j++) List[0][j] = List[0][j-1] + step;

  List[1].Resize(ListDim[1]);
  List[1][0] = 0.0;
  for (int j = 1; j < ListDim[1]; j++) List[1][j] = List[1][j-1] + step;
}

// ============================ Evaluate ==============================

void cBraninD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5/PI;
    double r = 6;
    double s = 10;
    double t = 1/(8*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    fobjs[0] = fobj;
}

// ============================= cHart3C ===============================

cHart3C :: cHart3C(void)
{
  NumVar = 3;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i < NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cHart3C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,3);
    cMatrix P(4,3);

    A[0][0] = 3.0;
    A[0][1] = 10.0;
    A[0][2] = 30.0;

    A[1][0] = 0.1;
    A[1][1] = 10.0;
    A[1][2] = 35.0;

    A[2][0] = 3.0;
    A[2][1] = 10.0;
    A[2][2] = 30.0;

    A[3][0] = 0.10;
    A[3][1] = 10.0;
    A[3][2] = 35.0;

    P[0][0] = pow(10,-4)*3689;
    P[0][1] = pow(10,-4)*1170;
    P[0][2] = pow(10,-4)*2673;

    P[1][0] =pow(10,-4)*4699;
    P[1][1] =pow(10,-4)*4387;
    P[1][2] =pow(10,-4)*7470;

    P[2][0] =pow(10,-4)*1091;
    P[2][1] =pow(10,-4)*8732;
    P[2][2] =pow(10,-4)*5547;

    P[3][0] =pow(10,-4)*381;
    P[3][1] =pow(10,-4)*5743;
    P[3][2] =pow(10,-4)*8828;

    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 3; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = - outer;
}

// ============================ Evaluate ==============================

void cHart3C :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    Evaluate(x, c, fobjs);
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double MA3 = 0.585 - 0.324*x1 - 0.379*x2 - 0.431*x3 - 0.208*x1*x2 + 0.326*x1*x3 + 0.193*x2*x3 + 0.225*x1*x1 + 0.263*x2*x2 + 0.274*x3*x3;

    fobjs[0] += 7.6*MA3;
}

// -------------------------------------------------------------------------
// Class cHart3D
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cHart3D ===============================

cHart3D :: cHart3D(void)
{
  NumVar = 3;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(1.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = 0.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================ Evaluate ==============================

void cHart3D :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,3);
    cMatrix P(4,3);

    A[0][0] = 3.0;
    A[0][1] = 10.0;
    A[0][2] = 30.0;

    A[1][0] = 0.1;
    A[1][1] = 10.0;
    A[1][2] = 35.0;

    A[2][0] = 3.0;
    A[2][1] = 10.0;
    A[2][2] = 30.0;

    A[3][0] = 0.10;
    A[3][1] = 10.0;
    A[3][2] = 35.0;

    P[0][0] = pow(10,-4)*3689;
    P[0][1] = pow(10,-4)*1170;
    P[0][2] = pow(10,-4)*2673;

    P[1][0] =pow(10,-4)*4699;
    P[1][1] =pow(10,-4)*4387;
    P[1][2] =pow(10,-4)*7470;

    P[2][0] =pow(10,-4)*1091;
    P[2][1] =pow(10,-4)*8732;
    P[2][2] =pow(10,-4)*5547;

    P[3][0] =pow(10,-4)*381;
    P[3][1] =pow(10,-4)*5743;
    P[3][2] =pow(10,-4)*8828;

    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 3; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = - outer;
}

// ============================= cHart6C ===============================

cHart6C :: cHart6C(void)
{
  NumVar = 6;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i < NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cHart6C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,6);
    cMatrix P(4,6);

    A[0][0] = 10;
    A[0][1] = 3;
    A[0][2] = 17;
    A[0][3] = 3.5;
    A[0][4] = 1.7;
    A[0][5] = 8;

    A[1][0] = 0.05;
    A[1][1] = 10;
    A[1][2] = 17;
    A[1][3] = 0.1;
    A[1][4] = 8;
    A[1][5] = 14;

    A[2][0] = 3;
    A[2][1] = 3.5;
    A[2][2] = 1.7;
    A[2][3] = 10;
    A[2][4] = 17;
    A[2][5] = 8;

    A[3][0] = 17;
    A[3][1] = 8;
    A[3][2] = 0.05;
    A[3][3] = 10;
    A[3][4] = 0.1;
    A[3][5] = 14;

    P[0][0] = pow(10,-4)*1312;
    P[0][1] = pow(10,-4)*1696;
    P[0][2] = pow(10,-4)*5569;
    P[0][3] = pow(10,-4)*124;
    P[0][4] = pow(10,-4)*8283;
    P[0][5] = pow(10,-4)*5886;

    P[1][0] =pow(10,-4)*2329;
    P[1][1] =pow(10,-4)*4135;
    P[1][2] =pow(10,-4)*8307;
    P[1][3] =pow(10,-4)*3736;
    P[1][4] =pow(10,-4)*1004;
    P[1][5] =pow(10,-4)*9991;

    P[2][0] =pow(10,-4)*2348;
    P[2][1] =pow(10,-4)*1451;
    P[2][2] =pow(10,-4)*3522;
    P[2][3] =pow(10,-4)*2883;
    P[2][4] =pow(10,-4)*3047;
    P[2][5] =pow(10,-4)*6650;

    P[3][0] =pow(10,-4)*4047;
    P[3][1] =pow(10,-4)*8828;
    P[3][2] =pow(10,-4)*8732;
    P[3][3] =pow(10,-4)*5743;
    P[3][4] =pow(10,-4)*1091;
    P[3][5] =pow(10,-4)*381;


    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 6; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = -outer;
}

// -------------------------------------------------------------------------
// Class cHart6D
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cHart6D ===============================

cHart6D :: cHart6D(void)
{
  NumVar = 6;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(1.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = 0.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================ Evaluate ==============================

void cHart6D :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    cVector alpha(4);
    alpha[0] = 1.0;
    alpha[1] = 1.2;
    alpha[2] = 3.0;
    alpha[3] = 3.2;

    cMatrix A(4,6);
    cMatrix P(4,6);

    A[0][0] = 10;
    A[0][1] = 3;
    A[0][2] = 17;
    A[0][3] = 3.5;
    A[0][4] = 1.7;
    A[0][5] = 8;

    A[1][0] = 0.05;
    A[1][1] = 10;
    A[1][2] = 17;
    A[1][3] = 0.1;
    A[1][4] = 8;
    A[1][5] = 14;

    A[2][0] = 3;
    A[2][1] = 3.5;
    A[2][2] = 1.7;
    A[2][3] = 10;
    A[2][4] = 17;
    A[2][5] = 8;

    A[3][0] = 17;
    A[3][1] = 8;
    A[3][2] = 0.05;
    A[3][3] = 10;
    A[3][4] = 0.1;
    A[3][5] = 14;

    P[0][0] = pow(10,-4)*1312;
    P[0][1] = pow(10,-4)*1696;
    P[0][2] = pow(10,-4)*5569;
    P[0][3] = pow(10,-4)*124;
    P[0][4] = pow(10,-4)*8283;
    P[0][5] = pow(10,-4)*5886;

    P[1][0] =pow(10,-4)*2329;
    P[1][1] =pow(10,-4)*4135;
    P[1][2] =pow(10,-4)*8307;
    P[1][3] =pow(10,-4)*3736;
    P[1][4] =pow(10,-4)*1004;
    P[1][5] =pow(10,-4)*9991;

    P[2][0] =pow(10,-4)*2348;
    P[2][1] =pow(10,-4)*1451;
    P[2][2] =pow(10,-4)*3522;
    P[2][3] =pow(10,-4)*2883;
    P[2][4] =pow(10,-4)*3047;
    P[2][5] =pow(10,-4)*6650;

    P[3][0] =pow(10,-4)*4047;
    P[3][1] =pow(10,-4)*8828;
    P[3][2] =pow(10,-4)*8732;
    P[3][3] =pow(10,-4)*5743;
    P[3][4] =pow(10,-4)*1091;
    P[3][5] =pow(10,-4)*381;


    double outer = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        double inner = 0;
        for (int jj = 0; jj< 6; jj++)
        {
            double xj = x[jj];
            double Aij = A[ii][jj];
            double Pij = P[ii][jj];
            inner += Aij*pow(xj-Pij,2);
        }
        double neww = alpha[ii]*exp(-inner);
        outer += neww;
    }

    fobjs[0] = -outer;
}

// -------------------------------------------------------------------------
// Class cRastriginC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cRastriginC ================================

cRastriginC :: cRastriginC(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  if (InpNumVar) NumVar = InpNumVar;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for(int i = 0; i < NumVar; i++)
  {
    Low[i] = -5.12;
    Upp[i] = 5.12;
  }
}

// ============================= Evaluate ================================

void cRastriginC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double f = 0;

  for(int i=0; i < NumVar; i++)
    f += x[i]*x[i] - 10*cos(2*PI*x[i]);

  f += 10*NumVar;

  fobjs[0] = f;
}

// -------------------------------------------------------------------------
// Class cRastriginD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cRastriginD ===============================

cRastriginD :: cRastriginD(void)
{
  NumVar = 2;
  NumConstr = 0;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(2*5.12/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = -5.12;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================= Evaluate ================================

void cRastriginD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Objective function evaluation.

  double f = 0;

  for(int i=0; i < NumVar; i++)
    f += x[i]*x[i] - 10*cos(2*PI*x[i]);

  f += 10*NumVar;

  fobjs[0] = f;
}

// ============================= cBraninC ===============================

cConstrainedBraninC :: cConstrainedBraninC(void)
{
  NumVar = 2;
  NumConstr = 1;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = -5;
  Upp[0] = 10;
  Low[1] = 0;
  Upp[1] = 15;
}

// ============================ Evaluate ==============================

void cConstrainedBraninC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5.0/PI;
    double r = 6.0;
    double s = 10.0;
    double t = 1.0/(8.0*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    c[0] = 1 - (x[0] + 5)*x[1]/45;
    fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cConstrainedBraninC :: EvalExactConstraint(int index, cVector& x, double &c)
{
  // Single constraint evaluation.
    if (index == 0){
        c = 1 - (x[0] + 5)*x[1]/45;
    }
    else{
        cout << "Invalid index in EvalExactConstraint!";
        exit(0);
    }
}

// -------------------------------------------------------------------------
// Class cConstrainedBraninD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================= cConstrainedBraninD ============================

cConstrainedBraninD :: cConstrainedBraninD(void)
{
  NumVar = 2;
  NumConstr = 1;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-3;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(15.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  List[0].Resize(ListDim[0]);
  List[0][0] = -5.0;
  for (int j = 1; j < ListDim[0]; j++) List[0][j] = List[0][j-1] + step;

  List[1].Resize(ListDim[1]);
  List[1][0] = 0.0;
  for (int j = 1; j < ListDim[1]; j++) List[1][j] = List[1][j-1] + step;
}

// ============================ Evaluate ==============================

void cConstrainedBraninD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Objective function evaluation.

    double PI = 3.14159265359;
    double a = 1.0;
    double b = 5.1/(4*pow(PI, 2.0));
    double d = 5/PI;
    double r = 6;
    double s = 10;
    double t = 1/(8*PI);

    double fobj = a*pow(x[1] - b*pow(x[0], 2) + d*x[0] - r, 2) + s*(1-t)*cos(x[0]) + s;

    c[0] = 1 - (x[0] + 5)*x[1]/45;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cGano2C
// -------------------------------------------------------------------------

// ============================= cGano2C ===============================

cGano2C :: cGano2C(void)
{
  NumVar = 2;
  NumConstr = 1;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0.10;
    Upp[i] = 10.0;
  }
}

// ============================ Evaluate ==============================

void cGano2C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // x[0] = 0.8846146; x[1] = 1.1500039;
  // Objective function evaluation.

  double fobj = 4*x[0]*x[0] + x[1]*x[1]*x[1] + x[0]*x[1];
  c[0] = 1/x[0] + 1/x[1] - 2.0;

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cGano2C :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
  // x[0] = 0.8846146; x[1] = 1.1500039;
  // Objective function evaluation.

  double fobj = 4*(x[0] + 0.1)*(x[0] + 0.1) + (x[1] + 0.1)*(x[1] + 0.1)*(x[1] + 0.1) + x[0]*x[1] + 0.1;
  c[0] = 1/x[0] + 1/(x[1] + 0.1) - 2.0 - 0.001;

  fobjs[0] = fobj;
}

// ============================ EvalExactFobj ==============================

void cGano2C :: EvalExactFobj(cVector &x, double &fobj)
{
  fobj = 4*x[0]*x[0] + x[1]*x[1]*x[1] + x[0]*x[1];
}

// ============================ EvalExactConstraint ==============================

void cGano2C :: EvalExactConstraint(int index, cVector& x, double &c)
{
  // Single constraint evaluation.
  if (index == 0){
      c = 1/x[0] + 1/x[1] - 2.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cGano2C :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
}

// -------------------------------------------------------------------------
// Class cG82C
// -------------------------------------------------------------------------

// ============================= cG82C ===============================

cG82C :: cG82C(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] =  0.0;
    Upp[i] = 10.0;
  }
}

// ============================ Evaluate ==============================

void cG82C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // x[0] = 1.2279713; x[1] = 4.2453733;
  // Objective function evaluation.

  double sin1 = sin(2*PI*x[0]);
  double sin2 = sin(2*PI*x[1]);

  double div = (x[0]*x[0]*x[0]*(x[0] + x[1]));

  if (div <= 1e-15) div = 1e-15;

  double fobj = -sin1*sin1*sin1*sin2/div;
  c[0] = x[0]*x[0] - x[1] + 1.0;
  c[1] = 1.0 - x[0] + pow((x[1] - 4.0), 2.0);

  // cout << "fobj = " << fobj << endl;
  // cout << "c1   = " << c[0] << endl;
  // cout << "c2   = " << c[1] << endl;

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cG82C :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // x[0] = 1.2279713; x[1] = 4.2453733;
    // Objective function evaluation.

    double sin1 = sin(2*PI*x[0] + 0.1);
    double sin2 = sin(2*PI*x[1] - 1);

    double div = (pow(x[0]+0.1, 3.0)*(x[0] + x[1] - 0.1));

    if (div <= 1e-15) div = 1e-15;

    double fobj = -sin1*sin1*sin1*sin2/div;
    c[0] = x[0]*x[0] - (x[1] - 0.1) + 1.0;
    c[1] = 1.0 - x[0] + pow((x[1] - 4.1), 2.0);

    // cout << "fobj = " << fobj << endl;
    // cout << "c1   = " << c[0] << endl;
    // cout << "c2   = " << c[1] << endl;

    fobjs[0] = fobj;
}

// ============================ EvalExactFobj ==============================

void cG82C :: EvalExactFobj(cVector &x, double &fobj)
{
  fobj = 4*x[0]*x[0] + x[1]*x[1]*x[1] + x[0]*x[1];
}

// ============================ EvalExactConstraint ==============================

void cG82C :: EvalExactConstraint(int index, cVector& x, double &c)
{
  // Single constraint evaluation.
  if (index == 0){
      c = 1/x[0] + 1/x[1] - 2.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cG82C :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
  approxc[1] = 1;
}

// -------------------------------------------------------------------------
// Class cKit5C
// -------------------------------------------------------------------------

// ============================= cKit5C ===============================

cKit5C :: cKit5C(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cKit5C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
  c[0] = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  c[1] = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  c[2] = (10*x[0] + x[1])/7 - 1;

  fobjs[0] = fobj;
}


// ============================ Evaluate ==============================

void cKit5C :: EvalExactFobj(cVector &x, double &fobj)
{
  fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
}

// ============================ Evaluate ==============================

void cKit5C :: EvalExactConstraint(int index, cVector& x, double &c)
{
  // Single constraint evaluation.
  if (index == 0){
      c = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  }
  else if (index == 1){
      c = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  }
  else if (index == 2){
      c = (10*x[0] + x[1])/7 - 1;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cKit5C :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
  approxc[1] = 1;
  approxc[2] = 1;
}

// -------------------------------------------------------------------------
// Class cKit5D
// -------------------------------------------------------------------------

// ============================= cKit5D ===============================

cKit5D :: cKit5D(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  ListDim = new int[NumVar];
  double step = 1.0e-2;
  for (int i = 0; i < NumVar; i++) ListDim[i] = round(1.0/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (int i = 0; i < NumVar; i++)
  {
    List[i].Resize(ListDim[i]);
    List[i][0] = 0.0;
    for (int j = 1; j < ListDim[i]; j++) List[i][j] = List[i][j-1] + step;
  }
}

// ============================ Evaluate ==============================

void cKit5D :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];
  // Objective function evaluation.

  double fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
  c[0] = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  c[1] = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  c[2] = (10*x[0] + x[1])/7 - 1;

  fobjs[0] = fobj;
}


// ============================ Evaluate ==============================

void cKit5D :: EvalExactFobj(int *algvar, double &fobj)
{
    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    fobj = -(x[0] - 1)*(x[0] - 1) - (x[1] - 0.5)*(x[1] - 0.5);
}

// ============================ Evaluate ==============================

void cKit5D :: EvalExactConstraint(int index, int *algvar, double &c)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Single constraint evaluation.
  if (index == 0){
      c = (((x[0] - 3)*(x[0] - 3) + (x[1] + 2)*(x[1] + 2))*exp(-(pow(x[1],7))))/12 - 1;
  }
  else if (index == 1){
      c = ((x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 0.5)*(x[1] - 0.5))/0.2 - 1;
  }
  else if (index == 2){
      c = (10*x[0] + x[1])/7 - 1;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cKit5D :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
  approxc[1] = 1;
  approxc[2] = 1;
}

// -------------------------------------------------------------------------
// Class c3BarTrussC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ c3BarTrussC ===============================
c3BarTrussC :: c3BarTrussC(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = 0.1;
  Upp[0] = 5.0;
  Low[1] = 0.1;
  Upp[1] = 5.0;
}

// ============================ Evaluate ===============================

void c3BarTrussC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Truss parameters

    double P = 20; // Nodal load

    double syc = 15; // Maximum permissible stress in compression
    double syt = 20;  // Maximum permissible stress in tension

    double A1 = x[0]; // Area of truss 1, which is the same as truss 3
    double A2 = x[1]; // Area of truss 2

    double w = 2*A1*sqrt(2) + A2; // Total weigth of the truss (divided by the constant \rho)

    double s1 = P*(A2 + sqrt(2)*A1)/(sqrt(2)*A1*A1 + 2*A1*A2);
    double s2 = P*(1/(A1 + sqrt(2)*A2));
    double s3 = -P*(A2/(sqrt(2)*A1*A1 + 2*A1*A2));

    // Constraints evaluation.

    c[0] = s1/syt - 1;
    c[1] = s2/syt - 1;
    c[2] = s3/syc - 1;

    // Objetive function evaluation

    double fobj = w;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class c3BarTrussC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ c3BarTrussC ===============================
c3BarTrussD :: c3BarTrussD(void)
{
  NumVar = 2;
  NumConstr = 3;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  double step = 1.0e-3;

  ListDim = new int[NumVar];
  int i,j;
  for (i = 0; i < NumVar; i++) ListDim[i] = round((5.0 - 0.1)/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (i = 0; i < NumVar; i++) List[i].Resize(ListDim[i]);

  for (i = 0; i < NumVar; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 0.1 + step*j;
}

// ============================ Evaluate ===============================

void c3BarTrussD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Truss parameters

    double P = 20; // Nodal load

    double syc = 15; // Maximum permissible stress in compression
    double syt = 20;  // Maximum permissible stress in tension

    double A1 = x[0]; // Area of truss 1, which is the same as truss 3
    double A2 = x[1]; // Area of truss 2

    double w = 2*A1*sqrt(2) + A2; // Total weigth of the truss (divided by the constant \rho)

    double s1 = P*(A2 + sqrt(2)*A1)/(sqrt(2)*A1*A1 + 2*A1*A2);
    double s2 = P*(1/(A1 + sqrt(2)*A2));
    double s3 = -P*(A2/(sqrt(2)*A1*A1 + 2*A1*A2));

    // Constraints evaluation.

    c[0] = s1/syt - 1;
    c[1] = s2/syt - 1;
    c[2] = s3/syc - 1;

    // Objetive function evaluation

    double fobj = w;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cNowackiBeamC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cNowackiBeamC ===============================
cNowackiBeamC :: cNowackiBeamC(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = 5;
  Upp[0] = 50;
  Low[1] = 50;
  Upp[1] = 250;
}

// ============================ Evaluate ===============================

void cNowackiBeamC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Beam parameters

    double l = 1500;   // Length of the beam, mm
    double b = x[0];   // Beam width, mm
    double h = x[1];   // Beam height, mm
    double F = 5000;   // Tip load, N

    double E  = 216620; // Young's modulus, MPa

    double sy   = 240;  // Yield stress, MPa
    double dmax = 5;    // Maximum displacement, mm

    // Auxiliary variables

    double Iy   = b*h*h*h/12; // Cross-sectional inertia
    double Area = b*h;        // Cross-sectional area

    // Constraints evaluation.

    c[0] = F*l*l*l/(3*E*Iy*dmax) - 1; // Maximum displacement constraint
    c[1] = 6*F*l/(b*h*h*sy) - 1;      // Maximum bending stress constraint

    // Objetive function evaluation

    double fobj = Area;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cNowackiBeamD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cNowackiBeamD ===============================
cNowackiBeamD :: cNowackiBeamD(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  double step = 1.0e-3;

  ListDim = new int[NumVar];
  int i,j;
  for (i = 0; i < 1; i++) ListDim[i] = round((50.0 - 5.0)/step) + 1;
  for (i = 1; i < 2; i++) ListDim[i] = round((250.0 - 50.0)/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (i = 0; i < NumVar; i++) List[i].Resize(ListDim[i]);

  for (i = 0; i < 1; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 5.0 + step*j;
  for (i = 1; i < 2; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 50.0 + step*j;
}

// ============================ Evaluate ===============================

void cNowackiBeamD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
    // Decodification of problem variables.

    cVector x(NumVar);
    for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

    // Beam parameters

    double l = 1500;   // Length of the beam, mm
    double b = x[0];   // Beam width, mm
    double h = x[1];   // Beam height, mm
    double F = 5000;   // Tip load, N

    double E  = 216620; // Young's modulus, MPa

    double sy   = 240;  // Yield stress, MPa
    double dmax = 5;    // Maximum displacement, mm

    // Auxiliary variables

    double Iy   = b*h*h*h/12; // Cross-sectional inertia
    double Area = b*h;        // Cross-sectional area

    // Constraints evaluation.

    c[0] = F*l*l*l/(3*E*Iy*dmax) - 1; // Maximum displacement constraint
    c[1] = 6*F*l/(b*h*h*sy) - 1;      // Maximum bending stress constraint

    // Objetive function evaluation

    double fobj = Area;

    fobjs[0] = fobj;
}

// -------------------------------------------------------------------------
// Class cBeamC
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBeamC ===============================

cBeamC :: cBeamC(void)
{
  NumVar = 10;
  NumConstr = 11;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<5; i++)
  {
    Low[i] = 35;
    Upp[i] = 65;
  }

  for (int i = 5; i<10; i++)
  {
    Low[i] = 1;
    Upp[i] = 4;
  }
}

// ============================ Evaluate ==============================

void cBeamC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Beam rotation.

  double teta[5];
  teta[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    teta[i] = ((P*l)/(E*I[i]))*(L + (l/2.0) - aux[i]) + teta[i-1];
  }

  // Beam displacement.

  double displ[5];
  displ[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    displ[i] = ((P*l*l)/(2.0*E*I[i])*(L-aux[i]+(2.0*l/3.0))) + (teta[i-1]*l) + displ[i-1];
  }

  // Stress constraints.

  int nc = 0;
  for (int i = 0; i < 5; i++)
  {
    c[nc++] = sig[i]/Sall - 1.0;
  }

  // Geometry constraints.

  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    c[nc++] = h/(20.0*b) - 1.0;
  }

  // Displacement constraint.

  c[nc++] = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;

  // Objetive function evaluation.

  fobjs[0] = (x[0]*x[5] + x[1]*x[6] + x[2]*x[7] + x[3]*x[8] + x[4]*x[9])*l;
}

// ============================ Evaluate ==============================

void cBeamC :: EvalExactConstraint(int index, cVector &x, double &c)
{
  // Decodification of problem variables.

  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Single constraint evaluation.
  if (index < 5){
      c = sig[index]/Sall - 1.0;
  }
  else if (index >=5 && index < 10){
      double h = x[index - 5];
      double b = x[index];
      c = h/(20.0*b) - 1.0;
  }
  else if (index == 10){
      c = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cBeamC :: GetApproxConstr(bool* approxc)
{
  approxc[0]  = 0;
  approxc[1]  = 0;
  approxc[2]  = 0;
  approxc[3]  = 0;
  approxc[4]  = 0;
  approxc[5]  = 0;
  approxc[6]  = 0;
  approxc[7]  = 0;
  approxc[8]  = 0;
  approxc[9]  = 0;
  approxc[10] = 0;
  approxc[11] = 0;
}

// -------------------------------------------------------------------------
// Class cBeamD
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cBeamD ===============================
cBeamD :: cBeamD(void)
{
  NumVar = 10;
  NumConstr = 11;
  NumObj = 1;

  // Create an array with the size of the list of each variable.

  double step = 1.0e-3;

  ListDim = new int[NumVar];
  int i,j;
  for (i = 0; i < 5; i++) ListDim[i] = round((65 - 35)/step) + 1;
  for (i = 5; i < 10; i++) ListDim[i] = round((4 - 1)/step) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];
  for (i = 0; i < NumVar; i++) List[i].Resize(ListDim[i]);

  for (i = 0; i < 5; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 35 + step*j;

  for (i = 5; i < 10; i++)
    for (j = 0; j < List[i].Dim(); j++) List[i][j] = 1 + step*j;

}
// ========================= Evaluate ================================

void cBeamD :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Beam rotation.

  double teta[5];
  teta[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    teta[i] = ((P*l)/(E*I[i]))*(L + (l/2.0) - aux[i]) + teta[i-1];
  }

  // Beam displacement.

  double displ[5];
  displ[0] = 0.0;
  for (int i = 1; i < 5; i++)
  {
    displ[i] = ((P*l*l)/(2.0*E*I[i])*(L-aux[i]+(2.0*l/3.0))) + (teta[i-1]*l) + displ[i-1];
  }

  // Stress constraints.

  int nc = 0;
  for (int i = 0; i < 5; i++)
  {
    c[nc++] = sig[i]/Sall - 1.0;
  }

  // Geometry constraints.

  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    c[nc++] = h/(20.0*b) - 1.0;
  }

  // Displacement constraint.

  c[nc++] = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;

  // Objetive function evaluation.

  double fobj = (x[0]*x[5] + x[1]*x[6] + x[2]*x[7] + x[3]*x[8] + x[4]*x[9])*l;

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cBeamD :: EvalExactConstraint(int index, int *algvar, double &c)
{
  // Decodification of problem variables.

  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) x[i] = List[i][algvar[i]];

  // Constraint evaluation.

  int P = 50000.0;            // Tip load (N)
  int L = 500.0;              // Total length (cm)
  int l = 100.0;              // Segment length (cm)
  double E = 2.0e7;           // Young's modulus (N/cm2)
  double Sall = 14000.0;      // Allowable stress (N/cm2)
  double vall = 2.7;          // Allowable displacement (cm)

  double aux[5];
  for (int i = 0; i < 5; i++)
  {
    aux[i] = l*(i+1);
  }

  // Bending moment.

  double M[5];
  for (int i = 0; i < 5; i++)
  {
    M[i] = P*(L + l - aux[i]);
  }

  // Inertia

  double I[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double b = x[i+5];
    I[i] = b*h*h*h/12.0;
  }

  // Stress.

  double sig[5];
  for (int i = 0; i < 5; i++)
  {
    double h = x[i];
    double W = I[i]/(h/2.0);
    sig[i] = M[i]/W;
  }

  // Single constraint evaluation.
  if (index < 5){
      c = sig[index]/Sall - 1.0;
  }
  else if (index >=5 && index < 10){
      double h = x[index - 5];
      double b = x[index];
      c = h/(20.0*b) - 1.0;
  }
  else if (index == 10){
      c = (P*l*l*l/(3*E))*(61/I[0] + 37/I[1] + 19/I[2] + 7/I[1] + 1/I[0])/vall - 1.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cBeamD :: GetApproxConstr(bool* approxc)
{
  approxc[0]  = 0;
  approxc[1]  = 0;
  approxc[2]  = 0;
  approxc[3]  = 0;
  approxc[4]  = 0;
  approxc[5]  = 0;
  approxc[6]  = 0;
  approxc[7]  = 0;
  approxc[8]  = 0;
  approxc[9]  = 0;
  approxc[10] = 0;
  approxc[11] = 0;
}

// ============================= cForresterC ===============================

cForresterAC :: cForresterAC(void)
{
  NumVar = 1;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cForresterAC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double a1   = 6.0*x[0] - 2.0;
  double fobj = a1*a1*sin(2*a1);

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cForresterAC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double A = 0.50;
  double B = 10.0;
  double C = -5.0;
  double D = 0.00;

  double a1   = A*pow(6.0*(x[0] + D) - 2.0, 2.0);
  double a2   = sin(12*(x[0] + D) - 4);
  double a3   = B*(x[0] + D - 0.5);
  double fobj = a1*a2 + a3 + C;

  fobjs[0] = fobj;
}

// ============================= cForresterBC ===============================

cForresterBC :: cForresterBC(void)
{
  NumVar = 1;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cForresterBC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double a1   = 6.0*x[0] - 2.0;
  double fobj = a1*a1*sin(2*a1);

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cForresterBC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double a1   = 6.0*x[0] - 2.0;
    double fobj = a1*a1*sin(2*a1);

    fobj = fobj - 5.0;

    fobjs[0] = fobj;
}

// ============================= cForresterCC ===============================

cForresterCC :: cForresterCC(void)
{
  NumVar = 1;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cForresterCC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double a1   = 6.0*x[0] - 2.0;
  double fobj = a1*a1*sin(2*a1);

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cForresterCC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double a1   = 6.0*(x[0] + 2.0) - 2.0;
    double fobj = a1*a1*sin(2*a1);

    fobjs[0] = fobj;
}

// ============================= cForresterCC ===============================

cForresterDC :: cForresterDC(void)
{
  NumVar = 1;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i<NumVar; i++)
  {
    Low[i] = 0;
    Upp[i] = 1;
  }
}

// ============================ Evaluate ==============================

void cForresterDC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  double a1   = 6.0*x[0] - 2.0;
  double fobj = a1*a1*sin(2*a1);

  fobjs[0] = fobj;
}

// ============================ Evaluate ==============================

void cForresterDC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

    double a1   = 6.0*x[0] - 2.0;
    double fobj = 3.0*a1*a1*sin(2*a1);

    fobjs[0] = fobj;
}

// ============================= cAckley5C ===============================

cAckley5C :: cAckley5C(void)
{
  NumVar = 5;
  NumConstr = 0;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0 ; i < NumVar; i++)
  {
    Low[i] = -2;
    Upp[i] =  2;
  }
}

// ============================ Evaluate ==============================

void cAckley5C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    double a, b, r;
    a = 20.0; b = 0.2; r = 2*PI;

    double aux1, aux2;
    aux1 = aux2 = 0;
    for (int i = 0; i < NumVar; i++)
    {
        aux1 += x[i]*x[i];
        aux2 += cos(r*x[i]);
    }

    double term1, term2, term3, term4;
    term1 = -a*exp(-b*sqrt((1.0/NumVar)*aux1));
    term2 = -1.0*exp((1.0/NumVar)*aux2);
    term3 = a;
    term4 = exp(1.0);

    fobjs[0] = term1 + term2 + term3 + term4;
}

// ============================ Evaluate ==============================

void cAckley5C :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation.

    Evaluate(x, c, fobjs);
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double x4 = x[3];
    double x5 = x[4];
    double MA5 = 0.585 - 0.00127*x1 - 0.00113*x2 - 0.00663*x3 - 0.0129*x4 - 0.00611*x5 + 0.00526*x1*x4 + 0.0106*x1*x5 - 0.000626*x2*x4 - 0.00310*x2*x5 - 0.00724*x4*x5 - 0.00096*x3*x3 - 0.0124*x4*x4 - 0.0101*x5*x5;

    fobjs[0] += 0.74*MA5;
}

// ============================= cBeamTimoshenko ===============================

cBeamTimoshenkoC :: cBeamTimoshenkoC(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = 400;  Low[1] = 1000;
  Upp[0] = 1000; Upp[1] = 2000;
}

// ============================ Evaluate ==============================

void cBeamTimoshenkoC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation

    double L = 5000;     // mm; Beam length
    double E = 2.1662E1; // Pa; Young modulus
    double nu = 0.27;    // Poisson's coefficient
    double q = 10;       // N/mm; Linearly distributed load

    double G = E/(2*(1 + nu)); // Pa; Shear modulus

    // Beam cross section

    double b = x[0];       // mm; Height
    double h = x[1];       // mm; Width
    double A = b*h;        // mm^2; Area
    double I = b*h*h*h/12; // mm^4; Inertia

    double wT = 5*q*L*L*L*L/(384*E*I) + q*L*L/(8*G*A); // mm; Displacement

    // Constraint evaluations

    double CostPerArea = 1.0;
    double Cost = CostPerArea*b/1000*h/1000;
    double CostMax = 0.8;
    c[0] = Cost/CostMax - 1.0;
    c[1] = h/(3*b) - 1.0;

    fobjs[0] = wT;
}

// ============================ Evaluate ==============================

void cBeamTimoshenkoC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation

    double L = 5000;     // mm; Beam length
    double E = 2.1662E1; // MPa; Young modulus
    double nu = 0.27;    // Poisson's coefficient
    double q = 10;       // N/mm; Linearly distributed load

    double G = E/(2*(1 + nu)); // Pa; Shear modulus

    // Beam cross section

    double b = x[0];       // mm; Height
    double h = x[1];       // mm; Width
    double A = b*h;        // mm^2; Area
    double I = b*h*h*h/12; // mm^4; Inertia

    double wE = 5*q*L*L*L*L/(384*E*I); // mm; Displacement

    // Constraint evaluations

    double CostPerArea = 1.0;
    double Cost = CostPerArea*b/1000*h/1000;
    double CostMax = 0.8;
    c[0] = Cost/CostMax - 1.0;
    c[1] = h/(3*b) - 1.0;

    fobjs[0] = wE;
}

// ============================ Evaluate ==============================

void cBeamTimoshenkoC :: EvalExactConstraint(int index, cVector& x, double &c)
{
  double b = x[0];
  double h = x[1];

  // Single constraint evaluation.
  if (index == 0){
      double CostPerArea = 1.0;
      double Cost = CostPerArea*b/1000*h/1000;
      double CostMax = 0.8;
      c = Cost/CostMax - 1.0;
  }
  else if (index == 1){
      c = h/(3*b) - 1.0;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cBeamTimoshenkoC :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 0;
  approxc[1] = 0;
}

// ============================= cBeamTimoshenko ===============================

cFGBeam :: cFGBeam(void)
{
  NumVar = 2;
  NumConstr = 2;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = 0.20; Low[1] = 0.00;
  Upp[0] = 1.00; Upp[1] = 10.0;
}

// ============================ Evaluate ==============================

void cFGBeam :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation

    double L  = 10;      // m; Beam length
    double Ec = 380e9;   // Pa; Young modulus
    double Em = 90e9;    // Pa; Young modulus
    double P  = 40e3;    // N; Tip load

    // Beam cross section

    double r = x[0];       // m; Radius
    double N = x[1];       // Power-law exponent

    double A = PI*r*r;        // m^2; Area
    double I = PI*r*r*r*r/2;  // m^4; Inertia

    // Displacement calculation (MCU using numerical integration)

    double ng = 10;
    cVector rvec, wvec;

    GaussPts1D(ng, rvec, wvec);

    double delta = 0.0;
    double J = L/2;

    for (int i = 0; i < ng; i++)
    {
        double xb = (rvec[i] + 1)*J;

        double Vm = pow(xb/L, N);
        double E  = Vm*Em + (1 - Vm)*Ec;

        double M  = P*(L - xb);
        double M1 = L - xb;

        delta += wvec[i]*M*M1/(E*I)*J;
    }

    // Constraint evaluations

    double VmFrac = 1.0/(N + 1.0);
    double VcFrac = 1.0 - VmFrac;
    double V = A*L;

    double CostM = 1.0;
    double CostC = 5.0;
    double CostMax = 50.0;

    double Cost = V*VmFrac*CostM + V*VcFrac*CostC;
    c[0] = Cost/CostMax - 1.0;

    double VcMin = 0.5;
    c[1] = 1.0 - VcFrac/VcMin;

    fobjs[0] = delta*1000;
}

// ============================ Evaluate ==============================

void cFGBeam :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation

    double L  = 10;      // m; Beam length
    double Ec = 380e9;   // Pa; Young modulus
    double Em = 90e9;    // Pa; Young modulus
    double P  = 40e3;    // N; Tip load

    // Beam cross section

    double r = x[0];       // m; Radius
    double N = x[1];       // Power-law exponent

    double A = PI*r*r;        // m^2; Area
    double I = PI*r*r*r*r/2;  // m^4; Inertia

    // Displacement calculation (One term Rayleigh-Ritz solution)

    double a1 = -(L*N + L)*P/(4*Ec*I*N + 4*Em*I);

    double w = -(a1*L*L)*1000; // mm; Displacement

    // Constraint evaluations

    double VmFrac = 1.0/(N + 1.0);
    double VcFrac = 1.0 - VmFrac;
    double V = A*L;

    double CostM = 1.0;
    double CostC = 5.0;
    double CostMax = 50.0;

    double Cost = V*VmFrac*CostM + V*VcFrac*CostC;
    c[0] = Cost/CostMax - 1.0;

    double VcMin = 0.5;
    c[1] = 1.0 - VcFrac/VcMin;

    fobjs[0] = w;
}

// ============================ Evaluate ==============================

void cFGBeam :: EvalExactConstraint(int index, cVector& x, double &c)
{
  double r = x[0];
  double N = x[1];

  double VmFrac = 1.0/(N + 1.0);
  double VcFrac = 1.0 - VmFrac;

  // Single constraint evaluation.
  if (index == 0){
      double L = 10;
      double A = PI*r*r;        // m^2; Area
      double V = A*L;

      double CostM = 1.0;
      double CostC = 5.0;
      double CostMax = 50.0;

      double Cost = V*VmFrac*CostM + V*VcFrac*CostC;
      c = Cost/CostMax - 1.0;
  }
  else if (index == 1){
      double VcMin = 0.5;
      c = 1.0 - VcFrac/VcMin;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ========================= GetApproxConstr ==========================

void cFGBeam :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 0;
  approxc[1] = 0;
}

// ================================ GaussPts1D =================================

void cFGBeam :: GaussPts1D(int npg, cVector &r, cVector &w)
{
    r.Resize(npg);
    w.Resize(npg);

    if (npg == 1)
    {
        r[0] = 0.000000000000000;
        w[0] = 2.000000000000000;
    }
    else if (npg == 2)
    {
        r[0] = -0.577350269189626;
        r[1] = 0.577350269189626;
        w[0] = 1.000000000000000;
        w[1] = 1.000000000000000;
    }
    else if (npg == 3)
    {
        r[0] = -0.774596669241483;
        r[1] = 0.000000000000000;
        r[2] = 0.774596669241483;
        w[0] = 0.555555555555556;
        w[1] = 0.888888888888889;
        w[2] = 0.555555555555556;
    }
    else if (npg == 4)
    {
        r[0] = -0.861136311594053;
        r[1] = -0.339981043584856;
        r[2] = 0.339981043584856;
        r[3] = 0.861136311594053;
        w[0] = 0.347854845137454;
        w[1] = 0.652145154862546;
        w[2] = 0.652145154862546;
        w[3] = 0.347854845137454;
    }
    else if (npg == 5)
    {
        r[0] = -0.906179845938664;
        r[1] = -0.538469310105683;
        r[2] = 0.000000000000000;
        r[3] = 0.538469310105683;
        r[4] = 0.906179845938664;
        w[0] = 0.23692688505618;
        w[1] = 0.478628670499367;
        w[2] = 0.568888888888889;
        w[3] = 0.478628670499367;
        w[4] = 0.236926885056189;
    }
    else if (npg == 6)
    {
        r[0] = -0.932469514203152;
        r[1] = -0.661209386466264;
        r[2] = -0.238619186083197;
        r[3] = 0.238619186083197;
        r[4] = 0.661209386466264;
        r[5] = 0.932469514203152;
        w[0] = 0.171324492379170;
        w[1] = 0.36076157304813;
        w[2] = 0.467913934572691;
        w[3] = 0.467913934572691;
        w[4] = 0.360761573048139;
        w[5] = 0.171324492379170;
    }
    else if (npg == 7)
    {
        r[0] = -0.949107912342758;
        r[1] = -0.741531185599394;
        r[2] = -0.405845151377397;
        r[3] = 0.000000000000000;
        r[4] = 0.405845151377397;
        r[5] = 0.741531185599394;
        r[6] = 0.949107912342758;
        w[0] = 0.129484966168870;
        w[1] = 0.279705391489277;
        w[2] = 0.381830050505119;
        w[3] = 0.417959183673469;
        w[4] = 0.381830050505119;
        w[5] = 0.279705391489277;
        w[6] = 0.129484966168870;
    }
    else if (npg == 8)
    {
        r[0] = -0.960289856497536;
        r[1] = -0.796666477413627;
        r[2] = -0.525532409916329;
        r[3] = -0.183434642495650;
        r[4] = 0.183434642495650;
        r[5] = 0.525532409916329;
        r[6] = 0.796666477413627;
        r[7] = 0.960289856497536;
        w[0] = 0.101228536290376;
        w[1] = 0.222381034453375;
        w[2] = 0.313706645877887;
        w[3] = 0.362683783378362;
        w[4] = 0.362683783378362;
        w[5] = 0.313706645877887;
        w[6] = 0.222381034453375;
        w[7] = 0.101228536290376;
    }
    else if (npg == 9)
    {
        r[0] = -0.968160239507626;
        r[1] = -0.836031107326636;
        r[2] = -0.613371432700590;
        r[3] = -0.324253423403809;
        r[4] = 0.000000000000000;
        r[5] = 0.324253423403809;
        r[6] = 0.613371432700590;
        r[7] = 0.836031107326636;
        r[8] = 0.968160239507626;
        w[0] = 0.081274388361574;
        w[1] = 0.180648160694857;
        w[2] = 0.260610696402935;
        w[3] = 0.312347077040003;
        w[4] = 0.330239355001260;
        w[5] = 0.312347077040003;
        w[6] = 0.260610696402935;
        w[7] = 0.180648160694857;
        w[8] = 0.081274388361574;
    }
    else if (npg == 10)
    {
        r[0] = -0.973906528517172;
        r[1] = -0.865063366688985;
        r[2] = -0.679409568299024;
        r[3] = -0.433395394129247;
        r[4] = -0.148874338981631;
        r[5] = 0.148874338981631;
        r[6] = 0.433395394129247;
        r[7] = 0.679409568299024;
        r[8] = 0.865063366688985;
        r[9] = 0.973906528517172;
        w[0] = 0.066671344308688;
        w[1] = 0.149451349150581;
        w[2] = 0.219086362515982;
        w[3] = 0.269266719309996;
        w[4] = 0.295524224714753;
        w[5] = 0.295524224714753;
        w[6] = 0.269266719309996;
        w[7] = 0.219086362515982;
        w[8] = 0.149451349150581;
        w[9] = 0.066671344308688;
    }
    else if (npg == 11)
    {
        r[0] = -0.269543155952345;
        r[1] =-0.519096129206812;
        r[2] = -0.730152005574049;
        r[3] = -0.887062599768095;
        r[4] = -0.978228658146057;
        r[5] = 0.000000000000000;
        r[6] = 0.269543155952345;
        r[7] = 0.519096129206812;
        r[8] = 0.730152005574049;
        r[9] = 0.887062599768095;
        r[10] = 0.978228658146057;
        w[0] = 0.262804544510247;
        w[1] = 0.233193764591990;
        w[2] = 0.186290210927734;
        w[3] = 0.125580369464905;
        w[4] = 0.055668567116175;
        w[5] = 0.272925086777901;
        w[6] = 0.262804544510247;
        w[7] = 0.233193764591990;
        w[8] = 0.186290210927734;
        w[9] = 0.125580369464905;
        w[10] = 0.055668567116174;
    }
    else if (npg == 12)
    {
        r[0] = -0.125233408511469;
        r[1] = -0.367831498998180;
        r[2] = -0.587317954286617;
        r[4] =  -0.769902674194305;
        r[4] = -0.904117256370475;
        r[5] = -0.981560634246719;
        r[6] = 0.125233408511469;
        r[7] = 0.367831498998180;
        r[8] = 0.587317954286617;
        r[9] =  0.769902674194305;
        r[10] = 0.904117256370475;
        r[11] = 0.981560634246719;
        w[0] = 0.249147045813403;
        w[1] = 0.233492536538355;
        w[2] = 0.203167426723066;
        w[3] = 0.160078328543346;
        w[4] = 0.106939325995318;
        w[5] = 0.047175336386512;
        w[6] = 0.249147045813403;
        w[7] = 0.233492536538355;
        w[8] = 0.203167426723066;
        w[9] = 0.160078328543346;
        w[10] = 0.106939325995318;
        w[11] = 0.047175336386512;
    }
    else if (npg == 25){
        r[0]  = -0.122864692610710;
        r[1]  = -0.243866883720988;
        r[2]  = -0.361172305809388;
        r[3]  = -0.473002731445715;
        r[4]  = -0.577662930241223;
        r[5]  = -0.673566368473468;
        r[6]  = -0.759259263037357;
        r[7]  = -0.833442628760834;
        r[8]  = -0.894991997878275;
        r[9]  = -0.942974571228974;
        r[10] = -0.976663921459517;
        r[11] = -0.995556969790498;
        r[12] =  0.000000000000000;
        r[13] =  0.122864692610710;
        r[14] =  0.243866883720988;
        r[15] =  0.361172305809388;
        r[16] =  0.473002731445715;
        r[17] =  0.577662930241223;
        r[18] =  0.673566368473468;
        r[19] =  0.759259263037357;
        r[20] =  0.833442628760834;
        r[21] =  0.894991997878275;
        r[22] =  0.942974571228974;
        r[23] =  0.976663921459517;
        r[24] =  0.995556969790498;
        w[0]  =  0.122242442990310;
        w[1]  =  0.119455763535785;
        w[2]  =  0.114858259145712;
        w[3]  =  0.108519624474264;
        w[4]  =  0.100535949067051;
        w[5]  =  0.091028261982964;
        w[6]  =  0.080140700335001;
        w[7]  =  0.068038333812357;
        w[8]  =  0.054904695975835;
        w[9]  =  0.040939156701306;
        w[10] =  0.026354986615032;
        w[11] =  0.011393798501026;
        w[12] =  0.123176053726715;
        w[13] =  0.122242442990310;
        w[14] =  0.119455763535785;
        w[15] =  0.114858259145712;
        w[16] =  0.108519624474264;
        w[17] =  0.100535949067051;
        w[18] =  0.091028261982964;
        w[19] =  0.080140700335001;
        w[20] =  0.068038333812357;
        w[21] =  0.054904695975835;
        w[22] =  0.040939156701306;
        w[23] =  0.026354986615032;
        w[24] =  0.011393798501026;
    }
    else
    {
        cout << "Maximum number of gauss points is 10." << endl;
        exit(0);
    }
}

// ============================= cBeamTimoshenko ===============================

cColumnBucklingRitzC :: cColumnBucklingRitzC(void)
{
  NumVar = 2;
  NumConstr = 1;
  NumObj = 1;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = 0.1; Low[1] = 0.1;
  Upp[0] = 1.0; Upp[1] = 1.0;
}

// ============================ Evaluate ==============================

void cColumnBucklingRitzC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    x[0] = 0.4;
    x[1] = 0.4;

    // Objective function evaluation

    double L = 10;       // m; Column length
    double E = 2.1662E7; // Pa; Young modulus

    double rb = x[0];
    double rt = x[1];

    double ndiv = 4;

    double Pcrit = EvalFiniteDifferences(L, E, rb, rt, ndiv);

    cout << "P = " << Pcrit << endl;

    exit(0);

    double I0 = x[0];    // m4; Inertia
    double f  = x[1];
    double k  = 1.0 + f*f*f*f;    // Inertia reduction factor - I = I0*[1 - x/(kL)]

    double coeffa = -4*pow(L, 15)/2625.0;
    double coeffb = -26*E*I0*pow(L, 13)/(175*k) + 36*E*I0*pow(L, 13)/175;
    double coeffc = 908*E*E*I0*I0*pow(L, 11)/(175*k) - 212*E*E*I0*I0*pow(L, 11)/(175*k*k) - 768*E*E*I0*I0*pow(L, 11)/175;
    double coeffd = -72*E*E*E*I0*I0*I0*pow(L, 9)/(5*k) + 144*E*E*E*I0*I0*I0*pow(L, 9)/(25*k*k) - 12*E*E*E*I0*I0*I0*pow(L, 9)/(25*k*k*k) + 48*E*E*E*I0*I0*I0*pow(L, 9)/5;

    cVector roots(3);
    SolvePolynome3(coeffa, coeffb, coeffc, coeffd, roots);

    double P = roots.Min();  // kN; Buckling load

    // Constraint evaluations

    double Volume = 2.0*PI*PI*L*k*pow(2.0, 0.5)/(3.0*I0)*(pow((I0/PI),(3.0/2.0)) - pow(((I0*k - I0)/(k*PI)), (3.0/2.0)));
    double VolumeMax = 1.0;
    c[0] = Volume/VolumeMax - 1.0;

    double If = I0*(1 - 1/k);
    double ri = pow((2*If/PI), 1.0/4.0);
    double rimin = 0.15;
    c[1] = 1 - ri/rimin;

    fobjs[0] = -P;
}

// ============================ Evaluate ==============================

void cColumnBucklingRitzC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation

    double L = 10;       // m; Column length
    double E = 2.1662E7; // Pa; Young modulus

    double I0 = x[0];    // m4; Inertia
    double f  = x[1];
    double k  = 1.0 + f*f*f*f;    // Inertia reduction factor - I = I0*[1 - x/(kL)]

    double P = 3*E*I0*(2*L*k - L)/(2*L*L*L*k);

    // Constraint evaluations

    double Volume = 2.0*PI*PI*L*k*pow(2.0, 0.5)/(3.0*I0)*(pow((I0/PI),(3.0/2.0)) - pow(((I0*k - I0)/(k*PI)), (3.0/2.0)));
    double VolumeMax = 1.0;
    c[0] = Volume/VolumeMax - 1.0;

    double If = I0*(1 - 1/k);
    double ri = pow((2*If/PI), 1.0/4.0);
    double rimin = 0.15;
    c[1] = 1 - ri/rimin;

    fobjs[0] = -P;
}

// ============================ Evaluate ==============================

void cColumnBucklingRitzC :: EvalExactConstraint(int index, cVector& x, double &c)
{
    double L = 10;       // m; Column length
    double I0 = x[0];    // m4; Inertia
    double f  = x[1];
    double k  = 1.0 + f*f*f*f;    // Inertia reduction factor - I = I0*[1 - x/(kL)]

  // Single constraint evaluation.
  if (index == 0)
  {
      double Volume = 2.0*PI*PI*L*k*pow(2.0, 0.5)/(3.0*I0)*(pow((I0/PI),(3.0/2.0)) - pow(((I0*k - I0)/(k*PI)), (3.0/2.0)));
      double VolumeMax = 1.0;
      c = Volume/VolumeMax - 1.0;
  }
  else if (index == 1)
  {
      double If = I0*(1 - 1/k);
      double ri = pow((2*If/PI), 1.0/4.0);
      double rimin = 0.15;
      c = 1 - ri/rimin;
  }
  else{
      cout << "Definition of an exact constraint missing!";
      exit(0);
  }
}

// ============================ Evaluate ==============================

void cColumnBucklingRitzC :: SolvePolynome3(double a0, double a1, double a2, double a3, cVector &x)
{
   double y[21];
   double z[21];

   for (int k=0; k<21; k++)
   {
       y[0]=1;
       y[1]=1;
       y[2]=1;

       /*  Bernoulli's algorithm */

       y[k+3]=-(((a1*y[k+2])+(a2*y[k+1])+(a3*y[k]))/a0);
       // cout<<y[10]/y[9]<<endl;
   }

   double alpha1=y[20]/y[19];

   double b0=a0;
   double b1=a1+alpha1*b0;
   double b2=-(a3/alpha1);
   //cout<<b0<<" " <<b1<<" "<<b2<<endl;

   for(int j=0;j<21;j++)
   {
      z[0]=0;
      z[1]=1;
      z[j+2]=-((b1*z[j+1]+b2*z[j])/b0);
   }

   double alpha2=(z[20]/z[19]);
   //cout<<" The second solution is alpha2=  " <<alpha2<<endl;

   double c0=b0;
   double c1=-(b2/alpha2);
   double alpha3=-(c1/c0);
   //cout<<c0<<" " <<c1<<endl;
   //cout<<" The third solution will be alpha3 = "<< alpha3<< "\n"<<endl;

   x[0] = alpha1; x[1] = alpha2; x[2] = alpha3;
}

// ============================ EvalFiniteDifferences ==============================

double cColumnBucklingRitzC :: EvalFiniteDifferences(double L, double E, double rb, double rt, double ndiv)
{
    double h      = L/ndiv;
    int    nn     = ndiv + 1;
    int    nnodes = nn + 2;

    cVector Nodes(nnodes); cVector Pos(nnodes);
    for (int i = 0; i < nnodes; i++)
    {
        Nodes[i] = (double) i;
        Pos[i]   = (Nodes[i] - 1.0)*h;
    }

    cMatrix A(ndiv, nn + 1); cMatrix B(ndiv, nn + 1);
    A.Zero( ); B.Zero( );

    for (int i = 1; i < nn; i++)
    {
        A[i - 1][i - 1] = 1.0;
        A[i - 1][i]     = -2.0;
        A[i - 1][i + 1] = 1.0;

        double I = GetInertiaX(rb, rt, Pos[i], L);
        B[i - 1][i]     = -1.0/I;
        B[i - 1][nn]    =  1.0/I;
    }

    for (int i = 0; i < ndiv; i++)
    {
        A[i][2] = A[i][0] + A[i][2];
        B[i][2] = B[i][0] + B[i][2];
    }

    cMatrix A0(ndiv, nn - 1); cMatrix B0(ndiv, nn - 1);
    A0.Zero( ); B0.Zero( );

    for (int i = 0; i < ndiv; i++)
    {
        for (int j = 0; j < (nn - 1); j++)
        {
            A0[i][j] = A[i][j + 2];
            B0[i][j] = B[i][j + 2];
        }
    }



    return 1;
}

// ============================ EvalFiniteDifferences ==============================

double cColumnBucklingRitzC :: GetInertiaX(double rb, double rt, double x, double L)
{
    double r = rb + (rt - rb)*x/L;
    double I = PI*r*r*r*r/2;

    return I;
}

// -------------------------------------------------------------------------
// Class cCONSTRC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cCONSTRC ===============================

cCONSTRC :: cCONSTRC(void)
{
  NumVar = 2;
  NumObj = 2;
  NumConstr = 2;

  Low = new double[NumVar];
  Upp = new double[NumVar];
  Low[0] = 0.1;
  Upp[0] = 1.0;
  Low[1] = 0.0;
  Upp[1] = 5.0;
 }

// ============================= Evaluate ================================

void cCONSTRC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Constraint evaluation.

  c[0] =  6-x[1]-(9*x[0]);
  c[1] =  1+x[1]-9*x[0];

  // Objective functions evaluation.

  fobjs[0] = x[0];
  fobjs[1] = (1+ x[1])/(x[0]);

 // cout << "Fobjs " << fobjs[0] << "  " << fobjs[1] << endl;

}

// -------------------------------------------------------------------------
// Class cTNKC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cTNKC ===============================

cTNKC :: cTNKC(void)
{
  NumVar = 2;
  NumObj = 2;
  NumConstr = 2;

  Low = new double[NumVar];
  Upp = new double[NumVar];
  Low[0] = 0;
  Upp[0] = PI;
  Low[1] = 0;
  Upp[1] = PI;

 }

// ============================= Evaluate ================================

void cTNKC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Constraint evaluation.

  c[0] =  -pow(x[0],2) -pow(x[1],2) + 1 + 0.1*cos(16*atan(x[0]/x[1]));
  c[1] =  pow(x[0] - 0.5, 2)/0.5 + pow(x[1]-0.5, 2)/0.5 - 1;

  // Objective functions evaluation.

  fobjs[0] = x[0];
  fobjs[1] = x[1];

//  cout << "Fobjs " << fobjs[0] << "  " << fobjs[1] << endl;

}

// -------------------------------------------------------------------------
// Class cZDT6C:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cZDT6C ===============================

cZDT6C :: cZDT6C(void)
{
  NumVar = 10;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
      Low[i] = 0.0;
      Upp[i] = 1.0;
  }
}

// ============================= Evaluate ================================

void cZDT6C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  fobjs[0] = x[0];

  double aux = 0.0;

  for (int i = 1; i < NumVar; i++)
  {
      aux += x[i];
  }

  aux = 1 + (aux*9)/(NumVar-1);

  fobjs[1] = aux*(1-pow((fobjs[0]/aux),2));

}

// -------------------------------------------------------------------------
// Class cSCHC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cSCHC ===============================

cSCHC :: cSCHC(void)
{
  NumVar = 1;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  Low[0] = -1000;
  Upp[0] = 1000;

}

// ============================= Evaluate ================================

void cSCHC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  fobjs[0] = pow(x[0], 2);

  fobjs[1] = pow((x[0]-2),2);
}

// -------------------------------------------------------------------------
// Class cZDT1C:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cZDT1C ===============================

cZDT1C :: cZDT1C(void)
{
  NumVar = 30;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
      Low[i] = 0.0;
      Upp[i] = 1.0;
  }
}

// ============================= Evaluate ================================

void cZDT1C :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  fobjs[0] = x[0];

  double aux = 0.0;

  for (int i = 1; i < NumVar; i++)
  {
      aux += x[i];
  }

  double aux2 = 0.0;

  aux2 = 1 + (9/(NumVar-1))*aux;

  fobjs[1] = aux2*(1-sqrt((fobjs[0]/aux2)));
}



// -------------------------------------------------------------------------
// Class cKURC:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cKURC ===============================

cKURC :: cKURC(void)
{
  NumVar = 3;
  NumObj = 2;
  NumConstr = 0;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
      Low[i] = -5;
      Upp[i] = 5;
  }
}

// ============================= Evaluate ================================

void cKURC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{

  // Objective functions evaluation.

  double aux1 = 0.0;
  double aux2 = 0.0;

  aux1 = -0.2*sqrt((pow(x[0],2)+pow(x[1],2)));
  aux2 = -0.2*sqrt((pow(x[1],2)+pow(x[2],2)));

  fobjs[0] = -10*(exp(aux1)+exp(aux2));

  double aux3 = 0.0;

  for (int i = 0; i < NumVar; i++)
  {
      aux3 += pow(fabs(x[i]),0.8)+5*sin(pow(x[i],3));
  }

  fobjs[1] = aux3;
}

// =========================== End of file =================================

