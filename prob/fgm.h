// -------------------------------------------------------------------------
// fgm.h - file containing the definition of the Functionally Graded
//         Material class.
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
// The class cFGM contains data and routines relevant to optimization
// of functionally graded structures.
// -------------------------------------------------------------------------

#ifndef _FGM_H
#define _FGM_H

#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <vector>

using namespace std;

#include "problem.h"
#include "vec.h"

// -------------------------------------------------------------------------
// FGM models:
//
typedef enum
{
  VOIGT,
  REUSS,
  MORI_TANAKA
} eFGMModel;

// -------------------------------------------------------------------------
// FGM volume fraction distribution:
//
typedef enum
{
  POWER_LAW,
  SIGMOYD,
  POLYNOMIAL,
  PCHIP,
  BSPLINE
} eFGMVolDist;

// -------------------------------------------------------------------------
// Definition of cFGM class:
//
class cFGM : public cProblem
{
 protected:
  static  eFGMModel   FGMModel;       // FGM model used to define eff. properties
  static  eFGMVolDist FGMVolDist;     // FGM volume fraction distribution
  static  int         NumFGMMat;      // Number of FGM materials
  static  cVector     FGMMat;            // Material list
  static  cVector     FGIDMat;    // Custom angle values
  static  double      dExp;          // Power-law exponent variation
  static  double      MinExp;        // Minimum power-law exponent
  static  double      MaxExp;        // Maximum power-law exponent
          int         NumCP;         // Number of control points; If power-law/sygmoide/exponential, NumCP = 1;
          double      dThk;          // Thickness variation
          double      MinThk;        // Minimum thickness
          double      MaxThk;        // Maximum thickness
          bool        ThkVar;        // Is thickness a variable?

  int     *ListDim; // Number of discrete values for each variable
  cVector *List;    // Lists of discrete values for each variable
  double  *Low;
  double  *Upp;


 protected:
          void     GaussPts1D(int, cVector &, cVector &);
          void     VolumeDist(eFGMVolDist, int, cVector, cVector &, double exp = 0);
          void     VolumeDist(eFGMVolDist, cVector, int, cVector, cVector &);
          void     PowerLaw(int, double, cVector, cVector &);
          void     Sigmoyd(int, double, cVector, cVector &);
          void     Polynomial(int, cVector, cVector &);
          void     PiecewiseCubicInterpolation(cVector, int, cVector, cVector &);
          void     Bspline(cVector, int, cVector, cVector &);
          int      FindSpan(int, int, double, cVector);
          void     BasisFuns(int, double, int, cVector, cVector&);
          void     CurvePoint(int, int, cVector, cVector, double, double&);
          void     CurveDerivs(int, int, cVector, cVector, double, cVector&);
          void     BSplineDer(cVector, int, cVector, cVector&);
          void     DersBasisFuns(int, double, int, int, cVector, cMatrix &);
          void     EffPropModel(eFGMModel, cVector, cVector &, cVector &, cVector &, cVector &, cVector &);
          void     Voigt(cVector, cVector &, cVector &, cVector &, cVector &);
          void     VoigtDensity(cVector, cVector &);
          void     Reuss(cVector, cVector &, cVector &, cVector &, cVector &);
          void     MoriTanaka(cVector, cVector &, cVector &, cVector &, cVector &);
          void     EvalVolumeRatio(double, double &);
          void     EvalVolumeRatio(cVector, double &);
          void     EvalDens(double, double &);
          void     EvalDens(cVector, double &);
          void     EvalCost(double, double &);
        double     EvalYieldStress(double Vck);
          void     CalcABDG(double, cVector, cVector, cVector, cVector,cMatrix &, cMatrix &, cMatrix &, cMatrix &, cMatrix &);
          void     CalcMb(double, cVector, cVector, cVector, cMatrix &);
          void     QMatrix(double, double, cMatrix &, cMatrix &);
 virtual  void     PrintResult(int**, ostream&);

 public:
          void     ReadFGMaterials(std::istream&);
          void     ReadFGMVolDist(std::istream&);
          void     ReadFGMModel(std::istream&);
          void     ReadFGMExpRange(std::istream&);
          void     ReadNumberOfControlPoints(std::istream&);
          void     ReadFGMThkRange(std::istream&);
          int      GetFGMNumMat(void) { return NumFGMMat; }

                   cFGM(void);
  virtual         ~cFGM(void);
  virtual void     LoadReadFunc(cInpMap&);
  virtual void     Init(void);
  virtual void     PrintVar(int**);
  virtual void     WriteVar(int**);
          void     GetBounds(int, int*, int*);
          void     GetDblBounds(double*, double*);
};

#endif
