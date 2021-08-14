// -------------------------------------------------------------------------
// rbf.h - file containing the definition of the cRBF class.
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
//
// The class cRBF contains data and methods relevant to create a surrogate
// model using Radial Basis Functions (RBF).
//
// Refs:
// FORRESTER, Alexander et al. Engineering design via surrogate modelling: a
// practical guide. John Wiley & Sons, 2008.
//
// KITAYAMA, Satoshi; ARAKAWA, Masao; YAMAZAKI, Koetsu. Sequential
// approximate optimization using radial basis function network for
// engineering  optimization. Optimization and Engineering,
// v. 12, n. 4, p. 535-557, 2011.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// void CreateModel(int nvar, int nout, int nsam, vector<cVector> &xs,
//                  vector<cVector> &ys)
//
//   nvar  -  number of problem variables                            (in)
//   nout  -  number of output responses                             (in)
//   nsam  -  number of samples                                      (in)
//   xs    -  array of sample points                                 (in)
//   ys    -  array with responses at each sample point              (in)
//
// This method creates a RBF model (R^nvar => R^nout) from a set of given
// samples.
// -------------------------------------------------------------------------
//
// void Evaluate(cVector &x, cVector &y)
//
//   x     -  given point (R^nvar)                                   (in)
//   y     -  responses (R^nout)                                     (out)
//
// This method evalutes the approximate responses of the RBF model
// at the given point.
//
// -------------------------------------------------------------------------
// Protected methods:
// -------------------------------------------------------------------------
//
// double CalcSigma(int nvar, int nsam, vector<cVector> xs)
//
//   nvar  -  number of problem variables                            (in)
//   nout  -  number of output responses                             (in)
//   nsam  -  number of samples                                      (in)
//
// This method returns the Gaussian parameter
// -------------------------------------------------------------------------

#ifndef _RBF_H
#define _RBF_H

#include <iostream>
#include <vector>
#include "optalg.h"
#include "vec.h"
#include "group.h"
#include "surr.h"

using namespace std;

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Definition of RBF class:
//
class cRBF : public cSURR
{
 protected:
  int              NumBasis;      // Number of basis
  double           Lambda;        // Regularization paramenter
  cMatrix          Sigma;         // Gaussian paramenter

  bool             c = true;

          void     CalcDistMax(cVector&);
          void     CalcDistMax(cVector&, double);
          void     CalcSigma(eSigmaType);
          void     CalcSigma(cMatrix, int);
          void     CalcSigmaNakayama(void);
          void     CalcSigmaKitayama(void);
          void     CalcSigmaKitayamaUniform(void);
          void     CalcSigmaKitayamaAdpScl(void);
          void     CalcSigmaLOOCV(void);
          void     CalcSigmakFOLDCV(void);

          void     CalcWeights(int, int, vector<cVector>, vector<cVector>);
          void     CalcH(vector<cVector>&);


 public:

          vector<cVector>  Weight;        // RBF Weights (NumOut x NumSample)
          cMatrix          Psi;           // Correlation matrix
          cMatrix          U;             // Inverted Correlation Matrix

          void             SetLambda(double l)  { Lambda = l; }
          void             DecompHmat(void);

                   cRBF(void);
                   cRBF(int, int, int, vector<cVector>, vector<double>);
  virtual          ~cRBF(void);
          void     GetSigma(cMatrix&) ;
          void     GetMu(double& mu, int out) { mu = Mu[out]; }
          void     GetSigmaSqr(double& sigmasqr, int out) { sigmasqr = SigmaSqr[out]; }
          void     CreateModel(const sSampData&, eSigmaType);
          void     CreateModel(const sSampData&,cMatrix);
          void     UpdateModel(eSigmaType, cVectorVec&,cVectorVec&);
          void     UpdateModel(eSigmaType, vector<cVector>, vector<cVector>, cVector, cVector);
          void     UpdateModel(cMatrix, vector<cVector>, vector<cVector>, cVector, cVector);
          void     Evaluate(cVector&, cVector&, vector<cVector> *vec = 0);
          void     MuSur(int);
          void     SigmaSqrSur(int);
          double   SSqrSur(cVector, int);
          double   CalcFi(double, int, int);
          double   EvalExpImp(cVector &,double);
          double   EvalLCB(cVector&);
          double   EvalProbImp(cVector &,double);
          double   EvalConstraintPF(cVector &, int, double){ }
          double   EvalInfPen(cVector&,int,double tol = 1e-6);           // LEO
          double   EvalProbFeas(cVector&,int,double tol = 1e-6);         // LEO
          double   EvalProbFeasTutum(cVector&,int,double tol = 1e-6);    // LEO
          double   EvalProbFeasBagheri(cVector&,int,double tol = 1e-6);  // LEO

 void              ReadSigmaType(std::istream&);

};

#endif
