// -------------------------------------------------------------------------
// KRG.h - file containing the definition of the cKRG class.
// -------------------------------------------------------------------------
// Copyright (c) 2018 LMCV/UFC
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
// The class cKRG contains data and methods relevant to create a surrogate
// model using Ordinary Kriging.
//
// Refs:
// JONES, D. R. A taxonomy of global optimization methods based on response
// surfaces. Journal of Global Optimization, v. 21, n. 4, p. 345–383, Dec 2001.
// ISSN 1573-2916.
//
// FORRESTER, Alexander et al. Engineering design via surrogate modelling: a
// practical guide. John Wiley & Sons, 2008.
//
// -------------------------------------------------------------------------
// Protected methods:
// -------------------------------------------------------------------------
//
// void MuSur(double &mu, int out)
//
//   mu  -  average of sample data evaluated by the Gaussian      (in/out)
//          Process
//   out -  number of the output of interest                      (in)
//
// This method returns the Gaussian average parameter considering the
// sample set of the specified output.
// -------------------------------------------------------------------------
//
// void SigmaSqrSur(double &mu, double sigmasqr, int out)
//
//   mu  -  mean of sample data evaluated by the Gaussian      (in)
//          Process
//   sigmasqrsqr - variance of the Gaussian Process               (in/out)
//   out -  number of the output of interest                      (in)
//
// This method returns the variance of the Gaussian Process considering the
// sample set of the specified output.
// -------------------------------------------------------------------------
//
// double Likelihood(vector<cVector> sampx, cVector Theta, int out)
//
//   sampx  -  sampling plan data                                 (in)
//   Theta  - Hyper-parameters                                    (in)
//   out    - number of the output of interest                      (in)
//
// This method returns the likelihood of the specified sampling plan to
// be represented by the hyper-parameters theta.
// -------------------------------------------------------------------------
//
// void MaxLikelihood(cVector &bestParticlevec)
//
//   bestParticlevec  -  sampling plan data                       (out)
//
// This method returns the likelihood of the specified sampling plan to
// be represented by the hyper-parameters theta.
// -------------------------------------------------------------------------

#ifndef _HKRG_H
#define _HKRG_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "group.h"

#include "surr.h"
#include "krg.h"
#include "benchmark.h"
#include "optalg.h"
#include "stdpso.h"
#include "stdde.h"
#include "problike.h"

using namespace std;

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Definition of KRG class:
//
class cHIERKRG : public cSURR
{
 protected:

//  int             Normalize;     // Normalization
  double            Lambda;        // Regularization paramenter
//  cVector          Xlow;         // X lower bound
//  cVector          Xupp;         // X upper bound
  int               SubPop;
  int               SubMaxGen;
  double            SubTolViol;
  double            SubMutProb;
  cOptAlgorithm     *SubAlgType;
  int               SubStallGen;
  eDifType          SubDifType;
  eSwaTopType       SubTopology;
  eSolType          SubSolType;
  bool*             ApproxC;
  eCorrelationType  CorrType;

 // vector<double>   SigmaSqr;      // Stores the SigmaSqr; used on statistical infill methods
 // vector<double>    Mu;            // Stores the Mean; used on statistical infill methods

 public:
           vector<double>    SSqr;
           vector<double>    SSqrl;
           vector<cVector>   BestTheta;
           vector<cVector>   BestThetad;
           vector<cVector>   BestThetaOld;
           vector<cVector>   BestThetadOld;
           vector<double>    BestRho;
           vector<double>    Prediction;
           vector<cVector>   Weight;
           vector<cVector>   Weightl;
           vector<cVector>   Weightd;
           vector<cVector>   Weightl1;
           bool              c = true;

           double            HPlow;
           double            HPupp;
           
           cMatrix           U1;
           cMatrix           U2;
           cMatrix           U3;
           cMatrix           U;              // Inverse of Psi
           cMatrix           Psi;            // Correlation matrix

           cMatrix*          Psil;           // LF model correlation matrix
           cMatrix*          Psildc;         // Inverse of Psil
           cMatrix*          Psid;           // Difference model correlation matrix
           cMatrix*          Psiddc;         // Inverse of Psid

           cVector           SigmaSqrl;      // Stores the SigmaSqr; used on statistical infill methods
           cVector           Mul;            // Stores the Mean; used on statistical infill methods
           cVector           SigmaSqrd;      // Stores the SigmaSqr; used on statistical infill methods
           cVector           Mud;            // Stores the Mean; used on statistical infill methods
           cVector*          fVec;           // Stores the difference vectors

                 cHIERKRG(void);
                 cHIERKRG(int, int, int, vector<cVector>, vector<double>);
  virtual        ~cHIERKRG(void);
          void   GetSigma(cVector&) ;
          cVector  GetBestThetad(int out) {return BestThetad[out];}
          cVector  GetBestTheta(int out)  {return BestTheta[out];}
          void   CreateModel(const sSampData&,eCorrelationType,double,double);
          void   UpdateModel(eCorrelationType,cVectorVec&,cVectorVec&,cVectorVec&,bool,bool);
          void   Evaluate(cVector&, cVector&, vector<cVector> *vec = 0);
          void   EvaluateLFM(cVector&, cVector&, vector<cVector> *vec = 0);
          void   Init(int, cMatrix &, cMatrix &, cMatrix &, double, double);
        double   Eval(cVector, int, eProbLikeType);
          void   EvalVel(int, int, int, int, cVector, cMatrix, cMatrix, cMatrix &);
          void   Mut(int, double, double, double, cMatrix &);
          void   UpdatePos(int, double, double, cMatrix &, cMatrix &);
          void   MaxLikelihood(cVector &, cHIERKRG*);
          void   MaxLikelihoodd(cVector &, cHIERKRG*);
        double   Likelihood(vector<cVector> &, cVector &, int);
        double   Likelihoodd(vector<cVector> &, cVector &, int);

        void     MuSur(double &, int, cMatrix &);
        void     MulSur(double &, int, cMatrix &);
        void     SigmaSqrlSur(double &, double &, int, cMatrix &);
        void     MudSur(double &, int, cMatrix &, cVector);
        void     SigmaSqrdSur(double &, double &, int, cMatrix &, cVector);
      double     SSqrSur(cVector, cVector &, cVector &, int);
      double     SSqrlSur(cVector, cVector &, cVector &, int);
      double     PredictionSur(cVector &, cVector &, cVector &, int);
      double     PredictionSurl(cVector &, cVector &, int);
        void     CalcCorrMat(vector<cVector>, cVector, cMatrix &, cMatrix &);
        void     EvalWeights(cVector &, cVector &, int);
        void     EvalWeightsl(cVector &, cVector &, int);

          double   EvalExpImp(cVector &,double);
          double   EvalLCB(cVector&);
          double   EvalProbImp(cVector &,double);

          double   EvalVFLCB(cVector&);
          double   EvalVFExpImp(cVector &,double,int&);
          double   EvalVFProbImp(cVector &,double,int&);

          double   EI1(cVector &,double);
          double   EI2(cVector &,double);
          double   LCB1(cVector&);
          double   LCB2(cVector&);
          double   PI1(cVector&, double);
          double   PI2(cVector &,double);

          double   EvalInfPen(cVector &x,int o,double tol = 1e-6){ return EvalInfPen(x, o, 2); }                // LEO
          double   EvalProbFeas(cVector&x,int o, double tol = 1e-6){ return EvalProbFeas(x, o, 2); }               // LEO
          double   EvalProbFeasTutum(cVector&x,int o, double tol = 1e-6){ return EvalProbFeasTutum(x, o, 2); }     // LEO
          double   EvalProbFeasBagheri(cVector&x,int o, double tol = 1e-6){ return EvalProbFeasBagheri(x, o, 2); } // LEO
          double   EvalProbFeasSohst(cVector&x,int o, double tol = 1e-6){ return EvalProbFeasSohst(x, o, 2); }     // LEO

          double   EvalConstraintPF(cVector &, int, int, double tol = 1e-6){ }
          double   EvalInfPen(cVector&,int, int,double tol = 1e-6);           // LEO
          double   EvalProbFeas(cVector&,int, int , double tol = 1e-6);         // LEO
          double   EvalProbFeasTutum(cVector&,int, int,double tol = 1e-6);    // LEO
          double   EvalProbFeasBagheri(cVector&,int, int,double tol = 1e-6);  // LEO
          double   EvalProbFeasSohst(cVector&,int, int,double tol = 1e-6);    // LEO
};

#endif
