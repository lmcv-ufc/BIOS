// -------------------------------------------------------------------------
// surr.h - file containing the definition of the cSURR class.
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
// The class cSURR contains data and methods relevant to create a surrogate
// model using Radial Basis Functions (RBF) or the KRGing Model (KRG).
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
//   sigtype  - type of calculation adopted to evaluate the sigmas
//              of RBF functions                                     (in)
//   xs   -  array of sample points                                  (in)
//   ys   -  array with responses at each sample point               (in)
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

#ifndef _SURR_H
#define _SURR_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "group.h"

using namespace std;

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;
class cRBF;

// -------------------------------------------------------------------------
// Infill Methods:
//
typedef enum
{
  GlobalBest,
  DensityFunc,
  ExpectedImprovement,
} eInfillType;

// -------------------------------------------------------------------------
// Sigma evaluation methods:
//
typedef enum
{
  NAKAYAMA,
  KITAYAMA,
  KITAYAMAUNIFORM,
  ADAPTIVE_SCALING,
  LOOCV,
  kFOLDCV
} eSigmaType;

// -------------------------------------------------------------------------
// Definition of sample data struct:
//
typedef std::vector<cVector> cVectorVec;

struct sSampData
{
  int          NumVar;       // Number of variables
  int          NumOut;       // Number of output functions
  int          NumSample;    // Number of samples
  int          NumValSample; // Number of validations samples
  cVectorVec   SampleX;      // Sample points (NumSample x NumVar)
  cVectorVec   SampleY;      // Sample responses (NumSample x NumOut)
  cVectorVec   ValSampleX;   // Validation Sample points (NumSample x NumVar)
  cVectorVec   ValSampleY;   // Validation Sample responses (NumSample x NumOut)

  int          NumApproxOut; // Number of approximated outputs.
  bool*        ApproxOut;    // Vector of approximated output flags.

  cVector      Xlow;         // X lower bound
  cVector      Xupp;         // X upper bound
  cVector      Ylow;         // Y lower bound
  cVector      Yupp;         // Y upper bound
  cVectorVec   ExactCons;    // Sample high-fidelity constraints. 
  cVectorVec   ExactFobj;    // Sample high-fidelity objective function.
 

  void Normalize(cVector&) const;
  void UndoNormalization(cVector&) const;
  void EvalYbounds(void);
};
istream& operator>> (istream&,sSampData&);

// -------------------------------------------------------------------------
// Definition of SURR class:
//
class cSURR
{
 protected:
  sSampData         sdata;
  cMatrix*          Hmat;
  cMatrix*          Hmatdc;

  vector<cVector>   NRMSE;
  vector<double>    Prediction;
  bool              c = true;
  static eSigmaType Sigtype;

 public:
          cVector  SigmaSqr;      // Stores the SigmaSqr; used on statistical infill methods
          cVector  Mu;            // Stores the Mean; used on statistical infill methods
          double   WEI;           // Stores the w factor, used for weighted Expected Improvement (WEI)
          double   Beta;          // Stores the beta factor, used for Lower Confidence Bound (LCB)


                   cSURR(void);
                   cSURR(int, int, int, vector<cVector>, vector<double>);
  virtual          ~cSURR(void);
   static void     ReadSigtype(void);
          void     Read(ifstream&, int&, int&, int&, cVector&, cVector&, vector<cVector> &, vector<cVector> &, int&, int&, vector<cVector> &, vector<cVector> &);
          void     GetData(int&, int&, int&, int &);
          void     GetApproxOut(bool* &ap)                { ap = sdata.ApproxOut;      }
          void     GetSampleX(vector<cVector> &sx)        { sx = sdata.SampleX;        }
          void     GetSampleY(vector<cVector> &sy)        { sy = sdata.SampleY;        }
          void     GetValSampleX(vector<cVector> &vsx)    { vsx = sdata.ValSampleX;    }
          void     GetValSampleY(vector<cVector> &vsy)    { vsy = sdata.ValSampleY;    }
          void     GetExactConst(vector<cVector> &cy)     { cy = sdata.ExactCons;      }
		  void     GetExactFobj(vector<cVector> &fobjsex) { fobjsex = sdata.ExactFobj; }
          void     GetNumSample(int &ns)                  { ns  = sdata.NumSample;     }
          void     GetNumOut(int &no)                     { no  = sdata.NumOut;        }

		  void     SetWEI(double w)  { WEI = w;       }
          void     SetBeta(double b) { Beta = b;      }

 static eSigmaType GetSigtype(void) { return Sigtype; }
		  bool     IsInSample(cVector &, double);

 const sSampData&  GetSampData(void) { return sdata; }
 virtual  void     Evaluate(cVector&, cVector&, vector<cVector> *vec = 0);
 virtual  double   EvalExpImp(cVector &,double);
 virtual  double   EvalProbImp(cVector &,double);
 virtual  double   EvalLCB(cVector &);
 virtual  double   Eval(cVector, int);
 virtual  double   EvalConstraintPF(cVector &, int, double);
 virtual  double   EvalInfPen(cVector &, int, double tol = 1e-6){ }
 virtual  double   EvalProbFeas(cVector &, int, double tol = 1e-6){ }
 virtual  double   EvalProbFeasTutum(cVector &, int, double tol = 1e-6){ }
 virtual  double   EvalProbFeasBagheri(cVector &, int, double tol = 1e-6){ }
 virtual  double   SSqrSur(cVector, int);
          void     InfillCriteria(eInfillType, cVector &, cRBF*);
};

#endif
