// -------------------------------------------------------------------------
// KRG.cpp - implementation of cKRG class.
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
// Created:      28-May-2019    Marina Alves Maia
//
// Modified:     23-Dec-2019    Marina Alves Maia
//                              Code refactoring
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//#include <bits/stdc++.h>

#ifdef _OMP_
#include "omp.h"
#endif

#include "krg.h"
#include "surr.h"
#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "matvec.h"
#include "gblvar.h"

#include "problike.h"
#include "penalty.h"

using namespace std;





// -------------------------------------------------------------------------
// Public methods:
//

// ================================= cKRG ==================================

cKRG :: cKRG(void)
{
    SubAlgType     = new cStandardPSO; // Method for maximization of likelihood
    SubPop         = 250;              // Population size (maximization of likelihood)
    SubMaxGen      = 100;              // Number of iterations (maximization of likelihood)
    SubTolViol     = 1e-5;             // Tolerance for violation (maximization of likelihood)
    SubMutProb     = 0.05;             // Mutation probability (maximization of likelihood)
    Hmat           = 0;                // Psi matrix
    Hmatdc         = 0;                // Inverse of Psi matrix
}

// ================================= ~cKRG =================================

cKRG :: ~cKRG()
{
}

// ============================ CreateModel ===============================

void cKRG :: CreateModel(const sSampData &sampdata, eCorrelationType corrtype,
                         double HyperParamLow, double HyperParamUpp)
{
  // Store samples points and values.

  sdata = sampdata;

    Mu.Resize(sdata.NumOut);
    SigmaSqr.Resize(sdata.NumOut);
    SSqr.resize(sdata.NumOut);

    // Store samples points and values.

    sdata = sampdata;

    CorrType = corrtype;
    HPlow = HyperParamLow;
    HPupp = HyperParamUpp;

    // Create virtual matrices

    if (Hmat == 0 && Hmatdc == 0)
    {
        Hmat = new cMatrix [sdata.NumOut];
        Hmatdc = new cMatrix [sdata.NumOut];
    }

    cVector bestlikelihood(sdata.NumOut);

    MaxLikelihood(bestlikelihood, this);        // Maximization of likelihood (definition of hyper-parameters)

    // EGO algorithm

    int bestind;
    cVector error(sdata.NumSample);
    double Errmin=10e30;

    cVectorVec SampleXbackup = sdata.SampleX;   //Creation of the samples vector (for variables and objective functions)
    cVectorVec SampleYbackup = sdata.SampleY;

    int NS = sdata.NumSample;

    for (int j = 0; j < NS; j++)                           //Index of the sample used as validation test on this iteration
    {
        sdata.SampleX = SampleXbackup;                     //Creation of the samples vector (for variables and objective functions)
        sdata.SampleY = SampleYbackup;

        sdata.SampleX.erase (sdata.SampleX.begin()+j);     //Erasing the sample used as validation
        sdata.SampleY.erase (sdata.SampleY.begin()+j);

        sdata.NumSample = sdata.SampleY.size();

        if (Weight.size() > 0) Weight.clear();

        for(int m = 0; m < sdata.NumOut; m++)
        {
            cVector tempmle1(sdata.NumSample);
            EvalWeights(tempmle1, BestTheta[m], m);
        }

        cVector yi(sdata.NumOut);
        Evaluate(SampleXbackup[j], yi);                                //Evaluating the validation sample

        for (int nt = 0; nt < sdata.NumOut; nt++)
        {
          SSqr[nt] = SSqrSur(SampleXbackup[j], BestTheta[nt], nt);
          double diff = SampleYbackup[j][nt]-yi[nt];
          double nerr = diff/sqrt(SSqr[nt]);
          if(isnan(nerr))
          {
              cout << "Could not compute error for validation (fit Kriging model)" << endl;
              exit(0);
          }
          error[j] += nerr;
        }

        if (abs(error[j]) < Errmin)
        {
            bestind = j;
            Errmin = abs(error[j]);
        }
      }

    if (abs(error.Max()) > 3.0|| abs(error.Min()) > 3.0)
    {
        cout << "\nWARNING!" << endl;
        cout << "The maximum error is greater than the allowed threshold (Jones et al., 1998). You can try to run again, changing the bounds for the hyperparameters." << endl;
        cout << "If a deterministic sampling technique is being used, you could try to change the sampling method.\n" << endl;
    }

      // Selecting the sample set with the lowest error

     sdata.SampleX = SampleXbackup;
     sdata.SampleY = SampleYbackup;

     sdata.NumSample = sdata.SampleX.size();

     cout << "\n Pos EGO " << endl;

     if (Weight.size() > 0) Weight.clear();

     for(int m = 0; m < sdata.NumOut; m++)
     {
         cVector tempmle1(sdata.NumSample);
         EvalWeights(tempmle1, BestTheta[m], m);

         cout << "Sample taken out: " << bestind << endl;
         cout << "Out " << m << endl;
         cout << "Mu: " << Mu[m] << endl;
         cout << "SigmaSqr: " << SigmaSqr[m] << endl;
     }

     if (SSqr.size() > 0) SSqr.clear();
}

// ============================ UpdateModel ===============================

void cKRG :: UpdateModel(eCorrelationType corrtype, cVectorVec &sx, cVectorVec &sy)
{
    int nn = sx.size();
    sdata.NumSample += nn;

    // Store sample points and values.

    for (int i = 0; i < nn; i++)
    {
      sdata.SampleX.push_back(sx[i]);
      sdata.SampleY.push_back(sy[i]);
    }

    // Train the surrogate model.

    if (BestTheta.size() > 0) BestTheta.clear();

    cVector bestlikelihood(sdata.NumOut);

    MaxLikelihood(bestlikelihood, this);

    if (Weight.size() > 0) Weight.clear();

    for(int m = 0; m < sdata.NumOut; m++)
    {
       cVector tempmle1(sdata.NumSample);
       EvalWeights(tempmle1, BestTheta[m], m);
    }

    if (SSqr.size() > 0) SSqr.clear();
}

// ================================ Evaluate ===============================

void cKRG :: Evaluate(cVector &x, cVector &y, vector<cVector> *sampx)
{
    int flag = 0;

    if (sampx == NULL) {flag = 0;}            // No sample has been specified
    else  {flag = 1;}                         // Otherwise

   // Evaluate the KRGing output.

   cVector thetatemp(sdata.NumVar);

   for (int i = 0; i < sdata.NumOut; i++)
   {
       thetatemp = BestTheta[i];

       if (flag == 1)                        // If sample is specified (Cross validation)
       {
           for (int j = 0; j < sdata.NumVar; j++) thetatemp[j] = pow(10, thetatemp[j]);
           double mu, sigmasqrsur;
           CalcCorrMat((*sampx), thetatemp, U);
           cMatrix Utemp = U;
           MuSur(mu, i, Utemp);
           SigmaSqrSur(mu, sigmasqrsur, i, Utemp);
           thetatemp = BestTheta[i];
       }

       y[i] = PredictionSur(x, thetatemp, i);
   }
}

// ========================== MuSur ==============================

void cKRG :: MuSur(double &mu, int out, cMatrix &Uaux)
{
    cVector v1(sdata.NumSample);
    cVector sy(sdata.NumSample);

    cMatrix Utemp;
    Utemp.Resize(sdata.NumSample, sdata.NumSample);

    Utemp = Uaux;

    for (int j = 0; j < sdata.NumSample; j++)                                // Vector of ones
    {
        v1[j] = 1.0;
    }

    for (int j = 0; j < sdata.NumSample; j++) sy[j] = sdata.SampleY[j][out];

    cVector mi11 = sy;
    Utemp.SolveLU(mi11);

    cVector mi21 = v1;
    Utemp.SolveLU(mi21);

    double mi1, mi2;
    mi1= v1*mi11;
    mi2 = v1*mi21;
    mu = mi1/mi2;

    if (isnan(mi11[0]))
    {
        cout << "Evaluation of Mu - cKRG :: MuSur( ) - is returning a NaN value. Please, check implementation." << endl;
        exit(0);
    }
}

// ========================== Likelihood ==============================

double cKRG :: Likelihood(vector<cVector> &sampx, cVector &Theta, int out)
{
    cVector teta(sdata.NumVar);
    for (int i = 0; i < sdata.NumVar; i++) teta[i] = pow(10, Theta[i]);

   double mu, sigmasqrsur, NegLnLike;

   cMatrix Uaux;

   CalcCorrMat(sampx, teta, Uaux);
   MuSur(mu, out, Uaux);
   SigmaSqrSur(mu, sigmasqrsur, out, Uaux);

   double LnDetPsi = 0.0;

   for (int i = 0; i < sdata.NumSample; i++) LnDetPsi += log(abs(Uaux[i][i]));

   double nsamp = sdata.NumSample;

   NegLnLike =  -1*(-(nsamp/2)*log(sigmasqrsur) - 0.5*LnDetPsi);

   return NegLnLike;
}

// ========================== Maximize Likelihood ==============================

void cKRG :: MaxLikelihood(cVector &bestParticlevec, cKRG *SurMod)
{

    cOptAlgorithm *alg = SubAlgType;
    cPenalty *pen = new cPenStatic;

    eSolType SubSolType = SOL_DBL_VEC;

    alg -> SetSolType(SubSolType);
    alg -> SetPopSize(SubPop);
    alg -> SetMaxGen(SubMaxGen);
    alg -> SetOptNum(1);
    alg -> SetFeedback(false);
    alg -> SetTolViol(SubTolViol);
    alg -> SetPenFunction(pen);
    alg -> SetMutProb(SubMutProb);

    int no;
    SurMod->GetNumOut(no);

    cout << "\nMaximization of the Likelihood Estimator . . . . . . . ." << endl;
    for (int i = 0; i < no; i++)
    {
        cProblem* probbest = new cProbLikelihood(SurMod, i);

        alg -> SetProblem(probbest);
        alg -> Init( );

        alg -> Solver();
        cOptSolution* best = alg->GetBest();
        bestParticlevec = best->GetObjFunc(0);

        BestTheta.push_back(best->GetBestVec());
    }
    cout << "Best thetas: " << endl;
    for (int i = 0; i < sdata.NumOut; i++)
    {
        cout << "Out: " << i+1 << "  Theta: ";
        BestTheta[i].Print();
    }
}

// ============================== GetSigma ================================

void cKRG :: EvalWeights(cVector &mle1, cVector &theta, int out)
{
    mle1.Resize(sdata.NumSample);

    cVector temp(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) temp[i] = pow(10, theta[i]);

    cMatrix Utemp;

    CalcCorrMat(sdata.SampleX, temp, Utemp);

    Hmatdc[out].Resize(sdata.NumSample, sdata.NumSample);

    Hmatdc[out] = Utemp;

    double mu, sigmasqrsur;
    MuSur(mu, out, Utemp);
    SigmaSqrSur(mu, sigmasqrsur, out, Utemp);
    Mu[out] = mu;
    SigmaSqr[out] = sigmasqrsur;

    cVector fi2(sdata.NumSample);

    for (int j = 0; j < sdata.NumSample; j++) fi2[j] = sdata.SampleY[j][out] - mu;

    mle1 = fi2;

    Utemp.SolveLU(mle1);

    Weight.push_back(mle1);
}


// ========================== PredictionSur ==============================

double cKRG :: PredictionSur(cVector &x, cVector &Theta, int out)
{
    cVector fi(sdata.NumSample);
    cVector teta(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) teta[i] = pow(10, Theta[i]);

    double p = 1.99;
    double d, dist;

    if(CorrType == GAUSS)
    {
        for (int i = 0; i < sdata.NumSample; i++)
        {
            double d = 0.0;

            for (int k = 0; k < sdata.NumVar; k++)
          d = d + teta[k]*(sdata.SampleX[i][k]-x[k])*(sdata.SampleX[i][k]-x[k]);

            fi[i] = exp(-d);
        }
    }
    else
    {
        for (int i = 0; i < sdata.NumSample; i++)
        {
            d = 1.0;

            for (int k = 0; k < sdata.NumVar; k++)
            {
                dist = abs(sdata.SampleX[i][k]-x[k]);
                d = d*exp(-(dist*sqrt(5.0)/teta[k]))*(1.0 + (dist*sqrt(5.0))/teta[k] + (5.0*dist*dist)/(3*teta[k]*teta[k]));
            }

            if (d < 0)
            {
                cout << "Computation of the Matérn basis function returned an invalid |x_i - x|" << endl;
                cout << "|x_i - x| was set to 0.0 (zero)." << endl;
                d = 0;
            }
            else if (d > 1.0)
            {
                cout << "Computation of the Matérn basis function returned an invalid |x_i - x|" << endl;
                cout << "|x_i - x| was set to 1.0 (one)." << endl;
                d = 1.0;
            }
            fi[i] = d;
        }
    }

    cVector mle1(sdata.NumSample);

    mle1 = Weight[out];

    double pred = Mu[out] + fi*mle1;

    return pred;
 }

// ========================== SSqrSur ==============================

double cKRG :: SSqrSur(cVector x, cVector &theta, int out)
{
    cVector fi(sdata.NumSample);
    cVector ones(sdata.NumSample);

    cVector teta(sdata.NumVar);

    for (int i = 0; i< sdata.NumVar; i++) teta[i] = pow(10, theta[i]);

    double sigmasqrsur;

    cMatrix Hm(sdata.NumSample, sdata.NumSample);

    Hm = Hmatdc[out];

    sigmasqrsur = SigmaSqr[out];

    double d, dist;

    if(CorrType == GAUSS)
    {
    for (int i = 0; i < sdata.NumSample; i++)
    {
        double d = 0.0;

        for (int k = 0; k < sdata.NumVar; k++) 
          d = d + teta[k]*(sdata.SampleX[i][k]-x[k])*(sdata.SampleX[i][k]-x[k]);

        fi[i] = exp(-d);

        ones[i] = 1.0;

    }
    }
    else
    {
        for (int i = 0; i < sdata.NumSample; i++)
        {
            d = 1.0;

            for (int k = 0; k < sdata.NumVar; k++)
            {
                dist = abs(sdata.SampleX[i][k]-x[k]);

                d = d*exp(-(dist*sqrt(5.0)/teta[k]))*(1.0 + (dist*sqrt(5.0))/teta[k] + (5.0*dist*dist)/(3*teta[k]*teta[k]));
            }

            if (d < 0 || d > 1.0)
            {
                cout << "Computation of the Matérn basis function returned an invalid |x_i - x|" << endl;
                cout << "|x_i - x| was set to 0 (zero)." << endl;
                d = 0;
            }

            fi[i] = d;

            ones[i] = 1.0;
        }
    }

    cVector s2c1 = fi;
    cVector onestemp = ones;
    Hm.SolveLU(s2c1);

    Hm.SolveLU(onestemp);

    double s2c = fi*s2c1;

    double ssqr = sigmasqrsur*(1.0 - s2c + (pow((1-ones*s2c1),2)/(ones*onestemp)));


    if (ssqr < 0)
    {
        cout << "Computation of the variance returned an invalid value (lower than 0.0): " << ssqr << endl;
        cout << "Variance was set to 0.0 (zero)." << endl;
        ssqr = 0;
    }

    return ssqr;
}


// ========================== SigmaSqrSur ==============================

void cKRG :: SigmaSqrSur(double &mu, double &sigmasqrsur, int out, cMatrix &Uaux)
{
    cVector sy(sdata.NumSample);
    cMatrix Psitemp(sdata.NumSample, sdata.NumSample);

    Psitemp = Uaux;

    for (int j = 0; j < sdata.NumSample; j++)
    {
        sy[j] = sdata.SampleY[j][out];
    }

    cVector sd1(sdata.NumSample);

    for (int j = 0; j < sdata.NumSample; j++)
    {
         sd1[j] = sy[j] - mu;
    }

     cVector var11 = sd1;
     Psitemp.SolveLU(var11);
     sigmasqrsur = (sd1*var11)/sdata.NumSample;
}

// ========================== Correlation Matrix ==============================

void cKRG :: CalcCorrMat(vector<cVector> sampx, cVector theta, cMatrix &Uaux)
{
    int nsamp = sampx.size();

    cMatrix Psi(nsamp, nsamp);
    Psi.Zero();

    // Evaluate the correlation matrix [Psi].

    double d, dist;

    if (CorrType == GAUSS)
    {

        for (int i = 0; i < nsamp; i++)
        {
            for (int j = i+1; j < nsamp; j++)
            {
                d = 0.0;

                for (int k = 0; k < sdata.NumVar; k++)  d = d + theta[k]*(sampx[i][k]-sampx[j][k])*(sampx[i][k]-sampx[j][k]);

                Psi[i][j] = exp(-d);
            }
        }
    }
    else
    {
        for (int i = 0; i < nsamp; i++)
        {
            for (int j = i+1; j < nsamp; j++)
            {
                d = 1.0;

                for (int k = 0; k < sdata.NumVar; k++)
                {
                    dist = abs(sampx[i][k]-sampx[j][k]);

                    d = d*exp(-(dist*sqrt(5.0)/theta[k]))*(1.0 + (dist*sqrt(5.0))/theta[k] + 5.0*dist*dist/(3*(theta[k]*theta[k])));
                }

                if (d < 0 || d > 1.0)
                {
                    d = 0;
                }

                Psi[i][j] = d;
            }
        }
    }

    cMatrix Psit(nsamp, nsamp);
    Psi.Transp(Psit);

    cMatrix eye(nsamp, nsamp);
    eye.Zero();
    for (int i = 0; i < nsamp; i++)
      for (int j = 0; j < nsamp; j++)
          {
              if (i == j)  eye[i][j] = 1.0;
          }

    cMatrix eps(nsamp, nsamp);
    eps.Zero();
    for (int i = 0; i < nsamp; i++)
    {
      for (int j = 0; j < nsamp; j++)
          {
              if (i == j)  eps[i][j] = 1.0e-8;
          }
    }

    Psi = Psi + Psit + eps +eye;

    Uaux.Resize(nsamp, nsamp);
    Uaux.Zero();
    Uaux = Psi;
    Uaux.DecompLU();
}


// ========================== EvalExpImp ===================================

double cKRG :: EvalExpImp(cVector &x, double ybest)
{
  double ei;

  double pi = atan(1)*4;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrSur(x, BestTheta[0], 0);
  Pred[0] = ytemp[0];
  s[0] = sqrt(ssqr[0]);

  bool isin = false;
  if (s[0] < 0.01)
  {
      isin = IsInSample(x, 1e-5);
  }

      if (s[0] <= 0 || isin == true)
      {
          ei = 0.0;
      }
      else
      {
          double EI1 = ybest - Pred[0];
          double erf1 = erf(EI1/(sqrt(2*ssqr[0])));
          double EI2 = 0.5 + 0.50*erf1;
          double EI3 = s[0]*(1/(sqrt(2*pi)));
          double EI4 = exp(-(EI1*EI1)/(2*ssqr[0]));
          ei = WEI*(EI1*EI2) + (1 - WEI)*(EI3*EI4);
      }

  return(ei);
}

// ========================== EvalProbImp ===================================

double cKRG :: EvalProbImp(cVector &x, double ybest)
{
  double PoI;

  double pi = atan(1)*4;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrSur(x, BestTheta[0], 0);
  Pred[0] = ytemp[0];
  s[0] = sqrt(ssqr[0]);

  bool isin = false;
  if (s[0] < 0.01)
  {
      isin = IsInSample(x, 1e-5);
  }

      if (s[0] <= 0 || isin == true)
      {
          PoI = 0.0;
      }
      else
      {
          double term = ybest - Pred[0];
          PoI = erf(term/(sqrt(2*ssqr[0])));
      }

  return(PoI);
}


// ========================== EvalProbImp ===================================

double cKRG :: EvalLCB(cVector &x)
{
  double LCB;

  double pi = atan(1)*4;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrSur(x, BestTheta[0], 0);
  Pred[0] = ytemp[0];
  s[0] = sqrt(ssqr[0]);

  bool isin = false;
  if (s[0] < 0.01)
  {
  isin = IsInSample(x, 1e-5);
  }

  LCB = Pred[0] - Beta*s[0];

  return(LCB);
}

// ========================== EvalInfPen ==============================

double cKRG :: EvalInfPen(cVector &x, int out, double tolviol)
{
  double pf;
  double PF = 0;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  Pred = ytemp[out];

  if (Pred > tolviol)
  {
      PF = 0;
  }
  else
  {
      PF = 1;
  }

  pf = PF;

  return pf;
}

// ========================== EvalProbFeas ==============================

double cKRG :: EvalProbFeas(cVector &x, int out, double tolviol)
{
  double pf;
  double PF = 0;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  ssqr = SSqrSur(x, BestTheta[out], out);

  double s;

  if(ssqr <= 0)
  {
      s = 0;
  }
  else
  {
    s  = sqrt(ssqr);
  }

  Pred = ytemp[out];

  bool isin = false;
  if (s < 10)
  {
  isin = IsInSample(x, 1e-5);
  }

  if (s <= tolviol || isin == true)
  {
      PF = 0;
  }
  else
  {
      double erfunc = erf((0 - Pred)/s);
      PF = 0.5 + 0.5*erfunc;
  }

  pf = PF;

  return pf;

}

// ========================== EvalProbFeasTutum ==============================

double cKRG :: EvalProbFeasTutum(cVector &x, int out, double tolviol)
{
    double pf;
    double PF = 0;

    double ssqr(sdata.NumOut);

    cVector ytemp(sdata.NumOut);
    Evaluate(x, ytemp);

    ssqr = SSqrSur(x, BestTheta[out], out);

    double s = sqrt(ssqr);

    if (ssqr <= 0)
    {
        s = 0;
    }
    else
    {
        s = sqrt(ssqr);
    }

    double Pred = ytemp[out];

    bool isin = false;

    if (s <= 0 || isin == true)
    {
        PF = 0;
    }
    else
    {
        double erfunc = erf((0 - Pred)/s);

        if (erfunc >= 1){
            PF = 0.5 + 0.5*erfunc;
        }
        else if(erfunc > 0){
            PF = 2 - erfunc;
        }
        else{
            PF = 0;
        }
    }
    pf = PF;

    return(pf);
}

// ========================== EvalProbFeasBagheri ==============================

double cKRG :: EvalProbFeasBagheri(cVector &x, int out, double tolviol)
{
  double pf;
  double PF = 0;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  ssqr = SSqrSur(x, BestTheta[out], out);

  double s;

  if(ssqr <= 0)
  {
      s = 0;
  }
  else
  {
    s  = sqrt(ssqr);
  }

  Pred = ytemp[out];

  bool isin = false;
  if (s < 10)
  {
  isin = IsInSample(x, 1e-5);
  }

  if (s <= tolviol || isin == true)
  {
      PF = 0;
  }
  else
  {
      double erfunc = erf((0 - Pred)/s);

      PF = 2*(0.5 + 0.5*erfunc);

      if (PF >= 1)
      {
          PF = 1;
      }
  }

  pf = PF;

  return pf;

}


// ============================== GetSigma ================================

void cKRG :: GetSigma(cVector &s)
{

}

// ============================== Init ================================

void cKRG :: Init(int i, cMatrix &xx, cMatrix &v, cMatrix &xp, double inf, double sup)
{
   // cout << "Inicializou indiv: " << i << endl;

    for(int j = 0; j < sdata.NumVar; j++)
    {
        // Initialize the particle position

        xx[i][j] = Utl::RandDouble(inf,sup);

        // Initialize the particle velocity

        v[i][j] = (sup - inf) * Utl::RandDouble(-1.0,1.0);

        // Assign the particle position to the best position so far

        xp[i][j] = xx[i][j];
    }
}

// ============================== Evaluate ================================

double cKRG :: Eval(cVector temp, int m)
{
    double fobj = Likelihood(sdata.SampleX, temp, m);

    return fobj;
}
