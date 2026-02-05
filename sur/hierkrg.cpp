// -------------------------------------------------------------------------
// cokrg.cpp - implementation of cCOKRG class.
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
#include <bits/stdc++.h>

#ifdef _OMP_
#include "omp.h"
#endif

#include "hierkrg.h"
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

// ================================= cCOKRG ==================================

cHIERKRG :: cHIERKRG(void)
{
    SubAlgType     = new cStandardPSO;  // Method for maximization of likelihood
    SubPop         = 50;               // Population size (maximization of likelihood)
    SubMaxGen      = 200;              // Number of iterations (maximization of likelihood)
    SubStallGen    = 50;
    SubDifType     = Loc2Best;
    SubTopology    = RING_TOPOLOGY;
    SubTolViol     = 1e-5;
    SubMutProb     = 0.05;
    Hmat = 0;      // C matrix
    Hmatdc = 0;    // Inverse of C matrix

    Psil   = 0;    // LF model correlation matrix
    Psildc = 0;    // Inverse of Psil
    Psid   = 0;    // Difference model correlation matrix
    Psiddc = 0;    // Inverse of Psid

    fVec   = 0;    // Matrix of difference vectors

    Update = 0;
}

// ================================= ~cCOKRG =================================

cHIERKRG :: ~cHIERKRG()
{
}

// ============================ CreateModel ===============================

void cHIERKRG :: CreateModel(const sSampData &sampdata, eCorrelationType corrtype,
                         double HyperParamLow, double HyperParamUpp)
{
    Update = 0;

    // Store samples points and values.

    sdata = sampdata;

    /*sdata.NumSampleLF = 7;
    sdata.NumSample = 4;

    cVector newx(1);
    cVector newylf(1);
    cVector newyhf(1);
    newx[0] = 1.0;
    newyhf[0] = 15.829731945974109;
    newylf[0] = 7.914865972987055;

    sdata.SampleX.push_back(newx);
    sdata.SampleXLF.push_back(newx);
    sdata.SampleY.push_back(newyhf);
    sdata.SampleYLF.push_back(newylf);*/

    Mu.Resize(sdata.NumOut);
    SigmaSqr.Resize(sdata.NumOut);
    SSqr.resize(sdata.NumOut);

    Mul.Resize(sdata.NumOut);
    SigmaSqrl.Resize(sdata.NumOut);

    Mud.Resize(sdata.NumOut);
    SigmaSqrd.Resize(sdata.NumOut);

    CorrType = corrtype;
    HPlow = HyperParamLow;
    HPupp = HyperParamUpp;

    /*cout << "HPlow: "  << HPlow << endl;
    cout << "HPupp: "  << HPupp << endl;
    cout << "rhoLow: " << rhoLow << endl;
    cout << "rhoUpp: " << rhoUpp << endl;

    cout << "NumSample = " << NumSample << endl;
    for (int i = 0; i < NumSample; i++) SampleX[i].Print( );
    //cout << "NumOut = " << NumOut << endl;
    for (int i = 0; i < NumSample; i++) SampleY[i].Print( );
    for (int i = 0; i < NumSample; i++) ExactFobj[i].Print( );*/

    // Get X bounds

  //  Xlow.Resize(nv);
 //   Xupp.Resize(nv);
 //   Xlow = low;
 //   Xupp = upp;

    // Get Y bounds.
//    sdata.EvalYbounds( );

    // Create virtual matrices

    if (Hmat == 0 && Hmatdc == 0 && Psil == 0 && Psildc == 0 && Psid == 0 && Psiddc == 0 && fVec == 0)
    {
        Hmat   = new cMatrix [sdata.NumOut];
        Hmatdc = new cMatrix [sdata.NumOut];
        Psil   = new cMatrix [sdata.NumOut];
        Psildc = new cMatrix [sdata.NumOut];
        Psid   = new cMatrix [sdata.NumOut];
        Psiddc = new cMatrix [sdata.NumOut];
        fVec   = new cVector [sdata.NumOut];
    }

    cVector bestlikelihood(sdata.NumOut);
    MaxLikelihood(bestlikelihood, this);

    if (Weightl.size() > 0) Weightl.clear();
    if (Weightl1.size() > 0) Weightl1.clear();

    for(int m = 0; m < sdata.NumOut; m++)
    {
      cVector tempmle1(sdata.NumSampleLF);
      EvalWeightsl(tempmle1, BestTheta[m], m);

      cout << "Out: " << m+1;
      cout << "  Mul: " << Mul[m];
      cout << "  SigmaSqrl: " << SigmaSqrl[m] << endl;
    }
    cVector bestlikelihoodd(sdata.NumOut);

    MaxLikelihoodd(bestlikelihoodd, this);

    if (Weight.size() > 0) Weight.clear();

    for(int m = 0; m < sdata.NumOut; m++)
    {
      // BestThetad[m] = -3.0; // APAGAR
      cVector tempmle1(sdata.NumSample);
      EvalWeights(tempmle1, BestThetad[m], m);

      cout << "Out: " << m+1;
      cout << "  SigmaSqr: " << SigmaSqr[m];
      cout << "  Mu: " << Mu[m] << endl;
    }

    if (SSqr.size() > 0) SSqr.clear();
    if (SSqrl.size() > 0) SSqrl.clear();
/*
    cVector xtest(1);
    cVector ytest(1);
    xtest[0] = 0.9;

    Evaluate(xtest, ytest);

    cout << "\nTEST EVALUATE" << endl;
    cout << "x = " << xtest[0] << "    y = " << ytest[0] << endl;*/

    /*
    cVector xtest(1);
    double s2 = 0;
    xtest[0] = 0.9;

    s2 = SSqrlSur(xtest, BestTheta[0], BestThetad[0], 0);

    cout << "\nTEST SSQR" << endl;
    cout << "x = " << xtest[0] << "    s = " << pow(s2, 0.5) << endl;*/
}

// ============================ UpdateModel ===============================

void cHIERKRG :: UpdateModel(eCorrelationType corrtype, cVectorVec &sx, cVectorVec &sy, cVectorVec &syl, bool addhf, bool addlf)
{
    Update = 0;

    int nn = sx.size();

    // Store sample points and values.

    for (int k = 0; k < nn; k++)
    {
        if (addhf == 1)
        {
            sdata.NumSample += 1;

            for (int i = 0; i < nn; i++)
            {
              sdata.SampleX.push_back(sx[i]);
              sdata.SampleY.push_back(sy[i]);
            }
        }

        if (addlf == 1)
        {
            sdata.NumSampleLF += 1;

            for (int i = 0; i < nn; i++)
            {
              sdata.SampleXLF.push_back(sx[i]);
              sdata.SampleYLF.push_back(syl[i]);
            }
        }
    }

    /*for (int k = 0; k < sdata.NumSample; k++)
    {
        cout << "AMOSTRA " << k+1 << endl;
        cout << "X =  ";
        sdata.SampleX[k].Print( );
        cout << "y =  ";
        sdata.SampleY[k].Print( );
    }

    for (int k = 0; k < sdata.NumSampleLF; k++)
    {
        cout << "AMOSTRA " << k+1 << endl;
        cout << "X =  ";
        sdata.SampleXLF[k].Print( );
        cout << "y =  ";
        sdata.SampleYLF[k].Print( );
    }*/

    // Train the surrogate model.

    if (BestTheta.size() > 0) BestTheta.clear();

    cVector bestlikelihood(sdata.NumOut);

    MaxLikelihood(bestlikelihood, this);

    if (Weightl.size() > 0) Weightl.clear();
    if (Weightl1.size() > 0) Weightl1.clear();

    for(int m = 0; m < sdata.NumOut; m++)
    {
      cVector tempmle1(sdata.NumSampleLF);
      EvalWeightsl(tempmle1, BestTheta[m], m);

      cout << "Out: " << m+1;
      cout << "  Mul: " << Mul[m];
      cout << "  SigmaSqrl: " << SigmaSqrl[m] << endl;
    }

    if (SSqr.size() > 0) SSqr.clear();
    if (SSqrl.size() > 0) SSqrl.clear();

    if (BestThetad.size() > 0) BestThetad.clear();

    cVector bestlikelihoodd(sdata.NumOut);

    MaxLikelihoodd(bestlikelihoodd, this);

    if (Weight.size() > 0) Weight.clear();

    for(int m = 0; m < sdata.NumOut; m++)
    {
      //BestThetad[m] = -3.0; // APAGAR
      cVector tempmle1(sdata.NumSample);
      EvalWeights(tempmle1, BestThetad[m], m);

      cout << "Out: " << m+1;
      cout << "  SigmaSqr: " << SigmaSqr[m];
      cout << "  Mu: " << Mu[m] << endl;
    }

    if (SSqr.size() > 0) SSqr.clear();
    if (SSqrl.size() > 0) SSqrl.clear();
}

// ================================ Evaluate ===============================

void cHIERKRG :: EvaluateLFM(cVector &x, cVector &y, vector<cVector> *sampx)
{
    int flag = 0;

    if (sampx == NULL) {flag = 0;}           // No sample has been specified
    else  {flag = 1;}                         // Otherwise

   // Evaluate the Co-KRGing output - LF model.

   cVector thetatemp(sdata.NumVar);

   for (int i = 0; i < sdata.NumOut; i++)
   {
       thetatemp = BestTheta[i];

       if (flag == 1)             // If sample is specified (Cross validation)
       {
           for (int j = 0; j < sdata.NumVar; j++) thetatemp[j] = pow(10, thetatemp[j]);
           double mu, sigmasqrsur;
           CalcCorrMat((*sampx), thetatemp, Psi, U);
           cMatrix Utemp = U;
           MulSur(mu, i, Utemp);
           SigmaSqrlSur(mu, sigmasqrsur, i, Utemp);
           thetatemp = BestTheta[i];
       }

       y[i] = PredictionSurl(x, thetatemp, i);
   }
}

// ================================ Evaluate ===============================

void cHIERKRG :: Evaluate(cVector &x, cVector &y, vector<cVector> *sampx)
{
    int flag = 0;

    if (sampx == NULL) {flag = 0;}           // No sample has been specified
    else  {flag = 1;}                         // Otherwise

   // Evaluate the KRGing output - LF model.

   cVector thetaltemp(sdata.NumVar);
   cVector thetadtemp(sdata.NumVar);

   for (int i = 0; i < sdata.NumOut; i++)
   {
       thetaltemp = BestTheta[i];
       thetadtemp = BestThetad[i];

       y[i] = PredictionSur(x, thetaltemp, thetadtemp, i);
   }
}

// ========================== MuSur ==============================

void cHIERKRG :: MuSur(double &mu, int out, cMatrix &Uaux)
{
    cVector v1(sdata.NumSample + sdata.NumSampleLF);
    cVector sy(sdata.NumSample + sdata.NumSampleLF);

    cMatrix Utemp;
    Utemp.Resize((sdata.NumSample + sdata.NumSampleLF), (sdata.NumSample + sdata.NumSampleLF));

    Utemp = Uaux;

    for (int j = 0; j < (sdata.NumSample + sdata.NumSampleLF); j++)                                // Vector of ones
    {
        v1[j] = 1.0;
    }

    for (int j = 0; j < (sdata.NumSample + sdata.NumSampleLF); j++)
    {
        if (j < sdata.NumSampleLF)
        {
            sy[j] = sdata.SampleYLF[j][out];
        }
        else
        {
            sy[j] = sdata.SampleY[j - sdata.NumSampleLF][out];
        }
    }

    cVector mi11 = sy;
    Utemp.SolveLU(mi11);

    if (isnan(mi11[out]))
    {
        cout << "mu" << endl;
        cout << "FLAAAAAAAAAAG  ==================" << endl;
        exit(0);
    }


    cVector mi21 = v1;
    Utemp.SolveLU(mi21);

    double mi1, mi2;
    mi1= v1*mi11;
    mi2 = v1*mi21;
    mu = mi1/mi2;

//   cout << "Mu " << mu << "MI1: " << mi1 << " MI2: " << mi2 <<  endl;
}

// ========================== MuSur ==============================

void cHIERKRG :: MulSur(double &mu, int out, cMatrix &Uaux)
{
    cVector v1(sdata.NumSampleLF);
    cVector sy(sdata.NumSampleLF);

    cMatrix Utemp;
    Utemp.Resize(sdata.NumSampleLF, sdata.NumSampleLF);

    Utemp = Uaux;

    for (int j = 0; j < sdata.NumSampleLF; j++)                                // Vector of ones
    {
        v1[j] = 1.0;
    }

    for (int j = 0; j < sdata.NumSampleLF; j++) sy[j] = sdata.SampleYLF[j][out];

    cVector mi11 = sy;
    Utemp.SolveLU(mi11);

    if (isnan(mi11[out]))
    {
        cout << "mul" << endl;
        cout << "FLAAAAAAAAAAG  ==================" << endl;
        exit(0);
    }


    cVector mi21 = v1;
    Utemp.SolveLU(mi21);

    double mi1, mi2;
    mi1= v1*mi11;
    mi2 = v1*mi21;
    mu = mi1/mi2;

//   cout << "Mu " << mu << "MI1: " << mi1 << " MI2: " << mi2 <<  endl;
}

// ========================== MuSur ==============================

void cHIERKRG :: MudSur(double &mu, int out, cMatrix &Uaux, cVector faux)
{
    cVector v1(sdata.NumSample);
    cVector sy(sdata.NumSample);

    cMatrix Utemp;
    Utemp.Resize(sdata.NumSample, sdata.NumSample);

    Utemp = Uaux;

    for (int j = 0; j < sdata.NumSample; j++) sy[j] = sdata.SampleY[j][out];

    cVector mi11 = sy;
    Utemp.SolveLU(mi11);

    if (isnan(mi11[out]))
    {
        sy.Print( );
        Utemp.Print( );
        cout << "mud" << endl;
        cout << "FLAAAAAAAAAAG  ==================" << endl;
        exit(0);
    }


    cVector mi21 = faux;
    Utemp.SolveLU(mi21);

    double mi1, mi2;
    mi1 = faux*mi11;
    mi2 = faux*mi21;
    mu = mi1/mi2;

//   cout << "Mu " << mu << "MI1: " << mi1 << " MI2: " << mi2 <<  endl;
}

// ========================== Likelihood ==============================

double cHIERKRG :: Likelihood(vector<cVector> &sampx, cVector &Theta, int out)
{
  /*  cout << "Calc likelihood" << endl;
    cout << "Amostra: " << sampx[1][0] << endl;*/
    cVector teta(sdata.NumVar);
    for (int i = 0; i < sdata.NumVar; i++) teta[i] = pow(10, Theta[i]);

   double mu, sigmasqrsur, NegLnLike;

   cMatrix Uaux;
   cMatrix Psiaux;

   CalcCorrMat(sampx, teta, Psiaux, Uaux);
   MulSur(mu, out, Uaux);
   SigmaSqrlSur(mu, sigmasqrsur, out, Uaux);

   //cMatrix Utemp(NumSample, NumSample);
   //Utemp = Uaux;

   double LnDetPsi = 0.0;

   for (int i = 0; i < sdata.NumSampleLF; i++) LnDetPsi += log(abs(Uaux[i][i]));

   double nsamp = sdata.NumSampleLF;

   NegLnLike =  -1*(-(nsamp/2)*log(sigmasqrsur) - 0.5*LnDetPsi);

/*   cout << "\nLIKELIHOOD" << endl;
   cout << "Likelihood: " << NegLnLike << endl;
   cout << "Mu: " << mu << endl;
   cout << "Sigmasqr: " << sigmasqrsur << endl;*/

   return NegLnLike;
}

// ========================== Likelihood ==============================

double cHIERKRG :: Likelihoodd(vector<cVector> &sampx, cVector &x, int out)
{
  /*  cout << "Calc likelihood" << endl;
    cout << "Amostra: " << sampx[1][0] << endl;*/
    cVector Theta(sdata.NumVar);

    for (int i = 0; i < (sdata.NumVar); i++) Theta[i] = x[i];

    cVector teta(sdata.NumVar);
    for (int i = 0; i < sdata.NumVar; i++) teta[i] = pow(10, Theta[i]);

   double mu, sigmasqrsur, NegLnLike;

   cVector faux(sdata.NumSample);

   for (int i = 0; i < sdata.NumSample; i++)
   {
       cVector predlf(sdata.NumOut);
       EvaluateLFM(sdata.SampleX[i], predlf);
       faux[i] = predlf[out];
   }

   cMatrix Psiaux;
   cMatrix Uaux;

   CalcCorrMat(sampx, teta, Psiaux, Uaux);

   MudSur(mu, out, Uaux, faux);
   SigmaSqrdSur(mu, sigmasqrsur, out, Uaux, faux);

   //cMatrix Utemp(NumSample, NumSample);
   //Utemp = Uaux;

   double LnDetPsi = 0.0;

   for (int i = 0; i < sdata.NumSample; i++) LnDetPsi += log(abs(Uaux[i][i]));

   double nsamp = sdata.NumSample;

   NegLnLike =  -1*(-(nsamp/2)*log(sigmasqrsur) - 0.5*LnDetPsi);

/*   cout << "\nLIKELIHOOD" << endl;
   cout << "Likelihood: " << NegLnLike << endl;
   cout << "Mu: " << mu << endl;
   cout << "Sigmasqr: " << sigmasqrsur << endl;*/

   return NegLnLike;
}

// ========================== Maximize Likelihood ==============================

void cHIERKRG :: MaxLikelihood(cVector &bestParticlevec, cHIERKRG *SurMod)
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

    alg -> SetStallGen(SubStallGen);

    if (alg->GetType() == STANDARD_PSO){
        alg -> SetSwarmTopology(SubTopology);
    }
    else if(alg->GetType() == STANDARD_DE){
        alg -> SetDifType(SubDifType);
    }

    int no;
    SurMod->GetNumOut(no);

    //no = 1; // Gambiarra
    for (int i = 0; i < no; i++)
    {
        cProblem* probbest = new cProbLikelihood(SurMod, i, BASIC);
        //probbest -> SetApproxC( ApproxC );

        if (Update == 1)
        {
            int NumInitSol = 1;
            sInpSol* SolVecInit = new sInpSol [NumInitSol];

            for(int k = 0; k < NumInitSol; k++)
            {
              SolVecInit[k].type = SubSolType;
              SolVecInit[k].CodVar.Resize(probbest -> GetNumVar( ));

              // Read each solution variable.
              for(int j = 0; j < probbest -> GetNumVar( ); ++j)
                SolVecInit[k].CodVar[j] = BestThetaOld[i][j];

              cout << "Initial Theta: ";
              SolVecInit[k].CodVar.Print( );
            }

            alg -> SetInpSolVec(SolVecInit);
            alg -> SetNumInpSol(NumInitSol);

            BestThetaOld.clear( );
        }

        alg -> SetProblem(probbest);
        alg -> Init( );
        alg -> Solver();
        cOptSolution* best = alg->GetBest();
        bestParticlevec = best->GetObjFunc(0);

        BestTheta.push_back(best->GetBestVec());
        BestThetaOld = BestTheta;
    }

    cout << "\nMAXIMUM LIKELIHOOD ESTIMATORS" << endl;
    for (int i = 0; i < sdata.NumOut; i++)
    {
        cout << "Out: " << i+1 << "  BestTheta = ";
        BestTheta[i].Print();
    }
}

// ========================== Maximize Likelihood ==============================

void cHIERKRG :: MaxLikelihoodd(cVector &bestParticlevec, cHIERKRG *SurMod)
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

    alg -> SetStallGen(SubStallGen);

    if (alg->GetType() == STANDARD_PSO){
        alg -> SetSwarmTopology(SubTopology);
    }
    else if(alg->GetType() == STANDARD_DE){
        alg -> SetDifType(SubDifType);
    }

    int no;
    SurMod->GetNumOut(no);

    //no = 1; // Gambiarra
    for (int i = 0; i < no; i++)
    {
        cProblem* probbest = new cProbLikelihood(SurMod, i, DIFFMODELFIT);
        //probbest -> SetApproxC( ApproxC );

        if (Update == 1)
        {
            int NumInitSol = 1;
            sInpSol* SolVecInit = new sInpSol [NumInitSol];

            for(int k = 0; k < NumInitSol; k++)
            {
              SolVecInit[k].type = SubSolType;
              SolVecInit[k].CodVar.Resize(probbest -> GetNumVar( ));

              // Read each solution variable.
              for(int j = 0; j < probbest -> GetNumVar( ); ++j)
                SolVecInit[k].CodVar[j] = BestThetadOld[i][j];

              cout << "Initial Theta: ";
              SolVecInit[k].CodVar.Print( );
            }

            alg -> SetInpSolVec(SolVecInit);
            alg -> SetNumInpSol(NumInitSol);

            BestThetadOld.clear( );
        }

        alg -> SetProblem(probbest);
        alg -> Init( );

        alg -> Solver();
        cOptSolution* best = alg->GetBest();
        bestParticlevec = best->GetObjFunc(0);

        cVector bestpart = best->GetBestVec();

        cVector besttheta(sdata.NumVar);
        for (int i = 0; i < sdata.NumVar; i++)
        {
            besttheta[i] = bestpart[i];
        }

        BestThetad.push_back(besttheta);
        BestThetadOld = BestThetad;
    }

    cout << "\nMAXIMUM LIKELIHOOD ESTIMATORS" << endl;
    for (int i = 0; i < sdata.NumOut; i++)
    {
        cout << "Out: " << i+1 << "  BestThetad = ";
        BestThetad[i].Print();
    }
}

// ============================== GetSigma ================================

void cHIERKRG :: EvalWeights(cVector &mle1, cVector &theta, int out)
{
//    cout << "\n\n\n\nENTROU AQUI ================ EVAL WEIGHTS" << endl;

    mle1.Resize(sdata.NumSample);
    cVector temp(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) temp[i] = pow(10, theta[i]);

    cMatrix Psitemp;
    cMatrix Utemp;

    CalcCorrMat(sdata.SampleX, temp, Psitemp, Utemp);

    Hmat[out].Resize(sdata.NumSample, sdata.NumSample);
    Hmat[out] = Psitemp;

    Hmatdc[out].Resize(sdata.NumSample, sdata.NumSample);
    Hmatdc[out] = Utemp;

    cVector faux(sdata.NumSample);

    for (int i = 0; i < sdata.NumSample; i++)
    {
        cVector predlf(sdata.NumOut);
        EvaluateLFM(sdata.SampleX[i], predlf);
        faux[i] = predlf[out];
    }

    fVec[out].Resize(sdata.NumSample);
    fVec[out] = faux;
    double mu, sigmasqrsur;
    MudSur(mu, out, Utemp, faux);
    SigmaSqrdSur(mu, sigmasqrsur, out, Utemp, faux);
    Mu[out] = mu;
    SigmaSqr[out] = sigmasqrsur;

    cVector sy(sdata.NumSample);

    for (int j = 0; j < (sdata.NumSample); j++)
    {
        sy[j] = sdata.SampleY[j][out];
        mle1[j] = sy[j] - Mu[out]*fVec[out][j];
    }

    Utemp.SolveLU(mle1);

    Weight.push_back(mle1);

    //cVector fi2(sdata.NumSample);

    // for (int j = 0; j < sdata.NumSample; j++) fi2[j] = sdata.SampleY[j][out] - mu;

    //mle1 = fi2;

    //Utemp.SolveLU(mle1);

    //Weight.push_back(mle1);
}

// ============================== GetSigma ================================

void cHIERKRG :: EvalWeightsl(cVector &mle1, cVector &theta, int out)
{
//    cout << "\n\n\n\nENTROU AQUI ================ EVAL WEIGHTS" << endl;

    mle1.Resize(sdata.NumSampleLF);

    cVector temp(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) temp[i] = pow(10, theta[i]);

    cMatrix Psitemp;
    cMatrix Utemp;

    CalcCorrMat(sdata.SampleXLF, temp, Psitemp, Utemp);

    Psil[out].Resize(sdata.NumSampleLF, sdata.NumSampleLF);
    Psil[out] = Psitemp;

    Psildc[out].Resize(sdata.NumSampleLF, sdata.NumSampleLF);
    Psildc[out] = Utemp;

    double mu, sigmasqrsur;
    MulSur(mu, out, Utemp);
    SigmaSqrlSur(mu, sigmasqrsur, out, Utemp);
    Mul[out] = mu;
    SigmaSqrl[out] = sigmasqrsur;

    cVector fi2(sdata.NumSampleLF);

    for (int j = 0; j < sdata.NumSampleLF; j++) fi2[j] = sdata.SampleYLF[j][out] - mu;

    mle1 = fi2;

    Utemp.SolveLU(mle1);

    Weightl.push_back(mle1);

    cVector ones(sdata.NumSampleLF);
    for (int i = 0; i < sdata.NumSampleLF; i++)
        ones[i] = 1.0;

    Utemp.SolveLU(ones);
    Weightl1.push_back(ones);
}

// ========================== PredictionSur ==============================

double cHIERKRG :: PredictionSur(cVector &x, cVector &Thetal, cVector &Thetad, int out)
{
    cVector fi(sdata.NumSample);
    cVector tetad(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) tetad[i] = pow(10, Thetad[i]);

  /*  teta[0] = pow(10, 1.3218);
    teta[1] = pow(10, 0.6637);*/

    double p = 1.99;
    double dd, dist;
    if(CorrType == GAUSS)
    {
    for (int i = 0; i < (sdata.NumSample); i++)
    {
        dd = 0.0;

        for (int k = 0; k < sdata.NumVar; k++)
         dd = dd + tetad[k]*(sdata.SampleX[i][k]-x[k])*(sdata.SampleX[i][k]-x[k]);

        fi[i] = exp(-dd);
       //  cout << "psi[" << i <<"] values: " << fi[i] << endl;
    }
    }
    else
    {
        for (int i = 0; i < (sdata.NumSample); i++)
        {
            dd = 1.0;

            for (int k = 0; k < sdata.NumVar; k++)
            {
                dist = abs(sdata.SampleX[i][k]-x[k]);
                dd = dd*exp(-(dist*sqrt(5.0)/tetad[k]))*(1.0 + (dist*sqrt(5.0))/tetad[k] + (5.0*dist*dist)/(3*tetad[k]*tetad[k]));
                if (dd < 0)
                {
                    dd = 0.0;
                }
                else if (dd > 1.0)
                {
                    dd = 1.0;
                }

                fi[i] = dd;
            }
      //      cout << "psi[" << i <<"] values: " << fi[i] << endl;
        }
    }

    cVector mle1(sdata.NumSample);

    mle1 = Weight[out];

    double ylf = PredictionSurl(x, Thetal, out);

    double pred = Mu[out]*ylf + fi*mle1;

/*  cout << "Mu " << out << ": " << Mu[out] << endl;
    cout << "Prediction " << out << ": " << pred << endl;*/

    return pred;
 }

// ========================== PredictionSur ==============================

double cHIERKRG :: PredictionSurl(cVector &x, cVector &Theta, int out)
{
    cVector fi(sdata.NumSampleLF);
    cVector teta(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) teta[i] = pow(10, Theta[i]);

  /*  teta[0] = pow(10, 1.3218);
    teta[1] = pow(10, 0.6637);*/

    double p = 1.99;
    double d, dist;

    if(CorrType == GAUSS)
    {
    for (int i = 0; i < sdata.NumSampleLF; i++)
    {
        double d = 0.0;

        for (int k = 0; k < sdata.NumVar; k++)
      d = d + teta[k]*(sdata.SampleXLF[i][k]-x[k])*(sdata.SampleXLF[i][k]-x[k]);

        fi[i] = exp(-d);

       //  cout << "psi[" << i <<"] values: " << fi[i] << endl;
    }
    }
    else
    {
        for (int i = 0; i < sdata.NumSample; i++)
        {
            d = 1.0;

            for (int k = 0; k < sdata.NumVar; k++)
            {
                dist = abs(sdata.SampleXLF[i][k]-x[k]);
                d = d*exp(-(dist*sqrt(5.0)/teta[k]))*(1.0 + (dist*sqrt(5.0))/teta[k] + (5.0*dist*dist)/(3*teta[k]*teta[k]));
            }

            if (d < 0)
            {
                cout << " ====aaa: " << d << endl;
                d = 0;
            }
            else if (d > 1.0)
            {
                cout << " ====aaa: " << d << endl;
                d = 1.0;
            }
            fi[i] = d;

      //      cout << "psi[" << i <<"] values: " << fi[i] << endl;
        }
    }

    cVector mle1(sdata.NumSampleLF);

    mle1 = Weightl[out];

    double pred = Mul[out] + fi*mle1;

/*  cout << "Mu " << out << ": " << Mu[out] << endl;
    cout << "Prediction " << out << ": " << pred << endl;*/

    return pred;
 }

// ========================== SSqrSur ==============================

double cHIERKRG :: SSqrSur(cVector x, cVector &Thetal, cVector &Thetad, int out)
{
    // c vector definition (fi)

    cVector fi(sdata.NumSample);
    cVector tetad(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) tetad[i] = pow(10, Thetad[i]);

  /*  teta[0] = pow(10, 1.3218);
    teta[1] = pow(10, 0.6637);*/

    double p = 1.99;
    double dd, dist;
    if(CorrType == GAUSS)
    {
    for (int i = 0; i < (sdata.NumSample); i++)
    {
        dd = 0.0;

        for (int k = 0; k < sdata.NumVar; k++)
         dd = dd + tetad[k]*(sdata.SampleX[i][k]-x[k])*(sdata.SampleX[i][k]-x[k]);

        fi[i] = exp(-dd);
       //  cout << "psi[" << i <<"] values: " << fi[i] << endl;
    }
    }
    else
    {
        for (int i = 0; i < (sdata.NumSample); i++)
        {
            dd = 1.0;

            for (int k = 0; k < sdata.NumVar; k++)
            {
                dist = abs(sdata.SampleX[i][k]-x[k]);
                dd = dd*exp(-(dist*sqrt(5.0)/tetad[k]))*(1.0 + (dist*sqrt(5.0))/tetad[k] + (5.0*dist*dist)/(3*tetad[k]*tetad[k]));
                if (dd < 0)
                {
                    dd = 0.0;
                }
                else if (dd > 1.0)
                {
                    dd = 1.0;
                }

                fi[i] = dd;
            }
      //      cout << "psi[" << i <<"] values: " << fi[i] << endl;
        }
    }

    // F vector definition

    cVector faux(sdata.NumSample);

    faux = fVec[out];

    double sigmasqrsur;

    cMatrix Hm(sdata.NumSample, sdata.NumSample);

    Hm = Hmatdc[out];

    sigmasqrsur  = SigmaSqr[out];

    cVector s2c1 = fi;
    cVector ftemp = faux;
    Hm.SolveLU(s2c1);

    Hm.SolveLU(ftemp);

    double s2c = fi*s2c1;

    double ylf = PredictionSurl(x, Thetal, out);

    double ssqr = sigmasqrsur*(1 - s2c + (pow((ylf-faux*s2c1),2)/(faux*ftemp)));//rho*rho*sigmasqrsurl + sigmasqrsurd - s2c + (pow((1-ones*s2c1),2)/(ones*onestemp)); // LAST TERM
    // double ssqr = sigmasqrsur*(1 - s2c);

    if (ssqr < 0)
    {
        cout << "ssqr: " << ssqr << " fi:";
        fi.Print();
        cout << "encerrou" << endl;
    }
   // cout << ssqr << endl;

    return ssqr;
}

// ========================== SSqrSur ==============================

double cHIERKRG :: SSqrlSur(cVector x, cVector &Thetal, cVector &Thetad, int out)
{
    // c vector definition (fi)

    cVector fi(sdata.NumSampleLF);
    cVector tetal(sdata.NumVar);

    for (int i = 0; i < sdata.NumVar; i++) tetal[i] = pow(10, Thetal[i]);

  /*  teta[0] = pow(10, 1.3218);
    teta[1] = pow(10, 0.6637);*/

    double p = 1.99;
    double dd, dist;
    if(CorrType == GAUSS)
    {
    for (int i = 0; i < (sdata.NumSampleLF); i++)
    {
        dd = 0.0;

        for (int k = 0; k < sdata.NumVar; k++)
         dd = dd + tetal[k]*(sdata.SampleXLF[i][k]-x[k])*(sdata.SampleXLF[i][k]-x[k]);

        fi[i] = exp(-dd);
       //  cout << "psi[" << i <<"] values: " << fi[i] << endl;
    }
    }
    else
    {
        for (int i = 0; i < (sdata.NumSampleLF); i++)
        {
            dd = 1.0;

            for (int k = 0; k < sdata.NumVar; k++)
            {
                dist = abs(sdata.SampleXLF[i][k]-x[k]);
                dd = dd*exp(-(dist*sqrt(5.0)/tetal[k]))*(1.0 + (dist*sqrt(5.0))/tetal[k] + (5.0*dist*dist)/(3*tetal[k]*tetal[k]));
                if (dd < 0)
                {
                    dd = 0.0;
                }
                else if (dd > 1.0)
                {
                    dd = 1.0;
                }

                fi[i] = dd;
            }
      //      cout << "psi[" << i <<"] values: " << fi[i] << endl;
        }
    }

    // F vector definition

    cVector onesaux(sdata.NumSampleLF);

    for (int i = 0; i < sdata.NumSampleLF; i++)
    {
        onesaux[i] = 1.0;
    }

    double sigmasqrsur;

    cMatrix Hm(sdata.NumSampleLF, sdata.NumSampleLF);

    Hm = Psildc[out];

    sigmasqrsur  = SigmaSqrl[out];

    cVector s2c1 = fi;
    cVector onestemp = onesaux;
    Hm.SolveLU(s2c1);

    double s2c = fi*s2c1;

    //double ssqr = sigmasqrsur*(1 - s2c + (pow((1-onesaux*s2c1),2)/(onesaux*Weightl1[out])));//rho*rho*sigmasqrsurl + sigmasqrsurd - s2c + (pow((1-ones*s2c1),2)/(ones*onestemp)); // LAST TERM
    double ssqr = sigmasqrsur*(1 - s2c);

    if (ssqr < 0)
    {
        cout << "ssqr: " << ssqr << " fi:";
        fi.Print();
        cout << "encerrou" << endl;
    }
   // cout << ssqr << endl;

    return ssqr;
}

// ========================== SigmaSqrSur ==============================

void cHIERKRG :: SigmaSqrlSur(double &mu, double &sigmasqrsur, int out, cMatrix &Uaux)
{
    cVector sy(sdata.NumSampleLF);
    cMatrix Psitemp(sdata.NumSampleLF, sdata.NumSampleLF);

    Psitemp = Uaux;

    for (int j = 0; j < sdata.NumSampleLF; j++)
    {
        sy[j] = sdata.SampleYLF[j][out];
    }

    cVector sd1(sdata.NumSampleLF);

    for (int j = 0; j < sdata.NumSampleLF; j++)
    {
         sd1[j] = sy[j] - mu;
    }

     cVector var11 = sd1;
     Psitemp.SolveLU(var11);
     sigmasqrsur = (sd1*var11)/sdata.NumSampleLF;

    if(sigmasqrsur < 0)
    {
     //   cout << "Sigma Sqr: " << sigmasqrsur << endl;
    }
}

// ========================== SigmaSqrSur ==============================

void cHIERKRG :: SigmaSqrdSur(double &mu, double &sigmasqrsur, int out, cMatrix &Uaux, cVector faux)
{
    cVector sy(sdata.NumSample);
    cMatrix Psitemp(sdata.NumSample, sdata.NumSample);

    Psitemp = Uaux;

    for (int j = 0; j < sdata.NumSample; j++) sy[j] = sdata.SampleY[j][out];

    cVector sd1(sdata.NumSample);

    for (int j = 0; j < sdata.NumSample; j++)
    {
         sd1[j] = sy[j] - mu*faux[j];
    }

     cVector var11 = sd1;
     Psitemp.SolveLU(var11);
     sigmasqrsur = (sd1*var11)/sdata.NumSample;

    if(sigmasqrsur < 0)
    {
     //   cout << "Sigma Sqr: " << sigmasqrsur << endl;
    }
}

// ========================== Correlation Matrix ==============================

void cHIERKRG :: CalcCorrMat(vector<cVector> sampx, cVector theta, cMatrix &Psi, cMatrix &Uaux)
{
    int nsamp = sampx.size();

    Psi.Resize(nsamp, nsamp);
    Psi.Zero();

    double p = 1.99;    // gaussian

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
                //      cout << "Psi[" << i <<"][" << j << "] values: " << Psi[i][j] << endl;
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
                                    cout << "===== aaa: " << d << endl;
                    d = 0;
                }

                Psi[i][j] = d;
                //      cout << "Psi[" << i <<"][" << j << "] values: " << Psi[i][j] << endl;
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

   //cout << "MATRIZ U " << U[1][0] << endl; //"  " << U[2][2] << setprecision(4) << endl;
}

// ========================== EvalExpImp ===================================

double cHIERKRG :: EvalVFExpImp(cVector &x, double ybest, int &fid)
{
  double ei, eil, eih;

  eil = EI1(x, ybest);
  eih = EI2(x, ybest);

  if (eil > eih)
  {
      ei = eil;
      fid = 1;
  }
  else
  {
      ei = eih;
      fid = 2;
  }

  return(ei);
}

// ========================== EvalExpImp ===================================

double cHIERKRG :: EvalExpImp(cVector &x, double ybest)
{
  double ei;

  ei = EI2(x, ybest);

  return(ei);
}

// ========================== EvalExpImp ===================================

double cHIERKRG :: EI1(cVector &x, double ybest)
{
  double ei;

  double pi = atan(1)*4;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrlSur(x, BestTheta[0], BestThetad[0], 0);
  Pred[0] = ytemp[0];
  s[0] = Mu[0]*sqrt(ssqr[0]);

  if (s[0] <= 1e-12)
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

// ========================== EvalExpImp ===================================

double cHIERKRG :: EI2(cVector &x, double ybest)
{
  double ei;

  double pi = atan(1)*4;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrSur(x, BestTheta[0], BestThetad[0], 0);
  Pred[0] = ytemp[0];
  s[0] = sqrt(ssqr[0]);

  if (s[0] <= 1e-12)
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

// ========================== EvalExpImp ===================================

double cHIERKRG :: EvalVFProbImp(cVector &x, double ybest, int &fid)
{
  double pi, pil, pih;

  pil = PI1(x, ybest);
  pih = PI2(x, ybest);

  if (pil > pih)
  {
      pi = pil;
      fid = 1;
  }
  else
  {
      pi = pih;
      fid = 2;
  }

  return(pi);
}

// ========================== EvalExpImp ===================================

double cHIERKRG :: EvalProbImp(cVector &x, double ybest)
{
  double pi;

  pi = PI2(x, ybest);

  return(pi);
}

// ========================== EvalProbImp ===================================

double cHIERKRG :: PI1(cVector &x, double ybest)
{
  double PoI;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrlSur(x, BestTheta[0], BestThetad[0], 0);
  Pred[0] = ytemp[0];
  s[0] = Mu[0]*sqrt(ssqr[0]);

  if (s[0] <= 1e-12)
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

double cHIERKRG :: PI2(cVector &x, double ybest)
{
  double PoI;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrSur(x, BestTheta[0], BestThetad[0], 0);
  Pred[0] = ytemp[0];
  s[0] = sqrt(ssqr[0]);

  if (s[0] <= 1e-12)
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

// ========================== EvalExpImp ===================================

double cHIERKRG :: EvalLCB(cVector &x)
{
double lcb;

lcb = LCB2(x);

return(lcb);
}

// ========================== EvalExpImp ===================================

double cHIERKRG :: EvalVFLCB(cVector &x)
{
double lcb, lcbl, lcbh;

lcbl = LCB1(x);
lcbh = LCB2(x);

if (lcbl < lcbh)
    lcb = lcbl;
else
    lcb = lcbh;

return(lcb);
}

// ========================== EvalProbImp ===================================

double cHIERKRG :: LCB1(cVector &x)
{
  double LCB;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrlSur(x, BestTheta[0], BestThetad[0], 0);
  Pred[0] = ytemp[0];
  s[0] = Mu[0]*sqrt(ssqr[0]);

  // LCB - Test
  LCB = Pred[0] - Beta*s[0];
  // LCB = (1 - Beta)*Pred[0] - Beta*s[0];

  return(LCB);
}

// ========================== EvalProbImp ===================================

double cHIERKRG :: LCB2(cVector &x)
{
  double LCB;

  cVector ytemp(sdata.NumOut);

  Evaluate(x, ytemp);

  cVector Pred(1);
  cVector s(1);
  cVector ssqr(1);

  ssqr[0] =  SSqrSur(x, BestTheta[0], BestThetad[0], 0);
  Pred[0] = ytemp[0];
  s[0] = sqrt(ssqr[0]);

  // LCB - Test
  LCB = Pred[0] - Beta*s[0];
  // LCB = (1 - Beta)*Pred[0] - Beta*s[0];

  return(LCB);
}

// ========================== EvalInfPen ==============================

double cHIERKRG :: EvalInfPen(cVector &x, int out, int fid, double tolviol)
{
  double pf;
  double PF = 0;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  Pred = ytemp[out];

  if (Pred > -tolviol)
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

double cHIERKRG :: EvalProbFeas(cVector &x, int out, int fid, double tolviol)
{
  double pf;
  double PF = 0;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  if (fid == 2)
    ssqr = SSqrSur(x, BestTheta[out], BestThetad[out], out);
  else if (fid == 1)
    ssqr = SSqrlSur(x, BestTheta[out], BestThetad[out], out);

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

  if (s <= tolviol)
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

double cHIERKRG :: EvalProbFeasTutum(cVector &x, int out, int fid, double tolviol)
{
    double pf;
    double PF = 0;

    double ssqr(sdata.NumOut);

    cVector ytemp(sdata.NumOut);
    Evaluate(x, ytemp);

    if (fid == 2)
      ssqr = SSqrSur(x, BestTheta[out], BestThetad[out], out);
    else if (fid == 1)
      ssqr = SSqrlSur(x, BestTheta[out], BestThetad[out], out);

    double s = sqrt(ssqr);

    if (ssqr <= 0)
    {
        s = 0;
    }
    else
    {
        s = sqrt(ssqr);
    }
     // cout << "CHEGOU AQUI EVAL PF \n\n\n\n" << endl;

    double Pred = ytemp[out];

    // Pred vai com valor trocado!

   //cout << "Pred: " << Pred << endl;

    if (s < tolviol)
    {
        PF = 0;
    }
    else
    {
        double erfunc = erf((0 - Pred)/s);
        //PF = 0.5 + 0.5*erfunc;

        if (erfunc >= 1){
            PF = 0.5 + 0.5*erfunc;
        }
        else if(erfunc > 0){
            PF = 2 - erfunc;
        }
        else{
            PF = 0;
        }
        // cout << "ei = " << ei << endl;
    }
    pf = PF;

    //cout << "CHEGOU AQUI EVAL PF \n\n\n\n" << endl;

    return(pf);
}

// ========================== EvalProbFeasBagheri ==============================

double cHIERKRG :: EvalProbFeasBagheri(cVector &x, int out, int fid, double tolviol)
{
  double pf;
  double PF = 0;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  if (fid == 2)
    ssqr = SSqrSur(x, BestTheta[out], BestThetad[out], out);
  else if (fid == 1)
    ssqr = SSqrlSur(x, BestTheta[out], BestThetad[out], out);

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

  if (s <= tolviol)
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

// ========================== EvalProbFeas ==============================

double cHIERKRG :: EvalProbFeasSohst(cVector &x, int out, int fid, double tolviol)
{
  double pf;
  double PF = 0;

  double pi = atan(1)*4;
  double n = nFacSohst;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  if (fid == 2)
    ssqr = SSqrSur(x, BestTheta[out], BestThetad[out], out);
  else if (fid == 1)
    ssqr = SSqrlSur(x, BestTheta[out], BestThetad[out], out);

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

  if (s <= tolviol)
  {
      PF = 0;
  }
  else
  {
      double erfunc = erf((0 - Pred)/s);
      PF = 0.5 + 0.5*erfunc;
  }

  pf = sin(PF*(pi/2.0));
  pf = pow(pf, n);

  return pf;
}


// ============================== GetSigma ================================

void cHIERKRG :: GetSigma(cVector &s)
{

}

// ============================== Init ================================

void cHIERKRG :: Init(int i, cMatrix &xx, cMatrix &v, cMatrix &xp, double inf, double sup)
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

double cHIERKRG :: Eval(cVector temp, int m, eProbLikeType plt)
{
    /*cout << "\nm: " << m << endl;
    temp.Print();*/

    double fobj;

    if (plt == BASIC)
    {
        fobj = Likelihood(sdata.SampleXLF, temp, m);
    }
    else if (plt == DIFFMODELFIT)
    {
        fobj = Likelihoodd(sdata.SampleX, temp, m);
    }
    else
    {
        cout << "Error: Evaluate not defined (problike)" <<endl;
    }

    return fobj;
}
