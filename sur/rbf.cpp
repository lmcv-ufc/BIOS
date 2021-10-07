// -------------------------------------------------------------------------
// rbf.cpp - implementation of cRBF class.
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
// Created:      06-Fev-2018    David Sena Balreira
//                              Luana Andreza Gomes Moura
//
// Modified:     19-Oct-2019    Leonardo Gonçalves Ribeiro
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//#include <bits/stdc++.h>

#include "rbf.h"
#include "surr.h"
#include "vec.h"
#include "mat.h"
#include "matvec.h"
#include "utl.h"
#include "matvec.h"

using namespace std;

// -------------------------------------------------------------------------
// Public methods:
//

// ================================= cRBF ==================================

cRBF :: cRBF(void)
{
  Hmat = 0;
}

// ================================= ~cRBF =================================

cRBF :: ~cRBF()
{
}

// ============================ CreateModel ===============================

void cRBF :: CreateModel(const sSampData &sampdata, eSigmaType sigtype)
{
 // Store samples points and values.

 sdata = sampdata;

 // Train the RBF model.

 Mu.Resize(sdata.NumOut);
 SigmaSqr.Resize(sdata.NumOut);

 Sigma.Resize(sdata.NumOut, sdata.NumSample);

  CalcSigma(sigtype);

  CalcWeights(sdata.NumSample, sdata.NumOut, sdata.SampleX, sdata.SampleY);

  // Decomp H Mat

  DecompHmat( );

  // Evaluate Mu e SigmaSqrSur
  for(int i = 0; i < sdata.NumOut; ++i)
  {
    MuSur(i);
    SigmaSqrSur(i);
  }
}

// ============================ CreateModel ===============================

void cRBF :: CreateModel(const sSampData &sampdata, cMatrix sig)
{
 // Store samples points and values.

 sdata = sampdata;

 // Get Y bounds.
 sdata.EvalYbounds( );

 // Train the RBF model.

 Sigma.Resize(sdata.NumOut, sdata.NumSample);
 CalcSigma(sig, sdata.NumSample);

 CalcWeights(sdata.NumSample, sdata.NumOut, sdata.SampleX, sdata.SampleY);

  // Decomp H Mat.
  DecompHmat( );

  // Calculando Mu e SigmaSqrSur
  for(int i = 0; i < sdata.NumOut; ++i)
  {
    MuSur(i);
    SigmaSqrSur(i);
  }
}

// ============================ UpdateModel ===============================

void cRBF :: UpdateModel(eSigmaType sigtype, cVectorVec& sx, cVectorVec &sy)
{
  int nn = sx.size();
  sdata.NumSample += nn;

  // Store samples points and values.

  for (int i = 0; i < nn; i++)
  {
    sdata.SampleX.push_back(sx[i]);
    sdata.SampleY.push_back(sy[i]);
  }

  // Train the RBF model.

  Sigma.Resize(sdata.NumOut, sdata.NumSample);

  CalcSigma(sigtype);

  Weight.clear( );  // Remove old weight vectors
  CalcWeights(sdata.NumSample, sdata.NumOut, sdata.SampleX, sdata.SampleY);   // Evaluate new weight vectors

  // Decomp H Mat.
  DecompHmat( );

  // Calculando Mu e SigmaSqrSur
  for(int i = 0; i < sdata.NumOut; ++i)
  {
    MuSur(i);
    SigmaSqrSur(i);
  }
}

void cRBF :: UpdateModel(eSigmaType sigtype, vector<cVector> sx, vector<cVector> sy, cVector cy, cVector fobjex)
{
  int nn = sx.size();
  sdata.NumSample += nn;

  // Store samples points and values.

  for (int i = 0; i < nn; i++)
  {
    sdata.SampleX.push_back(sx[i]);
    sdata.SampleY.push_back(sy[i]);

    if (cy.Dim( ) > 0)
    {
        sdata.ExactCons.push_back(cy);
    }

    if (fobjex.Dim( ) > 0)
    {
        sdata.ExactFobj.push_back(fobjex);
    }
  }

  // Update Y bounds.
  sdata.EvalYbounds( );

  // Train the RBF model.

  Sigma.Resize(sdata.NumOut, sdata.NumSample);

  CalcSigma(sigtype);

  Weight.clear( );  // Remove old weight vectors
  CalcWeights(sdata.NumSample, sdata.NumOut, sdata.SampleX, sdata.SampleY);   // Evaluate new weight vectors

  // Decomp H Mat.
  DecompHmat( );

  // Calculando Mu e SigmaSqrSur
  for(int i = 0; i < sdata.NumOut; ++i)
  {
    MuSur(i);
    SigmaSqrSur(i);
  }
}

// ============================ UpdateModel ===============================

void cRBF :: UpdateModel(cMatrix sig, vector<cVector> sx, vector<cVector> sy, cVector cy, cVector fobjex)
{
  int nn = sx.size();
  sdata.NumSample += nn;

  // Store samples points and values.

  for (int i = 0; i < nn; i++)
  {
    sdata.SampleX.push_back(sx[i]);
    sdata.SampleY.push_back(sy[i]);

    if (cy.Dim( ) > 0)
    {
        sdata.ExactCons.push_back(cy);
    }

    if (fobjex.Dim( ) > 0)
    {
        sdata.ExactFobj.push_back(fobjex);
    }
  }

   // Train the RBF model.

  Sigma.Resize(sdata.NumOut, sdata.NumSample);

  CalcSigma(sig, sdata.NumSample);

  Weight.clear( );  // Remove old weight vectors
  CalcWeights(sdata.NumSample, sdata.NumOut, sdata.SampleX, sdata.SampleY);   // Evaluate new weight vectors

  // Decomp H Mat.
  DecompHmat( );

  // Evaluate Mu e SigmaSqrSur
  for(int i = 0; i < sdata.NumOut; ++i)
  {
    MuSur(i);
    SigmaSqrSur(i);
  }

}

// ================================ Evaluate ===============================

void cRBF :: Evaluate(cVector &x, cVector &y, vector<cVector> *sampx)
{
    int flag = 0;

    if (sampx == NULL) { flag = 0;}           // No sample has been specified
    else  {flag = 1;}                         // Otherwise

    // Evaluate the distance from {x} to RBF centers.

    cVector dist(NumBasis);


   if (flag == 1)             // If sample is specified (Cross validation)
   {
      for (int i = 0; i < NumBasis; i++)
      {
        double r = 0.0;
        for (int j = 0; j < sdata.NumVar; j++)
        {
          double d = x[j] - (*sampx)[i][j];
          r += d*d;
        }

        dist[i] = r; // dist[i]= dist^2
      }
   }
   else                   // Nakayama and Kitayama
   {
       for (int i = 0; i < NumBasis; i++)
       {
         double r = 0.0;
         for (int j = 0; j < sdata.NumVar; j++)
         {
           double d = x[j] - sdata.SampleX[i][j];
           r += d*d;
         }
         dist[i] = r; // dist[i]= dist^2
       }
   }

   // Evaluate the RBF output.

   for (int i = 0; i < sdata.NumOut; i++)
   {
     double f = 0.0;
     for (int j = 0; j < NumBasis; j++){
       f += exp(-dist[j]/Sigma[i][j])*Weight[i][j];
     }
     y[i] = f;
   }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= CalcDistMax ===============================

void cRBF :: CalcDistMax(cVector &dmax)
{
  cVector dij(sdata.NumVar);
  dmax.Zero( );
  for (int i = 0; i < sdata.NumSample; i++)
  {
    for (int j = 0; j < sdata.NumSample; j++)
    {
      dij = sdata.SampleX[j] - sdata.SampleX[i];
      double dist = dij.Length( );
      dmax[i] = max(dmax[i], dist);
      dmax[j] = max(dmax[j], dist);
    }
  }
}

// ============================= CalcDistMax ===============================

void cRBF :: CalcDistMax(cVector &dmax, double sclfac)
{
  cVector dij(sdata.NumVar);
  dmax.Zero( );
  for (int i = 0; i < sdata.NumSample; i++)
  {
    for (int j = 0; j < sdata.NumSample; j++)
    {
      dij = sdata.SampleX[j] - sdata.SampleX[i];
      dij *= sclfac;
      double dist = dij.Length( );
      dmax[i] = max(dmax[i], dist);
      dmax[j] = max(dmax[j], dist);
    }
  }
}

// ============================== CalcSigma ================================

void cRBF :: CalcSigma(eSigmaType sigtype)
{
  if (sigtype == NAKAYAMA)
    CalcSigmaNakayama( );
  else if (sigtype == KITAYAMAUNIFORM)
    CalcSigmaKitayamaUniform( );
  else if (sigtype == KITAYAMA)
    CalcSigmaKitayama( );
  else if (sigtype == LOOCV)
    CalcSigmaLOOCV( );
  else if (sigtype == kFOLDCV)
    CalcSigmakFOLDCV( );
  else
    CalcSigmaKitayamaAdpScl( );
}

// ============================== GetSigma ================================

void cRBF :: GetSigma(cMatrix &s)
{
    s.Resize(sdata.NumOut, sdata.NumSample);

    int add = 0;

    if (sdata.ApproxOut[0] == 0)
    {
        add = 1;
    }

    for (int j = 0; j < sdata.NumOut; j++){
        if(sdata.ApproxOut[j+add] == 1){
          for (int i = 0; i < sdata.NumSample; i++)
          {
            s[j][i] = sqrt(Sigma[j][i]);
          }
        }
    }
}

// ============================== CalcSigma ================================

void cRBF :: CalcSigma(cMatrix signumber, int sampnum)
{
  for (int i = 0; i < sampnum; i++)
    for (int j = 0; j < sdata.NumOut; j++)
        Sigma[j][i] = signumber[j][i]*signumber[j][i];
}

// ======================= CalcSigmaKitayamaAdpScl =========================

void cRBF :: CalcSigmaKitayamaAdpScl(void)
{
  // Intialization.

  double a = sdata.NumSample - 1.0;
  double b = 1.0/(double)sdata.NumVar;
  cVector dmax(sdata.NumSample);

  // Adaptive scaling loop.

  int add = 0;

  if (sdata.ApproxOut[0] == 0)
  {
      add = 1;
  }

  double sclfac = 1.0;
  double alpha = 1.1;
  for (int iter = 0; iter < 200; iter++)
  {
    CalcDistMax(dmax, sclfac);

    double rmin = 1e30;
    for (int i = 0; i < sdata.NumSample; i++)
    {
        for (int j = 0; j < sdata.NumOut; j++){
          if (sdata.ApproxOut[j+add] == 1){
            double r = dmax[i]/(sqrt(sdata.NumVar)*pow(a, b));
            Sigma[j][i] = r*r;
            rmin = min(rmin, r);
          }
        }
    }

//    cout << "rmin = " << rmin << endl;

    if (rmin > 1.0) break;

    sclfac *= alpha;
  }
}

// ========================== CalcSigmaKitayama ============================

void cRBF :: CalcSigmaKitayama(void)
{
  cVector dmax(sdata.NumSample);
  CalcDistMax(dmax);

  double a = sdata.NumSample - 1.0;
  double b = 1.0/(double)sdata.NumVar;

  int add = 0;

  if (sdata.ApproxOut[0] == 0)
  {
      add = 1;
  }

  for (int i = 0; i < sdata.NumSample; i++)
  {
      for (int j = 0; j < sdata.NumOut; j++){
        if(sdata.ApproxOut[j+add] == 1){
          double r = dmax[i]/(sqrt(sdata.NumVar)*pow(a, b));
          Sigma[j][i] = r*r;
        }
      }
  }
}


// ========================== CalcSigmaKitayamaUniform ============================

void cRBF :: CalcSigmaKitayamaUniform(void)
{
  cVector dmax(sdata.NumSample);
  CalcDistMax(dmax);

  double b = 1.0/(double)sdata.NumVar;
  double Dmax = dmax.Max( );
  double r = Dmax/(sqrt(sdata.NumVar)*pow((double)sdata.NumSample, b));
  for (int i = 0; i < sdata.NumOut; i++)
  {
    if(sdata.ApproxOut[i] == 1){
      for (int j = 0; j < sdata.NumSample; j++){
          Sigma[i][j] = r*r;
      }
    }
  }
}


// ========================== CalcSigmaNakayama ============================

void cRBF :: CalcSigmaNakayama(void)
{
  cVector dmax(sdata.NumSample);
  CalcDistMax(dmax);

  double b = 1.0/(double)sdata.NumVar;
  double Dmax = dmax.Max( );
  double r = Dmax/(pow((double)sdata.NumVar*(double)sdata.NumSample, b));
  for (int i = 0; i < sdata.NumOut; i++)
  {
    if(sdata.ApproxOut[i] == 1){
      for (int j = 0; j < sdata.NumSample; j++){
          Sigma[i][j] = r*r;
      }
    }
  }
}

// ========================== CalcSigmaLOOCV ============================

void cRBF :: CalcSigmaLOOCV( )
{
  vector<double> sigmai(20);                                  //This array will contain all values possible (0.01 to 10 in a logarithmic scale)
  cVector Errmin(sdata.NumOut);                               //This will contain the minimum error achieved, and serve as comparation
  for (int i = 0; i < sdata.NumOut; i++){
      Errmin(i) = 10e30;
  }
  cVector Sigmac(sdata.NumOut);                               //This stores the current sigma with the minimum error
  cVector nrmse(sdata.NumOut);

  for (double i = 0; i < 20; i++)                             //Loop to store the possible sigma values to the sigmai list
  {
    double pot = -2+3*i/19;
    sigmai[i] = pow(10, pot) * sqrt(2);
  }

  for (int i = 0; i < 20; i ++)                               //Sigma to be tested on this iteration
  {
  double err = 0;

    for (int j = 0; j < sdata.NumSample; j++)                 //Index of the sample used as validation test on this iteration
    {

      vector<cVector> SampleXi = sdata.SampleX;               //Creation of the samples vector (for variables and objective functions)
      vector<cVector> SampleYi = sdata.SampleY;

      SampleXi.erase (SampleXi.begin()+j);                    //Erasing the sample used as validation
      SampleYi.erase (SampleYi.begin()+j);

      cMatrix sigvec(sdata.NumOut, sdata.NumSample-1);

      for (int jj = 0; jj < sdata.NumOut; jj++){
        for (int ii = 0; ii < sdata.NumSample-1; ii++)
        {
            sigvec[jj][ii] = sigmai[i];
        }
      }

      Sigma.Resize(sdata.NumOut, sdata.NumSample-1);

      CalcSigma(sigvec, sdata.NumSample-1);

      // Lambda = 1.0e-3;
      CalcWeights(sdata.NumSample-1, sdata.NumOut, SampleXi, SampleYi);

      cVector yi(sdata.NumOut);
      Evaluate(sdata.SampleX[j], yi, &SampleXi);                                //Evaluating the validation sample



      for (int nt = 0; nt < sdata.NumOut; nt++)
      {
        double diff = yi[nt]-sdata.SampleY[j][nt];
        double nerr = pow(pow(diff, 2)/1,0.5);
        if (isinf(nerr) || isnan(nerr))
        {
          err = abs(diff);
        }
        else
        {
          err = nerr;
        }
        nrmse[nt] = err;
      }
    }

    for (int nt = 0; nt < sdata.NumOut; nt++){
      if (nrmse[nt]<Errmin[nt])                                          //Comparing with previously established sigma and, if lower, define the new sigma
      {
        Sigmac[nt] = sigmai[i];
        Errmin[nt] = nrmse[nt];
      }
    }
  }

  cMatrix sigvec(sdata.NumOut, sdata.NumSample);
  for (int j = 0; j < sdata.NumOut; j++){
    for (int i = 0; i < sdata.NumSample; i++)
    {
        sigvec[j][i] = Sigmac[j];
    }
  }

  Sigma.Resize(sdata.NumOut, sdata.NumSample);

  CalcSigma(sigvec, sdata.NumSample);
}

// ========================== CalcSigmakFOLDCV ============================

void cRBF :: CalcSigmakFOLDCV( )
{
  vector<double> sigmai(20);                                    //This array will contain all values possible (0.01 to 10 in a logarithmic scale)*sqrt(2)
  cVector Errmin(sdata.NumOut);                                 //This will contain the minimum error achieved, and serve as comparation
  for (int i = 0; i < sdata.NumOut; i++){
      Errmin(i) = 10e30;
  }
  cVector Sigmac(sdata.NumOut);                                 //This stores the current sigma with the minimum error
  cVector nrmse(sdata.NumOut);
  int k = 5;
  if (sdata.NumSample < k)
  {
    k = sdata.NumSample;
  }

  cVector groupsize(k);

  for (double i = 0; i < 20; i++)                               //Loop to store the possible sigma values to the sigmai list
  {
    double pot = -2+3*i/19;
    sigmai[i] = pow(10, pot) * sqrt(2.0);
  }

  int n = sdata.NumSample/k;
  int N = k - sdata.NumSample % k;

  for (int i = 0; i < k; i++)
  {
      if (i < N)
      {
        groupsize[i] = n;
      }
      else
      {
        groupsize[i] = n+1;
      }
  }

  vector<int> sampindex(sdata.NumSample);
  cMatrix groups(k, n+1);

  for (int i = 0; i < sdata.NumSample; i++)
  {
      sampindex[i] = i;
  }

  random_shuffle(sampindex.begin(), sampindex.end());

  int s = 0;
  for (int i = 0; i < k; i++)
  {
      for (int j = 0; j < groupsize[i]; j++)
      {
        groups[i][j] = sampindex[s];
        s += 1;
      }
  }

  for (int i = 0; i < 20; i++)                               //Sigma to be tested on this iteration
  {
  double err = 0;

    for (int j = 0; j < k; j++)                              //Index of the group used as validation test on this iteration
    {

      vector<cVector>  SampleXi = sdata.SampleX;             //Creation of the samples vector (for variables and objective functions)
      vector<cVector>  SampleYi = sdata.SampleY;

      for (int g = 0; g < groupsize[j]; g++)
      {
        int gs = groupsize[j] - g;
        SampleXi.erase (SampleXi.begin()+gs);                //Erasing the group used as validation
        SampleYi.erase (SampleYi.begin()+gs);
      }

      cMatrix sigvec(sdata.NumOut, sdata.NumSample-groupsize[j]);

      for (int jj = 0; jj < sdata.NumOut; jj++){
        for (int ii = 0; ii < sdata.NumSample-groupsize[j]; ii++)
        {
          sigvec[jj][ii] = sigmai[i];
        }
      }

      Sigma.Resize(sdata.NumOut, sdata.NumSample-groupsize[j]);

      CalcSigma(sigvec, sdata.NumSample-groupsize[j]);
      // Lambda = 1.0e-3;
      CalcWeights(sdata.NumSample-groupsize[j], sdata.NumOut, SampleXi, SampleYi);

      cVector yi(sdata.NumOut);
      cMatrix Yi(groupsize[j], sdata.NumOut);

      for (int g = 0; g < groupsize[j]; g++)
      {
           Evaluate(sdata.SampleX[groups[j][g]], yi, &SampleXi); //Evaluating the validation sample
        for (int nt = 0; nt < sdata.NumOut; nt++)
        {
            Yi[g][nt] = yi[nt];
        }
      }

      for (int nt = 0; nt < sdata.NumOut; nt++)
      {
        double diff = 0;
        double nerr = 0;
        for (int g = 0; g < groupsize[j]; g++)
        {
          diff = (Yi[g][nt] - sdata.SampleY[groups[j][g]][nt]);
          nerr = pow(pow(diff,2)/groupsize[j],0.5);
          if (isinf(nerr))
          {
            err += abs(diff);
          }
          else
          {
            err += nerr;
          }
        }
        nrmse[nt] = err;
      }

    }
    for (int nt = 0; nt < sdata.NumOut; nt++){
      if (nrmse[nt] < Errmin[nt])                                          //Comparing with previously established sigma and, if lower, define the new sigma
      {
        Sigmac[nt] = sigmai[i];
        Errmin[nt] = nrmse[nt];
      }
    }
  }

  cMatrix sigvec(sdata.NumOut, sdata.NumSample);
  for (int j = 0; j < sdata.NumOut; j++){
    for (int i = 0; i < sdata.NumSample; i++)
    {
        sigvec[j][i] = Sigmac[j];
    }
  }
  Sigma.Resize(sdata.NumOut, sdata.NumSample);

  CalcSigma(sigvec, sdata.NumSample);
}

// ============================= CalcWeights ===============================

void cRBF :: CalcWeights(int nsamp, int nout, vector<cVector> sampx, vector<cVector> sampy)
{
  // Create virtual matrices.
  if (Hmat == 0) Hmat = new cMatrix [nout];

  for (int out = 0; out < sdata.NumOut; out ++){
      Hmat[out].Resize(nsamp, nsamp);
  }
  // Clear weight vector.

  if (Weight.size( ) > 0)
    Weight.clear( );

  // Define the number of centers (basis).

  NumBasis = nsamp;

  // Evaluate matrix [H].

  double d,s;

  cMatrix H(nsamp, nsamp);

  for(int out = 0; out < nout; out++)
    for (int i = 0; i < nsamp; i++)
      for (int j = 0; j < nsamp; j++)
      {
        s = 0.0;
        for (int k = 0; k < sdata.NumVar; k++)
	{
          d = sampx[i][k] - sampx[j][k];
          s += d*d;
        }
        Hmat[out][i][j] = exp(-s/Sigma[out][j]);
      }

  // Evaluate [A] = [H]t[H] + lambda*[I].

  for (int out = 0; out < nout; out++)
  {
    H = Hmat[out];

    cMatrix A(nsamp, nsamp);
    cMatrix Ht(nsamp, nsamp);
    H.Transp(Ht);
    A = Ht*H;
    for (int i = 0; i < nsamp; i++)
      A[i][i] += Lambda;

  // Evaluate the weight vectors {w}.

    cVector y(nsamp);
    cVector w(NumBasis);

    A.DecompLU( );

    for (int j = 0; j < nsamp; j++)
    {
      y[j] = sampy[j][out];
    }

    w = Ht*y;
    A.SolveLU(w); // Solve [A]{w} = [Ht]{y}

    Weight.push_back(w);  // Store in the array of weights
  }
}

// ============================= CalcHMatrix ===============================

void cRBF :: CalcH(vector<cVector> &H)
{
    //Define the number of centers (basis)
    NumBasis = sdata.NumSample;

    //Centers = sdata.SampleX

    double s;
    cMatrix U(sdata.NumSample, sdata.NumSample);
    for (int out = 0; out < sdata.NumOut; out++){
    for (int i = 0; i < sdata.NumSample; i++)
        for (int j = 0; j < sdata.NumSample; j++)
        {
            s = 0;
            for (int k = 0; k < sdata.NumVar; k++)
            {
                s += (sdata.SampleX[i][k] - sdata.SampleX[j][k]) * (sdata.SampleX[i][k] - sdata.SampleX[j][k]);
            }
                U[i][j]= s/Sigma[out][j];
                H[i][j] = exp(-U[i][j]);
        }
    }
}

// ============================= CalcFi ===============================

double cRBF :: CalcFi(double d, int out, int s)
{
    double fi = exp(-((d*d)/Sigma[out][s]));

    return fi;
}

// ========================== DecompHmat ==============================

void cRBF :: DecompHmat( )
{
    Hmatdc = Hmat;

    for (int out = 0; out < sdata.NumOut; out++){
      Hmatdc[out].DecompLU();
    }
}

// ========================== MuSur ==============================

void cRBF :: MuSur(int out)
{
        Mu[out] = 0;
}

// ========================== SigmaSqrSur ==============================

void cRBF :: SigmaSqrSur(int out)
{
        SigmaSqr[out] = 1.0;
}

// ========================== SSqrSur ==============================

double cRBF :: SSqrSur(cVector x, int out)
{
    double dist;
    cVector d(sdata.NumVar);
    cVector fi(sdata.NumSample);
    cVector v1(sdata.NumSample);

    for (int i = 0; i < sdata.NumSample; i++)
    {
        d = x - sdata.SampleX[i];
        dist = d.Length( );

        fi[i] = CalcFi(dist, out, i);
    }

    for (int j = 0; j < sdata.NumSample; j++)                                //Vector of ones
    {
        v1[j] = 1.0;
    }

    cVector s2c1 = fi;
    Hmatdc[out].SolveLU(s2c1);

    double S2c1 = (fi*s2c1);

    double ssqr = SigmaSqr[out]*(1 - S2c1); // + ((1 - S2c2)*(1 - S2c2))/S2c3);

    return ssqr;
}

// ========================== EvalExpImp ===================================

double cRBF :: EvalExpImp(cVector &x, double ybest)
{
  double ei;

  double pi = atan(1)*4;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  cVector Pred(1);

  ssqr = SSqrSur(x, 0);
  Pred[0] = ytemp[0];

  double s = sqrt(ssqr);

  bool isin = false;
  if (s < 100)
  {
  isin = IsInSample(x, 1e-5);
  }

  if (s <= 0.0 || isin == true)
  {
      ei = 0;
  }
  else
  {
      double EI1 = ybest - Pred[0];
      double erf1 = erf(EI1/(pow(2*ssqr, 0.5)));
      double EI2 = (0.50) + (0.50)*erf1;
      double EI3 = sqrt(ssqr)*(1/(sqrt(2*pi)));
      double EI4 = exp(-(EI1*EI1)/(2*ssqr));
      ei = WEI*EI1*EI2 + (1 - WEI)*EI3*EI4;
  }

  return(ei);
}

// ========================== EvalProbImp ===================================

double cRBF :: EvalProbImp(cVector &x, double ybest)
{
  double PoI;

  double pi = atan(1)*4;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  cVector Pred(1);

  ssqr = SSqrSur(x, 0);
  Pred[0] = ytemp[0];

  double s = sqrt(ssqr);

  bool isin = false;
  if (s < 100)
  {
  isin = IsInSample(x, 1e-5);
  }

  if (s <= 0.0 || isin == true)
  {
      PoI = 0;
  }
  else
  {
      double term = ybest - Pred[0];
      double erf1 = erf(term/(pow(2*ssqr, 0.5)));
      PoI = erf1;
  }

  return(PoI);
}



// ========================== EvalProbImp ===================================

double cRBF :: EvalLCB(cVector &x)
{
  double LCB;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  cVector Pred(1);

  ssqr = SSqrSur(x, 0);
  Pred[0] = ytemp[0];

  double s = sqrt(ssqr);

  bool isin = false;
  if (s < 100)
  {
  isin = IsInSample(x, 1e-5);
  }

  LCB = Pred[0] - Beta*s;

  return(LCB);
}

// ========================== EvalInfPen ==============================

double cRBF :: EvalInfPen(cVector &x, int out, double tolviol)
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

double cRBF :: EvalProbFeas(cVector &x, int out, double tolviol)
{
  double pf;
  double PF = 0;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  ssqr = SSqrSur(x, out);

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

double cRBF :: EvalProbFeasTutum(cVector &x, int out, double tolviol)
{
  double pf;
  double PF = 0;

 //x.Print();

 /* x[0] = 0;
  x[1] = 0;*/

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  ssqr = SSqrSur(x, out);

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

  if (s <= 0 || isin == true)
  {
      PF = 0;
  }
  else
  {
      double erfunc = erf((0 - Pred)/s);

      if (erfunc >= 1)
      {
          PF = 0.5 + 0.5*erfunc;
      }
      else if(erfunc > 0){
          PF = 2 - erfunc;
      }
      else
      {
          PF = 0;
      }
  }

  pf = PF;

  return(pf);
}

// ========================== EvalProbFeasBagheri ==============================

double cRBF :: EvalProbFeasBagheri(cVector &x, int out, double tolviol)
{
  double pf;
  double PF = 0;

  double ssqr;

  cVector ytemp(sdata.NumOut);
  Evaluate(x, ytemp);

  double Pred;

  ssqr = SSqrSur(x, out);

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

// ======================================================= End of file =====
