// -------------------------------------------------------------------------
// ego.cpp - implementation of cGREGO class.
// -------------------------------------------------------------------------
// Copyright (c) 2013 LMCV/UFC
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
// Created:      30-Jul-2019    Marina Alves Maia
//
// Modified:     17-Dec-2019    Marina Alves Maia
//                              Code refactoring.
// -------------------------------------------------------------------------

#include <string>
#include <vector>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "problem.h"
#include "sel.h"
#include "group.h"
#include "individual.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"
#include "ego.h"
#include "grego.h"
#include "krg.h"
#include "vec.h"
#include "surr.h"
#include "stdpso.h"
#include "stdga.h"
#include "probsurr.h"
#include "sao.h"
#include "samp.h"

// -------------------------------------------------------------------------
// Class EGO:
//

// -------------------------------------------------------------------------
// Public methods:

// =============================== cGREGO ===============================

cGREGO :: cGREGO(void)  : cSAO( )
{
  Type     = GREGO;
  Print    = false;
  WEI      = 0.5;
  ciclewei = 0;
  CorrType  = GAUSS;
  HyperParamLow = -1.0;
  HyperParamUpp = 2.0;
}

// ============================= ~cGREGO ==============================

cGREGO :: ~cGREGO(void)
{
}

// ================================ ReadWEI ================================

void cGREGO :: ReadWEI(istream &in)
{
  if (!(in >> WEI))
  {
    cout << "Error in the input of the expected improvement weight." << endl;
    exit(0);
  }
}

// =============================== ReadCyclicWei ================================

void cGREGO :: ReadCyclicWei(istream &in)
{
  if (!(in >> ciclewei))
  {
    cout << "Error in the definition of cyclic weighted expected improvement" << endl;
    exit(0);
  }
}

// ============================ ReadSigType =============================

void cGREGO :: ReadCorrType(istream &in)
{
    cout << "entrou aqui" << endl;
    char signame[100];

    if (!Utl::ReadString(in, signame))
    {
      cout << "Error in the input of the sigma definition method." << endl;
      exit(0);
    }

    if(string(signame) == "Gauss" || string(signame) == "GAUSS"
            || string(signame) == "gauss")
    {
        CorrType = GAUSS;
    }
    else if(string(signame) == "Matern" || string(signame) == "MATERN"
            || string(signame) == "matern")
    {
        CorrType = MATERN52;
    }
    else
    {
    cout << "Unknown correlation function definition: " << signame << endl;
    exit(0);
    }
}

// ======================== ReadHyperParam ==========================

void cGREGO :: ReadHyperParam(istream &in)
{
  if (!(in >> HyperParamLow) || !(in >> HyperParamUpp))
  {
    cout << "Error in the input of the hyperparameters." << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cGREGO :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cSAO :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("WEI.VALUE"                ,  makeReadObj(cGREGO, ReadWEI));
  im.Insert("USE.CYCLIC.WEI"           ,  makeReadObj(cGREGO, ReadCyclicWei));
  im.Insert("CORRELATION.TYPE"               ,  makeReadObj(cGREGO, ReadCorrType));
  im.Insert("KRG.HYPERPARAMETERS"               ,  makeReadObj(cGREGO, ReadHyperParam));
}

// =============================== Solver ==================================

void cGREGO :: Solver(void)
{
  vector<cVector> NRMSEopt;
  vector<cVector> RMAEopt;
  vector<cVector> NMAEopt;
  vector<cVector> ERRORopt;

  double w = WEI;

  vector<cVector> newsx;
  vector<cVector> newsy;
  int lastgen = MaxGen-1;

  int NumVar    = Prob->GetNumVar();      // Set number of variables
  int NumConstr = Prob->GetNumConstr();   // Set number of constraints
  int NumObj    = Prob->GetNumObj();      // Set number of objective functions

  // Set objective functions and constraints to be approximated by the surrogate model

  int nv, no;
  no = 0;

  SetApproxObj( );
  SetApproxConstr( );

  bool* ApproxOut = new bool[NumObj + NumConstr];
  for (int i = 0; i < (NumObj + NumConstr); i++ ){
      if ( i < NumObj ){
          ApproxOut[i] = ApproxObj[i];
      }
      else{
          ApproxOut[i] = ApproxC[i - NumObj];
      }
  }

  cout << "APROXIMACOES" << endl;
  for (int i = 0; i < NumObj + NumConstr; i++)
  {
      cout << ApproxOut[i] << endl;
  }

  for (int i = 0; i < NumObj+NumConstr; i++)
  {
      if (ApproxOut[i] == 1)
      {
          no += 1;
      }
  }

  cout << "NUMERO DE OUTPUTS: " << no << endl;

  nv = NumVar;
  //no = NumObj + NumConstr;

  cVector NRMSEopttemp(no);
  cVector RMAEopttemp(no);
  cVector NMAEopttemp(no);
  cVector ERRORopttemp(no);
  vector<int> GentoConvSur(OptNum);
  cVector BestY(no);
  cVector eimax(1);
  cVector ybestei(1);
  cVector cbest(NumConstr);

  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
   {
      GenStall = 0;
      bool ConvErrorFlag = false;
      int ConvError = 1e6;

     cout << "\n============================================================" << endl;
     cout << "\n%OPTIMIZATION.NUMBER\n" << opt + 1 << endl;

     // Init surrogate model

     cKRG *SurModel = new cKRG;

     int ns, flagvs, nvs;
     cVector low, upp;
     vector<cVector> sx;                 // Sampling plan
     vector<cVector> sy;
     vector<cVector> vsx;                // Validation points
     vector<cVector> vsy;
     vector<cVector> cy;
     vector<cVector> FobjExact;

     flagvs = FlagVS;
     nvs = Nvs;
     vsx = Vsx;
     vsy = Vsy;

     // Set initial sampling plan

     SetInitialSample(ns, nv, no, low, upp, sx, sy, ApproxOut, cy, FobjExact);

   /*   out << "\n\n" << endl;
     for (int i = 0; i < ns; i++)
     {
         out << sx[i][0] << " " << sx[i][1] << endl;
     }

     out << "\n\n" << endl;

     for (int i = 0; i < ns; i++)
     {
         out << sy[i][0] << endl;
     }

     exit (0);
*/

     cVector xb(nv);                    // Best individual
     cVector yb(no);

     cVector newpointx(nv);
     cVector newpointy(no);
     cVector ct;
     cVector yhfm(NumObj);


     eCorrelationType method = CorrType;

     // Create the surrogate model.

     cout << "Criar modelo. " << endl;

     SurModel->CreateModel(ns, nv, no, low, upp, method, HyperParamLow, HyperParamUpp, sx, sy, flagvs, vsx, vsy, ApproxOut, cy, FobjExact);
     Thetas.push_back(SurModel->BestTheta[0]);
     cout << "Theta: " << endl;
     SurModel->BestTheta[0].Print( );

     // Main loop

     for (int gen = 0; gen < MaxGen; gen++)
     {
         cout << "\nGEN: " << gen << endl;

         // Apply weighted expected improvement

         if (ciclewei)
         {
             double w[4] = {0.3, 0.5, 0.7, 0.5};
             int step = gen%4;
             SurModel->SetWEI(w[step]);
             cout << "w = " << w[step] << endl;
         }

         // Maximize Expected Improvement

         double ybei;
         double rand = Utl::RandDouble(0.0, 1.00);
         double egreedy = 0.10;

         if (rand < egreedy)
         {
             GenerateFeasibleSol(newpointx, ybei, SurModel);
             SurModel -> UndoNormalization(newpointx);
         }
         else
         {
             cout << "Minimizar Predicao: " << endl;
             MinimizeSurModelPrediction(newpointx, ybei, SurModel);
         }

         cout << "FINALIZOU MAX EI:" << ybei << endl;

         SurModel->GetBestFeasibleSample(xb, yb, 0, NumConstr, cbest);

    //     cout << "aTE AQUI  ok" << endl;

         eimax[0] = pow(10, -ybei);

    //     cout << "foooi" << endl;
         ybestei[0] = yb[0];

         EImax.push_back(eimax);
         Ybest.push_back(ybestei);

         if (NumConstr > 0)
         {
             Cbest.push_back(cbest);
         }

         cout << "Novo ponto da melhoria: " << endl;
         newpointx.Print();

         // Evaluate HFM

         ct.Resize(NumConstr);

         Prob->Evaluate(newpointx, ct, yhfm);

         for (int i = 0; i < no; i++)
         {
             // Store new point y, ys and ct

             if (ApproxObj[0] == 1)
             {
                 if (i == 0)
                 {
                     newpointy[i] = yhfm[i];
                 }
                 else
                 {
                     newpointy[i] = ct[i - NumObj];
                 }
             }
             else
             {
                 if (ApproxC[i] == 1)
                 {
                     newpointy[i] = ct[i];
                 }
             }
         }

         if (cy.empty( ))
         {
             ct.Resize(0);
         }

         // Normalize new sample point

         SurModel->Normalize(newpointx);

         newsx.push_back(newpointx);
         newsy.push_back(newpointy);

         // Update surrogate model

         cout << "Atualizar modelo" << endl;

         SurModel->UpdateModel(method, newsx, newsy, ct, yhfm);
         Thetas.push_back(SurModel->BestTheta[0]);

         newsx.clear();
         newsy.clear();

         // Get current best sample point

      //   cout << "Identificar melhor amostra" << endl;

         SurModel->GetBestFeasibleSample(xb, yb, 0, NumConstr, cbest);   // MARINA: GENERALIZAR PARA > 1 OUTPUT
         SurModel->UndoNormalization(xb);
         cout << "best = ";
         xb.Print();
         cout << "HFM = ";
         yb.Print();

         // Update Stallgen criterion

         if(gen == 0 || abs(BestY[0] - yb[0]) > 1e-5)
         {
             BestY = yb;
             GenStall = 0;
         }
         else
         {
             BestY = yb;
             GenStall += 1;
         }

         // Evaluate error measures

         if(FlagVS)
         {
             cVector vsytemp(no);
             vector<cVector> vsysur;

             for (int j = 0; j < nvs; j++)
             {
                 SurModel->Evaluate(vsx[j], vsytemp);
                 vsysur.push_back(vsytemp);
             }

             vsysur.push_back(yb);

             cVector temp, temp2, temp3, temp4;

             EvalErrorMeasures(vsysur, temp, temp2, temp3, temp4, no);

             vsysur.clear();

             // Update error metrics

             NRMSE.push_back(temp);
             RMAE.push_back(temp2);
             NMAE.push_back(temp3);
             ERROR.push_back(temp4);

             NRMSEopttemp = temp;
             RMAEopttemp = temp2;
             NMAEopttemp = temp3;
             ERRORopttemp = temp4;
         }

         double errormin = ERRORopttemp[0];

         int numsamp;
         SurModel->GetNumSample(numsamp);

         if (abs(errormin) < 0.01 && ConvErrorFlag == false)
         {
             ConvError = numsamp;
             ConvErrorFlag = true;
         }

         // Stopping criteria (Condition 01 - Reached max number of stall gen)

         if (GenStall == 10)
         {
                 lastgen = gen;
                 break;
         }
         else             //  (Condition 02 - Maximum number of infill points)
         {
             if (numsamp > Nmax)
             {
                 lastgen = gen;
                 break;
             }
         }
     }

     int sizeNRMSE = NRMSE.size();
     cout << "\n\n";

     for (int i = 0; i < sizeNRMSE; i++) cout << "NRMSE geracao " << i << " " << NRMSE[i][0] << endl;
     for (int i = 0; i < sizeNRMSE; i++) cout << "ERROR geracao " << i << " " << ERROR[i][0] << endl;

     out << "\n%SAMPLING.POINTS.TO.1PERC.ACCURACY" << endl;
     out << ConvError << endl;

     out << "\n%THETAS.GENERATION" << endl;
     out << opt+1 << endl;
     for (int j = 0; j < no; j++)
     {
         out << j << endl;
         for (int i = 0; i < Thetas.size( ); i++)
         {
             out << i << " ";
             for (int m = 0; m < nv; m++)
             {
                 out << Thetas[i][m] << " ";
             }
             out << endl;
         }
     }

     out << "\n%EI.GENERATION" << endl;
     out << opt+1 << endl;
     for (int i = 0; i < EImax.size( ); i++)
     {
         out << i << " ";
         out << EImax[i][0] << " " << Ybest[i][0] << endl;
     }

     out << "\n%NRMSE.GENERATION" << endl;
     out << opt+1 << endl;
     for (int j = 0; j < no; j++)
     {
         out << j << endl;
     for (int i = 0; i < sizeNRMSE; i++) out << i << " " << NRMSE[i][j] << endl;
     }

     out << "\n%ERROR.GENERATION" << endl;
     out << opt+1 << endl;
     for (int j = 0; j < no; j++)
     {
         out << j << endl;
     for (int i = 0; i < sizeNRMSE; i++) out << i << " " << ERROR[i][j] << endl;
     }

     NRMSE.clear();
     RMAE.clear();
     NMAE.clear();
     ERROR.clear();
     Thetas.clear();
     EImax.clear();
     Ybest.clear();

     // Store the best individual of the optimization.

     NRMSEopt.push_back(NRMSEopttemp);
     RMAEopt.push_back(RMAEopttemp);
     NMAEopt.push_back(NMAEopttemp);
     ERRORopt.push_back(ERRORopttemp);
     GentoConvSur[opt] = lastgen+1;

     // Print data in the output file.

     out << "\n%OBJECTIVE.FUNCTION" << endl;
     for (int i = 0; i < nv; i++) out << xb[i] << " ";
     out << endl;
     for (int i = 0; i < no; i++) out << yb[i] << " ";
     out << endl;


     int nsample;
      vector<cVector> sxfinal;
         vector<cVector> syfinal;
         vector<cVector> cyfinal;

         SurModel->GetNumSample(nsample);
         SurModel->GetSampleX(sxfinal);
         SurModel->GetSampleY(syfinal);
         SurModel->GetExactConst(cyfinal);

         out << "\n\n%NUMBER.OF.SAMPLES" << endl;
         out << nsample << endl;

         out << "\n%SAMPLE.X" << endl;

         for (int i = 0; i < nsample; i++)
         {
             for (int j = 0; j < nv; j++)
             {
                 out << sxfinal[i][j] << " ";
             }
             out << endl;
         }

         out << "\n\n%SAMPLE.Y" << endl;

         for (int i = 0; i < nsample; i++)
         {
             for (int j = 0; j < no; j++)
             {
                 out << syfinal[i][j] << " ";
             }
             out << endl;
         }

         if (!cyfinal.empty())
         {
         cout << "nsample: " << nsample << " vs " << cyfinal.size( ) << endl;

         out << "\n\n%EXACT.CONSTRAINTS" << endl;

         for (int i = 0; i < nsample; i++)
         {
             for (int j = 0; j < no; j++)
             {
                 out << cyfinal[i][j] << " ";
             }
             out << endl;
         }
         }
  }

  PostProcessingSur(OptNum, no, &NRMSEopt, &RMAEopt,&NMAEopt, &ERRORopt, &GentoConvSur);
}

// =============================== MinimizeSurModelPrediction ================================

void cGREGO :: MinimizeSurModelPrediction(cVector &xb, double &ybei, cKRG* Sur)
{
    cOptAlgorithm *alg = SubAlgType;
    cPenStatic *pen = new cPenStatic;
    pen -> SetFactor(10e15);

    alg -> SetSolType(SolType);
    alg -> SetPopSize(SubPop);
    alg -> SetMaxGen(SubMaxGen);
    alg -> SetOptNum(1);
    alg -> SetPrint(false);
    alg -> SetTolViol(SubTolViol);
    alg -> SetPenFunction(pen);
    alg -> SetMutProb(SubMutProb);

    cProblem* probbest = new cProbSurr(Sur, Prob, EVALUATE_SURROGATE);
    probbest -> SetApproxC( ApproxC );

    alg -> SetProblem(probbest);

    alg -> Init( );

    alg -> Solver();

    cOptSolution* best = alg -> GetBest();

    xb = best -> GetVec();
    ybei = best->GetObjFunc(0);
}

// =============================== GenerateFeasibleSol ================================

void cGREGO :: GenerateFeasibleSol(cVector &xb, double &ybei, cKRG* Sur)
{
    int nv = Sur->NumVar;
    int no;
    Sur->GetNumOut(no);
    double fobj = 10e20;

    cVector yout(no);

    int NumConstr = Prob->GetNumConstr();
    int NumObj = Prob->GetNumObj();

    cVector ct(NumConstr);

    bool* approxout = new bool[NumConstr+NumObj];
    Sur->GetApproxOut(approxout);

    cout << "========== entrou no aleatorio" << endl;

    while (fobj == 10e20)
    {

        for (int i = 0; i < nv; i++)         // Generating sampling point
        {
            xb[i] = Utl::RandDouble(0.0, 1.0);
        }

        Sur->Evaluate(xb, yout);            // Evaluate prediction

        if (NumConstr > 0)
        {
            if (approxout[0] == 1)
            {
                fobj = yout[0];

                for(int j = 0; j < NumConstr; j++)
                {
                    if (approxout[j+1] == 1)
                    {
                        ct[j] = yout[j+1];
                    }
                    else
                    {
                        double caux;
                        Prob->EvalExactConstraint(j, xb, caux);
                        ct[j] = caux;
                    }
                }
            }
            else
            {
                Prob->EvalExactFobj(xb, fobj);

                for(int j = 0; j < NumConstr; j++)
                {
                    if (approxout[j+1] == 1)
                    {
                        ct[j] = yout[j];
                    }
                    else
                    {
                        double caux;
                        Prob->EvalExactConstraint(j, xb, caux);
                        ct[j] = caux;
                    }
                }
            }
        }
        else
        {
            fobj = yout[0];
        }

    /*    cout << "fobj: " << fobj << endl;
        ct.Print();
        cout << "\n" << endl;*/

        if(NumConstr > 0)
        {
            if(ct.Max() > 0.0)
            {
              fobj = 10e20;
            }
        }
    }

    ybei = fobj;

}

// ======================================================= End of file =====

