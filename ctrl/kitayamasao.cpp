// -------------------------------------------------------------------------
// ego.cpp - implementation of cEGO class.
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
#include "kitayamasao.h"
#include "rbf.h"
#include "vec.h"
#include "surr.h"
#include "stdpso.h"
#include "probsurr.h"
#include "samp.h"

// -------------------------------------------------------------------------
// Class KitayamaSAO:
//

// -------------------------------------------------------------------------
// Static variables:
//

// -------------------------------------------------------------------------
// Set read functions labels:
//

/*
static bool ReadFuncRegister[] =
{
     CtrlMap( ).Insert("EXP.IMPROVEMENT.MAX", cEGO :: ReadExpImpMax),
};
*/

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cKitayamaSAO ==============================

cKitayamaSAO :: cKitayamaSAO(void) : cSAO( )
{
  Type    = KITSAO;
  Print   = false;
  SigType = kFOLDCV;
}

// ============================= ~cEGO ==============================

cKitayamaSAO :: ~cKitayamaSAO(void)
{
}

// ============================ ReadSigType =============================

void cKitayamaSAO :: ReadSigType(istream &in)
{
  cout << "entrou aqui" << endl;
  char signame[100];

  if (!Utl::ReadString(in, signame))
  {
    cout << "Error in the input of the sigma definition method." << endl;
    exit(0);
  }

  if(string(signame) == "KITAYAMA" || string(signame) == "Kitayama"
          || string(signame) == "KIT")
  {
      SigType = KITAYAMA;
  }
  else if(string(signame) == "UNIFORMKITAYAMA" || string(signame) == "UniformKitayama"
          || string(signame) == "UKIT")
  {
      SigType = KITAYAMAUNIFORM;
  }
  else if(string(signame) == "ADAPTIVESCALING" || string(signame) == "AdaptiveScaling"
          || string(signame) == "ASKIT")
  {
      SigType = ADAPTIVE_SCALING;
  }
  else if(string(signame) == "LEAVEONEOUTCROSSVALIDATION"
          || string(signame) == "LeaveOneOutCrossValidation" || string(signame) == "LOOCV")
  {
      SigType = LOOCV;
  }
  else if(string(signame) == "KFOLDCROSSVALIDATION" || string(signame) == "KFOLDCV"
          || string(signame) == "KFCV")
  {
      SigType = kFOLDCV;
  }
  else if(string(signame) == "NAKAYAMA" || string(signame) == "Nakayama"
          || string(signame) == "NAK")
  {
      SigType = NAKAYAMA;
  }
  else
  {
  cout << "Unknown sigma definition method: " << signame << endl;
  exit(0);
  }
}


// =============================== LoadReadFunc ============================

void cKitayamaSAO :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cSAO :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("SIGMA.TYPE"               ,  makeReadObj(cKitayamaSAO, ReadSigType));
}

// =============================== Solver ==================================

void cKitayamaSAO :: Solver(void)
{
  // Inicialização das variáveis globais do algoritmo

  vector<cVector> OptAllBestvar;  // Vetor que vai armazenar os melhores individuos de cada run
  vector<cVector> OptAllBestobj;  // Vetor que vai armazenar a melhor resposta de cada run
  cVector OptNRMSE(OptNum);       // Vetor que vai armazenar o menor NRMSE de cada run
  cVector OptGenToConv(OptNum);   // Vetor que vai armazenar o numero de gerações necessárias para convergência de cada run
  cVector OptSampToConv(OptNum);  // Vetor que vai armazenar o numero de amostras necessárias para convergência de cada run


  for (int opt = 0; opt < OptNum; opt++)
  {
    out << "\n%OPTIMIZATION.NUMBER\n" << opt + 1 << endl;
    // Conhecer o número de funções objetivo e restrições do problema

    int NumVar    = Prob -> GetNumVar();
    int NumConstr = Prob -> GetNumConstr();
    int NumObj    = Prob -> GetNumObj();

    // Gerar variáveis referentes à amostragem

    int nv, ns, no, flagvs, nvs;
    cVector low, upp;
    vector<cVector> sx;
    vector<cVector> sy;
    vector<cVector> vsx;
    vector<cVector> vsy;
    vector<cVector> cy;
    cVector cbest(NumConstr);
    vector<cVector> FobjExact;

    flagvs = FlagVS;
    nvs = Nvs;
    vsx = Vsx;
    vsy = Vsy;
    cout << flagvs << endl;
    cout << nvs << endl;
    cout << vsy[0][0] << endl;

    SetApproxObj( );
    SetApproxConstr( );

    cout << "\n" << endl;
    cout << "============================================================" << endl;
    cout << "\n%OPTIMIZATION.NUMBER " << opt + 1 << endl;
    out << "\n%OPTIMIZATION.NUMBER\n" << opt + 1 << endl;

    // Inicialização dos modelos substitutos que serão criados (Aproximação da função real e Função Densidade)

    cRBF *SurModel = new cRBF;
    cRBF *DensityFunction = new cRBF;

    SurModel -> SetLambda(1.0e-3);
    DensityFunction -> SetLambda(1.0e-3);

    nv = NumVar;
    no = NumObj + NumConstr;

    bool* ApproxOut = new bool[NumObj + NumConstr];
    for (int i = 0; i < (NumObj + NumConstr); i++ ){
        if ( i < NumObj ){
            ApproxOut[i] = ApproxObj[i];
        }
        else{
            ApproxOut[i] = ApproxC[i - NumObj];
        }
    }

    // Leitura do arquivo .smp

    SetInitialSample( ns, nv, no, low, upp, sx, sy, ApproxOut, cy, FobjExact);

    // Inicialização das variáveis locais (para cada otimização)

    double NRMSE = 0;
    int GenToConvSur = -1;
    int SampToConvSur = -1;
    bool Lowerror = false;

    // Construção do modelo inicial

    eSigmaType method = SigType;

    SurModel -> CreateModel(ns, nv, no, low, upp, method, sx, sy, flagvs, vsx, vsy, ApproxOut, cy, FobjExact);

    cVector allbestvar(nv);
    cVector allbestobj(NumObj);
    cVector xb(nv);
    cVector yb(no);
    SurModel -> GetBestFeasibleSample(xb, yb, 0, NumConstr, cbest);

    for (int i = 0; i < NumObj; i++){
        allbestobj[i] = yb[i];
    }

    // Construção do modelo inicial (Função Densidade)

    cVector df(1);
    vector<cVector> ones;
    for (int i = 0; i < ns; i++)
    {
        df[0] = 1.0;
        ones.push_back(df);
    }

    cMatrix sigmas(no, ns);

    SurModel -> GetSigma(sigmas);

    DensityFunction -> CreateModel(ns, nv, 1, low, upp, sigmas, sx, ones, ApproxOut, cy, FobjExact);

    // Colocar larguras no out
    out << "\n%SIGMAS" << endl;

    for (int run = 0; run < MaxGen; run++){
        cout << "------------------------------------------------------------" << endl;
        cout << "RUN = " << run + 1 << endl;

        // Inicialização de variáveis locais

        vector<cVector> newsx;
        vector<cVector> newsy;
        vector<cVector> newdf;
        cVector newpointx(nv);
        cVector newpointy(no);
        cVector yhfm(NumObj);
        cVector ct(NumConstr);
        int numsamp;


        cMatrix sigmas;
        SurModel -> GetSigma(sigmas);
        SurModel -> GetNumSample(numsamp);
        cout << "Sigma = " << sigmas[0][0] << endl;
        out << numsamp << "  " << sigmas[0][0] << endl;

        // Optimization of the SurModel surface

        cout << "best = ";
        MinimizeSurModelPrediction(newpointx, SurModel);    // O algoritmo de otimização é instanciado dentro do InfillCriteria
        newpointx.Print();

        // Is the best individual already on the sample?

        SurModel -> Normalize(newpointx);
        bool isinsamp = SurModel -> IsInSample(newpointx, 10e-6);
        SurModel -> UndoNormalization(newpointx);

        // If not, then add it to the sample

        if (!isinsamp){

        // Evaluate its HFM

        Prob -> Evaluate(newpointx, ct, yhfm);

        if (NumConstr > 0){
            if (ct.Max() <= TolViol){              // Se nenhuma restrição foi violada, compara-se yhfm com o melhor indivíduo até então. Se for melhor, atualiza esse melhor indivíduo
                if (yhfm[0] < allbestobj[0]){
                  allbestvar = newpointx;
                  allbestobj = yhfm;
                }
            }
        }
        else {
            if (yhfm[0] < allbestobj[0]){          // Se não houverem restrições, compara-se yhfm com o melhor indivíduo até então. Se for melhor, atualiza esse melhor indivíduo
              allbestvar = newpointx;
              allbestobj = yhfm;
            }
        }

        // Armazenando newpointy, concatenando yhfm e ct

        for (int i = 0; i < (NumObj + NumConstr); i++ ){
            if ( i < NumObj ){
                newpointy[i] = yhfm[i];
            }
            else{
                newpointy[i] = ct[i - NumObj];
            }
        }


        // Add it to the sample (of both the Surrogate Model and the Density Function)

        SurModel -> Normalize(newpointx);
        newsx.push_back(newpointx);
        newsy.push_back(newpointy);

        SurModel -> UpdateModel(method, newsx, newsy, ct, yhfm);

        // On the Density Function, the SampleY vector must be equal to 1.0

        SurModel -> GetSigma(sigmas);

        df[0] = 1.0;
        newdf.push_back(df);     

        DensityFunction -> UpdateModel(sigmas, newsx, newdf, ct, yhfm);

        // Clear the auxialiary vectors

        newsx.clear( );
        newsy.clear( );
        newdf.clear( );
        cout << "HFM sur = ";
        newpointy.Print();
        }

        // Optimization of the Density Function surface

        double count = nv/2;
        // double count = 1;
        for (int i = 1; i <= count; i++){
        SurModel -> GetNumSample(numsamp);
        if (numsamp >= Nmax){
            break;
        }
        out << numsamp << "  " << sigmas[0][0] << endl;

        cout << "df = ";
        MinimizeSurModelPrediction(newpointx, DensityFunction);
        newpointx.Print();

        // Evaluate the HFM of the best individual

        Prob -> Evaluate(newpointx, ct, yhfm);

        for (int i = 0; i < (NumObj + NumConstr); i++ ){
            if ( i < NumObj ){
                newpointy[i] = yhfm[i];
            }
            else{
                newpointy[i] = ct[i - NumObj];
            }
        }

        // Add it to the sample of the SurModel

        SurModel -> Normalize(newpointx);
        newsx.push_back(newpointx);
        newsy.push_back(newpointy);
        SurModel -> UpdateModel(method, newsx, newsy, ct, yhfm);

        // On the Density Function, the SampleY vector must be equal to 1.0

        SurModel -> GetSigma(sigmas);

        df[0] = 1.0;
        newdf.push_back(df);
        DensityFunction -> UpdateModel(sigmas, newsx, newdf, ct, yhfm);

        // Clear the auxiliary vectors

        newsx.clear( );
        newsy.clear( );
        newdf.clear( );

        cout << "HFM df = ";
        yhfm.Print();
        }

        /*
        vector<cVector>samplex;
        vector<cVector>sampley;
        int numsamp;
        DensityFunction -> GetSampleX(samplex);
        DensityFunction -> GetSampleY(sampley);
        DensityFunction -> GetNumSample(numsamp);
        cout << "SampleX Final = " << endl;
        for (int a = 0; a < numsamp; a++){
            samplex[a].Print();
        }*
        cout << "\n";

        cout << "SampleY Final = " << endl;
        for (int a = 0; a < numsamp; a++){
            sampley[a].Print();
        }*/

        SurModel -> GetBestFeasibleSample(xb, yb, 0, NumConstr, cbest);
        for (int i = 0; i < NumObj; i++){
            allbestobj[i] = yb[i];
        }
        cout << "Best individual = " << yb[0] << endl;

        if (FlagVS){
          double sqdiff = 0.0;
          double sqvsy = 0.0;
          double MSE = 0.0;

          for (int i = 0; i < nvs; i++){
            for (int j = 0; j < NumObj; j++){
              sqdiff += (vsy[i][j] - allbestobj[j])*(vsy[i][j] - allbestobj[j]);
              sqvsy += vsy[i][j]*vsy[i][j];
          }
          MSE = sqdiff/sqvsy;
          NRMSE = pow(MSE, 0.5);

          if (NRMSE < MinNRMSE && GenToConvSur == -1){
              Lowerror = true;
              GenToConvSur = run;
              SurModel -> GetNumSample(SampToConvSur);
              SampToConvSur = SampToConvSur - 2;
          }
        }
        }
        SurModel -> GetNumSample(numsamp);
        if (numsamp > Nmax){
            break;
        }

    }
    cout << "\n \n";
    cout << "============================================================" << endl;
    cout << "               OPTIMIZATION NUM "<< opt + 1 << "                 " << endl;
    cout << "============================================================" << endl;
    cout << "============================================================" << endl;
    cout << "                   OPTIMIZATION RESULTS                     " << endl;
    cout << "============================================================" << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "Best Individual: ";
    allbestvar.Print();
    cout << "------------------------------------------------------------" << endl;
    cout << "Fobj: ";
    allbestobj.Print();
    if (FlagVS){
    cout << "------------------------------------------------------------" << endl;
    cout << "NRMSE:   "  << NRMSE << endl;
    cout << "------------------------------------------------------------" << endl;

    if (Lowerror){
        cout << "------------------------------------------------------------" << endl;
        cout << "Iteration to achieve NRMSE < "<< MinNRMSE*100 <<"%:   "  << GenToConvSur << endl;
        cout << "Number of points needed        :   "  << SampToConvSur << endl;
        cout << "------------------------------------------------------------" << endl;
    }
    }

    OptAllBestvar.push_back(allbestvar);
    OptAllBestobj.push_back(allbestobj);
    if(FlagVS){
    OptNRMSE[opt] = NRMSE;
    OptGenToConv[opt] = GenToConvSur;
    OptSampToConv[opt] = SampToConvSur;
    }

    cout << "\n \n";
  }

  cout << "\n \n";
  cout << "============================================================" << endl;
  cout << "                           RESUMO                           " << endl;
  cout << "============================================================" << endl;
  cout << "============================================================" << endl;
  cout << "                   OPTIMIZATION RESULTS                     " << endl;
  cout << "============================================================" << endl;
  for (int i = 0; i < OptNum; i++){
  cout << "------------------------------------------------------------" << endl;
  cout << "               OPTIMIZATION NUM "<< i + 1 << "                 " << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Best Individual: ";
  OptAllBestvar[i].Print();
  cout << "------------------------------------------------------------" << endl;
  cout << "Fobj: ";
  OptAllBestobj[i].Print();
  if(FlagVS){
  cout << "------------------------------------------------------------" << endl;
  cout << "NRMSE:   "  << OptNRMSE[i] << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Iteration to achieve NRMSE < "<< MinNRMSE*100 <<"%:   "  << OptGenToConv[i] << endl;
  cout << "Number of points needed        :   "  << OptSampToConv[i] << endl;
  cout << "------------------------------------------------------------" << endl;
  }

  cout << "\n \n";
  }

  out << "\n%OBJECTIVE.FUNCTION" << endl;
  for (int i = 0; i < OptNum; i++){
      out << i + 1 << "  " << OptAllBestobj[i][0] << endl;
  }

  if (FlagVS){
  out << "\n%NRMSE" << endl;
  for (int i = 0; i < OptNum; i++){
      out << i + 1 << "  " << OptNRMSE[i] << endl;
  }

  out << "\n%NAPC" << endl;
  for (int i = 0; i < OptNum; i++){
      out << i + 1 << "  " << OptSampToConv[i] << endl;
  }
  }

}

// =============================== MinimizeSurModelPrediction ================================

void cKitayamaSAO :: MinimizeSurModelPrediction(cVector &xb, cRBF* Sur)
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
}

// ======================================================= End of file =====
