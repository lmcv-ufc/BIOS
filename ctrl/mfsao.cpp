// -------------------------------------------------------------------------
// mfsao.cpp - implementation of cMFSAO class.
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
// Created:      09-Jun-2021    Leonardo Gonçalves Ribeiro
// -------------------------------------------------------------------------

#include <string>
#include <vector>
#include <chrono>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "mfsao.h"
#include "sao.h"
#include "sampsao.h"
#include "problem.h"
#include "sel.h"
#include "group.h"
#include "individual.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"
#include "rbf.h"
#include "vec.h"
#include "surr.h"
#include "modpso.h"
#include "stdabc.h"
#include "stdga.h"
#include "stdais.h"
#include "stdde.h"
#include "lamga.h"
#include "modpso.h"
#include "rs.h"
#include "probsurr.h"
#include "samp.h"

// -------------------------------------------------------------------------
// Class SAO:
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

// ============================== cMFSAO ==============================

cMFSAO :: cMFSAO(void) : cOptAlgorithm( )
{
    SubPop         = 250;//100;
    SubMaxGen      = 100;//250;
    SubStallGen    = 100;//250;
    SubTolViol     = 1e-5;
    SubMutProb     = 0.02;
    Nmax           = 10e30;
    SubAlgType     = new cStandardPSO;
    MinNRMSE       = 0.01;
    SampleFileName = fname;
    InfillCriteria = EVALUATE_EXPECTED_IMPROVEMENT;

    WEI  = 0.5;
    Beta = 1.0;
    NFac = 0.15;

    CicleSize = 3;
    ListWEI.Resize(CicleSize);  ListWEI[0]  = 0.20; ListWEI[1]  = 0.35; ListWEI[2]  = 0.50;
    ListBeta.Resize(CicleSize); ListBeta[0] = 1.00; ListBeta[1] = 2.00; ListBeta[2] = 3.00;

    InputNS     = 0;
    InputNSLF   = 0;
    NestSamp    = 1;
    FlagVS      = 0;
    Nvs         = 0;
}

// =============================== LoadReadFunc ============================

void cMFSAO :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("SUB.POPULATION.SIZE"                  ,  makeReadObj(cMFSAO, ReadSubPop));
  im.Insert("SUB.MAXIMUM.GENERATIONS"              ,  makeReadObj(cMFSAO, ReadSubMaxGen));
  im.Insert("SUB.STALL.GEN"                        ,  makeReadObj(cMFSAO, ReadSubStallGen));
  im.Insert("SUB.CONSTRAINT.TOLERANCE"             ,  makeReadObj(cMFSAO, ReadSubTolViol));
  im.Insert("SUB.PSO.TOPOLOGY"                     ,  makeReadObj(cMFSAO, ReadSubPSOTopology));
  im.Insert("SUB.DE.TYPE"                          ,  makeReadObj(cMFSAO, ReadSubDEType));
  im.Insert("SUB.MUTATION.PROBABILITY"             ,  makeReadObj(cMFSAO, ReadSubMutProb));
  im.Insert("MAXIMUM.NUMBER.OF.POINTS"             ,  makeReadObj(cMFSAO, ReadNmax));
  im.Insert("MINIMUM.NRMSE"                        ,  makeReadObj(cMFSAO, ReadMinNRMSE));
  im.Insert("SUB.OPTIMIZATION.ALGORITHM"           ,  makeReadObj(cMFSAO, ReadSubAlgType));
  im.Insert("NUMBER.OF.INITIAL.SAMPLING.POINTS"    ,  makeReadObj(cMFSAO, ReadNumInitSamplingPoints));
  im.Insert("NUMBER.OF.INITIAL.LF.SAMPLING.POINTS" ,  makeReadObj(cMFSAO, ReadNumInitLFSamplingPoints));
  im.Insert("USE.NESTED.SAMPLE"                    ,  makeReadObj(cMFSAO, ReadUseNestedSample));
  im.Insert("SAMPLE.FILE.NAME"                     ,  makeReadObj(cMFSAO, ReadSampleFileName));
  im.Insert("VALIDATION.SAMPLES"                   ,  makeReadObj(cMFSAO, ReadValSamples));
  im.Insert("CONSTRAINT.HANDLING.METHOD"           ,  makeReadObj(cMFSAO, ReadConstrMethod));
  im.Insert("INFILL.CRITERIA"                      ,  makeReadObj(cMFSAO, ReadInfillCriteria));
  im.Insert("USE.CYCLIC.WEIGHTS"                   ,  makeReadObj(cMFSAO, ReadCicleWEI));
  im.Insert("WEI.VALUE"                            ,  makeReadObj(cMFSAO, ReadWEI));
  im.Insert("BETA.VALUE"                           ,  makeReadObj(cMFSAO, ReadBeta));
  im.Insert("N.FACTOR.SOHST"                       ,  makeReadObj(cMFSAO, ReadNFacSohst));
  im.Insert("CYCLIC.WEI.VALUE"                     ,  makeReadObj(cMFSAO, ReadCyclicWEI));
  im.Insert("CYCLIC.BETA.VALUE"                    ,  makeReadObj(cMFSAO, ReadCyclicBeta));
}

// =========================== ReadSubPop ===========================

void cMFSAO :: ReadSubPop(istream &in)
{
  if (!(in >> SubPop))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ========================= ReadSubMaxGen ==========================

void cMFSAO :: ReadSubMaxGen(istream &in)
{
  if (!(in >> SubMaxGen))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ========================= ReadSubMaxGen ==========================

void cMFSAO :: ReadSubStallGen(istream &in)
{
  if (!(in >> SubStallGen))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ======================== ReadSubTolViol ==========================

void cMFSAO :: ReadSubTolViol(istream &in)
{
  if (!(in >> SubTolViol))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ======================== ReadSubMutProb ==========================

void cMFSAO :: ReadSubMutProb(istream &in)
{
  if (!(in >> SubMutProb))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ======================== ReadSubMutProb ==========================

void cMFSAO :: ReadWEI(istream &in)
{
  if (!(in >> WEI))
  {
    cout << "Error in the input of the WEI weigth." << endl;
    exit(0);
  }
}

// ======================== ReadSubMutProb ==========================

void cMFSAO :: ReadBeta(istream &in)
{
  if (!(in >> Beta))
  {
    cout << "Error in the input of the beta factor for the Lower Confidence Bound criterion." << endl;
    exit(0);
  }
}

// ======================== ReadSubMutProb ==========================

void cMFSAO :: ReadNFacSohst(istream &in)
{
  if (!(in >> NFac))
  {
    cout << "Error in the input of the beta factor for the Lower Confidence Bound criterion." << endl;
    exit(0);
  }
}

// ======================== ReadWEI ==========================

void cMFSAO :: ReadCyclicWEI(istream &in)
{
  if (!(in >> CicleSize))
  {
    cout << "Error in the input of the size of the cyclic WEI." << endl;
    exit(0);
  }

  ListBeta.Resize(CicleSize); ListWEI.Resize(CicleSize);
  ListBeta.Zero( ); ListWEI.Zero( );

  for (int i = 0; i < CicleSize; i++)
  {
      if (!(in >> ListWEI[i]))
      {
        cout << "Error in the input of the WEI weigth (ID = " << i+1 << ")." << endl;
        exit(0);
      }
  }
}

// ======================== ReadBeta ==========================

void cMFSAO :: ReadCyclicBeta(istream &in)
{
    if (!(in >> CicleSize))
    {
      cout << "Error in the input of the size of the cyclic Beta." << endl;
      exit(0);
    }

    ListBeta.Resize(CicleSize); ListWEI.Resize(CicleSize);
    ListBeta.Zero( ); ListWEI.Zero( );

    for (int i = 0; i < CicleSize; i++)
    {
        if (!(in >> ListBeta[i]))
        {
          cout << "Error in the input of Beta (ID = " << i+1 << ")." << endl;
          exit(0);
        }
    }
}

// ======================== ReadSampleFileName ==========================

void cMFSAO :: ReadSampleFileName(istream &in)
{
  if (!(in >> SampleFileName))
  {
    cout << "Error in the input of the sample file name." << endl;
    exit(0);
  }
  SampleFileName = SampleFileName.substr(1, SampleFileName.size() - 2);
}

// ======================== ReadSubAlgType ==========================

void cMFSAO :: ReadSubAlgType(istream &in)
{
    // Read the algorithm label.

    char label[100];

    if (!Utl::ReadString(in, label))
    {
      cout << "Error in the input of the algorithm label." << endl;
      exit(0);
    }

    // Create apropriate optimization algorithm.
    if (string(label)=="StdGA" || string(label)=="stdga")
      SubAlgType = new cStandardGA( );
    else if (string(label)=="StdPSO" || string(label)=="stdpso")
      SubAlgType = new cStandardPSO( );
    else if (string(label)=="StdABC" || string(label)=="stdabc")
      SubAlgType = new cStandardABC( );
    else if (string(label)=="StdAIS" || string(label)=="stdais")
      SubAlgType = new cStandardAIS( );
    else if (string(label)=="RandSearch" || string(label)=="rs")
      SubAlgType = new cRandomSearch( );
    else if (string(label)=="StdDE" || string(label)=="stdde")
      SubAlgType = new cStandardDE( );
    else if (string(label)=="LamGA" || string(label)=="lamga")
      SubAlgType = new cLaminateGA;
    else if (string(label)=="LamPSO" || string(label)=="lampso")
      SubAlgType = new cLaminatePSO( );
    else
    {
        cout << "Unknown algorithm: " << label << endl;
        exit(0);
    }
}

// ======================== ReadSubPSOTopology ==========================

void cMFSAO :: ReadSubPSOTopology(istream &in)
{
    // Read the algorithm label.

    char label[100];

    if (!Utl::ReadString(in, label))
    {
      cout << "Error in the input of the algorithm label." << endl;
      exit(0);
    }

    // Create apropriate optimization algorithm.
    if (string(label)=="Ring" || string(label)=="ring")
      SubTopology = RING_TOPOLOGY;
    else if (string(label)=="Square" || string(label)=="Von Neumann" || string(label)=="Von-Neumann")
      SubTopology = SQUARE_TOPOLOGY;
    else if (string(label)=="GlobalBest" || string(label)=="GBest"  || string(label)=="Gbest" || string(label)=="Global")
      SubTopology = GBEST_TOPOLOGY;
    else
    {
        cout << "Unknown topology: " << label << endl;
        exit(0);
    }
}

// ======================== ReadConstrMethod ==========================

void cMFSAO :: ReadConstrMethod(istream &in)
{
    // Read the algorithm label.

    char label[100];

    if (!Utl::ReadString(in, label))
    {
      cout << "Error in the input of the algorithm label." << endl;
      exit(0);
    }

    // Create apropriate optimization algorithm.
    if (string(label)=="SOBESTER" || string(label)=="Sobester")
      ConstrMethod = INFINITE_PEN;
    else if (string(label)=="SCHONLAU" || string(label)=="Schonlau" || string(label)=="POF")
      ConstrMethod = POF_SCHONLAU;
    else if (string(label)=="TUTUM" || string(label)=="Tutum"  || string(label)=="FFT")
      ConstrMethod = POF_TUTUM;
    else if (string(label)=="BAGHERI" || string(label)=="Bagheri"  || string(label)=="FFB")
      ConstrMethod = POF_BAGHERI;
    else if (string(label)=="SOHST" || string(label)=="Sohst"  || string(label)=="FFS")
      ConstrMethod = POF_SOHST;
    else
    {
        cout << "Unknown constraint handling method: " << label << endl;
        exit(0);
    }
}

// ======================== ReadConstrMethod ==========================

void cMFSAO :: ReadInfillCriteria(istream &in)
{
    // Read the algorithm label.

    char label[100];

    if (!Utl::ReadString(in, label))
    {
      cout << "Error in the input of the algorithm label." << endl;
      exit(0);
    }

    // Create apropriate optimization algorithm.
    if (string(label)=="SURROGATEMINIMIZATION" || string(label)=="SURMIN" || string(label)=="SurrogateMinimization")
      InfillCriteria = EVALUATE_SURROGATE;
    else if (string(label)=="LOWERCONFIDENCEBOUND" || string(label)=="LowerConfidenceBound" || string(label)=="LCB")
      InfillCriteria = EVALUATE_LOWER_CONFIDENCE_BOUND;
    else if (string(label)=="PROBABILITYOFIMPROVEMENT" || string(label)=="ProbabilityOfImprovement"  || string(label)=="POI")
      InfillCriteria = EVALUATE_PROBABILITY_IMPROVEMENT;
    else if (string(label)=="EXPECTEDIMPROVEMENT" || string(label)=="ExpectedImprovement"  || string(label)=="EI")
      InfillCriteria = EVALUATE_EXPECTED_IMPROVEMENT;
    else if (string(label)=="VARIABLEFIDELITYLOWERCONFIDENCEBOUND" || string(label)=="VariableFidelityLowerConfidenceBound" || string(label)=="VFLCB")
      InfillCriteria = EVALUATE_VF_LOWER_CONFIDENCE_BOUND;
    else if (string(label)=="VARIABLEFIDELITYPROBABILITYOFIMPROVEMENT" || string(label)=="VariableFidelityProbabilityOfImprovement"  || string(label)=="VFPOI")
      InfillCriteria = EVALUATE_VF_PROBABILITY_IMPROVEMENT;
    else if (string(label)=="VARIABLEFIDELITYEXPECTEDIMPROVEMENT" || string(label)=="VariableFidelityExpectedImprovement"  || string(label)=="VFEI")
      InfillCriteria = EVALUATE_VF_EXPECTED_IMPROVEMENT;
    else
    {
        cout << "Unknown infill criteria: " << label << endl;
        exit(0);
    }
}

// ======================== ReadSubDEType ==========================

void cMFSAO :: ReadSubDEType(istream &in)
{
    // Read the algorithm label.

    char label[100];

    if (!Utl::ReadString(in, label))
    {
      cout << "Error in the input of the algorithm label." << endl;
      exit(0);
    }

    if (string(label) == "Rand1" || string(label) == "rand1" || string(label) == "Default" || string(label) == "def")
      SubDifType = Rand1;
    else if (string(label) == "LocalToBest" || string(label) == "Local2Best" || string(label) == "l2b" || string(label) == "local-to-best")
      SubDifType = Loc2Best;
    else if (string(label) == "BestWithJitter" || string(label) == "BestJitter" || string(label) == "BwJ" || string(label) == "best-with-jitter")
      SubDifType = BestJitter;
    else
    {
        cout << "Unknown topology: " << label << endl;
        exit(0);
    }
}

// =========================== ReadNmax =============================

void cMFSAO :: ReadNmax(istream &in)
{
  if (!(in >> Nmax))
  {
    cout << "Error in the input of the crossover rate." << endl;
    exit(0);
  }
}

// ======================== ReadSubMinNRMSE =========================

void cMFSAO :: ReadMinNRMSE(istream &in)
{
  if (!(in >> MinNRMSE))
  {
    cout << "Error in the input of the minimum NRMSE." << endl;
    exit(0);
  }
}

// ======================== ReadSubMinNRMSE =========================

void cMFSAO :: ReadCicleWEI(istream &in)
{
  if (!(in >> ciclewei))
  {
    cout << "Error in the input of the boolean related to the use of cyclic weights (should be 0 or 1)." << endl;
    exit(0);
  }
}

// ======================== ReadNumInitSamplingPoints =========================

void cMFSAO :: ReadNumInitSamplingPoints(istream &in)
{
  InputNS = 1;
  if (!(in >> NumInitSP))
  {
    cout << "Error in the input of the number of initial sampling points." << endl;
    exit(0);
  }
}

// ======================== ReadNumInitSamplingPoints =========================

void cMFSAO :: ReadNumInitLFSamplingPoints(istream &in)
{
  InputNSLF = 1;
  if (!(in >> NumInitSPLF))
  {
    cout << "Error in the input of the number of initial low-fidelity sampling points." << endl;
    exit(0);
  }
}

// ======================== ReadNumInitSamplingPoints =========================

void cMFSAO :: ReadUseNestedSample(istream &in)
{
  if (!(in >> NestSamp))
  {
    cout << "Error in the definition of the sample type. Use 1 for a nested sample, or 0 for a not nested sample." << endl;
    exit(0);
  }
}

// ======================== ReadNumInitSamplingPoints =========================

void cMFSAO :: ReadValSamples(istream &in)
{
  if (!(in >> Nvs))
  {
    cout << "Error in the input of the number of validation points." << endl;
    exit(0);
  }

  if (Nvs > 0) FlagVS = 1;

  int nv;
  if (!(in >> nv))
  {
    cout << "Error in the input of the number of the dimension of validation points." << endl;
    exit(0);
  }

  for (int i = 0; i < Nvs; i++)
  {
      cVector tempsx(nv);
      for (int j = 0; j < nv; j++)
      {
          if (!(in >> tempsx[j]))
          {
            cout << "Error in the input of the validation points." << endl;
            exit(0);
          }
      }

      Vsx.push_back(tempsx);
  }
  cVector tempsy(1);
  for (int i = 0; i < Nvs; i++)
  {
      if (!(in >> tempsy[0]))
      {
        cout << "Error in the input of the validation points." << endl;
        exit(0);
      }
      Vsy.push_back(tempsy);
  }
}


// ======================== SetApproxObj =========================

void cMFSAO :: SetApproxObj( )
{
/*
  int numobj = Prob -> GetNumObj( );
  bool* approxobj = new bool[numobj];
  approxobj = Prob -> GetApproxObj( );
  ApproxObj = approxobj;

  NumApproxObj = 0;

  for ( int i = 0; i < numobj; i++ ){
      if(ApproxObj[i] == 1){
          NumApproxObj += 1;
      }
  }
  */
}

// ======================== SetApproxConstr =========================

void cMFSAO :: SetApproxConstr( )
{
/*
    int numc = Prob -> GetNumConstr();
    bool* approxc = new bool[numc];
    approxc = Prob->GetApproxConstr( );
    ApproxC = approxc;
    NumApproxC = 0;
    for ( int i = 0; i < numc; i++ ){
        if(ApproxC[i] == 1){
            NumApproxC += 1;
        }
    }
*/
}

// ======================== SetInitialSample =========================

void cMFSAO :: SetInitialSample(sProbAppOut &appout, cSampSet* &set,sSampData &sdata, int &ev)
{

  // Sem arquivo smp...

  // Aqui as amostras são geradas pela classe Samp.
  ///
  // como seria feito?
  //
  // Seta o tamanho da amostra

    cout << "\n\nSample file not found!\n\n";
    cout << "Generating sample..." << endl;

    int nv = Prob -> VarNumEff( ); //Prob->VarNumRow( ) * Prob->VarNumCol( );
    int nslf = (InputNSLF) ? NumInitSPLF : 1.5*(nv + 1)*(nv + 2)/2.0;
    int ns   = (InputNS) ? NumInitSP : static_cast<int>(0.5*nslf)+1;
    set    = new cSampSet(nslf+ns,SolType,Prob,appout);

    // Create sample points.
    cout << "nv: " << nv << endl;
    cout << "nslf: " << nslf << endl;
    cout << "ns: " << ns << endl;

    if (NestSamp)
    {
        vector<cVector> sx;
        sx.reserve(nslf);
        cSamp InitialSample;

        bool check = 1;
        int count = 0;
        int countmax = 50;

        while (check == 1)
        {
            check = 0;
            sx.clear();
            InitialSample.InitSample(SampType, nv, nslf, sx);

            // Copy sample points to sample set - LF sample
            for (int i = 0; i < nslf; i++)
            {
              (*set)[i]->Init(sx[i]);
            }

            // Copy sample points to sample set - HF sample
            for (int i = 0; i < ns; i++)
            {
              // (*set)[i + nslf]->Init(sx[i]);
              (*set)[i + nslf]->Init(sx[nslf - ns + i]);
            }

            // Check if any two points are equal
            for (int m = 0; m < nslf; m++)
            {
                cVector sm;
                sm.Resize(nv);
                (*set)[m]->GetNormVar(sm);
                for (int n = 0; n < m; n++)
                {
                    double dist = 0;
                    cVector sn;
                    sn.Resize(nv);
                    (*set)[n]->GetNormVar(sn);
                    for (int k = 0; k < nv; k++)
                    {
                      dist += abs(sn[k] - sm[k]);
                    }
                    if (dist <= 1e-10)
                    {
                        check = 1;
                        count += 1;
                    }
                }
            }

            if (count > countmax)
            {
                cout << "\nBIOS was not able to generate a sample with enough different sampling points." << endl;
                cout << "Please, change the sampling technique or lower the number of initial sampling points\n" << endl;
                exit(0);
            }
        }

        double TLF, THF;
        TLF = THF = 0.0;
        // Evaluate each sampling point
        #pragma omp parallel for num_threads(omp_maxthread)
        for (int i = 0; i < nslf+ns; i++)
        {
          // cout << "i = " << i << endl;
          if (i < nslf)
          {
              //auto start = chrono::steady_clock::now();
              (*set)[i]->EvaluateLFP( );
              //auto end = chrono::steady_clock::now();
              //TLF += chrono::duration_cast<chrono::microseconds>(end - start).count();
          }
          else
          {
              //auto start = chrono::steady_clock::now();
              (*set)[i]->Evaluate( );
              //auto end = chrono::steady_clock::now();
              //THF += chrono::duration_cast<chrono::microseconds>(end - start).count();
          }
        }

        //cout << "MeanHF: " << THF/ns << endl;
        //cout << "MeanLF: " << TLF/ns << endl;
        // exit(0);

        ev+=ns;

        // Essa parte poderia estar depois...

        // Copy sample points to surrogate model, normalized and with correct number
        // of output (only approximated fobj + constraints).

        sdata.NumVar      = nv;
        sdata.NumSample   = ns;
        sdata.NumSampleLF = nslf;
        sdata.NumOut      = appout.GetNumAppOut( );

        sdata.SampleX.resize(ns);
        sdata.SampleY.resize(ns);

        sdata.SampleXLF.resize(nslf);
        sdata.SampleYLF.resize(nslf);
        for (int i = 0; i < nslf+ns; i++)
        {
          if (i < nslf)
          {
              // Input data.
              sdata.SampleXLF[i].Resize(nv);
              (*set)[i]->GetNormVar(sdata.SampleXLF[i]);

              // Output data.
              sdata.SampleYLF[i].Resize(sdata.NumOut);
              (*set)[i]->GetSurrOutRes(sdata.SampleYLF[i]);
          }
          else
          {
              // Input data.
              sdata.SampleX[i-nslf].Resize(nv);
              (*set)[i]->GetNormVar(sdata.SampleX[i-nslf]);

              // Output data.
              sdata.SampleY[i-nslf].Resize(sdata.NumOut);
              (*set)[i]->GetSurrOutRes(sdata.SampleY[i-nslf]);
          }
        }

        cout << "HIGH-FIDELITY SAMPLE" << endl;

        for (int i = 0; i < ns; i++)
        {
          cout << "X[" << i << "] = ";
          sdata.SampleX[i].Print( );
          cout << endl;
        }

        for (int i = 0; i < ns; i++)
        {
          cout << "Y[" << i << "] = ";
          sdata.SampleY[i].Print( );
          cout << endl;
        }

        cout << "LOW-FIDELITY SAMPLE" << endl;

        for (int i = 0; i < nslf; i++)
        {
          cout << "X[" << i << "] = ";
          sdata.SampleXLF[i].Print( );
          cout << endl;
        }

        for (int i = 0; i < nslf; i++)
        {
          cout << "Y[" << i << "] = ";
          sdata.SampleYLF[i].Print( );
          cout << endl;
        }
    }
    else
    {
        vector<cVector> sxlf;
        sxlf.reserve(nslf);
        vector<cVector> sx;
        sx.reserve(ns);
        cSamp InitialSample;

        bool check = 1;
        int count = 0;
        int countmax = 50;

        while (check == 1)
        {
            check = 0;
            sx.clear();
            InitialSample.InitSample(SampType, nv, nslf, sxlf);
            InitialSample.InitSample(SampType, nv, ns, sx);

            // Copy sample points to sample set - LF sample
            for (int i = 0; i < nslf; i++)
            {
              (*set)[i]->Init(sxlf[i]);
            }

            // Copy sample points to sample set - HF sample
            for (int i = 0; i < ns; i++)
            {
              // (*set)[i + nslf]->Init(sx[i]);
              (*set)[i + nslf]->Init(sx[i]);
            }

            // Check if any two points are equal - LF sample
            for (int m = 0; m < nslf; m++)
            {
                cVector sm;
                sm.Resize(nv);
                (*set)[m]->GetNormVar(sm);
                for (int n = 0; n < m; n++)
                {
                    double dist = 0;
                    cVector sn;
                    sn.Resize(nv);
                    (*set)[n]->GetNormVar(sn);
                    for (int k = 0; k < nv; k++)
                    {
                      dist += abs(sn[k] - sm[k]);
                    }
                    if (dist <= 1e-10)
                    {
                        check = 1;
                        count += 1;
                    }
                }
            }

            // Check if any two points are equal - HF sample
            for (int m = 0; m < ns; m++)
            {
                if (check == 1) break;
                cVector sm;
                sm.Resize(nv);
                (*set)[m+nslf]->GetNormVar(sm);
                for (int n = 0; n < m; n++)
                {
                    double dist = 0;
                    cVector sn;
                    sn.Resize(nv);
                    (*set)[n+nslf]->GetNormVar(sn);
                    for (int k = 0; k < nv; k++)
                    {
                      dist += abs(sn[k] - sm[k]);
                    }
                    if (dist <= 1e-10)
                    {
                        check = 1;
                        count += 1;
                    }
                }
            }

            if (count > countmax)
            {
                cout << "\nBIOS was not able to generate a sample with enough different sampling points." << endl;
                cout << "Please, change the sampling technique or lower the number of initial sampling points\n" << endl;
                exit(0);
            }
        }

        // Evaluate each sampling point
        #pragma omp parallel for num_threads(omp_maxthread)
        for (int i = 0; i < nslf+ns; i++)
        {
          if (i < nslf)
            (*set)[i]->EvaluateLFP( );
          else
            (*set)[i]->Evaluate( );
        }
        ev+=ns;

        // Essa parte poderia estar depois...

        // Copy sample points to surrogate model, normalized and with correct number
        // of output (only approximated fobj + constraints).

        sdata.NumVar      = nv;
        sdata.NumSample   = ns;
        sdata.NumSampleLF = nslf;
        sdata.NumOut      = appout.GetNumAppOut( );

        sdata.SampleX.resize(ns);
        sdata.SampleY.resize(ns);

        sdata.SampleXLF.resize(nslf);
        sdata.SampleYLF.resize(nslf);
        for (int i = 0; i < nslf+ns; i++)
        {
          if (i < nslf)
          {
              // Input data.
              sdata.SampleXLF[i].Resize(nv);
              (*set)[i]->GetNormVar(sdata.SampleXLF[i]);

              // Output data.
              sdata.SampleYLF[i].Resize(sdata.NumOut);
              (*set)[i]->GetSurrOutRes(sdata.SampleYLF[i]);
          }
          else
          {
              // Input data.
              sdata.SampleX[i-nslf].Resize(nv);
              (*set)[i]->GetNormVar(sdata.SampleX[i-nslf]);

              // Output data.
              sdata.SampleY[i-nslf].Resize(sdata.NumOut);
              (*set)[i]->GetSurrOutRes(sdata.SampleY[i-nslf]);
          }
        }

        cout << "HIGH-FIDELITY SAMPLE" << endl;

        for (int i = 0; i < ns; i++)
        {
          cout << "X[" << i << "] = ";
          sdata.SampleX[i].Print( );
          cout << endl;
        }

        for (int i = 0; i < ns; i++)
        {
          cout << "Y[" << i << "] = ";
          sdata.SampleY[i].Print( );
          cout << endl;
        }

        cout << "LOW-FIDELITY SAMPLE" << endl;

        for (int i = 0; i < nslf; i++)
        {
          cout << "X[" << i << "] = ";
          sdata.SampleXLF[i].Print( );
          cout << endl;
        }

        for (int i = 0; i < nslf; i++)
        {
          cout << "Y[" << i << "] = ";
          sdata.SampleYLF[i].Print( );
          cout << endl;
        }
    }
}

// ======================== ReadSampleFile =========================

void cMFSAO :: ReadSampleFile(ifstream &finp, int &ns, int &nv, int &no, cVector &low, cVector &upp,
                             vector<cVector> &sx, vector<cVector> &sy, vector<cVector> &cy, vector<cVector> &FobjExact)
{
    bool* readlabel = new bool[7];
    for (int i = 0; i < 7; i++)
    {
        readlabel[i] = 0;
    }

    string label;
    bool normalize = 0;

    while (finp >> label)
    {
        if (label == "%NUMBER.OF.SAMPLES")
        {
          finp >> ns;
          readlabel[0] = 1;
        }
        if (label == "%NUMBER.OF.VARIABLES")
        {
          finp >> nv;
          readlabel[1] = 1;
        }
        if (label == "%NORMALIZE")
        {
          finp >> normalize;
          readlabel[2] = 1;
        }
        if (label == "%LOWER.AND.UPPER.BOUNDS")
        {
          low.Resize(nv);
          upp.Resize(nv);
          for (int i = 0; i < nv; i++) finp >> low[i];
          for (int i = 0; i < nv; i++) finp >> upp[i];
          readlabel[3] = 1;
        }
        if (label == "%SAMPLE.X")
        {
          cVector xn(nv);
          for (int i = 0; i < ns; i++)
          {
            for (int j = 0; j < nv; j++)
            {
              finp >> xn[j];
              if (normalize)
                xn[j] = (xn[j] - low[j])/(upp[j] - low[j]);
            }
              sx.push_back(xn);
          }
          readlabel[4] = 1;
        }
        if (label == "%NUMBER.OF.OUTS")
        {
          finp >> no;
          readlabel[5] = 1;
        }
        if (label == "%SAMPLE.Y")
        {
          cVector yout(no);
          for (int i = 0; i < ns; i++)
          {
		   	  double aux;
              finp >> aux;
              if (aux == -123)
              {
                  return;
              }
            for (int j = 0; j < no; j++) finp >> yout[j];
            sy.push_back(yout);
          }
          readlabel[6] = 1;
        }
		if (label == "%EXACT.CONSTRAINTS")
        {
            int nexconst;
            finp >> nexconst;
            cVector exconst(nexconst);
            for (int i = 0; i < ns; i++)
            {
              for (int j = 0; j < nexconst; j++) finp >> exconst[j];
              cy.push_back(exconst);
            }
            readlabel[7] = 1;
        }
        if (label == "%EXACT.FOBJ")
        {
            cVector exfobj(1);
            for (int i = 0; i < ns; i++)
            {
              for (int j = 0; j < 1; j++) finp >> exfobj[j];
              FobjExact.push_back(exfobj);
            }
            readlabel[8] = 1;
        }
    }

    for (int i = 0; i < 7; i++)
    {
        bool stop = 0;
        if (readlabel[i] == 0)
        {
          cout << "The following label must be added to the sample file: ";
          if (i == 0) cout << "%NUMBER.OF.SAMPLES" << endl;
          if (i == 1) cout << "%NUMBER.OF.VARIABLES" << endl;
          if (i == 2) cout << "%NORMALIZE" << endl;
          if (i == 3) cout << "%LOWER.AND.UPPER.BOUNDS" << endl;
          if (i == 4) cout << "%SAMPLE.X" << endl;
          if (i == 5) cout << "%NUMBER.OF.OUTS" << endl;
          if (i == 6) cout << "%SAMPLE.Y" << endl;
          stop = 1;
        }
        if (stop)
        {
          cout << "Fill in the missing labels or remove the sample file to generate the sample automatically." << endl;
          exit(0);
        }
    }
}

// =============================== Solver ==================================

void cMFSAO :: Solver(void)
{
  //
  //  SAO Alg ( )
  //  {
  //     // (1) Fase de inicialização
  //      - Variáveis de pós-processamento.
  //      - Variáveis locais.
  //      - Iniciallização do Surrogate, avaliação do número de outputs. (Depende do problema de otimização)
  //
  //     // Para Opt de 1 até NumOpt
  //     // |
  //     // | - Inicialização do modeloSurrogate
  //     // | - Inicialização de variáveis de pós-processamento (por opt).
  //     // | - Avaliar Sample inicial.
  //     // | - 
  //     // | - Para pnt de 1 até NumPnt
  //     // |   |
  //     // |   |  - Escolher novo ponto pnew. (Infill)
  //     // |   |  - Atualizar modelo substituto.  SurrMod->Update(pnew). Aqui entra os Weight expected improvement.
  //     // |   |  - Stoppin criteria.
  //     // | - Pós-processamento por Opt.
  //
  //       - Pós-processamento global, de todas as optimizações..
  //  }
  //


  cSampSAO* newsmp;
  cVector   BestObj(OptNum);
  BestObj.Zero( );

  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Track number of individual evaluations.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Track the time spent in each phase

    double tev, tbuild, tinf;
    tinf = 0.0;

    // Create the population, mating pool and parent array.
    cSampSet    *smp;
    sProbAppOut appout(Prob);
    sSampData   sdata;

    auto start = chrono::steady_clock::now();
    SetInitialSample(appout,smp,sdata,EvalNum);
    auto end = chrono::steady_clock::now();
    tev = chrono::duration_cast<chrono::microseconds>(end - start).count();

    // Evaluate penalized objective function in samples.
    // Note: Pobj is not considered in infill procedure, but is used to select
    // the best feasible sample. It is adopted to handle situtations where no
    // feasible sample is avaliable.
    if (Pen)
      Pen->EvalPenObjFunc(smp, TolViol);

    // Track model hyperparameters in each iteration

    //ThetaIt.Resize(1, sdata.NumVar);
    //ThetaIt.Zero( );

    // Create the surrogate model.

    start = chrono::steady_clock::now();
    cSURR *SurModel = CreateSurrogate(sdata);
    end = chrono::steady_clock::now();
    tbuild = chrono::duration_cast<chrono::microseconds>(end - start).count();

    SurModel -> SetWEI(WEI);
    SurModel -> SetBeta(Beta);
    SurModel -> SetNFacSohst(NFac);

    // Store Hyperparameters

    cVector besttheta(sdata.NumVar);
    besttheta = SurModel -> GetBestThetad(0);

    PrintHyperPar(besttheta, sdata.NumVar, -1);

    // Evaluate initial sample points.

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    nHigFidSamp.Resize(MaxGen);
    nLowFidSamp.Resize(MaxGen);
    Tinf.Resize(MaxGen);
    Teval.Resize(MaxGen);
    Tbuild.Resize(MaxGen);
    UpdateTimeVar(opt, 0, sdata.NumSample, sdata.NumSampleLF, tinf, tev, tbuild);

    for (int step = 0; step < MaxGen; step++)
    {
      if ((step+1)%1 == 0 && Feedback) cout << "Step: " << step + 1 << endl;

      if (ciclewei)
      {
        int stepw = step % CicleSize;
        //ListWEI.Print( );
        SurModel -> SetWEI(ListWEI[stepw]);
        SurModel -> SetBeta(ListBeta[stepw]);
        if (Feedback) cout << "w    = " << ListWEI[stepw] << endl;
        if (Feedback) cout << "beta = " << ListBeta[stepw] << endl;
      }

      // Select new points.
      cVector  pntx(SurModel->GetSampData( ).NumVar);
      cVector npntx(SurModel->GetSampData( ).NumVar);
      cVector  pnty(SurModel->GetSampData( ).NumOut);

      cVector npntxl(SurModel->GetSampData( ).NumVar);
      cVector  pntyl(SurModel->GetSampData( ).NumOut);

      double  objf;
      start = chrono::steady_clock::now();
      EvalInfillCriteria(SurModel,smp,appout,pntx,objf);
      end = chrono::steady_clock::now();
      tinf = chrono::duration_cast<chrono::microseconds>(end - start).count();

      // Evaluate new point - Both the LF and HF sources are evaluated

      cVectorVec nx, ny, nyl;
      nx.push_back(pntx);
      bool addhf, addlf;

      double currbest = smp->BestSol( )->GetObjFunc(0);
      GetFidelityNewPoint(nx, addhf, addlf, currbest);

      start = chrono::steady_clock::now();
      if (addhf)
      {
          // High-Fidelity Eval
          newsmp = smp->PushBack(pntx);
          newsmp->Evaluate( );
          newsmp->GetSurrOutRes(pnty);
          ny.push_back(pnty);
          EvalNum++;
      }

      if (addlf)
      {
        // Low-Fidelity Eval
        newsmp = smp->PushBack(pntx);
        newsmp->EvaluateLFP( );
        newsmp->GetSurrOutRes(pntyl);
        nyl.push_back(pntyl);
      }
      end = chrono::steady_clock::now();
      tev = chrono::duration_cast<chrono::microseconds>(end - start).count();

      //SurModel->UpdateModel(method, nx, ny);
      start = chrono::steady_clock::now();
      UpdateSurrogate(nx,ny,nyl);
      end = chrono::steady_clock::now();
      tbuild = chrono::duration_cast<chrono::microseconds>(end - start).count();

      // Store Hyperparameters

      besttheta = SurModel -> GetBestThetad(0);
      PrintHyperPar(besttheta, sdata.NumVar, step);

      // Update variables related to PostProcessing.

      if (Pen) Pen->EvalPenObjFunc(smp, TolViol);
      UpdatePostVar(step, opt, lastBest,smp);

      int nh, nl;
      smp -> GetNumEval(nh, nl);
      if (Feedback) 
      {
        newsmp->Print( );
        cout << "Number of high-fidelity samples: " << nh << endl;
        cout << "Number of low-fidelity samples: " << nl << endl;
        cout << "Best Sample" << endl;
        smp->BestSol( )->Print( );
      }

      UpdateTimeVar(opt, step + 1, nh, nl, tinf, tev, tbuild);
    
      // Check conditions to stop optimization

      if (OptStopCrit(step, opt, lastBest, smp))
        break;   
    }

    // Store the best individual.
    best->Insert(smp->BestSol( ));

    // Print data in the output file.
    PrintPostVar(MaxGen, opt, EvalNum, smp);

    BestObj[opt] = smp->BestSol( )->GetObjFunc(0);

    PrintPostVarSur(MaxGen, opt, EvalNum, SurModel);

    int nhf, nlf;
    smp -> GetNumEval(nhf, nlf);

    PrintTimeVar( );
    NumberHFE.push_back(nhf);
    NumberLFE.push_back(nlf);
  }

  double meanhf = 0;
  double meanlf = 0;
  for (int opt = 0; opt < OptNum; opt++)
  {
      double nHF = NumberHFE[opt];
      double nLF = NumberLFE[opt];
      cout << "Opt Num " << opt+1 << "  nHFE = " << NumberHFE[opt] << "  nLFE = " << NumberLFE[opt] << "   BestObj = " << BestObj[opt] << endl;
      meanhf += nHF/OptNum;
      meanlf += nLF/OptNum;
  }
  cout << "\nMeanHF = " << meanhf << "  MeanLF = " << meanlf << endl;
}

// ============================ PrintHyperPar ============================

void cMFSAO :: PrintHyperPar(cVector theta, int nv, int step)
{
  if (!out) return;

  if (step == -1)
  {
      *out << "\n%HYPERPARAMETES.ITERATIONS\n";
  }


  for (int i = 0; i < nv; i++)
  {
      *out << theta[i] << "   ";
  }
  *out << endl;
}

// =============================== OptStopCrit =============================

bool cMFSAO :: OptStopCrit(int gen, int opt, double &lb, cGroup *mg)
{
  // Call standard optimization stopping criteria
  bool res = cOptAlgorithm :: OptStopCrit(gen,opt,lb,mg);

  int nhf, nlf;
  mg -> GetNumEval(nhf, nlf);

  // Maximum number of samples.
  if (!res && (nhf >= Nmax)) res = true;

  return res;
}

// =============================== PrintPostVar ============================

void cMFSAO :: PrintPostVarSur(int maxgen, int opt, double evnum, cSURR *SurMod)
{
  // Print sample

  if (!out) return;
  cout << "entrou print post var sur" << endl;

  int nh, nl, nv, no;
  SurMod -> GetNumSample(nh);
  SurMod -> GetNumSampleLF(nl);
  SurMod -> GetNumVar(nv);
  SurMod -> GetNumOut(no);

  vector<cVector> sx;
  vector<cVector> sy;
  vector<cVector> sxlf;
  vector<cVector> sylf;

  SurMod -> GetSampleX(sx);
  SurMod -> GetSampleY(sy);
  SurMod -> GetSampleXLF(sxlf);
  SurMod -> GetSampleYLF(sylf);

  *out << "\n%SAMPLE.X.HF\n";

  for (int i = 0; i < nh; i++)
  {
      *out << i+1 << "  ";
      for (int j = 0; j < nv; j++)
          *out << sx[i][j] << "  ";
      *out << endl;
  }

  *out << "\n%SAMPLE.Y.HF\n";

  for (int i = 0; i < nh; i++)
  {
      *out << i+1 << "  ";
      for (int j = 0; j < no; j++)
          *out << sy[i][j] << "  ";
      *out << endl;
  }

  *out << "\n%SAMPLE.X.LF\n";

  for (int i = 0; i < nl; i++)
  {
      *out << i+1 << "  ";
      for (int j = 0; j < nv; j++)
          *out << sxlf[i][j] << "  ";
      *out << endl;
  }

  *out << "\n%SAMPLE.Y.LF\n";

  for (int i = 0; i < nl; i++)
  {
      *out << i+1 << "  ";
      for (int j = 0; j < no; j++)
          *out << sylf[i][j] << "  ";
      *out << endl;
  }
}

// ============================ EvalErrorMeasures ===============================

void cMFSAO :: EvalErrorMeasures(vector<cVector> &vsysur, cVector &nrmse, cVector &rmae,
                                cVector &nmae, cVector &error, int no)
{
    // cout << "ENTROU NAS METRICAS" << endl;

    // Evaluate NRMSE and RMAE

     nrmse.Resize(no);
     nmae.Resize(no);
     rmae.Resize(no);
     error.Resize(no);

     cVector yave(no);
     cVector maxi(no);

     for (int i = 0; i < no; i++)
     {
         double a, b, sum, sumsamples = 0.0;

         for (int j = 0; j < Nvs; j++)
         {
          //   cout << "Surr predict: " << (vsysur[j][i]) << endl;
          //   cout << "HFM model: " << (Vsy[j][i]) << endl;
             a += (Vsy[j][i] - vsysur[j][i])*(Vsy[j][i] - vsysur[j][i]);
             b += Vsy[j][i]*Vsy[j][i];
             sum += abs(((Vsy[j][i])-(vsysur[j][i]))/((Vsy[j][i])));
             maxi[i] = max((Vsy[j][i])-abs((vsysur[j][i])), maxi[i]);
             sumsamples += (Vsy[j][i]);
         }

         yave[i] = sumsamples/Nvs;

         rmae[i] = (1/Nvs)*sum;
         nrmse[i] = pow(a/b, 0.5);

      //   cout << "NRMSE " << nrmse[i] << endl;
     //    cout << "RMAE " << rmae[i] << endl;
      }

     // Evaluate NMAE

     if (Nvs > 1)
     {
         for (int i = 0; i < no; i++)
         {
             double c = 0.0;
             for (int j = 0; j < Nvs; j++)
             {
                 c += pow((Vsy[j][i]) - yave[i],2);
             }
             nmae[i] = maxi[i]/pow((1/Nvs)*c,0.5);
         }
     }
     else
     {
         for (int i = 0; i < no; i++) nmae[i] = 0;  // NULL
     }

     // Evaluate error

     for (int i = 0; i < 1; i++)
     {
         cout << "Best predict HFM: " << (vsysur[Nvs][i]) << endl;
         cout << "HFM model: " << (Vsy[0][i]) << endl;
         error[i] = abs(((vsysur[Nvs][i]) - (Vsy[0][i]))/(Vsy[0][i]));
         cout << "ERROR " << error[i] << endl;
     }
 }

// =============================== PostProcessingSur ====================================

void cMFSAO :: PostProcessingSur(int nrun, int numout, std::vector<cVector> *NRMSE,
                                        std::vector<cVector> *RMAE, std::vector<cVector> *NMAE,
                                        std::vector<cVector> *ERROR, std::vector<int> *GentoConvSur)
{
/*
     cout << "ENTROU AQUI NO FINALZAO" << endl;
    cVector STDNRMSE(numout);
    cVector STDRMAE(numout);
    cVector STDNMAE(numout);
    cVector ERRORave(numout);

    cout << "ENTROU AQUI NO FINALZAO" << endl;

    if (Prob->GetNumObj() == 1)         // Only available for single-objective problems
    {
        // Evaluate mean NRMSE, RMAE and NMAE

        cVector NRMSEave(numout);
        cVector RMAEave(numout);
        cVector NMAEave(numout);

        for (int i = 0; i < numout; i++)
        {
            double sumNRMSE, sumRMAE, sumNMAE, sumERROR = 0.0;
            for (int j = 0; j < nrun; j++)
            {
                sumNRMSE += (*NRMSE)[j][i];
                sumRMAE += (*RMAE)[j][i];
                sumNMAE += (*NMAE)[j][i];
                sumERROR += (*ERROR)[j][i];
            }
            NRMSEave[i] = sumNRMSE/(nrun);
            RMAEave[i] = sumRMAE/(nrun);
            NMAEave[i] = sumNMAE/(nrun);
            ERRORave[i] = sumERROR/(nrun);
        }

        // Evaluate the average generation convergence of surrogate model

        double sumGENconv = 0.0;

        for (int j = 0; j < nrun; j++) sumGENconv += (*GentoConvSur)[j];

        double GentoConvSurave = sumGENconv/nrun;

        // Evaluate standard deviations

        for (int i = 0; i < numout; i++)
        {
            double tempstdNRMSE, tempstdRMAE, tempstdNMAE = 0.0;
            for (int j = 0; j < nrun; j++)
            {
                tempstdNRMSE += pow((*NRMSE)[j][i] - NRMSEave[i],2);
                tempstdRMAE += pow((*RMAE)[j][i] - RMAEave[i],2);
                tempstdNMAE += pow((*NMAE)[j][i] - NMAEave[i],2);
            }

            STDNRMSE[i] = pow(tempstdNRMSE/(nrun-1),0.5);
            STDRMAE[i] = pow(tempstdRMAE/(nrun-1),0.5);
            STDNMAE[i] = pow(tempstdNMAE/(nrun-1), 0.5);
        }

        out << "\n%NRMSE.OPTIMIZATION\n";
        for (int i = 0; i < numout; i++)
        {
            out << i+1 << endl;
            for (int j = 0; j < nrun; j++) out << (*NRMSE)[j][i] << endl;
        }

        out << "\n%RMAE.OPTIMIZATION\n";
        for (int i = 0; i < numout; i++)
        {
            out << i+1 << endl;
            for (int j = 0; j < nrun; j++) out << (*RMAE)[j][i] << endl;
        }

        out << "\n%ERROR.OPTIMIZATION\n";
        for (int i = 0; i < numout; i++)
        {
            out << i+1 << endl;
            for (int j = 0; j < nrun; j++) out << (*ERROR)[j][i] << endl;
        }

        out << "\n%STANDARD.DEVIATION.NRMSE\n";
        for (int i = 0; i < numout; i++) out << i+1 << "  " << STDNRMSE[i] << endl;

        if (RMAE != NULL)
        {
        out << "\n%STANDARD.DEVIATION.RMAE\n";
        for (int i = 0; i < numout; i++) out << i+1 << "  " << STDRMAE[i] << endl;
        }

        if (NMAE != NULL)
        {
        out << "\n%STANDARD.DEVIATION.NMAE\n";
        for (int i = 0; i < numout; i++) out << i+1 << "  " << STDNMAE[i] << endl;
        }

        out << "\n%AVERAGE.ERROR.CONV.SURROGATE\n";
         for (int i = 0; i < numout; i++) out << i+1 << "  " << ERRORave[i] << endl;

        out << "\n%GENERATIONS.CONV.SURROGATE\n";
        for (int i = 0; i < nrun; i++) out << i+1 << "  " << (*GentoConvSur)[i] << endl;

        out << "\n%AVERAGE.HFM.EVAL.CONV.SURROGATE\n";
        out << GentoConvSurave << endl;

        out << "\n%TOTAL.HFM.EVALUATIONS\n";
        out << sumGENconv << endl;
    }
    */
}

// ============================== MaxExpectImprov ==========================

void cMFSAO :: EvalInfillCriteria(cSURR* Sur, cSampSet *smp, sProbAppOut &appout,cVector &xb, double &yb)
{
    cOptAlgorithm *alg = SubAlgType;
    cPenalty *pen = new cPenStatic;

    alg -> SetSolType(SolType);
    alg -> SetPopSize(SubPop);
    alg -> SetMaxGen(SubMaxGen);
    alg -> SetOptNum(1);
    alg -> SetFeedback(false);
    alg -> SetTolViol(SubTolViol);
    alg -> SetPenFunction(pen);
    alg -> SetMutProb(SubMutProb);

    alg -> SetStallGen(SubStallGen);
    alg -> SetSampType(SampType);
    alg -> SetIntPopSamp(1);

    if (Pen)
        pen -> SetFactor(Pen->GetFactor( ));

    if (alg->GetType() == STANDARD_PSO){
        alg -> SetSwarmTopology(SubTopology);
    }
    else if(alg->GetType() == STANDARD_DE){
        alg -> SetDifType(SubDifType);
    }

    // Get best feasible sample objective function.
    double bestobj = smp->BestSol( )->GetObjFunc(0);

    cout << "BestObj = " << bestobj << endl;

    cProbSurr* probbest = new cProbSurr(Sur,Prob, &appout,bestobj,InfillCriteria);
    //cProbSurr* probbest = new cProbSurr(Sur,Prob,&appout,bestobj,EVALUATE_SURROGATE);
    probbest -> SetConstrMethod(ConstrMethod); // LEO

    alg -> SetProblem(probbest);
    alg -> Init( );
    alg -> Solver( );
    cOptSolution* best = alg -> GetBest();

    best -> GetNormVar(xb);
    yb = best -> GetObjFunc(0);

    if (Feedback)
    {
      // TODO: Printar de acordo com o tipo de infill criteria escolido.
      cout << "EI max: " << yb << endl;
      cout << yb << endl;
      cout << "pntX: ";
      xb.Print( );
    }
}

// ======================================================= End of file =====
