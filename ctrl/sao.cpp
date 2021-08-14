// -------------------------------------------------------------------------
// sao.cpp - implementation of cSAO class.
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
// Created:      30-Jul-2019    Leonardo Gonçalves Ribeiro
//
// Modified:     16-Nov-2020    Elias Saraiva Barroso
//               Moved ReadSampType to OptAlg.
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
// Auxiliary types:
//


sProbAppOut :: sProbAppOut( )
{
}

sProbAppOut :: sProbAppOut(cProblem* prob)
{
  // Get number of objective functions and constraints.
  numobjf  = prob->GetNumObj( );
  numconst = prob->GetNumConstr( );

  // Get approximated objective function data.
  ApproxOF = (numobjf)   ? new bool[numobjf]   : 0; 
  prob->GetApproxObj(ApproxOF);

  numappobjf = 0;
  for(int f = 0; f < numobjf; ++f) if (ApproxOF[f]) numappobjf++;

  // Get approximated constraints data.
  ApproxC  = (numconst) ? new bool[numconst] : 0; 
  prob->GetApproxConstr(ApproxC);

  numappconst = 0;
  for(int c = 0; c < numconst; ++c) if (ApproxC[c]) numappconst++;
}

sProbAppOut :: ~sProbAppOut( )
{
  if (ApproxOF) delete [] ApproxOF;
  if (ApproxC) delete  [] ApproxC;
}

int sProbAppOut :: GetAppConstrID(const int &id) const
{
  for(int i = 0, c = 0; c < numconst; ++c)
    if (ApproxC[c])
      if (i++ == id) return c;
}

int sProbAppOut :: GetExactConstrID(const int &id) const
{
  for(int i = 0, c = 0; c < numconst; ++c)
    if (!ApproxC[c])
      if (i++ == id) return c;
}

int sProbAppOut :: GetAppObjFuncID(const int &id)  const
{
  for(int i = 0, f = 0; f < numobjf; ++f)
    if (ApproxOF[f])
      if (i++ == id) return f;
}

int sProbAppOut :: GetAppConstrOutID(const int &id) const
{
  return id+numappobjf;
}

// -------------------------------------------------------------------------
// Class SAO:
//
// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cSAO ==============================

cSAO :: cSAO(void) : cOptAlgorithm( )
{
    SubPop         = 250;//100;
    SubMaxGen      = 100;//250;
    SubTolViol     = 1e-5;
    SubMutProb     = 0.02;
    Nmax           = 10e30;
    SubAlgType     = new cStandardPSO;
    MinNRMSE       = 0.01;
    SampleFileName = fname;
    InfillCriteria = EVALUATE_EXPECTED_IMPROVEMENT;

    WEI  = 0.5;
    Beta = 1.0;

    InputNS     = 0;
    FlagVS      = 0;
    Nvs         = 0;
}

// =============================== LoadReadFunc ============================

void cSAO :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("SUB.POPULATION.SIZE" ,  makeReadObj(cSAO, ReadSubPop));
  im.Insert("SUB.MAXIMUM.GENERATIONS" ,  makeReadObj(cSAO, ReadSubMaxGen));
  im.Insert("SUB.CONSTRAINT.TOLERANCE" ,  makeReadObj(cSAO, ReadSubTolViol));
  im.Insert("SUB.PSO.TOPOLOGY" ,  makeReadObj(cSAO, ReadSubPSOTopology));
  im.Insert("SUB.DE.TYPE" ,  makeReadObj(cSAO, ReadSubDEType));
  im.Insert("SUB.MUTATION.PROBABILITY" ,  makeReadObj(cSAO, ReadSubMutProb));
  im.Insert("MAXIMUM.NUMBER.OF.POINTS" ,  makeReadObj(cSAO, ReadNmax));
  im.Insert("MINIMUM.NRMSE" ,  makeReadObj(cSAO, ReadMinNRMSE));
  im.Insert("SUB.OPTIMIZATION.ALGORITHM" ,  makeReadObj(cSAO, ReadSubAlgType));
  im.Insert("NUMBER.OF.INITIAL.SAMPLING.POINTS" ,  makeReadObj(cSAO, ReadNumInitSamplingPoints));
  im.Insert("SAMPLE.FILE.NAME" ,  makeReadObj(cSAO, ReadSampleFileName));
  im.Insert("VALIDATION.SAMPLES" ,  makeReadObj(cSAO, ReadValSamples));
  im.Insert("CONSTRAINT.HANDLING.METHOD" ,  makeReadObj(cSAO, ReadConstrMethod));
  im.Insert("INFILL.CRITERIA" ,  makeReadObj(cSAO, ReadInfillCriteria));
  im.Insert("WEI.VALUE" ,  makeReadObj(cSAO, ReadWEI));
  im.Insert("BETA.VALUE" ,  makeReadObj(cSAO, ReadBeta));
}

// =========================== ReadSubPop ===========================

void cSAO :: ReadSubPop(istream &in)
{
  if (!(in >> SubPop))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ========================= ReadSubMaxGen ==========================

void cSAO :: ReadSubMaxGen(istream &in)
{
  if (!(in >> SubMaxGen))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ======================== ReadSubTolViol ==========================

void cSAO :: ReadSubTolViol(istream &in)
{
  if (!(in >> SubTolViol))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ======================== ReadSubMutProb ==========================

void cSAO :: ReadSubMutProb(istream &in)
{
  if (!(in >> SubMutProb))
  {
    cout << "Error in the input of the subproblem population size." << endl;
    exit(0);
  }
}

// ======================== ReadWEI ==========================

void cSAO :: ReadWEI(istream &in)
{
  if (!(in >> WEI))
  {
    cout << "Error in the input of the WEI weigth." << endl;
    exit(0);
  }
}

// ======================== ReadBeta ==========================

void cSAO :: ReadBeta(istream &in)
{
  if (!(in >> Beta))
  {
    cout << "Error in the input of the beta factor for the Lower Confidence Bound criterion." << endl;
    exit(0);
  }
}

// ======================== ReadSampleFileName ==========================

void cSAO :: ReadSampleFileName(istream &in)
{
  if (!(in >> SampleFileName))
  {
    cout << "Error in the input of the sample file name." << endl;
    exit(0);
  }
  SampleFileName = SampleFileName.substr(1, SampleFileName.size() - 2);
}

// ======================== ReadSubAlgType ==========================

void cSAO :: ReadSubAlgType(istream &in)
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

void cSAO :: ReadSubPSOTopology(istream &in)
{
    // Read the PSO topology label.

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

void cSAO :: ReadConstrMethod(istream &in)
{
    // Read the constrained method label.

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
    else
    {
        cout << "Unknown constraint handling method: " << label << endl;
        exit(0);
    }
}

// ======================== ReadInfillCriteria ==========================

void cSAO :: ReadInfillCriteria(istream &in)
{
    // Read the infill criteria label.

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
    else
    {
        cout << "Unknown infill criteria: " << label << endl;
        exit(0);
    }
}

// ======================== ReadSubDEType ==========================

void cSAO :: ReadSubDEType(istream &in)
{
    // Read the Differentiation label.

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

void cSAO :: ReadNmax(istream &in)
{
  if (!(in >> Nmax))
  {
    cout << "Error in the input of the crossover rate." << endl;
    exit(0);
  }
}

// ======================== ReadSubMinNRMSE =========================

void cSAO :: ReadMinNRMSE(istream &in)
{
  if (!(in >> MinNRMSE))
  {
    cout << "Error in the input of the minimum NRMSE." << endl;
    exit(0);
  }
}

// ======================== ReadNumInitSamplingPoints =========================

void cSAO :: ReadNumInitSamplingPoints(istream &in)
{
  InputNS = 1;
  if (!(in >> NumInitSP))
  {
    cout << "Error in the input of the number of initial sampling points." << endl;
    exit(0);
  }
}

// ======================== ReadValSamples =========================

void cSAO :: ReadValSamples(istream &in)
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

void cSAO :: SetApproxObj( )
{

}

// ======================== SetApproxConstr =========================

void cSAO :: SetApproxConstr( )
{

}

// ======================== SetInitialSample =========================

void cSAO :: SetInitialSample(sProbAppOut &appout, cSampSet* &set,sSampData &sdata, int &ev)
{
    cout << "\n\nSample file not found!\n\n";
    cout << "Generating sample..." << endl;
   
    int nv = Prob -> VarNumEff( ); //Prob->VarNumRow( ) * Prob->VarNumCol( );
    int ns = (InputNS) ? NumInitSP : 1.5*(nv + 1)*(nv + 2)/2.0;
    set    = new cSampSet(ns,SolType,Prob,appout);

    // Create sample points.
    cout << "nv: " << nv << endl;
    cout << "ns: " << ns << endl;
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
        InitialSample.InitSample(SampType, nv, ns, sx);

        // Copy sample points to sample set
        for (int i = 0; i < ns; i++)
        {
          (*set)[i]->Init(sx[i]);
        }

        // Check if any two points are equal
        for (int m = 0; m < ns; m++)
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
            cout << "Please, change the sampling technique or lower the number of initial sampling points.\n" << endl;
            exit(0);
        }
    }

    // Evaluate each sampling point
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < ns; i++)
    {
      (*set)[i]->Evaluate( );
    }
    ev+=ns;
   
    // Copy sample points to surrogate model, normalized and with correct number
    // of output (only approximated fobj + constraints). 
    
    sdata.NumVar    = nv;
    sdata.NumSample = ns; 
    sdata.NumOut    = appout.GetNumAppOut( );

    sdata.SampleX.resize(ns);
    sdata.SampleY.resize(ns);
    for (int i = 0; i < ns; i++)
    {
      // Input data.
      sdata.SampleX[i].Resize(nv);
      (*set)[i]->GetNormVar(sdata.SampleX[i]);

      // Output data.
      sdata.SampleY[i].Resize(sdata.NumOut);
      (*set)[i]->GetSurrOutRes(sdata.SampleY[i]);
    }

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
}

// ======================== ReadSampleFile =========================

void cSAO :: ReadSampleFile(ifstream &finp, int &ns, int &nv, int &no, cVector &low, cVector &upp,
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

void cSAO :: Solver(void)
{
  cSampSAO* newsmp;

  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Track number of individual evaluations.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Create the population, mating pool and parent array.
    cSampSet    *smp;
    sProbAppOut appout(Prob);
    sSampData sdata;
    SetInitialSample(appout,smp,sdata,EvalNum);

    // Evaluate penalized objective function in samples.
    // Note: Pobj is not considered in infill procedure, but is used to select
    // the best feasible sample. It is adopted to handle situtations where no
    // feasible sample is avaliable.

    if (Pen)
      Pen->EvalPenObjFunc(smp, TolViol);

    // Create the surrogate model.
    cSURR *SurModel = CreateSurrogate(sdata);

    SurModel -> SetWEI(WEI);
    SurModel -> SetBeta(Beta);

    // Evaluate initial sample points.

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    // Perform the GA iterations.

    for (int step = 0; step < MaxGen; step++)
    {
      if ((step+1)%1 == 0 && Feedback) cout << "Step: " << step + 1 << endl;

      if (ciclewei)
      {
        double w[3] = {0.2, 0.35, 0.5};
        double b[3] = {1.0, 2.0, 3.0};
        int stepw = step % 3;
        SurModel -> SetWEI(w[stepw]);
        SurModel -> SetBeta(b[stepw]);
        if (Feedback) cout << "w = " << w[stepw] << endl;
      }

      // Select new points.
      cVector  pntx(SurModel->GetSampData( ).NumVar);
      cVector npntx(SurModel->GetSampData( ).NumVar);
      cVector  pnty(SurModel->GetSampData( ).NumOut);
      double  objf;
      EvalInfillCriteria(SurModel,smp,appout,pntx,objf);


      // Evaluate new point.
      newsmp = smp->PushBack(pntx);
      newsmp->Evaluate( );
      newsmp->GetSurrOutRes(pnty);
      EvalNum++;

      // Update surrogate model.

      cVectorVec nx, ny;
      nx.push_back(pntx);
      ny.push_back(pnty);

      UpdateSurrogate(nx,ny);

      // Update variables related to PostProcessing.

      if (Pen) Pen->EvalPenObjFunc(smp, TolViol);

      UpdatePostVar(step, opt, lastBest,smp); 
      
      if (Feedback) 
      {
        cout << "\nNew Sample: " << endl;
        newsmp->Print( );
        cout << "Number of samples: " << smp->GetSize( ) << endl;
        cout << "Best Sample: " << endl;
        smp->BestSol( )->Print( );
        cout << endl;
      }
    
      // Check conditions to stop optimization

      if (OptStopCrit(step, opt, lastBest, smp))
        break;   
    }

    // Store the best individual.
    best->Insert(smp->BestSol( ));

    // Print data in the output file.

    PrintPostVar(MaxGen, opt, EvalNum, smp);
  }
}

// =============================== OptStopCrit =============================

bool cSAO :: OptStopCrit(int gen, int opt, double &lb, cGroup *mg)
{
  // Call standard optimization stopping criteria
  bool res = cOptAlgorithm :: OptStopCrit(gen,opt,lb,mg);

  // Maximum number of samples.
  if (!res && (mg->GetSize( ) >= Nmax)) res = true;

  return res;
}

// ============================ EvalErrorMeasures ===============================

void cSAO :: EvalErrorMeasures(vector<cVector> &vsysur, cVector &nrmse, cVector &rmae,
                                cVector &nmae, cVector &error, int no)
{
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
             a += (Vsy[j][i] - vsysur[j][i])*(Vsy[j][i] - vsysur[j][i]);
             b += Vsy[j][i]*Vsy[j][i];
             sum += abs(((Vsy[j][i])-(vsysur[j][i]))/((Vsy[j][i])));
             maxi[i] = max((Vsy[j][i])-abs((vsysur[j][i])), maxi[i]);
             sumsamples += (Vsy[j][i]);
         }

         yave[i] = sumsamples/Nvs;

         rmae[i] = (1/Nvs)*sum;
         nrmse[i] = pow(a/b, 0.5);
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

// ============================== MaxExpectImprov ==========================

void cSAO :: EvalInfillCriteria(cSURR* Sur, cSampSet *smp, sProbAppOut &appout,cVector &xb, double &yb)
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

    if (alg->GetType() == STANDARD_PSO){
        alg -> SetSwarmTopology(SubTopology);
    }
    else if(alg->GetType() == STANDARD_DE){
        alg -> SetDifType(SubDifType);
    }

    // Get best feasible sample objective function.
    double bestobj = smp->BestSol( )->GetObjFunc(0);

    cProbSurr* probbest = new cProbSurr(Sur,Prob, &appout,bestobj,InfillCriteria);

    probbest -> SetConstrMethod(ConstrMethod);

    alg -> SetProblem(probbest);
    alg -> Init( );
    alg -> Solver( );
    cOptSolution* best = alg -> GetBest();

    best -> GetNormVar(xb);
    yb = best -> GetObjFunc(0);

    if (Feedback)
    {
      cout << "Best acquisition function value: " << yb << endl;
    }
}

// ======================================================= End of file =====
