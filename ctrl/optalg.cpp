// -------------------------------------------------------------------------
// optalg.cpp - Implementation of control class.
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
// Created:      21-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     13-Mar-2013    Evandro Parente Junior
//               Renamed to optalg.cpp.
//
//               30-Nov-2013    Elias Saraiva Barroso
//               Added PSO, ABC and AIS optimization methods.
//
//               10-Nov-2015    Elias Saraiva Barroso
//               Implemented a new structure with pre and post-processing,
//               treating new algorithms (subclasses) abstractly.
//
//               03-Mar-2019    Marina Alves Maia
//               Added NSGA and LamNSGA optimizations methods.
//               Modified CreateOptAlg in order to associate modNSGA and
//               modLamNSGA to the multiobjective problems with up to two criteria.
//
//               16-Nov-2020    Elias Saraiva Barroso
//               Added ReadSampType method from SAO class.
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "optalg.h"
#include "input.h"
#include "group.h"
#include "optsolution.h"
#include "sel.h"
#include "penalty.h"
#include "problem.h"
#include "saokrg.h"
#include "grego.h"
#include "stdga.h"
#include "lamga.h"
#include "modnsgaII.h"
#include "modlamnsgaII.h"
#include "stdpso.h"
#include "modpso.h"
#include "stdabc.h"
#include "stdais.h"
#include "stdde.h"
#include "kitayamasao.h"
#include "saorbf.h"
#include "rs.h"
#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
//

// -------------------------------------------------------------------------
// Class cOptAlgReadEntry:
//

// ================================ Read ===================================

void cOptAlgReadEntry :: Read(istream &in)
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
    alg = new cStandardGA( );
  else if (string(label)=="LamGA" || string(label)=="lamga")
    alg = new cLaminateGA;
  else if (string(label)=="modNSGAII" || string(label)=="modnsgaII" || string(label)=="modnsgaii")
    alg = new cmodNSGAII( );
  else if (string(label)=="modLamNSGAII" || string(label)=="modlamnsgaII" || string(label)=="modlamnsgaii")
    alg = new cmodLaminateNSGAII( );
  else if (string(label)=="StdPSO" || string(label)=="stdpso")
    alg = new cStandardPSO( );
  else if (string(label)=="LamPSO" || string(label)=="lampso")
    alg = new cLaminatePSO( );
  else if (string(label)=="StdABC" || string(label)=="stdabc")
    alg = new cStandardABC( );
  else if (string(label)=="StdAIS" || string(label)=="stdais")
    alg = new cStandardAIS( );
  else if (string(label)=="StdDE" || string(label)=="stdde")
    alg = new cStandardDE( );
  else if (string(label)=="SAOKRG" || string(label)=="saokrg")
    alg = new cSAOKRG( );
  else if (string(label)=="SAORBF" || string(label)=="saobrf")
    alg = new cSAORBF( );
  else if (string(label)=="RandSearch" || string(label)=="rs")
    alg = new cRandomSearch( );
  else
  {
      cout << "Unknown algorithm: " << label << endl;
      exit(0);
  }

  // Load algorithm read functions.
  alg->LoadReadFunc(*inpmap);
}

// -------------------------------------------------------------------------
// Public methods:
//
// ============================= ReadOptNum ================================

void cOptAlgorithm :: ReadOptNum(istream &in)
{
  if (!(in >> OptNum))
  {
    cout << "Error in the input of the number of optimizations." << endl;
    exit(0);
  }
}

// ============================== ReadTolViol ==============================

void cOptAlgorithm :: ReadTolViol(istream &in)
{
  if (!(in >> TolViol))
  {
    cout << "Error in the input of the constraint tolerance." << endl;
    exit(0);
  }
}

// ============================== ReadTolSucRate ===========================

void cOptAlgorithm :: ReadTolSucRate(istream &in)
{
  if (!(in >> TolSucRate))
  {
    cout << "Error in the input of the success optimization tolerance." << endl;
    exit(0);
  }
}

// ============================== ReadPopSize ==============================

void cOptAlgorithm :: ReadPopSize(istream &in)
{
  if (!(in >> PopSize))
  {
    cout << "Error in the input of the population size." << endl;
    exit(0);
  }
}

// ============================= ReadSamptype =============================

void cOptAlgorithm :: ReadSampType(istream &in)
{
  if (!(in >> SampType))
  {
    cout << "Error in the input of the sampling method." << endl;
    exit(0);
  }
  IntPopSamp = true;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadCrossMethod =============================

void cOptAlgorithm :: ReadCrossMethod(istream &in)
{
  char cross[100];

  if (!Utl::ReadString(in, cross))
  {
    cout << "Error in the input of the crossover method." << endl;
    exit(0);
  }

  if (Type == STANDARD_GA || Type == LAMINATE_GA || Type == modNSGAII || Type == modLAMINATE_NSGAII)
  {
      if(string(cross) == "LinearCombination" || string(cross) == "linearcombination")
       CrossType = LINEAR_COMBINATION;
      else if(string(cross) == "Classical" || string(cross) == "classical")
       CrossType = CLASSICAL;
      else if(string(cross) == "BinDE" || string(cross) == "Binomial" || string(cross) == "BIN" || string(cross) == "Bin")
       CrossType = BIN;
      else
      {
        cout << "Unknown crossover method: " << cross << endl;
        exit(0);
      }
  }
  else
  {
  }
}

// ============================== ReadDifType ===============================

void cOptAlgorithm :: ReadDifType(istream &in)
{
  char label[100];

  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of selection method." << endl;
    exit(0);
  }

  if (string(label) == "Rand1" || string(label) == "rand1" || string(label) == "Default" || string(label) == "def")
    DifType = Rand1;
  else if (string(label) == "LocalToBest" || string(label) == "Local2Best" || string(label) == "l2b" || string(label) == "local-to-best")
    DifType = Loc2Best;
  else if (string(label) == "BestWithJitter" || string(label) == "BestJitter" || string(label) == "BwJ" || string(label) == "best-with-jitter")
    DifType = BestJitter;
  else
  {
    cout << "Unknown differentiation method: " << label << endl;
    exit(0);
  }
}

// ============================== ReadMaxGen ===============================

void cOptAlgorithm :: ReadMaxGen(istream &in)
{
  if (!(in >> MaxGen))
  {
    cout << "Error in the input of the maximum number of generations." << endl;
    exit(0);
  }
}

// =========================== ReadStallGen ================================

void cOptAlgorithm :: ReadStallGen(istream &in)
{
  if (!(in >> StallGen))
  {
    cout << "Error in the input of the stallgen number." << endl;
    exit(0);
  }
}

// ============================== ReadMigrGap ==============================

void cOptAlgorithm :: ReadMigrGap(istream &in)
{
  if (!(in >> MigrationGap))
  {
    cout << "Error in the input of the migration generation gap." << endl;
    exit(0);
  }
}

// ============================== ReadMigrNum ==============================

void cOptAlgorithm :: ReadMigrNum(istream &in)
{
  if (!(in >> MigrationNum))
  {
    cout << "Error in the input of the number of individuals for migration." << endl;
    exit(0);
  }
}

// ============================== ReadMutProb ==============================

void cOptAlgorithm :: ReadMutProb(istream &in)
{
  if (!(in >> MutProb))
  {
    cout << "Error in the input of the mutation rate." << endl;
    exit(0);
  }
}

// ============================== ReadMutRange =============================

void cOptAlgorithm :: ReadMutRange(istream &in)
{
  if (!(in >> MinMut) || !(in >> MaxMut))
  {
    cout << "Error in the input of the mutation range." << endl;
    exit(0);
  }
}

// ============================== ReadMaxThread ============================

void cOptAlgorithm :: ReadMaxThread(istream &in)
{
  if (!(in >> omp_maxthread) || (omp_maxthread < 1))
  {
    cout << "Error in the input of the maximum number of threads." << endl;
    exit(0);
  }
}

// ============================ ReadMinObjFunc =============================

void cOptAlgorithm :: ReadMinObjFunc(istream &in)
{
  if (!(in >> MinObjFunc))
    {
      cout << "Error in the input of the minimum objective function." << endl;
      exit(0);
    }

  MinObjFlag = true;
}

// ============================== ReadSeed =================================

void cOptAlgorithm :: ReadSeed(istream &in)
{
  unsigned int s;
  if (!(in >> s))
  {
    cout << "Error in the input of the seed." << endl;
    exit(0);
  }
  Utl :: SetSeed(s);
}

// ============================= ReadNumBestSol ============================

void cOptAlgorithm :: ReadNumBestSol(istream &in)
{
  if (!(in >> NumBestSol) || (NumBestSol < 1))
  {
    cout << "Error in the input of the number of output best solutions." << endl;
    exit(0);
  }
}

// ============================== ReadProblem ==============================

void cOptAlgorithm :: ReadProblem(istream &in)
{
  // Read the problem label.

  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of the problem type label." << endl;
    exit(0);
  }

  string problabel(label);
  eVarType  vartype = cIndividual     :: GetVarType(SolType);
  Prob              = cProblemFactory :: CreateProblem(problabel,(bool) vartype);
}

// ============================= ReadSelMethod =============================

void cOptAlgorithm :: ReadSelMethod(istream &in)
{
  char label[100];
  eSelType InpType;

  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of selection method." << endl;
    exit(0);
  }

  if (string(label) == "Ranking" || string(label) == "ranking")
    InpType = RANKING;
  else if (string(label) == "FitnessProportional" ||
           string(label) == "fitnessproportional")
    InpType = FITNESS_PROPORTIONAL;
  else if (string(label) == "Tournament" ||
           string(label) == "Tournament")
    InpType = TOURNAMENT;
  else if (string(label) == "TournamentNSGAII" ||
             string(label) == "TournamentnsgaII" ||
             string(label) == "tournamentnsgaii")
    InpType = TOURNAMENTNSGAII;
  else
  {
    cout << "Unknown selection method: " << label << endl;
    exit(0);
  }

  // Create selection method.
  if (Sel) delete Sel;
  Sel = cSelection :: CreateSelection(InpType);
}

// ============================== ReadPenalty ==============================

void cOptAlgorithm :: ReadPenalty(istream &in)
{
  char label[100];
  ePenType InpType;

  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of penalty method." << endl;
    exit(0);
  }

  if (string(label) == "Static" || string(label) == "static")
    InpType = STATIC;
  else if (string(label) == "Deb2000" || string(label) == "deb2000" ||
           string(label) == "deb")
    InpType = DEB2000;
  else if (string(label) == "Adaptive" || string(label) == "adaptive")
    InpType = ADAPTIVE;
  else if (string(label) == "Normalization"  || string (label) == "normalization")
      InpType = NORMALIZATION;
  else
  {
    cout << "Unknown penalty method: " << label << endl;
    exit(0);
  }

  // Create penalty method.
  if (Pen) delete Pen;
  Pen = cPenalty :: CreatePenalty(InpType);
  Pen->LoadReadFunc(*InpMap);
}

// ============================= ReadOptSolType ===========================

void cOptAlgorithm :: ReadOptSolType(istream &in)
{
  if (!(in >> SolType))
  {
    cout << "Error in the input of the optimization solution type label." << endl;
    exit(0);
  }
}

// ============================= ReadInpSol ================================

void cOptAlgorithm :: ReadInpSol(istream &in)
{
  eSolType type;
  int varsize;
  if (!(in >> NumInpSol) || !(in >> type) || !(in >> varsize))
  {
    cout << "Error in the input of the number of input solutions." << endl;
    exit(0);
  }

  if (NumInpSol < 0 || NumInpSol > PopSize)
  {
    cout << "Error: Invalid value for the number of input solutions." << endl;
    exit(0);
  }

  if (varsize < 0)
  {
    cout << "Error in the input of the variable size of the input solutions." << endl;
    exit(0);
  }
  InpSolVec = new sInpSol [NumInpSol];
  for(int i = 0; i < NumInpSol; ++i)
  {
    InpSolVec[i].type = type;
    InpSolVec[i].CodVar.Resize(varsize);

    // Read each solution variable.
    for(int j = 0; j < varsize; ++j)
      if (!(in >> InpSolVec[i].CodVar[j]))
      {
        cout << "Error in the input of the variable of the input solutions." << endl;
        exit(0);
      }
  }
}

// ================================= End ===================================

void cOptAlgorithm :: End(void)
{
  if (!out) return;

  // Print the end mark and the step labels.

  *out << "\n%END\n";
}

// =============================== LoadReadFunc ============================

void cOptAlgorithm :: LoadReadFunc(cInpMap &im)
{
  // Store input map.
  InpMap = &im;

  // Register read functions.
  im.Insert("OPTIMIZATION.NUMBER",        makeReadObj(cOptAlgorithm,ReadOptNum));
  im.Insert("POPULATION.SIZE",            makeReadObj(cOptAlgorithm,ReadPopSize));
  im.Insert("SWARM.SIZE",                 makeReadObj(cOptAlgorithm,ReadPopSize));
  im.Insert("CROSSOVER.METHOD",           makeReadObj(cOptAlgorithm,ReadCrossMethod));
  im.Insert("DIFFERENTIATION.METHOD",     makeReadObj(cOptAlgorithm,ReadDifType));
  im.Insert("MAXIMUM.GENERATIONS",        makeReadObj(cOptAlgorithm,ReadMaxGen));
  im.Insert("MAXIMUM.ITERATIONS",          makeReadObj(cOptAlgorithm,ReadMaxGen));
  im.Insert("STALL.GEN",                  makeReadObj(cOptAlgorithm,ReadStallGen));
  im.Insert("MIGRATION.GENERATION.GAP",   makeReadObj(cOptAlgorithm,ReadMigrGap));
  im.Insert("MIGRATION.INDIVIDUAL.NUMBER",makeReadObj(cOptAlgorithm,ReadMigrNum));
  im.Insert("MUTATION.PROBABILITY",       makeReadObj(cOptAlgorithm,ReadMutProb));
  im.Insert("MUTATION.RANGE",             makeReadObj(cOptAlgorithm,ReadMutRange));
  im.Insert("CONSTRAINT.TOLERANCE",       makeReadObj(cOptAlgorithm,ReadTolViol));
  im.Insert("SUCCESS.RATE.TOLERANCE",     makeReadObj(cOptAlgorithm,ReadTolSucRate));
  im.Insert("MAXIMUM.THREAD.NUMBER",      makeReadObj(cOptAlgorithm,ReadMaxThread));
  im.Insert("SEED",                       makeReadObj(cOptAlgorithm,ReadSeed));
  im.Insert("OUTPUT.SOLUTION.NUMBER",     makeReadObj(cOptAlgorithm,ReadNumBestSol));
  im.Insert("SELECTION.METHOD",           makeReadObj(cOptAlgorithm,ReadSelMethod));
  im.Insert("PENALTY.METHOD",             makeReadObj(cOptAlgorithm,ReadPenalty));
  im.Insert("MINIMUM.OBJECTIVE.FUNCTION", makeReadObj(cOptAlgorithm,ReadMinObjFunc));
  im.Insert("INDIVIDUAL.TYPE",            makeReadObj(cOptAlgorithm,ReadOptSolType));
  im.Insert("INPUT.SOLUTIONS",            makeReadObj(cOptAlgorithm,ReadInpSol));
  im.Insert("PROBLEM.TYPE",               makeReadObj(cOptAlgorithm,ReadProblem));
  im.Insert("SAMPLING.METHOD" ,           makeReadObj(cOptAlgorithm,ReadSampType));
}

// ============================== CreateOptAlg =============================

cOptAlgorithm* cOptAlgorithm :: CreateOptAlg(eOptAlgType type, cProblem* prob)
{
  cOptAlgorithm *alg = 0;
  eOptAlgType check;
  eOptAlgType check2;

  check = modNSGAII;
  check2 = modLAMINATE_NSGAII;

      switch(type)
      {
      case (STANDARD_GA):
        alg = new cStandardGA( );
      break;

      case (LAMINATE_GA):
        alg = new cLaminateGA( );
      break;

      case (STANDARD_PSO):
        alg = new cStandardPSO( );
      break;

      case (STANDARD_DE):
        alg = new cStandardDE( );
      break;

      case (LAMINATE_PSO):
        alg = new cLaminatePSO( );
      break;

      case (STANDARD_ABC):
        alg = new cStandardABC( );
      break;

      case (STANDARD_AIS):
        alg = new cStandardAIS( );
      break;

      case (modNSGAII):
         alg = new cmodNSGAII( );
      break;

      case (modLAMINATE_NSGAII):
         alg = new cmodLaminateNSGAII( );
      break;

      case (SAOKRG):
          alg = new cSAOKRG( );
      break;

//      case (KITSAO):
//'          alg = new cKitayamaSAO( );
//      break;

      case (SAORBF):
          alg = new cSAORBF( );
      break;
      }

      if (prob->GetNumObj( ) == 2)
      {
          if (type != check && type != check2)
          {
                  cout << "ERROR! The chosen algorithm is not supported for multiobjective opt. "
                          "The problem will be solved by the modNSGAII algorithm." << endl;
                  alg = new cmodNSGAII( );
          }
      }
      else if (prob->GetNumObj( ) > 2)
      {
          cout << "ERROR! The maximum number allowed is 2. Redefine your objective functions." << endl;
          exit(0);
      }
      else if (prob->GetNumObj( ) == 1)
      {
          if (type == check || type == check2)
          {
              cout << "ERROR! The chosen algorithm is not supported for single-objective optitmization. "
                      "The problem will be solved by the GA algorithm." << endl;
              alg = new cStandardGA( );
          }
      }

  return(alg);
}

// =============================== cControl ================================

cOptAlgorithm :: cOptAlgorithm(void)
{
  // Set default parameters.
  OptNum       = 1;
  TolViol      = 1.0e-3;
  TolSucRate   = 1.0e-4;
  PopSize      = 100;
  CrossType    = LINEAR_COMBINATION;
  cont         = 10;
  MaxGen       = 200;
  StallGen     = 200;
  MigrationGap = 1e6;
  MigrationNum = 0;
  MutProb      = 0.00;
  MaxMut       = 0.00;
  MinMut       = 0.00;
  NumBestSol   = 5;
  SampType     = NLHS;
  IntPopSamp   = false;
  out          = 0;

  MinObjFunc   = 0;
  MinObjFlag   = 0;


  // solution paameters.
  SolType      = SOL_DBL_VEC;
  NumInpSol    = 0;
  InpSolVec    = 0;

  Prob         = 0;
  Pen          = 0;
  Sel          = 0;

}

// =============================== ~cControl ===============================

cOptAlgorithm :: ~cOptAlgorithm(void)
{
}

// =============================== Init ====================================

void cOptAlgorithm :: Init(void)
{
  // Create selection method if needed.

  if (!Sel) Sel = cSelection :: CreateSelection(RANKING);

  // Create penalty method if needed.

  if (!Pen && Prob->GetNumConstr( ) > 0)
    Pen = cPenalty :: CreatePenalty(STATIC);

  if (Pen && Prob->GetNumConstr( ) == 0)
  {
    delete Pen;
    Pen = 0;
  }

  // Initialize problem.
  Prob->Init( );

  // Track the number of generations until convergence.

  if (Prob->GetNumObj( ) == 1)
  {
      GentoConv = new int[OptNum];
      GenBest   = new double[MaxGen];
      GenBestConstraints.Resize(MaxGen, Prob->GetNumConstr());
      GenBestConstraints.Zero();
     /* GenVarBest.Resize(MaxGen, Prob->GetNumVar());
      GenVarBest.Zero();*/
      GenAvg    = new double[MaxGen];
      MBestGen  = new double [MaxGen];
      best      = new cSolGroup(OptNum,Prob);

  // Initialize mean best solutions

  for(int gen = 0; gen < MaxGen; gen++) MBestGen[gen] = 0;
  }
}

// =============================== GetBestIndividuals ===========================

void cOptAlgorithm :: GetBestIndividuals(int gen, cVector &BestObjGen, cMatrix &BestIndGen, cGroup *mg)
{
    // Find the best individuals in a generation

    cVector aux(Prob->GetNumVar());

    aux  = mg->BestSol()->GetVec();

    for (int i = 0; i < Prob->GetNumVar(); i++)
    {
        BestIndGen[gen][i] = aux[i];
    }

    if (Pen)
    {
      BestObjGen[gen] = mg -> BestSol( )->GetPenObjFunc( );
    }
    else
    {
      BestObjGen[gen] = mg -> BestSol( )->GetObjFunc(0);
    }
}

// =============================== UpdatePostVar ===========================

void cOptAlgorithm :: UpdatePostVar(int gen, int opt, double &lb, cGroup *mg)
{
  if (Prob->GetNumObj( ) == 1) // Only available for single-objective problems
  {
  if (Pen)
  {
    GenBest[gen] = mg->BestSol( )->GetPenObjFunc( );
    cVector constr = mg->BestSol()->GetConstr();

    for (int i = 0; i < Prob->GetNumConstr( ); i++)
    {
        GenBestConstraints[gen][i] = constr[i];
    }
  }
  else
  {
    GenBest[gen] = mg->BestSol( )->GetObjFunc(0);
  }

  /* cVector aux(Prob->GetNumVar());

   aux  = mg->BestSol()->GetVec();

    for (int i = 0; i < Prob->GetNumVar(); i++)
    {
        GenVarBest[gen][i] = aux[i];
    }*/


  GenAvg[gen]  = mg->AvgObjFunc( );
  MBestGen[gen] += GenBest[gen];

  if (gen == 0)
  {
    lb = GenBest[gen];
    GentoConv[opt] = 1;
    StallGenCount = 0;
  }

  else if ( abs(GenBest[gen] - lb) > 1.0e-5)
  {
    GentoConv[opt] = gen + 1;
    lb = GenBest[gen];
    StallGenCount = 0;
  }

  else
    StallGenCount++;
  }
}

// =============================== OptStopCrit =============================

bool cOptAlgorithm :: OptStopCrit(int gen, int opt, double &lb, cGroup *mg)
{
  // Only available for single-objective problems
  // Multiobjective problems use the maximum generations as stopping criterion

  // Verify conditions to stop optimization

  bool StopCrit = false;

  if (Prob->GetNumObj( ) == 1)
  {

  // Condition 01 - max number of generations without improvement

  if (StallGenCount >= StallGen)
  {
      //cout << "Stall gen" << endl;
    StopCrit = true;
  }

  // Condition 02 - Reached the mínimum obj function

  if (MinObjFlag)
  {
    if (fabs(MinObjFunc - lb) <= TolSucRate)
      StopCrit = true;
  }

  // Fill the incomplete output data with the last values obtained

  if (StopCrit)
  {
    for(int newgen = gen + 1; newgen < MaxGen; newgen++)
      UpdatePostVar(newgen,opt,lb,mg);
  }
  }
  else
  {
      cout << "Maximum generations is currently being used as stopping criterion. " << endl;
  }

  return StopCrit;
}

// =============================== PrintPostVar ============================

void cOptAlgorithm :: PrintPostVar(int maxgen, int opt, double evnum, cGroup *mg)
{
  if (!out) return;

  cout << endl;

  *out << "\n%RESULT.BEST.INDIVIDUALS\n";

  for (int i = 0; i < maxgen; i++)
    *out << i+1 << "  " << GenBest[i] << endl;

  *out << "\n%RESULT.AVERAGE.INDIVIDUALS\n";

  for (int i = 0; i < maxgen; i++)
    *out << i+1 << "  " << GenAvg[i] << endl;

  // Store and print the best individual.

  if (Feedback)
  {
    cout << "Optimization Number: " << opt+1 << endl;
    best->GetSol(opt)->Print( );
    cout << "Number of individual evaluations: " << evnum << endl;
    cout << "Generations until convergence: " << GentoConv[opt] << endl;
  }

  *out << "\n%RESULT.INDIVIDUAL.EVALUATIONS\n";
  *out << evnum << "\n";

  // Write the best solutions 

  cSolGroup OptBestSol(NumBestSol,Prob);
  mg->Sort( );

  for (int sol_id = 0; sol_id < mg->GetSize( ); sol_id++)
  {
    if (sol_id == 0)
      OptBestSol.Insert(mg->GetSol(sol_id));
    else
    {
      bool AlreadyIn = false;
      for(int i = 0; i < OptBestSol.GetSize( ); i++)
        AlreadyIn = OptBestSol.GetSol(i)->CompVar(mg->GetSol(sol_id));

      if (AlreadyIn)
        continue;
      else
        OptBestSol.Insert(mg->GetSol(sol_id));
    }

    if (OptBestSol.GetSize( ) >= NumBestSol)
      break;
  }

  *out << "\n%RESULT.BEST.INDIVIDUALS\n";
  *out << OptBestSol.GetSize( ) << "\n\n";
  for(int i = 0; i < OptBestSol.GetSize( ); i++)
  {
    OptBestSol.GetSol(i)->Write(*out);
    OptBestSol.GetSol(i)->Print( );
    }
 // best->GetSol(opt)->Write( );
}

// =============================== PrintPostVar ============================

void cOptAlgorithm :: PrintPostVar(int maxgen, int opt, double evnum, cGroup *mg, std::vector<int> *SolRank)
{
  if (!out) return;

  cout << endl;

  if (Feedback)
  {
     cout << "Optimization Number: " << opt+1 << endl;
     for (int i = 0; i < mg->GetSize(); i++)
     {
        if ((*SolRank)[i] == 0)
        {
           cout << "\nSolution " << i << endl;
           mg->GetSol(i)->Print( );
           cout << "Number of individual evaluations: " << evnum << endl;
        }
      }
   }

   *out << "\n%GENERATION\n";
   *out << maxgen << endl;

   *out << "\n%RESULT.FEASIBLE.INDIVIDUALS\n";
   int solrankzero = 0;
   for (int k = 0; k < mg->GetSize(); k++)
   {
       if (mg->GetSol(k)->GetNormConst( ) == 0)
       {
         *out << (*SolRank)[k] << "   " << (mg->GetSol(k)->GetObjFunc(0)) << "   " << mg->GetSol(k)->GetObjFunc(1) << endl;
         if ((*SolRank)[k] == 0)
         {
           solrankzero += 1;
         }
       }
   }

   *out << "\n%RESULT.INDIVIDUALS.RANK.ZERO\n";
   *out << solrankzero << endl;

   *out << "\n%RESULT.INDIVIDUALS.PARETOSFRONT\n";
   for(int j = 0; j < mg->GetSize(); j++)
   {
     if ((*SolRank)[j] == 0 && mg->GetSol(j)->GetNormConst( ) == 0)
     {
       *out << "Solution " << j << endl;
       mg->GetSol(j)->Write(*out);
     }
   }

  /* *out << "\n%RESULT.INDIVIDUALS.PARETOSFRONT.CONSTR\n";
   for (int i = 0; i < Prob->GetNumConstr(); i++)
   {
     *out << "CONSTRAINT " << i << endl;
       for (int j = 0; j < PopSize; j++)
       {
          if ((*SolRank)[j] == 0 && mg->GetSol(j)->GetNormConst( ) == 0)
          {
            cVector constr = mg->GetSol(j)->GetConstr();
            *out << constr[i] << endl;
          }
       }
   }*/
}    

// =============================== GetBest ====================================

cOptSolution* cOptAlgorithm :: GetBest(void)
{
  cOptSolution *allBest = best->BestSol( );
  return allBest;
}

// =============================== PostProcessing ====================================

void cOptAlgorithm :: PostProcessing(void)
{
  if ( Prob->GetNumObj() == 1)
  {
  // Evaluate the success rate and print the best individual.

  cOptSolution *allBest = best->BestSol( );

  // Determine the success rate.

  int sucrate = 0;
  double diff,num;
  for (int i = 0; i < OptNum; i++)
  {
    if (Pen)
    {
      diff = best->GetSol(i)->GetPenObjFunc( ) - allBest->GetPenObjFunc( );
      num  = abs(allBest->GetPenObjFunc( ));
    }
    else
    {
      diff = best->GetSol(i)->GetObjFunc(0) - allBest->GetObjFunc(0);
      num  = abs(allBest->GetObjFunc(0));
    }
    num   = MAX(1.0, num);
    diff /= num;

    if (abs(diff) <= TolSucRate) sucrate++;
  }

  // Evalute the average number of iterations to convergence.

  double avgGentoConv = 0.0;
  for (int i = 0; i < OptNum; i++) avgGentoConv += GentoConv[i];
  avgGentoConv /= OptNum;

  // Evaluate Mean best and Standard Deviation.

  double MeanObjFunc = best->AvgObjFunc( );
  double StdDev = best->StdDevObjFunc( );

  double NMSRE;
  if (MinObjFlag) NMSRE = best->NormMeanSqrRootErr(MinObjFunc); 

  if (Feedback)
  { 
    cout << "Best solution in " << OptNum << " optimizations:" << endl;
    allBest->Print( );
    cout << "Success Rate is: " << double(100*sucrate/OptNum) << "%" << endl;
    cout << "Average Generations to Convergence: " << avgGentoConv << endl;
    cout << "Standard Deviation Best Individuals: " << StdDev << endl;
    cout << "Mean Best: " << MeanObjFunc << endl;
    if (MinObjFlag) cout << "Normalized Mean Square Root Error (Best Individuals): " << NMSRE << endl;
  }

  // Write mean best solutions results
  if (!out) return;

  *out << "\n%MEAN.BEST.OBJECTIVE.FUNCTION\n";

  for (int i = 0; i < MaxGen; i++)
    *out << i+1 << "  " << MBestGen[i]/OptNum << endl;

  *out << "\n%RESULT.MEAN.BEST\n";
  *out << MeanObjFunc << "\n";
  *out << "\n%RESULT.STANDARD.DEVIATION.BEST\n";
  *out << StdDev << "\n";

  if (MinObjFlag) 
    *out << "\n%RESULT.NORMALIZED.MEAN.SQUARE.ROOT.ERROR\n" << NMSRE << endl;

  *out << "\n%RESULT.BEST.INDIVIDUAL\n";
  allBest->Write(*out);
  *out << "\n%RESULT.SUCCESS.RATE\n";
  *out << double(100*sucrate)/OptNum << "\n";

  // Delete some variables.

  delete []GentoConv;
  delete []GenBest;
  delete []GenAvg;
  delete []MBestGen;
  delete best;
  }
}

// ======================================================= End of file =====
