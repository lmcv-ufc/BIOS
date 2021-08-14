// -------------------------------------------------------------------------
// metaopt.cpp - implementation of the benchmark class.
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
// Created:      16-Oct-2014    Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "metaopt.h"
#include "metaprob.h"
#include "stdpso.h"
#include "optalg.h"
#include "penalty.h"
#include "sel.h"
#include "utl.h"
#include "vec.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Register problems on the problem factory:
//
static const bool registeredProb[] =
{
  cProblemFactory :: Register("MetaOptRS"    , MakeProb<cMetaOptRSD>    ,".meta"),
  cProblemFactory :: Register("MetaOptGA"    , MakeProb<cMetaOptGAD>    ,".meta"),
  cProblemFactory :: Register("MetaOptPSO"   , MakeProb<cMetaOptPSOD>   ,".meta"),
  cProblemFactory :: Register("MetaOptABC"   , MakeProb<cMetaOptABCD>   ,".meta"),
  cProblemFactory :: Register("MetaOptAIS"   , MakeProb<cMetaOptAISD>   ,".meta"),
  cProblemFactory :: Register("MetaOptLamGA" , MakeProb<cMetaOptLamGAD> ,".meta"),
  cProblemFactory :: Register("MetaOptLamPSO", MakeProb<cMetaOptLamPSOD>,".meta")
};

// -------------------------------------------------------------------------
// Class cMetaOptContinuous:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cMetaOptContinuous ===========================

cMetaOptContinuous :: cMetaOptContinuous(void)
{
}
  
// ========================== ~cMetaOptContinuous ==========================

cMetaOptContinuous :: ~cMetaOptContinuous(void)
{
  delete [] Low;
  delete [] Upp;
}

// ============================== PrintVar =================================

void cMetaOptContinuous :: PrintVar(double *x)
{ 
  for (int i = 0; i < NumVar; i++)
    cout << x[i];
  cout << endl;
} 

// ============================= WriteVar ==================================

void cMetaOptContinuous :: WriteVar(double* x, ostream &out)
{
  for (int i = 0; i < NumVar; i++) 
    out << x[i];
  out << endl;
}

// ========================== GetBoundsDouble ==============================

void cMetaOptContinuous :: GetDblBounds(double *low,double *upp)
{
  for (int i = 0; i < NumVar; i++)
  {
    low[i] = Low[i];
    upp[i] = Upp[i];
  }
}

// -------------------------------------------------------------------------
// Class cMetaOptDiscrete:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//
// =========================== cMetaOptDiscrete ============================

cMetaOptDiscrete :: cMetaOptDiscrete(void)
{
  eMetaProbType type = cMetaProblem::GetInpType( );
  ProbSet = cMetaProblem :: CreateMetaProblem(type);

  // Default values.
  OptNumber  = 100;
  EvalNumber = 0;
  Topology = GBEST_TOPOLOGY;
}

// =========================== ~cMetaOptDiscrete ===========================

cMetaOptDiscrete :: ~cMetaOptDiscrete(void)
{
  delete [] ListDim;
  delete [] List;
}
 
// ============================== LoadReadFunc =============================

void cMetaOptDiscrete :: LoadReadFunc(cInpMap &im)
{
  // Register read functions.
  im.Insert("METAOPT.ALGORITHM.FLAG", makeReadObj(cMetaOptDiscrete,ReadAlgFlag));
  im.Insert("METAOPT.EVALUATE.NUMBER", makeReadObj(cMetaOptDiscrete,ReadEvalNumber));
  im.Insert("METAOPT.OPTIMIZATION.NUMBER", makeReadObj(cMetaOptDiscrete,ReadProbOptNum));
  im.Insert("METAOPT.SWARM.TOPOLOGY", makeReadObj(cMetaOptDiscrete,ReadSwarmTopology));
  im.Insert("METAPROB.TYPE", makeRead(cMetaProblem::ReadMetaProblem));
}

// ============================== PrintVar =================================

void cMetaOptDiscrete :: PrintVar(int *algvar)
{
  for (int i = 0; i < NumVar; i++) 
    cout << List[i][algvar[i]] << "  ";
  cout << endl;
}

// ============================= WriteVar ==================================

void cMetaOptDiscrete :: WriteVar(int *algvar, ostream &out)
{
  for (int i = 0; i < NumVar; i++) 
    out << List[i][algvar[i]] << "  ";
  out << endl;
}

// ============================= GetBounds =================================

void cMetaOptDiscrete :: GetBounds(int i, int *low, int *upp)
{
  *low = 0;
  *upp = ListDim[i]-1;
}

// ============================= ReadProbOptNum ============================

void cMetaOptDiscrete :: ReadProbOptNum(istream &in)
{
  // Read number of objective function evaluations requested.

  if (!(in >> OptNumber))
  {
    cout << "Error in the input of number of objective function evaluations." << endl;
    exit(0);
  }
}

// ============================= ReadEvalNumber ============================

void cMetaOptDiscrete :: ReadAlgFlag(istream &in)
{
  // Read the number of paramater flags
  
  int numflag;

  if (!(in >> numflag))
    Utl :: Exit("Error in the input of number of parameter flags.");

  // Read each flag

  bool flag;
  
  for(int i = 0; i < numflag; i++)
    if (!(in >> flag))
      Utl :: Exit("Error in the input of number of parameter flags.");
    else
      VarFlag.push_back(flag);
}

// ============================= ReadEvalNumber ============================

void cMetaOptDiscrete :: ReadEvalNumber(istream &in)
{
  // Read number of objective function evaluations requested.

  if (!(in >> EvalNumber))
  {
    cout << "Error in the input of number of objective function evaluations." << endl;
    exit(0);
  }
}

// ============================== ReadSwarmTopology ========================

void cMetaOptDiscrete :: ReadSwarmTopology(istream &in)
{
  // Read the swarm topology used on meta-optimization.
  
  char label[100];

  if (!Utl::ReadString(in, label))
    Utl::Exit("Error in the input of the algorithm label.");

  if (string(label) == "Ring")
    Topology = RING_TOPOLOGY;
  else if (string(label) == "Gbest")
    Topology = GBEST_TOPOLOGY;
  else if (string(label) == "Square")
    Topology = SQUARE_TOPOLOGY;
  else
    Utl::Exit("Unknown Swarm Topology type: " + string(label));
}


// ============================= Decode ====================================

void cMetaOptDiscrete :: Decode(cVector &algpar, cVector &Param)
{
  int ind_var = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
  {
    if (VarFlag[i]) 
    {
      Param[i] = algpar[ind_var];
      ind_var++;
    }
    else
      Param[i] = DefaultVal[i];
  }
}

// ============================= WriteControl ==============================

void cMetaOptDiscrete :: WriteControl(fstream &in, cVector &param)
{
  in << "%HEADER\n";
  in << "Created by MetaOpt-BIOS\n";
  in << "\n%MAXIMUM.THREAD.NUMBER\n1\n";
  in << "\n%OPTIMIZATION.NUMBER\n";
  in << (int) param[0] << endl;
  in << "\n%POPULATION.SIZE\n";
  in << (int) param[1] << endl;
  in << "\n%MAXIMUM.GENERATIONS\n";
  in << (int) param[2] << endl; 
  
  // Meta-optimization files will have the same selection and penalty values
  // from the main .opt file.
  
  in << "\n%CONSTRAINT.TOLERANCE\n";
  //in << cOptAlgorithm :: GetTolViol( ) << endl;

  in << "\n%PENALTY.METHOD\n";
  /*
  switch(cPenalty :: GetInpType( ))
  {
    case STATIC: 
      in << "'Static'\n";
    break;
      
    case DEB2000:
      in << "'Deb2000'\n";
    break;
      
    case ADAPTIVE:
      in << "'Adaptive'\n";
    break;
  }
  in << "\n%SELECTION.METHOD\n";
  switch(cSelection :: GetInpType( ))
  {
    case RANKING: 
      in << "'Ranking'\n";
    break;
      
    case FITNESS_PROPORTIONAL:
      in << "'FitnessProportional'\n";
    break;

    case TOURNAMENT:
      in << "'Tournament'\n";
    break;

    case TOURNAMENTNSGAII:
      in << "'TournamentNSGAII' \n";
    break;
  }
  */
}

// ============================= Evaluate ================================

void cMetaOptDiscrete :: Evaluate(int *algvar, cVector &c, cVector &fobjs)
{
  // Decodification of problem variables.
  
  cVector x(NumVar);
  for (int i = 0; i < NumVar; i++) 
    x[i] = List[i][algvar[i]];

 // Build the OPT input file.

  stringstream thread;
#ifdef _OMP_
  thread << omp_get_thread_num( );
#else
  thread << 0;
#endif
  
  string filename = "MetaOpt" + AlgorithmSymbol( ) + thread.str( );
  string optname = filename + ".opt";
  string outname = filename + ".out";

  fstream OptFile;
  fstream OutFile;

  double *MeanObjFunc = new double [ProbSet->GetNumProb( )];
  int ProbNum = ProbSet->GetNumProb( );

  for(int prob = 0; prob < ProbNum; prob++)
  {
    OptFile.open(optname.c_str( ),fstream::out | fstream::trunc);

    WriteAlgorithm(OptFile,x); 

    ProbSet->WriteProblem(filename,OptFile,prob);

    // Write end label
    
    OptFile << "\n%END";

    OptFile.close( );

    string cmd;
  
    // Run the optimization with BIOS.
    
    cmd = "./bios MetaOpt" + AlgorithmSymbol( ) + thread.str( ) + " -silent";
    
    int status = system(cmd.c_str( ));

    if (status)
      Utl::Exit("Error in the process of bios file.");

    // Open the .out file.
    
    OutFile.open(outname.c_str( ));
    
    if (!OutFile.is_open( ))
    {
      cout << "Error opening the out file for MetaOptimization." << endl;
      exit(0);
    }
    
    // Find Mean Best Obtained.

    string label;
    bool GotTag = false;

    while (OutFile >> label)
      if (label == "%RESULT.MEAN.BEST")
      {  
        OutFile >> MeanObjFunc[prob];
	GotTag = true;
        break;
      }
    
    OutFile.close( );

    if (GotTag == false)
    {
      delete [] MeanObjFunc;
      fobjs[0] = 1e5;
      return;
    }
 
    if (ProbSet->GetObjType ( ) == MAXIMIZATION)
      MeanObjFunc[prob] /= ProbSet->GetBestObjFunc(prob);
    else if (ProbSet->GetObjType ( ) == MINIMIZATION)
      MeanObjFunc[prob] = ProbSet->GetBestObjFunc(prob)/MeanObjFunc[prob];
  }

  // Display Constaints

  for(int i = 0; i < NumConstr; i++)
    if (MeanObjFunc[i] <= 0)
      c[i] = MeanObjFunc[i];
    else 
      c[i] = -MeanObjFunc[i];

  // Evaluate Objctive Function

  double fobjsum = 0;

  for(int prob = 0; prob < ProbSet->GetNumProb( ); prob++)
    fobjsum += MeanObjFunc[prob];

  double fobj = fobjsum/ProbNum;
 
  delete [] MeanObjFunc;

  fobjs[0] = -fobj;
}

// -------------------------------------------------------------------------
// Class cMetaOptRSD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptRSD ================================

cMetaOptRSD :: cMetaOptRSD(void) : cMetaOptDiscrete( )
{
  NumVar = 1;
  NumConstr = ProbSet->GetNumProb();

  ListDim = new int[NumVar];

  List = new cVector[NumVar];

  ListDim[0] = 1;
  List[0].Resize(ListDim[0]);
  List[0][0]= 1;
}

// ============================= WriteAlgorithm ============================

void cMetaOptRSD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Evaluate pop size.

  int PopSize=EvalNumber;

  int Gen = 0;

  // Write Control pameters.

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = PopSize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'StdPSO'\n";
}

// -------------------------------------------------------------------------
// Class cMetaOptGAD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptGAD ================================

cMetaOptGAD :: cMetaOptGAD(void) : cMetaOptDiscrete( )
{
  NumVar = 3;
  NumConstr = ProbSet->GetNumProb();

  // Default variables values with flag off.

  DefaultVal.Resize(NumVar);

  DefaultVal[0] = 1.0;  // Pop/Gen ratio
  DefaultVal[1] = 0.9;  // Crossover ration
  DefaultVal[2] = 0.05; // Mutation ration

  // Check active variables.

  if (VarFlag.empty( ))
    VarFlag.assign(3,1);
  else if (VarFlag.size( ) != 3)
    Utl::Exit("Error in the number of parameter flags given. Should be 3");

  // Removing Flagged off variables.

  NumVar = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
    if(VarFlag[i]) NumVar++;
  
  // Create an array with the size of the list of each variable.
  
  vector<double> step;

  if (VarFlag[0]) step.push_back(1.0e-2);
  if (VarFlag[1]) step.push_back(5.0e-2);
  if (VarFlag[2]) step.push_back(1.0e-2);

  vector<pair<double,double> > Ranges;

  if (VarFlag[0]) Ranges.push_back(make_pair(0.25,4.0)); // Pop/Gen rate 
  if (VarFlag[1]) Ranges.push_back(make_pair(0.4,1.0));  // Crossover rate
  if (VarFlag[2]) Ranges.push_back(make_pair(0.0,0.5));  // Mutate probability
  
  ListDim = new int[NumVar];
  
  for(int i = 0; i < NumVar; i++)
    ListDim[i] = round((Ranges[i].second-Ranges[i].first)/step[i]) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  for(int i = 0; i < NumVar; i++)
      List[i].Resize(ListDim[i]);
 
  for(int i = 0; i < NumVar; i++)
    for(int j = 0; j < List[i].Dim( ); j++)
      List[i][j] = Ranges[i].first + step[i]*j;
}

// ============================= WriteAlgorithm ============================

void cMetaOptGAD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Decode active parameter

  cVector Param(VarFlag.size());
  Decode(param, Param);

  // Evaluate population size

  int PopSize = round(sqrt(EvalNumber*Param[0]));
 
  if(PopSize < 1)
    PopSize = 1;
  
  int sonsize = round(PopSize*Param[1]);

  int Gen = round((double)(EvalNumber - PopSize)/sonsize);

  if(Gen < 0)
    Gen = 0;

  // Write Control pameters

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = PopSize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'StdGA'\n";

  // Write genetic algorithm parameters

  in << "\n%CROSSOVER.RATE\n";
  in << Param[1] << endl;
  in << "\n%MUTATION.PROBABILITY\n";
  in << Param[2] << endl;
}

// -------------------------------------------------------------------------
// Class cMetaOptLamGAD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptLamGAD ================================

cMetaOptLamGAD :: cMetaOptLamGAD(void) : cMetaOptDiscrete( )
{
  NumVar = 5;
  NumConstr = ProbSet->GetNumProb();

  // Default variables values with flag off

  DefaultVal.Resize(NumVar);

  DefaultVal[0] = 1.0;  // Pop/Gen ratio
  DefaultVal[1] = 0.9;  // Crossover ration
  DefaultVal[2] = 0.05; // Mutation ration
  DefaultVal[3] = 0.00; // Layer Swap Probability
  DefaultVal[4] = 0.00; // Layer Add/Dell Probability

  // Check active variables

  if (VarFlag.empty( ))
    VarFlag.assign(NumVar,1);
  else if (VarFlag.size( ) != 5)
    Utl::Exit("Error in the number of parameter flags given. Should be 5");

  // Removing Flagged off variables

  NumVar = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
    if(VarFlag[i]) NumVar++;
  
  // Create an array with the size of the list of each variable.

  vector<double> step;

  if (VarFlag[0]) step.push_back(1.0e-2);
  if (VarFlag[1]) step.push_back(5.0e-2);
  if (VarFlag[2]) step.push_back(1.0e-2);
  if (VarFlag[3]) step.push_back(1.0e-2);
  if (VarFlag[4]) step.push_back(1.0e-2);

  vector<pair<double,double> > Ranges;

  if (VarFlag[0]) Ranges.push_back(make_pair(0.25,4.0)); // Pop/Gen rate 
  if (VarFlag[1]) Ranges.push_back(make_pair(0.4,1.0));  // Crossover rate
  if (VarFlag[2]) Ranges.push_back(make_pair(0.0,0.5));  // Mutate probability
  if (VarFlag[3]) Ranges.push_back(make_pair(0.0,0.5));  // Swap probability
  if (VarFlag[4]) Ranges.push_back(make_pair(0.0,0.5));  // Add/Del probability
  
  ListDim = new int[NumVar];
  
  for(int i = 0; i < NumVar; i++)
    ListDim[i] = round((Ranges[i].second-Ranges[i].first)/step[i]) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  for(int i = 0; i < NumVar; i++)
      List[i].Resize(ListDim[i]);
 
  for(int i = 0; i < NumVar; i++)
    for(int j = 0; j < List[i].Dim( ); j++)
      List[i][j] = Ranges[i].first + step[i]*j;
}

// ============================= WriteAlgorithm ============================

void cMetaOptLamGAD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Decode active parameter

  cVector Param(VarFlag.size());
  Decode(param, Param);

  // Evaluate population size

  int PopSize = round(sqrt(EvalNumber*Param[0]));
 
  if(PopSize < 1)
    PopSize = 1;
  
  int sonsize = round(PopSize*Param[1]);

  int Gen = round((double)(EvalNumber - PopSize)/sonsize);

  if(Gen < 0)
    Gen = 0;

  // Write Control pameters

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = PopSize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'LamGA'\n";

  // Write genetic algorithm parameters

  in << "\n%CROSSOVER.RATE\n";
  in << Param[1] << endl;
  in << "\n%MUTATION.PROBABILITY\n";
  in << Param[2] << endl;
  in << "\n%LAYERSWAP.PROBABILITY\n";
  in << Param[3] << endl;
  in << "\n%LAYERADD.PROBABILITY\n";
  in << Param[4] << endl;
  in << "\n%LAYERDEL.PROBABILITY\n";
  in << Param[4] << endl;
}

// -------------------------------------------------------------------------
// Class cMetaOptPSOD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptPSOD ===============================

cMetaOptPSOD :: cMetaOptPSOD(void) : cMetaOptDiscrete( )
{
  NumVar = 5;
  NumConstr = ProbSet->GetNumProb();

  // Default variables values with flag off

  DefaultVal.Resize(NumVar);

  DefaultVal[0] = 1.0;  // Pop/Gen ratio
  DefaultVal[1] = 1.0;  // Inertia
  DefaultVal[2] = 0.00; // c1
  DefaultVal[3] = 0.00; // c2
  DefaultVal[4] = 0.00; // Mutation

  // Check active variables

  if (VarFlag.empty( ))
    VarFlag.assign(5,1);
  else if (VarFlag.size( ) != 5)
    Utl::Exit("Error in the number of parameter flags given. Should be 5");

  // Removing Flagged off variables

  NumVar = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
    if(VarFlag[i]) NumVar++;
  
  // Create an array with the size of the list of each variable.
  
  vector<double> step;

  if (VarFlag[0]) step.push_back(1.0e-2);
  if (VarFlag[1]) step.push_back(1.0e-2);
  if (VarFlag[2]) step.push_back(1.0e-2);
  if (VarFlag[3]) step.push_back(1.0e-2);
  if (VarFlag[4]) step.push_back(1.0e-2);

  vector<pair<double,double> > Ranges;

  if (VarFlag[0]) Ranges.push_back(make_pair(0.25,4.0)); // Pop/Gen rate 
  if (VarFlag[1]) Ranges.push_back(make_pair(0.0,2.0));  // Inertia
  if (VarFlag[2]) Ranges.push_back(make_pair(0.0,4.0));  // c1
  if (VarFlag[3]) Ranges.push_back(make_pair(0.0,4.0));  // c2
  if (VarFlag[4]) Ranges.push_back(make_pair(0.0,0.5));  // Mutation
  
  ListDim = new int[NumVar];
  
  for(int i = 0; i < NumVar; i++)
    ListDim[i] = round((Ranges[i].second-Ranges[i].first)/step[i]) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  for(int i = 0; i < NumVar; i++)
      List[i].Resize(ListDim[i]);
 
  for(int i = 0; i < NumVar; i++)
    for(int j = 0; j < List[i].Dim( ); j++)
      List[i][j] = Ranges[i].first + step[i]*j;
}

// ============================= WriteAlgorithm ============================

void cMetaOptPSOD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Decode active parameter

  cVector Param(VarFlag.size());
  Decode(param, Param);

  // Evaluate swarm size

  int SwarmSize = round(sqrt(Param[0]*EvalNumber));
  
  if(SwarmSize < 1)
    SwarmSize = 1;

  int Gen = round((double) (EvalNumber)/SwarmSize -1);

  if(Gen < 0)
    Gen = 0;

  // Write Control pameters

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = SwarmSize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'StdPSO'\n";
  
  // Write Particle Swarm optimizer parameters

  in << "\n%SWARM.TOPOLOGY\n";
  in << Topology << endl;

  in << "\n%PARTICLE.INERTIA\n";
  in << "'CONSTANT' " << Param[1] << endl;
  in << "\n%COGNITIVE.FACTOR\n";
  in << "'CONSTANT' " << Param[2] << endl;
  in << "\n%SOCIAL.FACTOR\n";
  in << "'CONSTANT' " << Param[3] << endl;
  in << "\n%MUTATION.PROBABILITY\n";
  in << Param[4] << endl;
}

// -------------------------------------------------------------------------
// Class cMetaOptLamPSOD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptLamPSOD ===============================

cMetaOptLamPSOD :: cMetaOptLamPSOD(void) : cMetaOptDiscrete( )
{
  eMetaProbType type = cMetaProblem::GetInpType( );
  ProbSet = cMetaProblem :: CreateMetaProblem(type);

  NumVar = 10;
  NumConstr = ProbSet->GetNumProb();

  // Default variables values with flag off

  DefaultVal.Resize(NumVar);

  DefaultVal[0] = 1.0;  // Pop/Gen ratio
  DefaultVal[1] = 0.0;  // Inertia
  DefaultVal[2] = 0.00; // c1
  DefaultVal[3] = 0.00; // c2
  DefaultVal[4] = 0.00; // Mutation probability
  DefaultVal[5] = 0.00; // LamMutation Thk probability
  DefaultVal[6] = 0.00; // LamMutation Ang probability
  DefaultVal[7] = 0.00; // LamMutation Mat probability
  DefaultVal[8] = 0.00; // Layer swap probability
  DefaultVal[9] = 0.00; // Layer Add/Del probability

  // Check active variables

  if (VarFlag.empty( ))
    VarFlag.assign(10,1);
  else if (VarFlag.size( ) != 10)
    Utl::Exit("Error in the number of parameter flags given. Should be 10!");

  // Removing Flagged off variables

  NumVar = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
    if(VarFlag[i]) NumVar++;
  
  // Create an array with the size of the list of each variable.
  
  vector<double> step;

  if (VarFlag[0]) step.push_back(1.0e-2);
  if (VarFlag[1]) step.push_back(1.0e-2);
  if (VarFlag[2]) step.push_back(1.0e-2);
  if (VarFlag[3]) step.push_back(1.0e-2);
  if (VarFlag[4]) step.push_back(1.0e-2);
  if (VarFlag[5]) step.push_back(1.0e-2);
  if (VarFlag[6]) step.push_back(1.0e-2);
  if (VarFlag[7]) step.push_back(1.0e-2);
  if (VarFlag[8]) step.push_back(1.0e-2);
  if (VarFlag[9]) step.push_back(1.0e-2);

  vector<pair<double,double> > Ranges;

  if (VarFlag[0]) Ranges.push_back(make_pair(0.25,4.0)); // Pop/Gen rate 
  if (VarFlag[1]) Ranges.push_back(make_pair(0.0,2.0));  // Inertia
  if (VarFlag[2]) Ranges.push_back(make_pair(0.0,4.0));  // c1
  if (VarFlag[3]) Ranges.push_back(make_pair(0.0,4.0));  // c2
  if (VarFlag[4]) Ranges.push_back(make_pair(0.0,0.5));  // Mutation
  if (VarFlag[5]) Ranges.push_back(make_pair(0.0,0.5));  // LamMutation Thk
  if (VarFlag[6]) Ranges.push_back(make_pair(0.0,0.5));  // LamMutation Ang
  if (VarFlag[7]) Ranges.push_back(make_pair(0.0,0.5));  // LamMutation Mat
  if (VarFlag[8]) Ranges.push_back(make_pair(0.0,0.5));  // Swap
  if (VarFlag[9]) Ranges.push_back(make_pair(0.0,0.5));  // Add / Del
  
  ListDim = new int[NumVar];
  
  for(int i = 0; i < NumVar; i++)
    ListDim[i] = round((Ranges[i].second-Ranges[i].first)/step[i]) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  for(int i = 0; i < NumVar; i++)
      List[i].Resize(ListDim[i]);
 
  for(int i = 0; i < NumVar; i++)
    for(int j = 0; j < List[i].Dim( ); j++)
      List[i][j] = Ranges[i].first + step[i]*j;
}

// ============================= WriteAlgorithm ============================

void cMetaOptLamPSOD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Decode active parameter

  cVector Param(VarFlag.size());
  Decode(param, Param);

  // Evaluate swarm size

  int SwarmSize = round(sqrt(Param[0]*EvalNumber));
  
  if(SwarmSize < 1)
    SwarmSize = 1;

  int Gen = round((double) (EvalNumber)/SwarmSize -1);

  if(Gen < 0)
    Gen = 0;

  // Write Control pameters

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = SwarmSize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'LamPSO'\n";
  
  // Write Particle Swarm optimizer parameters

  in << "\n%SWARM.TOPOLOGY\n";
  in << Topology << endl;

  in << "\n%PARTICLE.INERTIA\n";
  in << "'CONSTANT' " << Param[1] << endl;
  in << "\n%COGNITIVE.FACTOR\n";
  in << "'CONSTANT' " << Param[2] << endl;
  in << "\n%SOCIAL.FACTOR\n";
  in << "'CONSTANT' " << Param[3] << endl;
  in << "\n%MUTATION.PROBABILITY\n";
  in << Param[4] << endl;
  in << "\n%LAMINATE.MUTATION.PROBABILITY\n";
  in << Param[5] << " " << Param[6] << " " << Param[7] << endl;
  in << "\n%LAYERSWAP.PROBABILITY\n";
  in << Param[8] << endl;
  in << "\n%LAYERADD.PROBABILITY\n";
  in << Param[9] << endl;
  in << "\n%LAYERDEL.PROBABILITY\n";
  in << Param[9] << endl;
}

// -------------------------------------------------------------------------
// Class cMetaOptABCD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptABCD ===============================

cMetaOptABCD :: cMetaOptABCD(void) : cMetaOptDiscrete( )
{
  eMetaProbType type = cMetaProblem::GetInpType( );
  ProbSet = cMetaProblem :: CreateMetaProblem(type);

  NumVar = 4;
  NumConstr = ProbSet->GetNumProb();

  // Default variables values with flag off

  DefaultVal.Resize(NumVar);

  DefaultVal[0] =  1.0;  // Pop/Gen ratio
  DefaultVal[1] =  0.5;  // Food ratio
  DefaultVal[2] =   10;  // Max Trial Number
  DefaultVal[3] = 0.00; // ScoutBees Ratio

  // Check active variables

  if (VarFlag.empty( ))
    VarFlag.assign(4,1);
  else if (VarFlag.size( ) != 4)
    Utl::Exit("Error in the number of parameter flags given. Should be 4");

  // Removing Flagged off variables

  NumVar = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
    if(VarFlag[i]) NumVar++;
  
  // Create an array with the size of the list of each variable.
  
  vector<double> step;

  if (VarFlag[0]) step.push_back(1.0e-2);
  if (VarFlag[1]) step.push_back(1.0e-2);
  if (VarFlag[2]) step.push_back(1.0e-2);
  if (VarFlag[3]) step.push_back(1.0e-2);

  vector<pair<double,double> > Ranges;

  if (VarFlag[0]) Ranges.push_back(make_pair(0.25,4.0)); // Pop/Gen rate 
  if (VarFlag[1]) Ranges.push_back(make_pair(0.0,1.0));  // Food rate
  if (VarFlag[2]) Ranges.push_back(make_pair(0.0,4.0));  // Max Trial rate
  if (VarFlag[3]) Ranges.push_back(make_pair(0.0,0.5));  // ScoutBees rate
  
  ListDim = new int[NumVar];
  
  for(int i = 0; i < NumVar; i++)
    ListDim[i] = round((Ranges[i].second-Ranges[i].first)/step[i]) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  for(int i = 0; i < NumVar; i++)
      List[i].Resize(ListDim[i]);
 
  for(int i = 0; i < NumVar; i++)
    for(int j = 0; j < List[i].Dim( ); j++)
      List[i][j] = Ranges[i].first + step[i]*j;
}

// ============================= WriteAlgorithm ============================

void cMetaOptABCD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Decode active parameter

  cVector Param(VarFlag.size());
  Decode(param, Param);

  // Evaluate swarm size

  int SwarmSize = round(sqrt(Param[0]*EvalNumber));
  
  if(SwarmSize < 1)
    SwarmSize = 1;

  int Gen = round((double) (EvalNumber)/SwarmSize -1);

  if(Gen < 0)
    Gen = 0;

  // Write Control pameters

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = SwarmSize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'StdABC'\n";
  
  // Write Artificial Bee Colony optimizer parameters

  in << "\n%FOOD.RATE\n";
  in << Param[1] << endl;
  in << "\n%FOOD.MAXTRIAL.RATE\n";
  in << Param[2] << endl;
  in << "\n%SCOUTBEES.RATE\n";
  in << Param[3] << endl;
}

// -------------------------------------------------------------------------
// Class cMetaOptAISD:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cMetaOptAISD ===============================

cMetaOptAISD :: cMetaOptAISD(void) : cMetaOptDiscrete( )
{
  eMetaProbType type = cMetaProblem::GetInpType( );
  ProbSet = cMetaProblem :: CreateMetaProblem(type);

  NumVar = 4;
  NumConstr = ProbSet->GetNumProb();

  // Default variables values with flag off

  DefaultVal.Resize(NumVar);

  DefaultVal[0] =  1.0;  // Pop/Gen ratio
  DefaultVal[1] = 1.00;  // Clone rate
  DefaultVal[2] = 5.00;  // Hypermutate decay rate
  DefaultVal[3] = 0.00;  // Receptor editing rate

  // Check active variables

  if (VarFlag.empty( ))
    VarFlag.assign(4,1);
  else if (VarFlag.size( ) != 4)
    Utl::Exit("Error in the number of parameter flags given. Should be 4");

  // Removing Flagged off variables

  NumVar = 0;

  for(unsigned int i = 0; i < VarFlag.size( ); i++)
    if(VarFlag[i]) NumVar++;
  
  // Create an array with the size of the list of each variable.
  
  vector<double> step;

  if (VarFlag[0]) step.push_back(1.0e-2);
  if (VarFlag[1]) step.push_back(5.0e-2);
  if (VarFlag[2]) step.push_back(1.0);
  if (VarFlag[3]) step.push_back(1.0e-2);

  vector<pair<double,double> > Ranges;

  if (VarFlag[0]) Ranges.push_back(make_pair(0.25,4.0)); // Pop/Gen rate 
  if (VarFlag[1]) Ranges.push_back(make_pair(0.2,1.0));  // Clone rate
  if (VarFlag[2]) Ranges.push_back(make_pair(1.0,20.0));  // Decay rate
  if (VarFlag[3]) Ranges.push_back(make_pair(0.0,0.5));  // Receptor Edt rate
  
  ListDim = new int[NumVar];
  
  for(int i = 0; i < NumVar; i++)
    ListDim[i] = round((Ranges[i].second-Ranges[i].first)/step[i]) + 1;

  // Create the list of discrete values for each variable.

  List = new cVector[NumVar];

  for(int i = 0; i < NumVar; i++)
      List[i].Resize(ListDim[i]);
 
  for(int i = 0; i < NumVar; i++)
    for(int j = 0; j < List[i].Dim( ); j++)
      List[i][j] = Ranges[i].first + step[i]*j;
}

// ============================= WriteAlgorithm ============================

void cMetaOptAISD :: WriteAlgorithm(fstream &in, cVector &param)
{
  // Decode active parameter

  cVector Param(VarFlag.size());
  Decode(param, Param);

  // Evaluate memory cells size

  int memsize = round(sqrt(Param[0]*EvalNumber));
  
  if(memsize < 1)
    memsize = 1;

  int clonesize = round(memsize*Param[1]);
  int Rptsize = round(memsize*Param[3]);
 
  int Gen;

  if (clonesize+Rptsize > 0)
    Gen = round((double)(EvalNumber-memsize)/(clonesize + Rptsize));
  else
    Gen = 0;

  if(Gen < 0)
    Gen = 0;
  
  // Write Control pameters

  cVector ControlParam(3);

  ControlParam[0] = OptNumber; 
  ControlParam[1] = memsize; 
  ControlParam[2] = Gen;

  WriteControl(in,ControlParam);

  in << "\n%OPTIMIZATION.ALGORITHM\n";
  in << "'StdAIS'\n";
  
  // Write Artificial Immune System optimizer parameters

  in << "\n%CLONE.RATE\n";
  in << Param[1] << endl;
  in << "\n%HYPERMUTATE.DECAY.RATE\n";
  in << Param[2] << endl;
  in << "\n%RECEPTOR.EDITING.RATE\n";
  in << Param[3] << endl;
}


/*
// -------------------------------------------------------------------------
// Class cGSuite1C:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cGSuite1C ===============================

cGSuite1C :: cGSuite1C(void)
{
  Type = GSUITE1;
  NumVar = 13;
  NumConstr = 9;

  Low = new double[NumVar];
  Upp = new double[NumVar];
  Low[0] = 0;
  Upp[0] = 1;
  Low[1] = 0;
  Upp[1] = 1;
  Low[2] = 0;
  Upp[2] = 1;
  Low[3] = 0;
  Upp[3] = 1;
  Low[4] = 0;
  Upp[4] = 1;
  Low[5] = 0;
  Upp[5] = 1;
  Low[6] = 0;
  Upp[6] = 1;
  Low[7] = 0;
  Upp[7] = 1;
  Low[8] = 0;
  Upp[8] = 1;
  Low[9] = 0;
  Upp[9] = 100;
  Low[10] = 0;
  Upp[10] = 100;
  Low[11] = 0;
  Upp[11] = 100;
  Low[12] = 0;
  Upp[12] = 1;
 }

// ============================= Evaluate ================================

double cGSuite1C :: Evaluate(cVector &x, cVector &c)
{
  // Constraint evaluation.
  
  c[0] =  2.0*x[0] + 2.0*x[1] + x[9]  + x[10] - 10.0;
  c[1] =  2.0*x[0] + 2.0*x[2] + x[9]  + x[11] - 10.0;
  c[2] =  2.0*x[1] + 2.0*x[2] + x[10] + x[11] - 10.0;
  c[3] = -8.0*x[0] +     x[9];
  c[4] = -8.0*x[1] +     x[10];
  c[5] = -8.0*x[2] +     x[11];
  c[6] = -2.0*x[3] - 1.0*x[4] + x[9];
  c[7] = -2.0*x[5] - 1.0*x[6] + x[10];
  c[8] = -2.0*x[7] - 1.0*x[8] + x[11];
  
  // Objective function evaluation.
    
  double sum1 = 0.0;
  for (int i = 0; i < 4; i++) sum1 += x[i];
   
  double sum2 = 0.0;
  for (int i = 0; i < 4; i++) sum2 += x[i]*x[i];
    
  double sum3 = 0.0;
  for (int i = 4; i < 13; i++) sum3 += x[i];

  double fobj = 5.0*sum1 - 5.0*sum2 - sum3; 
  return(fobj);
}
*/  

// ======================================================= End of file =====
