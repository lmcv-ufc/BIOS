// -------------------------------------------------------------------------
// stdde.cpp - implementation of cStandardDE class.
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
// Created:      19-Jun-2020    Leonardo Gonçalves Ribeiro
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "stdde.h"
#include "problem.h"
#include "sel.h"
#include "group.h"
#include "individual.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"
#include "optalg.h"

// -------------------------------------------------------------------------
// Class StandardDE:
//

// -------------------------------------------------------------------------
// Public methods:
//
// ============================= ReadCrossRate =============================

void cStandardDE :: ReadCrossRate(istream &in)
{
  if (!(in >> CrossRate))
  {
    cout << "Error in the input of the crossover rate." << endl;
    exit(0);
  }
}

// ============================ ReadCrossRange =============================

void cStandardDE :: ReadCrossRange(istream &in)
{
  if (!(in >> MaxCross) || !(in >> MinCross))
  {
    cout << "Error in the input of the crossover range." << endl;
    exit(0);
  }
}

// ============================ ReadFWeight =============================

void cStandardDE :: ReadFWeight(istream &in)
{
  if (!(in >> CrossRate))
  {
    cout << "Error in the input of the weight." << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cStandardDE :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("CROSSOVER.RANGE",      makeReadObj(cStandardDE,ReadCrossRange));
  im.Insert("CROSSOVER.RATE",       makeReadObj(cStandardDE,ReadCrossRate));
  im.Insert("F.WEIGHT",             makeReadObj(cStandardDE,ReadFWeight));
}

// ============================== cStandardDE ==============================

cStandardDE :: cStandardDE(void) : cOptAlgorithm( )
{
  Type = STANDARD_DE;
  CrossRate = 0.80;
  MutProb = 0.0;
  FWeight = 0.85;
  DifType = Loc2Best;
  CrossType = BIN;
}

// ============================= ~cStandardDE ==============================

cStandardDE :: ~cStandardDE(void)
{
}

// =============================== Solver ==================================

void cStandardDE :: Solver(void)
{
   // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Check initial conditions

      if (PopSize < 4)
      {
          cout << "Population size should not be lower than 4" << endl;
          cout << "Population size was set to 4" << endl;
          PopSize = 4;
      }
      if (CrossRate < 0 || CrossRate > 1){
          cout << "Crossover rate should be in the range [0, 1]" << endl;
          cout << "Crossover rate was set to 0.8" << endl;
          CrossRate = 0.8;
      }
      if (FWeight < 0 || FWeight > 1){
          cout << "F should be in the range [0, 1]" << endl;
          cout << "F was set to 0.8" << endl;
          FWeight = 0.85;
      }

    // Track number of individual evaluations.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Randomize rates.

    RandomRates( );

    // Create the population, mating pool and parent array.

    cPopulation pop(PopSize, SolType, Prob);
    cPopulation popold(PopSize, SolType, Prob);

    // Evaluate initial sample points.

    vector<cVector> sx;
    int sizsamp = PopSize - NumInpSol;
    if (sizsamp < 1) IntPopSamp = 0;

    if (IntPopSamp)
    {
      sx.reserve(sizsamp);
      int nvar = Prob->VarNumRow( ) * Prob->VarNumCol( );

      cSamp InitialSample;
      InitialSample.InitSample(SampType,nvar,sizsamp,sx);
    }

    // Generate initial population.

    for (int i = 0; i < PopSize; i++)
    {
      // Initialize particles position and velocity.
      if (i < NumInpSol)
        pop[i]->Init(InpSolVec[i]); // Input values.
      else if (IntPopSamp)
        pop[i]->Init(sx[i]);        // Random values.
      else
        pop[i]->Init( );            // Random values.
    }

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < PopSize; i++)
    {
      pop[i]->Evaluate( );
      #pragma omp critical
      EvalNum++;
    }

    if (Pen)
      Pen->EvalPenObjFunc(&pop, TolViol);

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    // Perform the DE iterations.

    for (int gen = 0; gen < MaxGen; gen++)
    {
      if ((gen+1)%1 == 0 && Feedback) cout << "Generation: " << gen + 1 << endl;

      for (int i = 0; i < PopSize; i++)
      {
        popold[i]->Copy(pop[i]);
      }

      // Create the mating pool.
      cVector IndexPar(PopSize);

      Differentiation(pop, popold, IndexPar);

      // Crossover.

      Crossover(pop, popold, IndexPar);

      // Evaluate the objective function and constraints of offsprings.

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < PopSize; i++)
      {
        pop[i]->Evaluate( );

	#pragma omp critical
	EvalNum++;
      }

     // Evaluate the fitness function of the offspring population.

      if (Pen) Pen->EvalPenObjFunc(&pop, TolViol);

      Fitness(pop);

      for (int i = 0; i < PopSize; i++)
      {
          double newobj = 0;
          double oldobj = 0;
          if (Pen)
          {
              newobj = pop[i]->GetPenObjFunc( );
              oldobj = popold[i]->GetPenObjFunc( );
          }
          else
          {
              newobj = pop[i]->GetObjFunc(0);
              oldobj = popold[i]->GetObjFunc(0);
          }
          if (oldobj < newobj)
          {
              pop[i] -> Copy(popold[i]);
          }
      }

      // Migration.

      #ifdef _MPI_
      if (!((gen+1)%MigrationGap) && (gen+1) != MaxGen)
      {
        pop.Sort( );
        Migration(pop);
        if (Pen)
          Pen->EvalPenObjFunc(&pop, TolViol);
      }
      #endif

      // Update variables related to PostProcessing.

      UpdatePostVar(gen, opt, lastBest, &pop);
    
      // Check conditions to stop optimization

      if (OptStopCrit(gen, opt, lastBest, &pop))
        break;   
    }


    // Store the best individual.
    best->Insert(pop.BestSol( ));

    // Print data in the output file.

    PrintPostVar(MaxGen, opt, EvalNum, &pop);
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= Fitness ===================================

void cStandardDE :: Fitness(cPopulation &pop)
{
  // Compute the fitness function from the objective function.

  if (!Pen) // Unconstrained problems.
  {
    double bigval = abs(pop[0]->GetObjFunc(0));
    for (int i = 1; i < pop.GetSize( ); i++)
      bigval = MAX(bigval, pop[i]->GetObjFunc(0));

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < pop.GetSize( ); i++)  pop[i]->AssignFitFunc(1.10*bigval - pop[i]->GetObjFunc(0));
  }
  else      // Constrained problems.
  {
    double bigval = abs(pop[0]->GetPenObjFunc( ));
    for (int i = 1; i < pop.GetSize( ); i++)
      bigval = MAX(bigval, pop[i]->GetPenObjFunc( ));

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < pop.GetSize( ); i++)
      pop[i]->AssignFitFunc(1.10*bigval - pop[i]->GetPenObjFunc( ));
  }
}

// ============================ Differentiation ==================================

void cStandardDE :: Differentiation(cPopulation &pop, cPopulation &popold, cVector &IndexPar)
{
  // Auxiliary populations index
  vector<int> IndexPop1(PopSize);
  vector<int> IndexPop2(PopSize);
  vector<int> IndexPop3(PopSize);

  // Find the best individuals in a generation

  double AuxBestObj = 0;
  int BestIndex = 0;

  for (int i = 0; i < PopSize; i++)
  {
      if (i == 0){
          if (Pen)
          {
            AuxBestObj = popold.GetSol(i)->GetCurrPenObjFunc( );
          }
          else
          {
            AuxBestObj = popold.GetSol(i)->GetCurrObjFunc(0);
          }
          BestIndex = i;
      }

      if (Pen)
      {
        if (AuxBestObj >= popold.GetSol(i)->GetCurrPenObjFunc( ))
        {
            AuxBestObj = popold.GetSol(i)->GetCurrPenObjFunc( );
            BestIndex = i;
        }
      }
      else
      {
        if (AuxBestObj >= popold.GetSol(i)->GetCurrObjFunc(0))
        {
            AuxBestObj = popold.GetSol(i)->GetCurrObjFunc(0);
            BestIndex = i;
        }
      }
  }

  // Create the auxiliary populations (random index)

  for (int i = 0; i < PopSize; i++)
  {
    IndexPop1[i] = i;
    IndexPop2[i] = i;
    IndexPop3[i] = i;
  }

  random_shuffle(IndexPop1.begin(), IndexPop1.end());
  random_shuffle(IndexPop2.begin(), IndexPop2.end());
  random_shuffle(IndexPop3.begin(), IndexPop3.end());

  // Perform the differentiation operation.

  if(DifType == Rand1)
  {
      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < pop.GetSize( ); i = i+1)
      {
        pop[ i ]->Differentiation(DifType, FWeight, popold[BestIndex], popold[IndexPop1[i]], popold[IndexPop2[i]], popold[IndexPop3[i]]);
        IndexPar[i] = IndexPop3[i];
      }
  }
  else if(DifType == Loc2Best)
  {
      for (int i = 0; i < PopSize; i++)
      {
        IndexPop3[i] = i;
      }
      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < pop.GetSize( ); i = i+1)
      {
        pop[ i ]->Differentiation(DifType, FWeight, popold[BestIndex], popold[IndexPop1[i]], popold[IndexPop2[i]], popold[IndexPop3[i]]);
        IndexPar[i] = IndexPop3[i];
      }
  }
  else if(DifType == BestJitter)
  {
      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < pop.GetSize( ); i = i+1)
      {
        pop[ i ]->Differentiation(DifType, FWeight, popold[BestIndex], popold[IndexPop1[i]], popold[IndexPop2[i]], popold[IndexPop3[i]]);
        for (int i = 0; i < PopSize; i++)
        {
          IndexPop3[i] = BestIndex;
        }
        IndexPar[i] = IndexPop3[i];
      }
  }
}

// ============================ Crossover ==================================

void cStandardDE :: Crossover(cPopulation &son, cPopulation &parent, cVector IndexPar)
{
  // Perform the crossover operation.
  cPopulation auxpop(PopSize, SolType, Prob);

  for (int i = 0; i < PopSize; i++)
  {
    auxpop[i]->Copy(son[i]);
  }
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < son.GetSize( ); i = i+1)
  {
    //cout << "Individuo " << i << endl;
    son[ i ]->Crossover(CrossType, CrossRate, auxpop[i],parent[IndexPar[i]]);
    //cout << "\n" << endl;
  }
}

// ============================= Migration =================================

void cStandardDE :: Migration(cPopulation &pop)
{
#ifdef _MPI_

  // Get deme size.

  int popsize = pop.GetSize( );

  // Get size and rank.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );

  // Check population size constraint.

  if ((size-1)*MigrationNum >= pop.GetSize( ))
  {
    cout << "Number of received individuals in one deme is greater than the    original population size." << endl;
    exit(0);
  }

  // Send and receive individuals.

  for (int i = 0; i < MigrationNum; i++)
    for (int deme = 0; deme < size; deme++)
    {
      if (rank == deme)
        pop[i]->Send( );
      else
        pop[popsize-1-i]->Receive(deme);
    }

#endif
}

// ================================ Merge ==================================

void cStandardDE :: Merge(cPopulation &pop, cPopulation &son)
{
  // Get necessary data.

  int popsize = pop.GetSize( );
  int sonsize = son.GetSize( );
  int last = popsize - 1;
 // cout << "popsize = " << popsize << "  ";
 // cout << "sonsize = " << sonsize << "\n";

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < sonsize; i++)
  {
    pop[last-i]->Copy(son[i]);
  }
}

// ============================= RandomRates ==============================

void cStandardDE :: RandomRates(void)
{
  if (!MutProb && MaxMut)
    MutProb = Utl::RandDouble(MinMut, MaxMut);

  if (!CrossRate && MaxCross)
    CrossRate = Utl::RandDouble(MinCross, MaxCross);
}

// ======================================================= End of file =====
