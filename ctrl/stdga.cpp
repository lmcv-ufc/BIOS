// -------------------------------------------------------------------------
// stdga.cpp - implementation of cStandardGA class.
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
// Created:      22-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     15-Mar-2013    Evandro Parente Junior
//               Renamed from convga.cpp to stdga.cpp. The algorithm can
//               be used with different types of individuals.
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
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

#include "stdga.h"
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
// Class StandardGA:
//

// -------------------------------------------------------------------------
// Public methods:
//
// ============================= ReadCrossRate =============================

void cStandardGA :: ReadCrossRate(istream &in)
{
  if (!(in >> CrossRate))
  {
    cout << "Error in the input of the crossover rate." << endl;
    exit(0);
  }
}

// ============================ ReadCrossRange =============================

void cStandardGA :: ReadCrossRange(istream &in)
{
  if (!(in >> MaxCross) || !(in >> MinCross))
  {
    cout << "Error in the input of the crossover range." << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cStandardGA :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("CROSSOVER.RANGE", makeReadObj(cStandardGA,ReadCrossRange));
  im.Insert("CROSSOVER.RATE",  makeReadObj(cStandardGA,ReadCrossRate));
}

// ============================== cStandardGA ==============================

cStandardGA :: cStandardGA(void) : cOptAlgorithm( )
{
  Type = STANDARD_GA;
  CrossRate = 0.80;
  MutProb = 0.10;
}

// ============================= ~cStandardGA ==============================

cStandardGA :: ~cStandardGA(void)
{
}

// =============================== Solver ==================================

void cStandardGA :: Solver(void)
{
  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Track number of individual evaluations.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Randomize rates.

    RandomRates( );

    // Create the population, mating pool and parent array.

    int parsize = round(CrossRate*PopSize);
    int sonsize = parsize - parsize%2;
    cPopulation pop(PopSize, SolType, Prob);
    cPopulation son(sonsize, SolType, Prob);
    cPopulation matpool(PopSize, SolType, Prob);
    cPopulation parent(parsize, SolType, Prob);

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

    // Generate the initial population.

    for (int i = 0; i < PopSize; i++)
    {
      // Initialize genotypes.
      if (i < NumInpSol)
        pop[i]->Init(InpSolVec[i]); // Input values.
      else if (IntPopSamp)
        pop[i]->Init(sx[i]);        // Sample values.
      else
      {
        pop[i]->Init( );            // Random values.
      }
      EvalNum++;
    }

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < PopSize; i++)
    {
     // Evaluate the objective function and constraints of the initial pop.
      pop[i]->Evaluate( );
      #pragma omp critical
      EvalNum++;
    }

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    // Perform the GA iterations.

    for (int gen = 0; gen < MaxGen; gen++)
    {
      if ((gen+1)%1 == 0 && Feedback) cout << "Generation: " << gen + 1 << endl;
      // Evaluate the penalized objective function of the population.

      if (Pen)
        Pen->EvalPenObjFunc(&pop, TolViol);

      // Evaluate fitness function.

      Fitness(pop);

      // Create the mating pool.

      Sel->Select(1.0, &pop, &matpool);

      // Select parents for crossover.

      Sel->Select(CrossRate, &matpool, &parent);

      // Crossover.

      Crossover(parent, son);

      // Mutation.

      Mutation(son);

      // Evaluate the objective function and constraints of offsprings.

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < sonsize; i++)
      {
        son[i]->Evaluate( );

	#pragma omp critical
	EvalNum++;
      }

     // Evaluate the fitness function of the offspring population.

      if (Pen) Pen->EvalPenObjFunc(&son, TolViol);

      Fitness(son);

      // Merge the offspring population with the original one.

      pop.Sort( );
      Merge(pop, son);

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

void cStandardGA :: Fitness(cPopulation &pop)
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

// ============================ Crossover ==================================

void cStandardGA :: Crossover(cPopulation &parent, cPopulation &son)
{
  // Perform the crossover operation.

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < son.GetSize( ); i = i+2)
  {
    double r = Utl::RandDec( );
    son[ i ]->Crossover(CrossType, r, parent[i],parent[i+1]);
    son[i+1]->Crossover(CrossType, r, parent[i+1],parent[i]);
  }
}

// ============================== Mutation =================================

void cStandardGA :: Mutation(cPopulation &son)
{
  // Perform the mutation operation.

  if (MutProb)
  {
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < son.GetSize( ); i++)
      son[i]->Mutate(MutProb);
  }
}

// ============================= Migration =================================

void cStandardGA :: Migration(cPopulation &pop)
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

void cStandardGA :: Merge(cPopulation &pop, cPopulation &son)
{
  // Get necessary data.

  int popsize = pop.GetSize( );
  int sonsize = son.GetSize( );
  int last = popsize - 1;

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < sonsize; i++)
  {
    pop[last-i]->Copy(son[i]);
  }
}

// ============================= RandomRates ==============================

void cStandardGA :: RandomRates(void)
{
  if (!MutProb && MaxMut)
    MutProb = Utl::RandDouble(MinMut, MaxMut);

  if (!CrossRate && MaxCross)
    CrossRate = Utl::RandDouble(MinCross, MaxCross);
}

// ======================================================= End of file =====
