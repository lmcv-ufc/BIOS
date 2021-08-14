// -------------------------------------------------------------------------
// stdais.cpp - implementation of cStandardAIS class.
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
// Created:      29-Jun-2014    Elias Saraiva Barroso
//
// Modified:     
//  
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "stdais.h"
#include "problem.h"
#include "sel.h"
#include "group.h"
#include "optsolution.h"
#include "individual.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//
double cStandardAIS :: ReceptEditingRate = 0.05;
double cStandardAIS :: CloneRate         = 1.00;
double cStandardAIS :: Decay             = 5.00;

// -------------------------------------------------------------------------
// Set read functions labels:
//

static bool ReadFuncRegister[] =
{
  CtrlMap( ).Insert("CLONE.RATE", cStandardAIS :: ReadCloneRate),
  CtrlMap( ).Insert("HYPERMUTATE.DECAY.RATE", cStandardAIS :: ReadDecayRate),
  CtrlMap( ).Insert("RECEPTOR.EDITING.RATE", cStandardAIS :: ReadRptEdtRate)
};

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadCloneRate =============================

void cStandardAIS :: ReadCloneRate(void)
{
  if (!(in >> CloneRate))
  {
    cout << "Error in the input of the clone rate." << endl;
    exit(0);
  }
} 

// ============================= ReadRptEdtRate ============================

void cStandardAIS :: ReadRptEdtRate(void)
{
  if (!(in >> ReceptEditingRate))
  {
    cout << "Error in the input of the receptor editing rate." << endl;
    exit(0);
  }
}

// ============================= ReadDecayRate =============================

void cStandardAIS :: ReadDecayRate(void)
{
  if (!(in >> Decay))
  {
    cout << "Error in the input of the decay rate." << endl;
    exit(0);
  }
}

// ============================== cStandardAIS ==============================

cStandardAIS :: cStandardAIS(void)
{
  Type = STANDARD_AIS;
  Prob = 0;
  Pen  = 0;
  Sel  = 0;
}

// ============================= ~cStandardAIS ==============================

cStandardAIS :: ~cStandardAIS(void)
{
}

// =============================== Solver ==================================

void cStandardAIS :: Solver(void)
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
    
    // Create the memory and clone cells population.
    
    int CloneSize = round(CloneRate*PopSize);
    int ReceptEditingSize = round(PopSize * ReceptEditingRate);
    int MemoryRenewSize = PopSize - ReceptEditingSize;

    if (CloneSize <= 0)
      CloneSize = 1;

    if (MemoryRenewSize <= 0)
      MemoryRenewSize = 1;

    cPopulation MemoryCells(PopSize, SolType, Prob);
    cPopulation CloneCells(CloneSize, SolType, Prob);

    // Create pointer array to track all cells

    cOptSolution **AllCells = new cOptSolution* [CloneSize+PopSize];

    for(int i = 0; i < PopSize; i++)
      AllCells[i] = MemoryCells[i];

    for(int i = PopSize; i < (PopSize+CloneSize); i++)
      AllCells[i] = CloneCells[i-PopSize];

    // Generate the initial population.

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < PopSize; i++)
    {
      // Fill cells with random values.
      
      MemoryCells[i]->Init( );
      
      // Evaluate the objective function and constraints of the initial pop.
      
      MemoryCells[i]->Evaluate( );

      #pragma omp critical
      EvalNum++;
    }

    if (Pen)
      Pen->EvalPenObjFunc(&MemoryCells, TolViol);

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;
   
    if (Feedback) cout << "Generation: ";

    for (int gen = 0; gen < MaxGen; gen++)
    {
      if ((gen+1)%5 == 0 && Feedback) cout << "Generation: " << gen + 1 << endl;

      // Evaluate the penalized objective function of the population.
      
      if (Pen)
        Pen->EvalPenObjFunc(&MemoryCells, TolViol);
      
      // Evaluate fitness function.
      
      Fitness(MemoryCells);

      // Create the clone population.
      
      Sel->Select(CloneRate, &MemoryCells, &CloneCells);

      // Hypermutation.
      
      GlobalFitness(AllCells, (PopSize + CloneSize));
      Hypermutation(CloneCells, MemoryCells.BestSol( ), MemoryCells.WorstSol( ));
      
      // Evaluate the objective function and constraints of clones.
      
      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < CloneSize; i++)
      {
        CloneCells[i]->Evaluate( );
        
	#pragma omp critical
	EvalNum++;
      }

      // Evaluate the fitness function of the clone cells.
     
      if (Pen)
        Pen->EvalPenObjFunc(&CloneCells, TolViol);
      
      // Merge the clone cells with the memory cells.
      
      Merge(AllCells,(PopSize+CloneSize),MemoryCells, MemoryRenewSize);

      // Apply receptor editing operator.

      ReceptEditing(MemoryCells, ReceptEditingSize);

      // Evaluate new cells from receptor editing.
      
      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < ReceptEditingSize; i++)
      {
        MemoryCells[PopSize -1 - i]->Evaluate( );

	#pragma omp critical
        EvalNum++;
      }

      // Evaluate penalized objective function.
      
      if (Pen)
        Pen->EvalPenObjFunc(&MemoryCells, TolViol);
     
      // Migration.
      
      #ifdef _MPI_
      if (!((gen+1)%MigrationGap) && (gen+1) != MaxGen)
      {
        MemoryCells.Sort( );
        Migration(MemoryCells);
        if (Pen)
          Pen->EvalPenObjFunc(&MemoryCells, TolViol);
      }
      #endif 

      // Update variables related to PostProcessing

      UpdatePostVar(gen,opt,lastBest,&MemoryCells);
      
      // Check conditions to stop optimization

      if (OptStopCrit(gen,opt,lastBest,&MemoryCells))
        break;
    }
   
    // Store the best solution.

    best->Insert(MemoryCells.BestSol( ));
  
    // Print data in the output file.

    PrintPostVar(MaxGen,opt,EvalNum,&MemoryCells);

    // Delete some variables

    delete [] AllCells;
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= Fitness ===================================

void cStandardAIS :: Fitness(cPopulation &pop)
{
  // Compute the fitness function from the objective function.

  if (!Pen) // Unconstrained problems.
  {
    double bigval = abs(pop[0]->GetObjFunc(0));
    for (int i = 1; i < pop.GetSize( ); i++)
      bigval = MAX(bigval, pop[i]->GetObjFunc(0));

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < pop.GetSize( ); i++)
      pop[i]->AssignFitFunc(1.10*bigval - pop[i]->GetObjFunc(0));
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

// ============================= GlobalFitness =============================

void cStandardAIS :: GlobalFitness(cOptSolution** &allcells, int size)
{
  // Compute the fitness function from the objective function.

  if (!Pen) // Unconstrained problems.
  {
    double bigval = abs(allcells[0]->GetObjFunc(0));
    for (int i = 1; i < size; i++)
      bigval = MAX(bigval, allcells[i]->GetObjFunc(0));

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < size; i++)
      allcells[i]->AssignFitFunc(1.10*bigval - allcells[i]->GetObjFunc(0));
  }
  else      // Constrained problems.
  {
    double bigval = abs(allcells[0]->GetPenObjFunc( ));
    for (int i = 1; i < size; i++)
      bigval = MAX(bigval, allcells[i]->GetPenObjFunc( ));

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < size; i++)
      allcells[i]->AssignFitFunc(1.10*bigval - allcells[i]->GetPenObjFunc( ));
  }
}

// ============================ Hypermutation =============================

void cStandardAIS :: Hypermutation(cPopulation &clones, cOptSolution *bcell, cOptSolution *wcell)
{
  // Get best and worst fitness value.
  
  double bestfit = bcell->GetFitFunc( );
  double worstfit = wcell->GetFitFunc( );
  
  for (int i = 0; i < clones.GetSize( ); i++)
  {
   if (bestfit < clones[i]->GetFitFunc( ))
     bestfit = clones[i]->GetFitFunc( );
   if (worstfit > clones[i]->GetFitFunc( ))
     worstfit = clones[i]->GetFitFunc( );
  }

  // Perform the hypermutation operation.
  
  double r;

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < clones.GetSize( ); i++)
  {
    r = exp(-Decay*((clones[i]->GetFitFunc( )-worstfit)/(bestfit-worstfit)));
    clones[i]->Hypermutation(r);
  }
}

// ============================= Migration =================================

void cStandardAIS :: Migration(cPopulation &pop)
{
#ifdef _MPI_
	 return; 
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

void cStandardAIS :: Merge(cOptSolution** &allc, int s, cPopulation &mem, int r)
{
  // Sort all cells;

  sort(allc, allc + s, cOptSolution :: Compare);

  // Renew cells.

  #pragma omp parallel for num_threads(omp_maxthread)
  for(int i = 0; i < r; i++) 
    mem[i]->Copy(allc[i]);
}  

// ================================ ReceptEditing =========================

void cStandardAIS :: ReceptEditing(cPopulation &memory,int receptSize)
{
  // Get necessary data.

  int last = memory.GetSize() - 1;

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < receptSize; i++)
    memory[last - i]->Init( );
 }

// ============================= RandomRates ==============================

void cStandardAIS :: RandomRates(void)
{

}

// ======================================================= End of file =====
