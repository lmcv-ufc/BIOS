// -------------------------------------------------------------------------
// stdabc.cpp - implementation of cStandardABC class.
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
// Created:      02-Jul-2014    Elias Saraiva Barroso
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>  
#include <cmath>
#include <sstream> 
#include <map>
#include <algorithm>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "stdabc.h"
#include "problem.h"
#include "sel.h"
#include "group.h"
#include "food.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadScoutBeesRate ========================

void cStandardABC :: ReadScoutBeesRate(istream &in)
{
  if (!(in >> ScoutBeesRate)) 
  {
    cout << "Error in the input of the scout bees number." << endl;
    exit(0);
  }
}

// ============================== ReadFoodMaxTrialRate =====================

void cStandardABC :: ReadFoodMaxTrialRate(istream &in)
{
  if (!(in >> MaxTrialRate)) 
  {
    cout << "Error in the input of the max trial number (ABC)." << endl;
    exit(0);
  }
}

// ============================== ReadFoodSourceRate =======================

void cStandardABC :: ReadFoodSourceRate(istream &in)
{
  if (!(in >> FoodSourceRate)) 
  {
    cout << "Error in the input of the food source ratio." << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cStandardABC :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("FOOD.MAXTRIAL.RATE", makeReadObj(cStandardABC,ReadFoodMaxTrialRate));
  im.Insert("CROSSOVER.RATE",     makeReadObj(cStandardABC,ReadFoodSourceRate));
  im.Insert("SCOUTBEES.RATE",     makeReadObj(cStandardABC,ReadScoutBeesRate));
}

// ============================== cStandardABC =============================

cStandardABC :: cStandardABC(void)
{
  Type            = STANDARD_ABC;
  Prob            = 0;
  Pen             = 0;
  Sel             = 0;

  // Default parameters.
  MaxTrialRate    = 0.10;
  ScoutBeesRate   = 0.05;
  FoodSourceRate  = 0.50;
}

// ============================= ~cStandardABC =============================

cStandardABC :: ~cStandardABC(void)
{
}

// =============================== Solver ==================================

void cStandardABC :: Solver(void)
{
  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Track number of individual evaluations.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Randomize rates.

    if (!MutProb && MaxMut)
      MutProb = Utl::RandDouble(MinMut, MaxMut);

    // Set ABC parameters.

    int FoodsNumber = round(PopSize * FoodSourceRate);
    int OnlookerBeesNumber = round((PopSize - FoodsNumber)*(1 - ScoutBeesRate));
    int MaxTrial = round(MaxTrialRate*PopSize);
    int ScoutBeesNumber = PopSize - FoodsNumber - OnlookerBeesNumber;

    if (FoodsNumber <= 0) FoodsNumber = 1;
    if (OnlookerBeesNumber <= 0) OnlookerBeesNumber = 1;

    // Create populations for Foods, Solutions for EmployedBees and OnlookerBees.
   
    cFoodSource Foods(FoodsNumber, SolType, Prob);
    cFoodSource SolEmployedBees(FoodsNumber, SolType, Prob);
    cFoodSource SolOnlookerBees(OnlookerBeesNumber, SolType, Prob);

    int *SelIndices = new int [OnlookerBeesNumber];
    
    // Generate the initial population.
    
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < FoodsNumber; i++)
    {
      // Initialize foods.     
      if (i < NumInpSol)
        Foods[i]->Init(InpSolVec[i]); // Input values.
      else
        Foods[i]->Init( );            // Random values.
      
      // Evaluate the objective function and constraints of the initial foods source.
      
      Foods[i]->Evaluate( );

      #pragma omp critical
      EvalNum++;
    }

    if (Pen)
      Pen->EvalPenObjFunc(&Foods, TolViol);

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    for (int gen = 0; gen < MaxGen; gen++)
    {
      if ((gen+1)%5 == 0 && Feedback) cout << "Generation: " << gen + 1 << endl;

      // Phase 1 - Send Employed Bees to search for new solution.
      
      SendEmployedBees(Foods,SolEmployedBees);

      // Evaluate the mutate solutions.

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int j=0; j < FoodsNumber; j++)
      {
          SolEmployedBees[j]->Evaluate( );
          
	  #pragma omp critical
	  EvalNum++;
      }

      // Evaluate the penalized objective function of the population.

      if (Pen)
        Pen->EvalPenObjFunc(&SolEmployedBees, TolViol);

      // Update Food Position and FoodTrials.

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int j=0; j < FoodsNumber; j++)
      {
        if (*SolEmployedBees[j] < *Foods[j])
	{
          Foods[j]->Copy(SolEmployedBees[j]);
	  Foods[j]->ResetCounter( );
	}

	else
	  Foods[j]->IncrementCounter( );
      }
	  
      // Phase 2 - Send Onlooker Bees to search for new solutions.
     
      // Select solutions for the Onlookers Beers.
 
      Fitness(Foods);
      Sel->Select(1.0,&Foods,&SolOnlookerBees,SelIndices);

      SendOnlookerBees(Foods,SolOnlookerBees);

      // Evaluate the mutate solutions;

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int j=0; j < OnlookerBeesNumber; j++)
      {
          SolOnlookerBees[j]->Evaluate( );
          EvalNum++;
      }

      // Evaluate the penalized objective function of the new solutions.

      if (Pen)
        Pen->EvalPenObjFunc(&SolOnlookerBees, TolViol);

      // Update Food Position and FoodTrials.

      for (int j=0; j < OnlookerBeesNumber; j++)
      {
        for(int w=0; w < FoodsNumber; w++)       
	  if (*SolOnlookerBees[j] < *Foods[SelIndices[j]])
	  {
	    Foods[SelIndices[j]]->Copy(SolOnlookerBees[j]);
	    Foods[SelIndices[j]]->ResetCounter( );
	  }
	  else
	    Foods[SelIndices[j]]->IncrementCounter( );
      } 
      
      // Phase 3 - Send Scout Bees.

      // Select bad foods (violated the max trial)

      vector<cOptSolution*> VecScouts;  
      VecScouts.reserve(FoodsNumber);
      
      for (int j=0; j < FoodsNumber; j++)
        if (Foods[j]->GetCounter( ) > MaxTrial)
          VecScouts.push_back(Foods[j]);
     
      // Sort the selected bad foods.

      sort(VecScouts.begin(),VecScouts.end(),cOptSolution :: Compare);

      int scoutId = VecScouts.size();

      cOptSolution* BestSol = Foods.BestSol( );

      // Send ScoutBees to search for new food source.

      for(int j=0; j < ScoutBeesNumber; j++)
      {
        scoutId--;

        if(scoutId < 0)
          break;

        if (VecScouts[scoutId] != BestSol)
        {
	  VecScouts[scoutId]->Init( );
	  VecScouts[scoutId]->Evaluate( );

	  #pragma omp critical
	  EvalNum++;
        }
      }

      // Evaluate the penalized objective function of the population.

      if (Pen)
        Pen->EvalPenObjFunc(&Foods, TolViol);
 
      // Migration.
      
      #ifdef _MPI_
      if (!((gen+1)%MigrationGap) && (gen+1) != MaxGen)
      {
        Foods.Sort( );
        Migration(Foods);
        if (Pen)
          Pen->EvalPenObjFunc(Foods, TolViol);
      }
      #endif

      // Update variables related to PostProcessing

      UpdatePostVar(gen,opt,lastBest,&Foods);
      
      // Check conditions to stop optimization

      if (OptStopCrit(gen,opt,lastBest,&Foods))
        break;
    }
   
    // Store the best individual.

    best->Insert(Foods.BestSol( ));
  
    // Print data in the output file.

    PrintPostVar(MaxGen,opt,EvalNum,&Foods);

    // Delete some variables
    delete [] SelIndices; 
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================= Fitness ===================================

void cStandardABC :: Fitness(cFoodSource &pop)
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

// ============================ SendEmployedBees ==========================

void cStandardABC :: SendEmployedBees(cFoodSource &foods, cFoodSource &sol)
{
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < foods.GetSize( ); i++)
  {
    // Copy foods to solutions.
   
    sol[i]->Copy(foods[i]);

    // Evaluate ABC mutant solution.

    sol[i]->ABCmutate(foods,i);
  }
}

// ============================ SendOnlookerBees ==========================

void cStandardABC :: SendOnlookerBees(cFoodSource &foods, cFoodSource &sol)
{
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < sol.GetSize( ); i++)
  {
    // Evaluate ABC mutant solution.

    sol[i]->ABCmutate(foods,i);
  }
}

// ============================= Migration =================================

void cStandardABC :: Migration(cFoodSource &pop)
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

// ======================================================= End of file =====
