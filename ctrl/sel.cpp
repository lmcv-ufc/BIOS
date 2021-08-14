// -------------------------------------------------------------------------
// selection.cpp - implementation of the cSelection class.
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
// Created:      21-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     13-Mar-2013    Evandro Parente Junior
//               New definition of Select method.
// 
//               19-Oct-2015    Elias Saraiva Barroso
//               Define cSelection class methods considering the new cGroup class
//               hierarchy.
//
//               30-Aug-2017    Marina Alves Maia
//               Creation of new type of selection (TournamentNSGA).
// -------------------------------------------------------------------------

#include <iostream>
#include <cmath>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "sel.h"
#include "group.h"
#include "optsolution.h"
#include "input.h"
#include "utl.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ CreateSelection ============================

cSelection *cSelection :: CreateSelection(eSelType type)
{
  cSelection *sel = 0;
  
  switch(type)
  {
    case RANKING:
      sel = new cSelRank( );
    break;
      
    case FITNESS_PROPORTIONAL:
      sel = new cSelFitProp( );
    break;

    case TOURNAMENT:
      sel = new cSelTournament( );
    break;

    case TOURNAMENTNSGAII:
      sel = new cSelTournamentNSGAII( );
    break;
  }
  
  return(sel);
}

// ============================== cSelection ===============================

cSelection :: cSelection(void)
{
}

// ============================== ~cSelection ==============================

cSelection :: ~cSelection(void)
{
}

// -------------------------------------------------------------------------
// Class cSelRank:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cSelRank =================================

cSelRank :: cSelRank(void)
{
  Type = RANKING;
}

// ============================= ~cSelRank =================================

cSelRank :: ~cSelRank(void)
{
}

// =============================== Select ==================================

void cSelRank :: Select(double rate, cGroup *group, cGroup *sel, int *SelInds, std::vector<int> *vec, std::vector<double> *vec2)
{
  // Sort the original solution group.
  
  group->Sort( );
  
  // Calculate the sum of rankings.

  int selsize = sel->GetSize( );
  int groupsize = group->GetSize( );
  double sumrank = (groupsize+1)*groupsize/2.0;
  
  // Initialize roulette.
  
  cVector posrou(groupsize+1);
  posrou.Zero( );
  for (int i = 0; i < groupsize; i++)
    posrou[i+1] = posrou[i] + double(i + 1)/sumrank;

  // Select individuals until the new group is full.
  
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < selsize; i++)
  {
    double roll = Utl::RandDec( );

    bool ok = false;
    for (int j = 0; j < groupsize; j++)
    {
      if (roll >= posrou[groupsize-1-j] && roll <= posrou[groupsize-j])
      {
        sel->GetSol(i)->Copy(group->GetSol(j));
	ok = true;
	if (SelInds != NULL) // Elias (04-Nov-2014)
	  SelInds[i] = j; 
        break;
      }
    }
    if (!ok)  // Evandro (09-Mar-2013)
    {  
      cout << "SelRank" << endl;
      cout << "  groupsize = " << groupsize;
      cout << "  selsize = " << selsize << endl;
      int k = Utl::RandInt(0, groupsize-1);
      sel->GetSol(i)->Copy(group->GetSol(k));
      
      if (SelInds != NULL);
        SelInds[i] = k;     // Elias (04-Nov-2014)
    }	
  }
}

// -------------------------------------------------------------------------
// Class cSelFitProp:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSelFitProp ================================

cSelFitProp :: cSelFitProp(void)
{
  Type = FITNESS_PROPORTIONAL;
}

// ============================= ~cSelRank =================================

cSelFitProp :: ~cSelFitProp(void)
{
}

// =============================== Select ==================================

void cSelFitProp :: Select(double rate, cGroup *group, cGroup *sel, int *SelInds, std::vector<int> *vec, std::vector<double> *vec2)
{
  int selsize = sel->GetSize( );
  int groupsize = group->GetSize( );

  // Calculate the sum of fitness functions.
  
  double sumfit = 0.0;
  for (int i = 0; i < groupsize; i++)
    sumfit += group->GetSol(i)->GetFitFunc( );
  
  // Initialize roulette.
  
  cVector posrou(groupsize+1);
  posrou.Zero( );
  for (int i = 0; i < groupsize; i++)
    posrou[i+1] = posrou[i] + group->GetSol(i)->GetFitFunc( )/sumfit;
  
  // Select individuals until the new group is full.

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < selsize; i++)
  {
    double roll = Utl::RandDec( );
    
    bool ok = false;
    for (int j = 0; j < groupsize; j++)
    {
      if (roll >= posrou[j] && roll <= posrou[j+1])
      {
        sel->GetSol(i)->Copy(group->GetSol(j));
	ok = true;
	if (SelInds != 0) // Elias (04-Nov-2014)
	  SelInds[i] = j; 
        break;
      }
    }
    if (!ok)  // Evandro (09-Mar-2013)
    {  
      cout << "SelFitProp" << endl;
      cout << "  groupsize = " << groupsize;
      cout << "  selsize = " << selsize << endl;
      int k = Utl::RandInt(0, groupsize-1);
      sel->GetSol(i)->Copy(group->GetSol(k));
      
      if (SelInds != 0);
        SelInds[i] = k;     // Elias (04-Nov-2014)
    }	
  }
}

// -------------------------------------------------------------------------
// Class cSelTournament:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cSelTournament ==============================

cSelTournament :: cSelTournament(void)
{
  Type = TOURNAMENT;
}

// ============================= ~cSelTournament ==============================

cSelTournament :: ~cSelTournament(void)
{
}

// =============================== Select ==================================

void cSelTournament :: Select(double rate, cGroup *pop, cGroup *sel, int *SelInds, std::vector<int> *vec, std::vector<double> *vec2)
{
  int selsize = sel->GetSize( );
  int popsize = pop->GetSize( );

  //Define the number of players.

  int numplayers = 2;

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int j = 0; j < selsize; j++)
  {
    int winner = 0;
    double fit = -1.0;

    for (int l = 0; l < numplayers; l++)
    {
      double coord = Utl::RandDec( );
      double value = round((popsize-1)*coord);
      double aux = pop->GetSol(value)->GetFitFunc( );

      if (aux >= fit)
      {
        fit = aux;
        winner = value;
      }
    }

    sel->GetSol(j)->Copy(pop->GetSol(winner));


    if (SelInds != 0);
      SelInds[j] = winner;  // Elias (04-Nov-2014)
   }
}

// -------------------------------------------------------------------------
// Class cSelTournamentNSGAII:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cSelTournamentNSGA ==============================

cSelTournamentNSGAII :: cSelTournamentNSGAII(void)
{
  Type = TOURNAMENTNSGAII;
}

// ============================= ~cSelTournamentNSGA ==============================

cSelTournamentNSGAII :: ~cSelTournamentNSGAII(void)
{
}

// =============================== Select ==================================

void cSelTournamentNSGAII :: Select(double rate, cGroup *pop, cGroup *sel, int *SelInds, std::vector<int> *SolRank, std::vector<double> *CrowdDistanceValue)
{
  int selsize = sel->GetSize( );
  int popsize = pop->GetSize( );

  //Define the number of players.

  int numplayers = 2;
  int winner;
  int currentrank;
  double currentcrowddistance;

  #pragma omp parallel for num_threads(omp_maxthread)


  for (int j = 0; j < selsize; j++)
  {

    // Variables initialization

    int ranktemp = 100000;
    double crowddistancetemp = -100000;
    int valueanterior = -1000;

   //   cout << "\nSelecao entre: " << endl;

    for (int l = 0; l < numplayers; l++)
    {
      double coord = Utl::RandDec( );
      int value = round((popsize-1)*coord);

      while (value == valueanterior)
      {
          coord = Utl::RandDec( );
          value = round((popsize-1)*coord);                         // different solutions must be compared
      }

      currentrank = (*SolRank)[value];
      currentcrowddistance = (*CrowdDistanceValue)[value];

      if (currentrank < ranktemp)
      {
       ranktemp = currentrank;
       crowddistancetemp = currentcrowddistance;
       winner = value;
      }
      else if (currentrank == ranktemp)                         // if ranks are equal, compare crowding distance
      {
          if (currentcrowddistance > crowddistancetemp)
          {
              crowddistancetemp = currentcrowddistance;
              winner = value;
          }
          else if (currentcrowddistance == crowddistancetemp)  // if ranks and CDs are equal, compare the feasibility of the individuals
          {
              if (pop->GetSol(value)->GetNormConst( ) == 0 && pop->GetSol(valueanterior)->GetNormConst( ) == 0)
              {
                  if (pop->GetSol(value)->GetObjFunc(0) == pop->GetSol(valueanterior)->GetObjFunc(0) && pop->GetSol(value)->GetObjFunc(1) == pop->GetSol(valueanterior)->GetObjFunc(1))
                  {
                      if (pop->GetSol(value)->GetFeasibleConst( ) > pop->GetSol(valueanterior)->GetFeasibleConst( ))
                      {
                          winner = value;
                      }
                      else
                      {
                          winner = valueanterior;
                      }
                  }
                  else
                  {
                      winner = value;
                  }
              }
              else                                          // if ranks and CDs are equal and individuals are unfeasible,
              {                                             // choose winner randomly
                  int aux = Utl::RandInt(0, 1);
                  if (aux == 0)
                  {
                      winner = value;
                  }
                  else if (aux == 1)
                  {
                      winner = valueanterior;
                  }
              }
          }
      }

      valueanterior = value;
    }

    sel->GetSol(j)->Copy(pop->GetSol(winner));
   }

}

// ======================================================= End of file =====
