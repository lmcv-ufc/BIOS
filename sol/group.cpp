// -------------------------------------------------------------------------
// group.cpp - implementation of the cGroup class.
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
// Created:      14-Mar-2013    Evandro Parente Junior
//
// Modified:     03-Oct-2014    Elias Saraiva Barroso
//               Rename to cGroup.
// -------------------------------------------------------------------------

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cfloat>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "group.h"
#include "optsolution.h"
#include "individual.h"
#include "particle.h"
#include "food.h"
#include "sampsao.h"
#include "sao.h"
#include "problem.h"
#include "utl.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cGroup ===============================

cGroup :: cGroup(int n, cProblem *prob)
{
  GroupSize = n;
  Prob = prob;
  SolVec.reserve(GroupSize);
}

// ============================ ~cGroup ====================================

cGroup :: ~cGroup(void)
{
}

// =============================== Sort ====================================

void cGroup :: Sort(void)
{
  sort(SolVec.begin( ),SolVec.end( ),cOptSolution::Compare);
}

// ============================== BestSol ==================================

cOptSolution* cGroup :: BestSol(void)
{
  return(*min_element(SolVec.begin( ), SolVec.end( ), cOptSolution::Compare));
}

// ============================== WorstSol =================================

cOptSolution* cGroup :: WorstSol(void)
{
  cOptSolution *max = GetSol(0);

  for(int i = 1; i < GroupSize; i++)
    if (*max < *GetSol(i))
      max = GetSol(i);

  return max;
}

// ============================== GetSolId =================================

int cGroup :: GetSolId(cOptSolution *solptr)
{
  for(int i = 0; i < GroupSize; i++)
    if (GetSol(i) == solptr)
      return i;
  return -1;
}

// ============================= AvgObjFunc ================================

double cGroup :: AvgObjFunc(void)
{
  double sum = 0.0;

  if (Prob->GetNumConstr() == 0)
    for (int i = 0; i < GroupSize; i++)
      sum += GetSol(i)->GetObjFunc(0);

  else
    for (int i = 0; i < GroupSize; i++)
      sum += GetSol(i)->GetPenObjFunc( );

  return(sum/GroupSize);
}

// ============================= StdDeviation ==============================

double cGroup :: StdDevObjFunc(void)
{
  double sum = 0.0;
  double mean = AvgObjFunc( );

  if (Prob->GetNumConstr() == 0)
    for (int i = 0; i < GroupSize; i++)
      sum += pow((GetSol(i)->GetObjFunc(0) - mean),2);

  else
    for (int i = 0; i < GroupSize; i++)
      sum += pow((GetSol(i)->GetPenObjFunc( ) - mean),2);

  return(sqrt(sum/GroupSize));
}

// ========================== NormMeanSqrRootErr ===========================

double cGroup :: NormMeanSqrRootErr(double minobj)
{
  double sum = 0.0;

  if (Prob->GetNumConstr() == 0)
    for (int i = 0; i < GroupSize; i++)
      sum += pow((GetSol(i)->GetObjFunc(0) - minobj)/minobj,2);

  else
    for (int i = 0; i < GroupSize; i++)
      sum += pow((GetSol(i)->GetPenObjFunc( ) - minobj)/minobj,2);

  return(sqrt(sum/GroupSize));
}

// -------------------------------------------------------------------------
// cSolGroup class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSolGroup =================================

cSolGroup :: cSolGroup(int n,cProblem *p)
{
  Prob = p;
  SolVec.reserve(n);
  GroupSize = 0;
}

// ============================ ~cSolGroup =================================

cSolGroup :: ~cSolGroup(void)
{
  for (int i = 0; i < GroupSize; i++)
    delete SolVec[i];

  SolVec.clear( );
}

// ============================= Insert ====================================

void cSolGroup :: Insert(cOptSolution *p)
{
  SolVec.push_back(p->Clone( ));
  GroupSize++;
}

// -------------------------------------------------------------------------
// cPopulation class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPopulation ===============================

cPopulation :: cPopulation(int n, eSolType tS, cProblem *p) : cGroup(n,p)
{
  // Create Individuals.

  for (int i = 0; i < n; i++)
    SolVec.push_back(cIndividual :: CreateIndividual(tS,Prob));
}

// ============================ ~cPopulation ===============================

cPopulation :: ~cPopulation(void)
{
  for (int i = 0; i < GroupSize; i++)
    delete SolVec[i];

  SolVec.clear( );
}

// -------------------------------------------------------------------------
// cSwarm class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSwarm ====================================

cSwarm :: cSwarm(int n, eSolType tS, cProblem *p) : cGroup(n,p)
{
  // Create particles.

  for (int i = 0; i < n; i++)
    SolVec.push_back(cParticle :: CreateParticle(tS,Prob));
}

// ============================ ~cSwarm ====================================

cSwarm :: ~cSwarm(void)
{
  for (int i = 0; i < GroupSize; i++)
    delete SolVec[i];

  SolVec.clear( );
}

// -------------------------------------------------------------------------
// cFoodSource class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cFoodSource ===============================

cFoodSource :: cFoodSource(int n, eSolType tS, cProblem *p) : cGroup(n,p)
{
  // Create foods.

  for (int i = 0; i < n; i++)
    SolVec.push_back(cFood :: CreateFood(tS,Prob));
}

// ============================ ~cSwarm ====================================

cFoodSource :: ~cFoodSource(void)
{
  for (int i = 0; i < GroupSize; i++)
    delete SolVec[i];

  SolVec.clear( );
}

// -------------------------------------------------------------------------
// cSampSet class:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cSampSet ==================================

cSampSet :: cSampSet(int n, eSolType tS, cProblem *p, sProbAppOut &o) : cGroup(n,p)
{
  vtype      = tS;
  ProbAppOut = &o;

  // Create samples.

  for (int i = 0; i < n; i++)
    SolVec.push_back(cSampSAO :: CreateSamp(tS,Prob,o));
}

// ============================ ~cSampSet ==================================

cSampSet :: ~cSampSet(void)
{
  for (int i = 0; i < GroupSize; i++)
    delete SolVec[i];

  SolVec.clear( );
}

// ============================= BestFeasible ==============================

cSampSAO* cSampSet :: BestFeasible( )
{
  int nc;
  cSampSAO  *best = 0;
  double currbest = DBL_MAX;
  for (int i = 0; i < GroupSize; i++)
  {
    // Check if current sample violate any constraint.
    bool next = false;
    nc = SolVec[i]->GetProb( )->GetNumConstr( );
    for (int c = 0; c < nc; c++)
      if (SolVec[i]->GetConstr( )[c] > 0.0) next = true;

    if (next) continue;

    if (SolVec[i]->GetObjFunc(0) < currbest)
    {
      best     = (*this)[i];
      currbest = best->GetObjFunc(0);
    }
  }

  return best;
}

// ============================= PushBack ==================================

cSampSAO* cSampSet :: PushBack(cVector &nv)
{
  cSampSAO* newsmp = cSampSAO :: CreateSamp(vtype,Prob,*ProbAppOut);

  newsmp->Init(nv);
  SolVec.push_back(newsmp);
  GroupSize++;

  return newsmp;
}

// ======================================================= End of file =====
