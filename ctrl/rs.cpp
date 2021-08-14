// -------------------------------------------------------------------------
// rs.cpp - implementation of cRandomSearch class.
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
// Created:      12-Sep-2019    Elias Barroso
//
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

#include "rs.h"
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
// Class cRandomSearch:
//

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cRandomSearch ============================

cRandomSearch :: cRandomSearch(void) : cOptAlgorithm( )
{
  Type = RANDOM_SEARCH;
}

// ============================= ~cRandomSearch ============================

cRandomSearch :: ~cRandomSearch(void)
{
}

// =============================== Solver ==================================

void cRandomSearch :: Solver(void)
{
  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Initialize the surrogate model.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Create the population.
    cPopulation pop(PopSize, SolType, Prob);

    // Evaluate solutions.
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < PopSize; i++)
    {
      // Initialize genotypes.
      if (i < NumInpSol)
        pop[i]->Init(InpSolVec[i]); // Input values.
      else
        pop[i]->Init( );            // Random values.

     // Evaluate the objective function and constraints of the initial pop.

      pop[i]->Evaluate( );

      #pragma omp critical
      EvalNum++;
    }

    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    // Update variables related to PostProcessing.

    UpdatePostVar(0, opt, lastBest, &pop);

    // Store the best individual.

    best->Insert(pop.BestSol( ));

    // Print data in the output file.

    PrintPostVar(1, opt, EvalNum, &pop);
  }
}

// ======================================================= End of file =====
