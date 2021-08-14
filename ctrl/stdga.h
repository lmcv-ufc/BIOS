// -------------------------------------------------------------------------
// stdga.h - file containing the definition of the cStandardGA class.
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
//
// The cStandardGA class implements the standard Genetic Algorithm with
// selection, crossover and mutation operators. Constrained optimization 
// problems are solved using the penalty function approach.
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadCrossRate(void)
//  void ReadCrossRange(void)
//
// These methods read relevant data to the genetic algorithm and store them
// in protected variables.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// void Fitness(cGroup &pop)
//
//   pop - given population                                        (in/out)
//
// This method evaluates the fitness function value for all individuals
// of the given population.
// -------------------------------------------------------------------------
//
// void Mutation(cGroup &pop)
//
//   pop - given population                                        (in/out)
//
// This method applies the mutation genetic operator to all individuals
// of a given population. The mutation opertator is implemented by the 
// cIndividual class.
// -------------------------------------------------------------------------
//
// void Merge(cGroup &pop, cGroup &son)
//
//   pop     - original population / merged Population             (in)
//   son     - offspring population                                (out)
//
// This method exchanges the worst individuals from a given original popu-
// lation with individuals of an offspring population. The best individuals
// from the original population are kept in the resulting population 
// (elitism).
// -------------------------------------------------------------------------
//
// void Crossover(cGroup &parent, cGroup &son)
//
//   parent  - parent population                                   (in)
//
// This method performs the crossover operation generating a offspring 
// population. The crossover operator is implemented by the  cIndividual
// class.
// -------------------------------------------------------------------------
//
// void RandomRates(void)
//
// This method randomizes the genetic operator rates if ranges were
// specified instead of fixed values by the user.
// -------------------------------------------------------------------------

#ifndef _STDGA_H
#define _STDGA_H

#include "optalg.h"
#include "group.h"
#include "input.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cSelection;
class cIndividual;

// -------------------------------------------------------------------------
// Crossover Methods:
//
/*
typedef enum 
{
  LINEAR_COMBINATION,
  CLASSICAL
} eCrossType;*/

// -------------------------------------------------------------------------
// Definition of the Conventional GA class:
//
class cStandardGA : public cOptAlgorithm
{
 protected:

 // static eCrossType InpType;      // Crossover type read from input file
 
          double     CrossRate;    // Crossover rate
          double     MaxCross;     // Maximum crossover rate
          double     MinCross;     // Minimum crossover rate

  virtual void      Fitness(cPopulation&);
  virtual void      Mutation(cPopulation&);
//  virtual void      AddBestSample(cPopulation &, cVector, double, int);
  virtual void      Merge(cPopulation&,cPopulation&);
  virtual void      Crossover(cPopulation&,cPopulation&);
  virtual void      Migration(cPopulation&);
  virtual void      RandomRates(void);

 public:
  void              ReadCrossRate(std::istream&);
  void              ReadCrossRange(std::istream&);
  virtual void      LoadReadFunc(cInpMap&);

                    cStandardGA(void);
  virtual          ~cStandardGA(void);
  virtual void      Solver(void);
};

#endif
