// -------------------------------------------------------------------------
// stdabc.h - file containing the definition of the cStandardABC class.
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
// The cStandardABC class implements the standard Artificial Bee Colony 
// Algorithm presented in "AN IDEA BASED ON HONEY BEE SWARM FOR NUMERICAL
// OPTIMIZATION" (D. Karaboga,2005). This method has been developed for 
// continuous problems. Constrained optimization problems are solved using
// the penalty function approach.            
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadFoodSourceRatio(void)
//  void ReadFoodScoutBeesRatio(void)
//  void ReadFoodMaxTrial(void)
//
// These methods read relevant data to the artificial bee colony algorithm 
// and store them in protected variables.
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
// void SendEmployedBees(cGroup &food, cGroup& solution)
//
//   food - given food source population                           (in)
//   sol  - solutions generated.                                   (out)
//
// This method evaluate a population of solutions, generated from ABC mutate 
// operator, to all foods of a given food source population. The ABC mutate 
// opertator is implemented by the cIndividual class.
// -------------------------------------------------------------------------
//
// void SendOnlookerBees(cGroup &food, cGroup& solution)
//
//   food - given food source population                           (in)
//   sol  - solutions generated.                                   (out)
//
// This method evaluate a population of solutions, generated from ABC mutate 
// operator, to all foods of a given food source population. The ABC mutate 
// opertator is implemented by the cIndividual class.
// -------------------------------------------------------------------------

#ifndef _STDABC_H
#define _STDABC_H

#include "optalg.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cFoodSource;

// -------------------------------------------------------------------------
// Definition of the Conventional Artificial Bee Colony Optimization class:
//
class cStandardABC : public cOptAlgorithm
{
 protected:
          double   FoodSourceRate;
          double   ScoutBeesRate;
          double   MaxTrialRate;

  virtual void     Fitness(cFoodSource&);
  virtual void     SendEmployedBees(cFoodSource&,cFoodSource&);
  virtual void     SendOnlookerBees(cFoodSource&,cFoodSource&);
  virtual void     Migration(cFoodSource&);

 public:
          void    ReadFoodSourceRate(std::istream&);
          void    ReadFoodMaxTrialRate(std::istream&);
          void    ReadScoutBeesRate(std::istream&);
  virtual void    LoadReadFunc(cInpMap&);
 
                  cStandardABC(void);
  virtual        ~cStandardABC(void);
  virtual void    Solver(void);
};

#endif
