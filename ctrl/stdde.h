// -------------------------------------------------------------------------
// stdde.h - file containing the definition of the cStandardDE class.
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
// The cStandardDE class implements the standard Differential Evolution with
// differentiation, crossover and mutation operators. Constrained optimization
// problems are solved using the penalty function approach.
//
// -------------------------------------------------------------------------

#ifndef _STDDE_H
#define _STDDE_H

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
// Definition of the Conventional DE class:
//
class cStandardDE : public cOptAlgorithm
{
 protected:
 
          double     CrossRate;    // Crossover rate
          double     MaxCross;     // Maximum crossover rate
          double     MinCross;     // Minimum crossover rate
          double     FWeight;      // Scale factor

  virtual void      Fitness(cPopulation&);
  virtual void      Merge(cPopulation&,cPopulation&);
  virtual void      Differentiation(cPopulation&,cPopulation&, cVector&);
  virtual void      Crossover(cPopulation&,cPopulation&, cVector);
  virtual void      Migration(cPopulation&);
  virtual void      RandomRates(void);

 public:
  void              ReadCrossRate(std::istream&);
  void              ReadCrossRange(std::istream&);
  void              ReadFWeight(std::istream&);
  virtual void      LoadReadFunc(cInpMap&);
         void       SetDifType(eDifType d) {DifType = d;};

                    cStandardDE(void);
  virtual          ~cStandardDE(void);
  virtual void      Solver(void);
};

#endif
