// -------------------------------------------------------------------------
// modpso.h - file containing the definition modified PSO classes.
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
// The cLaminatedPSO class is an extension of the Standard Particle Swarm
// Optimization with operators specific to laminated structure problems.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// void Mutation(cPopulation &pop)
//
//   pop - given population                                        (in/out)
//
// This method applies the modified mutation operator to all individuals
// of a given population. This customized version contains not only the
// standard mutation, but also the layer swap, addition and deletion
// operators.
// -------------------------------------------------------------------------
//
// void RandomRates(void)
//
// This method generates randomized values for the  genetic operator
// rates and probabilities. This is also implemented here to randomize the
// probabilities of layer swap, addition and deletion.
// -------------------------------------------------------------------------

#ifndef _MODPSO_H
#define _MODPSO_H

#include "stdpso.h"
#include "lamalg.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cSelection;
class cIndividual;

// -------------------------------------------------------------------------
// Definition of the Laminate PSO class:
//
class cLaminatePSO : public cStandardPSO, public cLamAlg
{
 protected:
  virtual void    Mutation(cSwarm &);
  virtual void    RandomRates(void);

 public:
                  cLaminatePSO(void);
  virtual        ~cLaminatePSO(void);
};

#endif
