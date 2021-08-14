// -------------------------------------------------------------------------
// sel.h - file containing the definition of the cSelection class.
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
// The class cSelection provides a template for the implementation of diffe-
// rent selection methods.
//
// Selection
// |-- Ranking
// |-- Fitness Proportional (aka roulette wheel)
// |-- Tournament
// |-- NSGA Tournament
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
// void ReadSelMethod(void)
//
// This method reads the selection method from the input file.
// -------------------------------------------------------------------------
//
// eSelType GetInpType(void)
//
// This method returns the selection type read from the input file.
// -------------------------------------------------------------------------
//
// cSelection *CreateSelection(eSelType type)
//
//   type - Selection type                                         (in)
//
// This method creates a selection object based on the given selection type.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// eSelType GetType(void)
//
// This method returns the selection type.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// void Select(double rate, cGroup &pop, cGroup &sel, int SelInds,
//             vec<int> Vec, vec<double> Vec2  )
//
//   rate - percentage of the original population to be selected    (in)
//   pop  - original population                                     (in)
//   sel  - select individuals                                      (out)
//   SelInds - indices of selected individuals                      (in)
//   Vec - vector of integers which is empty by default, except for
//   the NSGATournament option, when Vec1 stores the rank of all
//   individuals. This parameter is only relevant in multiobjective
//   problems                                                       (in)
//   Vec2   - vector of double which is empty by default, except for
//   the NSGATournament option, when Vec2 stores the crowding distance of
//   all individuals. Again, this paramter is only relevant in
//   multiobjective problems)                                       (in)
//
// This method takes a population and selects a certain number of
// individuals from it, based on the given selection rate and the chosen
// selection method.
// -------------------------------------------------------------------------

#ifndef _SELECTION_H
#define _SELECTION_H

#include <vector>

// ------------------------------------------------------------------------
// Forward declarations:
//
class cIndividual;
class cGroup;

// -------------------------------------------------------------------------
// Fitness Function types:
//
typedef enum 
{
  RANKING,
  FITNESS_PROPORTIONAL,
  TOURNAMENT,
  TOURNAMENTNSGAII
} eSelType;

// -------------------------------------------------------------------------
// Definition of cSelection class:
//
class cSelection
{
 protected:
          eSelType Type;      // Selection type

 public:
  static  cSelection *CreateSelection(eSelType);

                      cSelection(void);
  virtual            ~cSelection(void);
          eSelType    GetType(void) { return Type; }
  virtual void        Select(double, cGroup *, cGroup *, int *SelInds=0, std::vector<int> *vec = 0, std::vector<double> *vec2 = 0) = 0;
};


// -------------------------------------------------------------------------
// Definition of Ranking Selection class:
//
class cSelRank : public cSelection
{
 public:
                cSelRank(void);
  virtual      ~cSelRank(void);
   void        Select(double, cGroup *, cGroup *, int *SelInds=0, std::vector<int> *vec = 0, std::vector<double> *vec2 = 0);
};


// -------------------------------------------------------------------------
// Definition of Fitness Proportional Selection class:
//
class cSelFitProp : public cSelection
{
 public:
               cSelFitProp(void);
  virtual     ~cSelFitProp(void);
    void      Select(double, cGroup *, cGroup *, int *SelInds=0, std::vector<int> *vec = 0, std::vector<double> *vec2 = 0);
};

// -------------------------------------------------------------------------
// Definition of Tournament Selection class:
//
class cSelTournament : public cSelection
{
 public:
        cSelTournament(void);
       ~cSelTournament(void);
  void       Select(double, cGroup *, cGroup *, int *SelInds=0, std::vector<int> *vec = 0, std::vector<double> *vec2 = 0);
};

// -------------------------------------------------------------------------
// Definition of TournamentNSGA Selection class:
//
class cSelTournamentNSGAII : public cSelection
{
 public:
        cSelTournamentNSGAII(void);
       ~cSelTournamentNSGAII(void);

  void  Select(double, cGroup *, cGroup *, int *SelInds=0, std::vector<int> *vec = 0, std::vector<double> *vec2 = 0);
};

#endif
