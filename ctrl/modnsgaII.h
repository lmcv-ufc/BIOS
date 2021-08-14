// -------------------------------------------------------------------------
// modNSGAII.h - file containing the definition of the cmodNSGAII class.
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
// The cmodNSGAII class implements a modified version of the Nondominated Sorting
// Genetic Algorithm II developed by Deb et al. (2002).
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadExploreBounds(void)
//
// This method reads the explore boundaries rate and store it in a protected
// variable.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// void Rank(cPopulation &pop)
//
//   pop     - given population                                    (in)
//   SolRank - vector with the rank of the solutions in pop        (out)
//
// This method evaluates the rank of all individuals of the given
// population (pop) and stores it in a vector of integers (SolRank).
// -------------------------------------------------------------------------
//
// void Fronts(cPopulation &pop, int &numfronts)
//
//   pop        - given population                                 (in)
//   numfronts  - number of fronts at the end of the sorting       (out)
//
// This method sorts the given population (pop) in fronts according to
// the rank of the individuals and provides the final number of fronts
// generated (numfronts).
// -------------------------------------------------------------------------
//
// void Merge(cPopulation &pop, cPopulation &son, cPopulation &newpop)
//
//   pop     - original population                                 (in)
//   son     - offspring population                                (in)
//   newpop  - merged population (pop + son)                       (out)
//
// This method merges two populations (pop and son)of size N each into one
// single population of size 2N.
// ------------------------------------------------------------------------
//
// void CrowdDistance(cPopulation &pop, int numfronts)
//
//   pop        - given population                                 (in)
//   numfronts  - number of fronts at the end of the sorting       (out)
//   CrowdingDistance - vector that stores the crowding distance
//                      of the solutions in the given population   (out)
//
// This method evaluates the crowding distance of the individuals of the
// given population based on their fronts classification.
// ------------------------------------------------------------------------
//
// void FillNextGeneration(cPopulation &newpop, cPopulation &pop)
//
//   newpop   - given population of size 2N                       (in)
//   pop      - next generation population of size N              (out)
//
// This method fills the population of the following generation by
// adding individuals with the lowest ranks and highest crowding distances
// until it reaches the maximum size of N.
// -------------------------------------------------------------------------
//
// void CrowdDistanceCut(cPopulation &pop, int frontid)
//
//   pop      - given population                                 (in)
//   frontid  - number of the front that surpasses the maximum
//              size of the population                           (in)
//
// This method evaluates the crowding distance of the individuals of the
// given population that belong to the specified front.
// -------------------------------------------------------------------------
//
// void Crossover(cGroup &parent, cGroup &son)
//
//   parent  - parent population                                   (in)
//   son     - son population                                      (out)
//
// This method performs the crossover operation generating the offspring
// population. The crossover operator is implemented by the cIndividual
// class. It is worth mentioning that the crossover operation suggested
// by the reference authors is recommended and designed for continuous
// search space, namely Simulated Binary Crossover (SBX),  which is not
// the case of the laminate problems (discrete variables). In this
// implementation, the crossover uses either a linear combination of the
// parents or the classical crossover. Further details are found in the
// cIndividual Class.
// ------------------------------------------------------------------------
//
// void Mutation(cPopulation &pop)
//
//   pop - given population                                        (in/out)
//
// This method applies the mutation genetic operator to all individuals
// of a given population. The mutation opertator is implemented by the
// cIndividual class. At this point, another modification is considered.
// The mutation operation in the original paper is done using a polynomial
// mutation, which is not considered to in this implementation.
// -------------------------------------------------------------------------
//
// void ExploreBounds(int gen, cPopulation &pop, int NumBoundsInd)
//
//   gen  - current generation                                    (in)
//   pop  - given population                                      (in/out)
//   NumBoundsInd - number of individuals that must be modified   (in)
//
// This method performs the explore boundaries operation generating
// individuals with their upper and lower limits evaluated by changing the
// material chromosome (default). The explore boundaries operator is
// implemented by the cIndividual class.
// -------------------------------------------------------------------------
//
// void RandomRates(void)
//
// This method randomizes the genetic operator rates if ranges were
// specified instead of fixed values by the user.
// -------------------------------------------------------------------------

#ifndef _modNSGAII_H
#define _modNSGAII_H

#include "optalg.h"
#include "group.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cSelection;
class cIndividual;

// -------------------------------------------------------------------------
// Auxiliary types:
//

typedef struct
{
    vector<int>    Front;          // Vector of solution addresses at each front
    int            FrontSize;      // Number of elements in each front (size)
    vector<double> CrowdDistance;  // Vector of crowding distance value of the solutions
}sFrontClass;

// -------------------------------------------------------------------------
// Definition of the NSGA class:
//
class cmodNSGAII : public cOptAlgorithm
{
 protected:
         double     ExpBoundsRate;   // Explore bounds rate
     sFrontClass*   FrontClass;      // Vector of Fronts

  virtual void      Mutation(cPopulation&);
  virtual void      Rank(cGroup *);
  virtual void      Fronts(cPopulation&, int &);
  virtual void      FillNextGen(cPopulation&, cPopulation&);
  virtual void      CrowdDistance(cPopulation&, int);
  virtual void      CrowdDistanceCut(cPopulation&, int);
  virtual void      Merge(cPopulation&, cPopulation&,cPopulation&);
  virtual void      Crossover(cPopulation&,cPopulation&);
  virtual void      ExploreBounds(int, cPopulation&, int);
  virtual void      Migration(cPopulation&);
  virtual void      RandomRates(void);
          int       checkdominance (cOptSolution *a, cOptSolution *b);
   vector<int>      SolRank;            // Vector of integers to store the rank of the individuals
   vector<double>   CrowdDistanceValue; // Vector of double to store the crowding distance values

 public:
          void      ReadExploreBounds(std::istream&);
  virtual void      LoadReadFunc(cInpMap&);

                    cmodNSGAII(void);
  virtual          ~cmodNSGAII(void);
  virtual void      Solver(void);
};

#endif
