// -------------------------------------------------------------------------
// group.h - file containing the definition of the cGroup class.
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
// The cGroup class stores the agents (ex: individuals or particles) used by the 
// optimization algorithms providing a simple interface to handle large 
// number of abstract objects (agents).
//  
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// int GetNumInd(void)
//
// This method returns the number of individuals of the population.
// -------------------------------------------------------------------------
//
// cIndividual* &operator[](int i)
//
//   i - index                                                     (in)
//    
// This method returns a pointer to the individual of index i.
// -------------------------------------------------------------------------
//
// void Sort(void)
//
// This method sorts the population based on the less than (<) operator of
// the cIndividual class. After sorting, the best individuals are placed on
// the first positions of the population.
// -------------------------------------------------------------------------
//
// cIndividual* BestInd(void)
//
// This method returns a pointer to the best individual of the population
// based on the less than (<) operator of the cIndividual class.
// -------------------------------------------------------------------------
//
// double AvgObjFunc(void)
//
// This method computes and returns the average value of the objective
// function of the individuals belonging to the population.
// -------------------------------------------------------------------------

#ifndef _GROUP_H
#define _GROUP_H

#include <vector>
using namespace std;

#include "optsolution.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class sProbAppOut;

// -------------------------------------------------------------------------
// Definition of cGroup class:
//
class cGroup
{
  public:
  vector<cOptSolution*> SolVec;
  int       GroupSize;          // Number of optimization solutions. (group size)
  cProblem* Prob;               // Problem
  eSolType  vtype; 


  public:

                cGroup(void) { }
                cGroup(int, cProblem *);
  virtual      ~cGroup(void);
  int           GetSize(void) { return GroupSize; }	      
  void          Sort(void);
  cOptSolution* BestSol(void);
  cOptSolution* WorstSol(void);
  cOptSolution* MeanSol(void);

  double        StdDevObjFunc(void);
  int           GetSolId(cOptSolution*);
  double        AvgObjFunc(void);
  double        NormMeanSqrRootErr(double);

  cOptSolution* GetSol(int i) {return SolVec[i];}
};

// -------------------------------------------------------------------------
// Definition of cSolGroup class:
//
class cSolGroup : public cGroup
{
 public:

                cSolGroup(int, cProblem*);
  virtual      ~cSolGroup(void);
  void          Insert(cOptSolution*);
};

// -------------------------------------------------------------------------
// Definition of cPopulation class:
//
class cPopulation : public cGroup
{
 public:
                cPopulation(int, eSolType, cProblem *);
  virtual      ~cPopulation(void);
  cIndividual*  operator[](int i) {return (cIndividual*) SolVec[i];}
};

// -------------------------------------------------------------------------
// Definition of cSwarm class:
//

class cSwarm : public cGroup
{
 public:
                cSwarm(int, eSolType, cProblem *);
               ~cSwarm(void);
  cParticle*   operator[](int i) {return (cParticle*) SolVec[i];}
  int        GetSolId(cOptSolution*);
};

// -------------------------------------------------------------------------
// Definition of cFoodSource class:
//

class cFoodSource : public cGroup
{
 public:
                cFoodSource(int, eSolType, cProblem *);
               ~cFoodSource(void);
  cFood*   operator[](int i) {return (cFood*) SolVec[i];}
};

// -------------------------------------------------------------------------
// Definition of cSampSet class:
//

class cSampSet : public cGroup
{
 protected:
  sProbAppOut* ProbAppOut;   

 public:
                cSampSet(int, eSolType, cProblem *, sProbAppOut&);
               ~cSampSet(void);
  cSampSAO*   operator[](int i) {return (cSampSAO*) SolVec[i];}
  cSampSAO*   BestFeasible(void);
  cSampSAO*   PushBack(cVector&);
};


#endif
