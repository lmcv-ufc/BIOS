// -------------------------------------------------------------------------
// stdpso.h - file containing the definition of the cStandardPSO class.
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
// The cStandardPSO class implements the standard Particle Swarm Algorithm
// with swarm topologies (gbest, ring, square) [1,2], dynamic parameters and 
// mutate operator [3] over particle position (after velocity evaluation and
// addition). Constrained optimization problems are solved using the penalty
// function approach. 
//
// References:
//
// Ref [1] - Bratton & Kennedy (2007): Defining a Standard for Particle Swarm
//                                     Optimization. Proceedings of the 2007
//                                     IEEE Swarm Intelligence Symposium (SIS
//                                     2007).
//                                      
// Ref [2] - Kennedy & Mendes (2006) : Neighborhood Topologies in Fully Informed
//                                     and Best-of-Neighborhood Particle
//                                     Swarms. IEEE TRANSACTIONS ON SYSTEMS,
//                                     MAN, AND CYBERNETICS—PART C: APPLICATIONS
//                                     AND REVIEWS, VOL. 36, NO. 4, JULY 2006.
//
// Ref [3] - Ratnaweera et al. (2004): Self-Organizing Hierarchical Particle
//                                     Swarm Optimizer With Time-Varying
//                                     Acceleration Coefficients. IEEE
//                                     TRANSACTIONS ON EVOLUTIONARY COMPUTATION,
//                                     VOL. 8, NO. 3, JUNE 2004
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadSwarmTopology(void)
//  void ReadInertia(void)
//  void ReadC1(void)
//  void ReadC2(void)
//  
// These methods read relevant data to the particle swarm optimization 
// algorithm and store them in protected variables.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// cParticle* EvalBestNeighboor(int id, cSwarm &swarm)
//
//   id    - given particle index                                  (in)
//   swarm - given swarm                                           (in)
//
// This method evaluate and return the reference to the best particle 
// (based on its best history position) on neighboorhood of the given 
// particle index and swarm. 
// -------------------------------------------------------------------------
//
// void InitParameters(void)
//
// This method initialize pso paramaters. Considering  the standard pso, 
// these are inertia, cognitive parameter (c1) and social parameter (c2).
// -------------------------------------------------------------------------
//
// void EvalParameters(double t)
//
//   t - given time variable (current iteration)                   (in)
//
// This method evaluate pso paramaters based on the current iteration number. 
// -------------------------------------------------------------------------
//
// void EvalVelocity(cSwarm &swarm)
//
//   swarm - given swarm                                           (in/out)
//
// This method evaluate velocity of each particle, considering the swarm
// neighborhooh lerning type (CANONICAL or FULLY_INFORMED). 
// -------------------------------------------------------------------------
//
// void Mutation(cSwarm &swarm)
//
//   swarm - given swarm                                           (in/out)
//
// This method applies the mutation genetic operator to all particles
// of a given swarm. The mutation opertator is implemented by the 
// cParticle class.
// -------------------------------------------------------------------------
//
// void GetNeighborhood(int id, cSwarm &swarm, vector<int> &Inds)
//
//   id    - particle index                                        (in)
//   swarm - given swarm                                           (in)
//   Inds  - neighbors indices                                     (out)
//
// This method find the neighborhood indexes of a given particle index and
// swarm, storing then at a given int vector.
// -------------------------------------------------------------------------
//
// int CorretSwarmSize(int swarmsize)
//
//   swarmsize - swarm size                                        (in)
//
// This method evaluate the swarm size based on the correct swarm matrix
// dimensions in square topology and return it. In other topologies this method
// does nothing. 
// -------------------------------------------------------------------------

#ifndef _STDPSO_H
#define _STDPSO_H

#include "optalg.h"
#include "algparam.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cSwarm;        

// Method to read topology from file

ostream& operator<<(ostream&,const eSwaTopType&);

// -------------------------------------------------------------------------
// Swarm Neighborhood Lerning:
//
typedef enum 
{
  CANONICAL,           // classical update formula. [2]
  FULLY_INFORMED       // fully informed formula.
} eSwaNbhLen;

// -------------------------------------------------------------------------
// Matrix index:
//
typedef struct
{
  double i;
  double j;
} sMatIndex;

// -------------------------------------------------------------------------
// Definition of the Conventional PSO class:
//
class cStandardPSO : public cOptAlgorithm
{
 protected:
          eSwaTopType Topology;
          eSwaNbhLen  NbhLenType;   // Neighborhood learning type
          cAlgParam   Inertia;     
          cAlgParam   c1;
          cAlgParam   c2;
        
          double      RedFactor;    // Redirect Factor
          sMatIndex   SquareMatSize;
          int         SwarmSize; 

  virtual cParticle*  EvalBestNeighboor(int,cSwarm&);
  virtual void        InitParameters(void);
  virtual void        EvalParameters(double);
  virtual void        EvalVelocity(cSwarm&);
  virtual void        Mutation(cSwarm&);
  virtual void        Migration(cSwarm&);
  virtual void        GetNeighborhood(int,cSwarm&,vector<int>&);
  virtual int         CorrectSwarmSize(int);  
  
 public:              
          void        ReadSwarmTopology(std::istream&);
          void        SetSwarmTopology(eSwaTopType t) { Topology        = t; }
          void        ReadInertia(std::istream&);
          void        ReadC1(std::istream&);
          void        ReadC2(std::istream&);
          void        LoadReadFunc(cInpMap&);
                      
                      cStandardPSO(void);
  virtual            ~cStandardPSO(void);
  virtual void        Solver(void);
};

#endif
