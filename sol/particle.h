// -------------------------------------------------------------------------
// particle.h - file containing the definition of the cParticle class.
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
// The cParticle define the particle swarm optimization agent.  The main
// difference of particle and its base class is the consideration of best history
// position. These positions are understood as design variables, and are
// definend in subclasses.
//
// Particle
// |-- VecInt
// |-- MatInt
// |-- VecDouble
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// void UpdateBestPos(bool force = false)
//
//   force - flag to skip the position verification                (in)
//
// this method update the best position if the best position objective function
// is worse than the current position objective function (consider penalized
// objective function in the case of constrained problems). The verification can
// be skipped using a true value on force flag.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// void UpdatePos(void)
//
// This method modify the particle position by adding particle velocity to it.
// At the end of the process, if the particle's position is outside the search
// space, the particle's position is changed to the search space limit and the
// particle's velocity is multiplied by -0.5 (only in variables that violated
// the search space limits).
// -------------------------------------------------------------------------
//
// void UpdateBestVar(void)
//
// This method update the particle best position variables. It is used by
// UpdateBestPos method. 
// -------------------------------------------------------------------------
//
// void EvalVelocity(double w, int num, double c[], cParticle *Pvec[])
//
//   w    - particle inertia                                       (in)
//   num  - number of vectors in the linear combination            (in)
//   c    - vector with the coeficients of linear combination      (in)
//   Pvec - array with pointer to particles used                   (in)
//
// This method evaluate particle velocity, as a linear combination of this
// particle velocity and current position, and particles best position (in
// Pvec): V = Vw + sum(c{i} * Pvec{i}.BestPos)  -> i in [ 0 , num-1 ]
// -------------------------------------------------------------------------
//
// void Mutate(double p)
//
//   p - mutation probability                                      (in)
//
// This method performs the mutation of the particle position according to the
// given probability.
// -------------------------------------------------------------------------

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>
using namespace std;

#include "vec.h"
#include "optsolution.h"
#include "stdpso.h"


// -------------------------------------------------------------------------
// Definition of cParticle class:
//
class cParticle : public cOptSolution
{
 protected:
         double    BestObjFuncVal;    // History best Objective function value
         double    BestPenObjFunc;    // History best Penalized objective function value
         double    BestFitFuncVal;    // History best Fitness function value
         cVector   BestConstr;        // History best Constraint values

 public:
  static void        ReadIndType(void);
  static cParticle*  CreateParticle(eSolType,cProblem*);
  static bool        Compare(cParticle*,cParticle*);


                     cParticle(cProblem*);
  virtual            ~cParticle(void);

          double    GetObjFunc(int id)  { return BestObjFuncVal; }
          double    GetPenObjFunc(void) { return BestPenObjFunc; }
          cVector&  GetConstr(void)     { return BestConstr;     }
       
  virtual void      Mutate(double)=0;
  virtual void      EvalVelocity(double,int,double[],cParticle*[]) = 0;
  virtual void      UpdatePos(void) = 0;
          void      UpdateBestPos(bool force = false);
  virtual void      UpdateBestVar(void) = 0;
};


// -------------------------------------------------------------------------
// Definition of cIndivIntVec class:
//
class cPartIntVec : public cParticle
{
 protected:
  int     *PosVar;       // Current optimization variables
  int     *BestPosVar;   // Best optimization variables
  cVector  VelVar;       // Optimization variables velocity

  cOptSolution*  NewObj(void) {return new cPartIntVec(Prob);}
 
 public:
                cPartIntVec(cProblem *);
               ~cPartIntVec(void);
  void          Init(void);
  void          Init(const sInpSol&);
  void          Init(const cVector&);
  void          Evaluate(void);
  void          Print(void);
  void          Write(std::ostream&);
  void          Copy(cOptSolution*);
  void          GetVar(int *);
  void          GetBestVar(int *);
  bool          CompVar(cOptSolution*);
  void          GetNormVar(cVector &xn);
  
  void          Mutate(double); 
  void          EvalVelocity(double,int,double[],cParticle*[]);
  void          UpdatePos(void);
  void          UpdateBestVar(void);

  void     Send(void);
  void     Receive(int);
};

// -------------------------------------------------------------------------
// Definition of cIndivIntMat class:
//
class cPartIntMat : public cParticle
{
 protected:
  int     **PosVar;       // Current optimization variables
  int     **BestPosVar;   // Best optimization variables
  cMatrix  VelVar;       // Optimization variables velocity
  
  cOptSolution*  NewObj(void) {return new cPartIntMat(Prob);}

 public:
                cPartIntMat(cProblem *);
               ~cPartIntMat(void);
  void          Init(void);
  void          Init(const sInpSol&);
  void          Init(const cVector&);
  void          Evaluate(void);
  void          Print(void);
  void          Write(std::ostream&);
  void          Copy(cOptSolution*);
  bool          CompVar(cOptSolution*);
  int **        GetVar(void) {return PosVar;}
  int **        GetBestVar(void) {return BestPosVar;}
  void          GetNormVar(cVector &xn);

  void          Mutate(double);
  void          LamMutate(double*);
  void          EvalVelocity(double,int,double[],cParticle*[]);
  void          UpdatePos(void);
  void          UpdateBestVar(void);
  void          Swap(double);
  void          Add(double);
  void          Delete(double);
  void          Send(void);
  void          Receive(int);
};

// -------------------------------------------------------------------------
// Definition of cPartDblVec class:
//
class cPartDblVec : public cParticle
{
 protected:
  cVector  PosVar;   // Current Optimization variables (vector of doubles)
  cVector  BestPosVar;   // Best optimization variables
  cVector  VelVar;       // Optimization variables velocity
  
  cOptSolution*  NewObj(void) {return new cPartDblVec(Prob);}
 
 public:
                cPartDblVec(cProblem *);
               ~cPartDblVec(void);
  void          Init(void);
  void          Init(const sInpSol&);
  void          Init(const cVector&);
  void          Evaluate(void);
  void          Print(void);
  void          Write(std::ostream&);
  void          Copy(cOptSolution *);
  bool          CompVar(cOptSolution*);
  void          GetNormVar(cVector &xn);
  cVector       GetVec(void) {return PosVar;}
  cVector       GetBestVec(void) {return BestPosVar;}

  void          Mutate(double);
                
  void          EvalVelocity(double,int,double[],cParticle*[]);
  void          UpdatePos(void);
  void          UpdateBestVar(void);

  void          Send(void);
  void          Receive(int);
};

#endif
