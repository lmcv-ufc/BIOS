// -------------------------------------------------------------------------
// individual.h - file containing the definition of the cIndividual class.
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
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// void Mutate(double p)
//
//   p - mutation probability                                      (in)
//
// This method performs the mutation of the individual according to the
// given probability.
// -------------------------------------------------------------------------
//
// void Crossover(eCrossType type, double r, cIndividual* par1, cIndividual* par2)
//
//   type - crossover method                                       (in)  
//   r    - coefficient of first parent                            (in)
//   par1 - pointer to the 1st parent                              (in)
//   par2 - pointer to the 2nd parent                              (in)
//
// This method performs the crossover operation computing the variables of
// the individual as a linear combination (ind = r*par1 + (1-r)*par2) or as
// the classical way (FALTA DESCREVER!)
// -------------------------------------------------------------------------

#ifndef _INDIVIDUAL_H
#define _INDIVIDUAL_H

#include <vector>
using namespace std;

#include "optsolution.h"
#include "stdga.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Definition of cIndividual class:
//
class cIndividual : public cOptSolution
{
protected:
        // Those  variable are used only in DE
        double    BestObjFuncVal;    // History best Objective function value
        double    BestPenObjFunc;    // History best Penalized objective function value
        double    BestFitFuncVal;    // History best Fitness function value
        cVector   BestConstr;        // History best Constraint values

 public:
  static cIndividual* CreateIndividual(eSolType,cProblem*);


                    cIndividual(cProblem*);
  virtual          ~cIndividual(void);
          bool      operator< (cIndividual&);
  virtual void      ExploreBounds(void) = 0;
  virtual void      ExploreUpperBounds(void) = 0;
   
  virtual void      Mutate(double) = 0;
  virtual void      AddBestSample(cVector, double) = 0;
  virtual void      GetBestSample(cVector &) = 0;

  virtual void      Hypermutation(double) = 0;
  virtual void      Crossover(eCrossType, double, cOptSolution*, cOptSolution*) = 0;
  virtual void      Differentiation(eDifType, double, cOptSolution*, cOptSolution*, cOptSolution*, cOptSolution*) = 0;
};

// -------------------------------------------------------------------------
// Definition of cIndivIntVec class:
//
class cIndivIntVec : public cIndividual
{
 protected:
  int     *Var;   // Optimization variables
  int     *Ib;

  cOptSolution*  NewObj(void) {return new cIndivIntVec(Prob);}

 public:
           cIndivIntVec(cProblem *);
          ~cIndivIntVec(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     Init(const cVector&);
  void     AddBestSample(cVector, double);
  void     GetBestSample(cVector &);
  void     Evaluate(void);
  void     Print(void);
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  void     GetNormVar(cVector &xn);
  bool     CompVar(cOptSolution*);
  void     GetVar(int *);
  void     Mutate(double);
  void     ExploreBounds(void);
  void     ExploreUpperBounds(void);
  void     Hypermutation(double);
  void     Crossover(eCrossType, double, cOptSolution*, cOptSolution*);
  void     Differentiation(eDifType, double, cOptSolution*, cOptSolution*, cOptSolution*, cOptSolution*);


  void     Send(void);
  void     Receive(int);
};


// -------------------------------------------------------------------------
// Definition of cIndivIntMat class:
//
class cIndivIntMat : public cIndividual
{
 protected:
  int    **Var;   // Optimization variables
  
  cOptSolution*  NewObj(void) {return new cIndivIntMat(Prob);}

 public:
           cPopulation* individual;
           cOptSolution* BestSol(void);
           cIndivIntMat(cProblem *);
          ~cIndivIntMat(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     Init(const cVector&);
  void     AddBestSample(cVector, double);
  void     GetBestSample(cVector &);
  void     Evaluate(void);
  void     Print(void);
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  bool     CompVar(cOptSolution*);
  int **   GetVar(void) {return Var;}
  void     WriteGetVar(int **);
  void     Mutate(double);
  void     ExploreBounds(void);
  void     ExploreUpperBounds(void);
  void     Hypermutation(double);
  void     Crossover(eCrossType, double, cOptSolution*, cOptSolution*);
  void     Differentiation(eDifType, double, cOptSolution*, cOptSolution*, cOptSolution*, cOptSolution*);
  void     GetNormVar(cVector &xn);
  
  void     LamMutate(double*);
  void     Swap(double);
  void     Add(double);
  void     Delete(double);
  void     Send(void);
  void     Receive(int);
};

// -------------------------------------------------------------------------
// Definition of cIndivDblVec class:
//
class cIndivDblVec : public cIndividual
{
 protected:
  cVector  Var;   // Optimization variables (vector of doubles)
  
  cOptSolution*  NewObj(void) {return new cIndivDblVec(Prob);}

 public:
           cIndivDblVec(cProblem *);
          ~cIndivDblVec(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     Init(const cVector&);
  void     AddBestSample(cVector, double);
  void     GetBestSample(cVector &);
  void     Evaluate(void);
  void     Print(void);
  void     GetVar(cVector&);
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  bool     CompVar(cOptSolution*);
  void     GetNormVar(cVector&);
  cVector  GetVec(void)     {return Var;}
  cVector  GetBestVec(void) {return Var;}

  void     Mutate(double);
  void     ExploreBounds(void);
  void     ExploreUpperBounds(void);
  void     Hypermutation(double);
  void     Crossover(eCrossType, double, cOptSolution*, cOptSolution*);
  void     Differentiation(eDifType, double, cOptSolution*, cOptSolution*, cOptSolution*, cOptSolution*);
  
  void     Send(void);
  void     Receive(int);
};

// -------------------------------------------------------------------------
// Definition of cIndivBinary class:
//
class cIndivBinVec : public cIndividual
{
 protected:
  int     *Var;     // Optimization variables
  int     *VarBin;  // Binary optimization variables
  int     *Scromo;
  int      Sind;
  
  cOptSolution*  NewObj(void) {return new cIndivBinVec(Prob);}

 public:
           cIndivBinVec(cProblem *);
          ~cIndivBinVec(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     AddBestSample(cVector, double);
  void     GetBestSample(cVector &);
  void     Evaluate(void);
  void     Print(void) ;
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  bool     CompVar(cOptSolution*);
  void     Mutate(double);
  void     ExploreBounds(void);
  void     ExploreUpperBounds(void);
  void     Hypermutation(double) {;}
  void     Crossover(eCrossType, double, cOptSolution*, cOptSolution*);
  void     Differentiation(eDifType, double, cOptSolution*, cOptSolution*, cOptSolution*, cOptSolution*);
  
  void     Decoding(void);
  void     Send(void);
  void     Receive(int);
};

#endif
