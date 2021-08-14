// -------------------------------------------------------------------------
// food.h - file containing the definition of the cIndividual class.
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
// The cFood class define the optimization agent used in Artificial Bee Colony
// method. The ABC mutation is the main operator used by the method. A counter
// of failed attempts to find a better solution is considered and used by the
// ABC optimization method. The ABC mutate is defined in subclasses.
//
// Food
// |-- VecInt
// |-- MatInt
// |-- VecDouble
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// void GetCounter(void)
//
// This method returns the counter of failed attempts solutions. 
//
// -------------------------------------------------------------------------
//
// void IncrementCounter(void)
//
// This method increment the counter of failed attempts solutions. 
//
// -------------------------------------------------------------------------
//
// void ResetCounter(void)
//
// This method set to zero the counter of failed attempts solutions. 
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// void ABCmutate(cFoodSource &foods, int FoodIndex);
//
//   foods     - The food source used in ABC mutate                (in)
//   FoodIndex - This food index in the food source                (in)
//
// This method apply ABC mutate on this food, doing a linear combination with
// this solution and a random solution in the foods.
// -------------------------------------------------------------------------

#ifndef _FOOD_H
#define _FOOD_H

#include <vector>
using namespace std;

#include "optsolution.h"
#include "vec.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cFoodSource;

// -------------------------------------------------------------------------
// Definition of cFood class:
//
class cFood : public cOptSolution
{
 protected:
         int FailCounter;     // failed Attempts to improvement
 
 public:
  static cFood*     CreateFood(eSolType,cProblem*);

                    cFood(cProblem*);
  virtual          ~cFood(void);
  
          int       GetCounter(void) {return FailCounter;}
	  void      IncrementCounter(void) {FailCounter++;}
	  void      ResetCounter(void) {FailCounter = 0;}
  virtual void      ABCmutate(cFoodSource&,int) = 0;
};


// -------------------------------------------------------------------------
// Definition of cFoodIntVec class:
//
class cFoodIntVec : public cFood
{
 protected:
  int     *Var;   // Optimization variables
  
  cOptSolution*  NewObj(void) {return new cFoodIntVec(Prob);}

 public:
           cFoodIntVec(cProblem *);
          ~cFoodIntVec(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     Evaluate(void);
  void     Print(void);
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  bool     CompVar(cOptSolution*);

  void     ABCmutate(cFoodSource&,int);

  void     Send(void);
  void     Receive(int);
};

// -------------------------------------------------------------------------
// Definition of cFoodIntMat class:
//
class cFoodIntMat : public cFood
{
 protected:
  int    **Var;   // Optimization variables
  
  cOptSolution*  NewObj(void) {return new cFoodIntMat(Prob);}


 public:
           cFoodIntMat(cProblem *);
          ~cFoodIntMat(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     Evaluate(void);
  void     Print(void);
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  bool     CompVar(cOptSolution*);
  
  void     ABCmutate(cFoodSource&,int);
  
  void     Swap(double);
  void     Add(double);
  void     Delete(double);
  
  void     Send(void);
  void     Receive(int);
};

// -------------------------------------------------------------------------
// Definition of cFoodDblVec class:
//
class cFoodDblVec : public cFood
{
 protected:
  cVector  Var;   // Optimization variables (vector of doubles)
  
  cOptSolution*  NewObj(void) {return new cFoodDblVec(Prob);}

 public:
           cFoodDblVec(cProblem *);
          ~cFoodDblVec(void);
  void     Init(void);
  void     Init(const sInpSol&);
  void     Evaluate(void);
  void     Print(void);
  void     Write(std::ostream&);
  void     Copy(cOptSolution *);
  bool     CompVar(cOptSolution*);
  
  void     ABCmutate(cFoodSource&,int);
  
  void     Send(void);
  void     Receive(int);
};

#endif
