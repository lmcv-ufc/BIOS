// -------------------------------------------------------------------------
// optsolution.h - file containing the definition of the cOptSolution class.
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
// The cOptSolution class defines the optimization solution described by the
// optimization agent of each method, as individual for genetic algorithms and
// particle for PSO. The generalized methods used for all agents are defined in
// the base class and the agent specified methods are defined in its subclasses.
//
// The class handles the storage of design variables, objective function and
// constraint vector. The design variables can be stored in
// different formats depending on the optimization algorithm and the problem
// to be solved. The optimization agents are generally created, stored and
// destroyed by the optimization algorithm.
//
// OptSolution
// |-- Individual
// |   |-- IndivIntVec
// |   |-- IndivIntMat
// |   |-- IndivDblVec
// |   |-- IndivBinVec
// |-- Particle
// |   |-- PartIntVec
// |   |-- PartIntMat
// |   |-- PartDblVec
// |-- Food
// |   |-- FoodIntVec
// |   |-- FoodIntMat
// |   |-- FoodDblVec
//
// -------------------------------------------------------------------------
// Virtual protected methods:
// -------------------------------------------------------------------------
//
// cOptSolution* NewObj(void)
//
// This method creates a object of same class and return it.
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// bool operator< (cOptSolution &sol)
//
//   sol - rhs solution                                            (in)
//
// This method compares the rhs solution with the lhs sol (this)
// based on the objective function (unconstrained problems) or the penalized
// objective function (constrained problems). It returns true if (f < sol.f)
// -------------------------------------------------------------------------
//
// cProblem* GetProb(void)
//
// This method returns a pointer to the optimization problem.
// -------------------------------------------------------------------------
//
// double GetCurrObjFunc(void)
//
// This method returns the current value of the objective function.
// -------------------------------------------------------------------------
//
// double GetCurrPenObjFunc(void)
//
// This method returns the current value of the penalized objective
// function. It is used in the solution of constrained problems.
// -------------------------------------------------------------------------
//
// double GetFitFunc(void)
//
// This method returns the current value of the fitness function.
// -------------------------------------------------------------------------
//
// const cVector& GetCurrConstr(void)
//
// This method returns the current values of the constraints.
// -------------------------------------------------------------------------
//
// void AssignObjFunc(double f)
//
//   f - value of the objective function                           (in)
//
// This method assigns the given value to the objective function of the
// solution.
// -------------------------------------------------------------------------
//
// void AssignPenObjFunc(double f)
//
//   f - value of the penalized objective function                 (in)
//
// This method assigns the given value to the penalized objective function
// of the solution. It is used in the solution of constrained problems.
// -------------------------------------------------------------------------
//
// void AssignFitFunc(double f)
//
//   f - value of the fitness function                             (in)
//
// This method assigns the given value to the fitness function of the
// solution. It is used by Genetic Algorithms.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// double GetObjFunc(void)
//
// This method is the same as GetCurrObjFunc, but in Particles subclass it
// return the history best objective function found.
// -------------------------------------------------------------------------
//
// double GetCurrPenObjFunc(void)
//
// This method is the same as GetCurrPenObjFunc, but in Particles subclass it
// return the history best penalized objective function found.
// -------------------------------------------------------------------------
//
// cVector& GetCurrConstr(void)
//
// This method is the same as GetCurrConstr, but in Particles subclass it
// return the history best constraints found.
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// cOptSolution* Clone(void)
//
// This methods creates a clone of this object and return it.
// -------------------------------------------------------------------------
//
// void Init(void)
//
// Uses the associated problem type to evaluate the initial values of the
// solution variables.
// -------------------------------------------------------------------------
//
// void Evaluate(void)
//
// This method uses the associated problem type to evaluate the objective
// function and the constraint vector of the solution.
// -------------------------------------------------------------------------
//
// void Print(void)
//
// Prints the solution data (objective function, constraint violations,
// design variable values) on the screen.
// -------------------------------------------------------------------------
//
// void Write(ostream &out)
//
//   out - output stream object                                  (in/out)
//
// Writes the solution data (objective function, constraint violations,
// design variable values) on the output file.
// -------------------------------------------------------------------------
//
// void Copy(cOptSolution* sol)
//
//   sol - pointer to the solution to be copied                  (in)
//
// This method copies the data of a given solution. Both solutions
// should be of same type, but need not to belong to the same group.
// -------------------------------------------------------------------------
//
// bool CompVar(cOptSolution* sol)
//
//   sol - pointer to the solution to be compared                (in)
//
// This method test if the variables of a given solution are equal to this
// object.
// -------------------------------------------------------------------------

#ifndef _OPTSOLUTION_H
#define _OPTSOLUTION_H

#include <vector>
#include <iosfwd>
using namespace std;

#include "vec.h"

// ------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cOptSolution;
class cIndividual;
class cParticle;
class cFood;
class cSampSAO;

// ------------------------------------------------------------------------
// Auxiliary types:
//
struct sOptSolCompare
{
  sOptSolCompare(bool t=false);
  bool operator() (cOptSolution*,cOptSolution*); // Comparison function
  bool HBest;    // History Best flag
}; 

// -------------------------------------------------------------------------
// Solution types:
//
typedef enum
{
  SOL_INT_VEC,
  SOL_INT_MAT,
  SOL_DBL_VEC,
  SOL_BIN_VEC,
} eSolType;

istream& operator>>(istream&,eSolType&);

// -------------------------------------------------------------------------
// Input Solutions:
//
typedef struct
{
  eSolType type;
  cVector  CodVar;  // Codified variables.
}  sInpSol;

// -------------------------------------------------------------------------
// Variable types:
//
typedef enum
{
  DISCRETE,
  CONTINUOUS
} eVarType;

// -------------------------------------------------------------------------
// Definition of cOptSolution class:
//
class cOptSolution
{
 protected:

          eSolType       OptSolType;    // Optimization Solution type
          cProblem*      Prob;          // Problem type
          double         PenObjFunc;    // Penalized objective function value
          double         NormConst;     // Normalized constraints value
          double         FeasibleConst; // Normalized nonviolated constraits values
          double         FitFuncVal;    // Fitness function value
          cVector        Constr;        // Constraint values
          cVector        Fobjs;         // Objective functions vector

  virtual cOptSolution*  NewObj(void) = 0;

 public:
  static  bool           Compare(cOptSolution *a,cOptSolution *b) {return *a < *b;}
  static  eVarType       GetVarType(const eSolType&);

                        cOptSolution(cProblem*);
  virtual              ~cOptSolution(void);
          bool          operator< (cOptSolution&);

          cProblem*     GetProb(void)      { return Prob; }
  virtual double        GetObjFunc(int id)   { return Fobjs[id]; }
          double        GetCurrObjFunc(int id) { return Fobjs[id]; }
  virtual double        GetPenObjFunc(void){ return PenObjFunc; }
          double        GetCurrPenObjFunc(void){ return PenObjFunc; }
          double        GetFitFunc(void)   { return FitFuncVal; }
          double        GetNormConst(void){ return NormConst; }
          double        GetFeasibleConst(void){ return FeasibleConst; }
  virtual int **        GetVar(void) { };
  virtual void          GetVar(cVector &) { };
  virtual cVector       GetVec(void) { };
  virtual cVector       GetBestVec(void) { };
  virtual void          GetVar(int*) { };

  virtual const cVector&  GetConstr(void) { return Constr; }
          const cVector&  GetCurrConstr(void) { return Constr; }
          void            AssignPenObjFunc(double f){ PenObjFunc = f; }
          void            AssignNormConst(double f){ NormConst = f; }
          void            AssignFeasibleConst(double g) {FeasibleConst = g;}
          void            AssignFitFunc(double f)   { FitFuncVal = f; }
  virtual void            GetNormVar(cVector&);
          cOptSolution*   Clone(void);
 
  virtual void          Init(void){ };
  virtual void          Init(int){ };
  virtual void          Init(const sInpSol&){ };
  virtual void          Init(const cVector&){ };
  virtual void          Evaluate(void) = 0;
  virtual void          Print(void) = 0;
  virtual void          Write(std::ostream&) = 0;
  virtual void          Copy(cOptSolution*) = 0;
  virtual bool          CompVar(cOptSolution*) = 0;

  // Laminate Operators.

  virtual void      LamMutate(double*);
  virtual void      Swap(double);
  virtual void      Add(double);
  virtual void      Delete(double);

  // These methods are important for MPI paralelization

  virtual void      Send(void) = 0;
  virtual void      Receive(int) = 0;
};

typedef std::vector<cOptSolution*> VecOptSol;

#endif
