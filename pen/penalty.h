// -------------------------------------------------------------------------
// penalty.h - file containing the definition of the cPenalty class.
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
// The class cPenalty provides the basic template for the implementation of
// different penalty functions methods, both adaptive and user-guided.
//
// Penalty
// |-- Static [1]
// |-- Deb [2]
// |-- Adaptive [3]
// |-- Constraints normalization [4]
//
// References:
// [1] Michalewicz and Schoenauer (1996)
// [2] Deb (2000)
// [3] Lemonge and Barbosa (2004)
// [4] Deb et al. (2002)
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
// void ReadPenalty(void)
//
// This method reads the type of the penalty function.
// -------------------------------------------------------------------------
//
// ePenType GetInpType(void)
//
// This method returns the type of the penalty function read from the input
// file.
// -------------------------------------------------------------------------
//
// cPenalty *CreatePenalty(ePenType type)
//
//   type - penalty type                                           (in)
//
// This method creates a penalty function object based on the given type 
// and returns a pointer to this object.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// void EvalPenObjFunc(cGroup *group, double tol)
//
//   group - group of solutions                                    (in/out)
//   tol   - tolerance for constraint violation                    (in)
//
// Evaluates the penalized objective function for the given population based
// on the chosen penalty method.
// -------------------------------------------------------------------------

#ifndef _PENALTY_H
#define _PENALTY_H

// -------------------------------------------------------------------------
// Forward declarations:
//
class cGroup;
class cIndividual;
class cInpMap;

// -------------------------------------------------------------------------
// Penalty Methods:
//
typedef enum 
{
  STATIC,
  DEB2000,
  ADAPTIVE,
  NORMALIZATION,
  FEASIBLECONSTRAINTS
} ePenType;

// -------------------------------------------------------------------------
// Definition of cPenalty class:
//
class cPenalty
{
 protected:
  static ePenType   InpType;   // Penalty type read from input file 
         ePenType   Type;      // Penalty type

 public:
  static  void      ReadPenalty(void);
  static  ePenType  GetInpType(void) { return InpType; }
  static  cPenalty *CreatePenalty(ePenType);

                    cPenalty(void);
  virtual          ~cPenalty(void);
          ePenType  GetType (void) { return Type; }
  virtual void      EvalPenObjFunc(cGroup *, double) = 0;
  virtual void      LoadReadFunc(cInpMap&) { }
         
};

// -------------------------------------------------------------------------
// Definition of cPenStatic class:
//
class cPenStatic : public cPenalty
{
 private:
          double K;  // Constant penalty parameter

 public:
          void  ReadFactor(std::istream&);
          void  SetFactor(double pf)  { K  = pf; }

                cPenStatic(void);
  virtual      ~cPenStatic(void);
          void  LoadReadFunc(cInpMap&);
          void  EvalPenObjFunc(cGroup *, double);
};

// -------------------------------------------------------------------------
// Definition of cPenDeb class:
//
class cPenDeb : public cPenalty
{
 public:
                cPenDeb(void);
  virtual      ~cPenDeb(void);
          void  EvalPenObjFunc(cGroup *, double);
};

// -------------------------------------------------------------------------
// Definition of cPenAdapt class:
//
class cPenAdaptive : public cPenalty
{
 public:
                cPenAdaptive(void);
  virtual      ~cPenAdaptive(void);
          void  EvalPenObjFunc(cGroup *, double);
};

// -------------------------------------------------------------------------
// Definition of cNormalization class:
//
class cNormalization : public cPenalty
{
 public:
                cNormalization(void);
  virtual      ~cNormalization(void);
          void  EvalPenObjFunc(cGroup *, double);
};

// -------------------------------------------------------------------------
// Definition of cFeasibleConstraints class:
//
class cFeasibleConstraints : public cPenalty
{
 public:
                cFeasibleConstraints(void);
  virtual      ~cFeasibleConstraints(void);
          void  EvalPenObjFunc(cGroup *, double);
};

#endif
