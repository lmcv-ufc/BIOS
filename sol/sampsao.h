// -------------------------------------------------------------------------
// sampsao.h - file containing the definition of the cSampSAO class.
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
// The cSampSAO defines each design from the surrogate data set. Each sample
// must have:
//   - Vector of design variables
//   - Vector of responses (Objective + Constraint functions)
// To provide an easy way of defining the best feasible designs, a penalized
// objective function is also computed.
//
// cSampSAO
// |-- VecInt
// |-- MatInt
// |-- VecDouble
//
// -------------------------------------------------------------------------

#ifndef _SAMPSAO_H
#define _SAMPSAO_H

#include <vector>
using namespace std;

#include "vec.h"
#include "optsolution.h"
#include "stdpso.h"

// -------------------------------------------------------------------------
// Definition of cParticle class:
//
class cSampSAO : public cOptSolution
{
 protected:
 sProbAppOut *appout;       // Approximated output data of opt. problem.

 public:
  static void        ReadIndType(void);
  static cSampSAO*   CreateSamp(eSolType,cProblem*,sProbAppOut&);

                     cSampSAO(cProblem*);
  virtual            ~cSampSAO(void);

  virtual void       GetSurrOutRes(cVector&) = 0;
};

// -------------------------------------------------------------------------
// Definition of cSampIntVec class:
//

class cSampIntVec : public cSampSAO
{
 protected:
  int     *Var;       // Current optimization variables

  cOptSolution*  NewObj(void) {return new cSampIntVec(Prob);}
 
 public:
                cSampIntVec(cProblem *);
               ~cSampIntVec(void);
  void          Init(void);
  void          Init(const sInpSol&);
  void          Init(const cVector&);
  void          Evaluate(void);
  void          Print(void);
  void          Write(std::ostream&);
  void          Copy(cOptSolution*);
  void          GetVar(int *);
  bool          CompVar(cOptSolution*);

  void          GetNormVar(cVector&);
  void          GetSurrOutRes(cVector&);

  void     Send(void) { }
  void     Receive(int) { }
};


// -------------------------------------------------------------------------
// Definition of cIndivIntMat class:
//

class cSampIntMat : public cSampSAO
{
 protected:
  int     **Var;       // Current optimization variables
  
  cOptSolution*  NewObj(void) {return new cSampIntMat(Prob);}

 public:
                cSampIntMat(cProblem *);
               ~cSampIntMat(void);
  void          Init(void);
  void          Init(const sInpSol&);
  void          Init(const cVector&);
  void          Evaluate(void);
  void          Print(void);
  void          Write(std::ostream&);
  void          Copy(cOptSolution*);
  bool          CompVar(cOptSolution*);
  int **        GetVar(void) {return Var;}

  void          GetNormVar(cVector&);
  void          GetSurrOutRes(cVector&);

  void          Send(void) { }
  void          Receive(int) { }
};


// -------------------------------------------------------------------------
// Definition of cPartDblVec class:
//
class cSampDblVec : public cSampSAO
{
 protected:
  cVector  Var;          // Sample optimization variables (vector of doubles)
  
  cOptSolution*  NewObj(void) {return new cSampDblVec(Prob);}
 
 public:
                cSampDblVec(cProblem *);
               ~cSampDblVec(void);
  void          Init(const cVector&);
  void          Evaluate(void);
  void          Print(void);
  void          Write(std::ostream&);
  void          Copy(cOptSolution*);
  bool          CompVar(cOptSolution*);
  cVector       GetVec(void) {return Var;}

  void          GetNormVar(cVector&);
  void          GetSurrOutRes(cVector&);

  void          Mutate(double);
                
  void          Send(void) { }
  void          Receive(int) { }
};

#endif
