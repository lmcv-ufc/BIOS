// -------------------------------------------------------------------------
// lamplt.h - file containing the definition of the Laminated Plate class.
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
// The class cLamPlate and its subclass implements analytical laminated
// composite plate problems.
// -------------------------------------------------------------------------

#ifndef _LAMPLT_H
#define _LAMPLT_H

#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <vector>
using namespace std;

#include "problem.h"
#include "lam.h"
#include "rbf.h"
#include "krg.h"

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;
class cMatrix;

// -------------------------------------------------------------------------
// Definition of cLamPlate class:
//
class cLamPlate : public cLaminated
{
 protected:
  static int            MaxSurGen;  // Max number of generations using SM
  static vector<double> Loads;      // Plate loads
  static vector<double> Dims;       // Plate dimensions 

 public:
          void    ReadLamMaxSurGen(std::istream&);
          void    ReadLamPltLoads(std::istream&);
          void    ReadLamPltDims(std::istream&);

                  cLamPlate(void);
  virtual        ~cLamPlate(void);
  virtual void    LoadReadFunc(cInpMap&);
  virtual void    FindNewSampPoint(cVector &, int, int);
  virtual double  EvalSampPoint(cVector &) { return 1.0; }
  virtual double  CalcBuckLoad(cMatrix &, double, double, double, double);          
};

// -------------------------------------------------------------------------
// Class for multi-objective weight and cost minimization of simply supported
// laminated plates with buckling constraint (Rao and Lakshmi, 2009):
//
class cLamPltBuckMOBJ: public cLamPlate
{
 public:

           cLamPltBuckMOBJ(void);
          ~cLamPltBuckMOBJ(void) { }
    void   Evaluate(int**, cVector &, cVector &);
};

// -------------------------------------------------------------------------

// Class for maximization of the load factor of simply supported laminated
// plates, considering Shear force - Liu et al. (2000):
//

class cLamPltLiu2000 : public cLamPlate
{
 public:
           cLamPltLiu2000(void);
          ~cLamPltLiu2000(void) { }

  void    Evaluate(int**, cVector &, cVector &);
};

// -------------------------------------------------------------------------
// Class for maximization of the load factor of simply supported laminated
// plates considering buckling and material failure - Kogiso et al. (1994):
//
class cLamPltLoadFactor : public cLamPlate
{
 protected:
  void     Analysis(cMatrix &, double &, double &);

 public:
           cLamPltLoadFactor(void);
          ~cLamPltLoadFactor(void) { }
  void    Evaluate(int**, cVector &, cVector &);
  void    EvalExactConstraint(int,int**,double &);
  void    GetApproxObj(bool*o) { o[0] = 1; }
  void    GetApproxConstr(bool*);
};

// -------------------------------------------------------------------------
// Class for cost minimization of simply supported laminated plates with
// failure criteria (Lopez 2009, ex 4.2):
//

class cLamPltMinLopez2009 : public cLamPlate
{
 public:
           cLamPltMinLopez2009(void);
          ~cLamPltMinLopez2009(void) { }
  void    Evaluate(int**, cVector &, cVector &);
};

#endif
