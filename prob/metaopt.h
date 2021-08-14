// -------------------------------------------------------------------------
// metaopt.h - file containing the definition cMetaOptimization class.
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
// The MetaOptimization class implements meta-optimization process as problem
// object. This special type of optimization problem is
// used to calibrate optimization variables for a set of problems. These
// problems are defined by cMetaProblem class. The Optimization method
// parameters are used as optimization variables in the meta-optimization
// problem. The variables can be defined as integer or continuous, and each
// subclass define a meta-optimization problem of speficied algorithm.
//
// MetaOptDiscrete
// |-- MetaGAD
// |-- MetaLamGAD
// |-- MetaPSOD
// |-- MetaLamPSOD
// |-- MetaLamPSOD
// |-- MetaABCD
// |-- MetaAISD
// MetaOptContinuous (*)
// |-- MetaGAC (*) 
// |-- MetaLamGAC (*) 
// |-- MetaPSOC (*)
// |-- MetaLamPSOC (*)
// |-- MetaABCC (*)
// |-- MetaAISC (*)
//
// (*) Not yet implemented.

// -------------------------------------------------------------------------

#ifndef _METAOPT_H
#define _METAOPT_H

#include "problem.h"
#include "stdpso.h"
#include <fstream>

using namespace std;

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cMetaProblem;

// ------------------------------------------------------------------------
// Definition of the discrete metaoptimization class:
//
class cMetaOptDiscrete : public cProblem
{
 protected:
  int          *ListDim;     // Number of discrete values for each variable
  cVector      *List;        // Lists of discrete values for each variable
  cMetaProblem *ProbSet;     // Set of problems to be metaoptimized
  cVector       DefaultVal;  // Default Variables values with off flag.

 public:
  int             EvalNumber;  // Number of evaluations requested.
  int             OptNumber;   // Number of opts to average mean values.
  vector<bool>    VarFlag;     // Flag to set algorithms variables on/off
  eSwaTopType     Topology;

  void            ReadProbOptNum(std::istream&); 
  void            ReadEvalNumber(std::istream&);
  void            ReadAlgFlag(std::istream&);
  void            ReadSwarmTopology(std::istream&);

                         cMetaOptDiscrete(void);
                         ~cMetaOptDiscrete(void);
         void            LoadReadFunc(cInpMap&);
          void           Decode(cVector&,cVector&);
          void           WriteControl(fstream&,cVector&);   
  virtual void           WriteAlgorithm(fstream&,cVector&) = 0;
  virtual string         AlgorithmSymbol(void) = 0;

  void     PrintVar(int*);
  void     WriteVar(int*,std::ostream&);
  void     GetBounds(int, int*, int*);
  void     Evaluate(int*,cVector&,cVector&);
};

//-------------------------------------------------------------------------
// Definition of the continuous metaoptimization class:
//
class cMetaOptContinuous : public cProblem
{
 protected:
  double *Low;
  double *Upp;

 public:
           cMetaOptContinuous(void);
          ~cMetaOptContinuous(void);
  void     PrintVar(double*);
  void     WriteVar(double*,std::ostream&);
  void     GetDblBounds(double*, double*);
};

// ------------------------------------------------------------------------
// Definition of MetaOptRSD class:
//
class cMetaOptRSD : public cMetaOptDiscrete
{
 public:
          cMetaOptRSD(void);
         ~cMetaOptRSD(void){ }
  void    WriteAlgorithm(fstream&,cVector&);
  string  AlgorithmSymbol(void) {return string("RS");}
};

// ------------------------------------------------------------------------
// Definition of MetaOptGAD class:
//
class cMetaOptGAD : public cMetaOptDiscrete
{
 public:
          cMetaOptGAD(void);
         ~cMetaOptGAD(void){ }
  void    WriteAlgorithm(fstream&,cVector&);
  string  AlgorithmSymbol(void) {return string("GA");}
};

// ------------------------------------------------------------------------
// Definition of MetaOptLamGAD class:
//
class cMetaOptLamGAD : public cMetaOptDiscrete
{
 public:
          cMetaOptLamGAD(void);
         ~cMetaOptLamGAD(void){ }
  void    WriteAlgorithm(fstream&,cVector&);
  string  AlgorithmSymbol(void) {return string("LamGA");}
};
// ------------------------------------------------------------------------
// Definition of MetaOptPSOD class:
//
class cMetaOptPSOD : public cMetaOptDiscrete
{
 public:
          cMetaOptPSOD(void);
         ~cMetaOptPSOD(void){ }
	 
  virtual void    WriteAlgorithm(fstream&,cVector&);
  virtual string  AlgorithmSymbol(void) {return string("PSO");}
};

// ------------------------------------------------------------------------
// Definition of MetaOptLamPSOD class:
//
class cMetaOptLamPSOD : public cMetaOptDiscrete
{
 public:
          cMetaOptLamPSOD(void);
         ~cMetaOptLamPSOD(void){ }
	 
  void    WriteAlgorithm(fstream&,cVector&);
  string  AlgorithmSymbol(void) {return string("LamPSO");}
};

// ------------------------------------------------------------------------
// Definition of MetaOptABCD class:
//
class cMetaOptABCD : public cMetaOptDiscrete
{
 public:
          cMetaOptABCD(void);
         ~cMetaOptABCD(void){ }
  void    WriteAlgorithm(fstream&,cVector&);
  string  AlgorithmSymbol(void) {return string("ABC");}
};

// ------------------------------------------------------------------------
// Definition of MetaOptAISD class:
//
class cMetaOptAISD : public cMetaOptDiscrete
{
 public:
          cMetaOptAISD(void);
         ~cMetaOptAISD(void){ }
  void    WriteAlgorithm(fstream&,cVector&);
  string  AlgorithmSymbol(void) {return string("AIS");}
};

#endif
