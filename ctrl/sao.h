// -------------------------------------------------------------------------
// sao.h - file containing the definition of the cSAO class.
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
// The cSAO class implements the general method for the solution of
// Sequential Approximate Optimization problems. It starts out by sampling
// in the problem design space, then fitting a metamodel and, finally, a
// new point is added in each iteration. Thus, the cSAO class provides
// important implementations such as:
//   - Solver( )
//   - SetInitialSample( )
//   - CreateSurrogate( )
//   - EvalInfillCriteria( )
//   - UpdateSurrogate( )
//
// -------------------------------------------------------------------------

#ifndef _SAO_H
#define _SAO_H

#include "optalg.h"
#include "penalty.h"
#include "group.h"
#include "surr.h"
#include "stdpso.h"
#include "samp.h"
#include "stdde.h"
#include "optsolution.h"
#include "probsurr.h"

// -------------------------------------------------------------------------
// Forward declarations:
//
class cProblem;
class cSelection;
class cIndividual;
class cSampSet;

// -------------------------------------------------------------------------
// Struct :
//
typedef struct
{  
  sInpSol var;
  cVector fvec;    // Vector of objective functions.
  cVector cvec;    // Vector of constraints. 
} sSAOSample;

struct sProbAppOut
{
 protected:
  bool *ApproxOF;   // Approximated objective function flag.
  bool *ApproxC;    // Approximated constraints flag.
  int   numobjf;
  int   numconst;
  int   numappconst;
  int   numappobjf;

 public:
  sProbAppOut(void);
  sProbAppOut(cProblem*);
 ~sProbAppOut(void);

  int    GetNumAppObjFunc( )  const { return numappobjf;             }
  int    GetNumAppConstr( )   const { return numappconst;            } 
  int    GetNumExactConstr( ) const { return numconst-numappconst;   } 
  int    GetNumAppOut( )      const { return numappobjf+numappconst; }

  int GetAppConstrID(const int&)    const;
  int GetExactConstrID(const int&)  const;
  int GetAppObjFuncID(const int&)   const;
  int GetAppConstrOutID(const int&) const;
};

// -------------------------------------------------------------------------
// Definition of the cSAO class:
//
class cSAO : public cOptAlgorithm
{
 protected:

    int             SubPop;
    int             SubMaxGen;
    int             GenStall;
    double          SubTolViol;
    double          SubMutProb;
    double          Nmax;
    cOptAlgorithm   *SubAlgType;
    int             NumInitSP;
    bool            InputNS;
    double          MinNRMSE;
    int             NumApproxObj;
    bool*           ApproxObj;
    int             NumApproxC;
    bool*           ApproxC;
    string          SampleFileName;
    eSwaTopType     SubTopology;
    eDifType        SubDifType;
    eConstrType     ConstrMethod; // LEO
    eProbSurrState  InfillCriteria;

    bool            FlagVS;
    int             Nvs;
    vector<cVector> Vsx;
    vector<cVector> Vsy;

    // Expect improvement data?
    double     WEI;
    double     Beta;
    bool       ciclewei;
    eSigmaType SigType;

 public:
                     cSAO(void);
  virtual           ~cSAO(void) {};

  virtual void      Solver(void);


          void      MaxExpectImprov(sProbAppOut&,cSURR*,cVector&,double&);

          void      MaxExpectImprov(cSURR*,cSampSet*,cVector&,double&);

          void      EvalInfillCriteria(cSURR*,cSampSet*,sProbAppOut&,cVector&,double&);

          double    GetBestFeasibleSample( );

  virtual cSURR*    CreateSurrogate(sSampData&) = 0;
  virtual void      UpdateSurrogate(cVectorVec&,cVectorVec&) = 0;

  virtual void      LoadReadFunc(cInpMap&);
  void              ReadSubPop(std::istream&);
  void              ReadSubMaxGen(std::istream&);
  void              ReadSubPSOTopology(std::istream&);
  void              ReadSubDEType(std::istream&);
  void              ReadSubTolViol(std::istream&);
  void              ReadSubMutProb(std::istream&);
  void              ReadSampleFileName(std::istream&);
  void              ReadNmax(std::istream&);
  void              ReadMinNRMSE(std::istream&);
  void              ReadValSamples(std::istream&);
  void              ReadSubAlgType(std::istream&);
  void              ReadNumInitSamplingPoints(std::istream&);
  void              ReadConstrMethod(std::istream&); // LEO
  void              ReadInfillCriteria(std::istream&); // LEO
  void              ReadWEI(std::istream&); // LEO
  void              ReadBeta(std::istream&); // LEO
  void              InitSample(sSampData&);
  void              SetApproxObj( );
  void              SetApproxConstr( );
  void              SetInitialSample(sProbAppOut&,cSampSet*&,sSampData&,int&);

  bool              OptStopCrit(int,int,double&,cGroup*);

  void              EvalErrorMeasures(vector<cVector> &, cVector &, cVector &, cVector &, cVector &, int);
  void              ReadSampleFile(ifstream&,int&,int&,int&,cVector&,cVector&,vector<cVector>&,vector<cVector>&, vector<cVector> &, vector<cVector> &);
};

#endif
