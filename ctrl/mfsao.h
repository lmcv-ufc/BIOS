// -------------------------------------------------------------------------
// mfsao.h - file containing the definition of the cMFSAO class.
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
// The cEGO class implements the standard Genetic Algorithm with
// selection, crossover and mutation operators. Constrained optimization
// problems are solved using the penalty function approach.
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadCrossRate(void)
//  void ReadCrossRange(void)
//
// These methods read relevant data to the genetic algorithm and store them
// in protected variables.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// void Fitness(cGroup &pop)
//
//   pop - given population                                        (in/out)
//
// This method evaluates the fitness function value for all individuals
// of the given population.
// -------------------------------------------------------------------------
//
// void Mutation(cGroup &pop)
//
//   pop - given population                                        (in/out)
//
// This method applies the mutation genetic operator to all individuals
// of a given population. The mutation opertator is implemented by the
// cIndividual class.
// -------------------------------------------------------------------------
//
// void Merge(cGroup &pop, cGroup &son)
//
//   pop     - original population / merged Population             (in)
//   son     - offspring population                                (out)
//
// This method exchanges the worst individuals from a given original popu-
// lation with individuals of an offspring population. The best individuals
// from the original population are kept in the resulting population
// (elitism).
// -------------------------------------------------------------------------
//
// void Crossover(cGroup &parent, cGroup &son)
//
//   parent  - parent population                                   (in)
//
// This method performs the crossover operation generating a offspring
// population. The crossover operator is implemented by the  cIndividual
// class.
// -------------------------------------------------------------------------
//
// void RandomRates(void)
//
// This method randomizes the genetic operator rates if ranges were
// specified instead of fixed values by the user.
// -------------------------------------------------------------------------

#ifndef _MFSAO_H
#define _MFSAO_H

#include "optalg.h"
#include "penalty.h"
#include "group.h"
#include "surr.h"
#include "stdpso.h"
#include "samp.h"
#include "stdde.h"
#include "optsolution.h"
#include "probsurr.h" // LEO

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
} sMFSAOSample;

// -------------------------------------------------------------------------
// Crossover Methods:
//
/*
typedef enum
{
  LINEAR_COMBINATION,
  CLASSICAL
} eCrossType;*/

// -------------------------------------------------------------------------
// Definition of the cSAO class:
//
class cMFSAO : public cOptAlgorithm
{
 protected:

    int             SubPop;
    int             SubMaxGen;
    int             SubStallGen;
    int             GenStall;
    double          SubTolViol;
    double          SubMutProb;
    double          Nmax;
    cOptAlgorithm   *SubAlgType;
    int             NumInitSP;
    bool            InputNS;
    int             NumInitSPLF;
    bool            InputNSLF;
    bool            NestSamp;
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

    vector<int>     NumberHFE;
    vector<int>     NumberLFE;

    bool            FlagVS;
    int             Nvs;
    vector<cVector> Vsx;
    vector<cVector> Vsy;

    // Expect improvement data?
    double     WEI;
    double     Beta;
    double     NFac;
    bool       ciclewei;
    eSigmaType SigType;

    int CicleSize;
    cVector ListWEI;
    cVector ListBeta;

 public:
                     cMFSAO(void);
  virtual           ~cMFSAO(void) {};

  virtual void      Solver(void);


          void      MaxExpectImprov(sProbAppOut&,cSURR*,cVector&,double&);

          void      MaxExpectImprov(cSURR*,cSampSet*,cVector&,double&);

          void      EvalInfillCriteria(cSURR*,cSampSet*,sProbAppOut&,cVector&,double&);

          double    GetBestFeasibleSample( );

  virtual cSURR*    CreateSurrogate(sSampData&) = 0;
  virtual void      UpdateSurrogate(cVectorVec&,cVectorVec&,cVectorVec&) = 0;
  virtual void      GetFidelityNewPoint(cVectorVec&, bool&, bool&, double) = 0;

  virtual void      LoadReadFunc(cInpMap&);
  void              ReadSubPop(std::istream&);
  void              ReadSubMaxGen(std::istream&);
  void              ReadSubStallGen(std::istream&);
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
  void              ReadNumInitLFSamplingPoints(std::istream&);
  void              ReadUseNestedSample(std::istream&);
  void              ReadConstrMethod(std::istream&); // LEO
  void              ReadInfillCriteria(std::istream&); // LEO
  void              ReadCicleWEI(std::istream&); // LEO
  void              ReadWEI(std::istream&); // LEO
  void              ReadBeta(std::istream&); // LEO
  void              ReadCyclicWEI(std::istream&); // LEO
  void              ReadCyclicBeta(std::istream&); // LEO
  void              ReadNFacSohst(std::istream&); // LEO
  void              InitSample(sSampData&);
  void              SetApproxObj( );
  void              SetApproxConstr( );
  void              SetInitialSample(sProbAppOut&,cSampSet*&,sSampData&,int&);
  void              PrintHyperPar(cVector,int,int);

  void              PrintPostVarSur(int,int,double,cSURR*);

  bool              OptStopCrit(int,int,double&,cGroup*);

  void              PenLFSamples(cGroup *);

  void              EvalErrorMeasures(vector<cVector> &, cVector &, cVector &, cVector &, cVector &, int);
  void              ReadSampleFile(ifstream&,int&,int&,int&,cVector&,cVector&,vector<cVector>&,vector<cVector>&, vector<cVector> &, vector<cVector> &);
  void              PostProcessingSur(int, int, std::vector<cVector> *nrmse = 0, std::vector<cVector> *rmae =0,
                                   std::vector<cVector> *nmae = 0, std::vector<cVector> *error = 0,
                                   std::vector<int> *gentoconv = 0);
};

#endif
