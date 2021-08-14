// -------------------------------------------------------------------------
// optalg.h - file containing the definition of the cOptAlgorithm class.
// -------------------------------------------------------------------------
// Copyright (c) 2013 LMCV/UFC
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
// The cOptAlgorithm class is responsible for the control of the optimization
// process. It provides both heuristic and bayesian optimization techniques.
//
// OptAlgorithm
// |-- StandardGA
// |---- LaminatedGA
// |-- NSGA
// |---- LaminatedNSGA
// |-- StandardPSO
// |---- LaminatedPSO
// |-- StandardDE
// |-- StandardABC
// |-- StandardAIS
// |-- RandomSearch
// |-- SAO
// |---- KRGSAO
// |---- RBFSAO
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadOptNum(void)
//  void ReadTolViol(void)
//  void ReadTolSucRate(void)
//  void ReadPopSize(void)
//  void ReadCrossMethod(void)
//  void ReadMaxGen(void)
//  void ReadStallGen(void)
//  void ReadMigrNum(void)
//  void ReadMigrGap(void)
//  void ReadMutRate(void)
//  void ReadMutRange(void)
//
// These methods read relevant data to the genetic algorithm and store them
// in protected variables.
// -------------------------------------------------------------------------
//
// void ReadAlg(void)
//
// This method reads the optimization algorithm from the input file.
// -------------------------------------------------------------------------
//
// void ReadMaxThread(void)
//
// This method reads the number of threads used on the program parallel
// version (OpenMP).
// -------------------------------------------------------------------------
//
// void End(void)
//
// This method finalizes the optimization process and prints the END label
// on the output file.
// -------------------------------------------------------------------------
//
// eOptAlgType GetInpType(void)
//
// This method returns the type of the algorithm read from the input file.
// -------------------------------------------------------------------------
//
// double GetTolViol(void)
//
// This method returns the tolerance for constraint violation.
// -------------------------------------------------------------------------
//
// cOptAlgorithm *CreateOptAlg(eOptAlgType)
//
//   type - algorithm type                                         (in)
//
// This method creates an optimization algorithm object according to the
// given type and returns a pointer to this object.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// void SetPenFunction(cPenalty *p)
//
//   p - pointer to penalty function                               (in)
//
// This method stores a pointer to the penalty function object to be used in
// the solution of contrained opmization problems.
// -------------------------------------------------------------------------
//
// void SetSelMethod(cSelection *s)
//
//   s - pointer to selection method                               (in)
//
// This method stores a pointer to the selection method to be used in
// Genetic Algorithm.
// -------------------------------------------------------------------------
//
// eOptAlgType GetType(void)
//
// This method returns the algorithm type.
// -------------------------------------------------------------------------
//
// cProblem *GetProblem(void)
//
// This method returns a pointer to the optimization problem being solved.
// -------------------------------------------------------------------------
//
// void Init(cProblem*)
//
// This method initialize some dynamic variables before start the solving
// operation.
// -------------------------------------------------------------------------
//
// void UpdatePostVar(int gen, int opt, double &lb, cGroup *mg)
//
//   gen - current generations number                              (in)
//   opt - current optimization number                             (in)
//   lb  - last best objective function                            (in/out)
//   mg  - solution group contaning the best solution              (in)
//
//
// This method update variables related to prost-processing. (e.g Mean best
// solutions over generations).
// -------------------------------------------------------------------------
//
// bool OptStopCrit(int gen, int opt, double &lb, cGroup *mg)
//
//   gen - current generations number                              (in)
//   opt - current optimization number                             (in)
//   lb  - last best objective function                            (in)
//   mg  - solution group contaning the best solution              (in)
//
// This method verify if conditions to stop optimization are satisfied. In the
// case of triggering, the reminder postprocessing data will be filled with the
// last available results, and return a true value. Return false otherwise.
// -------------------------------------------------------------------------
//
// void PrintPostVar(int maxgen, int opt, double evnum)
//
//   maxgen - maximum number of generations                        (in)
//   opt    - current optimization number                          (in)
//   evnum  - number of evaluations                                (in)
//
//
// This method print the output data of current optimization run on output
// file.
// -------------------------------------------------------------------------
//
// void PostProcessing(void)
//
// This method evaluate the postprocessing data related to all optimizations
// runs and print it in the output file.
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// virtual void Solver(void)
//
// This is the main method of the cOptAlgorithm class, where the necessary
// steps for solving the optimization problem are undertaken.
// -------------------------------------------------------------------------

#ifndef _OPTALGORITHM_H
#define _OPTALGORITHM_H

#include "input.h"
#include "mat.h"
#include "optsolution.h"
#include "samp.h"
#include <vector>
#include <iosfwd>

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cProblem;
class cPenalty;
class cSelection;
class cGroup;
class cSolGroup;
class cOptAlgorithm;

// -------------------------------------------------------------------------
// Control Algorithms:
//
typedef enum
{
  STANDARD_GA,            // Standard Genetic Algorithm
  LAMINATE_GA,	          // Genetic Algorithm with laminate operators
  modNSGAII,              // Nondominated Sorting Genetic Algorithm II
  modLAMINATE_NSGAII,     // Nondominated Sorting Genetic Algorithm with laminate operators
  STANDARD_PSO,           // Standard Particle Swarm Algorithm
  LAMINATE_PSO,           // Particle Swarm Algorithm with laminate operators
  STANDARD_ABC,           // Standard Artificial Bee Colony.
  STANDARD_AIS,           // Standard Artificial Immuny System.
  STANDARD_DE,            // Standard Differential Evolution.
  SAORBF,                 // Sequential Approximate Optimization using RBF models
  SAOKRG,                 // Sequential Approximate Optimization using KRG models
  RANDOM_SEARCH           // Random search.
} eOptAlgType;

// -------------------------------------------------------------------------
// Crossover Methods:
//
typedef enum
{
  LINEAR_COMBINATION,
  CLASSICAL,
  BIN
} eCrossType;

// -------------------------------------------------------------------------
// Differentiation Methods:
//

typedef enum
{
  Rand1,
  Loc2Best,
  BestJitter
} eDifType;

// -------------------------------------------------------------------------
// Swarm Topology:
//
typedef enum
{
  GBEST_TOPOLOGY,       // gbest population topology. [1,2]
  RING_TOPOLOGY,        // lbest (ring) population topology.
  SQUARE_TOPOLOGY       // square (von Neumann) population topology.
} eSwaTopType;

// -------------------------------------------------------------------------
// Crossover Methods:
//
class cOptAlgReadEntry : public cAbsReadEntry
{
 public:
  cOptAlgorithm *alg;
  cInpMap       *inpmap;

  void Read(std::istream&);
};

// -------------------------------------------------------------------------
// Definition of cOptAlgorithm class:
//
class cOptAlgorithm
{
 protected:
  cInpMap       *InpMap;        // Input map used to read input file 
  eCrossType     CrossType;     // Type of crossover read from input file
  eDifType       DifType;       // Type of differentiation read from input file (DE)
  int            OptNum;        // Number of optimizations
  int            PopSize;       // Population size
  int            cont;          // Population indice
  int            StallGen;      // Generations without improvement (stop. crit.)
  int            MaxGen;        // Maximum number of generations
  int            MigrationNum;  // Number of individuals to be migrated
  int            MigrationGap;  // Interval between migrations
  double         TolViol;       // Tolerance for constraint violation
  double         TolSucRate;    // Tolerance for success optimization
  double         MutProb;       // Mutation probability
  double         MaxMut;        // Maximum mutation probability
  double         MinMut;        // Minimum mutation probability
  int            NumBestSol;    // Number of best write on output file
  eSolType       SolType;       // Optimitzation Solution type read from the input file
  sInpSol*       InpSolVec;     // Vector of input solutions.
  int            NumInpSol;     // Number of input solutions.
  eSampType      SampType;      // Input sample type.
  bool           IntPopSamp;    // Flag to control sample code usage.
  bool           Feedback;      // Flag to print console feedback.
  std::ostream  *out;           // Output stream object. 

  double         MinObjFunc; // Minimun objective function
  bool           MinObjFlag; // Check if the minimum obj. func was provided

  vector<cVector>   NRMSE;         // Normalized Root Mean Squared Error
  vector<cVector>   RMAE;          // Relative Mean Average Error
  vector<cVector>   NMAE;          // Normalized Maximum Absolute Error
  vector<cVector>   ERROR;         // Error
  vector<cVector>   Thetas;
  vector<cVector>   EImax;         // Maximum EI of iterations
  vector<cVector>   Ybest;         // Best design
  vector<cVector>   Cbest;         // Best design constraints


          eOptAlgType    Type;          // Type of optimization algorithm
          cProblem      *Prob;          // Pointer to optimization problem
          cPenalty      *Pen;           // Pointer to penalty method
          cSelection    *Sel;           // Pointer to selection method

          int            StallGenCount; // Counter to StallGen
          int           *GentoConv;     // Number of genenations to converged
          double        *GenBest;       // Generations best solutions
          cMatrix        GenBestConstraints;  // Constraints related to the best individual throughout the generations
          cMatrix        GenVarBest;          // Best individual throughout the generations
          double        *GenAvg;        // Generations Average solutions
          double        *MBestGen;      // Generations mean best over optimizations
          cSolGroup     *best;          // Group of best solutions

 public:
          void           SetPenFunction(cPenalty* p)   { Pen          = p;    }
          void           SetSelMethod(cSelection* s)   { Sel          = s;    }
          void           SetOptNum(int n)              { OptNum       = n;    }
          void           SetPopSize(int s)             { PopSize      = s;    }
          void           SetStallGen(int sg)           { StallGen     = sg;   }
          void           SetMaxGen(int mg)             { MaxGen       = mg;   }
          void           SetMigrationNum(int mn)       { MigrationNum = mn;   }
          void           SetMigrationGap(int mg)       { MigrationGap = mg;   }
          void           SetTolViol(double tv)         { TolViol      = tv;   }
          void           SetTolSucRate(double tsr)     { TolSucRate   = tsr;  }
          void           SetMutProb(double mp)         { MutProb      = mp;   }
          void           SetMaxMut(double mm)          { MaxMut       = mm;   }
          void           SetMinMut(double mm)          { MinMut       = mm;   }
          void           SetProblem(cProblem* prob)    { Prob         = prob; }
          void           SetSolType(eSolType sol)      { SolType      = sol;  }
          void           SetOutStream(std::ostream &o) { out          = &o;   }
          void           SetFeedback(bool fb)          { Feedback     = fb;   }


          virtual  void  SetSwarmTopology(eSwaTopType){;}
          virtual  void  SetDifType(eDifType){;}


          void           ReadAlg(std::istream&);
          void           ReadOptNum(std::istream&);
          void           ReadTolViol(std::istream&);
          void           ReadTolSucRate(std::istream&);
          void           ReadPopSize(std::istream&);
          void           ReadSampType(std::istream&);
          void           ReadCrossMethod(std::istream&);
          void           ReadDifType(std::istream&);
          void           ReadMaxGen(std::istream&);
          void           ReadStallGen(std::istream&);
          void           ReadMigrNum(std::istream&);
          void           ReadMigrGap(std::istream&);
          void           ReadMutProb(std::istream&);
          void           ReadMutRange(std::istream&);
          void           ReadMaxThread(std::istream&);
          void           ReadSeed(std::istream&);
          void           ReadNumBestSol(std::istream&);
          void           ReadSelMethod(std::istream&);
          void           ReadPenalty(std::istream&);
          void           ReadOptSolType(std::istream&);
          void           ReadNumInpSol(std::istream&);
          void           ReadInpSol(std::istream&);
          void           ReadProblem(std::istream&);
          void           ReadMinObjFunc(std::istream&);
  virtual void           LoadReadFunc(cInpMap&);

          void           End(void);
          int            GetMaxGen(void)  { return MaxGen;  }
          int            GetPopSize(void) { return PopSize; }
          double         GetTolViol(void) { return TolViol; }
  static  cOptAlgorithm* CreateOptAlg(eOptAlgType, cProblem*);

                         cOptAlgorithm(void);
  virtual               ~cOptAlgorithm(void);
          eOptAlgType    GetType(void) { return Type; }
          cProblem*      GetProblem(void) { return Prob; }
          cOptSolution*  GetBest(void);
          void           Init(void);
          void           GetBestIndividuals(int,cVector&,cMatrix&,cGroup*);
          void           UpdatePostVar(int,int,double&,cGroup*);
          bool           OptStopCrit(int,int,double&,cGroup*);
          void           PrintPostVar(int,int,double,cGroup*);
          void           PrintPostVar(int,int,double,cGroup*, std::vector<int> *solrank);
          void           PostProcessing(void);

  virtual void           Solver(void) = 0;
};

#endif
