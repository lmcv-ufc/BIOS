// -------------------------------------------------------------------------
// problem.h - file containing the definition of the cProblem class.
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
// The class cProblem contains data and methods relevant to the optimization 
// problem, which can be defined as:
//
//    Find x that minimize f(x)
//      
//      subject to cj(x) <= 0,  j = 1..nconstr
//     
// The optimization algorithms handles both constrained (nconstr > 0) and
// unconstrained problems (nconstr = 0). The constraint functions can assume
// real (double) values.
//
// The variables x may assume integer or real (double) values. The
// problem indicates the type of variables to the optimization algorithm.
//
// Real variables should be used to solve continuous problems.
//
// Integer variables should be used to solve integer or discrete problems.
// In the latter case the variable represents the index of an array of
// possible values (e.g. a list of structural profiles or list with
// thicknesses of steel plates). If necessary the decodification of the
// variables should be performed by the problem.
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
// void ReadProblem(void)
//
// This method reads the problem type.
// -------------------------------------------------------------------------
//
// eProbType GetInpType(void)
//
// This method returns the type of the problem read from the input file.
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// eProbType GetType(void)
//
// This method returns the problem type.
// -------------------------------------------------------------------------
//
// int GetNumVar(void)
//
// This method returns the number of variables of the problem.
// -------------------------------------------------------------------------
//
// int GetNumConstr(void)
//
// This method returns the number of constraints of the problem.
//
// -------------------------------------------------------------------------
// Virtual methods:
// -------------------------------------------------------------------------
//
// void PrintVar(int *var)
//
//   var - algorithm variables                                     (in)
//
// This method prints the variables in a custom format. This method should
// be used when the problem variables are obtained from the optimization
// variable by a decoding process.
// -------------------------------------------------------------------------
//
// void WriteVar(int *var)
//
//   var - algorithm variables                                     (in)
//
// This method writes the variables in a custom format to the output file. 
// This method should be used when the problem variables are obtained from 
// the optimization variable by a decoding process.
// -------------------------------------------------------------------------
//
// void GetBounds(int i, int* low, int* upp)
//
//   i   - index of the variable                                   (in)
//   low - lower bound                                             (out)
//   upp - upper bound                                             (out)
//
// This method returns the bounds on the values of the given variable. It
// must be defined if the problem uses integer variables.
// -------------------------------------------------------------------------
//
// void GetBounds(int i, double* low, double* upp)
//
//   i   - index of the variable                                   (in)
//   low - lower bound                                             (out)
//   upp - upper bound                                             (out)
//
// This method returns the bounds on the values of the given variable. It
// must be defined if the problem uses continuous (double) variables.
// -------------------------------------------------------------------------
//
// void Evaluate(int *var, cVector &Constr, cVector &fobjs)
//
//   var    - algorithm variables                                  (in)
//   constr - constraint values                                    (out)
//   fobjs  - objective function(s) value(s)                       (out)
//
// This method evaluates the constraints and returns the value(s) of the
// objective function(s) for a given set of integer variables. If necessary
// the problem should decode the variables. This method must be defined if
// the problem uses integer variables.
// -------------------------------------------------------------------------
//
// void Evaluate(int **var, cVector &Constr, cVector &fobjs)
//
//   var    - algorithm variables                                  (in)
//   constr - constraint values                                    (out)
//   fobjs  - objective function(s) value(s)                       (out)
//
// This method evaluates the constraints and returns the value(s) of the
// objective function(s) for a given set of integer variables. If necessary
// the problem should decode the variables. This method must be defined if
// the problem uses integer variables.
// -------------------------------------------------------------------------
//
// void Evaluate(cVector &var, cVector &Constr, cVector &fobjs)
//
//   var    - algorithm variables                                  (in)
//   constr - constraint values                                    (out)
//   fobjs  - objective function(s) value(s)                       (out)
//
// This method evaluates the constraints and returns the value(s) of the
// objective function(s) for a given set of continuous (double) variables.
// It must be defined if the problem uses continuous (double) variables.
//
// -------------------------------------------------------------------------
// Pure virtual methods:
// -------------------------------------------------------------------------
//
// eIndType GetVarType(void)
//
// This method returns the type of variables each should be used in the
// solution of the optimization problem.
// -------------------------------------------------------------------------

#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "individual.h"
#include <map>
#include <string>

#include "mat.h"

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Definition of cProblem class:
//
class cProblem
{
 friend class cProblemFactory;

 protected:
  static string      InpLabel;   // Problem label read from input file
         string      Label;      // Problem label
         int         NumVar;     // Number of variables
         int         NumObj;     // Number of objective functions
         int         NumConstr;  // Number of constraints
  static double      MinObjFunc; // Minimun objective function
  static bool        MinObjFlag; // Check if the minimum obj. func was provided
         bool*       ApproxConst;
         bool*       ApproxC;

 public:
              
  static  void       ReadProblem(void);
  static  void       ReadMinObjFunc(void);
  static  string     GetInpLabel(void) { return InpLabel; }
  static  void       ReadProblemFile(string);
  static  void       Readinitialpopsize(void);

                     cProblem(void);
                     cProblem(string);
  virtual           ~cProblem(void);
          string     GetType(void)  { return Label; }
          int        GetNumVar(void){ return NumVar; }
          int        GetNumObj(void){ return NumObj; }
          int        GetNumConstr(void){ return NumConstr; }
          double     GetMinObjFunc(void){ return MinObjFunc; }
          bool       GetMinObjFlag(void){ return MinObjFlag; }
  virtual cProblem*  GetHFP(void){ return this; }
  virtual int        VarNumRow(void){ return NumVar; }
  virtual int        VarNumCol(void){ return 1; }
  virtual int        VarNumEff(void){ return VarNumRow( )*VarNumCol( ); }
  virtual void       Init(void){ }
  virtual void       LoadReadFunc(cInpMap&){ }
  virtual void       PrintVar(int*){ }
  virtual void       DecodeVar(int*, cVector &) { }
  virtual void       WriteVar(int*){ }
  virtual void       PrintVar(int**){ }
  virtual void       DecodeVar(int**, cVector &, cMatrix &) { }
  virtual void       WriteVar(int**){ }
  virtual void       GetBounds(int, int*, int*);
  virtual void       GetVarBounds(int, double &, double&);
  virtual void       GetDblBounds(double*, double*);
  virtual void       GetNormVar(cVector&,cVector&);
  virtual void       GetNormVar(int*,cVector&);
  virtual void       GetNormVar(int**,cVector&);
  virtual void       Evaluate(int*, cVector &, cVector &);
  virtual void       Evaluate(int**, cVector &, cVector &);
  virtual void       Evaluate(cVector &, cVector &, cVector &);
  virtual void       GetApproxConstr(bool*);
  virtual void       GetApproxObj(bool*);

  virtual void       EvalExactFobj(cVector&,double &);
  virtual void       EvalExactConstraint(int,cVector&,double &);

  virtual void       EvalExactFobj(int*,double &);
  virtual void       EvalExactConstraint(int,int*,double &);

  virtual void       EvalExactFobj(int**,double &);
  virtual void       EvalExactConstraint(int,int**,double &);

  virtual void       SetApproxC(bool*) { }  
};


// -------------------------------------------------------------------------
//
// The class cProblemFactory define a object factory used to create problems
// object. The main advantage of this approuch is to enable registration of
// new problems in a decentralized manner. A new problem can be implemented and
// registred in a new file, without having to modify the base class files. The
// class consider Singletorn pattern to ensure that a unique factory object is
// instantiated.
//
// -------------------------------------------------------------------------
// Local methods:
// -------------------------------------------------------------------------
//
// cProblemFactory& GetFactory(void)
//
// This method returns the single factory object (Singletorn pattern).
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
// bool Register(string label, ProbCreatorFunc func, string fext = "")
//
//   label - problem label                                           (in)
//   func  - function that creater the problem                       (in)
//   fext  - Problem file extension.                                 (in)
//
// This method register the problem label, the creational function and the
// problem file extension related to the given problem. The default value for
// fext means that the problem has not associated input file.
// -------------------------------------------------------------------------
//
// cProblem* CreateProblem(string label, bool IsCont)
//
//   label  - problem type                                         (in)
//   IsCont - flag used to choose continuos ou discrete problem    (in)
//
// This method creates a problem object according to the given label and boolean
// flag (used on benchmark problems to differentiate continuous problem version
// to discrete problem version) and returns a pointer to this object. If the
// problem has a file extention (added in Register function), the file will be
// read.
// -------------------------------------------------------------------------

typedef cProblem* (*ProbCreatorFunc)(bool);

class cProblemFactory
{
 private: 
  static  cProblemFactory& GetFactory(void);
 

  cProblemFactory(void) { }
 ~cProblemFactory(void) { }

  std::map<std::string, ProbCreatorFunc> CreationalMapFunc;
  std::map<std::string, std::string>          FileExtMap;

 public:
          
  static  bool             Register(std::string,ProbCreatorFunc,std::string fext = "");
  static  cProblem*        CreateProblem(string,bool);
  static  std::string      GetProbFileExt(cProblem*);
};

// Defining template methods to create the problem

template<typename ProbClass>
cProblem* MakeProb(bool IsCont) 
{
  return new ProbClass;
}

template<typename ContinuousProb, typename DiscreteProb>
cProblem* MakeProb(bool IsContinuous)
{
  if (IsContinuous) return new ContinuousProb; 
  else return new DiscreteProb;
}

#endif
