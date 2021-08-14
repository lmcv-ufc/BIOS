// -------------------------------------------------------------------------
// metaprob.h - file containing the definition cMetaProblem class.
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
// This class define the optimization problems used in the meta-optimization
// problem. This class stores the problema data, as input parameters and
// problem best objective values. Moreover, the class implements the
// WriteProblem method, used to write the problema input file.   
// -------------------------------------------------------------------------

#ifndef _METAPROB_H
#define _METAPROB_H

#include <fstream>

using namespace std;

// -------------------------------------------------------------------------
// MetaProblem types:
//
typedef enum 
{ 
  LAMINATE_PLATE_SET01, 
  LAMINATE_PLATE_SET02,
  LAMINATE_PLATE_SET_MIN_WEIGHT_01,
  LAMINATE_PLATE_SET_MIN_WEIGHT_02,
  LAMINATE_PLATE_SET_MIN_COST_01
} eMetaProbType;

// -------------------------------------------------------------------------
// Objective Function Type types:
//
typedef enum 
{ 
  MAXIMIZATION,
  MINIMIZATION
} eObjFuncType;

// ------------------------------------------------------------------------
// Definition of the metaoptimization problem set class:
//
class cMetaProblem
{
 protected:
  static eMetaProbType InpType;     // MetaProblems type read from input file
         eMetaProbType Type;        // MetaProb type 
	 eObjFuncType  ObjType;     // Objective Function Type
	 int           NumProb;     // Number of problems to e metaoptimized
	 double*       BestObjFunc; // Best values of each problem
  static int           DAng;

 public:
  static void          ReadMetaProblem(std::istream&);
  static eMetaProbType GetInpType(void) { return InpType; }
  static cMetaProblem* CreateMetaProblem(eMetaProbType);

         double        GetBestObjFunc(int id) {return BestObjFunc[id];}
                       cMetaProblem(void);
  virtual              ~cMetaProblem(void);

          int          GetNumProb(void) {return NumProb;} 
	  eObjFuncType GetObjType(void) {return ObjType;}
  virtual void         WriteProblem(string&,fstream&,int) = 0;
};

//-------------------------------------------------------------------------
// Definition of the laminate plate problems set01 class:
//

class cLamPltSet01 : public cMetaProblem
{
 public:
           cLamPltSet01(void);
           ~cLamPltSet01(void);
  void     WriteProblem(string&,fstream&,int);
};

//-------------------------------------------------------------------------
// Definition of the laminate plate problems set02 class:
//

class cLamPltSet02 : public cMetaProblem
{
 public:
           cLamPltSet02(void);
           ~cLamPltSet02(void);
  void     WriteProblem(string&,fstream&,int);
};

//-------------------------------------------------------------------------
// Definition of the laminate plate problems of weigth minimization 01 class:
//

class cLamPltSetMinWeight01 : public cMetaProblem
{
 public:
           cLamPltSetMinWeight01(void);
           ~cLamPltSetMinWeight01(void);
  void     WriteProblem(string&,fstream&,int);
};

//-------------------------------------------------------------------------
// Definition of the laminate plate problems of weigth minimization 02 class:
//

class cLamPltSetMinWeight02 : public cMetaProblem
{
 public:
           cLamPltSetMinWeight02(void);
           ~cLamPltSetMinWeight02(void);
  void     WriteProblem(string&,fstream&,int);
};

//-------------------------------------------------------------------------
// Definition of the laminate plate problems of cost minimization 01  class:
//

class cLamPltSetMinCost01 : public cMetaProblem
{
 public:
           cLamPltSetMinCost01(void);
           ~cLamPltSetMinCost01(void);
  void     WriteProblem(string&,fstream&,int);
};

#endif
