// -------------------------------------------------------------------------
// lam.h - file containing the definition of the Laminated class.
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
// The class cLaminated contains data and routines relevant to optimization
// of laminated composite structures.
// -------------------------------------------------------------------------

#ifndef _LAM_H
#define _LAM_H

#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <vector>

using namespace std;

#include "problem.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Laminate types:
//
typedef enum
{
  SYMMETRIC,
  BALANCED,
  SYMMETRIC_BALANCED,
  GENERAL
} eLamType;

// -------------------------------------------------------------------------
// Definition of cLaminated class:
//
class cLaminated : public cProblem
{
 protected:
  static  eLamType LamType;       // Laminate type
  static  double   dAng;          // Lamination angle variation
  static  double   MinAng;        // Minimum lamination angle
  static  double   MaxAng;        // Maximum lamination angle
  static  double   dThk;          // Thickness variation
  static  double   MinThk;        // Minimum thickness of each layer
  static  double   MaxThk;        // Maximum thickness of each layer
  static  int      MaxPlyNum;     // Maximum number of plies
  static  cVector  ThkVal;        // Custom thickness values
  static  cVector  AngVal;        // Custom angle values
  static  cVector  Thk;           // Thickness list
  static  cVector  Ang;           // Angle list
  static  int*     Mat;           // Material list
          int*     ListDim;       // List dimensions
          int      NumRow;        // Number of variable rows
          int      NumCol;        // Number of variable columns

 protected:
          bool     Decode(int**, cMatrix&);
          int      MaxContLay(cMatrix&);
          double   MaxContThk(cMatrix&);
          double   GetThickness(cMatrix&);
          double   GetAverageDensity(cMatrix&);
          double   GetWeight(cMatrix&);
          void     TMatrix(double, cMatrix&);
          void     TMatrix(double, cMatrix&, cMatrix&);
          void     QlMatrix(double*, cMatrix&);
          void     QlMatrix(double*, cMatrix&, int);	// Considers lamina degradation.
          void     QlMatrix(double*, cMatrix&, cMatrix&);
          void     QgMatrix(double, cMatrix&, cMatrix&);
          void     CalcA(cMatrix&, cMatrix&);
          void     CalcD(cMatrix&, cMatrix&);
          void     CalcAD(cMatrix&, cMatrix&, cMatrix&);
          void     CalcABD(cMatrix&, cMatrix&);
          void     CalcABD(cMatrix&, cMatrix&, cVector);
          void     MountABDG(cMatrix&, cMatrix&);
 virtual  void     PrintResult(int**, ostream&);

 public:
          void     ReadAngleRange(std::istream&);
          void     ReadAngleValues(std::istream&);
          void     ReadThickRange(std::istream&);
          void     ReadThickValues(std::istream&);
          void     ReadMaxPlyNum(std::istream&);
          void     ReadLamType(std::istream&);

                   cLaminated(void);
  virtual         ~cLaminated(void);
          int      VarNumRow(void) { return NumRow; }
          int      VarNumCol(void) { return NumCol; }
          int      VarNumEff(void);
  virtual void     Init(void);
  virtual void     LoadReadFunc(cInpMap&);
  virtual void     PrintVar(int**);
  virtual void     DecodeVar(int**, cVector &, cMatrix &);
  virtual void     WriteVar(int**,std::ostream&);
          void     GetBounds(int, int*, int*);
          double   GetThk(int id) {return Thk[id];}
};

#endif
