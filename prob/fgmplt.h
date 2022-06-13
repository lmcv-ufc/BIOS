// -------------------------------------------------------------------------
// fgmplt.h - file containing the definition of the FG Plate class.
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
// The class cFGMPlt and its subclass implements numerical functionally graded
// plate problems.
// -------------------------------------------------------------------------

#ifndef _FGMPLT_H
#define _FGMPLT_H

#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <vector>
using namespace std;

#include "problem.h"
#include "fgm.h"
#include "rbf.h"
#include "krg.h"

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;
class cMatrix;

// -------------------------------------------------------------------------
// Definition of FGMPlate class:
//
class cFGMPlate : public cFGM
{
 public:
                  cFGMPlate(void);
  virtual        ~cFGMPlate(void);
};

// -------------------------------------------------------------------------
// Definition of SquarePlateBuckFGM class:
//
// Maximization of the buckling critical load of a FG square plate, under
// a ceramic volume fraction-related constraint.
// This problem is described and solved in Ribeiro et al. [1] using a
// RBF-based SAO.
//
// Refs:
// [1] Ribeiro, L. G.; Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Surrogate based
//     optimization of functionally graded plates using radial basis functions.
//     Composite Structures, v. 252, 2020.
//
class cSquarePlateBuckFGM : public cFGMPlate
{
 protected:
  void     Analysis(cVector, double &);

 public:
           cSquarePlateBuckFGM(void);
          ~cSquarePlateBuckFGM(void) { }
  void     Evaluate(cVector & ,cVector &, cVector &);
  void     EvalExactFobj(cVector&,double &){ };
  void     EvalExactConstraint(int, cVector&, double &);
  void     GetApproxObj(bool*o) { o[0] = 1; }
  void     GetApproxConstr(bool*);
};

// -------------------------------------------------------------------------
// Definition of SquarePlateTridirBuckFGM class:
//
// Maximization of the buckling critical load of a FG square plate, under
// a ceramic volume fraction-related constraint. In this problem, the
// gradation is given in three directions [1].
//
// Refs:
// [1] Do, D.; Nguyen-Xuan, H.; Jaehong, L. Material optimization of
//     tri-directional functionally graded plates by using deep neural
//     network and isogeometric multimesh design approach. Applied
//     Mathematical Modelling, v. 87, 2020.
//
class cSquarePlateTridirBuckFGM : public cFGMPlate
{
 protected:
  void     Analysis(cVector, double &);

 public:
           cSquarePlateTridirBuckFGM(void);
          ~cSquarePlateTridirBuckFGM(void) { }
  void     Evaluate(cVector & ,cVector &, cVector &);
  void     EvalExactFobj(cVector&,double &){ };
  void     EvalExactConstraint(int, cVector&, double &);
  void     GetApproxObj(bool*o) { o[0] = 1; }
  void     GetApproxConstr(bool*);
};

// -------------------------------------------------------------------------
// Definition of SquarePlateTridirBuckFGM class:
//
// Maximization of the buckling critical load of a FG shell, under
// a ceramic volume fraction-related constraint. In this problem, the
// gradation is given in three directions.
//
class cShellTridirBuckFGM : public cFGMPlate
{
 protected:
  double   W_MObj;
  void     Analysis(cVector, double &);

 public:
  void     ReadW(std::istream&);

           cShellTridirBuckFGM(void);
          ~cShellTridirBuckFGM(void) { }
  void     Evaluate(cVector & ,cVector &, cVector &);
  void     EvalVolumeHole(cVector, double, cVector, double&);
  void     EvalExactFobj(cVector&,double &){ };
  void     EvalExactConstraint(int, cVector&, double &);
  void     GetApproxObj(bool*o) { o[0] = 1; }
  void     GetApproxConstr(bool*);
  void     LoadReadFunc(cInpMap&);
  void     Write(cVector&,ostream&);
};

// -------------------------------------------------------------------------
// Definition of SquarePlateFreqFGM class:
//
// Maximization of the fundamental frequency of a FG square plate, under
// a ceramic volume fraction-related constraint.
// This problem is described and solved in Ribeiro et al. [1] using a
// RBF-based SAO, and in Maia et al. [2] using a KRG-based SAO.
//
// Refs:
// [1] Ribeiro, L. G.; Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Surrogate based
//     optimization of functionally graded plates using radial basis functions.
//     Composite Structures, v. 252, 2020.
//
// [2] Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Kriging-based optimization of
//     functionally graded structures. Structural and Multidisciplinary Optimization.
//     2021.
//
class cSquarePlateFreqFGM : public cFGMPlate
{
 protected:
  void     Analysis(cVector, double &);

 public:
           cSquarePlateFreqFGM(void);
          ~cSquarePlateFreqFGM(void) { }
  void     Evaluate(cVector & ,cVector &, cVector &);
  void     EvalExactFobj(cVector&,double &){ };
  void     EvalExactConstraint(int, cVector&, double &);
  void     GetApproxObj(bool*o) { o[0] = 1; }
  void     GetApproxConstr(bool*);
};

// -------------------------------------------------------------------------
// Definition of cSquarePlateHoleBuckFGM class:
//
// Maximization of the critical buckling load of a FG square plate with a hole
// in its center, under ceramic volume fraction and mass-related constraints.
// This problem is described and solved in Ribeiro et al. [1] using a
// RBF-based SAO, and in Maia et al. [2] using a KRG-based SAO.
//
// Refs:
// [1] Ribeiro, L. G.; Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Surrogate based
//     optimization of functionally graded plates using radial basis functions.
//     Composite Structures, v. 252, 2020.
//
// [2] Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Kriging-based optimization of
//     functionally graded structures. Structural and Multidisciplinary Optimization.
//     2021.
//
class cSquarePlateHoleBuckFGM : public cFGMPlate
{
 protected:
  void     Analysis(cVector, double &);

 public:
           cSquarePlateHoleBuckFGM(void);
          ~cSquarePlateHoleBuckFGM(void) { }
  void     Evaluate(cVector & ,cVector &, cVector &);
  void     EvalExactFobj(cVector&,double &){ };
  void     EvalExactConstraint(int, cVector&, double &);
  void     GetApproxObj(bool*o) { o[0] = 1; }
  void     GetApproxConstr(bool*);
};

// -------------------------------------------------------------------------
// Definition of cSquarePlateFreqFrancoFGM class:
//
// Maximization of the fundamental frequency of a FG square plate, under
// fundamental frequency constraints.
// This problem is described and solved in Maia et al. [1] using a
// KRG-based SAO.
//
// Refs:
//
// [1] Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Kriging-based optimization of
//     functionally graded structures. Structural and Multidisciplinary Optimization.
//     2021.
//
class cSquarePlateFreqFrancoFGM : public cFGMPlate
{
 protected:
  void     Analysis(cVector, double &);

 public:
           cSquarePlateFreqFrancoFGM(void);
          ~cSquarePlateFreqFrancoFGM(void) { }
   void    Evaluate(cVector & ,cVector &, cVector &);
   void    EvalExactConstraint(int, cVector&, double &){ };
   void    GetApproxObj(bool*o) { o[0] = 1; }
   void    GetApproxConstr(bool*);
};

// -------------------------------------------------------------------------
// Definition of cScoordelisFGM class:
//
// Maximization of the mass of a shallow shell, under
// displacement constraints.
// This problem is described and solved in Maia et al. [1] using a
// KRG-based SAO.
//
// Refs:
//
// [1] Maia, M. A.; Parente Jr., E.; Melo, A. M. C. Kriging-based optimization of
//     functionally graded structures. Structural and Multidisciplinary Optimization.
//     2021.
//
class cScoordelisFGM : public cFGMPlate
{
 protected:
  void     AnalysisDispStress(cVector, double &, double &, double &);

 public:
           cScoordelisFGM(void);
          ~cScoordelisFGM(void) { }
  void    Evaluate(cVector & ,cVector &, cVector &);
  void    EvalExactConstraint(int, cVector&, double &){ };
  void    EvalExactFobj(cVector&, double &);
  void    GetApproxConstr(bool*);
  void    GetApproxObj(bool*);
};

// -------------------------------------------------------------------------
// Definition of cCircularPlateFreqFGM class:
//
// Maximization of the fundamental frequency of a FG circular plate,
// under a ceramic volume fraction and mass related constraints.
// This problem is described and solved in Barroso et al. [1] using a
// RBF and Kriging based SAO.
//
// Refs:
//
// [1] Barroso, E. S.; Ribeiro, L. G.; Maia, M. A.; Rocha, I. B. C. M.;
//     Parente Jr., E.; Melo, A. M. C. BIOS: An object-oriented framework
//     for Surrogate-Based Optimization using bio-inspired algorithms.
//     Structural and Multidisciplinary Optimization. Submitted for
//     publication.
//

class cCircularPlateFreqFGM : public cFGMPlate
{
 protected:
  void     Analysis(cVector, double &);

 public:
           cCircularPlateFreqFGM(void);
          ~cCircularPlateFreqFGM(void) { }
  void    Evaluate(cVector & ,cVector &, cVector &);
  void    EvalExactConstraint(int, cVector&, double &);
  void    GetApproxObj(bool*o) { o[0] = 1; }
  void    GetApproxConstr(bool*);
};

#endif
