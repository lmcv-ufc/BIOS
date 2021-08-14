// -------------------------------------------------------------------------
// benchmark.h - file containing the definition cBenchmark class.
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
// The aim of this class is to solve classical benchmark problems
// of numerical optimization and demonstrate how to implement a optimization
// problem using the BIOS system.
// -------------------------------------------------------------------------

#ifndef _BENCHMARK_H
#define _BENCHMARK_H

#include <vector>
#include "vec.h"
#include "mat.h"
#include "matvec.h"
#include "sysmat.h"
#include "problem.h"
#include "group.h"

using namespace std;

// ------------------------------------------------------------------------
// Definition of the discrete benchmark class:
//
class cBenchmark : public cProblem
{
  protected:
  static int InpNumVar;

  public:
  static void  ReadInpNumVar(void);

               cBenchmark(void);
  virtual     ~cBenchmark(void);
};

// ------------------------------------------------------------------------
// Definition of the discrete benchmark class:
//
class cBenchDiscrete : public cBenchmark
{
  protected:
  int     *ListDim; // Number of discrete values for each variable
  cVector *List;    // Lists of discrete values for each variable

  public:
           cBenchDiscrete(void);
  virtual ~cBenchDiscrete(void);
  void     PrintVar(int*);
  void     DecodeVar(int*, cVector &);
  void     WriteVar(int*,std::ostream&);
  void     GetBounds(int, int*, int*);
  void     GetVarBounds(int, double&, double&);

};

//-------------------------------------------------------------------------
// Definition of the continuous benchmark class:
//
class cBenchContinuous : public cBenchmark
{
  protected:
  double *Low;
  double *Upp;

  public:
           cBenchContinuous(void);
  virtual ~cBenchContinuous(void);
  void     GetDblBounds(double*, double*);
};

//--------------------------------------------------------------------------
// Definition of PeaksC class:
//
// Optimization problem Peaks. This problem deals with the optimization
// of a well-known function from base MATLAB [1], used for example pur-
// poses. In this class, the problem is solved using continuous variables.
//
// [1] https://www.mathworks.com/help/matlab/ref/peaks.html
//

class cPeaksC : public cBenchContinuous
{
  public:
         cPeaksC(void);
        ~cPeaksC(void){}
  void Evaluate(cVector & ,cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of PeaksD class:
//
// Optimization problem Peaks. This problem deals with the optimization
// of a well-known function from base MATLAB [1], used for example pur-
// poses. In this class, the problem is solved using discrete variables.
//
// [1] https://www.mathworks.com/help/matlab/ref/peaks.html
//

class cPeaksD : public cBenchDiscrete
{
  public:
          cPeaksD(void);
         ~cPeaksD(void){ }
  void  Evaluate(int *, cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of BraninC class:
//
// Optimization problem Branin. This problem is a very well-known optimization
// problem, usually employed for example purposes in SAO [1].  In this class,
// the problem is solved using continuous variables.
//
// [1] FORRESTER, Alexander et al. Engineering design via surrogate modelling: a
//     practical guide. John Wiley & Sons, 2008.
//

class cBraninC : public cBenchContinuous
{
  public:
         cBraninC(void);
        ~cBraninC(void){}
  void   Evaluate(cVector & ,cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of BraninD class:
//
// Optimization problem Branin. This problem is a very well-known optimization
// problem, usually employed for example purposes in SAO [1].  In this class,
// the problem is solved using discrete variables.
//
// [1] FORRESTER, Alexander et al. Engineering design via surrogate modelling: a
//     practical guide. John Wiley & Sons, 2008.
//

class cBraninD : public cBenchDiscrete
{
  public:
         cBraninD(void);
        ~cBraninD(void){}
  void   Evaluate(int *, cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of cHart3C class:
//
// Optimization problem Hartmann 3. This problem is a very well-known optimization
// problem, usually employed for example purposes in SAO [1].  In this class, the
// problem is solved using continuous variables.
//
// [1] Jones, D. et al. Efficient global optimization of expensive black box
//     functions. Journal of Global Optimization, v. 13, p. 455-492, 1998.
//

class cHart3C : public cBenchContinuous
{
  public:
         cHart3C(void);
        ~cHart3C(void){}
  void Evaluate(cVector & ,cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of cHart3D class:
//
// Optimization problem Hartmann 3. This problem is a very well-known optimization
// problem, usually employed for example purposes in SAO [1].  In this class, the
// problem is solved using discrete variables.
//
// [1] Jones, D. et al. Efficient global optimization of expensive black box
//     functions. Journal of Global Optimization, v. 13, p. 455-492, 1998.
//

class cHart3D : public cBenchDiscrete
{
  public:
          cHart3D(void);
         ~cHart3D(void){ }
  void  Evaluate(int *, cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of cHart6C class:
//
// Optimization problem Hartmann 6. This problem is a very well-known optimization
// problem, usually employed for example purposes in SAO [1].  In this class, the
// problem is solved using continuous variables.
//
// [1] Jones, D. et al. Efficient global optimization of expensive black box
//     functions. Journal of Global Optimization, v. 13, p. 455-492, 1998.
//

class cHart6C : public cBenchContinuous
{
  public:
         cHart6C(void);
        ~cHart6C(void){}
  void Evaluate(cVector & ,cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of cHart6D class:
//
// Optimization problem Hartmann 6. This problem is a very well-known optimization
// problem, usually employed for example purposes in SAO [1].  In this class, the
// problem is solved using discrete variables.
//
// [1] Jones, D. et al. Efficient global optimization of expensive black box
//     functions. Journal of Global Optimization, v. 13, p. 455-492, 1998.
//

class cHart6D : public cBenchDiscrete
{
  public:
          cHart6D(void);
         ~cHart6D(void){ }
  void  Evaluate(int *, cVector &, cVector &);
};

// ------------------------------------------------------------------------
// Definition of RastriginC class
//
// Optimization problem Rastrigin. This problem is a very well-known optimization
// problem, with multiple local minima [1]. In this class, the problem is solved
// using continuous variables.
//
// [1] Dieterich, J. M.; Hartke, B. Empirical review of standard benchmark functions
//     using evolutionary global optimization.
//     Acessed from: https://arxiv.org/pdf/1207.4318.pdf
//

class cRastriginC : public cBenchContinuous
{
  public:
          cRastriginC(void);
         ~cRastriginC(void){ }
  void  Evaluate(cVector &, cVector &, cVector &);
};

// ------------------------------------------------------------------------
// Definition of RastriginD class
//
// Optimization problem Rastrigin. This problem is a very well-known optimization
// problem, with multiple local minima [1]. In this class, the problem is solved
// using discrete variables.
//
// [1] Dieterich, J. M.; Hartke, B. Empirical review of standard benchmark functions
//     using evolutionary global optimization.
//     Acessed from: https://arxiv.org/pdf/1207.4318.pdf
//

class cRastriginD : public cBenchDiscrete
{
  public:
          cRastriginD(void);
         ~cRastriginD(void){ }
  void  Evaluate(int *, cVector &, cVector &);
};

//--------------------------------------------------------------------------
// Definition of ConstrainedBraninC class:
//
// Optimization problem Constrained Branin. This problem is a modification
// over a very well-known optimization problem, usually employed for example
// purposes in SAO [1]. In this class, the problem is solved using continuous
// variables.
//
// [1] Forrester, Alexander et al. Engineering design via surrogate modelling: a
//     practical guide. John Wiley & Sons, 2008.
//

class cConstrainedBraninC : public cBenchContinuous
{
  public:
         cConstrainedBraninC(void);
        ~cConstrainedBraninC(void){}
  void   Evaluate(cVector & ,cVector &, cVector &);
  void  EvalExactFobj(cVector&,double &){ };
  void  EvalExactConstraint(int, cVector&, double &);
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*c) { c[0] = 1; }
};

//--------------------------------------------------------------------------
// Definition of ConstrainedBraninD class:
//
// Optimization problem Constrained Branin. This problem is a modification
// over a very well-known optimization problem, usually employed for example
// purposes in SAO [1]. In this class, the problem is solved using discrete
// variables.
//
// [1] Forrester, Alexander et al. Engineering design via surrogate modelling: a
//     practical guide. John Wiley & Sons, 2008.
//

class cConstrainedBraninD : public cBenchDiscrete
{
  public:
         cConstrainedBraninD(void);
        ~cConstrainedBraninD(void){}
  void  Evaluate(int *, cVector &, cVector &);
  void  EvalExactFobj(int *,double &){ };
  void  EvalExactConstraint(int, int *, double &){ };
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*c){ c[0]  = 1; }
};

//--------------------------------------------------------------------------
// Definition of Kit5C class:
//
// Kitayama's constrained problem. This problem was solved by Kitayama et al. (2011) [1]
// using their own proposed SAO algorithm. In this class, the problem is
// solved using continuous variables.
//
// Refs:
// [1] Kitayama, Satoshi; Arakawa, Masao; Yamazaki, Koetsu. Sequential Approximate
//     Optimization using Radial Basis Function network for engineering optimization.
//     Optimization and Engineering, v. 12, p. 535-556, 2011.
//

class cKit5C : public cBenchContinuous
{
  public:
         cKit5C(void);
        ~cKit5C(void){}
  void  Evaluate(cVector & ,cVector &, cVector &);
  void  EvalExactFobj(cVector&,double &);
  void  EvalExactConstraint(int, cVector&, double &);
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*);
};

//--------------------------------------------------------------------------
// Definition of Kit5D class:
//
// Kitayama's constrained problem. This problem was solved by Kitayama et al. (2011) [1]
// using their own proposed SAO algorithm. In this class, the problem is
// solved using discrete variables.
//
// Refs:
// [1] Kitayama, Satoshi; Arakawa, Masao; Yamazaki, Koetsu. Sequential Approximate
//     Optimization using Radial Basis Function network for engineering optimization.
//     Optimization and Engineering, v. 12, p. 535-556, 2011.
//

class cKit5D : public cBenchDiscrete
{
  public:
         cKit5D(void);
        ~cKit5D(void){}
  void  Evaluate(int *, cVector &, cVector &);
  void  EvalExactFobj(int *,double &);
  void  EvalExactConstraint(int, int *, double &);
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*);
};

//--------------------------------------------------------------------------
// Definition of 3BarTrussC class:
//
// Optimization problem Three-Bar Truss. This problem is an engineering
// optimization problem, which deals with the mass optimization of a three-bar
// truss considering stress related constraints [1]. In this class, the problem
// is solved using continuous variables.
//
// [1] Rao, S. S. Engineering Optimization: Theory and Practice. John Wiley &
//     Sons, Inc. 4th Edition.
//

class c3BarTrussC : public cBenchContinuous
{
  public:
         c3BarTrussC(void);
        ~c3BarTrussC(void){}
  void   Evaluate(cVector & ,cVector &, cVector &);
  void   EvalExactFobj(cVector&,double &){ };
  void   EvalExactConstraint(int, cVector&, double &){ };
  void   GetApproxObj(bool*o) { o[0] = 1; }
  void   GetApproxConstr(bool*c) { c[0] = 1; c[1] = 1; c[2] = 1; }
};

//--------------------------------------------------------------------------
// Definition of 3BarTrussC class:
//
// Optimization problem Three-Bar Truss. This problem is an engineering
// optimization problem, which deals with the mass optimization of a three-bar
// truss considering stress related constraints [1]. In this class, the problem
// is solved using discrete variables.
//
// [1] Rao, S. S. Engineering Optimization: Theory and Practice. John Wiley &
//     Sons, Inc. 4th Edition.
//

class c3BarTrussD : public cBenchDiscrete
{
  public:
         c3BarTrussD(void);
        ~c3BarTrussD(void){}
  void  Evaluate(int* ,cVector &, cVector &);
  void  EvalExactFobj(int *,double &){ };
  void  EvalExactConstraint(int, int *, double &){ };
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*c){ c[0] = 1; c[1] = 1; c[2] = 1; }
};

//--------------------------------------------------------------------------
// Definition of NowackiBeamC class:
//
// Optimization problem Nowacki Beam. This problem is an engineering optimization
// problem which deals with the sectional area minimization of a cantilever
// beam, subjected to displacement and stress related constraints [1]. In this
// class, the problem is solved using continuous variables.
//
// [1] Nowacki, H. Modelling of design decisions for CAD, CAD Modelling, Systems
//     Engineering, CAD-Systems. Lecture notes in computer science.
//

class cNowackiBeamC : public cBenchContinuous
{
  public:
         cNowackiBeamC(void);
        ~cNowackiBeamC(void){}
  void   Evaluate(cVector & ,cVector &, cVector &);
  void   EvalExactFobj(cVector&,double &){ };
  void   EvalExactConstraint(int, cVector&, double &){ };
  void   GetApproxObj(bool*o) { o[0] = 1; }
  void   GetApproxConstr(bool*c) { c[0] = 1; c[1] = 1; }
};

//--------------------------------------------------------------------------
// Definition of NowackiBeamD class:
//
// Optimization problem Nowacki Beam. This problem is an engineering optimization
// problem which deals with the sectional area minimization of a cantilever
// beam, subjected to displacement and stress related constraints [1]. In this
// class, the problem is solved using discrete variables.
//
// [1] Nowacki, H. Modelling of design decisions for CAD, CAD Modelling, Systems
//     Engineering, CAD-Systems. Lecture notes in computer science.
//

class cNowackiBeamD : public cBenchDiscrete
{
  public:
         cNowackiBeamD(void);
        ~cNowackiBeamD(void){}
  void  Evaluate(int* ,cVector &, cVector &);
  void  EvalExactFobj(int *,double &){ };
  void  EvalExactConstraint(int, int *, double &){ };
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*c) { c[0] = 1; c[1] = 1; }
};

//--------------------------------------------------------------------------
// Definition of BeamC class:
//
// Optimization of a 5-segment cantilever beam. This problem is an engineering
// optimization problem which deals with the volume minimization of a cantilever
// beam with five different segments (with different width and height), subjected
// to displacement, geometry and stress related constraints [1]. In this class,
// the problem is solved using continuous variables.
//
// [1] Thanedar, P. B.; Vanderplaats, G. N. Survey of Discrete Variable Optimiza-
//     tion for Structural Design. Jornal of Structural Engineering, v. 121,
//     p. 301-306.
//

class cBeamC : public cBenchContinuous
{
  public:
         cBeamC(void);
        ~cBeamC(void){}
  void  Evaluate(cVector & ,cVector &, cVector &);
  void  EvalExactFobj(cVector&,double &){ };
  void  EvalExactConstraint(int, cVector&, double &);
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*);
};

//--------------------------------------------------------------------------
// Definition of BeamD class:
//
// Optimization of a 5-segment cantilever beam. This problem is an engineering
// optimization problem which deals with the volume minimization of a cantilever
// beam with five different segments (with different width and height), subjected
// to displacement, geometry and stress related constraints [1]. In this class,
// the problem is solved using discrete variables.
//
// [1] Thanedar, P. B.; Vanderplaats, G. N. Survey of Discrete Variable Optimiza-
//     tion for Structural Design. Jornal of Structural Engineering, v. 121,
//     p. 301-306.
//

class cBeamD : public cBenchDiscrete
{
  public:
         cBeamD(void);
        ~cBeamD(void){}
  void  Evaluate(int* ,cVector &, cVector &);
  void  EvalExactFobj(int *,double &){ };
  void  EvalExactConstraint(int, int *, double &);
  void  GetApproxObj(bool*o) { o[0] = 1; }
  void  GetApproxConstr(bool*c);
};

// ------------------------------------------------------------------------
// Definition of CONSTRC class:
//
class cCONSTRC : public cBenchContinuous
{
 public:
          cCONSTRC(void);
         ~cCONSTRC(void){ }
  void  Evaluate(cVector&, cVector &, cVector &);
};

// ------------------------------------------------------------------------
// Definition of TNKC class:
//
class cTNKC : public cBenchContinuous
{
 public:
          cTNKC(void);
         ~cTNKC(void){ }
  void  Evaluate(cVector&, cVector &, cVector &);
};

// ------------------------------------------------------------------------
// Definition of SCHC class:
//
class cSCHC : public cBenchContinuous
{
 public:
          cSCHC(void);
         ~cSCHC(void){ }
  void  Evaluate(cVector&, cVector &, cVector &);
};



// ------------------------------------------------------------------------
// Definition of ZDT1C class:
//
class cZDT1C : public cBenchContinuous
{
 public:
          cZDT1C(void);
         ~cZDT1C(void){ }
   void  Evaluate(cVector&, cVector &, cVector &);
};


// ------------------------------------------------------------------------
// Definition of ZDT6C class:
//
class cZDT6C : public cBenchContinuous
{
 public:
          cZDT6C(void);
         ~cZDT6C(void){ }
  void  Evaluate(cVector&, cVector &, cVector &);
};



// ------------------------------------------------------------------------
// Definition of KURC class:
//
class cKURC : public cBenchContinuous
{
 public:
          cKURC(void);
         ~cKURC(void){ }
   void  Evaluate(cVector&, cVector &, cVector &);
};

#endif
