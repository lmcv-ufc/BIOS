// -------------------------------------------------------------------------
// problike.h - file containing the definition of the cProbLikelihood class.
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
// The class cProbLikelihood defines the Maximum Likelihood Estimator (MLE)
// optimization problem. This is treated as a single-objective, continuous,
// unconstrained optimization problem.
//
// Refs:
// JONES, D. R. A taxonomy of global optimization methods based on response
// surfaces. Journal of Global Optimization, v. 21, n. 4, p. 345–383, Dec 2001.
// ISSN 1573-2916.
//
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
#ifndef _PROBLIKE_H
#define _PROBLIKE_H

#include <iostream>
#include <vector>
#include "vec.h"
#include "group.h"

#include "benchmark.h"
#include "optalg.h"
#include "stdpso.h"

using namespace std;

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;
class cKRG;

// ------------------------------------------------------------------------
// Definition of cProbLikelihood class:

class cProbLikelihood : public cBenchContinuous
{
  protected:
    cKRG             *Surr;
    int               Out;

  public:

            cProbLikelihood(cKRG *smod, int out);

         void    Evaluate(cVector &, cVector &, cVector &);
};

#endif
