// -------------------------------------------------------------------------
// problike.cpp - implementation of cProbLike class.
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
// Created:      23-Dec-2019    Marina Alves Maia
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//#include <bits/stdc++.h>

#ifdef _OMP_
#include "omp.h"
#endif

#include "krg.h"
#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "matvec.h"
#include "gblvar.h"
#include "problike.h"
#include "penalty.h"

using namespace std;

// ============================== ProbLikelihood ================================

cProbLikelihood :: cProbLikelihood(cKRG *SurMod, int out)
{
  NumVar    = SurMod->GetSampData( ).NumVar;
  NumConstr = 0;
  NumObj    = 1;

  Surr = SurMod;

  Low = new double[NumVar];
  Upp = new double[NumVar];

  for (int i = 0; i < NumVar; i++)
  {
    Low[i] = Surr->HPlow;
    Upp[i] = Surr->HPupp;
  }

  Out = out;
}


// ============================= Evaluate ================================

void cProbLikelihood :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation.

  fobjs[0] = Surr->Eval(x, Out);
}

// ======================================================= End of file =====
