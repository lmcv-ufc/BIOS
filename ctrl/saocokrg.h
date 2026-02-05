// -------------------------------------------------------------------------
// saocokrg.h - file containing the definition of the cSAOCOKRG class.
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
// The cSAOKRG class implements the standard Genetic Algorithm with
// selection, crossover and mutation operators. Constrained optimization
// problems are solved using the penalty function approach.
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
// void MaximizeExpectedImprovement(cVector &xb, cKRG* Sur)
//
//   xb -  vector of the design variables of the best       (out)
//         individual/particle
//   Sur - surrogate model                                  (in)
//
// -------------------------------------------------------------------------
//
// void PostProcessingSur(int nrun, int nout, std::vector<cVector> *nrmse = 0,
//      std::vector<cVector> *rmae =0, std::vector<cVector> *nmae = 0,
//      std::vector<cVector> *error = 0, std::vector<int> *gentooconvergence = 0)
//
//   nrun  - number of optimizations                        (in)
//   nout  - number of outputs                              (in)
//   nrmse - normalized root mean squared error             (out)
//   rmae  - root mean averaged error                       (out)
//   nmae  - normalized mean averaged error                 (out)
//   error - error between hfm opt and sur opt              (out)
//   gentoconv - number of generation until convergence     (out)
//
// This method print the output data of current optimization run on output
// file.

#ifndef _SAOCOKRG_H
#define _SAOCOKRG_H

#include "mfsao.h"
#include "cokrg.h"

// -------------------------------------------------------------------------
// Definition of the Conventional GA class:
//
class cSAOCOKRG : public cMFSAO
{
 protected:
    cCOKRG           *cokrg;
    eCorrelationType  CorrType;
    double            HyperParamLow;
    double            HyperParamUpp;
    double            RhoLow;
    double            RhoUpp;

    double            AddHF;
    double            AddLF;


 public:
                    cSAOCOKRG(void);
  virtual           ~cSAOCOKRG(void);

  virtual void      LoadReadFunc(cInpMap&);
  void              ReadCorrType(std::istream&);
  void              ReadHyperParam(std::istream&);
  void              ReadRho(std::istream&);

  cSURR*    CreateSurrogate(sSampData&);
  void      UpdateSurrogate(cVectorVec&,cVectorVec&,cVectorVec&);
  void      GetFidelityNewPoint(cVectorVec&, bool&, bool&, double);
};

#endif
