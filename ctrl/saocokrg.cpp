// -------------------------------------------------------------------------
// saocokrg.cpp - implementation of cSAOCOKRG class.
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
// Created:      30-Jul-2019    Marina Alves Maia
//
// Modified:     17-Dec-2019    Marina Alves Maia
//                              Code refactoring.
// -------------------------------------------------------------------------

#include <string>
#include <vector>

using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "problem.h"
#include "sel.h"
#include "group.h"
#include "individual.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"
#include "saocokrg.h"
#include "cokrg.h"
#include "vec.h"
#include "surr.h"
#include "stdpso.h"
#include "stdga.h"
#include "probsurr.h"
#include "mfsao.h"
#include "samp.h"

// -------------------------------------------------------------------------
// Class EGO:
//

// -------------------------------------------------------------------------
// Public methods:

// =============================== cSAOKRG =================================

cSAOCOKRG :: cSAOCOKRG(void)  : cMFSAO( )
{
  Type          = SAOCOKRG;
  CorrType      = GAUSS;

  HyperParamLow = -3.0;
  HyperParamUpp =  3.0;
  RhoLow        = -5.0;
  RhoUpp        =  5.0;
}

// ============================= ~cSAOKRG ==================================

cSAOCOKRG :: ~cSAOCOKRG(void)
{
}

// ============================ ReadCorrType ===============================

void cSAOCOKRG :: ReadCorrType(istream &in)
{
  char signame[100];

  if (!Utl::ReadString(in, signame))
  {
    cout << "Error in the input of the sigma definition method." << endl;
    exit(0);
  }

  if(string(signame) == "Gauss" || string(signame) == "GAUSS"
          || string(signame) == "gauss")
  {
    CorrType = GAUSS;
  }
  else if(string(signame) == "Matern" || string(signame) == "MATERN"
          || string(signame) == "matern")
  {
    CorrType = MATERN52;
  }
  else
  {
    cout << "Unknown correlation function definition: " << signame << endl;
    exit(0);
  }
}

// ======================== ReadHyperParam =================================

void cSAOCOKRG :: ReadHyperParam(istream &in)
{
  if (!(in >> HyperParamLow) || !(in >> HyperParamUpp))
  {
    cout << "Error in the input of the hyperparameters (lower and upper bounds)" << endl;
    exit(0);
  }
}

// ======================== ReadHyperParam =================================

void cSAOCOKRG :: ReadRho(istream &in)
{
  if (!(in >> RhoLow) || !(in >> RhoUpp))
  {
    cout << "Error in the input of rho (lower and upper bounds)" << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cSAOCOKRG :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cMFSAO :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("CORRELATION.TYPE"    ,  makeReadObj(cSAOCOKRG, ReadCorrType));
  im.Insert("KRG.HYPERPARAMETERS" ,  makeReadObj(cSAOCOKRG, ReadHyperParam));
  im.Insert("KRG.RHO"             ,  makeReadObj(cSAOCOKRG, ReadRho));
}

// =========================== CreateSurrogate =============================

cSURR* cSAOCOKRG ::  CreateSurrogate(sSampData &sdata)
{
  cokrg = new cCOKRG;
  cokrg->CreateModel(sdata,CorrType,HyperParamLow,HyperParamUpp, RhoLow, RhoUpp);

  return cokrg;
}

// =========================== UpdateSurrogate =============================

void cSAOCOKRG ::  UpdateSurrogate(cVectorVec &nx, cVectorVec &ny, cVectorVec &nylf)
{
  cokrg->UpdateModel(CorrType,nx,ny,nylf, AddHF, AddLF);
}

// =========================== UpdateSurrogate =============================

void cSAOCOKRG ::  GetFidelityNewPoint(cVectorVec &nx, bool &addhf, bool &addlf, double)
{
    AddHF = addhf = 1;
    AddLF = addlf = 1;
}

// ======================================================= End of file =====
