// -------------------------------------------------------------------------
// saorbf.cpp - implementation of cSAORBF class.
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
// Created:      30-Jul-2019    Leonardo Gonçaves Ribeiro
//
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

#include "saorbf.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"
#include "vec.h"
#include "rbf.h"
#include "surr.h"
#include "stdpso.h"
#include "probsurr.h"
#include "sao.h"
#include "samp.h"


// -------------------------------------------------------------------------
// Class SAORBF:
//

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cSAORBF ===============================

cSAORBF :: cSAORBF(void) : cSAO( )
{
  Type     = SAORBF;
  WEI      = 0.5;
  ciclewei = 0;
  SigType  = kFOLDCV;
  GenStall = 0;
}

// ============================== ~cSAORBF ===============================

cSAORBF :: ~cSAORBF(void)
{
}

// ============================ ReadSigType =============================

void cSAORBF :: ReadSigType(istream &in)
{
  char signame[100];

  if (!Utl::ReadString(in, signame))
  {
    cout << "Error in the input of the sigma definition method." << endl;
    exit(0);
  }

  if(string(signame) == "KITAYAMA" || string(signame) == "Kitayama"
          || string(signame) == "KIT")
  {
      SigType = KITAYAMA;
  }
  else if(string(signame) == "UNIFORMKITAYAMA" || string(signame) == "UniformKitayama"
          || string(signame) == "UKIT")
  {
      SigType = KITAYAMAUNIFORM;
  }
  else if(string(signame) == "ADAPTIVESCALING" || string(signame) == "AdaptiveScaling"
          || string(signame) == "ASKIT")
  {
      SigType = ADAPTIVE_SCALING;
  }
  else if(string(signame) == "LEAVEONEOUTCROSSVALIDATION"
          || string(signame) == "LeaveOneOutCrossValidation" || string(signame) == "LOOCV")
  {
      SigType = LOOCV;
  }
  else if(string(signame) == "KFOLDCROSSVALIDATION" || string(signame) == "KFOLDCV"
          || string(signame) == "KFCV")
  {
      SigType = kFOLDCV;
  }
  else if(string(signame) == "NAKAYAMA" || string(signame) == "Nakayama"
          || string(signame) == "NAK")
  {
      SigType = NAKAYAMA;
  }
  else
  {
  cout << "Unknown sigma definition method: " << signame << endl;
  exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cSAORBF :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cSAO :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("SIGMA.TYPE"               ,  makeReadObj(cSAORBF, ReadSigType));
}

// =========================== CreateSurrogate =============================

cSURR* cSAORBF ::  CreateSurrogate(sSampData &sdata) 
{
  rbf = new cRBF;
  rbf->SetLambda(0);
  rbf->CreateModel(sdata,SigType);

  return rbf;
}

// =========================== UpdateSurrogate =============================

void cSAORBF ::  UpdateSurrogate(cVectorVec &nx, cVectorVec &ny) 
{
  rbf->UpdateModel(SigType, nx, ny); 
}

// ======================================================= End of file =====
