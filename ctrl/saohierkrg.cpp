// -------------------------------------------------------------------------
// saohierkrg.cpp - implementation of cSAOHIERKRG class.
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
#include "saohierkrg.h"
#include "hierkrg.h"
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

cSAOHIERKRG :: cSAOHIERKRG(void)  : cMFSAO( )
{
  Type          = SAOHKRG;
  CorrType      = GAUSS;

  HyperParamLow = -3.0;
  HyperParamUpp =  3.0;
}

// ============================= ~cSAOKRG ==================================

cSAOHIERKRG :: ~cSAOHIERKRG(void)
{
}

// ============================ ReadCorrType ===============================

void cSAOHIERKRG :: ReadCorrType(istream &in)
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

void cSAOHIERKRG :: ReadHyperParam(istream &in)
{
  if (!(in >> HyperParamLow) || !(in >> HyperParamUpp))
  {
    cout << "Error in the input of the hyperparameters (lower and upper bounds)" << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cSAOHIERKRG :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cMFSAO :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("CORRELATION.TYPE"    ,  makeReadObj(cSAOHIERKRG, ReadCorrType));
  im.Insert("KRG.HYPERPARAMETERS" ,  makeReadObj(cSAOHIERKRG, ReadHyperParam));
}

// =========================== CreateSurrogate =============================

cSURR* cSAOHIERKRG ::  CreateSurrogate(sSampData &sdata)
{
  hkrg = new cHIERKRG;
  hkrg->CreateModel(sdata,CorrType,HyperParamLow,HyperParamUpp);

  return hkrg;
}

// =========================== UpdateSurrogate =============================

void cSAOHIERKRG ::  UpdateSurrogate(cVectorVec &nx, cVectorVec &ny, cVectorVec &nylf)
{
  hkrg->UpdateModel(CorrType,nx,ny,nylf, AddHF, AddLF);
}

// =========================== UpdateSurrogate =============================

void cSAOHIERKRG ::  GetFidelityNewPoint(cVectorVec &nx, bool &addhf, bool &addlf, double currbest)
{
    addhf = 1;
    addlf = 1;

    double ErrMax = 0.01;

    double inf1, inf2, error;

    if (InfillCriteria == EVALUATE_VF_EXPECTED_IMPROVEMENT)
    {
        inf1 = hkrg -> EI1(nx[0], currbest);
        inf2 = hkrg -> EI2(nx[0], currbest);

        if (inf1 != 0.0 || inf2 != 0.0)
        {
            if (inf1 == 0.0 && inf2 != 0.0)  addlf = 0;
            else if (inf2 == 0.0 && inf1 != 0.0)  addhf = 0;

            if (addhf == 1 && addlf == 1)
            {
                error = abs((inf1 - inf2)/inf2);
                if (error > ErrMax)
                {
                    if(inf1 > inf2) addhf = 0;
                    else addlf = 0;
                }
            }
        }
    }

    if (InfillCriteria == EVALUATE_VF_PROBABILITY_IMPROVEMENT)
    {
        inf1 = hkrg -> PI1(nx[0], currbest);
        inf2 = hkrg -> PI2(nx[0], currbest);

        if (inf1 != 0.0 || inf2 != 0.0)
        {
            if (inf1 == 0.0 && inf2 != 0)  addlf = 0;
            else if (inf2 == 0.0 && inf1 != 0)  addhf = 0;

            if (addhf == 1 && addlf == 1)
            {
                error = abs((inf1 - inf2)/inf2);
                if (error > ErrMax)
                {
                    if(inf1 > inf2) addhf = 0;
                    else addlf = 0;
                }
            }
        }
    }

    if (InfillCriteria == EVALUATE_VF_LOWER_CONFIDENCE_BOUND)
    {
        inf1 = hkrg -> LCB1(nx[0]);
        inf2 = hkrg -> LCB2(nx[0]);

        /*cVector t(5), td(5);

        t  = hkrg -> GetBestTheta(0);
        td = hkrg -> GetBestThetad(0);

        double s1 = hkrg -> SSqrlSur(nx[0], t, td, 0);
        double s2 = hkrg -> SSqrSur(nx[0], t, td, 0);

        cout << "inf1 = " << inf1 << endl;
        cout << "inf2 = " << inf2 << endl;

        cout << "s1 = " << s1 << endl;
        cout << "s2 = " << s2 << endl;*/

        if (inf1 != 0.0 || inf2 != 0.0)
        {
            if (addhf == 1 && addlf == 1)
            {
                error = abs((inf1 - inf2)/inf2);
                if (error > ErrMax)
                {
                    if(inf1 < inf2) addhf = 0;
                    else addlf = 0;
                }
            }
        }
        else
        {
            if(inf1 < inf2) addhf = 0;
            else addlf = 0;
        }
    }

    if (NestSamp) addlf = 1;

    if (addlf == 1) cout << "ADD LOW FID" << endl;
    if (addhf == 1) cout << "ADD HIG FID" << endl;

    AddLF = addlf;
    AddHF = addhf;
}

// ======================================================= End of file =====
