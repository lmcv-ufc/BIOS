// -------------------------------------------------------------------------
// surr.cpp - implementation of cSURR class.
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
// Created:      28-May-2019    Marina Alves Maia
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//#include <bits/stdc++.h>

#include "surr.h"
#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "matvec.h"
#include "input.h"
#include "gblvar.h"
#include "stdpso.h"
#include "benchmark.h"
#include "penalty.h"

using namespace std;

// -------------------------------------------------------------------------
// Static Variables:
//
eSigmaType cSURR :: Sigtype = NAKAYAMA;
//
// -------------------------------------------------------------------------
// Set read functions labels:
//
static bool ReadFuncRegister[] =
{
  CtrlMap( ).Insert("SIGMA.TYPE", cSURR::ReadSigtype)
};

// ============================ Normalize ==================================

void sSampData :: Normalize(cVector &x) const
{
  for (int i = 0; i < NumVar; i++) 
    x[i] = (x[i]-Xlow[i])/(Xupp[i] - Xlow[i]);
}

// ============================ UndoNormalization ==========================

void sSampData :: UndoNormalization(cVector &x) const
{
  for (int i = 0; i < NumVar; i++) 
    x[i] = (x[i]*(Xupp[i] - Xlow[i]))+Xlow[i];
}

// ============================ EvalYbounds ================================

void sSampData :: EvalYbounds( )
{
  Ylow.Resize(NumOut);
  Yupp.Resize(NumOut);
  Ylow = SampleY[0];
  Yupp = SampleY[0];
  for (int i = 1; i < NumSample; i++)
  {
    for (int j = 0; j < NumOut; j++)
    {
      Ylow[j] = min(Ylow[j], SampleY[i][j]);
      Yupp[j] = max(Yupp[j], SampleY[i][j]);
    }
  }
}

// ================================= operator>> ==================================

istream& operator>>(istream &in, sSampData &sd)
{
  // Read sample data.

  if (!(in >> sd.NumSample) || !(in >> sd.NumVar)) 
  {
    cout << "Error in the input of sample data." << endl;
    exit(0);
  }

  // Read normalization data.

  bool normalize;
  if (!(in >> normalize)) 
  {
    cout << "Error in the input of normalization flag." << endl;
    exit(0);
  }


  sd.Xlow.Resize(sd.NumVar);
  sd.Xupp.Resize(sd.NumVar);
  if (normalize)
  {
    for (int i = 0; i < sd.NumVar; i++) in >> sd.Xlow[i];
    for (int i = 0; i < sd.NumVar; i++) in >> sd.Xupp[i];
  }

  // Read sampling points and normalize to range [0, 1].

  cVector xn(sd.NumVar);
  for (int i = 0; i < sd.NumSample; i++)
  {
    for (int j = 0; j < sd.NumVar; j++)
    {
      if (!(in >> xn[j])) 
      {
	cout << "Error in the input of sample" << i  << "variable " << j << "." << endl;
	exit(0);
      }
      if (normalize)
	xn[j] = (xn[j] - sd.Xlow[j])/(sd.Xupp[j] - sd.Xlow[j]);
    }
    sd.SampleX.push_back(xn);
  }

  // Read sample values (i.o. outputs).

  if (!(in >> sd.NumOut)) 
  {
    cout << "Error in the input sample data." << endl;
    exit(0);
  }

  cVector yout(sd.NumOut);
  for (int i = 0; i < sd.NumSample; i++)
  {
    for (int j = 0; j < sd.NumOut; j++) in >> yout[j];
    sd.SampleY.push_back(yout);
  }

  // Read validation samples
  bool flagvs;  // Validation sample flag.

  in >> flagvs;

  if (flagvs)              // 1 means there are validations points
  {
      in >> sd.NumValSample;              // read the number of validations points

      // Read validation points and normalize to range [0, 1].

      cVector xtemp(sd.NumValSample);
      for (int i = 0; i < sd.NumValSample; i++)
      {
        for (int j = 0; j < sd.NumVar; j++)
        {
          in >> xtemp[j];
          if (normalize)
            xtemp[j] = (xtemp[j] - sd.Xlow[j])/(sd.Xupp[j] - sd.Xlow[j]);
        }
        sd.ValSampleY.push_back(xtemp);
      }

      // Read output values for validation points

      cVector youtemp(sd.NumOut);
      for (int i = 0; i < sd.NumValSample; i++)
      {
    for (int j = 0; j < sd.NumOut; j++) in >> youtemp[j];

	sd.ValSampleY.push_back(youtemp);
      }
  }
  else
  {
    cout << "No validation data has been provided." << endl;
  }
}

// -------------------------------------------------------------------------
// Public methods:
//
// ============================= ReadSigtype =============================

void cSURR :: ReadSigtype(void)
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
      Sigtype = KITAYAMA;
  }
  else if(string(signame) == "UNIFORMKITAYAMA" || string(signame) == "UniformKitayama"
	  || string(signame) == "UKIT")
  {
      Sigtype = KITAYAMAUNIFORM;
  }
  else if(string(signame) == "ADAPTIVESCALING" || string(signame) == "AdaptiveScaling"
	  || string(signame) == "ASKIT")
  {
      Sigtype = ADAPTIVE_SCALING;
  }
  else if(string(signame) == "LEAVEONEOUTCROSSVALIDATION"
	  || string(signame) == "LeaveOneOutCrossValidation" || string(signame) == "LOOCV")
  {
      Sigtype = LOOCV;
  }
  else if(string(signame) == "KFOLDCROSSVALIDATION" || string(signame) == "KFOLDCV"
	  || string(signame) == "KFCV")
  {
      Sigtype = kFOLDCV;
  }
  else if(string(signame) == "NAKAYAMA" || string(signame) == "Nakayama"
	  || string(signame) == "NAK")
  {
      Sigtype = NAKAYAMA;
  }
  else
  {
  cout << "Unknown sigma definition method: " << signame << endl;
  exit(0);
  }
}

// ================================= cSURR ==================================

cSURR :: cSURR(void)
{
}

// ================================= ~cSURR =================================

cSURR :: ~cSURR()
{
}


// ================================ GetData ================================

void cSURR :: GetData(int &ns, int &nv, int &no, int &nvs)
{
  ns  = sdata.NumSample;
  nv  = sdata.NumVar;
  no  = sdata.NumOut;
  nvs = sdata.NumValSample;
}

// ============================== IsInSample ===============================

bool cSURR :: IsInSample(cVector &x, double tol)
{
  cVector d(sdata.NumVar);
  for (int i = 0; i < sdata.NumSample; i++)
  {
    d = x - sdata.SampleX[i];
    if (d.Length( ) <= tol) return(true);
  }
  return(false);
}

// =========================== Virtual Methods =============================

// =============================== Evaluate ================================

void cSURR :: Evaluate(cVector &x, cVector &y, vector<cVector> *sampx)
{
  cout << "Error: Sur evaluate not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

double cSURR :: EvalExpImp(cVector &x, double ybest)
{
  cout << "Error: Evaluation of Expected Improvement not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

double cSURR :: EvalProbImp(cVector &x, double ybest)
{
  cout << "Error: Evaluation of Probability of Improvement not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

double cSURR :: EvalLCB(cVector &x)
{
  cout << "Error: Evaluation of Lower Confidence Bound not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

double cSURR :: Eval(cVector, int)
{
  cout << "Error: Evaluation of Likelihood not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

double cSURR :: EvalConstraintPF(cVector &, int, double)
{
  cout << "Error: Evaluation of a constraint's Probability of Feasibility not defined\n";
  exit(0);
}

// =============================== Evaluate ================================

double cSURR :: SSqrSur(cVector x, int out)
{
  cout << "Error: Evaluation of s2 not defined\n";
  exit(0);
}

// ======================================================= End of file =====
