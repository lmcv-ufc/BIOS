// -------------------------------------------------------------------------
// fgmplt.cpp - Implementation of the FG Plate problem class.
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
// Created:      27-Aug-2019    Marina Alves Maia
//
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string>
#include <vector>

#ifdef _OMP_
#include "omp.h"
#endif

#include "problem.h"
#include "fgm.h"
#include "fgmplt.h"
#include "input.h"
#include "utl.h"
#include "mat.h"
#include "vec.h"
#include "matvec.h"
#include "sysmat.h"
#include "material.h"
#include "optalg.h"
#include "rbf.h"
#include "gbldef.h"
#include "gblvar.h"

using namespace std;

// -------------------------------------------------------------------------
// Register problems on the problem factory:
//
static const bool registeredProb[] =
{
  cProblemFactory :: Register("SquarePlateBuckFGM"       , MakeProb<cSquarePlateBuckFGM>       ,".fgm"),
  cProblemFactory :: Register("SquarePlateFreqFGM"       , MakeProb<cSquarePlateFreqFGM>       ,".fgm"),
  cProblemFactory :: Register("SquarePlateHoleBuckFGM"   , MakeProb<cSquarePlateHoleBuckFGM>   ,".fgm"),
  cProblemFactory :: Register("SquarePlateFreqFrancoFGM" , MakeProb<cSquarePlateFreqFrancoFGM> ,".fgm"),
  cProblemFactory :: Register("ScoordelisDispFGM"        , MakeProb<cScoordelisFGM>            ,".fgm"),
  cProblemFactory :: Register("CircularPlateFreqFGM"     , MakeProb<cCircularPlateFreqFGM>     ,".fgm")
};

// -------------------------------------------------------------------------
// Public methods:
//


// ============================= cFGMPlate ================================

cFGMPlate :: cFGMPlate(void)
{
}

// ============================= ~cFGMPlate ================================

cFGMPlate :: ~cFGMPlate(void)
{
}

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateBuckFGM :: cSquarePlateBuckFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation

  // Ceramic volume percentage < Cmax
  // Cmax = 35%

  double vcratio;
  int numcp;

  if ((NumVar)%2 == 0){
      numcp = (NumVar)*2;
  }
  else{
      numcp = 2*NumVar - 1;
  }

  cVector Vcp(numcp);
  for (int i = 0; i < NumVar; i++){
      Vcp[i] = x[i];
      Vcp[numcp - i - 1] = x[i];
  }

  EvalVolumeRatio(Vcp, vcratio);

  c[0] = vcratio - 0.35;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateBuckFGM :: Analysis(cVector x, double &lbdb)
{
    // Evaluate volume fraction at gauss points according to a given distribution

    cVector Vcpg;

    int numcp;

    if ((NumVar)%2 == 0){
        numcp = (NumVar)*2;
    }
    else{
        numcp = 2*(NumVar) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 0; i < NumVar; i++){
        Vcp[i] = x[i];
        Vcp[numcp - i - 1] = x[i];
    }

    double thk = 1.0;

  int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string thread_number = thread.str();
  string cmd  = "del SqrPltBuck" + thread_number + ".dat";
  string cmd2 = "del SqrPltBuck" + thread_number + ".pos";
  string cmd3 = "rm SqrPltBuck" + thread_number + ".dat";
  string cmd4 = "rm SqrPltBuck" + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing SqrPltBuck.dat and plate.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing SqrPltBuck.dat and plate.pos files.\n";
#endif

  string begname = "datbegSqrPltBuck16x16.dat";
  string endname = "datendSqrPltBuck16x16.dat";

  string datname = "SqrPltBuck" + thread_number + ".dat";
  string posname = "SqrPltBuck" + thread_number + ".pos";

#ifdef _WIN32
  cmd = "type " + begname + " >> " + datname;
#else
  cmd = "cat " + begname + " >> " + datname;
#endif

  int status1 = system(cmd.c_str( ));
  int status2;

  if (status1)
  {
     cout << "Error in the copy of datbeg file.";
     lbdb = 0.0;
     return;
  }

  fstream dat;

  dat.open(datname.c_str( ));

  if (!dat.is_open( ))
  {
     cout << "Error opening the dat file for plate analysis." << endl;
     exit(0);
  }

  dat.seekp(0,ofstream::end);

  dat << "%SECTION.FGM.SHELL" << endl;
  dat << "1" << endl;
  dat << "1    1    " << thk << "    10    2    " << numcp << "    ";

  for (int i = 0; i < numcp; i++) dat << Vcp[i] << "    ";
  dat << endl;

  dat.close( );

#ifdef _WIN32
  cmd = "type " + endname + " >> " + datname;
#else
  cmd = "cat " + endname + " >> " + datname;
#endif

  status1 = system(cmd.c_str( ));

  if (status1)
  {
     cout << "Error in the copy of datend file.";
     lbdb = 0.0;
     exit(0);
  }

  // Run the analysis with FAST.

#ifdef _WIN32
  cmd = "fast.exe SqrPltBuck" + thread_number + " -silent";
#else
  cmd = "./fast SqrPltBuck" + thread_number + " -silent";
#endif

  status2 = system(cmd.c_str( ));

  if (status2)
  {
     cout << "Error in the analysis with fast.";
    #ifdef _WIN32
      cmd = "fast.exe SqrPltBuck" + thread_number + " -silent";
    #else
      cmd = "./fast SqrPltBuck" + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
          lbdb = 0.0;
      }
  }

  if (!status2)
  {
  // Open the pos file.

  ifstream pos;

  pos.open(posname.c_str( ));
  if (!pos.is_open( ))
  {
     cout << "Error opening the pos file for plate analysis." << endl;
     exit(0);
  }

  // Find buckling load factor

  string label;
  double buckfactor = 0;
  int mode;

  while (pos >> label)
  {

      if (label == "%RESULT.CASE.STEP.BUCKLING.FACTOR")
      {
          pos >> mode;
          pos >> buckfactor;
      }
   }

   if (buckfactor == 0)
   {
      cout << "Convergence not achieved in infill: " << endl;
   }

   // Push back the new targets Ybuck and Ystren
   lbdb = buckfactor;
  }
}

// ============================ Evaluate ==============================

void cSquarePlateBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    double vcratio;

    int numcp;

    if ((NumVar)%2 == 0){
        numcp = (NumVar)*2;
    }
    else{
        numcp = 2*NumVar - 1;
    }

    cVector Vcp(numcp);
    for (int i = 0; i < NumVar; i++){
        Vcp[i] = x[i];
        Vcp[numcp - i - 1] = x[i];
    }

    EvalVolumeRatio(Vcp, vcratio);

  // Single constraint evaluation.
    if (index == 0){
        c = vcratio - 0.35;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateBuckFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateFreqFGM :: cSquarePlateFreqFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateFreqFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj);  // Linearized Frequency Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation

  // Ceramic volume percentage < Cmax
  // Cmax = 50%

  double vcratio;
  int numcp;

  if ((NumVar)%2 == 0){
      numcp = (NumVar)*2;
  }
  else{
      numcp = 2*NumVar - 1;
  }

  cVector Vcp(numcp);
  for (int i = 0; i < NumVar; i++){
      Vcp[i] = x[i];
      Vcp[numcp - i - 1] = x[i];
  }

  EvalVolumeRatio(Vcp, vcratio);

  c[0] = vcratio - 0.50;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateFreqFGM :: Analysis(cVector x, double &lbdb)
{
    // Evaluate volume fraction at gauss points according to a given distribution

    cVector Vcpg;

    int numcp;

    if ((NumVar)%2 == 0){
        numcp = (NumVar)*2;
    }
    else{
        numcp = 2*(NumVar) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 0; i < NumVar; i++){
        Vcp[i] = x[i];
        Vcp[numcp - i - 1] = x[i];
    }

    double thk = 1.0;

  int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string thread_number = thread.str();
  string cmd = "del SqrPltFreq" + thread_number + ".dat";
  string cmd2 =  "del SqrPltFreq" + thread_number + ".pos";
  string cmd3 = "rm SqrPltFreq" + thread_number + ".dat";
  string cmd4 = "rm SqrPltFreq" + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing SqrPltFreq.dat and plate.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing SqrPltFreq.dat and plate.pos files.\n";
#endif

  string begname = "datbegSqrPltFreq16x16.dat";
  string endname = "datendSqrPltFreq16x16.dat";

  string datname = "SqrPltFreq" + thread_number + ".dat";
  string posname = "SqrPltFreq" + thread_number + ".pos";

#ifdef _WIN32
  cmd = "type " + begname + " >> " + datname;
#else
  cmd = "cat " + begname + " >> " + datname;
#endif

  int status1 = system(cmd.c_str( ));
  int status2;

  if (status1)
  {
     cout << "Error in the copy of datbeg file.";
     lbdb = 0.0;
     return;
  }

  fstream dat;

  dat.open(datname.c_str( ));

  if (!dat.is_open( ))
  {
     cout << "Error opening the dat file for plate analysis." << endl;
     exit(0);
  }

  dat.seekp(0,ofstream::end);

  dat << "%SECTION.FGM.SHELL" << endl;
  dat << "1" << endl;
  dat << "1    1    " << thk << "    10    2    " << numcp << "    ";

  for (int i = 0; i < numcp; i++) dat << Vcp[i] << "    ";
  dat << endl;

  dat.close( );

#ifdef _WIN32
  cmd = "type " + endname + " >> " + datname;
#else
  cmd = "cat " + endname + " >> " + datname;
#endif

  status1 = system(cmd.c_str( ));

  if (status1)
  {
     cout << "Error in the copy of datend file.";
     lbdb = 0.0;
     exit(0);
  }

  // Run the analysis with FAST.

#ifdef _WIN32
  cmd = "fast.exe SqrPltFreq" + thread_number + " -silent";
#else
  cmd = "./fast SqrPltFreq" + thread_number + " -silent";
#endif

  status2 = system(cmd.c_str( ));

  if (status2)
  {
     cout << "Error in the analysis with fast.";
    #ifdef _WIN32
      cmd = "fast.exe SqrPltFreq" + thread_number + " -silent";
    #else
      cmd = "./fast SqrPltFreq" + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
          lbdb = 0.0;
      }
  }

  // Open the pos file.

  ifstream pos;

  pos.open(posname.c_str( ));

  if (!pos.is_open( ))
  {
     cout << "Error opening the pos file for plate analysis." << endl;
     lbdb = 0;
     return;
  }

  // Find natural frequency

  string label;
  double vibfactor = 0;
  int mode;

  while (pos >> label)
  {

      if (label == "%RESULT.CASE.STEP.NATURAL.FREQUENCY")
      {
          pos >> mode;
          pos >> vibfactor;
      }
   }

  //cout << "                       \n buckfactor " << buckfactor << endl;
   if (vibfactor == 0)
   {
      cout << "Convergence not achieved in infill: " << endl;
//      exit(0);
   }

   // Push back the new targets Ybuck and Ystren
   lbdb = vibfactor;
}

// ============================ Evaluate ==============================

void cSquarePlateFreqFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    double vcratio;

    int numcp;

    if ((NumVar)%2 == 0){
        numcp = (NumVar)*2;
    }
    else{
        numcp = 2*NumVar - 1;
    }

    cVector Vcp(numcp);
    for (int i = 0; i < NumVar; i++){
        Vcp[i] = x[i];
        Vcp[numcp - i - 1] = x[i];
    }

    EvalVolumeRatio(Vcp, vcratio);

  // Single constraint evaluation.
    if (index == 0){
        c = vcratio - 0.50;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateFreqFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}

// ========================== cPltHoleFGM ===========================

cSquarePlateHoleBuckFGM :: cSquarePlateHoleBuckFGM(void)
{
  NumConstr = 2;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateHoleBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Evaluating the objective function

    double fobj;
    Analysis(x, fobj);
    fobjs[0] = -fobj;

    // Evaluating constraints

    double vcratio, rho, mass;
    int numcp;

    if ((NumVar - 1)%2 == 0){
        numcp = (NumVar - 1)*2;
    }
    else{
        numcp = 2*(NumVar - 1) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 1; i < NumVar; i++){
        Vcp[i - 1] = x[i];
        Vcp[numcp - i ] = x[i];
    }

    EvalVolumeRatio(Vcp, vcratio);

    rho = 0;
    EvalDens(Vcp, rho);
    double a = 0.72;
    double r = 0.072;
    mass = rho*(a*a - PI*r*r)*x[0];

    c[0] = vcratio - 0.50; // Ceramic volume fraction
    c[1] = mass - 100;     // Mass

    /*// Derivative assessment

    int nstep = 100;
    cVector xval(nstep);
    xval[0] = 0.0;
    for (int i = 1; i < nstep; i++){
        xval[i] = xval[i - 1] + 1.0/((double) nstep);
    }

    // tClock cputime = clock();

    cVector BSder(nstep);
    BSplineDer(Vcp, nstep, xval, BSder);

    // cputime = clock() - cputime;

    // cout << "time = " << (cputime/CLOCKS_PER_SEC) << endl;

    cVector ABSder(nstep);

    for (int i = 0; i < nstep; i++)
    {
        ABSder[i] = abs(BSder[i]);
    }

    double MaxDer;
    MaxDer = ABSder.Max();

    c[2] = MaxDer - tan(75*PI/180);*/
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateHoleBuckFGM :: Analysis(cVector x, double &lbdb)
{
    // Evaluate volume fraction at gauss points according to a given distribution

    cVector Vcpg;

    int numcp;

    double thk = x[0];




    if ((NumVar - 1)%2 == 0){
        numcp = (NumVar - 1)*2;
    }
    else{
        numcp = 2*(NumVar - 1) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 0; i < NumVar - 1; i++){
        Vcp[i] = x[i + 1];
        Vcp[numcp - i - 1] = x[i + 1];
    }

  int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string thread_number = thread.str();
  string cmd = "del SqrPltHoleBuck" + thread_number + ".dat";
  string cmd2 =  "del SqrPltHoleBuck" + thread_number + ".pos";
  string cmd3 = "rm SqrPltHoleBuck" + thread_number + ".dat";
  string cmd4 = "rm SqrPltHoleBuck" + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing SqrPltHoleBuck.dat and plate.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing SqrPltHoleBuck.dat and plate.pos files.\n";
#endif

  string begname = "datbegSqrPltHoleBuck.dat";
  string endname = "datendSqrPltHoleBuck.dat";

  string datname = "SqrPltHoleBuck" + thread_number + ".dat";
  string posname = "SqrPltHoleBuck" + thread_number + ".pos";

#ifdef _WIN32
  cmd = "type " + begname + " >> " + datname;
#else
  cmd = "cat " + begname + " >> " + datname;
#endif

  int status1 = system(cmd.c_str( ));
  int status2;

  if (status1)
  {
     cout << "Error in the copy of datbeg file.";
     lbdb = 0.0;
     return;
  }

  fstream dat;

  dat.open(datname.c_str( ));

  if (!dat.is_open( ))
  {
     cout << "Error opening the dat file for plate analysis." << endl;
     exit(0);
  }

  dat.seekp(0,ofstream::end);

  dat << "%SECTION.FGM.SHELL" << endl;
  dat << "1" << endl;
  dat << "1    1    " << thk << "    10    2    " << numcp << "    ";

  for (int i = 0; i < numcp; i++) dat << Vcp[i] << "    ";
  dat << endl;

  dat.close( );

#ifdef _WIN32
  cmd = "type " + endname + " >> " + datname;
#else
  cmd = "cat " + endname + " >> " + datname;
#endif

  status1 = system(cmd.c_str( ));

  if (status1)
  {
     cout << "Error in the copy of datend file.";
     lbdb = 0.0;
     exit(0);
  }

  // Run the analysis with FAST.

#ifdef _WIN32
  cmd = "fast.exe SqrPltHoleBuck" + thread_number + " -silent";
#else
  cmd = "./fast SqrPltHoleBuck" + thread_number + " -silent";
#endif

  status2 = system(cmd.c_str( ));

  if (status2)
  {
     cout << "Error in the analysis with fast.";
    #ifdef _WIN32
      cmd = "fast.exe SqrPltHoleBuck" + thread_number + " -silent";
    #else
      cmd = "./fast SqrPltHoleBuck" + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
          lbdb = 0.0;
      }
  }

  if (!status2)
  {
  // Open the pos file.

  ifstream pos;

  pos.open(posname.c_str( ));
  if (!pos.is_open( ))
  {
     cout << "Error opening the pos file for plate analysis." << endl;
     exit(0);
  }

  // Find buckling load factor

  string label;
  double buckfactor = 0;
  int mode;

  while (pos >> label)
  {

      if (label == "%RESULT.CASE.STEP.BUCKLING.FACTOR")
      {
          pos >> mode;
          pos >> buckfactor;
      }
   }

   if (buckfactor == 0)
   {
      cout << "Convergence not achieved in infill: " << endl;
   }

   // Push back the new targets Ybuck and Ystren
   lbdb = buckfactor;
  }
}

// ============================ Evaluate ==============================

void cSquarePlateHoleBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    double vcratio;

    int numcp;

    if ((NumVar - 1)%2 == 0){
        numcp = (NumVar - 1)*2;
    }
    else{
        numcp = 2*(NumVar - 1) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 1; i < NumVar; i++){
        Vcp[i - 1] = x[i];
        Vcp[numcp - i ] = x[i];
    }

    EvalVolumeRatio(Vcp, vcratio);

    // Exact constraint evaluation.

    if (index == 0){
        c = vcratio - 0.50;
    }
    else if (index == 1){
        double rho, mass;

        rho = 0;
        EvalDens(Vcp, rho);
        double a = 0.72;
        double r = 0.072;
        mass = rho*(a*a - PI*r*r)*x[0];

        c = mass - 100;
    }
    /*else if (index == 2){
        // Derivative assessment

        int nstep = 100;
        cVector xval(nstep);
        xval[0] = 0.0;
        for (int i = 1; i < nstep; i++){
            xval[i] = xval[i - 1] + 1/((double) nstep);
        }

        // tClock cputime = clock();

        cVector BSder(nstep);
        BSplineDer(Vcp, nstep, xval, BSder);

        // cputime = clock() - cputime;

        // cout << "time = " << (cputime/CLOCKS_PER_SEC) << endl;

        cVector ABSder(nstep);

        for (int i = 0; i < nstep; i++)
        {
            ABSder[i] = abs(BSder[i]);
        }

        double MaxDer;
        MaxDer = ABSder.Max();

        c = MaxDer - tan(75*PI/180);
    }*/
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateHoleBuckFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
  approxc[1] = 0;
}

// =================== cSquarePlateFreqFrancoFGM ======================

cSquarePlateFreqFrancoFGM :: cSquarePlateFreqFrancoFGM(void)
{
  NumConstr = 2;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateFreqFrancoFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Fundamental Frequency

  double freq;

  Analysis(x, freq);

  fobjs[0] = -freq;

  c[0] = 1.00 - freq/(3000.00);
  c[1] = freq/(8000.00) - 1.00;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateFreqFrancoFGM :: Analysis(cVector x, double &lbdb)
{
  double thk    = x[0];
  double pindex = x[1];

  int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string thread_number = thread.str();
  string cmd = "del SqrPltFrancoFreq" + thread_number + ".dat";
  string cmd2 =  "del SqrPltFrancoFreq" + thread_number + ".pos";
  string cmd3 = "rm SqrPltFrancoFreq" + thread_number + ".dat";
  string cmd4 = "rm SqrPltFrancoFreq" + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing SqrPltFrancoFreq.dat and SqrPltFrancoFreq.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing SqrPltFrancoFreq.dat and SqrPltFrancoFreq.pos files.\n";
#endif

  string begname = "datbegSqrPltFrancoFreq16x16.dat";
  string endname = "datendSqrPltFrancoFreq16x16.dat";

  string datname = "SqrPltFrancoFreq" + thread_number + ".dat";
  string posname = "SqrPltFrancoFreq" + thread_number + ".pos";

    #ifdef _WIN32
      cmd = "type " + begname + " >> " + datname;
    #else
      cmd = "cat " + begname + " >> " + datname;
    #endif

  int status = system(cmd.c_str( ));

  if (status)
  {
     cout << "Error in the copy of datbeg file.";
     lbdb = 0;
     return;
  }

  fstream dat;

  dat.open(datname.c_str( ));

  if (!dat.is_open( ))
  {
     cout << "Error opening the dat file for plate analysis." << endl;
     lbdb = 0;
     return;
  }

  dat.seekp(0,ofstream::end);

  dat << "%SECTION.FGM.SHELL" << endl;
  dat << "1" << endl;
  dat << "1    1    " << thk << "    10    1    1    " << pindex << endl;

  dat.close( );

#ifdef _WIN32
  cmd = "type " + endname + " >> " + datname;
#else
  cmd = "cat " + endname + " >> " + datname;
#endif

  status = system(cmd.c_str( ));

  if (status)
  {
     cout << "Error in the copy of datbeg file.";
     lbdb = 0.0;
     return;
  }

  // Run the analysis with FAST.

#ifdef _WIN32
  cmd = "fast.exe SqrPltFrancoFreq" + thread_number + " -silent";
#else
  cmd = "./fast SqrPltFrancoFreq" + thread_number + " -silent";
#endif

  status = system(cmd.c_str( ));

  if (status)
  {
     cout << "Error in the analysis with fast.";
     lbdb = 0;
     return;
  }

  // Open the pos file.

  ifstream pos;

  pos.open(posname.c_str( ));

  if (!pos.is_open( ))
  {
     cout << "Error opening the pos file for plate analysis." << endl;
     exit(0);
  }

  // Find fundamental frequency

  string label;
  double vibfactor = 0;
  int mode;

  while (pos >> label)
  {

      if (label == "%RESULT.CASE.STEP.NATURAL.FREQUENCY")
      {
          pos >> mode;
          pos >> vibfactor;
      }
   }

   if (vibfactor == 0)
   {
      cout << "Convergence not achieved in infill: " << endl;
   }

   // Push back the new targets Ybuck and Ystren
   lbdb = (vibfactor);
}

// ========================= GetApproxConstr ==========================

void cSquarePlateFreqFrancoFGM :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
  approxc[1] = 1;
}

// ========================== cScoordelisFGM ===========================

cScoordelisFGM :: cScoordelisFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cScoordelisFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  double thk    = x[0];
  double pindex = x[1];

  // Minimize mass

  double dens, mass;

  EvalDens(pindex, dens);

  double R = 2.54;

  double frac = 0.20/(2.0*PI);

  mass = frac*PI*((R+thk/2.0)*(R+thk/2.0)-(R-thk/2.0)*(R-thk/2.0))*dens*0.508;  // h*L*w*density

  // Evaluate displacements and stresses

 double disp, stress, yield;

  AnalysisDispStress(x, disp, stress, yield);

  // Set objective and constraint functions

  fobjs[0] = mass;

  double dispmax = 0.004;
  c[0] = abs(disp)/dispmax - 1.00;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cScoordelisFGM :: AnalysisDispStress(cVector x, double &disp,
                                            double &stress, double &yield)
{
    double thk    = x[0];
    double pindex = x[1];

    int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string thread_number = thread.str();
  string cmd = "del ScoordelisFGM" + thread_number + ".dat";
  string cmd2 =  "del ScoordelisFGM" + thread_number + ".pos";
  string cmd3 = "rm ScoordelisFGM" + thread_number + ".dat";
  string cmd4 = "rm ScoordelisFGM" + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing ScoordelisFGM.dat and ScoordelisFGM.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing ScoordelisFGM.dat and ScoordelisFGM.pos files.\n";
#endif

  string begname = "datbegScoordelisFGM.dat";
  string endname = "datendScoordelisFGM.dat";

  string datname = "ScoordelisFGM" + thread_number + ".dat";
  string posname = "ScoordelisFGM" + thread_number + ".pos";

#ifdef _WIN32
  cmd = "type " + begname + " >> " + datname;
#else
  cmd = "cat " + begname + " >> " + datname;
#endif

  int status = system(cmd.c_str( ));

  if (status)
  {
     cout << "Error in the copy of datbeg file.";
     disp = 1e6;
     stress = 1e20;
     return;
//     exit(EXIT_FAILURE);
  }

  fstream dat;

  dat.open(datname.c_str( ));

  if (!dat.is_open( ))
  {
     cout << "Error opening the dat file for plate analysis." << endl;
     exit(0);
  }

  dat.seekp(0,ofstream::end);

  dat << "%SECTION.FGM.SHELL" << endl;
  dat << "1" << endl;
  dat << "1    1    " << thk << "    10    1    1    " << pindex << endl;

  dat.close( );

#ifdef _WIN32
  cmd = "type " + endname + " >> " + datname;
#else
  cmd = "cat " + endname + " >> " + datname;
#endif

  status = system(cmd.c_str( ));

  if (status)
  {
     cout << "Error in the copy of datbeg file.";
     exit(EXIT_FAILURE);
  }

  // Run the analysis with FAST.

#ifdef _WIN32
  cmd = "fast.exe ScoordelisFGM" + thread_number + " -silent";
#else
  cmd = "./fast ScoordelisFGM" + thread_number + " -silent";
#endif

  status = system(cmd.c_str( ));

  if (status)
  {
     cout << "Error in the analysis with fast.";
  //   exit(EXIT_FAILURE);
     disp = -1e6;
     return;
  }

  // Open the pos file.

  ifstream pos;

  pos.open(posname.c_str( ));

  if (!pos.is_open( ))
  {
     cout << "Error opening the pos file for plate analysis." << endl;
     exit(0);
  }

  // Find buckling load factor

  string label;
  cVector genstress(8);
  cVector genstrain(8);
  cVector force(8);
  force.Zero();
  cVector desl(6);
  desl.Zero();
  double wcenter = -1.0e6;
  int step = 1;
  int nn, elmid;
  int numelm;

  while (pos >> label)
  {
      if (label == "%RESULT.CASE.STEP")
      {
          pos >> step;
          if (step == 6)
          {
              pos >> label;
              pos >> label;
              pos >> label;

              if (label == "%RESULT.CASE.STEP.NODAL.DISPLACEMENT")
              {
                 pos >> nn; // read number of nodes
                 pos >> label;

                     for (int j = 0; j < nn; j++)
                     {
                         pos >> elmid;
                         pos >> desl[0];
                         pos >> desl[1];
                         pos >> desl[2];
                         pos >> desl[3];
                         pos >> desl[4];
                         pos >> desl[5];

                         if (elmid == 181)
                         {
                             wcenter = desl[2];
                         }
                     }
              }
          }
      }
  }

  disp = wcenter;
}

// ============================ Evaluate ==============================

void cScoordelisFGM :: EvalExactFobj(cVector& x, double &fobj)
{
    double thk    = x[0];
    double pindex = x[1];

    // Minimize mass

    double dens, mass;

    EvalDens(pindex, dens);

    double R = 2.54;

    double frac = 0.20/(2.0*PI);

    mass = frac*PI*((R+thk/2.0)*(R+thk/2.0)-(R-thk/2.0)*(R-thk/2.0))*dens*0.508;  // h*L*w*density

    fobj = mass;
}

// ========================= GetApproxConstr ==========================

void cScoordelisFGM :: GetApproxObj(bool* approxobj)
{
  approxobj[0] = 0;
}

// ========================= GetApproxConstr ==========================

void cScoordelisFGM :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 1;
}

// ======================= cCircularPlateFreqFGM ===========================

cCircularPlateFreqFGM :: cCircularPlateFreqFGM(void)
{
  NumConstr = 2;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cCircularPlateFreqFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Fundamental frequency

    double fobj;
    Analysis(x, fobj);
    fobjs[0] = -fobj;

    // Constraint evaluation

    double vcratio, rho, mass, area;
    int numcp;

    if ((NumVar - 1)%2 == 0){
        numcp = (NumVar - 1)*2;
    }
    else{
        numcp = 2*(NumVar - 1) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 1; i < NumVar; i++){
        Vcp[i - 1] = x[i];
        Vcp[numcp - i] = x[i];
    }

    EvalVolumeRatio(Vcp, vcratio);

    rho = 0;
    EvalDens(Vcp, rho);
    area = PI*0.5*0.5;
    mass = rho*area*x[0];

    c[0] = (mass - 100)/100;
    c[1] = (vcratio - 0.35)/0.35;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cCircularPlateFreqFGM :: Analysis(cVector x, double &lbdb)
{
    // Evaluate volume fraction at gauss points according to a given distribution

    cVector Vcpg;

    int numcp;

    double thk = x[0];

    if ((NumVar - 1)%2 == 0){
        numcp = (NumVar - 1)*2;
    }
    else{
        numcp = 2*(NumVar - 1) - 1;
    }

    cVector Vcp(numcp);
    for (int i = 0; i < (NumVar - 1); i++){
        Vcp[i] = x[i + 1];
        Vcp[numcp - i - 1] = x[i + 1];
    }

  int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string thread_number = thread.str();
  string cmd = "del CircularPltFreq" + thread_number + ".dat";
  string cmd2 =  "del CircularPltFreq" + thread_number + ".pos";
  string cmd3 = "rm CircularPltFreq" + thread_number + ".dat";
  string cmd4 = "rm CircularPltFreq" + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing CircularPltFreq.dat and plate.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing CircularPltFreq.dat and plate.pos files.\n";
#endif

  string begname = "datbegCircularPltFreq.dat";
  string endname = "datendCircularPltFreq.dat";

  string datname = "CircularPltFreq" + thread_number + ".dat";
  string posname = "CircularPltFreq" + thread_number + ".pos";

#ifdef _WIN32
  cmd = "type " + begname + " >> " + datname;
#else
  cmd = "cat " + begname + " >> " + datname;
#endif

  int status1 = system(cmd.c_str( ));
  int status2;

  if (status1)
  {
     cout << "Error in the copy of datbeg file.";
     lbdb = 0.0;
     return;
  }

  fstream dat;

  dat.open(datname.c_str( ));

  if (!dat.is_open( ))
  {
     cout << "Error opening the dat file for plate analysis." << endl;
     exit(0);
  }

  dat.seekp(0,ofstream::end);

  dat << "%SECTION.FGM.SHELL" << endl;
  dat << "1" << endl;
  dat << "1    1    " << thk << "    10    2    " << numcp << "    ";

  for (int i = 0; i < numcp; i++) dat << Vcp[i] << "    ";
  dat << endl;

  dat.close( );

#ifdef _WIN32
  cmd = "type " + endname + " >> " + datname;
#else
  cmd = "cat " + endname + " >> " + datname;
#endif

  status1 = system(cmd.c_str( ));

  if (status1)
  {
     cout << "Error in the copy of datend file.";
     lbdb = 0.0;
     exit(0);
  }

  // Run the analysis with FAST.

#ifdef _WIN32
  cmd = "fast.exe CircularPltFreq" + thread_number + " -silent";
#else
  cmd = "./fast CircularPltFreq" + thread_number + " -silent";
#endif

  status2 = system(cmd.c_str( ));

  if (status2)
  {
     cout << "Error in the analysis with fast.";
    #ifdef _WIN32
      cmd = "fast.exe CircularPltFreq" + thread_number + " -silent";
    #else
      cmd = "./fast CircularPltFreq" + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
          lbdb = 0.0;
      }
  }

  // Open the pos file.

  ifstream pos;

  pos.open(posname.c_str( ));

  if (!pos.is_open( ))
  {
     cout << "Error opening the pos file for plate analysis." << endl;
     lbdb = 0;
     return;
  }

  // Find natural frequency

  string label;
  double vibfactor = 0;
  int mode;

  while (pos >> label)
  {

      if (label == "%RESULT.CASE.STEP.NATURAL.FREQUENCY")
      {
          pos >> mode;
          pos >> vibfactor;
      }
   }

   if (vibfactor == 0)
   {
      cout << "Convergence not achieved in infill: " << endl;
   }

   // Push back the new targets Ybuck and Ystren
   lbdb = vibfactor;
}

// ============================ Evaluate ==============================

void cCircularPlateFreqFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    double vcratio;

        int numcp;

        if ((NumVar - 1)%2 == 0){
            numcp = (NumVar - 1)*2;
        }
        else{
            numcp = 2*(NumVar - 1) - 1;
        }

        cVector Vcp(numcp);
        for (int i = 1; i < NumVar; i++){
            Vcp[i - 1] = x[i];
            Vcp[numcp - i ] = x[i];
        }

        EvalVolumeRatio(Vcp, vcratio);

      // Exact constraint evaluation.
        if (index == 0){
            double rho, mass;

            rho = 0;
            EvalDens(Vcp, rho);
            double area = PI*0.5*0.5;
            mass = rho*area*x[0];

            c = (mass - 100)/100;
        }
        else if (index == 1){
            c = (vcratio - 0.35)/0.35;
        }
        else{
            cout << "Definition of an exact constraint missing!";
            exit(0);
        }
}

// ========================= GetApproxConstr ==========================

void cCircularPlateFreqFGM :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 0;
  approxc[1] = 0;
}

// ======================================================= End of file =====
