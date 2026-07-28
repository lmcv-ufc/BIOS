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
#include <chrono>

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
  cProblemFactory :: Register("SquarePlateBuckFGM"       , MakeProb<cSquarePlateBuckFGM>         ,".fgm"),
  cProblemFactory :: Register("SquarePlateFreqFGM"       , MakeProb<cSquarePlateFreqFGM>         ,".fgm"),
  cProblemFactory :: Register("SquarePlateHoleBuckFGM"   , MakeProb<cSquarePlateHoleBuckFGM>     ,".fgm"),
  cProblemFactory :: Register("SquarePlateFreqFrancoFGM" , MakeProb<cSquarePlateFreqFrancoFGM>   ,".fgm"),
  cProblemFactory :: Register("ScoordelisDispFGM"        , MakeProb<cScoordelisFGM>              ,".fgm"),
  cProblemFactory :: Register("CircularPlateFreqFGM"     , MakeProb<cCircularPlateFreqFGM>       ,".fgm"),
  cProblemFactory :: Register("SquarePlateMFBuckFGM"     , MakeProb<cSquarePlateMFBuckFGM>       ,".fgm"),
  cProblemFactory :: Register("SquarePlateMFBuck3DirFGM" , MakeProb<cSquarePlateTriDirMFBuckFGM> ,".fgm"),
  cProblemFactory :: Register("SquarePlateCutOutFGM"     , MakeProb<cSquarePlateCutOutFGM>       ,".fgm"),
  cProblemFactory :: Register("ShallowShellTBuckFGM"     , MakeProb<cShallowShellMFThermBuckFGM> ,".fgm"),
  cProblemFactory :: Register("GuoVSCBuckMF"             , MakeProb<cSquarePlateMFBuckVSC>       ,".fgm")
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

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateMFBuckFGM :: cSquarePlateMFBuckFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
  CostMax = 0.65;
}

// ============================== Evaluate =================================

void cSquarePlateMFBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj, 2);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation

  // Ceramic volume percentage < Cmax
  // Cmax = 35%
  double Cmax = CostMax;

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

  c[0] = vcratio - Cmax;
}

// ============================== Evaluate =================================

void cSquarePlateMFBuckFGM :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    /*x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 0.4;
    x[3] = 0.0;
    x[4] = 0.0;*/

    /*x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    x[3] = 0.0;
    x[4] = 0.0;*/

    /*x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    x[3] = 0.45;
    x[4] = 0.0;*/

  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj, 1);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation

  // Ceramic volume percentage < Cmax
  // Cmax = 35%
  double Cmax = CostMax;

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

  c[0] = vcratio - Cmax;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateMFBuckFGM :: Analysis(cVector x, double &lbdb, int f)
{
    // Define number of gauss points

        int ngauss = 10;

        // Evaluate weights and abscissas of gauss points

        cVector r, w;
        GaussPts1D(ngauss, r, w);

        // Transform abscissas to thickness coordinate

        cVector t(ngauss);

        // for (int i = 0; i < ngauss; i++) t[i]  = (r[i] + 1)/2;  // t = 0 (bottom) and t = 1 (top)
        for (int i = 0; i < ngauss; i++) t[i]  = r[i]/2;  // t = -0.5 (bottom) and t = 0.5 (top)

        // Evaluate volume fraction at gauss points according to a given distribution

        cVector Vcpg;

        /*t.Resize(21);
        t[0] = -0.5;
        for (int i = 1; i < 21; i++) t[i]  = t[i-1] + (1.0/20.0);

        cVector V(5);
        V[0] = 1.0;
        V[1] = 0.84572;
        V[2] = 1.0;
        V[3] = 1.0;
        V[4] = 1.0;

        PiecewiseCubicInterpolation(V, 21, t, Vcpg);

        cout << "Vcpg = " << endl;
        Vcpg.Print();

        exit(0);*/

        /*cVector Vcpg;

        VolumeDist(FGMVolDist, ngauss, t, Vcpg, x[1]);*/

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

        /*cout << "\n\n ============ Matriz A ======= " << endl;
        A.Print();
        cout << "\n ============ Matriz B ======= " << endl;
        B.Print();
        cout << "\n ============ Matriz D ======= " << endl;
        D.Print();
        cout << "\n ============ Matriz G ======= " << endl;
        G.Print();
        cout << "\n ============ Matriz ABDG ======= " << endl;
        ABDG.Print();*/

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd, cmd2, cmd3, cmd4;
      if (f == 1)
      {
          cmd  = "del DoLee2D" + thread_number + ".dat";
          cmd2 = "del DoLee2D" + thread_number + ".pos";
          cmd3 = "rm DoLee2D" + thread_number + ".dat";
          cmd4 = "rm DoLee2D" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          cmd  = "del DoLee3D" + thread_number + ".dat";
          cmd2 = "del DoLee3D" + thread_number + ".pos";
          cmd3 = "rm DoLee3D" + thread_number + ".dat";
          cmd4 = "rm DoLee3D" + thread_number + ".pos";
      }


    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing DoLee.dat and plate.pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing DoLee.dat and plate.pos files.\n";
    #endif

      string begname, endname, datname, posname;

      if (f == 1)
      {
          begname = "datbegDoLeeSS22D.dat";
          endname = "datendDoLeeSS22D.dat";

          datname = "DoLee2D" + thread_number + ".dat";
          posname = "DoLee2D" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          begname = "datbegDoLee3D.dat";
          endname = "datendDoLee3D.dat";

          datname = "DoLee3D" + thread_number + ".dat";
          posname = "DoLee3D" + thread_number + ".pos";
      }

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

         //exit(EXIT_FAILURE);

         lbdb = 0.0;

         cout << "chega aqui" << endl;
         //exit(0);
         return;
      }

      fstream dat;

      dat.open(datname.c_str( ));

      if (!dat.is_open( ))
      {
         cout << "Error opening the dat file for plate analysis." << endl;
         exit(0);
      }

      if (f == 1)
      {

      dat.seekp(0,ofstream::end);

      dat << endl << endl << "%SECTION.FGM.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1    " << thk << "    10    2    " << numcp;

      for (int i = 0; i < numcp; i++)
      {
          dat << "    " << Vcp[i];
      }
      dat << endl;
      }
      else if (f == 2)
      {

      dat.seekp(0,ofstream::end);

      dat << endl << endl << "%SECTION.FGM.3D" << endl;
      dat << "1" << endl;
      dat << "1    1    2    " << numcp+2 << "    0.0    1.0";

      for (int i = 0; i < numcp; i++)
      {
          dat << "    " << Vcp[i];
      }
      dat << endl;
      }

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
      //   exit(EXIT_FAILURE);

         lbdb = 0.0;
         exit(0);
      }

      // Run the analysis with FAST.

    string fid;
    if (f == 1)
    {
        fid = "2D";
    }
    else
    {
        fid = "3D";
    }

    #ifdef _WIN32
      cmd = "fast.exe DoLee" + fid + thread_number + " -silent";
    #else
      cmd = "./fast DoLee" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
         //exit(EXIT_FAILURE);

    #ifdef _WIN32
      cmd = "fast.exe DoLee" + fid + thread_number + " -silent";
    #else
      cmd = "./fast DoLee" + fid + thread_number + " -silent";
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
      //string posname = "platehole.pos";
      //pos.open(posname.c_str());
      if (!pos.is_open( ))
      {
         cout << "Error opening the pos file for plate analysis." << endl;
         exit(0);
      }

      // Find buckling load factor

      string label;
      double buckfactor = 0;
    /*  cVector genstress(6);
      cVector genstrain(6);
      cVector force(8);
      force.Zero();*/
      int mode;

      while (pos >> label)
      {

          if (label == "%RESULT.CASE.STEP.BUCKLING.FACTOR")
          {
              pos >> mode;
              pos >> buckfactor;
          }
       }

      //cout << "                       \n buckfactor " << buckfactor << endl;
       if (buckfactor == 0)
       {
          cout << "Convergence not achieved in infill: " << endl;
          //exit(0);
       }

       // Push back the new targets Ybuck and Ystren
       lbdb = buckfactor;
      }
}

// ============================ Evaluate ==============================

void cSquarePlateMFBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    double vcratio;
    double Cmax = CostMax;

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
        c = vcratio - Cmax;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateMFBuckFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateTriDirMFBuckFGM :: cSquarePlateTriDirMFBuckFGM(void)
{
  NumConstr = 1;
  NumObj = 1;

  ResizeCP(6, 6, 4);
  CostMax = 0.7;
}

// ============================== Evaluate =================================

void cSquarePlateTriDirMFBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
    // Optimum (Cmax = 0.30)
    /*x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 0.0;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 0.0;
    x[6]  = 0.0;
    x[7]  = 0.5257;
    x[8]  = 1.0;
    x[9]  = 0.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 0.0;
    x[13] = 0.0;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 0.0;*/

    /*
    // Optimum BIOS (Cmax = 0.30)
    x[0]  = 1.0;
    x[1]  = 0.0;
    x[2]  = 0.0;
    x[3]  = 1.0;
    x[4]  = 0.3183;
    x[5]  = 1.0;
    x[6]  = 1.0;
    x[7]  = 0.0;
    x[8]  = 1.0;
    x[9]  = 0.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 0.0;
    x[13] = 0.0;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 0.0;
    */

    /*
    // Optimum (Cmax = 0.50)
    x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 0.5483;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 1.0;
    x[6]  = 0.9977;
    x[7]  = 1.0;
    x[8]  = 1.0;
    x[9]  = 1.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 0.0;
    x[13] = 0.0;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 0.0;
    */

    /*
    // Optimum BIOS (Cmax = 0.50)
    x[0]  = 1.0;
    x[1]  = 0.0;
    x[2]  = 1.0;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 1.0;
    x[6]  = 1.0;
    x[7]  = 1.0;
    x[8]  = 1.0;
    x[9]  = 1.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 0.0;
    x[13] = 0.0;
    x[14] = 0.0;
    x[15] = 0.615;
    x[16] = 0.0;
    x[17] = 0.0;
    */

    /*
    // Optimum (Cmax = 0.70)
    x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 1.0;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 1.0;
    x[6]  = 1.0;
    x[7]  = 1.0;
    x[8]  = 1.0;
    x[9]  = 1.0;
    x[10] = 1.0;
    x[11] = 0.0;
    x[12] = 1.0;
    x[13] = 0.2419;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 1.0;
    */

    /*
    // Optimum BIOS (Cmax = 0.70)
    x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 1.0;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 1.0;
    x[6]  = 1.0;
    x[7]  = 1.0;
    x[8]  = 1.0;
    x[9]  = 1.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 1.0;
    x[13] = 0.0;
    x[14] = 1.0;
    x[15] = 1.0;
    x[16] = 0.0;
    x[17] = 0.408;
    */

    // Objective function evaluation (using FAST)

    double fobj;
    Analysis(x, fobj, 2);  // Linearized Buckling Analysis
    fobjs[0] = -fobj;   // Maximization problem!
    cout << "fobj = " << fobj << endl;
    exit(0);

    // Constraint evaluation

    // Ceramic volume percentage < Cmax
    // Cmax = 65%

    double vcratio;
    double Cmax = CostMax;

    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    Vcp[0]  = Vcp[5]  = Vcp[30] = Vcp[35] = Vcp[108] = Vcp[113] = Vcp[138] = Vcp[143] = x[0];
    Vcp[1]  = Vcp[4]  = Vcp[31] = Vcp[34] = Vcp[109] = Vcp[112] = Vcp[139] = Vcp[142] = x[1];
    Vcp[2]  = Vcp[3]  = Vcp[32] = Vcp[33] = Vcp[110] = Vcp[111] = Vcp[140] = Vcp[141] = x[2];

    Vcp[6]  = Vcp[11] = Vcp[24] = Vcp[29] = Vcp[114] = Vcp[119] = Vcp[132] = Vcp[137] = x[3];
    Vcp[7]  = Vcp[10] = Vcp[25] = Vcp[28] = Vcp[115] = Vcp[118] = Vcp[133] = Vcp[136] = x[4];
    Vcp[8]  = Vcp[9]  = Vcp[26] = Vcp[27] = Vcp[116] = Vcp[117] = Vcp[134] = Vcp[135] = x[5];

    Vcp[12] = Vcp[17] = Vcp[18] = Vcp[23] = Vcp[120] = Vcp[125] = Vcp[126] = Vcp[131] = x[6];
    Vcp[13] = Vcp[16] = Vcp[19] = Vcp[22] = Vcp[121] = Vcp[124] = Vcp[127] = Vcp[130] = x[7];
    Vcp[14] = Vcp[15] = Vcp[20] = Vcp[21] = Vcp[122] = Vcp[123] = Vcp[128] = Vcp[129] = x[8];

    Vcp[36] = Vcp[41] = Vcp[66] = Vcp[71] = Vcp[72]  = Vcp[77]  = Vcp[102] = Vcp[107] = x[9];
    Vcp[37] = Vcp[40] = Vcp[67] = Vcp[70] = Vcp[73]  = Vcp[76]  = Vcp[103] = Vcp[106] = x[10];
    Vcp[38] = Vcp[39] = Vcp[68] = Vcp[69] = Vcp[74]  = Vcp[75]  = Vcp[104] = Vcp[105] = x[11];

    Vcp[42] = Vcp[47] = Vcp[60] = Vcp[65] = Vcp[78]  = Vcp[83]  = Vcp[96]  = Vcp[101] = x[12];
    Vcp[43] = Vcp[46] = Vcp[61] = Vcp[64] = Vcp[79]  = Vcp[82]  = Vcp[97]  = Vcp[100] = x[13];
    Vcp[44] = Vcp[45] = Vcp[62] = Vcp[63] = Vcp[80]  = Vcp[81]  = Vcp[98]  = Vcp[99]  = x[14];

    Vcp[48] = Vcp[53] = Vcp[54] = Vcp[59] = Vcp[84]  = Vcp[89]  = Vcp[90]  = Vcp[95]  = x[15];
    Vcp[49] = Vcp[52] = Vcp[55] = Vcp[58] = Vcp[85]  = Vcp[88]  = Vcp[91]  = Vcp[94]  = x[16];
    Vcp[50] = Vcp[51] = Vcp[56] = Vcp[57] = Vcp[86]  = Vcp[87]  = Vcp[92]  = Vcp[93]  = x[17];

    //cout << "1" << endl;
    EvalVolumeRatio3D(Vcp, vcratio, 6, 6, 4);
    //cout << "2" << endl;

    c[0] = vcratio - Cmax;
}

// ============================== Evaluate =================================

void cSquarePlateTriDirMFBuckFGM :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{

    /*// Optimum (Cmax = 0.30)
    x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 0.0;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 0.0;
    x[6]  = 0.0;
    x[7]  = 0.5257;
    x[8]  = 1.0;
    x[9]  = 0.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 0.0;
    x[13] = 0.0;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 0.0;



    // Optimum (Cmax = 0.50)
    x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 0.5483;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 1.0;
    x[6]  = 0.9977;
    x[7]  = 1.0;
    x[8]  = 1.0;
    x[9]  = 1.0;
    x[10] = 0.0;
    x[11] = 0.0;
    x[12] = 0.0;
    x[13] = 0.0;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 0.0;



    // Optimum (Cmax = 0.70)
    x[0]  = 1.0;
    x[1]  = 1.0;
    x[2]  = 1.0;
    x[3]  = 1.0;
    x[4]  = 1.0;
    x[5]  = 1.0;
    x[6]  = 1.0;
    x[7]  = 1.0;
    x[8]  = 1.0;
    x[9]  = 1.0;
    x[10] = 1.0;
    x[11] = 0.0;
    x[12] = 1.0;
    x[13] = 0.2419;
    x[14] = 0.0;
    x[15] = 0.0;
    x[16] = 0.0;
    x[17] = 1.0;*/


    // Objective function evaluation (using FAST)

    double fobj;
    Analysis(x, fobj, 1);  // Linearized Buckling Analysis
    fobjs[0] = -fobj;   // Maximization problem!

    // Constraint evaluation

    // Ceramic volume percentage < Cmax
    // Cmax = 50%

    double vcratio;
    double Cmax = CostMax;

    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    Vcp[0]  = Vcp[5]  = Vcp[30] = Vcp[35] = Vcp[108] = Vcp[113] = Vcp[138] = Vcp[143] = x[0];
    Vcp[1]  = Vcp[4]  = Vcp[31] = Vcp[34] = Vcp[109] = Vcp[112] = Vcp[139] = Vcp[142] = x[1];
    Vcp[2]  = Vcp[3]  = Vcp[32] = Vcp[33] = Vcp[110] = Vcp[111] = Vcp[140] = Vcp[141] = x[2];

    Vcp[6]  = Vcp[11] = Vcp[24] = Vcp[29] = Vcp[114] = Vcp[119] = Vcp[132] = Vcp[137] = x[3];
    Vcp[7]  = Vcp[10] = Vcp[25] = Vcp[28] = Vcp[115] = Vcp[118] = Vcp[133] = Vcp[136] = x[4];
    Vcp[8]  = Vcp[9]  = Vcp[26] = Vcp[27] = Vcp[116] = Vcp[117] = Vcp[134] = Vcp[135] = x[5];

    Vcp[12] = Vcp[17] = Vcp[18] = Vcp[23] = Vcp[120] = Vcp[125] = Vcp[126] = Vcp[131] = x[6];
    Vcp[13] = Vcp[16] = Vcp[19] = Vcp[22] = Vcp[121] = Vcp[124] = Vcp[127] = Vcp[130] = x[7];
    Vcp[14] = Vcp[15] = Vcp[20] = Vcp[21] = Vcp[122] = Vcp[123] = Vcp[128] = Vcp[129] = x[8];

    Vcp[36] = Vcp[41] = Vcp[66] = Vcp[71] = Vcp[72]  = Vcp[77]  = Vcp[102] = Vcp[107] = x[9];
    Vcp[37] = Vcp[40] = Vcp[67] = Vcp[70] = Vcp[73]  = Vcp[76]  = Vcp[103] = Vcp[106] = x[10];
    Vcp[38] = Vcp[39] = Vcp[68] = Vcp[69] = Vcp[74]  = Vcp[75]  = Vcp[104] = Vcp[105] = x[11];

    Vcp[42] = Vcp[47] = Vcp[60] = Vcp[65] = Vcp[78]  = Vcp[83]  = Vcp[96]  = Vcp[101] = x[12];
    Vcp[43] = Vcp[46] = Vcp[61] = Vcp[64] = Vcp[79]  = Vcp[82]  = Vcp[97]  = Vcp[100] = x[13];
    Vcp[44] = Vcp[45] = Vcp[62] = Vcp[63] = Vcp[80]  = Vcp[81]  = Vcp[98]  = Vcp[99]  = x[14];

    Vcp[48] = Vcp[53] = Vcp[54] = Vcp[59] = Vcp[84]  = Vcp[89]  = Vcp[90]  = Vcp[95]  = x[15];
    Vcp[49] = Vcp[52] = Vcp[55] = Vcp[58] = Vcp[85]  = Vcp[88]  = Vcp[91]  = Vcp[94]  = x[16];
    Vcp[50] = Vcp[51] = Vcp[56] = Vcp[57] = Vcp[86]  = Vcp[87]  = Vcp[92]  = Vcp[93]  = x[17];

    //cout << "1" << endl;
    EvalVolumeRatio3D(Vcp, vcratio, 6, 6, 4);
    //cout << "2" << endl;

    c[0] = vcratio - Cmax;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateTriDirMFBuckFGM :: Analysis(cVector x, double &lbdb, int f)
{
    // Evaluate volume fraction at gauss points according to a given distribution

    cVector Vcpg;

    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    Vcp[0]  = Vcp[5]  = Vcp[30] = Vcp[35] = Vcp[108] = Vcp[113] = Vcp[138] = Vcp[143] = x[0];
    Vcp[1]  = Vcp[4]  = Vcp[31] = Vcp[34] = Vcp[109] = Vcp[112] = Vcp[139] = Vcp[142] = x[1];
    Vcp[2]  = Vcp[3]  = Vcp[32] = Vcp[33] = Vcp[110] = Vcp[111] = Vcp[140] = Vcp[141] = x[2];

    Vcp[6]  = Vcp[11] = Vcp[24] = Vcp[29] = Vcp[114] = Vcp[119] = Vcp[132] = Vcp[137] = x[3];
    Vcp[7]  = Vcp[10] = Vcp[25] = Vcp[28] = Vcp[115] = Vcp[118] = Vcp[133] = Vcp[136] = x[4];
    Vcp[8]  = Vcp[9]  = Vcp[26] = Vcp[27] = Vcp[116] = Vcp[117] = Vcp[134] = Vcp[135] = x[5];

    Vcp[12] = Vcp[17] = Vcp[18] = Vcp[23] = Vcp[120] = Vcp[125] = Vcp[126] = Vcp[131] = x[6];
    Vcp[13] = Vcp[16] = Vcp[19] = Vcp[22] = Vcp[121] = Vcp[124] = Vcp[127] = Vcp[130] = x[7];
    Vcp[14] = Vcp[15] = Vcp[20] = Vcp[21] = Vcp[122] = Vcp[123] = Vcp[128] = Vcp[129] = x[8];

    Vcp[36] = Vcp[41] = Vcp[66] = Vcp[71] = Vcp[72]  = Vcp[77]  = Vcp[102] = Vcp[107] = x[9];
    Vcp[37] = Vcp[40] = Vcp[67] = Vcp[70] = Vcp[73]  = Vcp[76]  = Vcp[103] = Vcp[106] = x[10];
    Vcp[38] = Vcp[39] = Vcp[68] = Vcp[69] = Vcp[74]  = Vcp[75]  = Vcp[104] = Vcp[105] = x[11];

    Vcp[42] = Vcp[47] = Vcp[60] = Vcp[65] = Vcp[78]  = Vcp[83]  = Vcp[96]  = Vcp[101] = x[12];
    Vcp[43] = Vcp[46] = Vcp[61] = Vcp[64] = Vcp[79]  = Vcp[82]  = Vcp[97]  = Vcp[100] = x[13];
    Vcp[44] = Vcp[45] = Vcp[62] = Vcp[63] = Vcp[80]  = Vcp[81]  = Vcp[98]  = Vcp[99]  = x[14];

    Vcp[48] = Vcp[53] = Vcp[54] = Vcp[59] = Vcp[84]  = Vcp[89]  = Vcp[90]  = Vcp[95]  = x[15];
    Vcp[49] = Vcp[52] = Vcp[55] = Vcp[58] = Vcp[85]  = Vcp[88]  = Vcp[91]  = Vcp[94]  = x[16];
    Vcp[50] = Vcp[51] = Vcp[56] = Vcp[57] = Vcp[86]  = Vcp[87]  = Vcp[92]  = Vcp[93]  = x[17];

    double thk = 1.0;

  int num_thread = 0;
#ifdef _OMP_
  num_thread = omp_get_thread_num( );
#endif

  stringstream thread;
  thread << num_thread;

  string fid;
  if (f == 1)
      fid = "2D";
  else
      fid = "3D";


  string thread_number = thread.str();
  string cmd  = "del SqrPltBuck" + fid + thread_number + ".dat";
  string cmd2 = "del SqrPltBuck" + fid + thread_number + ".pos";
  string cmd3 = "rm SqrPltBuck" + fid + thread_number + ".dat";
  string cmd4 = "rm SqrPltBuck" + fid + thread_number + ".pos";

#ifdef _WIN32
  if (system(cmd.c_str()) || system(cmd2.c_str()))
     cout << "Problem on removing SqrPltBuck.dat and plate.pos files.\n";
#else
  if (system(cmd3.c_str()) || system(cmd4.c_str()))
     cout << "Problem on removing SqrPltBuck.dat and plate.pos files.\n";
#endif

  string begname, endname, datname, posname;
  if (f == 1)
  {
      begname = "datbegSqrPltBuck3DirSS22D.dat";
      endname = "datendSqrPltBuck3DirSS22D.dat";

      datname = "SqrPltBuck2D" + thread_number + ".dat";
      posname = "SqrPltBuck2D" + thread_number + ".pos";
  }
  else
  {
      begname = "datbegSqrPltBuck3Dir3D.dat";
      endname = "datendSqrPltBuck3Dir3D.dat";

      datname = "SqrPltBuck3D" + thread_number + ".dat";
      posname = "SqrPltBuck3D" + thread_number + ".pos";
  }

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

  if (f == 1)
  {
      dat << "%SECTION.FGM.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1    " << thk << "    10    3    " << numcp+10 << endl;
      dat << "-5.0  5.0  -5.0  5.0" << endl; // lxlow; lxupp; lylow; lyupp;
      dat << "6  6  4" << endl; // ncp_x; ncp_y; ncp_z;
      dat << "3  3  3" << endl; // Cubic in all coordinates
      for (int i = 0; i < numcp/6; i++) dat << Vcp[i*6 + 0] << "  " << Vcp[i*6 + 1] << "  " << Vcp[i*6 + 2] << "  " << Vcp[i*6 + 3] << "  " << Vcp[i*6 + 4] << "  " << Vcp[i*6 + 5] << endl;
      dat << endl;
  }
  else
  {
      dat << "%SECTION.FGM.3D" << endl;
      dat << "1" << endl;
      dat << "1    1    4    " << numcp+12 << endl;
      dat << "0.0  10.0  0.0  10.0  0.0  1.0" << endl; // lxlow; lxupp; lylow; lyupp;
      dat << "6  6  4" << endl; // ncp_x; ncp_y; ncp_z;
      dat << "3  3  3" << endl; // Cubic in all coordinates
      for (int i = 0; i < numcp/6; i++) dat << Vcp[i*6 + 0] << "  " << Vcp[i*6 + 1] << "  " << Vcp[i*6 + 2] << "  " << Vcp[i*6 + 3] << "  " << Vcp[i*6 + 4] << "  " << Vcp[i*6 + 5] << endl;
      dat << endl;

  }

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
  cmd = "fast.exe SqrPltBuck" + fid + thread_number + " -silent";
#else
  cmd = "./fast SqrPltBuck" + fid + thread_number + " -silent";
#endif
  status2 = system(cmd.c_str( ));

  if (status2)
  {
     cout << "Error in the analysis with fast.";
    #ifdef _WIN32
      cmd = "fast.exe SqrPltBuck" + fid + thread_number + " -silent";
    #else
      cmd = "./fast SqrPltBuck" + fid + thread_number + " -silent";
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

void cSquarePlateTriDirMFBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    // Ceramic volume percentage < Cmax
    // Cmax = 30%

    double vcratio;
    double Cmax = CostMax;

    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    Vcp[0]  = Vcp[5]  = Vcp[30] = Vcp[35] = Vcp[108] = Vcp[113] = Vcp[138] = Vcp[143] = x[0];
    Vcp[1]  = Vcp[4]  = Vcp[31] = Vcp[34] = Vcp[109] = Vcp[112] = Vcp[139] = Vcp[142] = x[1];
    Vcp[2]  = Vcp[3]  = Vcp[32] = Vcp[33] = Vcp[110] = Vcp[111] = Vcp[140] = Vcp[141] = x[2];

    Vcp[6]  = Vcp[11] = Vcp[24] = Vcp[29] = Vcp[114] = Vcp[119] = Vcp[132] = Vcp[137] = x[3];
    Vcp[7]  = Vcp[10] = Vcp[25] = Vcp[28] = Vcp[115] = Vcp[118] = Vcp[133] = Vcp[136] = x[4];
    Vcp[8]  = Vcp[9]  = Vcp[26] = Vcp[27] = Vcp[116] = Vcp[117] = Vcp[134] = Vcp[135] = x[5];

    Vcp[12] = Vcp[17] = Vcp[18] = Vcp[23] = Vcp[120] = Vcp[125] = Vcp[126] = Vcp[131] = x[6];
    Vcp[13] = Vcp[16] = Vcp[19] = Vcp[22] = Vcp[121] = Vcp[124] = Vcp[127] = Vcp[130] = x[7];
    Vcp[14] = Vcp[15] = Vcp[20] = Vcp[21] = Vcp[122] = Vcp[123] = Vcp[128] = Vcp[129] = x[8];

    Vcp[36] = Vcp[41] = Vcp[66] = Vcp[71] = Vcp[72]  = Vcp[77]  = Vcp[102] = Vcp[107] = x[9];
    Vcp[37] = Vcp[40] = Vcp[67] = Vcp[70] = Vcp[73]  = Vcp[76]  = Vcp[103] = Vcp[106] = x[10];
    Vcp[38] = Vcp[39] = Vcp[68] = Vcp[69] = Vcp[74]  = Vcp[75]  = Vcp[104] = Vcp[105] = x[11];

    Vcp[42] = Vcp[47] = Vcp[60] = Vcp[65] = Vcp[78]  = Vcp[83]  = Vcp[96]  = Vcp[101] = x[12];
    Vcp[43] = Vcp[46] = Vcp[61] = Vcp[64] = Vcp[79]  = Vcp[82]  = Vcp[97]  = Vcp[100] = x[13];
    Vcp[44] = Vcp[45] = Vcp[62] = Vcp[63] = Vcp[80]  = Vcp[81]  = Vcp[98]  = Vcp[99]  = x[14];

    Vcp[48] = Vcp[53] = Vcp[54] = Vcp[59] = Vcp[84]  = Vcp[89]  = Vcp[90]  = Vcp[95]  = x[15];
    Vcp[49] = Vcp[52] = Vcp[55] = Vcp[58] = Vcp[85]  = Vcp[88]  = Vcp[91]  = Vcp[94]  = x[16];
    Vcp[50] = Vcp[51] = Vcp[56] = Vcp[57] = Vcp[86]  = Vcp[87]  = Vcp[92]  = Vcp[93]  = x[17];

    //cout << "1" << endl;
    EvalVolumeRatio3D(Vcp, vcratio, 6, 6, 4);
    //cout << "2" << endl;

  // Single constraint evaluation.
    if (index == 0){
        c = vcratio - Cmax;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateTriDirMFBuckFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateCutOutFGM :: cSquarePlateCutOutFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
  FreqNatNormMin = 0.009;
  FreqNatNormMax = 0.015;
}

// ============================== Evaluate =================================

void cSquarePlateCutOutFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // for (int i = 0; i < 7; i++) x[i] = 0.0;
  // x[0] = 1.0; x[1] = 1.0; x[2] = 1.0; x[3] = 1.0; x[4] = 1.0; x[5] = 0.475; x[6] = 0.0;

  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj, 2);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation (using FAST)

  double Fmax = FreqNatNormMax; double Fmin = FreqNatNormMin;

  double ce;
  AnalysisC(x, ce, 2);    // Natural frequency

  // Normalization
  double h = 0.20; double rhoc = 2370; double Ec = 348.43e9; double nuc = 0.30; double Gc = Ec/(2*(1 + nuc));
  double NormFac = h*pow(rhoc/Gc, 0.5);
  double F = ce*NormFac;

  // cout << "F = " << F << "  NormFac = " << NormFac << "    ce = " << ce <<  endl;
  // exit(0);

  double cmax = F/Fmax - 1;
  double cmin = 1 - F/Fmin;

  if (cmax > cmin)
      c[0] = cmax;
  else
      c[0] = cmin;

}

// ============================== Evaluate =================================

void cSquarePlateCutOutFGM :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Objective function evaluation (using FAST)

    double fobj;
    Analysis(x, fobj, 1);  // Linearized Buckling Analysis
    fobjs[0] = -fobj;   // Maximization problem!

    // Constraint evaluation (using FAST)

    double Fmax = FreqNatNormMax; double Fmin = FreqNatNormMin;

    double ce;
    AnalysisC(x, ce, 1);    // Natural frequency

    // Normalization
    double h = 0.20; double rhoc = 2370; double Ec = 348.43e9; double nuc = 0.30; double Gc = Ec/(2*(1 + nuc));
    double NormFac = h*pow(rhoc/Gc, 0.5);
    double F = ce*NormFac;

    double cmax = F/Fmax - 1;
    double cmin = 1 - F/Fmin;

    if (cmax > cmin)
        c[0] = cmax;
    else
        c[0] = cmin;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateCutOutFGM :: Analysis(cVector x, double &lbdb, int f)
{
    // Define number of gauss points

        int ngauss = 10;

        // Evaluate weights and abscissas of gauss points

        cVector r, w;
        GaussPts1D(ngauss, r, w);

        // Transform abscissas to thickness coordinate

        cVector t(ngauss);

        // for (int i = 0; i < ngauss; i++) t[i]  = (r[i] + 1)/2;  // t = 0 (bottom) and t = 1 (top)
        for (int i = 0; i < ngauss; i++) t[i]  = r[i]/2;  // t = -0.5 (bottom) and t = 0.5 (top)

        // Evaluate volume fraction at gauss points according to a given distribution

        cVector Vcpg;

        /*t.Resize(21);
        t[0] = -0.5;
        for (int i = 1; i < 21; i++) t[i]  = t[i-1] + (1.0/20.0);

        cVector V(5);
        V[0] = 1.0;
        V[1] = 0.84572;
        V[2] = 1.0;
        V[3] = 1.0;
        V[4] = 1.0;

        PiecewiseCubicInterpolation(V, 21, t, Vcpg);

        cout << "Vcpg = " << endl;
        Vcpg.Print();

        exit(0);*/

        /*cVector Vcpg;

        VolumeDist(FGMVolDist, ngauss, t, Vcpg, x[1]);*/

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

        double thk = 0.2;

        /*cout << "\n\n ============ Matriz A ======= " << endl;
        A.Print();
        cout << "\n ============ Matriz B ======= " << endl;
        B.Print();
        cout << "\n ============ Matriz D ======= " << endl;
        D.Print();
        cout << "\n ============ Matriz G ======= " << endl;
        G.Print();
        cout << "\n ============ Matriz ABDG ======= " << endl;
        ABDG.Print();*/

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd, cmd2, cmd3, cmd4;
      if (f == 1)
      {
          cmd  = "del PlateCutOutLF" + thread_number + ".dat";
          cmd2 = "del PlateCutOutLF" + thread_number + ".pos";
          cmd3 = "rm PlateCutOutLF" + thread_number + ".dat";
          cmd4 = "rm PlateCutOutLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          cmd  = "del PlateCutOutHF" + thread_number + ".dat";
          cmd2 = "del PlateCutOutHF" + thread_number + ".pos";
          cmd3 = "rm PlateCutOutHF" + thread_number + ".dat";
          cmd4 = "rm PlateCutOutHF" + thread_number + ".pos";
      }


    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing dat and pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing dat and pos files.\n";
    #endif

      string begname, endname, datname, posname;

      int NumElm;

      if (f == 1)
      {
          begname = "datbegPlateCutoutThermBuck_n4p3.dat";
          endname = "datendPlateCutoutThermBuck_n4p3.dat";

          NumElm = 32;

          datname = "PlateCutOutLF" + thread_number + ".dat";
          posname = "PlateCutOutLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          begname = "datbegPlateCutoutThermBuck_n32p3.dat";
          endname = "datendPlateCutoutThermBuck_n32p3.dat";

          NumElm = 2048;

          datname = "PlateCutOutHF" + thread_number + ".dat";
          posname = "PlateCutOutHF" + thread_number + ".pos";
      }

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

         //exit(EXIT_FAILURE);

         lbdb = 0.0;

         cout << "chega aqui" << endl;
         //exit(0);
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

      dat << endl << endl << "%SECTION.FGM.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1    " << thk << "    10    2    " << numcp;

      for (int i = 0; i < numcp; i++)
      {
          dat << "    " << Vcp[i];
      }
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
      //   exit(EXIT_FAILURE);

         lbdb = 0.0;
         exit(0);
      }

      // Write the element temperature

      double Telm = 1;

      fstream dat2;

      dat2.open(datname.c_str( ));

      if (!dat2.is_open( ))
      {
         cout << "Error opening the dat file for plate analysis." << endl;
         exit(0);
      }

      dat2.seekp(0,ofstream::end);

      dat2 << endl << endl << "%LOAD.CASE.ELEMENT.TEMPERATURE" << endl;
      dat2 << NumElm << endl;
      for (int i = 0; i < NumElm; i++)
      {
          // dat2 << i+1 << "    " << p <<  "    " << nt << "    ";
          // for (int j = 0; j < nt; j++) dat2 << T[j] << "    ";
          // dat2 << endl;
          dat2 << i+1 << "    " << 1 <<  "    " << 2 << "    ";
          for (int j = 0; j < 2; j++) dat2 << Telm << "    ";
          dat2 << endl;
      }
      dat2 << endl;
      dat2 << "%END";

      dat2.close( );

      // Run the analysis with FAST.

    string fid;
    if (f == 1)
    {
        fid = "LF";
    }
    else
    {
        fid = "HF";
    }

    #ifdef _WIN32
      cmd = "fast.exe PlateCutOut" + fid + thread_number + " -silent";
    #else
      cmd = "./fast PlateCutOut" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
         //exit(EXIT_FAILURE);

    #ifdef _WIN32
      cmd = "fast.exe PlateCutOut" + fid + thread_number + " -silent";
    #else
      cmd = "./fast PlateCutOut" + fid + thread_number + " -silent";
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
      //string posname = "platehole.pos";
      //pos.open(posname.c_str());
      if (!pos.is_open( ))
      {
         cout << "Error opening the pos file for plate analysis." << endl;
         exit(0);
      }

      // Find buckling load factor

      string label;
      double buckfactor = 0;
    /*  cVector genstress(6);
      cVector genstrain(6);
      cVector force(8);
      force.Zero();*/
      int mode;

      // Maximum failure index.

        /*double S11, S22, S33, T12, T13, T23, term1, term2, term3, term4, Svm;

        double fi;
        double maxFI = 0.0;

        cVector GenStress(8);
        cVector GenStrain(8);

        cVector BendingStrain(3);
        cVector ShearStrain(2);

        cVector BendingStress(3);
        cVector ShearStress(2);

        cVector Stress(3);

        cMatrix C(8,8);
        cMatrix lQb(3,3);
        cMatrix lQs(2,2);
        cMatrix Tb(3,3);
        cMatrix Ts(2,2);

        C.Zero( );

        C = ABDG;

        cMatrix S(8,8);
        C.CompInverse(S);

        int nsteps = 10;
        cVector CoordZ(nsteps + 1);
        cVector CoordZn(nsteps + 1);
        CoordZ.Zero( ); CoordZn.Zero( );

        CoordZ[0]  = -thk/2.0;
        CoordZn[0] = -1.0;

        for (int lam = 0; lam < nsteps; lam++)
          CoordZ[lam+1] = CoordZ[lam] + thk/((double)nsteps);

        for (int lam = 0; lam < nsteps; lam++)
          CoordZn[lam+1] = CoordZn[lam] + 2.0/((double)nsteps);

        cVector VcSteps(nsteps+1);
        VolumeDist(FGMVolDist, VcSteps, nsteps+1, CoordZ, VcSteps);

        // Evaluate the effective properties at gauss points

        cVector En, Nun, Kn, Gn, Rhon;

        EffPropModel(FGMModel, VcSteps, En, Nun, Kn, Gn, Rhon);*/

        while (pos >> label)
        {
          int numelm, elmid, npg;

          if (label == "%RESULT.CASE.STEP.BUCKLING.FACTOR")
          {
              pos >> mode;
              pos >> buckfactor;
          }

          /*if (label == "%RESULT.CASE.STEP.ELEMENT.GAUSS.SCALAR.DATA")
          {
            pos >> numelm;

            for (int i = 0; i < numelm; i++)
            {
              pos >> elmid;
              pos >> npg;

              for (int j = 0; j < npg; j++)
              {
                GenStress.Zero( );
                pos >> GenStress[0] >> GenStress[1] >> GenStress[2] >> GenStress[3] >> GenStress[4] >> GenStress[5] >> GenStress[6] >> GenStress[7];

                GenStress[3] *= -1.0;
                GenStress[4] *= -1.0;
                GenStress[5] *= -1.0;
                GenStress[6] *= -1.0;
                GenStress[7] *= -1.0;

                GenStrain.Zero( );
                GenStrain = S*GenStress;

                for (int lam = 0; lam < nsteps; lam++)
                {
                  // Get the constitutive and transformation matrices.

                  BendingStrain.Zero( );
                  ShearStrain.Zero( );
                  Stress.Zero( );
                  lQb.Zero( );
                  lQs.Zero( );
                  Tb.Zero( );
                  Ts.Zero( );

                  QMatrix(En[lam], Nun[lam], lQb, lQs);

                  // Lower Point.

                  BendingStrain[0] = GenStrain[0] + CoordZ[lam]*GenStrain[3];
                  BendingStrain[1] = GenStrain[1] + CoordZ[lam]*GenStrain[4];
                  BendingStrain[2] = GenStrain[2] + CoordZ[lam]*GenStrain[5];

                  ShearStrain[0] = GenStrain[6];
                  ShearStrain[1] = GenStrain[7];

                  BendingStress = lQb*BendingStrain;  // S_XX S_YY T_XY
                  ShearStress = lQs*ShearStrain;      // T_XZ T_YZ

                  S11 = BendingStress[0];
                  S22 = BendingStress[1];
                  S33 = 0.0;
                  T12 = BendingStress[2];
                  T13 = ShearStress[0];
                  T23 = ShearStress[1];

                  term1 = pow((S11 - S22), 2);
                  term2 = pow((S22 - S33), 2);
                  term3 = pow((S22 - S11), 2);
                  term4 = 6*(T23*T23 + T13*T13 + T12*T12);

                  Svm = sqrt(0.5*(term1 + term2 + term3 + term4));

                  double Em = FGMMat[0];
                  double Ec = FGMMat[3];

                  double q   = 90.0e9;   // Stress transfer parameter
                  double Sym = 493.7e6; // Yield stress
                  double Sy  = Sym*((1 - VcSteps[lam]) + (q + Em)/(q + Ec)*(Ec/Em)*VcSteps[lam]);

                  fi = Svm/Sy;

                  if (Svm > maxFI) maxFI = fi;
                }
              }
            }
          }*/
        }

      //cout << "                       \n buckfactor " << buckfactor << endl;
       if (buckfactor == 0)
       {
          cout << "Convergence not achieved in infill: " << endl;
          //exit(0);
       }

       // Push back the new targets Ybuck and Ystren
       lbdb = buckfactor;
      }
}

// ============================== Analysis =================================

void cSquarePlateCutOutFGM :: AnalysisC(cVector x, double &vib, int f)
{
    // Define number of gauss points

        int ngauss = 10;

        // Evaluate weights and abscissas of gauss points

        cVector r, w;
        GaussPts1D(ngauss, r, w);

        // Transform abscissas to thickness coordinate

        cVector t(ngauss);

        // for (int i = 0; i < ngauss; i++) t[i]  = (r[i] + 1)/2;  // t = 0 (bottom) and t = 1 (top)
        for (int i = 0; i < ngauss; i++) t[i]  = r[i]/2;  // t = -0.5 (bottom) and t = 0.5 (top)

        // Evaluate volume fraction at gauss points according to a given distribution

        cVector Vcpg;

        /*t.Resize(21);
        t[0] = -0.5;
        for (int i = 1; i < 21; i++) t[i]  = t[i-1] + (1.0/20.0);

        cVector V(5);
        V[0] = 1.0;
        V[1] = 0.84572;
        V[2] = 1.0;
        V[3] = 1.0;
        V[4] = 1.0;

        PiecewiseCubicInterpolation(V, 21, t, Vcpg);

        cout << "Vcpg = " << endl;
        Vcpg.Print();

        exit(0);*/

        /*cVector Vcpg;

        VolumeDist(FGMVolDist, ngauss, t, Vcpg, x[1]);*/

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

        double thk = 0.2;

        /*cout << "\n\n ============ Matriz A ======= " << endl;
        A.Print();
        cout << "\n ============ Matriz B ======= " << endl;
        B.Print();
        cout << "\n ============ Matriz D ======= " << endl;
        D.Print();
        cout << "\n ============ Matriz G ======= " << endl;
        G.Print();
        cout << "\n ============ Matriz ABDG ======= " << endl;
        ABDG.Print();*/

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd, cmd2, cmd3, cmd4;
      if (f == 1)
      {
          cmd  = "del PlateCutOutLF" + thread_number + ".dat";
          cmd2 = "del PlateCutOutLF" + thread_number + ".pos";
          cmd3 = "rm PlateCutOutLF" + thread_number + ".dat";
          cmd4 = "rm PlateCutOutLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          cmd  = "del PlateCutOutHF" + thread_number + ".dat";
          cmd2 = "del PlateCutOutHF" + thread_number + ".pos";
          cmd3 = "rm PlateCutOutHF" + thread_number + ".dat";
          cmd4 = "rm PlateCutOutHF" + thread_number + ".pos";
      }


    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing dat and pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing dat and pos files.\n";
    #endif

      string begname, endname, datname, posname;

      if (f == 1)
      {
          begname = "datbegPlateCutoutVib_n4p3.dat";
          endname = "datendPlateCutoutVib_n4p3.dat";

          datname = "PlateCutOutLF" + thread_number + ".dat";
          posname = "PlateCutOutLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          begname = "datbegPlateCutoutVib_n32p3.dat";
          endname = "datendPlateCutoutVib_n32p3.dat";

          datname = "PlateCutOutHF" + thread_number + ".dat";
          posname = "PlateCutOutHF" + thread_number + ".pos";
      }

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

         //exit(EXIT_FAILURE);

         vib = 0.0;

         cout << "chega aqui" << endl;
         //exit(0);
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

      dat << endl << endl << "%SECTION.FGM.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1    " << thk << "    10    2    " << numcp;

      for (int i = 0; i < numcp; i++)
      {
          dat << "    " << Vcp[i];
      }
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
      //   exit(EXIT_FAILURE);

         vib = 0.0;
         exit(0);
      }

      // Run the analysis with FAST.

    string fid;
    if (f == 1)
    {
        fid = "LF";
    }
    else
    {
        fid = "HF";
    }

    #ifdef _WIN32
      cmd = "fast.exe PlateCutOut" + fid + thread_number + " -silent";
    #else
      cmd = "./fast PlateCutOut" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
         //exit(EXIT_FAILURE);

    #ifdef _WIN32
      cmd = "fast.exe PlateCutOut" + fid + thread_number + " -silent";
    #else
      cmd = "./fast PlateCutOut" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
          vib = 0.0;
      }

      }

      if (!status2)
      {
      // Open the pos file.

      ifstream pos;

      pos.open(posname.c_str( ));
      //string posname = "platehole.pos";
      //pos.open(posname.c_str());
      if (!pos.is_open( ))
      {
         cout << "Error opening the pos file for plate analysis." << endl;
         exit(0);
      }

      // Find buckling load factor

      string label;
      double vibfactor = 0;
    /*  cVector genstress(6);
      cVector genstrain(6);
      cVector force(8);
      force.Zero();*/
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
          //exit(0);
       }

       // Push back the new targets Ybuck and Ystren
       vib = vibfactor;
      }
}

// ============================== Analysis =================================

void cSquarePlateCutOutFGM :: AnalysisStress(cVector x, double &lbdb, int f)
{
    // Define number of gauss points

        int ngauss = 10;

        // Evaluate weights and abscissas of gauss points

        cVector r, w;
        GaussPts1D(ngauss, r, w);

        // Transform abscissas to thickness coordinate

        cVector t(ngauss);

        // for (int i = 0; i < ngauss; i++) t[i]  = (r[i] + 1)/2;  // t = 0 (bottom) and t = 1 (top)
        for (int i = 0; i < ngauss; i++) t[i]  = r[i]/2;  // t = -0.5 (bottom) and t = 0.5 (top)

        // Evaluate volume fraction at gauss points according to a given distribution

        cVector Vcpg;

        /*t.Resize(21);
        t[0] = -0.5;
        for (int i = 1; i < 21; i++) t[i]  = t[i-1] + (1.0/20.0);

        cVector V(5);
        V[0] = 1.0;
        V[1] = 0.84572;
        V[2] = 1.0;
        V[3] = 1.0;
        V[4] = 1.0;

        PiecewiseCubicInterpolation(V, 21, t, Vcpg);

        cout << "Vcpg = " << endl;
        Vcpg.Print();

        exit(0);*/

        /*cVector Vcpg;

        VolumeDist(FGMVolDist, ngauss, t, Vcpg, x[1]);*/

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
        cout << "Vcp: ";
        Vcp.Print( );

        double thk = 0.1;

        VolumeDist(FGMVolDist, Vcp, ngauss, t, Vcpg);

        // Evaluate the effective properties at gauss points

        cVector Epg, Nupg, Kpg, Gpg, Rhopg;

        EffPropModel(FGMModel, Vcpg, Epg, Nupg, Kpg, Gpg, Rhopg);

        // Evaluate the ABDG matrices

        cMatrix A, B, D, G, ABDG;

        CalcABDG(thk, r, w, Epg, Nupg, A, B, D, G, ABDG);

        /*cout << "\n\n ============ Matriz A ======= " << endl;
        A.Print();
        cout << "\n ============ Matriz B ======= " << endl;
        B.Print();
        cout << "\n ============ Matriz D ======= " << endl;
        D.Print();
        cout << "\n ============ Matriz G ======= " << endl;
        G.Print();
        cout << "\n ============ Matriz ABDG ======= " << endl;
        ABDG.Print();*/

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd, cmd2, cmd3, cmd4;
      if (f == 1)
      {
          cmd  = "del PlateCutOutLF" + thread_number + ".dat";
          cmd2 = "del PlateCutOutLF" + thread_number + ".pos";
          cmd3 = "rm PlateCutOutLF" + thread_number + ".dat";
          cmd4 = "rm PlateCutOutLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          cmd  = "del PlateCutOutHF" + thread_number + ".dat";
          cmd2 = "del PlateCutOutHF" + thread_number + ".pos";
          cmd3 = "rm PlateCutOutHF" + thread_number + ".dat";
          cmd4 = "rm PlateCutOutHF" + thread_number + ".pos";
      }


    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing dat and pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing dat and pos files.\n";
    #endif

      string begname, endname, datname, posname;

      int NumElm;

      if (f == 1)
      {
          begname = "datbegPlateCutoutBuck_n4p3.dat";
          endname = "datendPlateCutoutBuck_n4p3.dat";

          NumElm = 32;

          datname = "PlateCutOutLF" + thread_number + ".dat";
          posname = "PlateCutOutLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          begname = "datbegPlateCutoutBuck_n32p3.dat";
          endname = "datendPlateCutoutBuck_n32p3.dat";

          NumElm = 2048;

          datname = "PlateCutOutHF" + thread_number + ".dat";
          posname = "PlateCutOutHF" + thread_number + ".pos";
      }

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

         //exit(EXIT_FAILURE);

         lbdb = 0.0;

         cout << "chega aqui" << endl;
         //exit(0);
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

      dat << endl << endl << "%SECTION.FGM.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1    " << thk << "    10    2    " << numcp;

      for (int i = 0; i < numcp; i++)
      {
          dat << "    " << Vcp[i];
      }
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
      //   exit(EXIT_FAILURE);

         lbdb = 0.0;
         exit(0);
      }

      // Write the element temperature

      double Telm = 300;

      fstream dat2;

      dat2.open(datname.c_str( ));

      if (!dat2.is_open( ))
      {
         cout << "Error opening the dat file for plate analysis." << endl;
         exit(0);
      }

      dat2.seekp(0,ofstream::end);

      dat2 << endl << endl << "%LOAD.CASE.ELEMENT.TEMPERATURE" << endl;
      dat2 << NumElm << endl;
      for (int i = 0; i < NumElm; i++)
      {
          // dat2 << i+1 << "    " << p <<  "    " << nt << "    ";
          // for (int j = 0; j < nt; j++) dat2 << T[j] << "    ";
          // dat2 << endl;
          dat2 << i+1 << "    " << 1 <<  "    " << 2 << "    ";
          for (int j = 0; j < 2; j++) dat2 << Telm << "    ";
          dat2 << endl;
      }
      dat2 << endl;
      dat2 << "%END";

      dat2.close( );

      // Run the analysis with FAST.

    string fid;
    if (f == 1)
    {
        fid = "LF";
    }
    else
    {
        fid = "HF";
    }

    #ifdef _WIN32
      cmd = "fast.exe PlateCutOut" + fid + thread_number + " -silent";
    #else
      cmd = "./fast PlateCutOut" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
         //exit(EXIT_FAILURE);

    #ifdef _WIN32
      cmd = "fast.exe PlateCutOut" + fid + thread_number + " -silent";
    #else
      cmd = "./fast PlateCutOut" + fid + thread_number + " -silent";
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
      //string posname = "platehole.pos";
      //pos.open(posname.c_str());
      if (!pos.is_open( ))
      {
         cout << "Error opening the pos file for plate analysis." << endl;
         exit(0);
      }

      // Find buckling load factor

      string label;
      double buckfactor = 0;
    /*  cVector genstress(6);
      cVector genstrain(6);
      cVector force(8);
      force.Zero();*/
      int mode;

      // Maximum failure index.

        double S11, S22, S33, T12, T13, T23, term1, term2, term3, term4, Svm;

        double fi;
        double maxFI = 0.0;

        cVector GenStress(8);
        cVector GenStrain(8);

        cVector BendingStrain(3);
        cVector ShearStrain(2);

        cVector BendingStress(3);
        cVector ShearStress(2);

        cVector Stress(3);

        cMatrix C(8,8);
        cMatrix lQb(3,3);
        cMatrix lQs(2,2);
        cMatrix Tb(3,3);
        cMatrix Ts(2,2);

        C.Zero( );

        C = ABDG;

        cMatrix S(8,8);
        C.CompInverse(S);

        int nsteps = 10;
        cVector CoordZ(nsteps + 1);
        cVector CoordZn(nsteps + 1);
        CoordZ.Zero( ); CoordZn.Zero( );

        CoordZ[0]  = -thk/2.0;
        CoordZn[0] = -1.0;

        for (int lam = 0; lam < nsteps; lam++)
          CoordZ[lam+1] = CoordZ[lam] + thk/((double)nsteps);

        for (int lam = 0; lam < nsteps; lam++)
          CoordZn[lam+1] = CoordZn[lam] + 2.0/((double)nsteps);

        cVector VcSteps(nsteps+1);
        VolumeDist(FGMVolDist, VcSteps, nsteps+1, CoordZ, VcSteps);

        // Evaluate the effective properties at gauss points

        cVector En, Nun, Kn, Gn, Rhon;

        EffPropModel(FGMModel, VcSteps, En, Nun, Kn, Gn, Rhon);

        while (pos >> label)
        {
          int numelm, elmid, npg;

          if (label == "%RESULT.CASE.STEP.BUCKLING.FACTOR")
          {
              pos >> mode;
              pos >> buckfactor;
          }

          if (label == "%RESULT.CASE.STEP.ELEMENT.GAUSS.SCALAR.DATA")
          {
            pos >> numelm;

            for (int i = 0; i < numelm; i++)
            {
              pos >> elmid;
              pos >> npg;

              for (int j = 0; j < npg; j++)
              {
                GenStress.Zero( );
                pos >> GenStress[0] >> GenStress[1] >> GenStress[2] >> GenStress[3] >> GenStress[4] >> GenStress[5] >> GenStress[6] >> GenStress[7];

                GenStress[3] *= -1.0;
                GenStress[4] *= -1.0;
                GenStress[5] *= -1.0;
                GenStress[6] *= -1.0;
                GenStress[7] *= -1.0;

                GenStrain.Zero( );
                GenStrain = S*GenStress;

                for (int lam = 0; lam < nsteps; lam++)
                {
                  // Get the constitutive and transformation matrices.

                  BendingStrain.Zero( );
                  ShearStrain.Zero( );
                  Stress.Zero( );
                  lQb.Zero( );
                  lQs.Zero( );
                  Tb.Zero( );
                  Ts.Zero( );

                  QMatrix(En[lam], Nun[lam], lQb, lQs);

                  // Lower Point.

                  BendingStrain[0] = GenStrain[0] + CoordZ[lam]*GenStrain[3];
                  BendingStrain[1] = GenStrain[1] + CoordZ[lam]*GenStrain[4];
                  BendingStrain[2] = GenStrain[2] + CoordZ[lam]*GenStrain[5];

                  ShearStrain[0] = GenStrain[6];
                  ShearStrain[1] = GenStrain[7];

                  BendingStress = lQb*BendingStrain;  // S_XX S_YY T_XY
                  ShearStress = lQs*ShearStrain;      // T_XZ T_YZ

                  S11 = BendingStress[0];
                  S22 = BendingStress[1];
                  S33 = 0.0;
                  T12 = BendingStress[2];
                  T13 = ShearStress[0];
                  T23 = ShearStress[1];

                  term1 = pow((S11 - S22), 2);
                  term2 = pow((S22 - S33), 2);
                  term3 = pow((S22 - S11), 2);
                  term4 = 6*(T23*T23 + T13*T13 + T12*T12);

                  Svm = sqrt(0.5*(term1 + term2 + term3 + term4));

                  double Em = FGMMat[0];
                  double Ec = FGMMat[3];

                  double q   = 90.0e9;   // Stress transfer parameter
                  double Sym = 493.7e6; // Yield stress
                  double Sy  = Sym*((1 - VcSteps[lam]) + (q + Em)/(q + Ec)*(Ec/Em)*VcSteps[lam]);

                  fi = Svm/Sy;

                  if (Svm > maxFI) maxFI = fi;
                }
              }
            }
          }
        }

        cout << "maxFI = " << maxFI << endl;

      //cout << "                       \n buckfactor " << buckfactor << endl;
       if (buckfactor == 0)
       {
          cout << "Convergence not achieved in infill: " << endl;
          //exit(0);
       }

       // Push back the new targets Ybuck and Ystren
       lbdb = buckfactor;
      }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateCutOutFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 1;
}

// ========================== cSquarePlateBuckFGM ===========================

cShallowShellMFThermBuckFGM :: cShallowShellMFThermBuckFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
  CostMax = 0.7;
}

// ============================== Evaluate =================================

void cShallowShellMFThermBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // x[0] = 1.0; x[1] = 1.0; x[2] = 1.0; x[3] = 0.0; x[4] = 0.0; x[5] = 0.0; x[6] = 0.0; x[7] = 0.162; x[8] = 1.0;

  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj, 2);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation

  // Cost < Cmax
  // Cmax = 0.7
  double Cmax = CostMax;

  double vcratio; int numcp;

  numcp = NumVar;
  cVector Vcp(numcp);
  for (int i = 0; i < NumVar; i++){
      Vcp[i] = x[i];
  }

  EvalVolumeRatio(Vcp, vcratio);

  double R = 2.54; double t = 0.0127; double L = 0.508;
  double CostM = 1; double CostC = 20;
  double rm = R - t/2; double rM = R + t/2; double Volume = PI*(rM*rM - rm*rm)*L;

  double Cost = vcratio*Volume*CostC + (1 - vcratio)*Volume*CostM;

  c[0] = Cost - Cmax;
}

// ============================== Evaluate =================================

void cShallowShellMFThermBuckFGM :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    /*x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 0.4;
    x[3] = 0.0;
    x[4] = 0.0;*/

    /*x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    x[3] = 0.0;
    x[4] = 0.0;*/

    /*x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    x[3] = 0.45;
    x[4] = 0.0;*/

  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(x, fobj, 1);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!

  // Constraint evaluation

  // Cost < Cmax
  // Cmax = 0.7
  double Cmax = CostMax;

  double vcratio; int numcp;

  numcp = NumVar;
  cVector Vcp(numcp);
  for (int i = 0; i < NumVar; i++){
      Vcp[i] = x[i];
  }

  EvalVolumeRatio(Vcp, vcratio);

  double R = 2.54; double t = 0.0127; double L = 0.508;
  double CostM = 1; double CostC = 20;
  double rm = R - t/2; double rM = R + t/2; double Volume = PI*(rM*rM - rm*rm)*L;

  double Cost = vcratio*Volume*CostC + (1 - vcratio)*Volume*CostM;

  c[0] = Cost - Cmax;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cShallowShellMFThermBuckFGM :: Analysis(cVector x, double &lbdb, int f)
{
    // Define number of gauss points

        int ngauss = 10;

        // Evaluate weights and abscissas of gauss points

        cVector r, w;
        GaussPts1D(ngauss, r, w);

        // Transform abscissas to thickness coordinate

        cVector t(ngauss);

        // for (int i = 0; i < ngauss; i++) t[i]  = (r[i] + 1)/2;  // t = 0 (bottom) and t = 1 (top)
        for (int i = 0; i < ngauss; i++) t[i]  = r[i]/2;  // t = -0.5 (bottom) and t = 0.5 (top)

        // Evaluate volume fraction at gauss points according to a given distribution

        cVector Vcpg;

        int numcp;

        /*if ((NumVar)%2 == 0){
            numcp = (NumVar)*2;
        }
        else{
            numcp = 2*(NumVar) - 1;
        }

        cVector Vcp(numcp);
        for (int i = 0; i < NumVar; i++){
            Vcp[i] = x[i];
            Vcp[numcp - i - 1] = x[i];
        }*/

        numcp = NumVar;
        cVector Vcp(numcp);
        for (int i = 0; i < NumVar; i++){
            Vcp[i] = x[i];
        }

        double thk = 0.0127;

        /*
        // Heat Conduction

        double Ti = 1.0; double Ts = 0.5;
        int p  = 1; int ne = 10; int nt = p*ne + 1;

        cVector T(nt); cVector tHC(nt), VcHC(nt);
        T.Zero( ); tHC.Zero( ); VcHC.Zero( );
        tHC[0] = -0.5;
        for (int i = 1; i < nt; i++) tHC[i] = tHC[i - 1] + 1.0/((double)nt - 1.0);

        cVector FEMParam(2); FEMParam[0] = ne; FEMParam[1] = p;
        VolumeDist(FGMVolDist, Vcp, nt, tHC, VcHC);

        HeatConductionFEM(thk, Ti, Ts, nt, tHC, VcHC, FEMParam, T);
        */

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd, cmd2, cmd3, cmd4;
      if (f == 1)
      {
          cmd  = "del ShallowShellLF" + thread_number + ".dat";
          cmd2 = "del ShallowShellLF" + thread_number + ".pos";
          cmd3 = "rm ShallowShellLF" + thread_number + ".dat";
          cmd4 = "rm ShallowShellLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          cmd  = "del ShallowShellHF" + thread_number + ".dat";
          cmd2 = "del ShallowShellHF" + thread_number + ".pos";
          cmd3 = "rm ShallowShellHF" + thread_number + ".dat";
          cmd4 = "rm ShallowShellHF" + thread_number + ".pos";
      }


    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing DoLee.dat and plate.pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing DoLee.dat and plate.pos files.\n";
    #endif

      string begname, endname, datname, posname;
      int NumElm;

      if (f == 1)
      {
          begname = "datbegShallowShelln8p2.dat";
          endname = "datendShallowShelln8p2.dat";
          // NumElm = 4;
          // NumElm = 16;
          NumElm = 64;
          // NumElm = 4096;

          datname = "ShallowShellLF" + thread_number + ".dat";
          posname = "ShallowShellLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          begname = "datbegShallowShelln64p3.dat";
          endname = "datendShallowShelln64p3.dat";
          // NumElm = 4;
          // NumElm = 16;
          // NumElm = 64;
          NumElm = 4096;

          datname = "ShallowShellHF" + thread_number + ".dat";
          posname = "ShallowShellHF" + thread_number + ".pos";
      }

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

         //exit(EXIT_FAILURE);

         lbdb = 0.0;

         cout << "chega aqui" << endl;
         //exit(0);
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

      dat << endl << endl << "%SECTION.FGM.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1    " << thk << "    10    2    " << numcp;

      for (int i = 0; i < numcp; i++)
      {
          dat << "    " << Vcp[i];
      }
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
      //   exit(EXIT_FAILURE);

         lbdb = 0.0;
         exit(0);
      }

      /*// Write the element temperature

      fstream dat2;

      dat2.open(datname.c_str( ));

      if (!dat2.is_open( ))
      {
         cout << "Error opening the dat file for plate analysis." << endl;
         exit(0);
      }

      dat2.seekp(0,ofstream::end);

      dat2 << endl << endl << "%LOAD.CASE.ELEMENT.TEMPERATURE" << endl;
      dat2 << NumElm << endl;
      for (int i = 0; i < NumElm; i++)
      {
          dat2 << i+1 << "    " << p <<  "    " << nt << "    ";
          for (int j = 0; j < nt; j++) dat2 << T[j] << "    ";
          dat2 << endl;
      }
      dat2 << endl;
      dat2 << "%END";

      dat2.close( );
      */

      // Run the analysis with FAST.

    string fid;
    if (f == 1)
    {
        fid = "LF";
    }
    else
    {
        fid = "HF";
    }

    #ifdef _WIN32
      cmd = "fast.exe ShallowShell" + fid + thread_number + " -silent";
    #else
      cmd = "./fast ShallowShell" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
         //exit(EXIT_FAILURE);

    #ifdef _WIN32
      cmd = "fast.exe ShallowShell" + fid + thread_number + " -silent";
    #else
      cmd = "./fast ShallowShell" + fid + thread_number + " -silent";
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
      //string posname = "platehole.pos";
      //pos.open(posname.c_str());
      if (!pos.is_open( ))
      {
         cout << "Error opening the pos file for plate analysis." << endl;
         exit(0);
      }

      // Find buckling load factor

      string label;
      double buckfactor = 0;
    /*  cVector genstress(6);
      cVector genstrain(6);
      cVector force(8);
      force.Zero();*/
      int mode;

      while (pos >> label)
      {

          if (label == "%RESULT.CASE.STEP.NATURAL.FREQUENCY")
          {
              pos >> mode;
              pos >> buckfactor;
          }
       }

      //cout << "                       \n buckfactor " << buckfactor << endl;
       if (buckfactor == 0)
       {
          cout << "Convergence not achieved in infill: " << endl;
          //exit(0);
       }

       // Push back the new targets Ybuck and Ystren
       lbdb = buckfactor;
      }
}

// ============================ Evaluate ==============================

void cShallowShellMFThermBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{

    // Constraint evaluation

    // Cost < Cmax
    // Cmax = 0.7
    double Cmax = CostMax;

    double vcratio; int numcp;

    numcp = NumVar;
    cVector Vcp(numcp);
    for (int i = 0; i < NumVar; i++){
        Vcp[i] = x[i];
    }

    EvalVolumeRatio(Vcp, vcratio);

    double R = 2.54; double t = 0.0127; double L = 0.508;
    double CostM = 1; double CostC = 20;
    double rm = R - t/2; double rM = R + t/2; double Volume = PI*(rM*rM - rm*rm)*L;

    double Cost = vcratio*Volume*CostC + (1 - vcratio)*Volume*CostM;

  // Single constraint evaluation.
    if (index == 0){
        c = Cost - Cmax;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cShallowShellMFThermBuckFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateMFBuckVSC :: cSquarePlateMFBuckVSC(void)
{
  NumConstr = 0;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateMFBuckVSC :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
  // x[0] = 0.253556; x[1] = 0.735111;
  // x[2] = 0.253556; x[3] = 0.735111;
  // x[4] = 0.253556; x[5] = 0.735111;
  // x[6] = 0.253556; x[7] = 0.735111;

  // Decode variable vector

  // x.Print( );

  int NumLam = 2*NumVar;

  cMatrix Ang(NumLam, 2);
  Ang.Zero( );

  for (int i = 0; i < x.Dim( )/2; i++)
  {
      Ang[2*i][0]   = 90.0*x[2*i];
      Ang[2*i][1]   = 90.0*x[2*i + 1];

      Ang[2*i + 1][0]   = -90.0*x[2*i];
      Ang[2*i + 1][1]   = -90.0*x[2*i + 1];

      Ang[NumLam - (2*i+1) - 1][0] = -90.0*x[2*i];
      Ang[NumLam - (2*i+1) - 1][1] = -90.0*x[2*i + 1];

      Ang[NumLam - (2*i) - 1][0] = 90.0*x[2*i];
      Ang[NumLam - (2*i) - 1][1] = 90.0*x[2*i + 1];
  }

  // cout << "Ang = " << endl;
  // Ang.Print( );

  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(Ang, fobj, 2);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;      // Maximization problem!
  // cout << "fobj = " << fobj << endl;
}

// ============================== Evaluate =================================

void cSquarePlateMFBuckVSC :: EvaluateLFP(cVector &x, cVector &c, cVector &fobjs)
{
    // Decode variable vector

    int NumLam = 2*NumVar;

    cMatrix Ang(NumLam, 2);
    Ang.Zero( );

    for (int i = 0; i < x.Length( ); i++)
    {
        Ang[2*i][0]   = 90.0*x[i];
        Ang[2*i+1][1] = 90.0*x[i];

        Ang[NumLam - (2*i+1) - 1][1] = 90.0*x[i];
        Ang[NumLam - (2*i) - 1][0]   = 90.0*x[i];
    }
  // Objective function evaluation (using FAST)

  double fobj;
  Analysis(Ang, fobj, 1);  // Linearized Buckling Analysis
  fobjs[0] = -fobj;   // Maximization problem!
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateMFBuckVSC :: Analysis(cMatrix Ang, double &lbdb, int f)
{
        double thklam = 0.000127;

        /*cout << "\n\n ============ Matriz A ======= " << endl;
        A.Print();
        cout << "\n ============ Matriz B ======= " << endl;
        B.Print();
        cout << "\n ============ Matriz D ======= " << endl;
        D.Print();
        cout << "\n ============ Matriz G ======= " << endl;
        G.Print();
        cout << "\n ============ Matriz ABDG ======= " << endl;
        ABDG.Print();*/

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd, cmd2, cmd3, cmd4;
      if (f == 1)
      {
          cmd  = "del GuoVSCLF" + thread_number + ".dat";
          cmd2 = "del GuoVSCLF" + thread_number + ".pos";
          cmd3 = "rm GuoVSCLF" + thread_number + ".dat";
          cmd4 = "rm GuoVSCLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          cmd  = "del GuoVSCHF" + thread_number + ".dat";
          cmd2 = "del GuoVSCHF" + thread_number + ".pos";
          cmd3 = "rm GuoVSCHF" + thread_number + ".dat";
          cmd4 = "rm GuoVSCHF" + thread_number + ".pos";
      }


    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing DoLee.dat and plate.pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing DoLee.dat and plate.pos files.\n";
    #endif

      string begname, endname, datname, posname;

      if (f == 1)
      {
          begname = "datbegGuoVSCLF.dat";
          endname = "datendGuoVSCLF.dat";

          datname = "GuoVSCLF" + thread_number + ".dat";
          posname = "GuoVSCLF" + thread_number + ".pos";
      }
      else if (f == 2)
      {
          begname = "datbegGuoVSCHF.dat";
          endname = "datendGuoVSCHF.dat";

          datname = "GuoVSCHF" + thread_number + ".dat";
          posname = "GuoVSCHF" + thread_number + ".pos";
      }

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

         //exit(EXIT_FAILURE);

         lbdb = 0.0;

         cout << "chega aqui" << endl;
         //exit(0);
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

      dat << endl << endl << "%SECTION.VSC.SHELL" << endl;
      dat << "1" << endl;
      dat << "1    1.0    0.0    0.0    1    0.1    16 " << endl;

      for (int i = 0; i < 2*NumVar; i++)
      {
          dat << "1    " << thklam << "    " << Ang[i][0] << "    " << Ang[i][1] << endl;
      }

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
      //   exit(EXIT_FAILURE);

         lbdb = 0.0;
         exit(0);
      }

      // Run the analysis with FAST.

    string fid;
    if (f == 1)
    {
        fid = "LF";
    }
    else
    {
        fid = "HF";
    }

    #ifdef _WIN32
      cmd = "fast.exe GuoVSC" + fid + thread_number + " -silent";
    #else
      cmd = "./fast GuoVSC" + fid + thread_number + " -silent";
    #endif

      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
         //exit(EXIT_FAILURE);

    #ifdef _WIN32
      cmd = "fast.exe DoLee" + fid + thread_number + " -silent";
    #else
      cmd = "./fast DoLee" + fid + thread_number + " -silent";
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
      //string posname = "platehole.pos";
      //pos.open(posname.c_str());
      if (!pos.is_open( ))
      {
         cout << "Error opening the pos file for plate analysis." << endl;
         exit(0);
      }

      // Find buckling load factor

      string label;
      double buckfactor = 0;
    /*  cVector genstress(6);
      cVector genstrain(6);
      cVector force(8);
      force.Zero();*/
      int mode;

      while (pos >> label)
      {

          if (label == "%RESULT.CASE.STEP.BUCKLING.FACTOR")
          {
              pos >> mode;
              pos >> buckfactor;
          }
       }

      //cout << "                       \n buckfactor " << buckfactor << endl;
       if (buckfactor == 0)
       {
          cout << "Convergence not achieved in infill: " << endl;
          //exit(0);
       }

       // Push back the new targets Ybuck and Ystren
       lbdb = buckfactor;
      }
}

// ======================================================= End of file =====
