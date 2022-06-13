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
  cProblemFactory :: Register("SquarePlateTridirBuckFGM" , MakeProb<cSquarePlateTridirBuckFGM> ,".fgm"),
  cProblemFactory :: Register("ShellTridirBuckFGM"       , MakeProb<cShellTridirBuckFGM>       ,".fgm"),
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

// -------------------------------------------------------------------------
// Public methods:
//

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateTridirBuckFGM :: cSquarePlateTridirBuckFGM(void)
{
  NumConstr = 1;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateTridirBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
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

    // Objective function evaluation (using FAST)

    double fobj;
    Analysis(x, fobj);  // Linearized Buckling Analysis
    fobjs[0] = -fobj;   // Maximization problem!

    // Constraint evaluation

    // Ceramic volume percentage < Cmax
    // Cmax = 50%

    double vcratio;
    double Cmax = 0.30;

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

void cSquarePlateTridirBuckFGM :: Analysis(cVector x, double &lbdb)
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
      fid = "2D";


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

      begname = "datbegDoLeeTridir16x16.dat";
      endname = "datendDoLeeTridir16x16.dat";

      datname = "SqrPltBuck2D" + thread_number + ".dat";
      posname = "SqrPltBuck2D" + thread_number + ".pos";

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
      dat << "1    1    " << thk << "    10    3    " << numcp+10 << endl;
      dat << "-5.0  5.0  -5.0  5.0" << endl; // lxlow; lxupp; lylow; lyupp;
      dat << "6  6  4" << endl; // ncp_x; ncp_y; ncp_z;
      dat << "3  3  3" << endl; // Cubic in all coordinates
      for (int i = 0; i < numcp/6; i++) dat << Vcp[i*6 + 0] << "  " << Vcp[i*6 + 1] << "  " << Vcp[i*6 + 2] << "  " << Vcp[i*6 + 3] << "  " << Vcp[i*6 + 4] << "  " << Vcp[i*6 + 5] << endl;
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

void cSquarePlateTridirBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    // Ceramic volume percentage < Cmax
    // Cmax = 50%

    double vcratio;
    double Cmax = 0.30;

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

void cSquarePlateTridirBuckFGM :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadW ======================================

void cShellTridirBuckFGM :: ReadW(istream &in)
{
  if (!(in >> W_MObj))
  {
    cout << "Error in the input of the laminate maximum number of plies." << endl;
    exit(0);
  }
}

// ============================== LoadReadFunc =============================

void cShellTridirBuckFGM :: LoadReadFunc(cInpMap &im)
{
  // Load base class functions.
  cFGM :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("MULTIOBJECTIVE.WEIGHT",makeReadObj(cShellTridirBuckFGM,ReadW));
}

// ========================== cShellTridirBuckFGM ==========================

cShellTridirBuckFGM :: cShellTridirBuckFGM(void)
{
  NumConstr = 0;
  NumObj    = 1;
  W_MObj    = 1;
}

// ============================== Evaluate =================================

void cShellTridirBuckFGM :: Evaluate(cVector &x, cVector &c, cVector &fobjs)
{
// Test solutions.
/*
    for (int i = 0; i < 32; i++)
        x[i] = 1.0;
*/

    // w = 0.6, Vc = 33%
    /*
    x[0] = 1;
    x[1] = 1;
    x[2] = 1;
    x[3] = 1;
    x[4] = 0;
    x[5] = 1;
    x[6] = 0.8384;
    x[7] = 0;
    x[8] = 0;
    x[9] = 1;
    x[10] = 0;
    x[11] = 0;
    x[12] = 1;
    x[13] = 1;
    x[14] = 1;
    x[15] = 1;
    x[16] = 0;
    x[17] = 0;
    x[18] = 0;
    x[19] = 0;
    x[20] = 0;
    x[21] = 0;
    x[22] = 0;
    x[23] = 0;
    x[24] = 0;
    x[25] = 0;
    x[26] = 0;
    x[27] = 0;
    x[28] = 0;
    x[29] = 0;
    x[30] = 0;
    x[31] = 0;


    // w = 0.3, Vc = 70%
    x[0] = 1;
    x[1] = 1;
    x[2] = 1;
    x[3] = 1;
    x[4] = 1;
    x[5] = 1;
    x[6] = 1;
    x[7] = 1;
    x[8] = 1;
    x[9] = 1;
    x[10] = 1;
    x[11] = 1;
    x[12] = 1;
    x[13] = 1;
    x[14] = 1;
    x[15] = 1;
    x[16] = 1;
    x[17] = 1;
    x[18] = 1;
    x[19] = 1;
    x[20] = 0;
    x[21] = 1;
    x[22] = 0;
    x[23] = 0;
    x[24] = 0;
    x[25] = 0.2792;
    x[26] = 0;
    x[27] = 0;
    x[28] = 1;
    x[29] = 1;
    x[30] = 1;
    x[31] = 1;

    // w = 0.5, Vc = 50%
    x[0] = 1;
    x[1] = 1;
    x[2] = 1;
    x[3] = 1;
    x[4] = 1;
    x[5] = 1;
    x[6] = 1;
    x[7] = 1;
    x[8] = 1;
    x[9] = 1;
    x[10] = 1;
    x[11] = 1;
    x[12] = 1;
    x[13] = 1;
    x[14] = 1;
    x[15] = 1;
    x[16] = 0;
    x[17] = 0;
    x[18] = 0;
    x[19] = 0;
    x[20] = 0;
    x[21] = 0;
    x[22] = 0;
    x[23] = 0;
    x[24] = 0;
    x[25] = 0;
    x[26] = 0;
    x[27] = 0;
    x[28] = 0;
    x[29] = 0;
    x[30] = 0;
    x[31] = 0;
   
   */

    // Objective function evaluation (using FAST)

    double fobj;
    Analysis(x, fobj);  // Linearized Buckling Analysis

    // Normalized buckling load

    double a  = 1, pi = atan(1.0)*4.0;
    double Ec = 348.43e9, nuc = 0.24;
    double h  = 0.02;
    double Dc = Ec*h*h*h/(12*(1 - nuc*nuc));

    double lbdn = fobj*a*a/(pi*pi*Dc);

    // Evaluation of the Ceramic Volume Fraction

    double vcratio;
    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    int Num0[8]  = {0 , 7 , 56 , 63 , 192, 199, 248, 255};
    int Num1[8]  = {1 , 6 , 57 , 62 , 193, 198, 249, 254};
    int Num2[8]  = {2 , 5 , 58 , 61 , 194, 197, 250, 253};
    int Num3[8]  = {3 , 4 , 59 , 60 , 195, 196, 251, 252};
    int Num4[8]  = {8 , 15, 48 , 55 , 200, 207, 240, 247};
    int Num5[8]  = {9 , 14, 49 , 54 , 201, 206, 241, 246};
    int Num6[8]  = {10, 13, 50 , 53 , 202, 205, 242, 245};
    int Num7[8]  = {11, 12, 51 , 52 , 203, 204, 243, 244};
    int Num8[8]  = {16, 23, 40 , 47 , 208, 215, 232, 239};
    int Num9[8]  = {17, 22, 41 , 46 , 209, 214, 233, 238};
    int Num10[8] = {18, 21, 42 , 45 , 210, 213, 234, 237};
    int Num11[8] = {19, 20, 43 , 44 , 211, 212, 235, 236};
    int Num12[8] = {24, 31, 32 , 39 , 216, 223, 224, 231};
    int Num13[8] = {25, 30, 33 , 38 , 217, 222, 225, 230};
    int Num14[8] = {26, 29, 34 , 37 , 218, 221, 226, 229};
    int Num15[8] = {27, 28, 35 , 36 , 219, 220, 227, 228};
    int Num16[8] = {64, 71, 120, 127, 128, 135, 184, 191};
    int Num17[8] = {65, 70, 121, 126, 129, 134, 185, 190};
    int Num18[8] = {66, 69, 122, 125, 130, 133, 186, 189};
    int Num19[8] = {67, 68, 123, 124, 131, 132, 187, 188};
    int Num20[8] = {72, 79, 112, 119, 136, 143, 176, 183};
    int Num21[8] = {73, 78, 113, 118, 137, 142, 177, 182};
    int Num22[8] = {74, 77, 114, 117, 138, 141, 178, 181};
    int Num23[8] = {75, 76, 115, 116, 139, 140, 179, 180};
    int Num24[8] = {80, 87, 104, 111, 144, 151, 168, 175};
    int Num25[8] = {81, 86, 105, 110, 145, 150, 169, 174};
    int Num26[8] = {82, 85, 106, 109, 146, 149, 170, 173};
    int Num27[8] = {83, 84, 107, 108, 147, 148, 171, 172};
    int Num28[8] = {88, 95, 96 , 103, 152, 159, 160, 167};
    int Num29[8] = {89, 94, 97 , 102, 153, 158, 161, 166};
    int Num30[8] = {90, 93, 98 , 101, 154, 157, 162, 165};
    int Num31[8] = {91, 92, 99 , 100, 155, 156, 163, 164};

    for (int i = 0; i < 8; i++)
    {
        Vcp[Num0[i]]  = x[0];
        Vcp[Num1[i]]  = x[1];
        Vcp[Num2[i]]  = x[2];
        Vcp[Num3[i]]  = x[3];
        Vcp[Num4[i]]  = x[4];
        Vcp[Num5[i]]  = x[5];
        Vcp[Num6[i]]  = x[6];
        Vcp[Num7[i]]  = x[7];
        Vcp[Num8[i]]  = x[8];
        Vcp[Num9[i]]  = x[9];
        Vcp[Num10[i]] = x[10];
        Vcp[Num11[i]] = x[11];
        Vcp[Num12[i]] = x[12];
        Vcp[Num13[i]] = x[13];
        Vcp[Num14[i]] = x[14];
        Vcp[Num15[i]] = x[15];
        Vcp[Num16[i]] = x[16];
        Vcp[Num17[i]] = x[17];
        Vcp[Num18[i]] = x[18];
        Vcp[Num19[i]] = x[19];
        Vcp[Num20[i]] = x[20];
        Vcp[Num21[i]] = x[21];
        Vcp[Num22[i]] = x[22];
        Vcp[Num23[i]] = x[23];
        Vcp[Num24[i]] = x[24];
        Vcp[Num25[i]] = x[25];
        Vcp[Num26[i]] = x[26];
        Vcp[Num27[i]] = x[27];
        Vcp[Num28[i]] = x[28];
        Vcp[Num29[i]] = x[29];
        Vcp[Num30[i]] = x[30];
        Vcp[Num31[i]] = x[31];
    }

    EvalVolumeRatio3D(Vcp, vcratio, 8, 8, 4);

    cVector center(2); center[0] = 0.0; center[1] = 0.0; // parametric coordinates
    double radius = 0.2;                                 // parametric coordinates
    double vcratiohole;

    EvalVolumeHole(center, radius, Vcp, vcratiohole);    // center of the hole ; radius ; Vcp ; vcratiohole

    // Plate volume

    double Area   = 1.0;
    double thk    = 0.02;
    double Vplate = Area*thk;

    // Hole volume

    double rad      = 0.1;
    double AreaHole = pi*rad*rad;
    double Vhole    = AreaHole*thk;

    // Total volume

    double Vol   = Vplate - Vhole;
    double VcTot = vcratio*Vplate - vcratiohole*Vhole;
    double Vcrtot = VcTot/Vol;

    // Cost

    // 0.966920

    double rhom = 8000, rhoc = 2730;
    double Cm   = 3   , Cc   = 50;

    double TotalCost = 0.966920*thk*(Vcrtot*rhoc*Cc + (1 - Vcrtot)*rhom*Cm);

    // buck: lbdn ; cost: TotalCost

    // Stores the value for each objective function

    double w = W_MObj;
    double m = 2.0;

    double CostMin = 464.1216;
    double CostMax = 2639.6916;
    double BuckMin = -2.140551;
    double BuckMax = -1.169951;

    fobjs[0] = pow(w * (TotalCost-CostMin)/(CostMax-CostMin),m) + pow( (1.0-w)*(-lbdn-BuckMin)/(BuckMax-BuckMin),m);
    
    // Feedback
    //cout << "fobj: " <<fobjs[0] << endl;
    //cout <<setprecision(10);
    //cout << scientific << setprecision(6);
    //cout << "Total Cost " << TotalCost << endl;
    //cout << "-lbdn " << -lbdn << endl;
}

// ============================ Evaluate ==============================

void cShellTridirBuckFGM :: EvalVolumeHole(cVector c, double r, cVector Vcp, double &v)
{
    int Nbx = 8; int Nby = 8; int Nbz = 4;

    cMatrix *CP;

    CP = new cMatrix[Nbz];

    for (int m = 0; m < Nbz; m++){
        cMatrix auxM(Nbx,Nby);
        for (int j = 0; j < Nby; j++){
            for (int i = 0; i < Nbx; i++)
            {
              auxM[j][i] = Vcp[m*Nbx*Nby + j*Nbx + i];
            }
        }
        CP[m].Resize(Nbx, Nby);
        CP[m] = auxM;
    }

    int px, py, pz;
    px = py = pz = 3;

    int mx = Nbx + px + 1;
    int my = Nby + py + 1;
    int mz = Nbz + pz + 1;

    // Knot vectors

    double lxupp, lxlow, lyupp, lylow, lzupp, lzlow;

    lxupp = lyupp = lzupp =  1.0;
    lxlow = lylow = lzlow = -1.0;

    // x axis
    cVector Ux( mx );
    for (int i = 0; i < mx; i++){
        if (i < (px + 1)) Ux[i] = 0;
        else if (i >= Nbx) Ux[i] = 1;
        else{
            Ux[i] = Ux[i - 1] + 1.0/(Nbx - px);
        }
    }

    for (int i = 0; i < mx; i++){
        Ux[i] = Ux[i]*(lxupp - lxlow) + lxlow; // Knot vector defined between lxlow and lxupp
    }

    // y axis
    cVector Uy( my );
    for (int i = 0; i < my; i++){
        if (i < (py + 1)) Uy[i] = 0;
        else if (i >= Nby) Uy[i] = 1;
        else{
            Uy[i] = Uy[i - 1] + 1.0/(Nby - py);
        }
    }

    for (int i = 0; i < my; i++){
        Uy[i] = Uy[i]*(lyupp - lylow) + lylow; // Knot vector defined between lylow and lyupp
    }

    // z axis
    cVector Uz( mz );
    for (int i = 0; i < mz; i++){
        if (i < (pz + 1)) Uz[i] = 0;
        else if (i >= Nbz) Uz[i] = 1;
        else{
            Uz[i] = Uz[i - 1] + 1.0/(Nbz - pz); // Knot vector defined between 0 and 1
        }
    }

    for (int i = 0; i < mz; i++){
        Uz[i] = Uz[i]*(lzupp - lzlow) + lzlow; // Knot vector defined between lzlow and lzupp
    }

    int n1 = 5; int n2 = 5; int n3 = 5;
    //int n1 = 10; int n2 = 10; int n3 = 10;
    cVector x1(n1), x2(n2), x3(n3);

    for (int i = 0; i < n1; i++)
    {
        x1[i] = c[0] - r*((double)n1 - i - 1.0)/((double)n1 - 1.0);
        x2[i] = c[1] - r*((double)n1 - i - 1.0)/((double)n1 - 1.0);
        x3[i] =      1.0*((double)n1 - i - 1.0)/((double)n1 - 1.0);
    }

    cVector V2(1);
    int ContCircle = 0; v = 0;

    for (int i = 0; i < (n1 - 1); i++)
    {
        for (int j = 0; j < (n2 - 1); j++)
        {
            for (int k = 0; k < (n3 - 1); k++)
            {
                cMatrix cdnt(3,1);
                cdnt[0][0] = (x1[i] + x1[i + 1])/2.0;
                cdnt[1][0] = (x2[j] + x2[j + 1])/2.0;
                cdnt[2][0] = (x3[k] + x3[k + 1])/2.0;

                // cout << cdnt[0][0] << "  "  << cdnt[1][0] << "  "  << cdnt[2][0] << endl;

                // check if inside the circle

                double dist = pow(cdnt[0][0]*cdnt[0][0] + cdnt[1][0]*cdnt[1][0], 0.5);

                // if (dist > r) cout << "NOT ON CIRCLE" << endl;
                if (dist > r) continue;

                // evaluate volume fraction
                //cout << "in hole" << endl;

                BsplineSol(CP, 1, cdnt, V2, Nbx, Nby, Nbz, px, py, pz, Ux, Uy, Uz);
                v += V2[0];
                ContCircle += 1;
            }
        }
    }

    // cout << endl;
    // cout << "ContCircle = " << ContCircle << endl;

    v = v/(double)ContCircle;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cShellTridirBuckFGM :: Analysis(cVector x, double &lbdb)
{
    // Evaluate volume fraction at gauss points according to a given distribution

    cVector Vcpg;

    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    // Filling control points vector

    int Num0[8]  = {0 , 7 , 56 , 63 , 192, 199, 248, 255};
    int Num1[8]  = {1 , 6 , 57 , 62 , 193, 198, 249, 254};
    int Num2[8]  = {2 , 5 , 58 , 61 , 194, 197, 250, 253};
    int Num3[8]  = {3 , 4 , 59 , 60 , 195, 196, 251, 252};
    int Num4[8]  = {8 , 15, 48 , 55 , 200, 207, 240, 247};
    int Num5[8]  = {9 , 14, 49 , 54 , 201, 206, 241, 246};
    int Num6[8]  = {10, 13, 50 , 53 , 202, 205, 242, 245};
    int Num7[8]  = {11, 12, 51 , 52 , 203, 204, 243, 244};
    int Num8[8]  = {16, 23, 40 , 47 , 208, 215, 232, 239};
    int Num9[8]  = {17, 22, 41 , 46 , 209, 214, 233, 238};
    int Num10[8] = {18, 21, 42 , 45 , 210, 213, 234, 237};
    int Num11[8] = {19, 20, 43 , 44 , 211, 212, 235, 236};
    int Num12[8] = {24, 31, 32 , 39 , 216, 223, 224, 231};
    int Num13[8] = {25, 30, 33 , 38 , 217, 222, 225, 230};
    int Num14[8] = {26, 29, 34 , 37 , 218, 221, 226, 229};
    int Num15[8] = {27, 28, 35 , 36 , 219, 220, 227, 228};
    int Num16[8] = {64, 71, 120, 127, 128, 135, 184, 191};
    int Num17[8] = {65, 70, 121, 126, 129, 134, 185, 190};
    int Num18[8] = {66, 69, 122, 125, 130, 133, 186, 189};
    int Num19[8] = {67, 68, 123, 124, 131, 132, 187, 188};
    int Num20[8] = {72, 79, 112, 119, 136, 143, 176, 183};
    int Num21[8] = {73, 78, 113, 118, 137, 142, 177, 182};
    int Num22[8] = {74, 77, 114, 117, 138, 141, 178, 181};
    int Num23[8] = {75, 76, 115, 116, 139, 140, 179, 180};
    int Num24[8] = {80, 87, 104, 111, 144, 151, 168, 175};
    int Num25[8] = {81, 86, 105, 110, 145, 150, 169, 174};
    int Num26[8] = {82, 85, 106, 109, 146, 149, 170, 173};
    int Num27[8] = {83, 84, 107, 108, 147, 148, 171, 172};
    int Num28[8] = {88, 95, 96 , 103, 152, 159, 160, 167};
    int Num29[8] = {89, 94, 97 , 102, 153, 158, 161, 166};
    int Num30[8] = {90, 93, 98 , 101, 154, 157, 162, 165};
    int Num31[8] = {91, 92, 99 , 100, 155, 156, 163, 164};

    for (int i = 0; i < 8; i++)
    {
        Vcp[Num0[i]]  = x[0];
        Vcp[Num1[i]]  = x[1];
        Vcp[Num2[i]]  = x[2];
        Vcp[Num3[i]]  = x[3];
        Vcp[Num4[i]]  = x[4];
        Vcp[Num5[i]]  = x[5];
        Vcp[Num6[i]]  = x[6];
        Vcp[Num7[i]]  = x[7];
        Vcp[Num8[i]]  = x[8];
        Vcp[Num9[i]]  = x[9];
        Vcp[Num10[i]] = x[10];
        Vcp[Num11[i]] = x[11];
        Vcp[Num12[i]] = x[12];
        Vcp[Num13[i]] = x[13];
        Vcp[Num14[i]] = x[14];
        Vcp[Num15[i]] = x[15];
        Vcp[Num16[i]] = x[16];
        Vcp[Num17[i]] = x[17];
        Vcp[Num18[i]] = x[18];
        Vcp[Num19[i]] = x[19];
        Vcp[Num20[i]] = x[20];
        Vcp[Num21[i]] = x[21];
        Vcp[Num22[i]] = x[22];
        Vcp[Num23[i]] = x[23];
        Vcp[Num24[i]] = x[24];
        Vcp[Num25[i]] = x[25];
        Vcp[Num26[i]] = x[26];
        Vcp[Num27[i]] = x[27];
        Vcp[Num28[i]] = x[28];
        Vcp[Num29[i]] = x[29];
        Vcp[Num30[i]] = x[30];
        Vcp[Num31[i]] = x[31];
    }

    double thk = 0.02;

      int num_thread = 0;
    #ifdef _OMP_
      num_thread = omp_get_thread_num( );
    #endif

      stringstream thread;
      thread << num_thread;

      string thread_number = thread.str();
      string cmd  = "del ShellBuck" + thread_number + ".dat";
      string cmd2 = "del ShellBuck" + thread_number + ".pos";
      string cmd3 = "rm ShellBuck" + thread_number + ".dat";
      string cmd4 = "rm ShellBuck" + thread_number + ".pos";

    #ifdef _WIN32
      if (system(cmd.c_str()) || system(cmd2.c_str()))
         cout << "Problem on removing SqrPltBuck.dat and plate.pos files.\n";
    #else
      if (system(cmd3.c_str()) || system(cmd4.c_str()))
         cout << "Problem on removing SqrPltBuck.dat and plate.pos files.\n";
    #endif

      string begname, endname, datname, posname;

      begname = "datbegShellHoleTridir.txt";
      endname = "datendShellHoleTridir.txt";

      datname = "ShellBuck" + thread_number + ".dat";
      posname = "ShellBuck" + thread_number + ".pos";

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
      dat << "1    1    " << thk << "    10    3    " << numcp+10 << endl;
      dat << "0  1  -0.5  0.5" << endl;                                // lxlow; lxupp; lylow; lyupp;
      dat << "8  8  4" << endl;                                             // ncp_x; ncp_y; ncp_z;
      dat << "3  3  3" << endl;                                             // Cubic in all coordinates
      for (int i = 0; i < numcp/8; i++) dat << Vcp[i*8 + 0] << "  " << Vcp[i*8 + 1] << "  " << Vcp[i*8 + 2] << "  " << Vcp[i*8 + 3] << "  " << Vcp[i*8 + 4] << "  " << Vcp[i*8 + 5] << "  " << Vcp[i*8 + 6] << "  " << Vcp[i*8 + 7] << endl;
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
        cmd = "./fast ShellBuck" + thread_number + " -silent";
      #endif
      status2 = system(cmd.c_str( ));

      if (status2)
      {
         cout << "Error in the analysis with fast.";
        #ifdef _WIN32
          cmd = "fast.exe SqrPltBuck" + thread_number + " -silent";
        #else
          cmd = "./fast ShellBuck" + thread_number + " -silent";
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
      pos.close( );

      if (buckfactor == 0)
      {
         cout << "Convergence not achieved in infill: " << endl;
      }

      // Push back the new targets Ybuck and Ystren
      lbdb = buckfactor;
    }
}

// ============================ Evaluate ==============================

void cShellTridirBuckFGM :: EvalExactConstraint(int index, cVector& x, double &c)
{
    // Ceramic volume percentage < Cmax
    // Cmax = 50%

    double vcratio;
    double Cmax = 0.30;

    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    int Num0[8]  = {0 , 7 , 56 , 63 , 192, 199, 248, 255};
    int Num1[8]  = {1 , 6 , 57 , 62 , 193, 198, 249, 254};
    int Num2[8]  = {2 , 5 , 58 , 61 , 194, 197, 250, 253};
    int Num3[8]  = {3 , 4 , 59 , 60 , 195, 196, 251, 252};
    int Num4[8]  = {8 , 15, 48 , 55 , 200, 207, 240, 247};
    int Num5[8]  = {9 , 14, 49 , 54 , 201, 206, 241, 246};
    int Num6[8]  = {10, 13, 50 , 53 , 202, 205, 242, 245};
    int Num7[8]  = {11, 12, 51 , 52 , 203, 204, 243, 244};
    int Num8[8]  = {16, 23, 40 , 47 , 208, 215, 232, 239};
    int Num9[8]  = {17, 22, 41 , 46 , 209, 214, 233, 238};
    int Num10[8] = {18, 21, 42 , 45 , 210, 213, 234, 237};
    int Num11[8] = {19, 20, 43 , 44 , 211, 212, 235, 236};
    int Num12[8] = {24, 31, 32 , 39 , 216, 223, 224, 231};
    int Num13[8] = {25, 30, 33 , 38 , 217, 222, 225, 230};
    int Num14[8] = {26, 29, 34 , 37 , 218, 221, 226, 229};
    int Num15[8] = {27, 28, 35 , 36 , 219, 220, 227, 228};
    int Num16[8] = {64, 71, 120, 127, 128, 135, 184, 191};
    int Num17[8] = {65, 70, 121, 126, 129, 134, 185, 190};
    int Num18[8] = {66, 69, 122, 125, 130, 133, 186, 189};
    int Num19[8] = {67, 68, 123, 124, 131, 132, 187, 188};
    int Num20[8] = {72, 79, 112, 119, 136, 143, 176, 183};
    int Num21[8] = {73, 78, 113, 118, 137, 142, 177, 182};
    int Num22[8] = {74, 77, 114, 117, 138, 141, 178, 181};
    int Num23[8] = {75, 76, 115, 116, 139, 140, 179, 180};
    int Num24[8] = {80, 87, 104, 111, 144, 151, 168, 175};
    int Num25[8] = {81, 86, 105, 110, 145, 150, 169, 174};
    int Num26[8] = {82, 85, 106, 109, 146, 149, 170, 173};
    int Num27[8] = {83, 84, 107, 108, 147, 148, 171, 172};
    int Num28[8] = {88, 95, 96 , 103, 152, 159, 160, 167};
    int Num29[8] = {89, 94, 97 , 102, 153, 158, 161, 166};
    int Num30[8] = {90, 93, 98 , 101, 154, 157, 162, 165};
    int Num31[8] = {91, 92, 99 , 100, 155, 156, 163, 164};

    for (int i = 0; i < 8; i++)
    {
        Vcp[Num0[i]]  = x[0];
        Vcp[Num1[i]]  = x[1];
        Vcp[Num2[i]]  = x[2];
        Vcp[Num3[i]]  = x[3];
        Vcp[Num4[i]]  = x[4];
        Vcp[Num5[i]]  = x[5];
        Vcp[Num6[i]]  = x[6];
        Vcp[Num7[i]]  = x[7];
        Vcp[Num8[i]]  = x[8];
        Vcp[Num9[i]]  = x[9];
        Vcp[Num10[i]] = x[10];
        Vcp[Num11[i]] = x[11];
        Vcp[Num12[i]] = x[12];
        Vcp[Num13[i]] = x[13];
        Vcp[Num14[i]] = x[14];
        Vcp[Num15[i]] = x[15];
        Vcp[Num16[i]] = x[16];
        Vcp[Num17[i]] = x[17];
        Vcp[Num18[i]] = x[18];
        Vcp[Num19[i]] = x[19];
        Vcp[Num20[i]] = x[20];
        Vcp[Num21[i]] = x[21];
        Vcp[Num22[i]] = x[22];
        Vcp[Num23[i]] = x[23];
        Vcp[Num24[i]] = x[24];
        Vcp[Num25[i]] = x[25];
        Vcp[Num26[i]] = x[26];
        Vcp[Num27[i]] = x[27];
        Vcp[Num28[i]] = x[28];
        Vcp[Num29[i]] = x[29];
        Vcp[Num30[i]] = x[30];
        Vcp[Num31[i]] = x[31];
    }

    //cout << "1" << endl;
    EvalVolumeRatio3D(Vcp, vcratio, 8, 8, 4);
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

void cShellTridirBuckFGM :: GetApproxConstr(bool *approxc)
{

}

// ========================= Write ====================================

void cShellTridirBuckFGM :: Write(cVector &x, ostream &out)
{
    // Evaluation of the Ceramic Volume Fraction

    double vcratio;
    int numcp = NumVar*8;

    cVector Vcp(numcp);

    // Filling control points vector

    int Num0[8]  = {0 , 7 , 56 , 63 , 192, 199, 248, 255};
    int Num1[8]  = {1 , 6 , 57 , 62 , 193, 198, 249, 254};
    int Num2[8]  = {2 , 5 , 58 , 61 , 194, 197, 250, 253};
    int Num3[8]  = {3 , 4 , 59 , 60 , 195, 196, 251, 252};
    int Num4[8]  = {8 , 15, 48 , 55 , 200, 207, 240, 247};
    int Num5[8]  = {9 , 14, 49 , 54 , 201, 206, 241, 246};
    int Num6[8]  = {10, 13, 50 , 53 , 202, 205, 242, 245};
    int Num7[8]  = {11, 12, 51 , 52 , 203, 204, 243, 244};
    int Num8[8]  = {16, 23, 40 , 47 , 208, 215, 232, 239};
    int Num9[8]  = {17, 22, 41 , 46 , 209, 214, 233, 238};
    int Num10[8] = {18, 21, 42 , 45 , 210, 213, 234, 237};
    int Num11[8] = {19, 20, 43 , 44 , 211, 212, 235, 236};
    int Num12[8] = {24, 31, 32 , 39 , 216, 223, 224, 231};
    int Num13[8] = {25, 30, 33 , 38 , 217, 222, 225, 230};
    int Num14[8] = {26, 29, 34 , 37 , 218, 221, 226, 229};
    int Num15[8] = {27, 28, 35 , 36 , 219, 220, 227, 228};
    int Num16[8] = {64, 71, 120, 127, 128, 135, 184, 191};
    int Num17[8] = {65, 70, 121, 126, 129, 134, 185, 190};
    int Num18[8] = {66, 69, 122, 125, 130, 133, 186, 189};
    int Num19[8] = {67, 68, 123, 124, 131, 132, 187, 188};
    int Num20[8] = {72, 79, 112, 119, 136, 143, 176, 183};
    int Num21[8] = {73, 78, 113, 118, 137, 142, 177, 182};
    int Num22[8] = {74, 77, 114, 117, 138, 141, 178, 181};
    int Num23[8] = {75, 76, 115, 116, 139, 140, 179, 180};
    int Num24[8] = {80, 87, 104, 111, 144, 151, 168, 175};
    int Num25[8] = {81, 86, 105, 110, 145, 150, 169, 174};
    int Num26[8] = {82, 85, 106, 109, 146, 149, 170, 173};
    int Num27[8] = {83, 84, 107, 108, 147, 148, 171, 172};
    int Num28[8] = {88, 95, 96 , 103, 152, 159, 160, 167};
    int Num29[8] = {89, 94, 97 , 102, 153, 158, 161, 166};
    int Num30[8] = {90, 93, 98 , 101, 154, 157, 162, 165};
    int Num31[8] = {91, 92, 99 , 100, 155, 156, 163, 164};

    for (int i = 0; i < 8; i++)
    {
        Vcp[Num0[i]]  = x[0];
        Vcp[Num1[i]]  = x[1];
        Vcp[Num2[i]]  = x[2];
        Vcp[Num3[i]]  = x[3];
        Vcp[Num4[i]]  = x[4];
        Vcp[Num5[i]]  = x[5];
        Vcp[Num6[i]]  = x[6];
        Vcp[Num7[i]]  = x[7];
        Vcp[Num8[i]]  = x[8];
        Vcp[Num9[i]]  = x[9];
        Vcp[Num10[i]] = x[10];
        Vcp[Num11[i]] = x[11];
        Vcp[Num12[i]] = x[12];
        Vcp[Num13[i]] = x[13];
        Vcp[Num14[i]] = x[14];
        Vcp[Num15[i]] = x[15];
        Vcp[Num16[i]] = x[16];
        Vcp[Num17[i]] = x[17];
        Vcp[Num18[i]] = x[18];
        Vcp[Num19[i]] = x[19];
        Vcp[Num20[i]] = x[20];
        Vcp[Num21[i]] = x[21];
        Vcp[Num22[i]] = x[22];
        Vcp[Num23[i]] = x[23];
        Vcp[Num24[i]] = x[24];
        Vcp[Num25[i]] = x[25];
        Vcp[Num26[i]] = x[26];
        Vcp[Num27[i]] = x[27];
        Vcp[Num28[i]] = x[28];
        Vcp[Num29[i]] = x[29];
        Vcp[Num30[i]] = x[30];
        Vcp[Num31[i]] = x[31];
    }

    EvalVolumeRatio3D(Vcp, vcratio, 8, 8, 4);

    cVector center(2); center[0] = 0.0; center[1] = 0.0; // parametric coordinates
    double radius = 0.2;                                 // parametric coordinates
    double vcratiohole;

    EvalVolumeHole(center, radius, Vcp, vcratiohole);    // center of the hole ; radius ; Vcp ; vcratiohole

    // Plate volume

    double Area   = 1.0;
    double thk    = 0.02;
    double Vplate = Area*thk;

    // Hole volume

    double pi = atan(1.0)*4.0;
    double rad      = 0.1;
    double AreaHole = pi*rad*rad;
    double Vhole    = AreaHole*thk;

    // Total volume

    double Vol   = Vplate - Vhole;
    double VcTot = vcratio*Vplate - vcratiohole*Vhole;
    double Vcrtot = VcTot/Vol;

    // Cost

    // 0.966920

    double rhom = 8000, rhoc = 2730;
    double Cm   = 3   , Cc   = 50;

    double TotalCost = 0.966920*thk*(Vcrtot*rhoc*Cc + (1 - Vcrtot)*rhom*Cm);

    // buck: lbdn ; cost: TotalCost

    // Stores the value for each objective function

    double w = W_MObj;
    double m = 2.0;

    double CostMin = 464.1216;
    double CostMax = 2639.6916;
    double BuckMin = -2.154576253;
    double BuckMax = -1.177744365;

    out << "'CeramicFraction' 'TotalCost'\n" ;
    out << scientific << Vcrtot << " " << TotalCost << "\n\n";
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
