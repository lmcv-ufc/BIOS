// -------------------------------------------------------------------------
// metaprob.cpp - implementation of the metaproblem class.
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
// Created:      16-Aug-2014    Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <sstream>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "metaprob.h"
#include "fstream"
#include "utl.h"
#include "vec.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Static variables:
//

eMetaProbType cMetaProblem :: InpType;

// -------------------------------------------------------------------------
// Class cMetaProblem:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= CreateMetaProblem =========================

cMetaProblem* cMetaProblem :: CreateMetaProblem(eMetaProbType type)
{
  cMetaProblem *metaprob = 0;

  switch(type)
  {
    case LAMINATE_PLATE_SET01:
      metaprob = new cLamPltSet01;
    break;
    
    case LAMINATE_PLATE_SET02:
      metaprob = new cLamPltSet02;
    break;

    case LAMINATE_PLATE_SET_MIN_WEIGHT_01:
      metaprob = new cLamPltSetMinWeight01;
    break;

    case LAMINATE_PLATE_SET_MIN_WEIGHT_02:
      metaprob = new cLamPltSetMinWeight02;
    break;

    case LAMINATE_PLATE_SET_MIN_COST_01:
      metaprob = new cLamPltSetMinCost01;
    break;
  }

  return(metaprob);
}

// ============================== ReadMetaProblem ==============================

void cMetaProblem :: ReadMetaProblem(istream &in)
{
  // Read the metaproblem label.

  char label[100];
  if (!Utl::ReadString(in, label))
    Utl :: Exit("Error in the input of the problem type label.");

  // Create the appropriate metaproblem.

  if (string(label) == "LamPltSet01")
    InpType = LAMINATE_PLATE_SET01;
  else if (string(label) == "LamPltSet02")
    InpType = LAMINATE_PLATE_SET02;
  else if (string(label) == "LamPltSetMinWeight01")
    InpType = LAMINATE_PLATE_SET_MIN_WEIGHT_01;
  else if (string(label) == "LamPltSetMinWeight02")
    InpType = LAMINATE_PLATE_SET_MIN_WEIGHT_02;
  else if (string(label) == "LamPltSetMinCost01")
    InpType = LAMINATE_PLATE_SET_MIN_COST_01;
  else
    Utl :: Exit("Unknown metaproblem type: "+string(label));
}

// ============================== MetaProblem =============================

cMetaProblem :: cMetaProblem(void)
{
  BestObjFunc = 0;
}

// ============================== ~MetaProblem ============================

cMetaProblem :: ~cMetaProblem(void)
{
  if (BestObjFunc)
    delete [] BestObjFunc;
}

// -------------------------------------------------------------------------
// Class cLamPltSet01:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cLamPltSet01 ================================

cLamPltSet01 :: cLamPltSet01(void) : cMetaProblem( )
{
  Type = LAMINATE_PLATE_SET01;
  ObjType = MAXIMIZATION;
  NumProb = 28;

  BestObjFunc = new double [28];

  BestObjFunc[0]  = -12690.7;
  BestObjFunc[1]  = -10008.8;
  BestObjFunc[2]  = -7967.54;
  BestObjFunc[3]  = -6617.98;
  BestObjFunc[4]  = -3447.57;
  BestObjFunc[5]  = -2550.72;
  BestObjFunc[6]  = -2024.23;
  BestObjFunc[7]  = -1677.01;
  BestObjFunc[8]  = -1517.08;
  BestObjFunc[9]  = -1160.12;
  BestObjFunc[10] = -939.145;
  BestObjFunc[11] = -788.882;
  BestObjFunc[12] = -936.839;
  BestObjFunc[13] = -780.669;
  BestObjFunc[14] = -669.171;
  BestObjFunc[15] = -585.525;
  BestObjFunc[16] = -494.263;
  BestObjFunc[17] = -466.815;
  BestObjFunc[18] = -442.246;
  BestObjFunc[19] = -419.251;
  BestObjFunc[20] = -465.200;
  BestObjFunc[21] = -452.958;
  BestObjFunc[22] = -441.344;
  BestObjFunc[23] = -413.043;
  BestObjFunc[24] = -455.328;
  BestObjFunc[25] = -448.430;
  BestObjFunc[26] = -441.737;
  BestObjFunc[27] = -413.624;
}

// =========================== ~cLamPltSet01 ===============================

cLamPltSet01 :: ~cLamPltSet01(void)
{
  delete [] BestObjFunc;
}
  
// ============================= WriteProblem ==============================

void cLamPltSet01 :: WriteProblem(string &fname, fstream &optfile, int id)
{
  // Open the problem file 

  string probfname = fname + ".lam";
  fstream in(probfname.c_str( ), std::fstream::out | fstream::trunc);

  static const int Col = 4;
  static const int Row = 7;

  static double LoadNY[Row][Col] = {{43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0}};

  static double DimB[Row][Col] = {{0.127,0.127,0.127,0.127},
                                  {0.254,0.254,0.254,0.254},
                                  {0.381,0.381,0.381,0.381},
                                  {0.508,0.508,0.508,0.508},
                                  {1.016,1.016,1.016,1.016},
                                  {1.524,1.524,1.524,1.524},
                                  {2.032,2.032,2.032,2.032}};

  // Write Individual type in .opt file

  optfile << "\n%INDIVIDUAL.TYPE\n";
  optfile << "'IntegerMatrix'\n";

  // Write problem name in .opt file

  optfile << "\n%PROBLEM.TYPE\n";
  optfile << "'LamPltLoadFactor'\n";

  // Write the Header

  in << "\n%HEADER\nThis file was created by BIOSv3\n";

  // Write problem condictions

  int i,j;

  i = id/Col;
  j = id%Col;

  in << "\n%LAMPLT.LOADS\n";
  in << "175.0 " << LoadNY[i][j] << " 0 0 0 0" << endl;

  in << "\n%LAMPLT.DIMENSIONS\n";
  in << "0.508 " << DimB[i][j] << endl;

  // Write material

  in << "\n%MATERIAL\n";
  in << "1\n";
  in << "\n%MATERIAL.ORTHOTROPIC\n";
  in << "1\n";
  in << "1 127.59e9 13.03e9 9e9 0.3 0.3 0.49 6.41e9 7.1e9 6.2e9 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10\n";
  in << "\n%MATERIAL.DENSITY\n";
  in << "1\n";
  in << "1 1605\n";
  
  // Laminate properties

  in << "\n%LAMINATE.ANGLE.RANGE\n";
  in << "0 45 90\n";
  in << "\n%LAMINATE.THICKNESS.VALUES\n";
  in << "1\n";
  in << "1.27e-4\n";
  in << "\n%LAMINATE.TYPE\n";
  in << "'Symmetric_Balanced'\n";
  in << "\n%LAMINATE.MAXIMUM.PLY.NUMBER\n";
  in << "48\n";  
  
  in << "\n%END\n";  

  // Close problem file

  in.close( );
}

// -------------------------------------------------------------------------
// Class cLamPltSet02:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cLamPltSet02 ================================

cLamPltSet02 :: cLamPltSet02(void) : cMetaProblem( )
{
  Type = LAMINATE_PLATE_SET02;
  ObjType = MAXIMIZATION;
  NumProb = 28;

  BestObjFunc = new double [28];

  BestObjFunc[0]  = -13479.8;
  BestObjFunc[1]  = -10491.4;
  BestObjFunc[2]  = -8339.42;
  BestObjFunc[3]  = -6882.55;
  BestObjFunc[4]  = -3564.13;
  BestObjFunc[5]  = -2699.46;
  BestObjFunc[6]  = -2114.29;
  BestObjFunc[7]  = -1737.56;
  BestObjFunc[8]  = -1569.56;
  BestObjFunc[9]  = -1200.25;
  BestObjFunc[10] = -971.631;
  BestObjFunc[11] = -816.170;
  BestObjFunc[12] = -936.839;
  BestObjFunc[13] = -780.699;
  BestObjFunc[14] = -669.171;
  BestObjFunc[15] = -585.525;
  BestObjFunc[16] = -531.549;
  BestObjFunc[17] = -501.845;
  BestObjFunc[18] = -466.904;
  BestObjFunc[19] = -434.391;
  BestObjFunc[20] = -522.398;
  BestObjFunc[21] = -499.792;
  BestObjFunc[22] = -464.401;
  BestObjFunc[23] = -430.641;
  BestObjFunc[24] = -519.543;
  BestObjFunc[25] = -499.766;
  BestObjFunc[26] = -464.814;
  BestObjFunc[27] = -431.275;
}

// =========================== ~cLamPltSet02 ===============================

cLamPltSet02 :: ~cLamPltSet02(void)
{
  delete [] BestObjFunc;
}
  
// ============================= WriteProblem ==============================

void cLamPltSet02 :: WriteProblem(string &fname, fstream &optfile, int id)
{
  // Open the problem file 

  string probfname = fname + ".lam";
  fstream in(probfname.c_str( ), std::fstream::out | fstream::trunc);

  static const int Col = 4;
  static const int Row = 7;

  static double LoadNY[Row][Col] = {{43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0},
                                    {43.75,87.5,131.25,175.0}};

  static double DimB[Row][Col] = {{0.127,0.127,0.127,0.127},
                                  {0.254,0.254,0.254,0.254},
                                  {0.381,0.381,0.381,0.381},
                                  {0.508,0.508,0.508,0.508},
                                  {1.016,1.016,1.016,1.016},
                                  {1.524,1.524,1.524,1.524},
                                  {2.032,2.032,2.032,2.032}};

  // Write Individual type in .opt file

  optfile << "\n%INDIVIDUAL.TYPE\n";
  optfile << "'IntegerMatrix'\n";

  // Write problem name in .opt file

  optfile << "\n%PROBLEM.TYPE\n";
  optfile << "'LamPltLoadFactor'\n";

  // Write the Header

  in << "\n%HEADER\nThis file was created by BIOSv3\n";
  
  // Write problem condictions

  int i,j;

  i = id/Col;
  j = id%Col;

  in << "\n%LAMPLT.LOADS\n";
  in << "175.0 " << LoadNY[i][j] << " 0 0 0 0" << endl;

  in << "\n%LAMPLT.DIMENSIONS\n";
  in << "0.508 " << DimB[i][j] << endl;

  // Write material

  in << "\n%MATERIAL\n";
  in << "1\n";
  in << "\n%MATERIAL.ORTHOTROPIC\n";
  in << "1\n";
  in << "1 127.59e9 13.03e9 9e9 0.3 0.3 0.49 6.41e9 7.1e9 6.2e9 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10 1e10\n";
  in << "\n%MATERIAL.DENSITY\n";
  in << "1\n";
  in << "1 1605\n";
  
  // Laminate properties

  in << "\n%LAMINATE.ANGLE.RANGE\n";
  in << "0 15 90\n";
  in << "\n%LAMINATE.THICKNESS.VALUES\n";
  in << "1\n";
  in << "1.27e-4\n";
  in << "\n%LAMINATE.TYPE\n";
  in << "'Symmetric_Balanced'\n";
  in << "\n%LAMINATE.MAXIMUM.PLY.NUMBER\n";
  in << "48\n";  
  
  in << "\n%END\n";  
  
  // Close problem file

  in.close( );
}

// -------------------------------------------------------------------------
// Class cLamPltSetMinWeight01:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cLamPltSetMinWeight01 ========================

cLamPltSetMinWeight01 :: cLamPltSetMinWeight01(void) : cMetaProblem( )
{
  Type = LAMINATE_PLATE_SET_MIN_WEIGHT_01;
  ObjType = MINIMIZATION;
  NumProb = 9;

  BestObjFunc = new double [9];

  BestObjFunc[0]  = 76.2696;
  BestObjFunc[1]  = 82.6254;
  BestObjFunc[2]  = 95.337;
  BestObjFunc[3]  = 63.558;
  BestObjFunc[4]  = 82.6254;
  BestObjFunc[5]  = 101.693;
  BestObjFunc[6]  = 31.779;
  BestObjFunc[7]  = 38.1348;
  BestObjFunc[8]  = 44.4906;
}

// =========================== ~cLamPltSetMinWeight01 ===============================

cLamPltSetMinWeight01 :: ~cLamPltSetMinWeight01(void)
{
  delete [] BestObjFunc;
}
  
// ============================= WriteProblem ==============================

void cLamPltSetMinWeight01 :: WriteProblem(string &fname, fstream &optfile, int id)
{
  // Open the problem file 

  string probfname = fname + ".lam";
  fstream in(probfname.c_str( ), std::fstream::out | fstream::trunc);

  static double LoadNX[9] = {2e6,2e6,2e6,
                             2e6,2e6,2e6,
                            -2e6,-2e6,-2e6};

  static double LoadNY[9] = {2e6,2e6,2e6,
                            -2e6,-2e6,-2e6,
                            -2e6,-2e6,-2e6};

  static double LoadNXY[9] = {0.0,5e5,10e5,
                             0.0,5e5,10e5,
                             0.0,5e5,10e5};

  // Write Individual type in .opt file

  optfile << "\n%INDIVIDUAL.TYPE\n";
  optfile << "'IntegerMatrix'\n";

  // Write problem name in .opt file

  optfile << "\n%PROBLEM.TYPE\n";
  optfile << "'LamPltMinLopez2009'\n";

  // Write the Header

  in << "\n%HEADER\nThis file was created by BIOSv3\n";
  
  // Write problem condictions
  
  in << "\n%LAMPLT.LOADS\n";
  in << LoadNX[id]  << " " << LoadNY[id] << " " << LoadNXY[id] << " 0 0 0\n";

  in << "\n%LAMPLT.DIMENSIONS\n";
  in << "1.0 1.0\n" << endl;

  // Write material

  in << "\n%MATERIAL\n";
  in << "1\n";
  in << "\n%MATERIAL.ORTHOTROPIC\n";
  in << "1\n";
  in << "1 116.6e9 7.673e9 9e9 0.27 0.27 0.49 4.173e9 7.1e9 6.2e9 2062e6 1701e6 70e6 240e6 230e9 230e9 105e6 105e6 105e6\n";
  in << "\n%MATERIAL.DENSITY\n";
  in << "1\n";
  in << "1 1605\n";
  
  // Laminate properties

  in << "\n%ORTHOTROPIC.FAILURE.CRITERION\n";
  in << "'TsaiWu'\n";
  in << "\n%LAMINATE.ANGLE.RANGE\n";
  in << "0 45 90\n";
  in << "\n%LAMINATE.THICKNESS.VALUES\n";
  in << "2\n";
  in << "1.0e-4\n";
  in << "0.0\n";
  in << "\n%LAMINATE.TYPE\n";
  in << "'Symmetric_Balanced'\n";
  in << "\n%LAMINATE.MAXIMUM.PLY.NUMBER\n";
  in << "200\n";  
  
  in << "\n%END\n";  
  
  // Close problem file

  in.close( );
}

// -------------------------------------------------------------------------
// Class cLamPltSetMinWeight02:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cLamPltSetMinWeight02 ========================

cLamPltSetMinWeight02 :: cLamPltSetMinWeight02(void) : cMetaProblem( )
{
  Type = LAMINATE_PLATE_SET_MIN_WEIGHT_02;
  ObjType = MINIMIZATION;
  NumProb = 9;

  BestObjFunc = new double [9];

  BestObjFunc[0]  = 76.2696;
  BestObjFunc[1]  = 82.6254;
  BestObjFunc[2]  = 95.337;
  BestObjFunc[3]  = 63.558;
  BestObjFunc[4]  = 82.6254;
  BestObjFunc[5]  = 101.693;
  BestObjFunc[6]  = 31.779;
  BestObjFunc[7]  = 38.1348;
  BestObjFunc[8]  = 44.4906;
}

// =========================== ~cLamPltSetMinWeight02 ===============================

cLamPltSetMinWeight02 :: ~cLamPltSetMinWeight02(void)
{
  delete [] BestObjFunc;
}
  
// ============================= WriteProblem ==============================

void cLamPltSetMinWeight02 :: WriteProblem(string &fname, fstream &optfile, int id)
{
  // Open the problem file 

  string probfname = fname + ".lam";
  fstream in(probfname.c_str( ), std::fstream::out | fstream::trunc);

  static double LoadNX[9] = {2e6,2e6,2e6,
                             2e6,2e6,2e6,
                            -2e6,-2e6,-2e6};

  static double LoadNY[9] = {2e6,2e6,2e6,
                            -2e6,-2e6,-2e6,
                            -2e6,-2e6,-2e6};

  static double LoadNXY[9] = {0.0,5e5,10e5,
                             0.0,5e5,10e5,
                             0.0,5e5,10e5};

  // Write Individual type in .opt file

  optfile << "\n%INDIVIDUAL.TYPE\n";
  optfile << "'IntegerMatrix'\n";

  // Write problem name in .opt file

  optfile << "\n%PROBLEM.TYPE\n";
  optfile << "'LamPltMinLopez2009'\n";

  // Write the Header

  in << "\n%HEADER\nThis file was created by BIOSv3\n";

  // Write problem condictions
  
  in << "\n%LAMPLT.LOADS\n";
  in << LoadNX[id]  << " " << LoadNY[id] << " " << LoadNXY[id] << " 0 0 0\n";

  in << "\n%LAMPLT.DIMENSIONS\n";
  in << "1.0 1.0\n" << endl;

  // Write material

  in << "\n%MATERIAL\n";
  in << "1\n";
  in << "\n%MATERIAL.ORTHOTROPIC\n";
  in << "1\n";
  in << "1 116.6e9 7.673e9 9e9 0.27 0.27 0.49 4.173e9 7.1e9 6.2e9 2062e6 1701e6 70e6 240e6 230e9 230e9 105e6 105e6 105e6\n";
  in << "\n%MATERIAL.DENSITY\n";
  in << "1\n";
  in << "1 1605\n";
  
  // Laminate properties

  in << "\n%ORTHOTROPIC.FAILURE.CRITERION\n";
  in << "'TsaiWu'\n";
  in << "\n%LAMINATE.ANGLE.RANGE\n";
  in << "0 15 90\n";
  in << "\n%LAMINATE.THICKNESS.VALUES\n";
  in << "2\n";
  in << "1.0e-4\n";
  in << "0.0\n";
  in << "\n%LAMINATE.TYPE\n";
  in << "'Symmetric_Balanced'\n";
  in << "\n%LAMINATE.MAXIMUM.PLY.NUMBER\n";
  in << "200\n";  
  
  in << "\n%END\n";  

  // Close problem file

  in.close( );
}

// -------------------------------------------------------------------------
// Class cLamPltSetMinCost01:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== cLamPltSetMinCost01 ================================

cLamPltSetMinCost01 :: cLamPltSetMinCost01(void) : cMetaProblem( )
{
  Type = LAMINATE_PLATE_SET_MIN_COST_01;
  ObjType = MINIMIZATION;
  NumProb = 3;

  BestObjFunc = new double [3];

  BestObjFunc[0]  = -13479.8;
  BestObjFunc[1]  = -10491.4;
  BestObjFunc[2]  = -8339.42;
}

// =========================== ~cLamPltSet02 ===============================

cLamPltSetMinCost01 :: ~cLamPltSetMinCost01(void)
{
  delete [] BestObjFunc;
}
  
// ============================= WriteProblem ==============================

void cLamPltSetMinCost01 :: WriteProblem(string &fname, fstream &optfile, int id)
{
  // Open the problem file 

  string probfname = fname + ".lam";
  fstream in(probfname.c_str( ), std::fstream::out | fstream::trunc);

  static double LoadNX[3] = {2e6,-2e6,-2e6};
  static double LoadNY[3] = {2e6,-2e6,2e6};

  // Write Individual type in .opt file

  optfile << "\n%INDIVIDUAL.TYPE\n";
  optfile << "'IntegerMatrix'\n";

  // Write problem name in .opt file

  optfile << "\n%PROBLEM.TYPE\n";
  optfile << "'LamPltMinCost'\n";

  // Write the Header

  in << "\n%HEADER\nThis file was created by BIOSv3\n";

  // Write problem condictions

  in << "\n%LAMPLT.LOADS\n";
  in << LoadNX[id] << " " << LoadNY[id] << " 0 0 0 0" << endl;

  in << "\n%LAMPLT.DIMENSIONS\n";
  in << "1.0 1.0" << endl;

  // Write material

  in << "\n%MATERIAL\n";
  in << "2\n";
  in << "\n%MATERIAL.ORTHOTROPIC\n";
  in << "2\n";
  in << "1 116.6e9 7.673e9 9e9 0.27 0.27 0.49 4.173e9 7.1e9 6.2e9 2062e6 1701e6 70e6 240e6 230e9 230e9 105e6 105e6 105e6\n";
  in << "2 37.6e9 9.584e9 9.584e9 0.26 0.26 0.49 4.081e8 4.081e9 4.081e9 1134e6 1031e6 54e6 150e6 1e10 1e10 75e6 75e6 75e6\n";
  in << "\n%MATERIAL.DENSITY\n";
  in << "2\n";
  in << "1 1605\n";
  in << "2 1903\n";

  in << "\n%ORTHOTROPIC.FAILURE.CRITERION\n";
  in << "'TsaiWu'\n";
  
  // Laminate properties

  in << "\n%LAMINATE.ANGLE.RANGE\n";
  in << "0 45 90\n";
  in << "\n%LAMINATE.THICKNESS.VALUES\n";
  in << "1\n";
  in << "1.0e-4\n";
  in << "\n%LAMINATE.TYPE\n";
  in << "'Symmetric_Balanced'\n";
  in << "\n%LAMINATE.MAXIMUM.PLY.NUMBER\n";
  in << "100\n";  
  
  in << "\n%END\n";  

  // Close problem file

  in.close( );
}

// ======================================================= End of file =====
