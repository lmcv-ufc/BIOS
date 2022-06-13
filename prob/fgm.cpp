// -------------------------------------------------------------------------
// fgm.cpp - Implementation of the Functionally Graded Material problem class.
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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#ifdef _OMP_
#include "omp.h"
#endif

#include "fgm.h"
#include "input.h"
#include "gbldef.h"
#include "gblvar.h"
#include "utl.h"
#include "mat.h"
#include "vec.h"
#include "material.h"
#include <queue>

using namespace std;

// -------------------------------------------------------------------------
// Static variables:
//

eFGMModel   cFGM :: FGMModel     = VOIGT;
eFGMVolDist cFGM :: FGMVolDist   = POLYNOMIAL;
int         cFGM :: NumFGMMat    = 2;
double      cFGM :: dExp         = 1.00;
double      cFGM :: MinExp       = 0.00;
double      cFGM :: MaxExp       = 10.00;
cVector     cFGM :: FGIDMat;
cVector     cFGM :: FGMMat;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ Init =======================================

void cFGM :: Init( )
{    // Create an array with the size of the list of each variable.

    int index = 0;

    if(ThkVar){
        index = 1;
    }

    NumVar = NumCP + index;

    cout << "NumCp: " << NumCP << " numero de variaveis: " << NumVar << endl;

    Low = new double[NumVar];
    Upp = new double[NumVar];

    if(ThkVar){
        Low[0] = MinThk;
        Upp[0] = MaxThk;

        cout << "MinThk: " << Low[0] << " UpperThk: " << Upp[0] << endl;
    }

    if(FGMVolDist == POWER_LAW || FGMVolDist == SIGMOYD || FGMVolDist == POLYNOMIAL){
        if (NumCP > 0)
        {
        Low[index] = MinExp;
        Upp[index] = MaxExp;
        cout << "MinExp: " << Low[index] << " UpperExp: " << Upp[index] << endl;
        }
    }

    if(FGMVolDist == PCHIP || FGMVolDist == BSPLINE){
        for (int i = 0; i < NumVar; i++){
            Low[i + index] = 0.0;
            Upp[i + index] = 1.0;
        }
     }
}

// ============================== LoadReadFunc =============================

void cFGM :: LoadReadFunc(cInpMap &im)
{
  // Register read functions.
  im.Insert("FGM.MATERIALS",makeReadObj(cFGM, ReadFGMaterials));
  im.Insert("FGM.VOLUME.DISTRIBUTION",makeReadObj(cFGM, ReadFGMVolDist));
  im.Insert("FGM.MODEL",makeReadObj(cFGM, ReadFGMModel));
  im.Insert("FGM.EXPONENT.RANGE",makeReadObj(cFGM, ReadFGMExpRange));
  im.Insert("FGM.NUMBER.OF.CONTROL.POINTS",makeReadObj(cFGM, ReadNumberOfControlPoints));
  im.Insert("FGM.THICKNESS.RANGE",makeReadObj(cFGM, ReadFGMThkRange));

  // Register material read functions.
  im.Insert("MATERIAL",makeRead(cMaterial :: ReadNumMat));
  im.Insert("MATERIAL.DENSITY",makeRead(cMaterial :: ReadDensity));
  im.Insert("MATERIAL.COST",makeRead(cMaterial :: ReadCost));
  im.Insert("MATERIAL.ISOTROPIC",makeRead(cMaterial :: ReadIso));
  im.Insert("MATERIAL.ORTHOTROPIC",makeRead(cMaterial :: ReadOrtho));
  im.Insert("ISOTROPIC.FAILURE.CRITERION",makeRead(cMatIso :: ReadFailType));
  im.Insert("ORTHOTROPIC.FAILURE.CRITERION",makeRead(cMatOrtho :: ReadFailType));
}

// =========================== ReadFGMaterials =============================

void cFGM :: ReadFGMaterials(istream &in)
{
  if (!(in >> NumFGMMat) || NumFGMMat == 0 || NumFGMMat > 2)
  {
    cout << "Error in the input of the number of fg materials." << endl;
    exit(0);
  }

  FGIDMat.Resize(NumFGMMat);

  for (int i = 0; i < NumFGMMat; i++)
  {
    if (!(in >> FGIDMat[i]))
    {
      cout << "Error in the input of the possible functionally graded materials." << endl;
      exit(0);
    }
  }

  //FGIDMat.Print();

  // Fill vector with material info

  FGMMat.Resize(NumFGMMat*3);

  for (int i = 0; i < NumFGMMat; i++)
  {
    cMaterial *mat = cMaterial::GetMaterial(FGIDMat[i]);
    double *paramtemp = new double[mat->NumParam( )];
    double dens       = mat -> GetDensity();

    mat->GetParam(paramtemp);

    eMatType matid  = mat->GetType( );

   if (matid == MAT_ISOTROPIC)
   {
       FGMMat[i*3]   = paramtemp[0];    // E1
       FGMMat[i*3+1] = paramtemp[1];    // nu
       FGMMat[i*3+2] = dens;            // density
   }
   else
   {
       FGMMat[i*3]   = paramtemp[0];    // E1
       FGMMat[i*3+1] = paramtemp[3];    // nu
       FGMMat[i*3+2] = dens;            // density
   }

   delete []paramtemp;
  }

  //cout << "PROPRIEDADES" << endl;
  //FGMMat.Print();
}

// ============================== ReadFGMModel ==============================

void cFGM :: ReadFGMModel(istream &in)
{
  char label[100];

  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of FGM model." << endl;
    exit(0);
  }

  if (string(label) == "Voigt" || string(label) == "voigt")
    FGMModel = VOIGT;

  else if (string(label) == "Reuss" || string(label) == "reuss")
    FGMModel = REUSS;

  else if (string(label) == "Mori-tanaka" || string(label) == "Mori-Tanaka"
           || string(label) == "MoriTanaka" || string(label) == "Moritanaka")
    FGMModel = MORI_TANAKA;

  else
  {
    cout << "Unknown fgm model type: " << label << endl;
    exit(0);
  }
}

// ============================== ReadFGMVolDist ==============================

void cFGM :: ReadFGMVolDist(istream &in)
{
  char label[100];

  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of FGM volume fraction." << endl;
    exit(0);
  }

  if (string(label) == "Power-Law" || string(label) == "Power-law"
          || string(label) == "Powerlaw" || string(label) == "PowerLaw"){
    FGMVolDist = POWER_LAW;
    NumCP = 1;
  }

  else if (string(label) == "Sigmoyd" || string(label) == "sigmoyd"){
    FGMVolDist = SIGMOYD;
    NumCP = 1;
  }
  else if (string(label) == "Polynomial" || string(label) == "polynomial"){
    FGMVolDist = POLYNOMIAL;
    NumCP = 1;
  }
  else if (string(label) == "PiecewiseCubicInterpolation" || string(label) == "HermiteCubicInterpolation" || string(label) == "PCHIP"){
    FGMVolDist = PCHIP;
  }
  else if (string(label) == "BSpline" || string(label) == "B-spline" || string(label) == "B-Spline"){
    FGMVolDist = BSPLINE;
  }

  else
  {
    cout << "Unknown fgm volume fraction distribution type: " << label << endl;
    exit(0);
  }
}

// ============================== ReadFGMExpRange ==============================

void cFGM :: ReadFGMExpRange(istream &in)
{
  if (!(in >> MinExp) || !(in >> MaxExp))
  {
    cout << "Error in the input of the exponent range." << endl;
    exit(0);
  }
}

// ============================== ReadFGMExpRange ==============================

void cFGM :: ReadNumberOfControlPoints(istream &in)
{
  if (!(in >> NumCP))
  {
    cout << "Error in the input of the number of control points" << endl;
    exit(0);
  }
}

// ============================== ReadFGMThkRange ==============================

void cFGM :: ReadFGMThkRange(istream &in)
{
  ThkVar = true;
  if (!(in >> MinThk) || !(in >> MaxThk))
  {
    cout << "Error in the input of the thickness range." << endl;
    exit(0);
  }
}

// ============================= cFGM ================================

cFGM :: cFGM(void)
{
    ThkVar = false;
}

// ============================ ~cFGM ================================

cFGM :: ~cFGM(void)
{

}

// ============================== PrintVar =================================

void cFGM :: PrintVar(int **algvar)
{
  // Decode the variables.

}

// ============================== PrintResult ======================

void cFGM :: PrintResult(int **algvar, ostream &res)
{
  /*cVector countAng(Ang.Dim());
  ostringstream currline, lamline;
  queue<pair<string,string> > OutPutQueue;

  countAng.Zero();

  // Abort if Feedback is active

  if (!Feedback)
    return;

  // Decode the variables.

  cMatrix laminate;
  Decode(algvar, laminate);

  // Write the objective function.

  res << fixed;
  res << "-------------------------------------------------------------\n";
  res << "Optimization results\n";
  res << "-------------------------------------------------------------\n";
  res << "Lamina   Thkn    Angle   Mat\n";
  res << "-------------------------------------------------------------\n";

  int lam = 1;

  for (int j = 0; j < laminate.NCol( ); j++)
  {

    lamline << setw(6) << setprecision(1) << left << lam++ << "  ";
    currline << setw(6) << setprecision(3) << right << laminate[0][j] << "  ";
    currline << setw(5) << setprecision(3) << laminate[1][j] << "  ";

    if (LamType == BALANCED || LamType == SYMMETRIC_BALANCED)
      currline << setw(4) << setprecision(0) << laminate[2][j] << " (balanced)";

    else
      currline << setw(4) << setprecision(0) << laminate[2][j];

    // Add the output to queue.

    OutPutQueue.push(make_pair<string,string>(lamline.str(),currline.str()));

    // Clear stream objects.

    lamline.str("");
    currline.str("");

    if (LamType == BALANCED || LamType == SYMMETRIC_BALANCED)
    {
      j++;
      lam++;
    }

    if (LamType == SYMMETRIC || LamType == SYMMETRIC_BALANCED)
      if (j+1 == laminate.NCol( )/2)
        break;
  }

  // Write the laminate on the screen.

  while (!OutPutQueue.empty( ))
  {
    // Get current string.

    string currlam  = OutPutQueue.front( ).first;
    string currtext = OutPutQueue.front( ).second;

    // Track how many times this string appear.

    int mult = 1;
    OutPutQueue.pop( );

    while (!OutPutQueue.empty( ))
    {
      if (currtext == OutPutQueue.front( ).second)
        mult++;

      else
        break;

      OutPutQueue.pop( );
    }

    if (mult > 1)
      res << currlam << currtext << " [" << mult << "]\n";

    else
      res << currlam << currtext << endl;
  }

  res << "-------------------------------------------------------------\n";
  res << "Laminate thickness   = " << setprecision(3) << GetThickness(laminate) << endl;

  // Print the number and percentage of each angle

  res << "\n=============== Angle number and percentages ======================= " << endl;

  for (int i=0; i < laminate.NCol(); i++)
    for (int j=0; j < Ang.Dim(); j++)
    {
      if (laminate[1][i] == Ang[j] || laminate[1][i] == -Ang[j]) countAng[j] += 1;
      if (i == (laminate.NCol()-1)) res << "Ang " << setw(2) << Ang[j] << "deg.: " << setw(2) << countAng[j] << " ..... " << (100*countAng[j]/laminate.NCol()) << "\%" << endl;
    }

  res << "======================================\n\n";*/
}

// ============================== WriteVar =================================

void cFGM :: WriteVar(int **algvar)
{
  // Decode the variables.

/*  cMatrix layup;
  Decode(algvar, layup);*/

  PrintResult(algvar, cout);

  // Print the laminate layup.

  /*for (int i = 0; i < layup.NRow( ); i++)
  {
    for (int j = 0; j < layup.NCol( ); j++)
      out << layup[i][j] << " ";

    out << endl;
  }*/
}

// =============================== GetBounds ===============================

void cFGM :: GetBounds(int i, int *low, int *upp)
{
  *low = 0;
  *upp = ListDim[i]-1;
}

// ============================ GetDblBounds ================================

void cFGM :: GetDblBounds(double *low, double *upp)
{
  for (int i = 0; i < NumVar; i++)
  {
    low[i] = Low[i];
    upp[i] = Upp[i];
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ GaussPts1D =================================

void cFGM :: GaussPts1D(int npg, cVector &r, cVector &w)
{
    r.Resize(npg);
    w.Resize(npg);

    if (npg == 1)
    {
        r[0] = 0.000000000000000;
        w[0] = 2.000000000000000;
    }
    else if (npg == 2)
    {
        r[0] = -0.577350269189626;
        r[1] = 0.577350269189626;
        w[0] = 1.000000000000000;
        w[1] = 1.000000000000000;
    }
    else if (npg == 3)
    {
        r[0] = -0.774596669241483;
        r[1] = 0.000000000000000;
        r[2] = 0.774596669241483;
        w[0] = 0.555555555555556;
        w[1] = 0.888888888888889;
        w[2] = 0.555555555555556;
    }
    else if (npg == 4)
    {
        r[0] = -0.861136311594053;
        r[1] = -0.339981043584856;
        r[2] = 0.339981043584856;
        r[3] = 0.861136311594053;
        w[0] = 0.347854845137454;
        w[1] = 0.652145154862546;
        w[2] = 0.652145154862546;
        w[3] = 0.347854845137454;
    }
    else if (npg == 5)
    {
        r[0] = -0.906179845938664;
        r[1] = -0.538469310105683;
        r[2] = 0.000000000000000;
        r[3] = 0.538469310105683;
        r[4] = 0.906179845938664;
        w[0] = 0.23692688505618;
        w[1] = 0.478628670499367;
        w[2] = 0.568888888888889;
        w[3] = 0.478628670499367;
        w[4] = 0.236926885056189;
    }
    else if (npg == 6)
    {
        r[0] = -0.932469514203152;
        r[1] = -0.661209386466264;
        r[2] = -0.238619186083197;
        r[3] = 0.238619186083197;
        r[4] = 0.661209386466264;
        r[5] = 0.932469514203152;
        w[0] = 0.171324492379170;
        w[1] = 0.36076157304813;
        w[2] = 0.467913934572691;
        w[3] = 0.467913934572691;
        w[4] = 0.360761573048139;
        w[5] = 0.171324492379170;
    }
    else if (npg == 7)
    {
        r[0] = -0.949107912342758;
        r[1] = -0.741531185599394;
        r[2] = -0.405845151377397;
        r[3] = 0.000000000000000;
        r[4] = 0.405845151377397;
        r[5] = 0.741531185599394;
        r[6] = 0.949107912342758;
        w[0] = 0.129484966168870;
        w[1] = 0.279705391489277;
        w[2] = 0.381830050505119;
        w[3] = 0.417959183673469;
        w[4] = 0.381830050505119;
        w[5] = 0.279705391489277;
        w[6] = 0.129484966168870;
    }
    else if (npg == 8)
    {
        r[0] = -0.960289856497536;
        r[1] = -0.796666477413627;
        r[2] = -0.525532409916329;
        r[3] = -0.183434642495650;
        r[4] = 0.183434642495650;
        r[5] = 0.525532409916329;
        r[6] = 0.796666477413627;
        r[7] = 0.960289856497536;
        w[0] = 0.101228536290376;
        w[1] = 0.222381034453375;
        w[2] = 0.313706645877887;
        w[3] = 0.362683783378362;
        w[4] = 0.362683783378362;
        w[5] = 0.313706645877887;
        w[6] = 0.222381034453375;
        w[7] = 0.101228536290376;
    }
    else if (npg == 9)
    {
        r[0] = -0.968160239507626;
        r[1] = -0.836031107326636;
        r[2] = -0.613371432700590;
        r[3] = -0.324253423403809;
        r[4] = 0.000000000000000;
        r[5] = 0.324253423403809;
        r[6] = 0.613371432700590;
        r[7] = 0.836031107326636;
        r[8] = 0.968160239507626;
        w[0] = 0.081274388361574;
        w[1] = 0.180648160694857;
        w[2] = 0.260610696402935;
        w[3] = 0.312347077040003;
        w[4] = 0.330239355001260;
        w[5] = 0.312347077040003;
        w[6] = 0.260610696402935;
        w[7] = 0.180648160694857;
        w[8] = 0.081274388361574;
    }
    else if (npg == 10)
    {
        r[0] = -0.973906528517172;
        r[1] = -0.865063366688985;
        r[2] = -0.679409568299024;
        r[3] = -0.433395394129247;
        r[4] = -0.148874338981631;
        r[5] = 0.148874338981631;
        r[6] = 0.433395394129247;
        r[7] = 0.679409568299024;
        r[8] = 0.865063366688985;
        r[9] = 0.973906528517172;
        w[0] = 0.066671344308688;
        w[1] = 0.149451349150581;
        w[2] = 0.219086362515982;
        w[3] = 0.269266719309996;
        w[4] = 0.295524224714753;
        w[5] = 0.295524224714753;
        w[6] = 0.269266719309996;
        w[7] = 0.219086362515982;
        w[8] = 0.149451349150581;
        w[9] = 0.066671344308688;
    }
    else if (npg == 11)
    {
        r[0] = -0.269543155952345;
        r[1] =-0.519096129206812;
        r[2] = -0.730152005574049;
        r[3] = -0.887062599768095;
        r[4] = -0.978228658146057;
        r[5] = 0.000000000000000;
        r[6] = 0.269543155952345;
        r[7] = 0.519096129206812;
        r[8] = 0.730152005574049;
        r[9] = 0.887062599768095;
        r[10] = 0.978228658146057;
        w[0] = 0.262804544510247;
        w[1] = 0.233193764591990;
        w[2] = 0.186290210927734;
        w[3] = 0.125580369464905;
        w[4] = 0.055668567116175;
        w[5] = 0.272925086777901;
        w[6] = 0.262804544510247;
        w[7] = 0.233193764591990;
        w[8] = 0.186290210927734;
        w[9] = 0.125580369464905;
        w[10] = 0.055668567116174;
    }
    else if (npg == 12)
    {
        r[0] = -0.125233408511469;
        r[1] = -0.367831498998180;
        r[2] = -0.587317954286617;
        r[4] =  -0.769902674194305;
        r[4] = -0.904117256370475;
        r[5] = -0.981560634246719;
        r[6] = 0.125233408511469;
        r[7] = 0.367831498998180;
        r[8] = 0.587317954286617;
        r[9] =  0.769902674194305;
        r[10] = 0.904117256370475;
        r[11] = 0.981560634246719;
        w[0] = 0.249147045813403;
        w[1] = 0.233492536538355;
        w[2] = 0.203167426723066;
        w[3] = 0.160078328543346;
        w[4] = 0.106939325995318;
        w[5] = 0.047175336386512;
        w[6] = 0.249147045813403;
        w[7] = 0.233492536538355;
        w[8] = 0.203167426723066;
        w[9] = 0.160078328543346;
        w[10] = 0.106939325995318;
        w[11] = 0.047175336386512;
    }
    else if (npg == 25){
        r[0]  = -0.122864692610710;
        r[1]  = -0.243866883720988;
        r[2]  = -0.361172305809388;
        r[3]  = -0.473002731445715;
        r[4]  = -0.577662930241223;
        r[5]  = -0.673566368473468;
        r[6]  = -0.759259263037357;
        r[7]  = -0.833442628760834;
        r[8]  = -0.894991997878275;
        r[9]  = -0.942974571228974;
        r[10] = -0.976663921459517;
        r[11] = -0.995556969790498;
        r[12] =  0.000000000000000;
        r[13] =  0.122864692610710;
        r[14] =  0.243866883720988;
        r[15] =  0.361172305809388;
        r[16] =  0.473002731445715;
        r[17] =  0.577662930241223;
        r[18] =  0.673566368473468;
        r[19] =  0.759259263037357;
        r[20] =  0.833442628760834;
        r[21] =  0.894991997878275;
        r[22] =  0.942974571228974;
        r[23] =  0.976663921459517;
        r[24] =  0.995556969790498;
        w[0]  =  0.122242442990310;
        w[1]  =  0.119455763535785;
        w[2]  =  0.114858259145712;
        w[3]  =  0.108519624474264;
        w[4]  =  0.100535949067051;
        w[5]  =  0.091028261982964;
        w[6]  =  0.080140700335001;
        w[7]  =  0.068038333812357;
        w[8]  =  0.054904695975835;
        w[9]  =  0.040939156701306;
        w[10] =  0.026354986615032;
        w[11] =  0.011393798501026;
        w[12] =  0.123176053726715;
        w[13] =  0.122242442990310;
        w[14] =  0.119455763535785;
        w[15] =  0.114858259145712;
        w[16] =  0.108519624474264;
        w[17] =  0.100535949067051;
        w[18] =  0.091028261982964;
        w[19] =  0.080140700335001;
        w[20] =  0.068038333812357;
        w[21] =  0.054904695975835;
        w[22] =  0.040939156701306;
        w[23] =  0.026354986615032;
        w[24] =  0.011393798501026;
    }
    /*elseif (npg == 13)
        r = [-0.230458315955135,-0.448492751036447,-0.642349339440340, ...
             -0.801578090733310,-0.917598399222978,-0.984183054718588, ...
              0.000000000000000, 0.230458315955135, 0.448492751036447, ...
              0.642349339440340, 0.801578090733310, 0.917598399222978, ...
              0.984183054718588];
        w = [ 0.226283180262897, 0.207816047536888, 0.178145980761946, ...
              0.138873510219787, 0.092121499837728, 0.040484004765316, ...
              0.232551553230874, 0.226283180262897, 0.207816047536888, ...
              0.178145980761946, 0.138873510219787, 0.092121499837728, ...
              0.040484004765316];
      elseif (npg == 14)
        r = [-0.108054948707344,-0.319112368927890,-0.515248636358154, ...
             -0.687292904811685,-0.827201315069765,-0.928434883663573, ...
             -0.986283808696812, 0.108054948707344, 0.319112368927890, ...
              0.515248636358154, 0.687292904811685, 0.827201315069765, ...
              0.928434883663573, 0.986283808696812];
        w = [ 0.215263853463158, 0.205198463721296, 0.185538397477938, ...
              0.157203167158193, 0.121518570687903, 0.080158087159760, ...
              0.035119460331752, 0.215263853463158, 0.205198463721296, ...
              0.185538397477938, 0.157203167158193, 0.121518570687903, ...
              0.080158087159760, 0.035119460331752];
      elseif (npg == 15)
        r = [-0.201194093997434,-0.394151347077563,-0.570972172608539, ...
             -0.724417731360170,-0.848206583410427,-0.937273392400706, ...
             -0.987992518020485, 0.000000000000000, 0.201194093997434, ...
              0.394151347077563, 0.570972172608539, 0.724417731360170, ...
              0.848206583410427, 0.937273392400706, 0.987992518020485];
        w = [ 0.198431485327112, 0.186161000015562, 0.166269205816994, ...
              0.139570677926154, 0.107159220467172, 0.070366047488108, ...
              0.030753241996117, 0.202578241925561, 0.198431485327112, ...
              0.186161000015562, 0.166269205816994, 0.139570677926154, ...
              0.107159220467172, 0.070366047488108, 0.030753241996117];
      elseif (npg == 16)
        r = [-0.095012509837637,-0.281603550779259,-0.458016777657227, ...
             -0.617876244402644,-0.755404408355003,-0.865631202387832, ...
             -0.944575023073233,-0.989400934991650, 0.095012509837637, ...
              0.281603550779259, 0.458016777657227, 0.617876244402644, ...
              0.755404408355003, 0.865631202387832, 0.944575023073233, ...
              0.98940093499165];
        w = [ 0.189450610455068, 0.182603415044924, 0.169156519395002, ...
              0.149595988816577, 0.124628971255534, 0.095158511682493, ...
              0.062253523938648, 0.027152459411754, 0.189450610455068, ...
              0.182603415044924, 0.169156519395002, 0.149595988816577, ...
              0.124628971255534, 0.095158511682493, 0.062253523938648, ...
              0.027152459411754];
      elseif (npg == 17)
        r = [-0.178484181495848,-0.351231763453876,-0.512690537086477, ...
             -0.657671159216691,-0.657671159216691,-0.880239153726986, ...
             -0.950675521768768,-0.990575475314417, 0.000000000000000, ...
              0.178484181495848, 0.351231763453876, 0.512690537086477, ...
              0.657671159216691, 0.781514003896801, 0.880239153726986, ...
              0.950675521768768, 0.990575475314417];
        w = [ 0.176562705366993, 0.168004102156450, 0.154045761076810, ...
              0.135136368468525, 0.111883847193404, 0.085036148317179, ...
              0.055459529373987, 0.024148302868548, 0.179446470356206, ...
              0.176562705366993, 0.168004102156450, 0.154045761076810, ...
              0.135136368468525, 0.111883847193404, 0.085036148317179, ...
              0.055459529373987, 0.024148302868548];
      elseif (npg == 18)
        r = [-0.084775013041735,-0.251886225691505,-0.411751161462843, ...
             -0.559770831073947,-0.691687043060353,-0.803704958972523, ...
             -0.892602466497556,-0.955823949571398,-0.991565168420931, ...
              0.084775013041735, 0.251886225691505, 0.411751161462843, ...
              0.559770831073947, 0.691687043060353, 0.803704958972523, ...
              0.892602466497556, 0.955823949571398, 0.991565168420931];
        w = [ 0.16914238296314,  0.164276483745833, 0.154684675126265, ...
              0.140642914670651, 0.122555206711478, 0.100942044106287, ...
              0.076425730254889, 0.049714548894970, 0.021616013526483, ...
              0.169142382963144, 0.164276483745833, 0.154684675126265, ...
              0.140642914670651, 0.122555206711478, 0.100942044106287, ...
              0.076425730254889, 0.049714548894970, 0.021616013526483];
      elseif (npg == 19)
        r = [-0.160358645640225,-0.316564099963630,-0.464570741375961, ...
             -0.600545304661681,-0.720966177335229,-0.822714656537143, ...
             -0.903155903614818,-0.960208152134830,-0.992406843843584, ...
              0.000000000000000, 0.160358645640225, 0.316564099963631, ...
              0.464570741375961, 0.600545304661681, 0.720966177335229, ...
              0.822714656537143, 0.903155903614818, 0.960208152134830, ...
              0.992406843843584];
        w = [ 0.158968843393954, 0.152766042065860, 0.142606702173607, ...
              0.128753962539336, 0.111566645547334, 0.091490021622450, ...
              0.069044542737641, 0.044814226765700, 0.019461788229726, ...
              0.161054449848784, 0.158968843393954, 0.152766042065860, ...
              0.142606702173607, 0.128753962539336, 0.111566645547334, ...
              0.091490021622450, 0.069044542737641, 0.044814226765700, ...
              0.019461788229726];
      elseif (npg == 20)
        r = [-0.076526521133497,-0.227785851141645,-0.373706088715419, ...
             -0.510867001950827,-0.636053680726515,-0.746331906460151, ...
             -0.839116971822219,-0.912234428251326,-0.963971927277914, ...
             -0.993128599185095, 0.076526521133497, 0.227785851141645, ...
              0.373706088715419, 0.510867001950827, 0.636053680726515, ...
              0.746331906460151, 0.839116971822219, 0.912234428251326, ...
              0.963971927277914, 0.993128599185095];
        w = [ 0.152753387130726, 0.149172986472604, 0.142096109318382, ...
              0.131688638449177, 0.118194531961518, 0.101930119817240, ...
              0.083276741576705, 0.062672048334109, 0.040601429800387, ...
              0.017614007139152, 0.152753387130726, 0.149172986472604, ...
              0.142096109318382, 0.131688638449177, 0.118194531961518, ...
              0.101930119817240, 0.083276741576705, 0.062672048334109, ...
              0.040601429800387, 0.017614007139152];
      elseif (npg == 21)
        r = [-0.145561854160895,-0.288021316802401,-0.424342120207439, ...
             -0.551618835887220,-0.667138804197412,-0.768439963475678, ...
             -0.853363364583317,-0.920099334150401,-0.967226838566306, ...
             -0.993752170620389, 0.000000000000000, 0.145561854160895, ...
              0.288021316802401, 0.424342120207439, 0.551618835887220, ...
              0.667138804197412, 0.768439963475678, 0.853363364583317, ...
              0.920099334150401, 0.967226838566306, 0.993752170620389];
        w = [ 0.144524403989970, 0.139887394791073, 0.132268938633337, ...
              0.121831416053728, 0.108797299167148, 0.093444423456034, ...
              0.076100113628379, 0.057134425426857, 0.036953789770852, ...
              0.016017228257774, 0.146081133649690, 0.144524403989970, ...
              0.139887394791073, 0.132268938633337, 0.121831416053728, ...
              0.108797299167148, 0.093444423456034, 0.076100113628379, ...
              0.057134425426857, 0.036953789770852, 0.016017228257774];
      elseif (npg == 22)
        r = [-0.069739273319722,-0.207860426688221,-0.341935820892084, ...
             -0.469355837986757,-0.587640403506912,-0.694487263186683, ...
             -0.787816805979208,-0.865812577720300,-0.926956772187174, ...
             -0.970060497835429,-0.994294585482399, 0.069739273319722, ...
              0.207860426688221, 0.341935820892084, 0.469355837986757, ...
              0.587640403506912, 0.694487263186683, 0.787816805979208, ...
              0.865812577720300, 0.926956772187174, 0.970060497835429, ...
              0.994294585482399];
        w = [ 0.139251872855632, 0.136541498346015, 0.131173504787062, ...
              0.123252376810512, 0.112932296080539, 0.100414144442881, ...
              0.085941606217068, 0.069796468424520, 0.052293335152683, ...
              0.033774901584814, 0.014627995298272, 0.139251872855632, ...
              0.136541498346015, 0.131173504787062, 0.123252376810512, ...
              0.112932296080539, 0.100414144442881, 0.085941606217068, ...
              0.069796468424520, 0.052293335152683, 0.033774901584814, ...
              0.014627995298272];
      elseif (npg == 23)
        r = [-0.133256824298466,-0.264135680970345,-0.390301038030291, ...
             -0.509501477846007,-0.619609875763646,-0.718661363131950, ...
             -0.804888401618840,-0.876752358270442,-0.932971086826016, ...
             -0.972542471218115,-0.994769334997552, 0.000000000000000, ...
              0.133256824298466, 0.264135680970345, 0.390301038030291, ...
              0.509501477846007, 0.619609875763646, 0.718661363131950, ...
              0.804888401618840, 0.876752358270442, 0.932971086826016, ...
              0.972542471218115, 0.994769334997552];
        w = [ 0.132462039404697, 0.128905722188082, 0.123049084306729, ...
              0.114996640222411, 0.104892091464541, 0.092915766060035, ...
              0.079281411776719, 0.064232421408526, 0.048037671731085, ...
              0.030988005856979, 0.013411859487142, 0.133654572186106, ...
              0.132462039404697, 0.128905722188082, 0.123049084306729, ...
              0.114996640222411, 0.104892091464541, 0.092915766060035, ...
              0.079281411776719, 0.064232421408526, 0.048037671731085, ...
              0.030988005856979, 0.013411859487142];
      elseif (npg == 24)
        r = [-0.064056892862606,-0.191118867473616,-0.315042679696163, ...
             -0.433793507626045,-0.545421471388840,-0.648093651936975, ...
             -0.740124191578554,-0.820001985973903,-0.886415527004401, ...
             -0.938274552002733,-0.974728555971309,-0.995187219997021, ...
              0.064056892862606, 0.191118867473616, 0.315042679696163, ...
              0.433793507626045, 0.545421471388840, 0.648093651936975, ...
              0.740124191578554, 0.820001985973903, 0.886415527004401, ...
              0.938274552002733, 0.974728555971309, 0.995187219997021];
        w = [ 0.127938195346752, 0.125837456346828, 0.121670472927803, ...
              0.115505668053726, 0.107444270115966, 0.097618652104114, ...
              0.086190161531953, 0.073346481411080, 0.059298584915437, ...
              0.044277438817420, 0.028531388628934, 0.012341229799987, ...
              0.127938195346752, 0.125837456346828, 0.121670472927803, ...
              0.115505668053726, 0.107444270115966, 0.097618652104114, ...
              0.086190161531953, 0.073346481411080, 0.059298584915437, ...
              0.044277438817420, 0.028531388628934, 0.012341229799987];
      elseif (npg == 25)
        r = [-0.122864692610710,-0.243866883720988,-0.361172305809388, ...
             -0.473002731445715,-0.577662930241223,-0.673566368473468, ...
             -0.759259263037357,-0.833442628760834,-0.894991997878275, ...
             -0.942974571228974,-0.976663921459517,-0.995556969790498, ...
              0.000000000000000, 0.122864692610710, 0.243866883720988, ...
              0.361172305809388, 0.473002731445715, 0.577662930241223, ...
              0.673566368473468, 0.759259263037358, 0.833442628760834, ...
              0.894991997878275, 0.942974571228974, 0.976663921459517, ...
              0.995556969790498];
        w = [ 0.122242442990310, 0.119455763535785, 0.114858259145712, ...
              0.108519624474264, 0.100535949067051, 0.091028261982964, ...
              0.080140700335001, 0.068038333812357, 0.054904695975835, ...
              0.040939156701306, 0.026354986615032, 0.011393798501026, ...
              0.123176053726715, 0.122242442990310, 0.119455763535785, ...
              0.114858259145712, 0.108519624474264, 0.100535949067051, ...
              0.091028261982964, 0.080140700335001, 0.068038333812357, ...
              0.054904695975835, 0.040939156701306, 0.026354986615032, ...
              0.011393798501026];*/
    else
    {
        cout << "Maximum number of gauss points is 10." << endl;
        exit(0);
    }
}

// ============================== VolumeDist ===============================

void cFGM :: VolumeDist(eFGMVolDist VolDistType, int ngauss, cVector t, cVector &Vc, double pindex)
{
  if (VolDistType == POWER_LAW)
  {
      PowerLaw(ngauss, pindex, t, Vc);
  }
  else if (VolDistType == SIGMOYD)
  {
      Sigmoyd(ngauss, pindex, t, Vc);
  }
  else if (VolDistType == POLYNOMIAL)
  {
      Polynomial(ngauss, t, Vc);
  }
  else
  {
      cout << "Unknown fgm volume fraction distribution type." << endl;
      exit(0);
  }
}

// ============================== VolumeDist ===============================

void cFGM :: VolumeDist(eFGMVolDist VolDistType, cVector V, int ngauss, cVector t, cVector &Vc)
{

  if (VolDistType == PCHIP)
  {
      PiecewiseCubicInterpolation(V, ngauss, t, Vc);
  }
  else if (VolDistType == BSPLINE)
  {
      Bspline(V, ngauss, t, Vc);
  }
  else{
      cout << "Unknown fgm volume fraction distribution type." << endl;
      exit(0);
  }
}

// ============================== PowerLaw ===============================

void cFGM :: PowerLaw(int ngauss, double pindex, cVector t, cVector &Vc)
{
    Vc.Resize(ngauss);

    double h = 1.0;

    cVector ttemp(ngauss);

    for(int i = 0; i < ngauss; i++)
    {
        ttemp[i] = t[i]    ;                     // Bottom: -0.5 Top: 0.5
        Vc[i] = 1*pow((2*t[i]+h)/(2*h), pindex);
    }
}

// ============================== Sigmoyd ===============================

void cFGM :: Sigmoyd(int ngauss, double pindex, cVector t, cVector &Vc)
{
    Vc.Resize(ngauss);

    double h = 1.0;

    cVector ttemp(ngauss);

    double aux;

    for(int i = 0; i < ngauss; i++)
    {
        ttemp[i] = t[i]    ;                     // Bottom: -0.5 Top: 0.5
        if (ttemp[i] <= 0.0)
        {
            Vc[i] = (1.0/2.0)*pow(((h/2.0)+ttemp[i])/(h/2.0), pindex);
          //  cout << "t: " << ttemp[i] << " Vc: " << Vc[i] << endl;
        }
        else
        {
            if (ttemp[i] > 0.498)
            {
                Vc[i] = 1.0;
            }
            else
            {
                    aux = 1.0-(1.0/2.0)*pow(((h/2.0)-ttemp[i])/(h/2.0), pindex);
                    Vc[i] = aux;
                    }

            //    cout << "tu: " << ttemp[i] << " Vc: " << aux << endl;
        }
    }
}

// ============================== Polynomial ===============================

void cFGM :: Polynomial(int ngauss, cVector t, cVector &Vc)
{
    Vc.Resize(ngauss);
    double a = -1.00;
    double b = -1.00;
    double c = 0.00;
    double n = 1.00;

    for(int i = 0; i < ngauss; i++) Vc[i] = pow(1 - a*t[i] + b*pow(t[i],c), n);
}

// ============================== PCHIP ===============================

void cFGM :: PiecewiseCubicInterpolation(cVector V, int np, cVector t, cVector& Vc)
{
    int ncp = V.Dim();                                               // Number of control points
    double h   = 1/((double)ncp - 1);                                // Thickness between each control point

    cVector X(ncp), S(ncp);
    X[0] = -0.5;
    for (int i = 0; i < ncp; i++){
        if (i != 0){
            X[i] = X[i-1] + h;                                       // Position of control points (should be evenly spaced)
        }

        if ( i == 0){
            S[i] = (-3*V[i] + 4*V[i + 1] - V[i + 2])/(2*h);          // If S1
        }
        else if(i == (ncp - 1)){
            S[i] = ( 3*V[i] - 4*V[i - 1] + V[i - 2])/(2*h);          // If Sn
        }
        else{
            S[i] = (V[i + 1] - V[i - 1])/(2*h);                      // Otherwise
        }
    }

    // Fitting the derivatives to its lower and upper bounds

    cVector lowS(ncp);
    cVector uppS(ncp);

    lowS[0]       = -3*V[0]/h;
    lowS[ncp - 1] = -3*(1 - V[ncp - 1])/h;
    for (int i = 1; i < (ncp - 1); i++){
        cVector lowtemp(2);
        lowtemp[0] = -3*V[i]/h;
        lowtemp[1] = -3*(1 - V[i])/h;
        lowS[i] = lowtemp.Max();
    }

    uppS[0]       = 3*(1 - V[0])/h;
    uppS[ncp - 1] = 3*V[ncp - 1]/h;
    for (int i = 1; i < (ncp - 1); i++){
        cVector upptemp(2);
        upptemp[0] = 3*(1 - V[i])/h;
        upptemp[1] = 3*V[i]/h;
        uppS[i] = upptemp.Min();
    }

    for (int i = 0; i < ncp; i++){
        if(lowS[i] > S[i]){
            S[i] = lowS[i];
        }
        if(uppS[i] < S[i]){
            S[i] = uppS[i];
        }
    }

    // Assessment of volume fraction on evaluated points

    Vc.Resize(np);

    for (int i = 0; i < np; i++){
        int l = 0;
        int r = 0;
        for (int j = 1; j < ncp; j++){
            if(t[i] < X[j]){                                         // Defining which control points will be used
                l = j - 1;
                r = j    ;
                break;
            }
        }

        double dist = t[i] - X[l];
        double f1 = V[l];
        double f2 = V[r];
        double d1 = S[l];
        double d2 = S[r];

        double df = (f2 - f1)/h;
        double c2 = - ( 2.0 * d1 - 3.0 * df + d2 ) / h;
        double c3 =   (       d1 - 2.0 * df + d2 ) / h / h;

        Vc[i] = f1 + ( dist ) * ( d1 + ( dist ) * ( c2+ ( dist ) *   c3 ) );
    }
}

// ============================== Bspline ===============================

void cFGM :: Bspline(cVector V, int nt, cVector t, cVector& Vc)
{
    for (int i = 0; i < nt; i++){
        t[i] = t[i] + 0.5;
    }

    int n = V.Dim();            // Number of Control Points
    int p = 3;                  // Deegree = 3 (Cubic)
    int m = n + p + 1;

    cVector U( m );     // Knot vector
    for (int i = 0; i < m; i++){
        if (i < (p + 1)) U[i] = 0;
        else if (i >= n) U[i] = 1;
        else{
            U[i] = U[i - 1] + 1.0/(n - p);
        }
    }

    Vc.Resize(nt);
    cVector C(nt);

    for (int i = 0; i < nt; i++){
        CurvePoint( n, p, U, V, t[i], C[i] );
        Vc[i] = C[i];
    }

   // cout << "Vc: " << Vc[0] << endl;
}

// ============================== FindSpan ===============================

int cFGM :: FindSpan(int n, int p, double u, cVector U)
{
    if ( u == U[n+1])  return(n);
    int low = p;
    int high = n + 1;
    int mid = (low + high)/2;
    while ( u < U[mid] || u >= U[mid + 1] ){
        if (u < U[mid]) high = mid;
        else            low  = mid;
        mid = (low + high)/2;
    }
    return (mid);
}

// ============================== FindSpan ===============================

void cFGM :: BasisFuns(int i, double u, int p, cVector U, cVector &N)
{
    N[0] = 1.0;
    cVector left(p + 1), right(p + 1);
    for (int j = 1; j <= p; j++){
        left[j]  = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for ( int r = 0; r < j; r++){
            double temp = N[r]/(right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1]*temp;
            saved = left[j - r]*temp;
        }
        N[j] = saved;
    }
}

// ============================== FindSpan ===============================

void cFGM :: CurvePoint(int n, int p, cVector U, cVector P, double u, double& C)
{
    int span = FindSpan( n, p, u, U );
    cVector N(n);
    BasisFuns( span, u, p, U, N );
    C = 0.0;
    for ( int i = 0; i <= p; i++ ){
        C = C + N[i]*P[span - p + i];
    }
}

// ============================= BSplineDer =============================

void cFGM :: BSplineDer(cVector V, int nt, cVector t, cVector& BSder)
{
    int n = V.Dim();            // Number of Control Points
    int p = 3;                  // Deegree = 3 (Cubic)
    int m = n + p + 1;

    cVector U( m );     // Knot vector
    for (int i = 0; i < m; i++){
        if (i < (p + 1)) U[i] = 0;
        else if (i >= n) U[i] = 1;
        else{
            U[i] = U[i - 1] + 1.0/(n - p);
        }
    }

    cVector ders(2);

    for (int i = 0; i < nt; i++){
        CurveDerivs( n, p, U, V, t[i], ders );
        /*cout << "for " << xval[i] << "    ";
        ders.Print();
        cout << "----------------------------------" << endl;*/
        BSder[i] = ders[1];
    }
}

// ============================== FindSpan ===============================

void cFGM :: CurveDerivs(int n, int p, cVector U, cVector P, double u, cVector&CK)
{
    int d = 1;
    int du = 0;
    if (d < p){
        du = d;
    }
    else {
        du = p;
    }
    for (int k = p+1; k <= d; k++) CK[k] = 0.0;

    int span = FindSpan( n, p, u, U );

    cMatrix nders(du+1,p+1);
    DersBasisFuns(span,u,p,du,U,nders);

    for (int k=0; k<=du; k++)
    {
    CK[k] = 0.0;
      for (int j=0; j<=p; j++){
          CK[k] = CK[k] + nders[k][j]*P[span-p+j];
      }
    }
}

// ============================== FindSpan ===============================

void cFGM :: DersBasisFuns(int i, double u, int p, int n, cVector U, cMatrix &ders)
{
    double ndu[p+1][p+1];
    ndu[0][0] = 1.0;
    cVector left(p + 1), right(p + 1);

    for (int j = 1; j <= p; j++)
    {
        left[j] = u-U[i+1-j];
        right[j] = U[i+j]-u;
        double saved = 0.0;
        for (int r = 0; r < j; r++)
        {
          ndu[j][r] = right[r+1] + left[j-r] ;
          double temp = ndu[r][j-1]/ndu[j][r];

          ndu[r][j] = saved + right[r + 1]*temp;
          saved = left[j-r]*temp;
        }
        ndu[j][j] = saved;
    }

    for (int j=0; j<=p; j++){
        ders[0][j] = ndu[j][p];
    }

    double a[2][p+1];
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < (p+1); j++){
            a[i][j] = 0.0;
        }
    }

    for (int r = 0; r <= p; r++)
    {
        int s1 = 0;
        int s2 = 1;

        a[0][0] = 1.0;

        for (int k = 1; k <= n; k++)
        {
            double d = 0.0;
            int rk = r - k;
            int pk = p - k;
            if(r >= k){
                a[s2][0] = a[s1][0]/ndu[pk+1][rk];
                d = a[s2][0]*ndu[rk][pk];
            }

            int j1, j2;

            if (rk >= -1) {
                j1 = 1;
            }
            else j1 = -rk;

            if ( (r-1) <= pk) {
                j2 = k-1;
            }
            else j2 = p-r;

            for (int j=j1; j<=j2; j++)
            {
              a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
              d = d + a[s2][j]*ndu[rk+j][pk];
            }

            if (r <= pk)
            {
            a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
            d = d + a[s2][k]*ndu[r][pk];
            }
            ders[k][r] = d;
            int j=s1;
            s1=s2;
            s2=j;
        }
        int varr = p;
        for (int k = 1; k <= n; k++) {
            //cout << "k = " << k << endl;
            for (int j = r; j<= p; j++){
                //cout << "j = " << j << endl;
                //cout << "derskj antes = " << ders[k][j] << endl;
                ders[k][j] *= varr;
                //cout << "derskj depois = " << ders[k][j] << endl;
            }
            varr *= (p-k);
        }
    }

    /*cout << "ders = " << endl;
    for (int i = 0; i <= 1; i++){
        for (int j = 0; j <= p; j++)
        {
             cout << ders[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;*/
}

// ============================== FGMModel ===============================

void cFGM :: EffPropModel(eFGMModel ModelType, cVector Vcpg, cVector &Epg,
                      cVector &Nupg, cVector &Kpg, cVector &Gpg, cVector &Rhopg)
{
  if (ModelType == VOIGT)    // Rule of Mixtures
  {
      Voigt(Vcpg, Epg, Nupg, Kpg, Gpg);
  }
  else if (ModelType == REUSS)
  {
      Reuss(Vcpg, Epg, Nupg, Kpg, Gpg);
  }
  else if (ModelType == MORI_TANAKA)
  {
      MoriTanaka(Vcpg, Epg, Nupg, Kpg, Gpg);
  }
  else
  {
      cout << "Unknown fgm model type." << endl;
      exit(0);
  }

  VoigtDensity(Vcpg, Rhopg);
}

// ============================== Voigt ===============================

void cFGM :: Voigt(cVector Vcpg, cVector &Epg, cVector &Nupg, cVector &Kpg,
                   cVector &Gpg)
{
    // cout << "ENTROU AQUI VOIGT" << endl;

    int ngauss = Vcpg.Dim();
    Epg.Resize(ngauss);
    Nupg.Resize(ngauss);
    Kpg.Resize(ngauss);
    Gpg.Resize(ngauss);

    /*cout << "Em = " << FGMMat[0] << endl;
    cout << "Ec = " << FGMMat[2] << endl;
    cout << "Num = " << FGMMat[1] << endl;
    cout << "Nuc = " << FGMMat[3] << endl;*/

    double Em = FGMMat[0];
    double Ec = FGMMat[3];
    double Num = FGMMat[1];
    double Nuc = FGMMat[4];

    //FGMMat.Print();

    for (int i = 0; i < ngauss; i++)
    {
        Epg[i] = Em + (Ec - Em)*Vcpg[i];
        Nupg[i] = Num + (Nuc - Num)*Vcpg[i];
        Kpg[i] = Epg[i]/(3*(1 - 2*Nupg[i]));
        Gpg[i] = Epg[i]/(2*(1+Nupg[i]));
    }

  //  cout << "Epg: " << Epg[0] << endl;
}

// ========================== VoigtDensity ============================

void cFGM :: VoigtDensity(cVector Vcpg, cVector &Rhopg)
{
    //cout << "DENSITY" << endl;

    int ngauss = Vcpg.Dim();
    Rhopg.Resize(ngauss);

    double Rhom = FGMMat[2];
    double Rhoc = FGMMat[5];

    //FGMMat.Print();

    for (int i = 0; i < ngauss; i++)
    {
        Rhopg[i] = Rhom + (Rhoc - Rhom)*Vcpg[i];
    }
}

// ============================== Reuss ===============================

void cFGM :: Reuss(cVector Vcpg, cVector &Epg, cVector &Nupg, cVector &Kpg,
                   cVector &Gpg)
{
    int ngauss = Vcpg.Dim();
    Epg.Resize(ngauss);
    Nupg.Resize(ngauss);
    Kpg.Resize(ngauss);
    Gpg.Resize(ngauss);

    double Em = FGMMat[0];
    double Ec = FGMMat[3];
    double Num = FGMMat[1];
    double Nuc = FGMMat[4];

    for (int i = 0; i < ngauss; i++)
    {
        Epg[i] = Em*Ec/(Em*Vcpg[i] + Ec*(1 - Vcpg[i]));
        Nupg[i] = Num*Nuc/(Num*Vcpg[i] + Nuc*(1 - Vcpg[i]));
        Kpg[i] = Epg[i]/(3*(1 - 2*Nupg[i]));
        Gpg[i] = Epg[i]/(2*(1+Nupg[i]));
    }
}

// ============================== MoriTanaka ===============================

void cFGM :: MoriTanaka(cVector Vcpg, cVector &Epg, cVector &Nupg, cVector &Kpg,
                   cVector &Gpg)
{
    int ngauss = Vcpg.Dim();
    Epg.Resize(ngauss);
    Nupg.Resize(ngauss);
    Kpg.Resize(ngauss);
    Gpg.Resize(ngauss);

    double Em = FGMMat[0];
    double Ec = FGMMat[3];
    double Num = FGMMat[1];
    double Nuc = FGMMat[4];

    double Km = Em/(3*(1 - 2*Num));
    double Kc = Ec/(3*(1 - 2*Nuc));
    double Gm = Em/(2*(1 + Num));
    double Gc = Ec/(2*(1 + Nuc));

    cVector Vmpg(ngauss);

    for (int i = 0; i < ngauss; i++) Vmpg[i] = 1 - Vcpg[i];

    double fm = Gm*(9*Km + 8*Gm)/(6*(Km + 2*Gm));

    for (int i = 0; i < ngauss; i++)
    {
        Kpg[i]  = Km + Vcpg[i]/(1/(Kc - Km) + Vmpg[i]/(Km + 4*Gm/3));
        Gpg[i]  = Gm + Vcpg[i]/(1/(Gc - Gm) + Vmpg[i]/(Gm + fm));

        Epg[i] = (9*Kpg[i]*Gpg[i])/(3*Kpg[i] + Gpg[i]);
        Nupg[i] = (3*Kpg[i] - 2*Gpg[i])/(2*(3*Kpg[i] + Gpg[i]));
    }
}

// ============================== EvalVolumeRatio =================================

void cFGM :: EvalVolumeRatio(double exp, double &v)
{
    cVector t(201);

    double lowt = -0.5;
    double uppt =  0.5;

    t[0] = lowt;
    for (int i = 1; i < 201; i++) t[i]  = t[i - 1] + (uppt - lowt)/200;

    cVector Vcpg;

    VolumeDist(FGMVolDist, 201, t, Vcpg, exp);

    double vtemp = 0;

    for (int i = 0; i < 201; i++){
        vtemp += Vcpg[i]/201;
    }

    v = vtemp;


}

// ============================== EvalVolumeRatio =================================

void cFGM :: EvalVolumeRatio(cVector Vcp, double &v)
{
    int ngauss = 25;
    cVector r, w;
    GaussPts1D(ngauss, r, w);

    cVector t(ngauss);
    for (int i = 0; i < ngauss; i++) t[i]  = r[i]/2;

    cVector Vcpg;

    VolumeDist(FGMVolDist, Vcp, ngauss, t, Vcpg);

    double vtemp = 0;

    for (int i = 0; i < ngauss; i++){
        vtemp += Vcpg[i]*w[i];
    }

    vtemp = vtemp/2;

    v = vtemp;
}

// ============================== EvalVolumeRatio3D =================================

void cFGM :: EvalVolumeRatio3D(cVector Vcp, double &v, int Nbx, int Nby, int Nbz)
{
    int ngauss = 5;

    int ngaussthk = 2;

    cVector r, w;
    GaussPts1D(ngauss, r, w);

    cVector rthk, wthk;
    GaussPts1D(ngaussthk, rthk, wthk);

    //r.Print();
    //w.Print();

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

    int n    = ngauss;
    int nthk = ngaussthk;

    cVector V2(1);
    V2[0] = 0.0;
    int NumTot = n*n*n;
    cVector Vcpg(NumTot);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < nthk; k++)
            {
                cMatrix cdnt(3,1);
                cdnt[0][0] = r[i];
                cdnt[1][0] = r[j];
                cdnt[2][0] = rthk[k];
                BsplineSol(CP, 1, cdnt, V2, Nbx, Nby, Nbz, px, py, pz, Ux, Uy, Uz);
                Vcpg[i*n*n + j*n + k] = V2[0];
            }
        }
    }

    v = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < nthk; k++)
            {
                v += Vcpg[i*n*n + j*n + k]*w[i]*w[j]*wthk[k];
            }
        }
    }

    v = v/8.0;

    delete [] CP;
}

// ===============================================================

void cFGM :: EvalDens(double exp, double &m)
{
    double vcratio;

    EvalVolumeRatio(exp, vcratio);

    m = (vcratio)*FGMMat[5] + (1.0-vcratio)*FGMMat[2];
}

// ============================== EvalDens =================================

void cFGM :: EvalDens(cVector Vcp, double &m)
{
    double vcratio;
    EvalVolumeRatio(Vcp, vcratio);

    m = vcratio*FGMMat[5] + (1 - vcratio)*FGMMat[2];
}

// ============================== EvalCost =================================

void cFGM :: EvalCost(double exp, double &c)
{
    double vcratio;

    EvalVolumeRatio(exp, vcratio);

    c = (vcratio)*FGMMat[5]*50.0 + (1.0-vcratio)*FGMMat[2]*3.0;
}


// ============================== EvalYieldStress =================================

double cFGM :: EvalYieldStress(double Vck)
{
    double q = 120.0e9; // Given
    double sigmaym = 400.0e6;  // Stainless steel

    double Em = FGMMat[0];
    double Ec = FGMMat[3];

    double Vmk = 1.0 - Vck;

    double sigmayk = sigmaym*(Vmk + (((q+Em)/(q+Ec))*(Ec/Em))*Vck);

    return sigmayk;
}

// =============================== CalcABDG =================================

void cFGM :: CalcABDG(double t, cVector r, cVector w, cVector Epg, cVector Nupg,
                      cMatrix &A, cMatrix &B, cMatrix &D, cMatrix &G, cMatrix &ABDG)
{
  int ngauss = Epg.Dim();

  // Initialize the given matrices.

  A.Resize(3, 3);
  B.Resize(3, 3);
  D.Resize(3, 3);
  G.Resize(2, 2);
  ABDG.Resize(8, 8);

  A.Zero();
  B.Zero();
  D.Zero();
  G.Zero();
  ABDG.Zero();

  cMatrix Q(3, 3);
  cMatrix Qs(2, 2);

  // Compute the thickness of the FG structure.

  double thk = t;

  // Evaluate the ADBG matrices

  double J = thk/2;              // Jacobian
  double z;

  for (int i = 0; i < ngauss; i++)
  {
      z = r[i]*thk/2;
      QMatrix(Epg[i], Nupg[i], Q, Qs);

      for (int j = 0; j < 3; j++)
      {
          for (int k = 0; k < 3; k++)
          {
              A[j][k] = A[j][k] + J*w[i]*Q[j][k];
              B[j][k] = B[j][k] + J*w[i]*z*Q[j][k];
              D[j][k] = D[j][k] + J*w[i]*z*z*Q[j][k];

              if (j < 2 && k < 2) G[j][k] = G[j][k] + (5.00/6.00)*J*w[i]*Qs[j][k];
          }
      }
  }

  // Mount ABDG matrix

  for (int i = 0; i < 3; i++)
  {
      for (int j = 0; j < 3; j++)
      {
          ABDG[i][j] = A[i][j];
      }
  }

  for (int i = 0; i < 3; i++)
  {
      for (int j = 3; j < 6; j++)
      {
          ABDG[i][j] = B[i][j-3];
      }
  }

  for (int i = 3; i < 6; i++)
  {
      for (int j = 0; j < 3; j++)
      {
          ABDG[i][j] = B[i-3][j];
      }
  }

  for (int i = 3; i < 6; i++)
  {
      for (int j = 3; j < 6; j++)
      {
          ABDG[i][j] = D[i-3][j-3];
      }
  }

  for (int i = 6; i < 8; i++)
  {
      for (int j = 6; j < 8; j++)
      {
          ABDG[i][j] = G[i-6][j-6];
      }
  }
}

// =============================== QMatrix =================================

void cFGM :: QMatrix(double E, double Nu, cMatrix &Q, cMatrix &Qs)
{
  // Matrix [Q] (membrane).

   Q.Zero();
   Qs.Zero();

   Q[0][0] = E/(1 - Nu*Nu);
   Q[1][1] = Q[0][0];
   Q[0][1] = Nu*E/(1 - Nu*Nu);
   Q[1][0] = Q[0][1];
   Q[2][2] = E/(2*(1 + Nu));

   // Matrix [Qs] (transverse shear).

   Qs[1][1] = Q[2][2];
   Qs[0][0] = Q[2][2];
}

// =============================== CalcABDG =================================

void cFGM :: CalcMb(double t, cVector r, cVector w, cVector Rhopg, cMatrix &Mb)
{
  int ngauss = Rhopg.Dim();
  //Rhopg.Print();

  // Initialize the Mb matrix.

  Mb.Resize(5, 5);

  Mb.Zero();

  // Compute the thickness of the FG structure.

  double thk = t;

  // Evaluate the ADBG matrices

  double J = thk/2;              // Jacobian
  double z;
  double I0 = 0;
  double I1 = 0;
  double I2 = 0;


  for (int i = 0; i < ngauss; i++)
  {
      z = r[i]*thk/2;

      I0 = I0 + J*w[i]*Rhopg[i];
      I1 = I1 + J*w[i]*Rhopg[i]*z;
      I2 = I2 + J*w[i]*Rhopg[i]*z*z;
  }

  // Mount Mb matrix

  Mb[0][0] = I0;
  Mb[1][1] = I0;
  Mb[2][2] = I0;
  Mb[3][3] = I2;
  Mb[4][4] = I2;
  Mb[0][4] = I1;
  Mb[1][3] = -I1;
  Mb[3][1] = -I1;
  Mb[4][0] = I1;
}

// ============================== BsplineSol ===============================

void cFGM :: BsplineSol(cMatrix *CP, int nt, cMatrix coord, cVector& Vc, int nbx, int nby, int nbz, int px, int py, int pz, cVector Ux, cVector Uy, cVector Uz)
{
    Vc.Resize(nt);
    cVector C(nt);

    for (int i = 0; i < nt; i++){
        SolidPoint( nbx, nby, nbz, px, py, pz, Ux, Uy, Uz, CP, coord, C[i] );
        Vc[i] = C[i];
    }
}

// ============================== SolidPoint ===============================

void cFGM :: SolidPoint(int nbx, int nby, int nbz, int px, int py, int pz, cVector Ux, cVector Uy, cVector Uz, cMatrix *CP, cMatrix coord, double &C)
{
    int spanx = FindSpan( nbx, px, coord[0][0], Ux);
    cVector Nx(nbx);
    BasisFuns( spanx, coord[0][0], px, Ux, Nx );

    int spany = FindSpan( nby, py, coord[1][0], Uy);
    cVector Ny(nby);
    BasisFuns( spany, coord[1][0], py, Uy, Ny );

    int spanz = FindSpan( nbz, pz, coord[2][0], Uz);
    cVector Nz(nbz);
    BasisFuns( spanz, coord[2][0], pz, Uz, Nz );

    C = 0.0;
    for ( int i = 0; i <= px; i++ ){
        for (int j = 0; j <= py; j++ ){
            for (int m = 0; m <= pz; m++){
                C = C + Nx[i]*Ny[j]*Nz[m]*CP[m][i][j];
            }
        }
    }
}


// ======================================================= End of file =====
