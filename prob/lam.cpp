// -------------------------------------------------------------------------
// lam.cpp - Implementation of the Laminated problem class.
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
// Created:      21-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     08-Apr-2014    Elias Saraiva Barroso
//               Creation of cLamPltBuck and cLamPltLoadFactor classes.
//
// Modified:     10-Apr-2014    Evandro Parente Junior
//               ??
//
// Modified:	 Jul-2014	    Joao Paulo Bernhardt
//		         Added Progressive Failure problem classes.
//
// Modified:     03-Mar-2019    Marina Alves Maia
//               Added Average Density function.
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

#include "lam.h"
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
eLamType cLaminated :: LamType    = GENERAL;
double   cLaminated :: MinAng     = 0.0;
double   cLaminated :: MaxAng     = 0.0;
double   cLaminated :: dAng       = 0.0;
double   cLaminated :: dThk       = 0.0;
double   cLaminated :: MinThk     = 0.0;
double   cLaminated :: MaxThk     = 0.0;
int      cLaminated :: MaxPlyNum  = 0;
cVector  cLaminated :: ThkVal;
cVector  cLaminated :: AngVal;
cVector  cLaminated :: Thk;
cVector  cLaminated :: Ang;
int*     cLaminated :: Mat;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ ReadAngleRange =============================

void cLaminated :: ReadAngleRange(istream &in)
{
  if (!(in >> MinAng) || !(in >> dAng) || !(in >> MaxAng))
  {
    cout << "Error in the input of the lamination angle range." << endl;
    exit(0);
  }
}

// =========================== ReadAngleValues =============================

void cLaminated :: ReadAngleValues(istream &in)
{
  int numang = 0;

  if (!(in >> numang) || numang == 0)
  {
    cout << "Error in the input of the number of possible angle values." << endl;
    exit(0);
  }

  AngVal.Resize(numang);

  for (int i = 0; i < numang; i++)
  {
    if (!(in >> AngVal[i]))
    {
      cout << "Error in the input of the possible angle values." << endl;
      exit(0);
    }
  }
}

// ============================ ReadThickRange =============================

void cLaminated :: ReadThickRange(istream &in)
{
  if (!(in >> MinThk) || !(in >> dThk) || !(in >> MaxThk))
  {
    cout << "Error in the input of the laminate thickness range." << endl;
    exit(0);
  }
}

// =========================== ReadThickValues =============================

void cLaminated :: ReadThickValues(istream &in)
{
  int numthk = 0;

  if (!(in >> numthk) || numthk == 0)
  {
    cout << "Error in the input of the number of possible thickness values." << endl;
    exit(0);
  }

  ThkVal.Resize(numthk);

  for (int i = 0; i < numthk; i++)
  {
    if (!(in >> ThkVal[i]))
    {
      cout << "Error in the input of the possible thickness values." << endl;
      exit(0);
    }
  }
}

// ============================ ReadMaxPlyNum ==============================

void cLaminated :: ReadMaxPlyNum(istream &in)
{
  if (!(in >> MaxPlyNum))
  {
    cout << "Error in the input of the laminate maximum number of plies." << endl;
    exit(0);
  }
}

// ============================== ReadLamType ==============================

void cLaminated :: ReadLamType(istream &in)
{
  char label[100];
  
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of laminate type." << endl;
    exit(0);
  }

  if (string(label) == "Symmetric" || string(label) == "symmetric")
    LamType = SYMMETRIC;
  
  else if (string(label) == "Balanced" || string(label) == "balanced")
    LamType = BALANCED;
  
  else if (string(label) == "Symmetric_Balanced" || string(label) == "symmetric_balanced")
    LamType = SYMMETRIC_BALANCED;
  
  else if (string(label) == "General" || string(label) == "general")
    LamType = GENERAL;
  
  else
  {
    cout << "Unknown laminate type: " << label << endl;
    exit(0);
  }
}

// ============================= cLaminated ================================

cLaminated :: cLaminated( )
{
}

// ============================ ~cLaminated ================================

cLaminated :: ~cLaminated(void)
{
  delete[]ListDim;
  delete[]Mat;
}

// ============================ Init =======================================

void cLaminated :: Init( )
{
  // Add the problem reader functions tags

  NumRow = 3;

  switch (LamType)
  {
    case SYMMETRIC:
    case BALANCED:
      NumCol = (int) MaxPlyNum/2;
      break;
    case SYMMETRIC_BALANCED:
      NumCol = (int) MaxPlyNum/4;
      break;
    case GENERAL:
      NumCol = MaxPlyNum;
      break;
  }

  ListDim = new int[NumRow];

  if (dThk && ThkVal.Dim( ))
    ListDim[0] = 1 + ((MaxThk - MinThk) / dThk) + ThkVal.Dim( );
  
  else if (dThk && !ThkVal.Dim( ))
    ListDim[0] = 1 + (MaxThk - MinThk) / dThk;
  
  else
    ListDim[0] = ThkVal.Dim( );

  if (dAng && AngVal.Dim( ))
    ListDim[1] = 1 + ((MaxAng - MinAng) / dAng) + AngVal.Dim( );
  
  else if (dAng && !AngVal.Dim( ))
    ListDim[1] = 1 + (MaxAng - MinAng) / dAng;
  
  else
    ListDim[1] = AngVal.Dim( );

  ListDim[2] = cMaterial::GetNumOrtho( );

  Thk.Resize(ListDim[0]);
  Ang.Resize(ListDim[1]);
  Mat = new int[cMaterial::GetNumOrtho()];

  if (dThk)
  {
    for (int i = 0; i < 1 + ((MaxThk - MinThk)/dThk); i++)
      Thk[i] = MinThk + i*dThk;
  }

  if (ThkVal.Dim( ))
  {
    for (int i = 0; i < ThkVal.Dim( ); i++)
      Thk[ListDim[0] - 1 - i] = ThkVal[i];
  }
 
  if (dAng)
  {
    for (int i = 0; i < 1 + ((MaxAng - MinAng)/dAng); i++)
      Ang[i] = MinAng + i*dAng;
  }
 
  if (AngVal.Dim( ))
  {
    for (int i = 0; i < AngVal.Dim( ); i++)
      Ang[ListDim[1] - 1 - i] = AngVal[i];
  }
  
  int pos = 0;

  for (int i = 0; i < cMaterial::GetNumMat( ); i++)
  {
    cMaterial *mat = cMaterial::GetMaterial(i+1);
  
    if (mat->GetType( ) == MAT_ORTHOTROPIC)
    {
      Mat[pos] = mat->GetLabel( );
      pos++;
    }
  }

  // Sorting the Thickness vector

  Thk.Sort(ASCENDING); 
}

// ============================== LoadReadFunc =============================

void cLaminated :: LoadReadFunc(cInpMap &im)
{
  // Register read functions.
  im.Insert("LAMINATE.ANGLE.RANGE",makeReadObj(cLaminated,ReadAngleRange));
  im.Insert("LAMINATE.ANGLE.VALUES",makeReadObj(cLaminated,ReadAngleValues));
  im.Insert("LAMINATE.MAXIMUM.PLY.NUMBER",makeReadObj(cLaminated,ReadMaxPlyNum));
  im.Insert("LAMINATE.THICKNESS.RANGE",makeReadObj(cLaminated,ReadThickRange));
  im.Insert("LAMINATE.THICKNESS.VALUES",makeReadObj(cLaminated,ReadThickValues));
  im.Insert("LAMINATE.TYPE",makeReadObj(cLaminated,ReadLamType));

  // Register material read functions.
  im.Insert("MATERIAL",makeRead(cMaterial :: ReadNumMat));
  im.Insert("MATERIAL.DENSITY",makeRead(cMaterial :: ReadDensity));
  im.Insert("MATERIAL.COST",makeRead(cMaterial :: ReadCost));
  im.Insert("MATERIAL.ISOTROPIC",makeRead(cMaterial :: ReadIso));
  im.Insert("MATERIAL.ORTHOTROPIC",makeRead(cMaterial :: ReadOrtho));
  im.Insert("ISOTROPIC.FAILURE.CRITERION",makeRead(cMatIso :: ReadFailType));
  im.Insert("ORTHOTROPIC.FAILURE.CRITERION",makeRead(cMatOrtho :: ReadFailType));
}

// ============================== PrintVar =================================

void cLaminated :: PrintVar(int **algvar)
{
  // Decode the variables.

  cMatrix layup;
  Decode(algvar, layup);

  PrintResult(algvar, cout); 

  // Print the laminate layup.

  for (int i = 0; i < layup.NRow( ); i++)
  {
    for (int j = 0; j < layup.NCol( ); j++)
      cout << layup[i][j] << " ";
    
    cout << endl;
  }
}

// ============================== PrintVar =================================

void cLaminated :: DecodeVar(int **algvar, cVector &xsamp, cMatrix &layup)
{
  // Decode the variables.

  Decode(algvar, layup);

  xsamp.Resize(NumVar);

  for (int i = 0; i < NumVar; i++)
  {
      if (algvar[1][i] == 0)
          xsamp[i] = 0.0;
      else if (algvar[1][i] == 1)
          xsamp[i] = 0.5;
      else
          xsamp[i] = 1.0;
  }
}

// ============================== PrintVar =================================

int cLaminated :: VarNumEff(void)
{
  int low,upp;
  int n = VarNumRow( );
  int m = VarNumCol( );

  double eff = 0;

  for (int i = 0; i < n; i++)
  {
    GetBounds(i, &low, &upp);
    if (low != upp)
    {
        eff += 1.0;
    }
  }
  return n*m*eff/3.0;
}

// ============================== PrintResult ======================

void cLaminated :: PrintResult(int **algvar, ostream &res)
{
  cVector countAng(Ang.Dim());
  ostringstream currline, lamline;
  queue<pair<string,string> > OutPutQueue;

  countAng.Zero();

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

  res << "======================================\n\n";
}

// ============================== WriteVar =================================

void cLaminated :: WriteVar(int **algvar, ostream &out)
{
  // Decode the variables.

  cMatrix layup;
  Decode(algvar, layup);

  //PrintResult(algvar, cout); 

  // Print the laminate layup.

  for (int i = 0; i < layup.NRow( ); i++)
  {
    for (int j = 0; j < layup.NCol( ); j++)
      out << layup[i][j] << " ";
 
    out << endl;
  }
}

// =============================== GetBounds ===============================

void cLaminated :: GetBounds(int i, int *low, int *upp)
{
  *low = 0;
  *upp = ListDim[i]-1;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ================================ Decode =================================

bool cLaminated :: Decode(int **algvar, cMatrix &layup)
{
  int j;

  // Check bounds.

  int MaxMat = cMaterial::GetNumMat( ) - 1;

/*  MaxThk = Thk[Thk.Dim()-1];

  for (j = 0; j < NumCol; j++)
  {
    if (Thk[algvar[0][j]] < -1 || Thk[algvar[0][j]] > MaxThk)
    {
      cout << "Error: thk = " << Thk[algvar[0][j]] << " MaxThk = " << MaxThk << endl;
      return(false);
    }
    if (Ang[algvar[1][j]] < 0 || Ang[algvar[1][j]] > MaxAng)
    {
      cout << "Error: ang = " << Ang[algvar[1][j]] << " MaxAng = " << MaxAng << endl;
      return(false);
    }
    if (algvar[2][j] < 0 || algvar[2][j] > MaxMat)
    {
      cout << "Error:  mat = " << Mat[algvar[2][j]] << " MaxMat = " << MaxMat << endl;
      return(false);
    }
  }*/  // Rewrite this part

  // Count the number of layers with thickness > 0.

  int nvarlay = 0;
  for (j = 0; j < NumCol; j++)
    if (algvar[0][j] != -1 && Thk[algvar[0][j]] != 0.00) nvarlay++;

  if (nvarlay == 0) return(0);  // Empty layup

  // Remove layers with zero thickness.

  int **lamvar = new int*[3];
 
  for (int i = 0; i < 3; i++) lamvar[i] = new int[nvarlay];

  int pos = 0;
 
  for (j = 0; j < NumCol; j++)
    if (algvar[0][j] != -1 && Thk[algvar[0][j]] != 0.00)
    {
      lamvar[0][pos] = algvar[0][j];
      lamvar[1][pos] = algvar[1][j];
      lamvar[2][pos] = algvar[2][j];
      pos++;
    }

  // Get the laminate layup from the optimization variables.

  if (LamType == SYMMETRIC)
  {
    layup.Resize(3, 2*nvarlay);

    for (j = 0; j < nvarlay; j++)
    {
      int jsup = 2*nvarlay - j - 1; // Symmetric of layer j
      layup[0][j] = layup[0][jsup] = Thk[lamvar[0][j]];
      layup[1][j] = layup[1][jsup] = Ang[lamvar[1][j]];
      layup[2][j] = layup[2][jsup] = Mat[lamvar[2][j]];
    }
  }
  else if (LamType == BALANCED)
  {
    layup.Resize(3, 2*nvarlay);

    for (j = 0; j < nvarlay; j++)
    {
      int j1 = 2*j;     //  Angle layer
      int j2 = j1 + 1;  // -Angle layer
      layup[0][j1] = layup[0][j2] = Thk[lamvar[0][j]];
      layup[1][j1] =  Ang[lamvar[1][j]];
      layup[1][j2] = -Ang[lamvar[1][j]];
      layup[2][j1] = layup[2][j2] = Mat[lamvar[2][j]];
    }
  }
  else if (LamType == SYMMETRIC_BALANCED)
  {
    layup.Resize(3, 4*nvarlay);

    for (j = 0; j < nvarlay; j++)
    {
      if (Ang[lamvar[1][j]] == 0 || Ang[lamvar[1][j]] == 90)
      {
      int jinf1 = 2*j;                 //  Angle layer
      int jinf2 = jinf1 + 1;           // -Angle layer
      int jsup1 = 4*nvarlay - 2*j - 1; // Symmetric of layer jinf1
      int jsup2 = jsup1 - 1;           // Symmetric of layer jinf2
      layup[0][jinf1] = layup[0][jinf2] = layup[0][jsup1] = layup[0][jsup2] = (Thk[lamvar[0][j]]);
      layup[1][jinf1] = layup[1][jsup1] =  Ang[lamvar[1][j]];
      layup[1][jinf2] = layup[1][jsup2] = -Ang[lamvar[1][j]];
      layup[2][jinf1] = layup[2][jinf2] = layup[2][jsup1] = layup[2][jsup2] = Mat[lamvar[2][j]];
      }
      else
      {
      int jinf1 = 2*j;                 //  Angle layer
      int jinf2 = jinf1 + 1;           // -Angle layer
      int jsup1 = 4*nvarlay - 2*j - 1; // Symmetric of layer jinf1
      int jsup2 = jsup1 - 1;           // Symmetric of layer jinf2
      layup[0][jinf1] = layup[0][jinf2] = layup[0][jsup1] = layup[0][jsup2] = Thk[lamvar[0][j]];
      layup[1][jinf1] = layup[1][jsup1] =  Ang[lamvar[1][j]];
      layup[1][jinf2] = layup[1][jsup2] = -Ang[lamvar[1][j]];
      layup[2][jinf1] = layup[2][jinf2] = layup[2][jsup1] = layup[2][jsup2] = Mat[lamvar[2][j]];
      }

    }
  }
  else if (LamType == GENERAL)
  {
    layup.Resize(3, nvarlay);

    for (j = 0; j < nvarlay; j++)
    {
      layup[0][j] = Thk[lamvar[0][j]];
      layup[1][j] = Ang[lamvar[1][j]];
      layup[2][j] = Mat[lamvar[2][j]];
    }
  }

  // Release memory.

  for (int i = 0; i < 3; i++) delete []lamvar[i];
  delete []lamvar;

  return(1);
}

// ============================== MaxContLay ===============================

int cLaminated :: MaxContLay(cMatrix &layup)
{
  // Compute the maximum number of contiguous layers with the same angle.

  int nlam = layup.NCol( );
  int numcont = 1;
  int maxcont = numcont;
  double angle = layup[1][0];
 
  if (angle == -90) angle = 90;
 
  for (int lam = 1; lam < nlam; lam++)
  {
    double nxtang = layup[1][lam];
 
    if (nxtang == -90) nxtang = 90;

    if (nxtang == angle)
    {
      numcont++;
      if (numcont > maxcont) maxcont = numcont;
    }
 
    else
    {
      numcont = 1;
      angle = nxtang;
    }
  }

  return(maxcont);
}

// ============================== MaxContThk ===============================

double cLaminated :: MaxContThk(cMatrix &layup)
{
  // Compute the thickness of contiguous layers with the same angle.

  int nlam = layup.NCol( );
  double thkcont = layup[0][0];  // Thickness of first layer
  double maxcont = thkcont;
  double angle = layup[1][0];
 
  if (angle == -90) angle = 90;
 
  for (int lam = 1; lam < nlam; lam++)
  {
    double nxtang = layup[1][lam];
 
    if (nxtang == -90) nxtang = 90;

    if (nxtang == angle)
    {
      thkcont += layup[0][lam];  // Update the thickness of cont. layers
 
      if (thkcont > maxcont) maxcont = thkcont;
    }
 
    else
    {
      thkcont = layup[0][lam];   // Reset the thickness of cont. layers
      angle = nxtang;
    }
  }

  return(maxcont);
}

// ============================= GetThickness ==============================

double cLaminated :: GetThickness(cMatrix &layup)
{
  int nlam = layup.NCol( );

  double tt = 0.0;
 
  for (int lam = 0; lam < nlam; lam++) tt += layup[0][lam];

  return(tt);
}

// ============================== GetAverageDensity ========================

double cLaminated :: GetAverageDensity(cMatrix &layup)
{
    int nlam = layup.NCol( );

    double dens = 0.0;
    double avdens = 0.0;
    int lanmat;

    for (int lam = 0; lam < nlam; lam++)
    {
        lanmat = (int)layup[2][lam];
        dens = dens + cMaterial::GetMaterial(lanmat)->GetDensity( );
    }

    avdens = dens/nlam;

//    cout << "Densidade media " << avdens << endl;

    return(avdens);
}

// =============================== TMatrix =================================

void cLaminated :: TMatrix(double theta, cMatrix &T)
{
  // Trigonometric constants.

  double c  = cos(theta);
  double s  = sin(theta);
  double c2 = c*c;
  double s2 = s*s;

  // Transformation matrix (in-plane).

  T[0][0] = T[1][1] = c2;
  T[0][1] = T[1][0] = s2;
  T[0][2] = c*s;
  T[1][2] = -c*s;
  T[2][0] = -2.0*c*s;
  T[2][1] = 2.0*c*s;
  T[2][2] = c2 - s2;
}

// =============================== TMatrix =================================

void cLaminated :: TMatrix(double theta, cMatrix &Tm, cMatrix &Ts)
{
  // Trigonometric constants.

  double c  = cos(theta);
  double s  = sin(theta);
  double c2 = c*c;
  double s2 = s*s;

  // Membrane (in-plane) transformation matrix.

  Tm[0][0] = Tm[1][1] = c2;
  Tm[0][1] = Tm[1][0] = s2;
  Tm[0][2] = c*s;
  Tm[1][2] = -c*s;
  Tm[2][0] = -2.0*c*s;
  Tm[2][1] = 2.0*c*s;
  Tm[2][2] = c2 - s2;

  // Shear transformation matrix.

  Ts[0][0] = Ts[1][1] = c;
  Ts[0][1] = -s;
  Ts[1][0] = s;
}

// =============================== QlMatrix ================================

void cLaminated :: QlMatrix(double *param, cMatrix &Ql)
{
  // Get material parameters.

  double E1   = param[0];
  double E2   = param[1];
  double Nu12 = param[3];
  double Nu21 = Nu12*E2/E1;
  double G12  = param[6];

  // Evaluate the local constitutive matrix (in-plane).

  Ql[0][0] = E1/(1.0 - Nu12*Nu21);
  Ql[0][1] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Ql[0][2] = 0.0;

  Ql[1][0] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Ql[1][1] = E2/(1.0 - Nu12*Nu21);
  Ql[1][2] = 0.0;

  Ql[2][0] = 0.0;
  Ql[2][1] = 0.0;
  Ql[2][2] = G12;
}

// =============================== QlMatrix with Failure Index ================================

void cLaminated :: QlMatrix(double *param, cMatrix &Ql, int FLi)
{
  // Get material parameters.
 
  double df = 0.0001; 	// degradation factor
  double E1;
  double E2 ;
  double Nu12;
  double Nu21;
  double G12;

  switch(FLi)
  {case 1:	//Fiber Failure
	   E1   = param[0]*df;
	   E2   = param[1]*df;
	   Nu12 = param[3]*df;
	   Nu21 = Nu12*E2/E1;
	   G12  = param[6]*df;
      break;
	case 2:	//Matrix Failure
	   E1   = param[0];
	   E2   = param[1]*df;
	   Nu12 = param[3]*df;
	   Nu21 = Nu12*E2/E1;
	   G12  = param[6]*df;
   	break;
	case 3:	//Shear Failure
	   E1   = param[0];
	   E2   = param[1]*df;
	   Nu12 = param[3]*df;
	   Nu21 = Nu12*E2/E1;
	   G12  = param[6]*df;
		break;
	default:	//No Failure
	   E1   = param[0];
	   E2   = param[1];
	   Nu12 = param[3];
	   Nu21 = Nu12*E2/E1;
	   G12  = param[6];
		break;
  }
 
  // Evaluate the local constitutive matrix (in-plane).

  Ql[0][0] = E1/(1.0 - Nu12*Nu21);
  Ql[0][1] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Ql[0][2] = 0.0;

  Ql[1][0] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Ql[1][1] = E2/(1.0 - Nu12*Nu21);
  Ql[1][2] = 0.0;

  Ql[2][0] = 0.0;
  Ql[2][1] = 0.0;
  Ql[2][2] = G12;
}

// =============================== QlMatrix ================================

void cLaminated :: QlMatrix(double *param, cMatrix &Qlm, cMatrix &Qls)
{
  // Get material parameters.

  double E1   = param[0];
  double E2   = param[1];
  double Nu12 = param[3];
  double Nu21 = Nu12*E2/E1;
  double G12  = param[6];
  double G13  = param[7];
  double G23  = param[8];

  // Evaluate the local membrane constitutive matrix (Qlm).

  Qlm[0][0] = E1/(1.0 - Nu12*Nu21);
  Qlm[0][1] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Qlm[0][2] = 0.0;

  Qlm[1][0] = (Nu12*E2)/(1.0 - Nu12*Nu21);
  Qlm[1][1] = E2/(1.0 - Nu12*Nu21);
  Qlm[1][2] = 0.0;

  Qlm[2][0] = 0.0;
  Qlm[2][1] = 0.0;
  Qlm[2][2] = G12;

  // Evaluate the local shear constitutive matrix (Qls).

  Qls[0][0] = G23;
  Qls[0][1] = 0.0;

  Qls[1][0] = 0.0;
  Qls[1][1] = G13;
}

// =============================== QgMatrix ================================

void cLaminated :: QgMatrix(double theta, cMatrix &Ql, cMatrix &Qg)
{
  // Evaluate auxiliary constants.

  double c  = cos(theta);
  double s  = sin(theta);
  double c2 = c*c;
  double s2 = s*s;
  double c3 = c*c2;
  double s3 = s*s2;
  double c4 = c2*c2;
  double s4 = s2*s2;

  // Evaluate the global constitutive matrix [Qg].

  Qg[0][0] = c4*Ql[0][0] + 2*s2*c2*(Ql[0][1] + 2*Ql[2][2]) + s4*Ql[1][1];
  Qg[0][1] = s2*c2*(Ql[0][0] + Ql[1][1] - 4*Ql[2][2]) + (c4 + s4)*Ql[0][1];
  Qg[0][2] = s*c3*(Ql[0][0] - Ql[0][1] - 2*Ql[2][2]) + s3*c*(Ql[0][1] - Ql[1][1] + 2*Ql[2][2]);

  Qg[1][0] = Qg[0][1];
  Qg[1][1] = s4*Ql[0][0] + 2*s2*c2*(Ql[0][1] + 2*Ql[2][2]) + c4*Ql[1][1];
  Qg[1][2] = s3*c*(Ql[0][0] - Ql[0][1] - 2*Ql[2][2]) + s*c3*(Ql[0][1] - Ql[1][1] + 2*Ql[2][2]);

  Qg[2][0] = Qg[0][2];
  Qg[2][1] = Qg[1][2];
  Qg[2][2] = s2*c2*(Ql[0][0] + Ql[1][1] - 2*Ql[0][1] - 2*Ql[2][2]) + (c4 + s4)*Ql[2][2];
}

// ================================ CalcA ==================================

void cLaminated :: CalcA(cMatrix &layup, cMatrix &A)
{
  // Initialize the membrane stiffness matrix.

  A.Zero( );

  // Evaluate the contribution of each ply.

  int nlam = layup.NCol( );
  cMaterial *prevmat = 0;        // Previous material
  cMatrix Ql(3,3);         // Local constitutive matrix
  cMatrix Qg(3,3);         // Global constitutive matrix
 
  for (int lam = 0; lam < nlam; lam++)
  {
    // Evaluate the local constitutive matrix [Ql]. This matrix is computed
    // only when the layer material is different from the previous one.

    cMaterial *mat = cMaterial::GetMaterial((int)layup[2][lam]);

    if (mat != prevmat)
    {
      double *param = new double[9];

       eMatType matid  = mat->GetType( );

      if (matid == MAT_ISOTROPIC)
      {
          double *paramtemp = new double[mat->NumParam( )];
          mat->GetParam(paramtemp);

          param[0] = paramtemp[0];      // E1
          param[1] = paramtemp[0];      // E2
          param[3] = paramtemp[1];      // v
          param[6] = paramtemp[0]/(2*(1+paramtemp[1]));  // G

      //    cout << "\nE1 " << param[0] << " E2 " << param[1] << " v " << param[3] << " G12 " << param[6] << endl;

          delete []paramtemp;
      }
      else
      {
          mat->GetParam(param);
    //      cout << "E1 " << param[0] << " E2 " << param[1] << " v " << param[3] << " G12 " << param[6] << endl;
      }

      QlMatrix(param, Ql);
      delete []param;
      prevmat = mat;
    }

    // Evaluate the global constitutive matrix [Qg].

    double theta = layup[1][lam]*PI/180.0;  // deg to rad
    QgMatrix(theta, Ql, Qg);

    // Compute [A] matrix.

    double dz = layup[0][lam];  // Layer thickness
 
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) A[i][j] += Qg[i][j]*dz;
  }
}

// ================================ CalcD ==================================

void cLaminated :: CalcD(cMatrix &layup, cMatrix &D)
{
  // Initialize [D] matrix.

  D.Zero( );

  // Compute the laminate thickness.

  double lamthk = GetThickness(layup);

  // Evaluate the contribution of each ply.

  int nlam = layup.NCol( );
  int prevmat = -1;        // Previous material
  double zb = -lamthk/2.0; // Bottom coordinate (layer)
  double zt;               // Top coordinate (layer)
  cMatrix Ql(3,3);         // Local constitutive matrix
  cMatrix Qg(3,3);         // Global constitutive matrix
 
  for (int lam = 0; lam < nlam; lam++)
  {
    // Evaluate the local constitutive matrix [Ql]. This matrix is computed
    // only when the layer material is different from the previous one.

    int idmat = (int)layup[2][lam];
    if (idmat != prevmat)
    {
      cMaterial *mat = cMaterial::GetMaterial(idmat);
      double *param = new double[mat->NumParam( )];
      mat->GetParam(param);
      QlMatrix(param, Ql);
      delete []param;
      prevmat = idmat;
    }

    // Evaluate the global constitutive matrix [Qg].

    double theta = layup[1][lam]*PI/180.0;  // deg to rad
    QgMatrix(theta, Ql, Qg);

    // Compute [D] matrix.

    zt = zb + layup[0][lam];
    double d3z = zt*zt*zt - zb*zb*zb;
 
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) D[i][j] += Qg[i][j]*d3z/3.0;

    // Update the bottom coordinate.

    zb = zt;
  }
}

// =============================== CalcAD ==================================

void cLaminated :: CalcAD(cMatrix &layup, cMatrix &A, cMatrix &D)
{
  // Initialize given matrices.

  A.Zero( );
  D.Zero( );

  // Compute the laminate thickness.

  double lamthk = GetThickness(layup);

  // Evaluate the contribution of each ply.

  int nlam = layup.NCol( );
  int prevmat = -1;        // Previous material
  double zb = -lamthk/2.0; // Bottom coordinate (layer)
  double zt;               // Top coordinate (layer)
  cMatrix Ql(3,3);         // Local constitutive matrix
  cMatrix Qg(3,3);         // Global constitutive matrix
 
  for (int lam = 0; lam < nlam; lam++)
  {
    // Evaluate the local constitutive matrix [Ql]. This matrix is computed
    // only when the layer material is different from the previous one.

    int idmat = (int)layup[2][lam];
 
    if (idmat != prevmat)
    {
      cMaterial *mat = cMaterial::GetMaterial(idmat);
      double *param = new double[mat->NumParam( )];
      mat->GetParam(param);
      QlMatrix(param, Ql);
      delete []param;
      prevmat = idmat;
    }

    // Evaluate the global constitutive matrix [Qg].

    double theta = layup[1][lam]*PI/180.0;  // deg to rad
    QgMatrix(theta, Ql, Qg);

    // Compute [A] and [D] matrices.

    zt = zb + layup[0][lam];
    double dz  = zt - zb;
    double d3z = zt*zt*zt - zb*zb*zb;
 
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        A[i][j] += Qg[i][j]*dz;
        D[i][j] += Qg[i][j]*d3z/3.0;
      }
    }

    // Update the bottom coordinate.

    zb = zt;
  }
}

// =============================== CalcABD =================================

void cLaminated :: CalcABD(cMatrix &layup, cMatrix &C)
{
  // Initialize the given matrix.

  C.Zero( );

  // Compute the laminate thickness.

  double lamthk = GetThickness(layup);

  // Evaluate the contribution of each ply.

  int nlam = layup.NCol( );
  int prevmat = -1;        // Previous material
  double zb = -lamthk/2.0; // Bottom coordinate (layer)
  double zt;               // Top coordinate (layer)
  cMatrix Ql(3,3);         // Local constitutive matrix
  cMatrix Qg(3,3);         // Global constitutive matrix
 
  for (int lam = 0; lam < nlam; lam++)
  {
    // Evaluate the local constitutive matrix [Ql]. This matrix is computed
    // only when the layer material is different from the previous one.

    int idmat = (int)layup[2][lam];
 
    if (idmat != prevmat)
    {
      cMaterial *mat = cMaterial::GetMaterial(idmat);
      double *param = new double[mat->NumParam( )];
      mat->GetParam(param);
      QlMatrix(param, Ql);
      delete []param;
      prevmat = idmat;
    }

    // Evaluate the global constitutive matrix [Qg].

    double theta = layup[1][lam]*PI/180.0;  // deg to rad
    QgMatrix(theta, Ql, Qg);

    // Compute the [C] matrix.

    zt = zb + layup[0][lam];
    double dz  = zt - zb;
    double d2z = zt*zt - zb*zb;
    double d3z = zt*zt*zt - zb*zb*zb;
 
    for (int i = 0; i < 3; i++)
    {
      int i3 = i + 3;
 
      for (int j = 0; j < 3; j++)
      {
        int j3 = j + 3;
        C[i ][j ] += Qg[i][j]*dz;          // [A]
        C[i ][j3] += Qg[i][j]*d2z/2.0;     // [B]
        C[i3][j ] += Qg[i][j]*d2z/2.0;     // [B]
        C[i3][j3] += Qg[i][j]*d3z/3.0;     // [D]
      }
    }

    // Update the bottom coordinate.

    zb = zt;
  }
}

// =============================== CalcABD with Failure Index List =================================

void cLaminated :: CalcABD(cMatrix &layup, cMatrix &C, cVector FL)
{
  // Initialize the given matrix.

  C.Zero( );

  // Compute the laminate thickness.

  double lamthk = GetThickness(layup);

  // Evaluate the contribution of each ply.

  int nlam = layup.NCol( );
  int FLi;						// Failure Index
  double zb = -lamthk/2.0; // Bottom coordinate (layer)
  double zt;               // Top coordinate (layer)
  cMatrix Ql(3,3);         // Local constitutive matrix
  cMatrix Qg(3,3);         // Global constitutive matrix
 
  for (int lam = 0; lam < nlam; lam++)
  {
    // Evaluate the local constitutive matrix [Ql].
      int idmat = (int)layup[2][lam];
      FLi = FL[lam];
      cMaterial *mat = cMaterial::GetMaterial(idmat);
      double *param = new double[mat->NumParam( )];
      mat->GetParam(param);
      QlMatrix(param, Ql, FLi);
      delete []param;

    // Evaluate the global constitutive matrix [Qg].

    double theta = layup[1][lam]*PI/180.0;  // deg to rad
    QgMatrix(theta, Ql, Qg);

    // Compute the [C] matrix.

    zt = zb + layup[0][lam];
    double dz  = zt - zb;
    double d2z = zt*zt - zb*zb;
    double d3z = zt*zt*zt - zb*zb*zb;
   
    for (int i = 0; i < 3; i++)
    {
      int i3 = i + 3;
   
      for (int j = 0; j < 3; j++)
      {
        int j3 = j + 3;
        C[i ][j ] += Qg[i][j]*dz;          // [A]
        C[i ][j3] += Qg[i][j]*d2z/2.0;     // [B]
        C[i3][j ] += Qg[i][j]*d2z/2.0;     // [B]
        C[i3][j3] += Qg[i][j]*d3z/3.0;     // [D]
      }
    }

    // Update the bottom coordinate.

    zb = zt;
  }
}


// ============================== MountABDG ================================

void cLaminated :: MountABDG(cMatrix &layup, cMatrix &C)
{
  // Initialize matrices.

  cMatrix lQb(3,3);  // Local bending constitutive matrix.
  cMatrix lQs(2,2);  // Local shear constitutive matrix.
  cMatrix gQb(3,3);  // Global bending constitutive matrix.
  cMatrix gQs(2,2);  // Global shear constitutive matrix.

  cMatrix A(3,3);    // A matrix.
  cMatrix B(3,3);    // B matrix.
  cMatrix D(3,3);    // D matrix.
  cMatrix G(2,2);    // G matrix.

  A.Zero( );
  B.Zero( );
  D.Zero( );
  G.Zero( );

  // Total laminate thickness and vertical coordinates.

  int nlam = layup.NCol( );
  cVector CoordZ(nlam+1);
  CoordZ.Zero( );

  double thksum = 0.0;
 
  for (int lam = 0; lam < nlam; lam++)
    thksum += layup[0][lam];

  CoordZ[0] = -thksum/2.0;

  for (int lam = 0; lam < nlam; lam++)
    CoordZ[lam+1] = CoordZ[lam] + layup[0][lam];

  // Evaluate the contribution of each ply.

  for (int lam = 0; lam < nlam; lam++)
  {
    // Clean existing matrices.

    lQb.Zero( );
    lQs.Zero( );
    gQb.Zero( );
    gQs.Zero( );

    // Define some constants.

    double theta = layup[1][lam]*PI/180.0;
    double c  = cos(theta);
    double s  = sin(theta);
    double c2 = c*c;
    double s2 = s*s;
    double c3 = c*c2;
    double s3 = s*s2;
    double c4 = c2*c2;
    double s4 = s2*s2;

    // Get the local constitutive matrices.

    cMaterial *mat = cMaterial::GetMaterial((int)layup[2][lam]);

    double *param = new double[mat->NumParam( )];
    mat->GetParam(param);

    QlMatrix(param, lQb, lQs);

    // Evaluate the global bending constitutive matrix (gQb)

    gQb[0][0] = c4*lQb[0][0] + 2*s2*c2*(lQb[0][1] + 2*lQb[2][2]) + s4*lQb[1][1];
    gQb[0][1] = s2*c2*(lQb[0][0] + lQb[1][1] - 4*lQb[2][2]) + (c4 + s4)*lQb[0][1];
    gQb[0][2] = s*c3*(lQb[0][0] - lQb[0][1] - 2*lQb[2][2]) + s3*c*(lQb[0][1] - lQb[1][1] + 2*lQb[2][2]);

    gQb[1][0] = s2*c2*(lQb[0][0] + lQb[1][1] - 4*lQb[2][2]) + (c4 + s4)*lQb[0][1];
    gQb[1][1] = s4*lQb[0][0] + 2*s2*c2*(lQb[0][1] + 2*lQb[2][2]) + c4*lQb[1][1];
    gQb[1][2] = s3*c*(lQb[0][0] - lQb[0][1] - 2*lQb[2][2]) + s*c3*(lQb[0][1] - lQb[1][1] + 2*lQb[2][2]);

    gQb[2][0] = s*c3*(lQb[0][0] - lQb[0][1] - 2*lQb[2][2]) + s3*c*(lQb[0][1] - lQb[1][1] + 2*lQb[2][2]);
    gQb[2][1] = s3*c*(lQb[0][0] - lQb[0][1] - 2*lQb[2][2]) + s*c3*(lQb[0][1] - lQb[1][1] + 2*lQb[2][2]);
    gQb[2][2] = s2*c2*(lQb[0][0] + lQb[1][1] - 2*lQb[0][1] - 2*lQb[2][2]) + (c4 + s4)*lQb[2][2];

    // Evaluate the global shear constitutive matrix (gQs)

    gQs[0][0] = c2*lQs[0][0] + s2*lQs[1][1];
    gQs[0][1] = s*c*(lQs[1][1] - lQs[0][0]);

    gQs[1][0] = s*c*(lQs[1][1] - lQs[0][0]);
    gQs[1][1] = s2*lQs[0][0] + c2*lQs[1][1];

    // Evaluate A, B, D, G matrices.

    double dz  = CoordZ[lam+1] - CoordZ[lam];
    double d2z = CoordZ[lam+1]*CoordZ[lam+1] - CoordZ[lam]*CoordZ[lam];
    double d3z = CoordZ[lam+1]*CoordZ[lam+1]*CoordZ[lam+1] - CoordZ[lam]*CoordZ[lam]*CoordZ[lam];

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
        A[i][j] += gQb[i][j]*dz;
        B[i][j] += gQb[i][j]*d2z/2.0;
        D[i][j] += gQb[i][j]*d3z/3.0;
      }

    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        G[i][j] += (5.0/6.0)*gQs[i][j]*dz;

    delete[] param;
  }

  // Mount the complete C Matrix.

  C.Zero( );

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      C[i][j] = A[i][j];
      C[3+i][j] = C[i][3+j] = B[i][j];
      C[3+i][3+j] = D[i][j];
    }

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      C[6+i][6+j] = G[i][j];
}

// ======================================================= End of file =====
