// -------------------------------------------------------------------------
// material.cpp - implementation of the Material class.
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
// Created:      13-Nov-2012     Iuri Barcelos Rocha
//               Based on cMaterial class of FAST system.
//
// Modified:     12-Apr-2014     Evandro Parente Junior
//               Implementation of 2D failure criteria.
//
//		 Jul-2014	 Joao Paulo Bernhardt
//		 Implementation of new failure criteria and progressive failure support.
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include "material.h"
#include "gblvar.h"
#include "input.h"
#include "utl.h"
#include "vec.h"

using namespace std;

// -------------------------------------------------------------------------
// Static variables:
//
int         cMaterial :: NumMat   = 0;
int         cMaterial :: NumIso   = 0;
int         cMaterial :: NumOrtho = 0;
cMaterial** cMaterial :: VecMat   = 0;
static double f12param            = -0.5;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadNumMat ===============================

void cMaterial :: ReadNumMat(istream &in)
{
  // Read the number of materials.

  if (!(in >> NumMat) || (NumMat <= 0))
  {
    cout << "Error in the input of the number of materials!" << endl;
    exit(0);
  }

  // Alloc the array of materials.

  VecMat = new cMaterial*[NumMat];
}

// ================================ ReadIso ================================

void cMaterial :: ReadIso(istream &in)
{
  // Read the number of elastic isotropic materials.

  int n;

  if (!(in >> n))
  {
    cout << "Error in the input of the number of isotropic materials!" << endl;
    exit(0);
  }

  NumIso += n;

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of elastic isotropic material " << i+1 << "!" << endl;
      exit(0);
    }
    cMatIso *mat = new cMatIso(id);
    mat->Read( );
  }
}

// ================================ ReadOrtho ==============================

void cMaterial :: ReadOrtho(istream &in)
{
  // Read the number of elastic orthotropic materials.

  int n;

  if (!(in >> n))
  {
    cout << "Error in the input of the number of orthotropic materials!" << endl;
    exit(0);
  }

  NumOrtho += n;

  // Read each material.

  int id;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> id) || id < 1 || id > NumMat)
    {
      cout << "Error in the input of elastic orthotropic material " << i+1 << "!" << endl;
      exit(0);
    }
    cMatOrtho *mat = new cMatOrtho(id);
    mat->Read( );
  }
}

// ================================ Destroy ================================

void cMaterial :: Destroy(void)
{
  // Destroy each material.

  for (int i = 0; i < NumMat; i++) delete VecMat[i];

  // Release the array of materials.

  delete []VecMat;
}

// =============================== ReadDensity =============================

void cMaterial :: ReadDensity(istream &in)
{
  int n;

  if (!(in >> n))
  {
    cout << "Error on reading number of densities!" << endl;
    exit(0);
  }

  int label;
  double dens;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> label) || !(in >> dens))
    {
      cout << "Error on reading density " << i+1 << "!" << endl;
      exit(0);
    }

    cMaterial *pcMat = GetMaterial(label);
    if (!pcMat)
    {
      cout << "Error on reading density " << i+1 << "!" << endl;
      exit(0);
    }
    pcMat->Density = dens;
 }
}

// =============================== ReadCost ================================

void cMaterial :: ReadCost(istream &in)
{
  int n;

  if (!(in >> n))
  {
    cout << "Error on reading number of costs!" << endl;
    exit(0);
  }

  int label;
  double cost;
  for (int i = 0; i < n; i++)
  {
    if (!(in >> label) || !(in >> cost))
    {
      cout << "Error on reading cost " << i+1 << "!" << endl;
      exit(0);
    }

    cMaterial *pcMat = GetMaterial(label);
    if (!pcMat)
    {
      cout << "Error on reading cost " << i+1 << " (material " << i+1;
      cout << " do not exist!)." << endl;
      exit(0);
    }
    pcMat->Cost = cost;
 }
}

// ============================== GetMaterial ==============================

cMaterial *cMaterial :: GetMaterial(int id)
{
  if (id > 0 && id <= NumMat) return(VecMat[id-1]);
  return(0);
}

// ============================== cMaterial ================================

cMaterial :: cMaterial(int id)
{
  Density = 1.0;
  Label = id;           // Store the material label
  VecMat[id-1] = this;  // Add the material to the material vector
}

// ============================== ~cMaterial ===============================

cMaterial :: ~cMaterial(void)
{
}

// -------------------------------------------------------------------------
// Class cMatIso:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Static variables:
//
eFailType cMatIso :: FailType = VON_MISES;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadFailType ==============================

void cMatIso :: ReadFailType(istream &in)
{
  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of isotropic failure criteria." << endl;
    exit(0);
  }

  if (string(label) == "VonMises" || string(label) == "vonmises")
    FailType = VON_MISES;
  else
  {
    cout << "Unknown isotropic failure criterium: " << label << endl;
    exit(0);
  }
}

// ================================ cMatIso ================================

cMatIso :: cMatIso(int id) : cMaterial(id)
{
  Type = MAT_ISOTROPIC;
  E  = 1.0;
  Nu = 0.0;
  SMYS = 0.0;
  SMTS = 0.0;
}

// ============================== ~cMatIso =================================

cMatIso :: ~cMatIso(void)
{
}

// ================================= Read ==================================

void cMatIso :: Read(void)
{
  if (!(in >> E) || !(in >> Nu) || !(in >> SMYS) || !(in >> SMTS))
  {
    cout << "Error in the input of material " << Label << "!" << endl;
    exit(0);
  }
}

// ============================== SafetyFactor =============================

double cMatIso :: SafetyFactor(cVector &stress)
{
  double sf = 0.0;  // Safety factor

  if (FailType == VON_MISES)
  {
    if (stress.Dim( ) == 3)
      sf = VonMises2D(stress);
    else
      sf = VonMises3D(stress);
  }
  else
  {
    cout << "Invalid failure criterion for isotropic materials." << endl;
    exit(0);
  }

  return(sf);
}

double cMatIso :: SafetyFactor(cVector &stress, int &FLind)
{
  double sf = 0.0;  // Safety factor

  if (FailType == VON_MISES)
  {
    if (stress.Dim( ) == 3)
      sf = VonMises2D(stress);
    else
      sf = VonMises3D(stress);
  }
  else
  {
    cout << "Invalid failure criterion for isotropic materials." << endl;
    exit(0);
  }

  return(sf);
}

// =============================== GetParam ================================

void cMatIso :: GetParam(double *param)
{
  param[0] = E;
  param[1] = Nu;
  param[2] = SMYS;
  param[3] = SMTS;
}

// ============================= VonMises2D ================================

double cMatIso :: VonMises2D(cVector &stress)
{
  double aux1 = stress[0]*stress[0];
  double aux2 = stress[0]*stress[1];
  double aux3 = stress[1]*stress[1];
  double aux4 = stress[2]*stress[2];

  double SigVM = sqrt(aux1 - aux2 + aux3 + 3.0*aux4);
  double sf = SMYS/SigVM;

 // cout << "SF " << sf << endl;

  return(sf);
}

// ============================= VonMises3D ================================

double cMatIso :: VonMises3D(cVector &stress)
{
  double aux1 = (stress[0] - stress[1])*(stress[0] - stress[1]);
  double aux2 = (stress[0] - stress[2])*(stress[0] - stress[2]);
  double aux3 = (stress[1] - stress[2])*(stress[1] - stress[2]);
  double aux4 = stress[3]*stress[3];
  double aux5 = stress[4]*stress[4];
  double aux6 = stress[5]*stress[5];

  double SigVM = sqrt(0.5*(aux1 + aux2 + aux3) + 3.0*(aux4 + aux5 + aux6));
  double sf = SMYS/SigVM;


  return(sf);
}

// -------------------------------------------------------------------------
// Class cMatOrtho:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Static variables:
//
eFailType cMatOrtho :: FailType = TSAI_WU;

// -------------------------------------------------------------------------
// Public methods:
//

// ========================= TsaiWuF12Param ================================

void cMatOrtho :: TsaiWuF12Param(double param)
{ 
  f12param = param; 
}

// ============================== ReadFailType =============================

void cMatOrtho :: ReadFailType(istream &in)
{
  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of orthotropic failure criteria." << endl;
    exit(0);
  }

  if (string(label) == "MaxStress" || string(label) == "maxstress" ||
      string(label) == "MaximumStress" || string(label) == "maximumstress")
    FailType = MAX_STRESS;
  else if (string(label) == "MaxStrain" || string(label) == "maxstrain" ||
      string(label) == "MaximumStrain" || string(label) == "maximumstrain")
    FailType = MAX_STRAIN;
  else if (string(label) == "Tsai-Hill" || string(label) == "tsai-hill" ||
           string(label) == "TsaiHill" || string(label) == "tsaihill")
    FailType = TSAI_HILL;
  else if (string(label) == "Tsai-Wu" || string(label) == "tsai-wu" ||
           string(label) == "TsaiWu" || string(label) == "tsaiwu")
    FailType = TSAI_WU;
  else if (string(label) == "Tsai-Wu3D" || string(label) == "tsai-wu3d" ||
           string(label) == "TsaiWu3D" || string(label) == "tsaiwu3d")
    FailType = TSAI_WU;
  else
  {
    cout << "Unknown orthotropic failure criterium: " << label << endl;
    exit(0);
  }
}

// =============================== cMatOrtho ===============================

cMatOrtho :: cMatOrtho(int id) : cMaterial(id)
{
  Type = MAT_ORTHOTROPIC;
  E1   = 1.0;
  E2   = 1.0;
  E3   = 1.0;
  Nu12 = 0.0;
  Nu13 = 0.0;
  Nu23 = 0.0;
  G12  = 0.0;
  G13  = 0.0;
  G23  = 0.0;
  Xt   = 0.0;
  Xc   = 0.0;
  Yt   = 0.0;
  Yc   = 0.0;
  Zt   = 0.0;
  Zc   = 0.0;
  S6   = 0.0;
  S4   = 0.0;
  S5   = 0.0;
}

// ============================== ~cMatOrtho ===============================

cMatOrtho :: ~cMatOrtho(void)
{
}

// ================================== Read =================================

void cMatOrtho :: Read(void)
{
  if (!(in >> E1)   || !(in >> E2)   || !(in >> E3)  || !(in >> Nu12) ||
      !(in >> Nu13) || !(in >> Nu23) || !(in >> G12) || !(in >> G13)  ||
      !(in >> G23)  || !(in >> Xt)   || !(in >> Xc)  || !(in >> Yt)   ||
      !(in >> Yc)   || !(in >> Zt)   || !(in >> Zc)  || !(in >> S6)   ||
      !(in >> S4)   || !(in >> S5))
  {
    cout << "Error in the input of material " << Label << "!" << endl;
    exit(0);
  }
}

// ================================ GetParam ===============================

void cMatOrtho :: GetParam(double *param)
{
  param[0] = E1;
  param[1] = E2;
  param[2] = E3;
  param[3] = Nu12;
  param[4] = Nu13;
  param[5] = Nu23;
  param[6] = G12;
  param[7] = G13;
  param[8] = G23;
}

// ============================== SafetyFactor =============================

double cMatOrtho :: SafetyFactor(cVector &stress)

// For MAX_STRAIN, the &stress vector should actually receive strain (eps1, eps2..) parameters.

{
  double sf = 0.0;  // Safety factor

  if (FailType == MAX_STRESS)
  {
    if (stress.Dim( ) == 3)
      sf = MaxStress2D(stress);
    else
      sf = MaxStress3D(stress);
  }
  else if (FailType == MAX_STRAIN)
  {
    if (stress.Dim( ) == 3)
      sf = MaxStrain2D(stress);
    else
      sf = MaxStrain3D(stress);
  }
  else if (FailType == TSAI_HILL)
  {
    if (stress.Dim( ) == 3)
      sf = TsaiHill(stress);
    else
    {
      cout << "Invalid failure criterium Tsai-Hill 3D!\n";
      exit(0);
    }
  }
  else if (FailType == TSAI_WU)
  {
    if (stress.Dim( ) == 3)
      sf = TsaiWu2D(stress);
    else
      sf = TsaiWu3D(stress);
  }
  else
  {
    cout << "Invalid failure criterion for orthotropic materials." << endl;
    exit(0);
  }

  return(sf);
}

// ============================== SafetyFactor Progressive =============================

// cVector &stress actually is [sig1,sig2,sig3,eps1,eps2,eps3]
// int &FLind identifies the failure mode.

double cMatOrtho :: SafetyFactor(cVector &stress, int &FLind)
{
  double sf = 0.0;  // Safety factor

  if (FailType == MAX_STRESS)
  {
    if (stress.Dim( ) == 6)
      sf = MaxStress2D(stress, FLind);
    else
      sf = MaxStress3D(stress);
  }
  else if (FailType == MAX_STRAIN)
  {
    if (stress.Dim( ) == 6)
      sf = MaxStrain2D(stress, FLind);
    else
      sf = MaxStrain3D(stress);
  }
  else if (FailType == TSAI_HILL)
  {
    if (stress.Dim( ) == 3)
      sf = TsaiHill(stress);
    else
    {
      cout << "Invalid failure criterium Tsai-Hill 3D!\n";
      exit(0);
    }
  }
  else if (FailType == TSAI_WU)
  {
    if (stress.Dim( ) == 6)
      sf = TsaiWu2D(stress, FLind);
    else
      sf = TsaiWu3D(stress);
  }
  else
  {
    cout << "Invalid failure criterion for orthotropic materials." << endl;
    exit(0);
  }

  return(sf);
}

// =============================== MaxStrain2D =============================

double cMatOrtho :: MaxStrain2D(cVector &strain)
{
  double TOL = 1.0e-12;
  double sf  = 1.0/TOL;
  double aux;

  // Normal strain.

  if (fabs(strain[3]) > TOL)
  {
    if (strain[3] > 0.0)
      aux =  (Xt/E1)/strain[3];
    else
      aux = (-Xc/E1)/strain[3];
    sf = Utl::Min(sf, aux);
  }

  if (fabs(strain[4]) > TOL)
  {
    if (strain[4] > 0.0)
      aux =  (Yt/E2)/strain[4];
    else
      aux = (-Yc/E2)/strain[4];
    sf = Utl::Min(sf, aux);
  }

  // Shear strain.

  if (fabs(strain[5]) > TOL)
  {
    aux = (S6/G12)/fabs(strain[5]);
    sf  = Utl::Min(sf, aux);
  }

  return(sf);
}

// =============================== MaxStrain2D Progressive =============================

double cMatOrtho :: MaxStrain2D(cVector &strain, int &FLind)
{
  double TOL = 1.0e-12;
  double sf  = 1.0/TOL;
  double aux;

  // Normal strain.

  if (fabs(strain[3]) > TOL)
  {
    if (strain[3] > 0.0)
      aux =  (Xt/E1)/strain[3];
    else
      aux = (-Xc/E1)/strain[3];
    sf = Utl::Min(sf, aux);
    FLind = 1;
  }

  if (fabs(strain[4]) > TOL)
  {
    if (strain[4] > 0.0)
      aux =  (Yt/E2)/strain[4];
    else
      aux = (-Yc/E2)/strain[4];
    sf = Utl::Min(sf, aux);
    FLind = 2;
  }

  // Shear strain.

  if (fabs(strain[5]) > TOL)
  {
    aux = (S6/G12)/fabs(strain[5]);
    sf  = Utl::Min(sf, aux);
    FLind = 3;
  }
  return(sf);
}

// =============================== MaxStrain3D =============================

double cMatOrtho :: MaxStrain3D(cVector &strain)
{
  return(0); 
}

// =============================== MaxStress2D =============================

double cMatOrtho :: MaxStress2D(cVector &stress)
{
  double TOL = 1.0e-12;
  double sf  = 1.0/TOL;
  double aux;

  // Normal stresses.

  if (fabs(stress[0]) > TOL)
  {
    if (stress[0] > 0.0)
      aux =  Xt/stress[0];
    else
      aux = -Xc/stress[0];
    sf = Utl::Min(sf, aux);
  }

  if (fabs(stress[1]) > TOL)
  {
    if (stress[1] > 0.0)
      aux =  Yt/stress[1];
    else
      aux = -Yc/stress[1];
    sf = Utl::Min(sf, aux);
  }

  // Shear stress.

  if (fabs(stress[2]) > TOL)
  {
    aux = S6/fabs(stress[2]);
    sf  = Utl::Min(sf, aux);
  }

  return(sf);
}

// =============================== MaxStress2D Progressive =============================

double cMatOrtho :: MaxStress2D(cVector &stress, int &FLind)
{
  double TOL = 1.0e-12;
  double sf  = 1.0/TOL;
  double aux;

  // Normal stresses.

  if (fabs(stress[0]) > TOL)
  {
    if (stress[0] > 0.0)
      aux =  Xt/stress[0];
    else
      aux = -Xc/stress[0];
    sf = Utl::Min(sf, aux);
    FLind = 1;
  }

  if (fabs(stress[1]) > TOL)
  {
    if (stress[1] > 0.0)
      aux =  Yt/stress[1];
    else
      aux = -Yc/stress[1];
    sf = Utl::Min(sf, aux);
    FLind = 2;
  }

  // Shear stress.

  if (fabs(stress[2]) > TOL)
  {
    aux = S6/fabs(stress[2]);
    sf  = Utl::Min(sf, aux);
    FLind = 3;
  }

  return(sf);
}

// =============================== MaxStress3D =============================

double cMatOrtho :: MaxStress3D(cVector &stress)
{
  double TOL = 1.0e-12;
  double sf  = 1.0/TOL;
  double aux;

  // Normal stresses.

  if (fabs(stress[0]) > TOL)
  {
    if (stress[0] > 0.0)
      aux =  Xt/stress[0];
    else
      aux = -Xc/stress[0];
    sf = Utl::Min(sf, aux);
  }

  if (fabs(stress[1]) > TOL)
  {
    if (stress[1] > 0.0)
      aux =  Yt/stress[1];
    else
      aux = -Yc/stress[1];
    sf = Utl::Min(sf, aux);
  }

  if (fabs(stress[2]) > TOL)
  {
    if (stress[2] > 0.0)
      aux =  Zt/stress[2];
    else
      aux = -Zc/stress[2];
    sf = Utl::Min(sf, aux);
  }

  // Shear stresses.

  if (fabs(stress[5]) > TOL)
  {
    aux = S4/fabs(stress[5]);
    sf  = Utl::Min(sf, aux);
  }

  if (fabs(stress[4]) > TOL)
  {
    aux = S5/fabs(stress[4]);
    sf  = Utl::Min(sf, aux);
  }

  if (fabs(stress[3]) > TOL)
  {
    aux = S6/fabs(stress[3]);
    sf  = Utl::Min(sf, aux);
  }

  return(sf);
}

// ================================ TsaiHill ===============================

double cMatOrtho :: TsaiHill(cVector &stress)
{
  double F1;
  double A;
    if (stress[0] > 0.0)
      F1 = Xt;
    else
      F1 = Xc;
    A = 1.0/(F1*F1);

  double F2;
  double B;
    if (stress[1] > 0.0)
      F2 = Yt;
    else
      F2 = Yc;
    B = 1.0/(F2*F2);

  double D;
  D  = 1.0/(S6*S6);

  double C = -1.0*A;

  double aux1 = stress[0]*stress[0];
  double aux2 = stress[1]*stress[1];
  double aux3 = stress[0]*stress[1];
  double aux4 = stress[3]*stress[3];

  double fi = A*aux1 + B*aux2 + C*aux3 + D*aux4;

  return(fi);
}

// ================================ TsaiWu2D ===============================

double cMatOrtho :: TsaiWu2D(cVector &stress)
{
 /* cout << Xt << Xc << Yt << Yc << S6 << endl;
  cout << stress[0] << stress[1] << stress[2] << endl;*/

  double f1  = 1.0/Xt - 1.0/Xc;
  double f2  = 1.0/Yt - 1.0/Yc;
  double f11 = 1.0/(Xt*Xc);
  double f22 = 1.0/(Yt*Yc);
  double f12 = f12param * sqrt(f11*f22);
  double f66 = 1.0/(S6*S6);

  double aux1 = stress[0]*stress[0];
  double aux2 = stress[1]*stress[1];
  double aux3 = stress[0]*stress[1];
  double aux6 = stress[2]*stress[2];

  double a  = f11*aux1 + f22*aux2 + 2.0*f12*aux3 + f66*aux6;
  double b  = f1*stress[0] + f2*stress[1];
  double sf = (-b + sqrt(b*b + 4.0*a))/(2.0*a);

  return(sf);
}

// ================================ TsaiWu2D Progressive ===============================

double cMatOrtho :: TsaiWu2D(cVector &stress, int &FLind)
{
  double f1  = 1.0/Xt - 1.0/Xc;
  double f2  = 1.0/Yt - 1.0/Yc;
  double f11 = 1.0/(Xt*Xc);
  double f22 = 1.0/(Yt*Yc);
  double f12 = f12param * sqrt(f11*f22);
  double f66 = 1.0/(S6*S6);

  double aux1 = stress[0]*stress[0];
  double aux2 = stress[1]*stress[1];
  double aux3 = stress[0]*stress[1];
  double aux6 = stress[2]*stress[2];

  double a  = f11*aux1 + f22*aux2 + 2.0*f12*aux3 + f66*aux6;
  double b  = f1*stress[0] + f2*stress[1];
  double sf = (-b + sqrt(b*b + 4.0*a))/(2.0*a);

  double H1 = fabs(f1*stress[0] + f11*aux1);
  double H2 = fabs(f2*stress[1] + f22*aux2);
  double H6 = fabs(f66*aux6);

  cVector H(3);
  H[0] = H1;
  H[1] = H2;
  H[2] = H6;
  if (H.Max() == H1) FLind = 1;
  else if (H.Max() == H2) FLind = 2;
  else if (H.Max() == H6) FLind = 3;

  return(sf);
}

// =============================== TsaiWu3D ================================

double cMatOrtho :: TsaiWu3D(cVector &stress)
{
  double f1  = 1.0/Xt - 1.0/Xc;
  double f2  = 1.0/Yt - 1.0/Yc;
  double f3  = 1.0/Zt - 1.0/Zc;
  double f11 = 1.0/(Xt*Xc);
  double f22 = 1.0/(Yt*Yc);
  double f33 = 1.0/(Zt*Zc);
  double f44 = 1.0/(S5*S5);
  double f55 = 1.0/(S4*S4);
  double f66 = 1.0/(S6*S6);
  double f12 = -0.5*sqrt(f11*f22);
  double f13 = -0.5*sqrt(f11*f33);
  double f23 = -0.5*sqrt(f22*f33);

  double aux1 = stress[0]*stress[0];
  double aux2 = stress[1]*stress[1];
  double aux3 = stress[2]*stress[2];
  double aux4 = stress[4]*stress[4];
  double aux5 = stress[5]*stress[5];
  double aux6 = stress[3]*stress[3];

  double aux12 = stress[0]*stress[1]; //aux3
  double aux13 = stress[0]*stress[2];
  double aux23 = stress[1]*stress[2];

  double a  = f11*aux1 + f22*aux2 + f33*aux3 + f44*aux4 + f55*aux5 + f66*aux6 + 2.0*f12*aux12 + 2.0*f13*aux13 + 2.0*f23*aux23;
  double b  = f1*stress[0] + f2*stress[1] + f3*stress[2];
  double sf = (-b + sqrt(b*b + 4.0*a))/(2.0*a);

  return(sf);
}

// ======================================================= End of file =====
