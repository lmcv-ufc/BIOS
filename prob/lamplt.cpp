// -------------------------------------------------------------------------
// lamplt.cpp - Implementation of the Laminated Plate problem class.
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
// Modified:     Jul-2014    Joao Paulo Bernhardt
//               Added Progressive Failure problem classes
//
// Modified:	 04-Feb-2015	Elias Saraiva Barroso
//               Creation of cLamPlate class.
//
// Modified:	 28-Jun-2016	Elias Saraiva Barroso
//		         Creation of the multi-objective cLamPltBuckMOBJ class.
//
// Modified:     04-Dec-2017    Marina Alves Maia
//               Creation of the multiobjective cLamPltFrequency class.
//
// Modified:     30-Jul-2018	Evandro Parente and David Sena
//               Implementation of Surrogate Based Optimization.
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
#include "lam.h"
#include "lamplt.h"
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
// Static variables:
//
int            cLamPlate :: MaxSurGen = -1;
vector<double> cLamPlate :: Loads;
vector<double> cLamPlate :: Dims;

// -------------------------------------------------------------------------
// Register problems on the problem factory:
//
static const bool registeredProb[] =
{
  cProblemFactory :: Register("LamPltBuckMOBJ"         , MakeProb<cLamPltBuckMOBJ>              ,".lam"),
  cProblemFactory :: Register("LamPltLiu2000"          , MakeProb<cLamPltLiu2000>               ,".lam"),
  cProblemFactory :: Register("LamPltLoadFactor"       , MakeProb<cLamPltLoadFactor>            ,".lam"),
  cProblemFactory :: Register("LamPltMinLopez2009"     , MakeProb<cLamPltMinLopez2009>          ,".lam"),
  cProblemFactory :: Register("LamPltMFBuck"           , MakeProb<cSquarePlateMFBuckLam>        ,".lam"),
  cProblemFactory :: Register("LamPltBuckKeshtegarMF"  , MakeProb<cSquarePlateKeshtegarBuckLam> ,".lam")
};

// -------------------------------------------------------------------------
// Public methods:
//

// =========================== ReadLamMaxSurGen ============================

void cLamPlate :: ReadLamMaxSurGen(istream &in)
{
  if (!(in >> MaxSurGen))
    Utl::Exit("Error in the input of the LamPlate MaxSurGen.");
  else
    cout << "MaxSurGen = " << MaxSurGen << endl;
}

// ============================ ReadLamPltLoads ============================

void cLamPlate :: ReadLamPltLoads(istream &in)
{
  double load;

  for(int i = 0; i < 6; i++)
  {
    if (!(in >> load))
      Utl::Exit("Error in the input of the LamPlate Loads.");
    else
      Loads.push_back(load);
  }
}

// =============================== ReadLamPltDims ==========================

void cLamPlate :: ReadLamPltDims(istream &in)
{
  double dim;

  for (int i = 0; i < 2; i++)
  {
    if (!(in >> dim))
      Utl::Exit("Error in the input of the LamPlate dimensions.");
    else
      Dims.push_back(dim);
  }
}

// ============================= cLamPlate =================================

cLamPlate :: cLamPlate(void)
{
}

// ============================= ~cLamPlate ================================

cLamPlate :: ~cLamPlate(void)
{
}

// ============================== LoadReadFunc =============================

void cLamPlate :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cLaminated :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("LAMPLT.MAX.SUR.GENERATIONS",makeReadObj(cLamPlate,ReadLamMaxSurGen));
  im.Insert("LAMPLT.LOADS",makeReadObj(cLamPlate,ReadLamPltLoads));
  im.Insert("LAMPLT.DIMENSIONS",makeReadObj(cLamPlate,ReadLamPltDims));
}

// ======================== FindNewSampPoint  ==============================

void cLamPlate :: FindNewSampPoint(cVector &bestPosition, int SwarmSize, int MaxGen)
{
 // cout << "FindNewSampPoint\n";
 
  // PSO
 
  double best = 0;
  double w = 0.7;
  double c1 = 2.5;
  double c2 = 2.5;
  double tempo = 0;
  cVector fobj(SwarmSize);
  cVector fobjp(SwarmSize);
  cMatrix xx(SwarmSize,NumVar);
  cMatrix xn(SwarmSize,NumVar);
  cMatrix v(SwarmSize,NumVar);
  cMatrix xp(SwarmSize,NumVar);
  cVector xg(NumVar);
  double bestParticle = 0;
  bestPosition(NumVar);
  int inferior = 0;
  int superior = 1;
 
  for (int i = 0; i < SwarmSize; i++)
  {
    for (int j = 0; j < NumVar; j++)
    {
      // Initialize the particle position
 
      xx[i][j] = Utl::RandDouble(inferior, superior);
 
      // Initialize the particle velocity
         
      v[i][j] = (superior - inferior) * Utl::RandDouble(-1,1);
 
      // Assign the particle position to the best position so far
 
      xp[i][j] = xx[i][j];
    }
 
    cVector temp(NumVar);
    for (int k = 0; k < NumVar; k++)
      temp[k] = xx[i][k];
 
    fobj[i] = EvalSampPoint(temp);        
    fobjp[i] = fobj[i];
 
    // Find the best neighbor
 
    for (int k = 0; k < SwarmSize; k++)
    {
      if (k == 0)
      {
        best = fobjp[k];
        for (int q = 0; q < NumVar; q++)
          xg[q] = xx[k][q];
      }
      else if (k > 0)
      {
        if (fobjp[k] < best)
        {
          best = fobjp[k];
          for (int q = 0; q < NumVar; q++)
            xg[q] = xx[k][q];
        }
      }
    }
  }
 
  // Run the generations
 
  for (int gen = 0; gen < MaxGen; gen++)
  {
    // End find the best neighbor
 
    for (int i = 0; i < SwarmSize; i++)
    {
      // Calculate the rates r1 and r2 randomly
 
      double r1 = Utl::RandDouble(0,1);
      double r2 = Utl::RandDouble(0,1);
 
      // Calculate the velocity of the particle
 
      for (int r = 0; r < NumVar; r++)
        v[i][r] = w*v[i][r]+c1*r1*(xp[i][r]-xx[i][r])+c2*r2*(xg[r]-xx[i][r]);
 
      // Actualize the position of the particle
 
      for (int s = 0; s < NumVar; s++)
      {
        xx[i][s] = xx[i][s] + v[i][s];
        if(xx[i][s] < inferior)
        {
          xx[i][s] = inferior;
          v[i][s] = -0.5 * v[i][s];
        }
        if (xx[i][s] > superior)
        {
          xx[i][s] = superior;
          v[i][s] = -0.5 * v[i][s];
        }
      }
 
      // Temporary variable of the objective function
 
      cVector temporary(NumVar);
      for (int k = 0; k < NumVar; k++)
        temporary[k] = xx[i][k];
 
      tempo = EvalSampPoint(temporary);        
                           
      // Actualize or not actualize the best position of the particle
 
      if (tempo < fobjp[i])
      {
        fobj[i] = tempo;
        fobjp[i] = tempo;
 
        // Actualize the best position of the particle
 
        for (int p = 0; p < NumVar; p++)
          xp[i][p] = xx[i][p];
      }
      else
      {
        // Do not actualize the best position of the particle
 
        fobj[i] = tempo;
      }
      if (fobjp[i] < best)
      {
        best = fobjp[i];
        for (int q = 0; q < NumVar; q++)
          xg[q] = xp[i][q];
      }    
    }
  }

  // Find the best particle position and objective function
 
  for (int i = 0; i < SwarmSize; i++)
  {
    if (i == 0)
    {
      bestParticle = fobjp[i];
      for (int j = 0; j < NumVar; j++)
        bestPosition[j] = xp[i][j];
    }
    else if (i > 0)
    {
      if (fobjp[i] < bestParticle)
      {
        bestParticle = fobjp[i];
        for (int j = 0; j < NumVar; j++)
          bestPosition[j] = xp[i][j];
      }
    }
  }
}

// ============================ CalcBuckLoad  ==============================

double cLamPlate :: CalcBuckLoad(cMatrix &layup, double a, double b, 
		                 double Nx, double Ny)
{
  // Get stiffness data.

  cMatrix D(3,3);
  CalcD(layup, D);
  double D11 = D[0][0];
  double D12 = D[0][1];
  double D22 = D[1][1];
  double D66 = D[2][2];

  // Evaluate the minimum buckling load.

  int pmax = 21;
  int qmax = 21;
  double minlbd = 0.0;
  for (int p = 1; p < pmax; p++)
  {
    for (int q = 1; q < qmax; q++)
    {
      double pa2 = (p/a)*(p/a);
      double pa4 = pa2*pa2;
      double qb2 = (q/b)*(q/b);
      double qb4 = qb2*qb2;
      double lbd = (PI*PI*(D11*pa4 + 2.0*(D12 + 2.0*D66)*pa2*qb2 + D22*qb4))/(pa2*Nx + qb2*Ny);

      if (p == 1 && q == 1)
        minlbd = lbd;
      else if (lbd < minlbd)
        minlbd = lbd;
    }
  }

  return(minlbd);
}

// -------------------------------------------------------------------------
// Public methods:
//

// -------------------------------------------------------------------------
// Class cLamPltBuckMOBJ:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cLamPltBuck ===============================

cLamPltBuckMOBJ :: cLamPltBuckMOBJ(void)
{
  NumObj = 1;                  // Number of objective functions
  NumConstr = 1;               // Number of constraints
  NumConstr += 2;              // Display constraints
}

// ============================== Evaluate =================================

void cLamPltBuckMOBJ :: Evaluate(int** algvar, cVector &c, cVector &fobjs)
{
  // Problem data.

  double limlbd = 100.0; // Threshold buckling factor

  // Default load and dimensions values

  double Nx = 175.0;   // Compressive load in the x-direction
  double Ny = 175.0;   // Compressive load in the y-direction

  double a  = 0.9144; //0.92; // Plate length
  double b  = 0.7620; //0.75; // Plate width

  if (!Loads.empty( ))
  {
    Nx  = Loads[0];
    Ny  = Loads[1];
  }

  if (!Dims.empty( ))
  {
    a  = Dims[0];
    b  = Dims[1];
  }

  // Decode the variables.

  cMatrix layup;
  if (!Decode(algvar, layup))
  {
    cout << "Decoding process failure." << endl;
    fobjs[0] = 10e6;
    //fobjs[1] = 10e6;
  }

  // Get stiffness data.

  cMatrix D(3,3);
  CalcD(layup, D);
  double D11 = D[0][0];
  double D12 = D[0][1];
  double D22 = D[1][1];
  double D66 = D[2][2];

  // Evaluate the minimum buckling load.

  int pmax = 21;
  int qmax = 21;
  double minlbd = 0.0;
  for (int p = 1; p < pmax; p++)
  {
    for (int q = 1; q < qmax; q++)
    {
      double pa2 = (p/a)*(p/a);
      double pa4 = pa2*pa2;
      double qb2 = (q/b)*(q/b);
      double qb4 = qb2*qb2;
      double lbd = (PI*PI*(D11*pa4 + 2.0*(D12 + 2.0*D66)*pa2*qb2 + D22*qb4))/(pa2*Nx + qb2*Ny);

      if (p == 1 && q == 1)
        minlbd = lbd;
      else if (lbd < minlbd)
        minlbd = lbd;
    }
  }

 // cout << "Carga de flambagem: " << minlbd << endl;

  // Buckling constraint.

  c[0] = limlbd/minlbd-1;

  // Plate weight and cost.

  int nlam = layup.NCol( );

 // cout << "Numero de laminas " << nlam << endl;

  int lammat;
  double W = 0.0;
  double C = 0.0;
  double dens = 0.0;
  double relcost = 0.0;
  double area = a*b;
  for (int lam = 0; lam < nlam; lam++)
  {
    lammat = (int)layup[2][lam];
    dens = cMaterial::GetMaterial(lammat)->GetDensity( );
    relcost = cMaterial::GetMaterial(lammat)->GetCost( );
    W += 9.807*area*dens*layup[0][lam];
    C += area*dens*layup[0][lam]*relcost;
  }

  // Stores the value for each objective function

  double w = 1.0;
  double m = 2.0;

  double Wmin = 55.729;
  double Wmax = 89.935;
  double Cmin = 9.171;
  double Cmax = 45.461;

  fobjs[0] = pow(w * (W-Wmin)/(Wmax-Wmin),m) + pow( (1.0-w)*(C-Cmin)/(Cmax-Cmin),m);

  // Display constraints.
  c[1] = -W;
  c[2] = -C;
}

// -------------------------------------------------------------------------
// Class cLamPltLiu2000:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cLamPltLiu2000 ===========================

cLamPltLiu2000 :: cLamPltLiu2000(void)
{
  NumConstr = 6;
  NumObj = 1;
}
// ============================== Evaluate =================================

void cLamPltLiu2000 :: Evaluate(int **algvar, cVector &c, cVector &fobjs)
{
  // Problem data.

  int maxcont = 4;        // Max. num. of contiguous layers with same angle

  // Default load and dimensions values.

  double Nx  = 20000.0;   // Axial load in the x-direction
  double Ny  = 2000.0;    // Axial load in the y-direction
  double Nxy = 1000;      // Shear load
  double Mx = 0;
  double My = 0;
  double Mxy = 0;

  double a = 24;          // Plate length
  double b = 24;          // Plate width

  if (!Loads.empty( ))
  {
    Nx  = Loads[0];
    Ny  = Loads[1];
    Nxy = Loads[2];
  }

  if (!Dims.empty( ))
  {
    a  = Dims[0];
    b  = Dims[1];
  }

  // Decode the variables.

  cMatrix layup;
  if (!Decode(algvar, layup))
  {
    cout << "Decoding process failure." << endl;
    fobjs[0] = 1.0e6;
  }

  // Compute the ABD matrix.

  cMatrix C(6,6);
  CalcABD(layup, C);

  // Get stiffness data.

  double D11 = C[3][3];
  double D12 = C[3][4];
  double D22 = C[4][4];
  double D66 = C[5][5];

  // Evaluate the minimum buckling load.

  int pmax = 21;
  int qmax = 21;
  double minlbd = 0.0;
  for (int p = 1; p < pmax; p++)
  {
    for (int q = 1; q < qmax; q++)
    {
      double pa2 = (p/a)*(p/a);
      double pa4 = pa2*pa2;
      double qb2 = (q/b)*(q/b);
      double qb4 = qb2*qb2;
      double lbd = (PI*PI*(D11*pa4 + 2.0*(D12 + 2.0*D66)*pa2*qb2 + D22*qb4))/(pa2*Nx + qb2*Ny);

      if (p == 1 && q == 1)
        minlbd = lbd;
      else if (lbd < minlbd)
        minlbd = lbd;
    }
  }

  static double r[][2] = {{0.00,11.71},
                          {0.20,11.80},
                          {0.50,12.20},
                          {1.00,13.17},
                          {2.00,10.80},
                          {3.00,9.950},
                          {5.00,9.250},
                          {10.0,8.700},
                          {20.0,8.400},
                          {40.0,8.250},
                          {10e6,8.130}};

  double B1=0;
  double R = sqrt(D11*D22)/(D12+2*D66);

  for (int i = 1; i <= 10; i++)
  {
    if (R >= r[i-1][0] && R<= r[i][0])
    {
      double dr = r[i][0] - r[i-1][0];
      double db = r[i][1] - r[i-1][1];

      B1 = ((R-r[i-1][0])*db + r[i-1][1]*dr)/dr;
      break;
    }
  }

  double shearlbd;

  if (R > 1)
    shearlbd = 4*B1*sqrt(sqrt((D11*D22*D22*D22)))/(b*b*Nxy);
  else
    shearlbd = 4*B1*sqrt((D22*(D12+2*D66)))/(b*b*Nxy);

  double comb_minlbd = 1/(1/minlbd + 1/(shearlbd*shearlbd));

  // Constraint in the maximum contiguous layers with the same angle.

  int numcont = MaxContLay(layup);
  c[0] = double(numcont)/double(maxcont) - 1.0;

  // Constraint in the minimun percentual of 90 and 0 degree plies.

  int count90 = 0;
  int count0 = 0;

  for(int lam = 0; lam < layup.NCol( ); lam++)
  {
    if (layup[1][lam] == 90.0)
      count90++;

    if (layup[1][lam] == 0.0)
      count0++;
  }

  c[1] =0;
  c[2] =0;

  // Additional display constraints.

  c[3] = -fabs(minlbd);        // Buckling factor
  c[4] = -fabs(comb_minlbd);   // Combined Buckling factor
  c[5] = -shearlbd;            // Shear buckling factor

  // Objective function: maximization of the load factor.

  if (fabs(shearlbd) < comb_minlbd)
    fobjs[0] = -fabs(shearlbd);
  else
    fobjs[0] = -comb_minlbd;
}

// -------------------------------------------------------------------------
// Class cLamPltLoadFactor:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cLamPltLoadFactor ==========================

cLamPltLoadFactor :: cLamPltLoadFactor(void)
{
    NumObj = 1;
    NumConstr = 1;
}

// ============================== Evaluate =================================

void cLamPltLoadFactor :: Evaluate(int** algvar, cVector &c, cVector &fobjs)
{
  // Decode the variables.

  cMatrix layup;
  if (!Decode(algvar, layup))
  {
    cout << "Decoding process failure." << endl;
    exit(0);
  }

  // Evaluate the load factors.

  double lbdbck, lbdstr;
  Analysis(layup, lbdbck, lbdstr);

  // Constraint in the maximum contiguous layers with the same angle.

  int maxcont = 4;        // Max. num. of contiguous layers with same angle
  int numcont = MaxContLay(layup);
  c[0] = double(numcont)/double(maxcont) - 1.0;

  // Additional display constraints.

  // c[1] = -lbdbck;   // Buckling factor
  // c[2] = -lbdstr;   // Strength factor

  // Objective function: maximization of the load factor.

  if (lbdstr < lbdbck)
    fobjs[0] = -lbdstr;
  else
    fobjs[0] = -lbdbck;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cLamPltLoadFactor :: Analysis(cMatrix &layup, double &lbdb, double &lbds)
{
  // Problem data.

  double lim_e1  = 0.008; // Threshold e1 strain
  double lim_e2  = 0.029; // Threshold e2 strain
  double lim_e12 = 0.015; // Threshold e12 strain
  double Sf = 1.5;        // Strain safety factor

  // Default load and dimensions values

  double Nx  = 175.0;     // Axial load in the x-direction
  double Ny  = 21.875;    // Axial load in the y-direction
  double Nxy = 0.0;       // Shear load
  double Mx  = 0.0;       // Moment in the x-direction
  double My  = 0.0;       // Moment in the y-direction
  double Mxy = 0.0;       // Moment in the xy-direction

  double a = 0.508;       // Plate length
  double b = 0.127;       // Plate width

  if (!Loads.empty( ))
  {
    Nx  = Loads[0];
    Ny  = Loads[1];
    Nxy = Loads[2];
    Mx  = Loads[3];
    My  = Loads[4];
    Mxy = Loads[5];
  }

  if (!Dims.empty( ))
  {
    a = Dims[0];
    b = Dims[1];
  }

  // Compute the ABD matrix.

  cMatrix C(6,6);
  CalcABD(layup, C);

  // Get stiffness data.

  double D11 = C[3][3];
  double D12 = C[3][4];
  double D22 = C[4][4];
  double D66 = C[5][5];

  // Evaluate the minimum buckling load.

  int pmax = 21;
  int qmax = 21;
  double minlbd = 0.0;
  for (int p = 1; p < pmax; p++)
  {
    for (int q = 1; q < qmax; q++)
    {
      double pa2 = (p/a)*(p/a);
      double pa4 = pa2*pa2;
      double qb2 = (q/b)*(q/b);
      double qb4 = qb2*qb2;
      double lbd = (PI*PI*(D11*pa4 + 2.0*(D12 + 2.0*D66)*pa2*qb2 + D22*qb4))/(pa2*Nx + qb2*Ny);

      if (p == 1 && q == 1)
        minlbd = lbd;
      else if (lbd < minlbd)
        minlbd = lbd;
    }
  }

  // Solve [C]{genstrain} = {genstress}.

  cVector genstress(6);
  cVector genstrain(6);
  genstress[0] = Nx;
  genstress[1] = Ny;
  genstress[2] = Nxy;
  genstress[3] = Mx;
  genstress[4] = My;
  genstress[5] = Mxy;
  C.Solve(genstress, genstrain);

  // Get the membrane strains and curvatures.

  cVector memstrain(3);
  cVector curvature(3);
  memstrain[0] = genstrain[0];
  memstrain[1] = genstrain[1];
  memstrain[2] = genstrain[2];
  curvature[0] = genstrain[3];
  curvature[1] = genstrain[4];
  curvature[2] = genstrain[5];

  // Evaluate the laminate thickness.

  double lamthk = GetThickness(layup);

  // Evaluate strain criterion for each ply.

  int nlam = layup.NCol( );
  double z = -lamthk/2.0;
  double minstr;
  cVector epsl(3);
  cVector epsg(3);
  cMatrix T(3,3);
  for (int lam = 0; lam < nlam; lam++)
  {
    // Compute the global strains at the bottom of the layer.

    epsg[0] = memstrain[0] + z*curvature[0];
    epsg[1] = memstrain[1] + z*curvature[1];
    epsg[2] = memstrain[2] + z*curvature[2];

    // Compute the local strains.

    double theta = layup[1][lam]*PI/180;

    TMatrix(theta, T);

    epsl = T*epsg;

    // Store the smallest factor.

    if (lam == 0) minstr = lim_e1/(fabs(epsl[0])*Sf);
    if (minstr > lim_e1/(fabs(epsl[0])*Sf))  minstr = lim_e1/(fabs(epsl[0])*Sf);
    if (minstr > lim_e2/(fabs(epsl[1])*Sf))  minstr = lim_e2/(fabs(epsl[1])*Sf);
    if (minstr > lim_e12/(fabs(epsl[2])*Sf)) minstr = lim_e12/(fabs(epsl[2])*Sf);

    // Compute the global strains at the top of the layer.

    z += layup[0][lam];
    epsg[0] = memstrain[0] + z*curvature[0];
    epsg[1] = memstrain[1] + z*curvature[1];
    epsg[2] = memstrain[2] + z*curvature[2];

    // Compute the local strains.

    epsl = T*epsg;

    // Store the smallest factor.

    if (minstr > lim_e1/(fabs(epsl[0])*Sf))  minstr = lim_e1/(fabs(epsl[0])*Sf);
    if (minstr > lim_e2/(fabs(epsl[1])*Sf))  minstr = lim_e2/(fabs(epsl[1])*Sf);
    if (minstr > lim_e12/(fabs(epsl[2])*Sf)) minstr = lim_e12/(fabs(epsl[2])*Sf);
  }

  // Return load factors.

  lbdb = minlbd;   // Buckling factor
  lbds = minstr;   // Strength factor
}

// ============================ Evaluate ==============================

void cLamPltLoadFactor :: EvalExactConstraint(int index, int** algvar, double &c)
{
  cMatrix layup;
  if (!Decode(algvar, layup))
  {
    cout << "Decoding process failure." << endl;
    exit(0);
  }
  // Exact constraint evaluation.
    if (index == 0){
        int maxcont = 4;        // Max. num. of contiguous layers with same angle
        int numcont = MaxContLay(layup);
        c = double(numcont)/double(maxcont) - 1.0;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cLamPltLoadFactor :: GetApproxConstr(bool* approxc)
{
  approxc[0] = 0;
}

// -------------------------------------------------------------------------
// Class cLamPltMinLopez2009:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================ cLamPltMinLopez2009 ========================

cLamPltMinLopez2009 :: cLamPltMinLopez2009(void)
{
  NumConstr = 2;
  NumObj = 1;

  // F12 = 0 is adopted in TsaiWu criterio safety factor computation.
  cMatOrtho :: TsaiWuF12Param(0.0);
}

// ============================== Evaluate =================================

void cLamPltMinLopez2009 :: Evaluate(int** algvar, cVector &c, cVector &fobjs)
{
  // Problem data.

  double Sf      = 1.0;   // Safety Factor

  // Default load and dimensions values (case 01).

  double Nx  = 10000e3;   // Axial load in the x-direction
  double Ny  = 0.0;       // Axial load in the y-direction
  double Nxy = 0.0;       // Shear load
  double Mx  = 0.0;       // Moment in the x-direction
  double My  = 0.0;       // Moment in the y-direction
  double Mxy = 0.0;       // Moment in the xy-direction

  double a = 1;       // Plate length
  double b = 1;       // Plate width

  if (!Loads.empty( ))
  {
    Nx  = Loads[0];
    Ny  = Loads[1];
    Nxy = Loads[2];
    Mx  = Loads[3];
    My  = Loads[4];
    Mxy = Loads[5];
  }

  if (!Dims.empty( ))
  {
    a  = Dims[0];
    b  = Dims[1];
  }

  // Decode the variables.

  cMatrix layup;
  if (!Decode(algvar, layup))
    fobjs = 1e6;

  // Compute the ABD matrix.

  cMatrix C(6,6);
  CalcABD(layup, C);

  // Solve [C]{genstrain} = {genstress}.

  cVector genstress(6);
  cVector genstrain(6);
  genstress[0] = Nx;
  genstress[1] = Ny;
  genstress[2] = Nxy;
  genstress[3] = Mx;
  genstress[4] = My;
  genstress[5] = Mxy;
  C.Solve(genstress, genstrain);

  // Get the membrane strains and curvatures.

  cVector memstrain(3);
  cVector curvature(3);
  memstrain[0] = genstrain[0];
  memstrain[1] = genstrain[1];
  memstrain[2] = genstrain[2];
  curvature[0] = genstrain[3];
  curvature[1] = genstrain[4];
  curvature[2] = genstrain[5];

  // Evaluate the laminate thickness.

  double lamthk = GetThickness(layup);

  // Evaluate strain criterion for each ply.

  int nlam = layup.NCol( );
  int prevmat = -1;
  double z = -lamthk/2.0;
  double critSFinf;
  double critSFsup;
  double mincritSF = 1e12;
  cVector epsl(3);
  cVector epsg(3);
  cVector sigl(3);
  cMatrix T(3,3);
  cMatrix Ql(3,3);

  for (int lam = 0; lam < nlam; lam++)
  {
    // Find this layer's material ID.

    int idmat = (int)layup[2][lam];

    // Compute the global strains at the bottom of the layer.

    epsg[0] = memstrain[0] + z*curvature[0];
    epsg[1] = memstrain[1] + z*curvature[1];
    epsg[2] = memstrain[2] + z*curvature[2];

    // Evaluate the local constitutive matrix [Ql]. This matrix is computed
    // only when the layer material is different from the previous one.

    if (idmat != prevmat)
    {
      cMaterial *mat = cMaterial::GetMaterial(idmat);
      double *param = new double[mat->NumParam( )];
      mat->GetParam(param);
      QlMatrix(param, Ql);
      delete []param;
      prevmat = idmat;
    }

    // Compute the local strains and stresses.

    TMatrix(layup[1][lam]*PI/180.0, T);
    epsl = T*epsg;
    sigl = Ql*epsl;

    // Criteria Safety Factor.

    critSFinf = (cMaterial::GetMaterial(idmat)->SafetyFactor(sigl))/Sf;

    if (lam == 0) mincritSF = critSFinf;
    if (mincritSF > critSFinf) mincritSF = critSFinf;

    // Compute the global strains at the top of the layer.

    z += layup[0][lam];
    epsg[0] = memstrain[0] + z*curvature[0];
    epsg[1] = memstrain[1] + z*curvature[1];
    epsg[2] = memstrain[2] + z*curvature[2];

    // Compute the local strains and stresses.

    epsl = T*epsg;
    sigl = Ql*epsl;

    // Criteria Safety Factor.

    critSFsup = (cMaterial::GetMaterial(idmat)->SafetyFactor(sigl))/Sf;

    if (critSFsup < mincritSF) mincritSF = critSFsup;
  }

  // Plate weight.

  int lammat;
  double W = 0.0;
  double dens = 0.0;
  double area = a*b;
  for (int lam = 0; lam < nlam; lam++)
  {
    lammat = (int)layup[2][lam];
    dens = cMaterial::GetMaterial(lammat)->GetDensity( );
    W += 9.81*area*dens*layup[0][lam];
  }

  // Constraint in criteria safety factor.

  c[0] = 1.0 - mincritSF;
  c[1] = -mincritSF;

  // Objective function: weight reduction.

  fobjs[0] = W;
}

// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateMFBuckLam :: cSquarePlateMFBuckLam(void)
{
  NumConstr = 1;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateMFBuckLam :: Evaluate(int** algvar, cVector &c, cVector &fobjs)
{
    // Decode the variables.

    cMatrix layup;
    if (!Decode(algvar, layup))
    {
      cout << "Decoding process failure." << endl;
      exit(0);
    }

    // Evaluate the load factors.

    double lbdbck;
    Analysis(2, layup, lbdbck);

    // Constraint in the maximum contiguous layers with the same angle.

    int maxcont = 4;        // Max. num. of contiguous layers with same angle
    int numcont = MaxContLay(layup);
    c[0] = double(numcont)/double(maxcont) - 1.0;

    // Objective function: maximization of the load factor.

    fobjs[0] = -lbdbck;
}

// ============================== Evaluate =================================

void cSquarePlateMFBuckLam :: EvaluateLFP(int** algvar, cVector &c, cVector &fobjs)
{
    // Decode the variables.

    cMatrix layup;
    if (!Decode(algvar, layup))
    {
      cout << "Decoding process failure." << endl;
      exit(0);
    }

    // Evaluate the load factors.

    double lbdbck;
    Analysis(1, layup, lbdbck);

    // Constraint in the maximum contiguous layers with the same angle.

    int maxcont = 4;        // Max. num. of contiguous layers with same angle
    int numcont = MaxContLay(layup);
    c[0] = double(numcont)/double(maxcont) - 1.0;

    // Objective function: maximization of the load factor.

    fobjs[0] = -lbdbck;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Analysis =================================

void cSquarePlateMFBuckLam :: Analysis(int f, cMatrix &layup, double &lbdb)
{
    int num_thread = 0;
  #ifdef _OMP_
    num_thread = omp_get_thread_num( );
  #endif

    stringstream thread;
    thread << num_thread;

    string thread_number = thread.str();
    string cmd  = "del LamSqrPltBuck" + thread_number + ".dat";
    string cmd2 = "del LamSqrPltBuck" + thread_number + ".pos";
    string cmd3 = "rm LamSqrPltBuck" + thread_number + ".dat";
    string cmd4 = "rm LamSqrPltBuck" + thread_number + ".pos";

  #ifdef _WIN32
    if (system(cmd.c_str()) || system(cmd2.c_str()))
       cout << "Problem on removing plate.dat and plate.pos files.\n";
  #else
    if (system(cmd3.c_str()) || system(cmd4.c_str()))
       cout << "Problem on removing plate.dat and plate.pos files.\n";
  #endif

    string begname, endname;

    if (f == 1)
    {
        begname = "datbegSqrPltCCLamLF_4x4.dat";
        endname = "datendSqrPltLamLF_4x4.dat";
    }
    else if (f == 2)
    {
        begname = "datbegSqrPltCCLamHF.dat";
        endname = "datendSqrPltLamHF.dat";
    }

    string datname = "LamSqrPltBuck" + thread_number + ".dat";
    string posname = "LamSqrPltBuck" + thread_number + ".pos";


    #ifdef _WIN32
      cmd = "type " + begname + " >> " + datname;
    #else
      cmd = "cat " + begname + " >> " + datname;
    #endif

    int status = system(cmd.c_str( ));

    if (status)
    {
       cout << "Error in the copy of datbeg file.";
       exit(EXIT_FAILURE);
    }

    fstream dat;

    dat.open(datname.c_str( ));

    if (!dat.is_open( ))
    {
       cout << "Error opening the dat file for plate analysis." << endl;
       exit(0);
    }

    dat.seekp(0,ofstream::end);

    dat << endl << endl << "%SECTION.LAMINATED.SHELL" << endl;
    dat << "1" << endl;
    dat << "1 1.0 0.0 0.0 " << layup.NCol( ) << endl;

    for (int i = 0; i < layup.NCol(); i++)
        dat << "1   0.01    " << layup[1][i] << endl;

    dat << endl;

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
    cmd = "fast.exe LamSqrPltBuck" + thread_number + " -silent";
  #else
    cmd = "./fast LamSqrPltBuck" + thread_number + " -silent";
  #endif

    status = system(cmd.c_str( ));

    if (status)
    {
       cout << "Error in the analysis with fast.";
       exit(EXIT_FAILURE);
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

// ============================ Evaluate ==============================

void cSquarePlateMFBuckLam :: EvalExactConstraint(int index, int** algvar, double &c)
{
    // Decode the variables.

    cMatrix layup;
    if (!Decode(algvar, layup))
    {
      cout << "Decoding process failure." << endl;
      exit(0);
    }

  // Single constraint evaluation.
    if (index == 0)
    {
        int maxcont = 4;        // Max. num. of contiguous layers with same angle
        int numcont = MaxContLay(layup);
        c = double(numcont)/double(maxcont) - 1.0;
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateMFBuckLam :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}



// ========================== cSquarePlateBuckFGM ===========================

cSquarePlateKeshtegarBuckLam :: cSquarePlateKeshtegarBuckLam(void)
{
  NumConstr = 1;
  NumObj = 1;
}

// ============================== Evaluate =================================

void cSquarePlateKeshtegarBuckLam :: Evaluate(int** algvar, cVector &c, cVector &fobjs)
{
    int tb = MinThk*100; int tu = MaxThk*100; int tMax = 100;
    //cout << "StartRepair" << endl;
    RepairAlgorithm(algvar, tb, tu, tMax);
    //cout << "EndRepair" << endl;

    // Decode the variables.

    cMatrix layup;
    if (!Decode(algvar, layup))
    {
      cout << "Decoding process failure." << endl;
      exit(0);
    }

    // layup[0][0] = 0.11; layup[0][1] = 0.07; layup[0][2] = 0.12; layup[0][3] = 0.10; layup[0][4] = 0.10;
    // layup[0][9] = 0.11; layup[0][8] = 0.07; layup[0][7] = 0.12; layup[0][6] = 0.10; layup[0][5] = 0.10;
    // layup[1][0] = 35;   layup[1][1] = -15;  layup[1][2] = -47;  layup[1][3] = -38;  layup[1][4] = -2;
    // layup[1][9] = 35;   layup[1][8] = -15;  layup[1][7] = -47;  layup[1][6] = -38;  layup[1][5] = -2;

    // Evaluate the load factors.

    double lbdbck;
    Analysis(2, layup, lbdbck);

    // Constraint in the maximum total thickness.

    double MaxThk = 1.0;
    double thk = 0;
    for (int i = 0; i < layup.NCol( ); i++) thk += layup[0][i];
    c[0] = thk/MaxThk - 1.0;

    if (abs(c[0]) > 1e-6)
    {
        cout << "c[0] = " << c[0] << endl;
        cout << "thk = " << thk << endl;
        exit(0);
    }


    // Objective function: maximization of the load factor.

    fobjs[0] = -lbdbck;
}

// ============================== Evaluate =================================

void cSquarePlateKeshtegarBuckLam :: EvaluateLFP(int** algvar, cVector &c, cVector &fobjs)
{
    int tb = MinThk*100; int tu = MaxThk*100; int tMax = 100;
    RepairAlgorithm(algvar, tb, tu, tMax);

    // Decode the variables.

    cMatrix layup;
    if (!Decode(algvar, layup))
    {
      cout << "Decoding process failure." << endl;
      exit(0);
    }

    // Evaluate the load factors.

    double lbdbck;
    Analysis(1, layup, lbdbck);

    // Constraint in the maximum total thickness.

    double MaxThk = 1.0;
    double thk = 0;
    for (int i = 0; i < layup.NCol( ); i++) thk += layup[0][i];
    c[0] = thk/MaxThk - 1.0;

    if (abs(c[0]) > 1e-6)
    {
        cout << "c[0] = " << c[0] << endl;
        cout << "thk = " << thk << endl;
        exit(0);
    }

    // Objective function: maximization of the load factor.

    fobjs[0] = -lbdbck;
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ========================== RepairAlgorithm ===========================

void cSquarePlateKeshtegarBuckLam :: RepairAlgorithm(int** algvar, int tb, int tu, int tMax)
{
    if (LamType == SYMMETRIC)
    {
        //cout << "0" << endl;
        int SizeThk = NumCol;
        //cout << "0.1" << endl;
        cVector thk(SizeThk);
        //cout << "0.2" << endl;
        for (int i = 0; i < SizeThk; i++) thk[i] = algvar[0][i];
        //cout << "1" << endl;

        // thk in cm!
        for (int i = 0; i < SizeThk; i++) thk[i] = thk[i] + tb;
        //cout << "2" << endl;

        // Algorithm: List of priority

        cVector ListPrior(SizeThk);
        for (int i = 0; i < SizeThk; i++)
        {
            ListPrior[i] = i;
        }

        //cout << "3" << endl;
        // Redefine thickness so that sum < tMax

        int Sum = 0;
        for (int i = 0; i < (SizeThk - 1); i++) Sum += thk[ListPrior[i]];

        thk[ListPrior[SizeThk - 1]] = tMax/2 - Sum;

        //cout << "4" << endl;
        // Repair so that thickness for each ply is between the desired values (tb < t < tu)

        int id = 1;
        while (thk[ListPrior[SizeThk - 1]] < tb)
        {
            if (thk[ListPrior[SizeThk - 1 - id]] == tb)
            {
                id += 1;
                if (id == SizeThk) id = 1;
                continue;
            }
            thk[ListPrior[SizeThk - 1]] += 1;
            thk[ListPrior[SizeThk - 1 - id]] -= 1;

            id += 1;
            if (id == SizeThk) id = 1;
            continue;
        }

        //cout << "5" << endl;

        id = 1;
        while (thk[ListPrior[SizeThk - 1]] > tu)
        {
            if (thk[ListPrior[id - 1]] == tu)
            {
                id += 1;
                if (id == SizeThk) id = 1;
                continue;
            }
            thk[ListPrior[SizeThk - 1]] -= 1;
            thk[ListPrior[id - 1]]     += 1;

            id += 1;
            if (id == SizeThk) id = 1;
            continue;
        }

        //cout << "6" << endl;
        // Define algvar again
        for (int i = 0; i < SizeThk; i++) thk[i] = thk[i] - tb;

        for (int i = 0; i < SizeThk; i++)
        {
            algvar[0][i]                = thk[i];
            algvar[0][2*NumCol - i - 1] = thk[i];
        }

        //cout << "7" << endl;
    }
    else
    {
        cVector thk(NumCol);
        for (int i = 0; i < NumCol; i++) thk[i] = algvar[0][i];

        // thk in cm!
        for (int i = 0; i < NumCol; i++) thk[i] = thk[i] + tb;

        // Algorithm: List of priority
        int IndMid = NumCol/2;

        cVector ListPrior(NumCol);
        int Cont1 = 0; int Cont2 = 0;
        for (int i = 0; i < NumCol; i++)
        {
            int r = i%2;
            if (r == 0)
            {
                ListPrior[i] = Cont1; Cont1 += 1;
            }
            else
            {
                ListPrior[i] = NumCol - Cont2 - 1; Cont2 = Cont2 + 1;
            }
        }

        // Redefine thickness so that sum < tMax

        int Sum = 0;
        for (int i = 0; i < (NumCol - 1); i++) Sum += thk[ListPrior[i]];

        thk[ListPrior[NumCol - 1]] = tMax - Sum;

        // Repair so that thickness for each ply is between the desired values (tb < t < tu)

        int id = 1;
        while (thk[ListPrior[NumCol - 1]] < tb)
        {
            if (thk[ListPrior[NumCol - 1 - id]] == tb)
            {
                id += 1;
                if (id == NumCol) id = 1;
                continue;
            }
            thk[ListPrior[NumCol - 1]] += 1;
            thk[ListPrior[NumCol - 1 - id]] -= 1;

            id += 1;
            if (id == NumCol) id = 1;
            continue;
        }


        id = 1;
        while (thk[ListPrior[NumCol - 1]] > tu)
        {
            if (thk[ListPrior[id - 1]] == tu)
            {
                id += 1;
                if (id == NumCol) id = 1;
                continue;
            }
            thk[ListPrior[NumCol - 1]] -= 1;
            thk[ListPrior[id - 1]]     += 1;

            id += 1;
            if (id == NumCol) id = 1;
            continue;
        }

        // Define algvar again
        for (int i = 0; i < NumCol; i++) thk[i] = thk[i] - tb;

        for (int i = 0; i < NumCol; i++) algvar[0][i] = thk[i];
    }
}

// ============================== Analysis =================================

void cSquarePlateKeshtegarBuckLam :: Analysis(int f, cMatrix &layup, double &lbdb)
{
    int num_thread = 0;
  #ifdef _OMP_
    num_thread = omp_get_thread_num( );
  #endif

    stringstream thread;
    thread << num_thread;

    string thread_number = thread.str();
    string cmd  = "del LamSqrPltBuck" + thread_number + ".dat";
    string cmd2 = "del LamSqrPltBuck" + thread_number + ".pos";
    string cmd3 = "rm LamSqrPltBuck" + thread_number + ".dat";
    string cmd4 = "rm LamSqrPltBuck" + thread_number + ".pos";

  #ifdef _WIN32
    if (system(cmd.c_str()) || system(cmd2.c_str()))
       cout << "Problem on removing plate.dat and plate.pos files.\n";
  #else
    if (system(cmd3.c_str()) || system(cmd4.c_str()))
       cout << "Problem on removing plate.dat and plate.pos files.\n";
  #endif

    string begname, endname;

    if (f == 1)
    {
        begname = "datbegSqrPltKeshtegarLamLF_4x4.dat";
        endname = "datendSqrPltKeshtegarLamLF_4x4.dat";
    }
    else if (f == 2)
    {
        begname = "datbegSqrPltKeshtegarLamHF.dat";
        endname = "datendSqrPltKeshtegarLamHF.dat";
    }

    string datname = "LamSqrPltBuck" + thread_number + ".dat";
    string posname = "LamSqrPltBuck" + thread_number + ".pos";


    #ifdef _WIN32
      cmd = "type " + begname + " >> " + datname;
    #else
      cmd = "cat " + begname + " >> " + datname;
    #endif

    int status = system(cmd.c_str( ));

    if (status)
    {
       cout << "Error in the copy of datbeg file.";
       exit(EXIT_FAILURE);
    }

    fstream dat;

    dat.open(datname.c_str( ));

    if (!dat.is_open( ))
    {
       cout << "Error opening the dat file for plate analysis." << endl;
       exit(0);
    }

    dat.seekp(0,ofstream::end);

    dat << endl << endl << "%SECTION.LAMINATED.SHELL" << endl;
    dat << "1" << endl;
    dat << "1 1.0 0.0 0.0 " << layup.NCol( ) << endl;

    for (int i = 0; i < layup.NCol(); i++)
        dat << "1   " << layup[0][i] << "    " << layup[1][i] << endl;

    dat << endl;

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
    cmd = "fast.exe LamSqrPltBuck" + thread_number + " -silent";
  #else
    cmd = "./fast LamSqrPltBuck" + thread_number + " -silent";
  #endif

    status = system(cmd.c_str( ));

    if (status)
    {
       cout << "Error in the analysis with fast.";
       exit(EXIT_FAILURE);
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

// ============================ Evaluate ==============================

void cSquarePlateKeshtegarBuckLam :: EvalExactConstraint(int index, int** algvar, double &c)
{
    int tb = MinThk*100; int tu = MaxThk*100; int tMax = 100;
    RepairAlgorithm(algvar, tb, tu, tMax);

    // Decode the variables.

    cMatrix layup;
    if (!Decode(algvar, layup))
    {
      cout << "Decoding process failure." << endl;
      exit(0);
    }

  // Single constraint evaluation.
    if (index == 0)
    {
        // Constraint in the maximum total thickness.

        double MaxThk = 1.0;
        double thk = 0;
        for (int i = 0; i < layup.NCol( ); i++) thk += layup[0][i];
        c = thk/MaxThk - 1.0;

        if (abs(c) > 1e-6)
        {
            cout << "c = " << c << endl;
            cout << "thk = " << thk << endl;
            exit(0);
        }
    }
    else{
        cout << "Definition of an exact constraint missing!";
        exit(0);
    }
}

// ========================= GetApproxConstr ==========================

void cSquarePlateKeshtegarBuckLam :: GetApproxConstr(bool *approxc)
{
  approxc[0] = 0;
}


// ======================================================= End of file =====
