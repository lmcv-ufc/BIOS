// ------------------------------------------------------------------------
// material.h - file containing the definition of the Material class.
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
// -----------------------------------------------------------------------
//
// The Material class basically reads the material parameters from the
// input file. These parameters are stored and can be queried by others
// when these parameters are necessary.
//
// Material
// |-- Isotropic
// |-- Orthotropic
//
// -------------------------------------------------------------------------

#ifndef _MATERIAL_H
#define _MATERIAL_H

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Failure criterium types:
//
typedef enum
{
  VON_MISES,
  MAX_STRAIN,
  MAX_STRESS,
  TSAI_HILL,
  TSAI_WU
} eFailType;

// -------------------------------------------------------------------------
// Material types:
//
typedef enum
{
  MAT_ISOTROPIC,       // Elastic isotropic material
  MAT_ORTHOTROPIC      // Elastic orthotropic material
} eMatType;

// ------------------------------------------------------------------------
// Definition of the Material class:
//
class cMaterial
{
 private:
  static int         NumMat;
  static int         NumIso;
  static int         NumOrtho;
  static cMaterial **VecMat;

 protected:
  eMatType  Type;      // Material type
  int       Label;     // Material label
  double    Density;   // Mass density
  double    Cost;      // Material cost

 public:
  static void       ReadNumMat(std::istream&);
  static void       ReadIso(std::istream&);
  static void       ReadOrtho(std::istream&);
  static void       ReadDensity(std::istream&);
  static void       ReadCost(std::istream&);
  static void       Destroy(void);
  static int        GetNumMat(void) { return NumMat; }
  static int        GetNumIso(void) { return NumIso; }
  static int        GetNumOrtho(void) { return NumOrtho; }
  static cMaterial *GetMaterial(int);

                    cMaterial(int);
  virtual          ~cMaterial(void);
          eMatType  GetType(void) { return Type; }
          int       GetLabel(void) { return Label; }
          double    GetDensity(void) {return Density; }
          double    GetCost(void) {return Cost; }
  virtual int       NumParam(void) = 0;
  virtual void      GetParam(double *) = 0;
  virtual void      Read(void) = 0;
  virtual double    SafetyFactor(cVector&) = 0;
  virtual double    SafetyFactor(cVector&, int&) = 0;
};

// ------------------------------------------------------------------------
// Definition of the Elastic Isotropic Material class:
//
class cMatIso : public cMaterial
{
 protected:
  static eFailType FailType; // Failure criterium
         double    E;        // Young's modulus
         double    Nu;       // Poisson's ratio
         double    SMYS;     // Specified minimum yield stress at room temperature based on the engineering stress-strain curve
         double    SMTS;     // Specified minimum tensile strength at room temperature based on the engineering stress-strain curve

 public:
  static   void   ReadFailType(std::istream&);

                  cMatIso(int);
  virtual        ~cMatIso(void);
           void   Read(void);
           void   GetParam(double *);
           int    NumParam(void) { return 4; } //E, Nu, SMYS, SMTS
           double SafetyFactor(cVector&);
           double SafetyFactor(cVector&, int&);
           double VonMises2D(cVector&);
           double VonMises3D(cVector&);
};

// ------------------------------------------------------------------------
// Definition of the Elastic Orthotropic Material class:
//
class cMatOrtho : public cMaterial
{
 protected:
  static eFailType FailType;                   // Failure criterium
         double    E1,E2,E3;                   // Young's moduli
         double    Nu12,Nu13,Nu23;             // Poisson's ratios
         double    G12,G13,G23;                // Shear moduli
         double    Xt,Xc,Yt,Yc,Zt,Zc,S4,S5,S6; // Strength parameters

 public:
  static   void   ReadFailType(std::istream&);
  static   void   TsaiWuF12Param(double);

                  cMatOrtho(int);
  virtual        ~cMatOrtho(void);
           void   Read(void);
           void   GetParam(double *);
           int    NumParam(void) { return 9; }
           double SafetyFactor(cVector&);
           double SafetyFactor(cVector&, int&);    // SF with failure type identification (for progressive failure)
	   double MaxStrain2D(cVector&);
	   double MaxStrain2D(cVector&, int&);
	   double MaxStrain3D(cVector&);
           double MaxStress2D(cVector&);
           double MaxStress2D(cVector&, int&);
           double MaxStress3D(cVector&);
           double TsaiHill(cVector&);
           double TsaiWu2D(cVector&);
           double TsaiWu2D(cVector&, int&);
           double TsaiWu3D(cVector&);
};

#endif
