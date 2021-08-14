// -------------------------------------------------------------------------
// samp.h - file containing the definition of the cSamp class.
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
//
// The class cSamp contains different methods to define the initial sample
// used in Surrogate Models. Most routines employed were taken from:
//           https://people.sc.fsu.edu/~jburkardt/cpp_src/
//
// cSamp
// |- Hammersley Sequence Sampling
// |- Latin Hypercube Sampling
// |- N-Latin Hypercube Sampling
// |- Random Sampling
// |- Sobol Sequence Sampling
// |- Centroidal Voronoi Tesselation
// |- Latinized Centroidal Voronoi Tesselation
// |- Hammersley
// |- Hammersley
//
// Refs:
// FORRESTER, Alexander et al. Engineering design via surrogate modelling: a
// practical guide. John Wiley & Sons, 2008.
//
// STEPONAVICE, Ingrida; SHIRAZI-MANESH, Mojdeh; HYNDMAN, Rob J.; SMITH-MILES,
// Kate; VILLANOVA, Laura. On sampling methods for costly multi-objective
// black-box optimization. Optimization and Its Applications. v. 107, p. 273-296, 2016.
//
// CHO, Inyong; LEE, Yongbin; RYU, Dongheum; CHOI, Dong-Hoon. Comparison study of
// sampling methods for computer experiments using various performance measures.
// v. 55, p. 221-235, 2017.
//
// -------------------------------------------------------------------------

#ifndef _SAMP_H
#define _SAMP_H

#include <iostream>
#include <vector>
#include <cmath>

#include "vec.h"
#include "mat.h"
#include "utl.h"
#include "matvec.h"
#include "group.h"

using namespace std;

// -------------------------------------------------------------------------
// Forward Declarations:
//
class cVector;

// -------------------------------------------------------------------------
// Sampling methods:
//
typedef enum
{
  HAMMERSLEY,
  LHS,
  NLHS,
  RANDOM,
  SOBOL,
  CVT,
  LCVT
} eSampType;
istream &operator>>(istream&,eSampType&);

// -------------------------------------------------------------------------
// Definition of Sampling class:
//
class cSamp
{
 protected:
  int              NumVar;        // Number of variables
  int              NumSample;     // Number of samples
  cVector          Xlow;          // X lower bound
  cVector          Xupp;          // X upper bound

 public:
                   cSamp(void);
  virtual          ~cSamp(void);
          void     Read(ifstream&, int&, int&, int&, vector<cVector> &, vector<cVector> &);
          void     InitSample(eSampType, int, int, vector<cVector> &);
          void     CreateSample(eSampType, vector<cVector> &);
          void     CalcSampHammersley(int, int, vector<cVector>  &);
          void     CalcSampLHS(int, int, int&,vector<cVector>  &);
          void     CalcSampRandom(int, int, int&,vector<cVector>  &);
          void     CalcSampSobol( int, int, int, vector <cVector> &);
          void     CalcSampNLHS(int, int, vector<cVector>  &);
          void     CalcSampCVT (int, int, int, int, int, int, int, int, int* , double* ,int *, double *, double *, vector<cVector>  & );
          void     CalcSampLCVT (int, int, int, int, int, int, int, int, int* , double* ,int *, double *, double *, vector<cVector>  & );

          int           i4vec_sum ( int, cVector );
          int           prime ( int n );
          double       *r8mat_uniform_01_new ( int m, int n, int &seed );
          int          *perm_uniform_new ( int n, int &seed );
          int           i4_uniform_ab ( int ilo, int ihi, int &seed );
          int           get_seed ( );
          void          i8_sobol ( int dim_num, long long int *seed, double quasi[ ] );
          int           i8_bit_lo0 ( long long int n );
          void          cvt_sample ( int dim_num, int n, int n_now, int sample, bool initialize, int *seed, double r[] );
          void          cvt_iterate ( int dim_num, int n, int batch, int sample, bool initialize, int sample_num, int *seed, double r[], double *it_diff, double *energy );
          int           i4_min ( int i1, int i2 );
          void          find_closest ( int dim_num, int n, int sample_num, double s[], double r[], int nearest[] );
          unsigned long random_initialize ( int seed );
          void          r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
          void          i4_to_halton_sequence ( int dim_num, int n, int step, int seed[], int leap[], int base[], double r[] );
          void          tuple_next_fast ( int m, int n, int rank, int x[] );
          void          user ( int dim_num, int n, int *seed, double r[] );
          double        r8_huge ( void );
          bool          halham_leap_check ( int dim_num, int leap[] );
          bool          halham_n_check ( int n );
          bool          halham_dim_num_check ( int dim_num );
          bool          halham_seed_check ( int dim_num, int seed[] );
          bool          halham_step_check ( int step );
          bool          halton_base_check ( int dim_num, int base[] );
          void          r8mat_latinize ( int m, int n, double table[] );
          int          *r8vec_sort_heap_index_a ( int n, double a[] );
};

#endif
