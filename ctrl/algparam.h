// -------------------------------------------------------------------------
// algparam.h - file containing the definition of the cAlgParam class.
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
// The cAlgParam class handle algorithm parameters that varies with 
// iteration/generation.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// cAlgParam(eParamType t, double vi, double vf)
//
//   t  - The type of algorithm parameter                          (in)
//   vi - Initial value                                            (in)
//   vf - Final Value                                              (in)
//
// This constructor stores variables used by the class.
// -------------------------------------------------------------------------
//
// void setLimits(double ti, double tf)
//
//   ti - Initial time variable (iteration/generation)             (in)
//   tf - Final time variable (iteration/generation)               (in)
//
// This method stores values for the initial and final iteration/generation.
// -------------------------------------------------------------------------
//
// void Evaluate(double t)
//
//   t - Current time variable (iteration/generation)              (in)
//
// This method compute the value of the parameter and hold it on Value 
// variable.
// -------------------------------------------------------------------------
//
// double GetValue(void)
//
// This method returns the current value of the parameter.
// -------------------------------------------------------------------------
//
// ifstream& operator>>(std::ifstream &in, cAlgParam &par) 
//
//   in  - stream object                                           (in)
//   par - Algorithm parameter object                              (out)
//
// This method read a given returns the problem type.
// -------------------------------------------------------------------------

#ifndef _ALGPARAM_H
#define _ALGPARAM_H

#include <istream>

using namespace std;

// -------------------------------------------------------------------------
// Algorithms parameters types:
//
typedef enum 
{
  CONSTANT,
  LINEAR, 
  EXPONENTIAL,
  PARABOLIC,
  PROBABILISTIC
} eParamType;

// -------------------------------------------------------------------------
// Definition of cAlgParam class:
//
class cAlgParam 
{
 protected:
  double           lim_inf,lim_sup; 
  double           t_inf,t_sup;
  double           Value;
  eParamType       Type;

 public:
  
                   cAlgParam(void);
                   cAlgParam(eParamType, double, double = 0);
                  ~cAlgParam(void);
         void      setLimits(double,double);
         void      Evaluate(double);
        
         double    GetValue(void) {return Value;}
    
  friend istream& operator>>(istream&,cAlgParam&);
};

#endif
