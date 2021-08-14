// -------------------------------------------------------------------------
// algparam.cpp - implementation of algorithm parameter class.
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
// Created:      20-Jul-2014    Elias Saraiva Barroso
//
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

#include <string>
#include "algparam.h"
#include <stdlib.h>

// -------------------------------------------------------------------------
// Public methods:
//

// ================================ cAlgParam ==============================

cAlgParam :: cAlgParam(void)
{
}

// ================================ ~cAlgParam =============================

cAlgParam :: ~cAlgParam(void)
{
}

// ================================ cAlgParam ==============================

cAlgParam :: cAlgParam(eParamType t, double vi, double vf)
{
  Type = t;
  Value = 0;

  switch (t)
  {
    case CONSTANT:
      lim_inf = vi;
    break;

    case LINEAR:
    case EXPONENTIAL:
    case PARABOLIC:
    case PROBABILISTIC:
      lim_inf = vi;
      lim_sup = vf;
    break;
  }
}

// ================================ Evaluate ===============================

void cAlgParam :: Evaluate(double t)
{
  switch (Type)
  {
    case CONSTANT:
      Value = lim_inf;
    break;

    case LINEAR:
      Value = (lim_inf + ((lim_sup-lim_inf)/(t_sup-t_inf))*t);
    break;

    case EXPONENTIAL:
      Value = (lim_inf*pow((lim_inf/lim_sup),(t-t_inf)/(t_inf-t_sup)));
    break;

    case PARABOLIC:
      Value = (lim_inf - lim_sup) * pow((t/t_sup),2) - 2*(lim_inf - lim_sup)*(t/t_sup) + lim_inf;
    break;

    case PROBABILISTIC:
      return;
    break;
  }
}

// ================================ setParam  ==============================

void cAlgParam :: setLimits(double ti, double tf)
{
  switch (Type)
  {
    case CONSTANT:
    break;

    case LINEAR:
    case EXPONENTIAL:
    case PARABOLIC:
    case PROBABILISTIC:
      t_inf = ti;
      t_sup = tf;
    break;
  }
}

// ================================ operator>>  ==============================

istream&  operator>> (std::istream &in, cAlgParam &par)
{ 
  string str;

  in >> str; 

  if (str == "'CONSTANT'")
  {
    par.Type = CONSTANT;
    in >> par.lim_inf;
  }
  else if (str == "'LINEAR'")
  {
    par.Type = LINEAR;
    in >> par.lim_inf;
    in >> par.lim_sup;
  } 
  else if (str == "'EXPONENTIAL'")
  {
    par.Type = EXPONENTIAL;
    in >> par.lim_inf;
    in >> par.lim_sup;
  }
  
  else if (str == "'PARABOLIC'")
  {  
    par.Type = PARABOLIC;
    in >> par.lim_inf;
    in >> par.lim_sup;
  }
  else if (str == "'PROBABILISTIC'")
  {  
    par.Type = PARABOLIC;
    in >> par.lim_inf;
    in >> par.lim_sup;
  }
  else
  {
    cout << "Unknown AlgParameter type: " << str << endl;
    exit(0);
  }

  return in;
}

// ======================================================= End of file =====
