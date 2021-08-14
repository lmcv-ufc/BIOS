// -------------------------------------------------------------------------
// utl.c - utilitary functions.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>

#include "utl.h"

static unsigned int seed = 1; 

// -------------------------------------------------------------------------
// Public functions:
//

// ============================ UtlNextLabel ============================

int Utl :: NextLabel(ifstream &fp, char *s)
{
  while (fp.get() != '%')
    if (fp.eof())
      return(0);

  fp >> s;
  return(1);
}

// ==============================  ReadString  ==========================

int Utl :: ReadString(istream &fp, char *s)
{
  int i,c;
 
  while ((c = fp.get()) != '\'')
    if (fp.eof())
      return(0);

  for (i = 0; (c = fp.get()) != '\''; i++)
  {
    if (fp.eof()) return(0);
    s[i] = c;
  }

  s[i] = '\0';
  return(1);
}

// ============================== RandInt =============================

int Utl :: RandInt(int floor, int ceiling)
{
  return floor + rand() % (ceiling-floor+1);
}

// ============================== RandDec =============================

double Utl :: RandDec(void)
{
  return rand() / (RAND_MAX+1.0);
}

// ============================ RandDouble ============================

double Utl :: RandDouble(double floor, double ceiling)
{
  double f = (double)rand() / RAND_MAX;
  return floor + f * (ceiling - floor);
}

// ================================ Min ===============================

double Utl :: Min(double val1, double val2)
{
  if (val1 < val2)
    return val1;
  else
    return val2;
}

// ================================ Max ===============================

double Utl :: Max(double val1, double val2)
{
  if (val1 > val2)
    return val1;
  else
    return val2;
}

// ================================ Exit ==============================

void Utl :: Exit(string messeger)
{
  cout << endl << messeger << endl;
  exit(0);
}

// ================================ GetSeed ===========================

int  Utl :: GetSeed(void)
{
  return seed; 
}

// ================================ SetSeed ===========================

void Utl :: SetSeed(unsigned int s)
{
  seed = s;
  srand(seed);
}

// =============================================== End of file ========
