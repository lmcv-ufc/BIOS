// -------------------------------------------------------------------------
// lamalg.cpp - implementation of cLamAlg class.
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
// Created:      20-Aug-2014    Elias Saraiva Barroso
//
// Modified:     
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

#include "lamalg.h"
#include "input.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Static variables:
//
double cLamAlg :: LamMutProb[3] = {0.00, 0.00, 0.00};
double cLamAlg :: SwapProb      = 0.05;
double cLamAlg :: MaxSwap       = 0.00;
double cLamAlg :: MinSwap       = 0.00;
double cLamAlg :: AddProb       = 0.00;
double cLamAlg :: MaxAdd        = 0.00;
double cLamAlg :: MinAdd        = 0.00;
double cLamAlg :: DelProb       = 0.00;
double cLamAlg :: MaxDel        = 0.00;
double cLamAlg :: MinDel        = 0.00;

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= ReadLamMutProb ============================

void cLamAlg :: ReadLamMutProb(void)
{
  if (!(in >> LamMutProb[0]) || !(in >> LamMutProb[1]) || 
      !(in >> LamMutProb[2]))
  {
    cout << "Error in the input of the laminate mutation." << endl;
    exit(0);
  }
}

// ============================= ReadSwapProb ==============================

void cLamAlg :: ReadSwapProb(void)
{
  if (!(in >> SwapProb))
  {
    cout << "Error in the input of the layer swap probability." << endl;
    exit(0);
  }
}

// ============================= ReadSwapRange =============================

void cLamAlg :: ReadSwapRange(void)
{
  if (!(in >> MinSwap) || !(in >> MaxSwap))
  {
    cout << "Error in the input of the layer swap range." << endl;
    exit(0);
  }
}

// ============================== ReadAddProb ==============================

void cLamAlg :: ReadAddProb(void)
{
  if (!(in >> AddProb))
  {
    cout << "Error in the input of the layer addition probability." << endl;
    exit(0);
  }
}

// ============================== ReadAddRange =============================

void cLamAlg :: ReadAddRange(void)
{
  if (!(in >> MinAdd) || !(in >> MaxAdd))
  {
    cout << "Error in the input of the layer addition range." << endl;
    exit(0);
  }
}

// ============================== ReadDelProb ==============================

void cLamAlg :: ReadDelProb(void)
{
  if (!(in >> DelProb))
  {
    cout << "Error in the input of the layer deletion probability." << endl;
    exit(0);
  }
}

// ============================== ReadDelRange =============================

void cLamAlg :: ReadDelRange(void)
{
  if (!(in >> MinDel) || !(in >> MaxDel))
  {
    cout << "Error in the input of the layer deletion range." << endl;
    exit(0);
  }
}

// ============================= cLamAlg ===================================

cLamAlg :: cLamAlg(void)
{
}

// ============================= ~cLamAlg ==================================

cLamAlg :: ~cLamAlg(void)
{
}

// ======================================================= End of file =====
