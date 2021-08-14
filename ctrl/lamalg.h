// -------------------------------------------------------------------------
// lamalg.h - file containing the definition of the cLaminateAlgorithm class.
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
// The Laminated Algorithm class store all the input data used on the 
// optimization algorithms which use specific laminated operators
//
// -------------------------------------------------------------------------
// Static methods:
// -------------------------------------------------------------------------
//
//  void ReadSwapRate(void)
//  void ReadSwapRange(void)
//  void ReadAddRate(void)
//  void ReadAddRange(void)
//  void ReadDelRate(void)
//  void ReadDelRange(void)
// -------------------------------------------------------------------------

#ifndef _LAMALG_H
#define _LAMALG_H

// -------------------------------------------------------------------------
// Definition of the Laminate Algorithm class:
//
class cLamAlg
{
 protected:
  static double LamMutProb[3];        // Laminate mutate probability
  static double SwapProb;             // Layer swap probability
  static double MaxSwap;              // Maximum layer swap probability
  static double MinSwap;              // Minimum layer swap probability
  static double AddProb;              // Layer addition probability
  static double MaxAdd;               // Maximum layer addition probability
  static double MinAdd;               // Minimum layer addition probability
  static double DelProb;              // Layer deletion probability
  static double MaxDel;               // Maximum layer deletion probability
  static double MinDel;               // Minimum layer deletion probability

 public:
  static  void    ReadLamMutProb(void);
  static  void    ReadSwapProb(void);
  static  void    ReadSwapRange(void);
  static  void    ReadAddProb(void);
  static  void    ReadAddRange(void);
  static  void    ReadDelProb(void);
  static  void    ReadDelRange(void);
 
                   cLamAlg(void);
                  ~cLamAlg(void);
};

#endif
