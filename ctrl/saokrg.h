// -------------------------------------------------------------------------
// saokrg.h - file containing the definition of the cSAOKRG class.
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
// The cSAOKRG class implements routines related to performing SAO using
// Kriging, such as:
//   - CreateSurrogate( )
//   - UpdateSurrogate( )
//
// -------------------------------------------------------------------------

#ifndef _EGO_H
#define _EGO_H

#include "sao.h"
#include "krg.h"

// -------------------------------------------------------------------------
// Definition of the SAOKRG class:
//
class cSAOKRG : public cSAO
{
 protected:
    cKRG             *krg;
    eCorrelationType  CorrType;
    double            HyperParamLow;
    double            HyperParamUpp;


 public:
                    cSAOKRG(void);
  virtual           ~cSAOKRG(void);

  virtual void      LoadReadFunc(cInpMap&);
  void              ReadCorrType(std::istream&);
  void              ReadHyperParam(std::istream&);

  cSURR*    CreateSurrogate(sSampData&);
  void      UpdateSurrogate(cVectorVec&,cVectorVec&);
};

#endif
