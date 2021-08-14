// ------------------------------------------------------------------------
// gbldef.h - global definitions.
// ------------------------------------------------------------------------
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

#ifndef _GBLDEF_H
#define _GBLDEF_H

// ------------------------------------------------------------------------
// Useful constants:
//
const double PI = 3.141592653589793;
const double GR = 9.80665; 

// ------------------------------------------------------------------------
// Useful functions:
//
template <typename T>
inline T MAX(T x, T y)
{
  if (x < y)
    return y;
  else
    return x;
}

template <typename T>
inline T MIN(T x, T y)
{
  if (x > y)
    return y;
  else
    return x;
}

inline double rad2deg(double r)                                                      
{                                                                             
  return(180.0*r/PI);                                                         
}                                                                             
                                                                              
inline double deg2rad(double d)                                                      
{                                                                             
  return(PI*d/180.0);                                                         
}

#endif
