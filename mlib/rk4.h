// -------------------------------------------------------------------------
// rk4.h - Solution of IVPs using the Runge-Kutta Method of 4th-order.
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
// The class RK4 is a implementation of the Runge-Kutta Method of 4th-order
// to the solution of Initial Value Problems (IVP's), defined as:
//
//                             dx/dt = f(t, x)
//                             x(t0) = x0
//
// where t is an scalar (e.g. time), x(t) is the unknown function defined
// in Rn, and f(t, x) is a known vector function in Rn. Thus, the IVP
// is a system of n first-order differential equations.
//
// -------------------------------------------------------------------------
// Client functions:
// -------------------------------------------------------------------------
//
// typedef void (FuncDer)(void *data, double t, cVector &x, cVector &f)
//
//   data - problem data                                               (in)
//   t    - given time                                                 (in)
//   x    - vector of variables at t                                   (in)
//   f    - vector containing the first derivates                     (out)
//
// This method should be implemented by the clients of the RK4 class to
// compute the first derivatives at a given point f(t, x) = dx/dt. The
// variable data may contain additional parameters used in the computation
// of the first derivatives.
//
// -------------------------------------------------------------------------
// Public methods:
// -------------------------------------------------------------------------
//
// RK4(FuncDer *f, int n)
//
//   f  -  function to compute the first derivatives                   (in)
//   n  -  problem dimension                                           (in)
//
// Constructor method.
// -------------------------------------------------------------------------
//
// void Step(void *data, double t0, double t, cVector &x0, cVector &x)
//
//   data - problem data                                               (in)
//   t0   - initial time                                               (in)
//   t    - final time                                                 (in)
//   x0   - vector of variables at t0                                  (in)
//   x    - vector of variables at t                                  (out)
//
// This method integrates the system of first-order diffential equation from
// t to t0 using one step of the RK4 method (step size h = t - t0). It can
// be called several times by the client for increasing t values in order
// to obtain the time history of x(t).
// -------------------------------------------------------------------------

#ifndef _RK4_H
#define _RK4_H

#include "vec.h"

// -------------------------------------------------------------------------
// Auxiliary types:
//
typedef void (FuncDer)(void *, double, cVector &, cVector &);

// -------------------------------------------------------------------------
// Definition of RK4 class:
//
class RK4
{
 protected:
  FuncDer *Fder; // Function to compute the derivates
  int      n;    // Problem dimension
  cVector  f;    // Vector of derivatives
  cVector  k1;   // RK4 vector
  cVector  k2;   // RK4 vector
  cVector  k3;   // RK4 vector
  cVector  k4;   // RK4 vector

 public:
                RK4(FuncDer *, int);
  virtual      ~RK4(void);

          void  Step(void *, double, double, cVector &, cVector &);
};

#endif
