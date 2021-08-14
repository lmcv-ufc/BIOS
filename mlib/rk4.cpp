// -------------------------------------------------------------------------
// rk4.cpp - implementation of the 4th order Runge-Kutta Method to solve
//           Initial Value Problems (IVPs)
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
// Created:      10-Jun-2016    Evandro Parente Junior
//
// Modified:
// -------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "rk4.h"
#include "vec.h"

// -------------------------------------------------------------------------
// Public methods:
//

// ================================== RK4 ==================================

RK4 :: RK4(FuncDer *df, int dim)
{
  Fder = df;
  n = dim;
  f.Resize(n);
  k1.Resize(n);
  k2.Resize(n);
  k3.Resize(n);
  k4.Resize(n);
}

// ================================= ~RK4 ==================================

RK4 :: ~RK4(void)
{
}

// ================================= Step ==================================

void RK4 :: Step(void *data, double t0, double t, cVector &x0, cVector &x)

{
  // Compute the time step (increment).

  double dt = t - t0;

  // Evaluate k1 = dt*f(t0, x0).

  Fder(data, t0, x0, f);
  k1 = dt*f;

  // Evaluate k2 = dt*f(t0 + 0.5*dt, x0 + 0.5*k1).

  x = x0 + 0.5*k1;
  Fder(data, t0 + 0.5*dt, x, f);
  k2 = dt*f;

  // Evaluate k3 = dt*f(t0 + 0.5*dt, x0 + 0.5*k2).

  x = x0 + 0.5*k2;
  Fder(data, t0 + 0.5*dt, x, f);
  k3 = dt*f;

  // Evaluate k4 = dt*f(t0 + dt, x0 + k3).

  x = x0 + k3;
  Fder(data, t0 + dt, x, f);
  k4 = dt*f;

  // Evaluate the new point.

  for (int i = 0; i < n; i++)
  {
    x[i] = x0[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
  }
}

// ================================================ End of file ============
