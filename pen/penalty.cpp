// -------------------------------------------------------------------------
// penalty.cpp - implementation of cPenalty class.
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
// Created:      21-Apr-2012    Iuri Barcelos Rocha
//
// Modified:     15-Mar-2013    Evandro Parente Junior
//               Creation of cPopulation class.
//
// Modified:     06-Nov-2014    Elias Saraiva Barroso
//               Define cPenalty class methods considering the new cGroup class
//               hierarchy.
//
// Modified:     30-Aug-2017    Marina Alves Maia
//               Creation of new type of penalty (Normalization), which must
//               be used only on the NSGA or LamNSGA algorithms.
// -------------------------------------------------------------------------

#include <cmath>
#include <iostream>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#include "penalty.h"
#include "group.h"
#include "optsolution.h"
#include "utl.h"
#include "input.h"
#include "gbldef.h"
#include "gblvar.h"

// -------------------------------------------------------------------------
// Public methods:
//


// ============================= CreatePenalty =============================

cPenalty *cPenalty :: CreatePenalty(ePenType type)
{
  cPenalty *pen = 0;

  switch(type)
  {
    case STATIC:
      pen = new cPenStatic( );
    break;
      
    case DEB2000:
      pen = new cPenDeb( );
    break;
      
    case ADAPTIVE:
      pen = new cPenAdaptive( );
    break;

    case NORMALIZATION:
       pen = new cNormalization( );
    break;
  }
  
  return(pen);
}

// =============================== cPenalty ================================

cPenalty :: cPenalty(void)
{
}

// ============================== ~cPenalty ================================

cPenalty :: ~cPenalty(void)
{
}

// -------------------------------------------------------------------------
// Class cPenStatic:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== cPenStatic ===============================

cPenStatic :: cPenStatic(void)
{
  Type = STATIC;
  K    = 1.0e6;
}

// ============================= ~cPenStatic ===============================

cPenStatic :: ~cPenStatic(void)
{
}

// ============================== ReadFactor ===============================

void cPenStatic :: ReadFactor(istream &in)
{
  if (!(in >> K))
  {
    cout << "Error in the input of the static penalty factor." << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cPenStatic :: LoadReadFunc(cInpMap &im)
{
  // Register read functions.
  im.Insert("PENALTY.CONSTANT.FACTOR", makeReadObj(cPenStatic,ReadFactor));
}

// ============================ EvalPenObjFunc =============================

void cPenStatic :: EvalPenObjFunc(cGroup *group, double tol)
{
  int groupsize = group->GetSize( );

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < groupsize; i++)
  {                                       
    double fobj = group->GetSol(i)->GetCurrObjFunc(0);
    cVector constr = group->GetSol(i)->GetCurrConstr( );       
                                                    
    if (constr.Max( ) <= tol)        // Feasible individual
      group->GetSol(i)->AssignPenObjFunc(fobj);
    else                             // Unfeasible individual
    {
      double fpen = 0.0;
      for (int j = 0; j < constr.Dim( ); j++)
        if (constr[j] > 0.0) fpen += K*constr[j];

      group->GetSol(i)->AssignPenObjFunc(fobj + fpen);
    }
  }
}

// -------------------------------------------------------------------------
// Class cPenDeb:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// =============================== cPenDeb =================================

cPenDeb :: cPenDeb(void)
{
  Type = DEB2000;
}

// ============================== ~cPenDeb =================================

cPenDeb :: ~cPenDeb(void)
{
}

// ============================ EvalPenObjFunc =============================

void cPenDeb :: EvalPenObjFunc(cGroup *group, double tol)
{
  int groupsize = group->GetSize( );
  double maxobj = group->GetSol(0)->GetCurrObjFunc(0);
  for (int i = 1; i < groupsize; i++)
    maxobj = MAX(maxobj, group->GetSol(i)->GetCurrObjFunc(0));
  
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < groupsize; i++)
  {
    double fobj = group->GetSol(i)->GetCurrObjFunc(0);
    cVector constr = group->GetSol(i)->GetCurrConstr( );
    
    if (constr.Max( ) <= tol)        // Feasible individual
      group->GetSol(i)->AssignPenObjFunc(fobj);
    else                             // Unfeasible individual
    {
      double fpen = 0.0;
      for (int j = 0; j < constr.Dim( ); j++)
        if (constr[j] > 0) fpen += constr[j];

      group->GetSol(i)->AssignPenObjFunc(maxobj + fpen);
    }
  }
}

// -------------------------------------------------------------------------
// Class cPenAdaptive:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cPenAdaptive ==============================

cPenAdaptive :: cPenAdaptive(void)
{
  Type = ADAPTIVE;
}

// ============================ ~cPenAdaptive ==============================

cPenAdaptive :: ~cPenAdaptive(void)
{
}

// ============================ EvalPenObjFunc =============================

void cPenAdaptive :: EvalPenObjFunc(cGroup *group, double tol)
{
  // Evaluate objective function and constraint violation sums.
  
  int groupsize = group->GetSize( );
  cVector constr = group->GetSol(0)->GetCurrConstr( );
  int nc = constr.Dim( );
  cVector sumviol(nc);
  sumviol.Zero( );
  double sumobj = 0.0;
  for (int i = 0; i < groupsize; i++)
  {
    sumobj += group->GetSol(i)->GetCurrObjFunc(0);

    constr = group->GetSol(i)->GetCurrConstr( );
    for (int j = 0; j < nc; j++)
      if (constr[j] > 0) sumviol[j] += constr[j];
  }
  
  // Evaluate the square sum of constraint violations.
  
  double sqrviol = 0.0;
  for (int i = 0; i < nc; i++) sqrviol += sumviol[i]*sumviol[i];
  
  // Evaluate the penalty parameters.
  
  cVector k(nc);
  k.Zero( );
  
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < nc; i++)
    if (abs(sumobj)*sumviol[i]/sqrviol) 
      k[i] = abs(sumobj)*sumviol[i]/sqrviol;
  
  // Evaluate the penalized objective function.
  
  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < groupsize; i++)
  {
    double  fobj = group->GetSol(i)->GetCurrObjFunc(0);
    cVector constr = group->GetSol(i)->GetCurrConstr( );
    
    if (constr.Max( ) <= tol)        // Feasible individual
      group->GetSol(i)->AssignPenObjFunc(fobj);
    else                             // Unfeasible individual
    {
      double fpen = 0.0;
      for (int j = 0; j < constr.Dim( ); j++)
        if (constr[j] > 0) fpen += k[j]*constr[j];
         
      if (fobj > sumobj/groupsize)
        group->GetSol(i)->AssignPenObjFunc(fobj + fpen);
      else
        group->GetSol(i)->AssignPenObjFunc((sumobj/groupsize) + fpen);
    }
  }
}

// -------------------------------------------------------------------------
// Class cNormalization:
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Public methods:
//

// ============================= cNormalization ==============================

cNormalization :: cNormalization(void)
{
  Type = NORMALIZATION;
}

// ============================ ~cNormalization ==============================

cNormalization :: ~cNormalization(void)
{
}

// ============================ EvalNormalization =============================

void cNormalization :: EvalPenObjFunc(cGroup *group, double tol)
{
  // Evaluate objective function and constraint violation sums.

  int groupsize = group->GetSize( );
  cVector constr = group->GetSol(0)->GetCurrConstr( );
  int nc = constr.Dim( );
  cVector maxconst(nc);
  maxconst.Zero( );

  // Evaluate the maximum values of each constraint

  for (int j = 0; j < nc; j++)
  {
      for (int i = 0; i < groupsize; i++)
      {
          constr = group->GetSol(i)->GetCurrConstr( );

//          cout << "Restricao " << j << " do individuo " << i << ": " << constr[j] << endl;

          if (constr[j] > 0)
          {
              if (constr[j] > maxconst[j])
              {
                  maxconst[j] = constr[j];
              }
              else
              {
                  maxconst[j] = maxconst[j];
              }
          }
      }
  }


  // Normalizes all constraints according to their maximum values

   for (int i = 0; i < groupsize; i++)
   {
       double normalizedconst = 0.0;

       constr = group->GetSol(i)->GetCurrConstr( );

       for (int j = 0; j < nc; j++)
       {
           if (maxconst[j]==0)
           {normalizedconst+=0;}
           else if (constr[j] > 0)
           {normalizedconst=normalizedconst+constr[j]/maxconst[j];}
       }
       group->GetSol(i)->AssignNormConst(normalizedconst);
   }

   // Evaluate objective function and constraint sums.

  cVector maxconstfeasible(nc);
  maxconstfeasible.Zero( );

  // Evaluate the maximum values of each constraint

  for (int j = 0; j < nc; j++)
  {
      for (int i = 0; i < groupsize; i++)
      {
          constr = group->GetSol(i)->GetCurrConstr( );

          if (constr[j] < maxconstfeasible[j])
          { maxconstfeasible[j] = constr[j]; }
          else
          { maxconstfeasible[j] = maxconstfeasible[j];}
      }
  }

  // Normalizes all constraints according to their maximum values

   for (int i = 0; i < groupsize; i++)
   {
       double normalizedconstfeasible = 0.0;

       constr = group->GetSol(i)->GetCurrConstr( );

       for (int j = 0; j < nc; j++)
       {
           if (maxconstfeasible[j]==0)
           {normalizedconstfeasible+=0;}
           else if (constr[j] <= 0)
           {normalizedconstfeasible=normalizedconstfeasible+constr[j]/maxconstfeasible[j];}
       }
       group->GetSol(i)->AssignFeasibleConst(normalizedconstfeasible);
   }
}

// ======================================================= End of file =====
