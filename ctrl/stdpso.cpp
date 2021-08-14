// -------------------------------------------------------------------------
// stdpso.cpp - implementation of cStandardPSO class.
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
// Created:      26-Nov-2014    Elias Saraiva Barroso
//
// Modified:
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <math.h>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "stdpso.h"
#include "algparam.h"
#include "problem.h"
#include "sel.h"
#include "group.h"
#include "particle.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Auxiliary methods:
//

// ====================== operator<< (eSwaTopType) =========================

ostream& operator<<(ostream &out ,const eSwaTopType &type)
{
  switch(type)
  {
    case GBEST_TOPOLOGY:
      out << "'Gbest'";
    break;

    case RING_TOPOLOGY:
      out << "'Ring'";
    break;

    case SQUARE_TOPOLOGY:
      out << "'Square'";
    break;
  }

  return out;
}

// -------------------------------------------------------------------------
// Public methods:
//

// ============================== ReadSwarmTopology ========================

void cStandardPSO :: ReadSwarmTopology(istream &in)
{
  // Read the algorithm label.
  
  char label[100];

  if (!Utl::ReadString(in, label))
    Utl::Exit("Error in the input of the algorithm label.");

  if (string(label) == "Ring")
    Topology = RING_TOPOLOGY;
  else if (string(label) == "Gbest")
    Topology = GBEST_TOPOLOGY;
  else if (string(label) == "Square")
    Topology = SQUARE_TOPOLOGY;
  else
    Utl::Exit("Unknown Swarm Topology type: " + string(label));
}

// ============================== ReadInertia ==============================

void cStandardPSO :: ReadInertia(istream &in)
{
  if (!(in >> Inertia)) 
  {
    cout << "Error in the input of the inertia." << endl;
    exit(0);
  }
}

// ============================== ReadC1 ===================================

void cStandardPSO :: ReadC1(istream &in)
{
  if (!(in >> c1)) 
  {
    cout << "Error in the input of the cognitive factor." << endl;
    exit(0);
  }
}

// ============================== ReadC2 ===================================

void cStandardPSO :: ReadC2(istream &in)
{
  if (!(in >> c2)) 
  {
    cout << "Error in the input of the social factor." << endl;
    exit(0);
  }
}

// =============================== LoadReadFunc ============================

void cStandardPSO :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.
  im.Insert("PARTICLE.INERTIA", makeReadObj(cStandardPSO,ReadInertia));
  im.Insert("COGNITIVE.FACTOR", makeReadObj(cStandardPSO,ReadC1));
  im.Insert("SOCIAL.FACTOR",    makeReadObj(cStandardPSO,ReadC2));
  im.Insert("SWARM.TOPOLOGY",   makeReadObj(cStandardPSO,ReadSwarmTopology));
}

// ============================== cStandardPSO =============================

cStandardPSO :: cStandardPSO(void)
: Inertia(CONSTANT,0.72),c1(CONSTANT,1.50),c2(CONSTANT,1.50)
{
  Type       = STANDARD_PSO;
  Prob       = 0;
  Pen        = 0;
  SwarmSize  = 0;
  Topology   = GBEST_TOPOLOGY;
  NbhLenType = CANONICAL; 
  RedFactor  = -0.50;
}

// ============================= ~cStandardPSO ==============================

cStandardPSO :: ~cStandardPSO(void)
{
}

// =============================== Solver ==================================

void cStandardPSO :: Solver(void)
{
    // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Track number of individual evaluations.

    int EvalNum = 0;

    // Track the best objective function.

    double lastBest = 0.0;

    // Randomize rates.
    if (!MutProb && MaxMut)
      MutProb = Utl::RandDouble(MinMut, MaxMut);

    // Correct swarm size in the case of square topology.

    if (!SwarmSize)
      SwarmSize = CorrectSwarmSize(PopSize);

    // Create the swarm.

    cSwarm Swarm(SwarmSize, SolType, Prob);

    // Evaluate initial sample points.
    
    vector<cVector> sx;
    int sizsamp = SwarmSize - NumInpSol;
    if (sizsamp < 1) IntPopSamp = 0;

    if (IntPopSamp)
    {
      sx.reserve(sizsamp);
      int nvar = Prob->VarNumRow( ) * Prob->VarNumCol( );

      cSamp InitialSample;
      InitialSample.InitSample(SampType,nvar,sizsamp,sx);
    }

    // Generate initial swarm.

    for (int i = 0; i < SwarmSize; i++)
    {
      // Initialize particles position and velocity.
      if (i < NumInpSol)
        Swarm[i]->Init(InpSolVec[i]); // Input values.
      else if (IntPopSamp)  
        Swarm[i]->Init(sx[i]);        // Random values.
      else 
        Swarm[i]->Init( );            // Random values.
    }

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < SwarmSize; i++)
    {
      Swarm[i]->Evaluate( );
      #pragma omp critical
      EvalNum++;
    }

    if (Pen)
      Pen->EvalPenObjFunc(&Swarm, TolViol);

    // Update the best position.

    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < SwarmSize; i++)  // Update Best position without
      Swarm[i]->UpdateBestPos(true);     // checking condictions.       

    // Initialize the algorithm parameters.
  
    InitParameters( );
  
    if (Feedback) cout << "Optimization: " << opt + 1 << endl;

    // Perform the PSO iterations.

    for (int gen = 0; gen < MaxGen; gen++)
    {
      if ((gen+1)%5 == 0 && Feedback) cout << "Generation: " << gen + 1 << endl;

      // Evaluate Algorithm Parameters

      EvalParameters(gen);
      
      // Evaluate Velocity

      EvalVelocity(Swarm);

      // Mutation.

      Mutation(Swarm);
      
      // Update particles position (x = x + v)

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int i = 0; i < SwarmSize; i++)
        Swarm[i]->UpdatePos( );

      // Evaluate objective function.

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int j = 0; j < SwarmSize; j++)
      {
        Swarm[j]->Evaluate( );
	
        #pragma omp critical
	EvalNum++;
      }

      // Evaluate the penalized objective function.

      if (Pen)
        Pen->EvalPenObjFunc(&Swarm, TolViol);

      // Update particle best values (Xp).

      #pragma omp parallel for num_threads(omp_maxthread)
      for (int j = 0; j < SwarmSize; j++) 
        Swarm[j]->UpdateBestPos( );                

      // Migration.

      #ifdef _MPI_
      if (!((gen+1)%MigrationGap) && (gen+1) != MaxGen)
      {
        Migration(Swarm);
          
        if (Pen)
        Pen->EvalPenObjFunc(&Swarm, TolViol); 
      }
      #endif

      // Update variables related to PostProcessing.

      UpdatePostVar(gen, opt,lastBest, &Swarm);

      // Check conditions to stop optimization

      if (OptStopCrit(gen,opt,lastBest,&Swarm))
        break;
    }
   
    // Store the best individual.
    best->Insert(Swarm.BestSol( ));

    // Print data in the output file.
    PrintPostVar(MaxGen,opt,EvalNum,&Swarm);
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================ CorrectSwarmSize ===========================

int cStandardPSO :: CorrectSwarmSize(int swarmsize)
{
  int size = swarmsize;

  if (Topology == SQUARE_TOPOLOGY)
  {
    // found the best matrix format
    
    int i=1,j, last_j, diff, lastdiff, diff_ij, lastdiff_ij;

    lastdiff_ij = abs(1-swarmsize);
    last_j = swarmsize;   

    while (1)
    {
      i++;
      j = 1;
      lastdiff = abs(i*j - swarmsize);

      while(1)
      {
        j++;
        diff = abs(i*j - swarmsize);
        if (diff >= lastdiff)
        {
          j--;
          break;
        }
           
        lastdiff = diff;
      }
           
      diff_ij = abs(i-j);

      if (diff_ij >= lastdiff_ij)
      {
        i--;
        j = last_j;
        break;
      }
 
      lastdiff_ij = diff_ij;
      last_j = j;
    }

    int newsize = i*j;
    
    if (newsize != swarmsize)
    {
      cout << "The swarm size is set to " << newsize; 
      cout << " because square topology! (" << i << "," << j << ")\n";
    }
    else 
      cout << " Square Topology Matrix [" << i << "][" << j << "]\n";
    
    size = newsize;
 
    SquareMatSize.i = i;
    SquareMatSize.j = j;
    
  }

  return size;
}

// ============================ EvalBestNeighboor ==========================

cParticle* cStandardPSO :: EvalBestNeighboor(int id, cSwarm &swarm)
{ 
  cOptSolution *local_best = 0;

  if (Topology == GBEST_TOPOLOGY)
    local_best = swarm.BestSol( );
  else
  {
    vector<cOptSolution*> local;
    
    // Get particle's neighborhood indices.
    
    vector<int> NbhInds;
    GetNeighborhood(id,swarm,NbhInds);
    
    for(unsigned int i = 0; i < NbhInds.size( ); i++)
      local.push_back(swarm[NbhInds[i]]);
  
    local_best = *min_element(local.begin( ),local.end( ),cOptSolution::Compare); 
  }

  cParticle *part = dynamic_cast<cParticle*>(local_best);

  if (part)
    return part; 
  else
    Utl :: Exit("Invalid downcast on cStandardPSO :: EvalBEstNeighboor");
    
  return 0;
}

// ============================ GetNeighborhood ============================

void cStandardPSO :: GetNeighborhood(int id, cSwarm &swarm, vector<int> &Inds)
{
  // Required variables.
  
  int i,j;                              // Particle matrix index   
  int endId, nrow, ncol, i_abv, i_bel, j_lef, j_rig;

  switch(Topology)
  {  
    case SQUARE_TOPOLOGY:
      endId = swarm.GetSize( ) - 1;
      
      // Get row and column size.
      
      nrow = SquareMatSize.i;
      ncol = SquareMatSize.j;
         
      // Get particle matrix index.
      
      j = id%ncol;
      i = (id-j)/ncol;
 
      // Get neighbors index

      if (i-1 < 0)      // Above
        i_abv = nrow-1;
      else
        i_abv = i-1; 
      
      if (i+1 > nrow-1) // Below
        i_bel = 0;
      else i_bel = i+1;

      if (j-1 < 0)
        j_lef = ncol-1; // Left
      else
        j_lef = j-1; 
      
      if (j+1 > ncol-1) // Right
        j_rig = 0;
      else j_rig = j+1;
     
      Inds.push_back(id);
      Inds.push_back(i_abv*ncol + j);
      Inds.push_back(i_bel*ncol + j);
      Inds.push_back(i*ncol + j_lef);
      Inds.push_back(i*ncol + j_rig);
    break;

    case RING_TOPOLOGY:
    
      // Left, center and right particles.
      
      endId = swarm.GetSize( ) - 1;

      if(id > 0 && id < endId)
      { 
        Inds.push_back(id-1);
        Inds.push_back(id);
        Inds.push_back(id+1);
      }
      else
      {
        Inds.push_back(0);
        Inds.push_back(endId);
        
        if(id == 0)
          Inds.push_back(1);
        else
          Inds.push_back(endId-1);
      }

    break;
  
    case GBEST_TOPOLOGY:
      for(int k = 0; k < swarm.GetSize( ); k++)
        Inds.push_back(k);
    break;
  }
}

// ============================ InitParameters =============================

void cStandardPSO :: InitParameters()
{
  Inertia.setLimits(0,MaxGen); 
  c1.setLimits(0,MaxGen); 
  c2.setLimits(0,MaxGen); 
}

// ============================ EvalParameters =============================

void cStandardPSO :: EvalParameters(double time)
{
  Inertia.Evaluate(time);
  c1.Evaluate(time);
  c2.Evaluate(time);
}

// ============================ EvalVelocity ===============================

void cStandardPSO :: EvalVelocity(cSwarm &swarm)
{
  // Compute the velocity for each particle.

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < swarm.GetSize( ); i++)
  {
    if (NbhLenType == CANONICAL)
    {
      double r1 = Utl::RandDec( );
      double r2 = Utl::RandDec( );
      
      cParticle *Xg = EvalBestNeighboor(i,swarm);

      double Coef[] = {c1.GetValue( )*r1, c2.GetValue( )*r2};

      cParticle* Parts[] = {swarm[i],Xg};

      swarm[i]->EvalVelocity(Inertia.GetValue( ),2,Coef,Parts);
 
    }
    else if (NbhLenType == FULLY_INFORMED)
    {
      double r1 = Utl::RandDec( );

      vector<int> NbhInds;
      GetNeighborhood(i,swarm,NbhInds);

      int size = NbhInds.size( ); // Neighborhood size

      double *Coef = new double [size];
      cParticle** Parts = new cParticle* [size];

      for(int j = 0; j < size; j++)
      {
        Coef[j]  = (c1.GetValue( )*r1*(1.0/size));
        Parts[j] = swarm[NbhInds[j]];
      }

      swarm[i]->EvalVelocity(Inertia.GetValue( ),size,Coef,Parts);

      // Delete dynamics variables
      
      delete [] Coef;
      delete [] Parts;
    }
  }
}

// ============================== Mutation =================================

void cStandardPSO :: Mutation(cSwarm &swarm)
{
  // Perform the mutation operation.

  if (MutProb)
  { 
    #pragma omp parallel for num_threads(omp_maxthread)
    for (int i = 0; i < swarm.GetSize( ); i++) 
      swarm[i]->Mutate(MutProb);
  }
}

// ============================= Migration =================================

void cStandardPSO :: Migration(cSwarm &swarm)
{
#ifdef _MPI_
  
  // Get deme size.

  int popsize = swarm.GetSize( );
  
  // Get size and rank.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );

  // Check population size constraint.

  if ((size-1)*MigrationNum >= xp.GetSize( ))
  {
    cout << "Number of received individuals in one deme is greater than the original population size." << endl;
    exit(0);
  }

  // Send and receive individuals.
  
  for (int i = 0; i < MigrationNum; i++)
  {
    // Get best and worst particle ids from swarm.

    int BestPartId = xp.GetSolId(xp.BestSol( ));
    int WorstPartId = xp.GetSolId(xp.WorstSol( ));

    if (BestPartId ==-1) exit(0);  

    // Make the Migration process on the best particles position.

    for (int deme = 0; deme < size; deme++)
    {
      if (rank == deme)
        swarm[BestPartId]->Send( );
      else
        swarm[WorstPartId]->Receive(deme);
    }
  }

#endif
}

// ======================================================= End of file =====
