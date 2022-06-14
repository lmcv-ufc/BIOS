// -------------------------------------------------------------------------
// main.cpp - main driver of BIOS.
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
//               Building of problem, algorithm, selection method and
//               penalty function objects.
//
//               16-Jul-2014    Elias Sariava Barroso
//               Update Version label do 3.0.
//
//               13-Nov-2015    Elias Sariava Barroso
//               Created global variable to set the number of threads used in
//               the parallel (OpenMP) program version.
//
//               26-Out-2018    Leonardo Gonçalves Ribeiro
//               Consertado Evaluate( ) das funções para englobar o parâmetro gen()
//
// -------------------------------------------------------------------------

#include <time.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

#ifdef _MPI_
#include <mpi.h>
#endif

#ifdef _OMP_
#include "omp.h"
#endif

#include "input.h"
#include "problem.h"
#include "individual.h"
#include "optalg.h"
#include "sel.h"
#include "penalty.h"
#include "utl.h"

typedef long double tClock;
const string Version("4.0.0 - Feb/2022");

// -------------------------------------------------------------------------
// Global variables:
//
ifstream in;   // input file
string fname;  // input file name
int omp_maxthread = 1;
unsigned int seed;

// -------------------------------------------------------------------------
// main function:
//
int main(int argc, char **argv)
{
  // Get the maximum number of threads from OpenMP.

#ifdef _OMP_
  #pragma omp master
    omp_maxthread = omp_get_max_threads( );
#endif

  // Start MPI.

#ifdef _MPI_
  MPI::Init(argc, argv);
#endif

  // Set rand() seed.
  unsigned int seed;

#ifdef _MPI_
  seed = (unsigned)time(0)+1+MPI::COMM_WORLD.Get_rank( );
#else
  seed = (unsigned)time(0);
#endif
  Utl :: SetSeed(seed);

  // Get program options (feedback).
  bool Feedback = true;
  if (argc > 2 && string(argv[2]) == "-silent") Feedback = false;

  if (Feedback)
  {
    cout << endl;
    cout << "\t===================================================================" << endl;
    cout << "\t==                                                               ==" << endl;
    cout << "\t==        Biologically Inspired Optimization System - BIOS       ==" << endl;
    cout << "\t==                                                               ==" << endl;
    cout << "\t==  Laboratorio de Mecanica Computacional e Visualizacao - LMCV  ==" << endl;
    cout << "\t==                                                               ==" << endl;
    cout << "\t==            Universidade Federal do Ceara - UFC                ==" << endl;
    cout << "\t==                                                               ==" << endl;
    cout << "\t==                   Version " << left << setw(36) << Version <<"==" << endl;
    cout << "\t===================================================================" << endl << endl;
  }

  // Get the name of the input file.

  if (argc == 1)
  {
    cout << endl << "\tEnter the input file name [.opt]......... ";
    cin >> fname;
  }
  else
    fname = string(argv[1]);

  // Start time counter.

  tClock cputime = clock( );
  time_t start = time(NULL);

  // Open the input file.

  string aux = fname + ".opt";
  in.open(aux.c_str());
  if (!in.is_open())
  {
    cout << "Invalid input file." << endl;
    return(1);
  }

  // Open the output file.
  ofstream out;

  aux = fname + ".out";
  out.open(aux.c_str());
  if (!out.is_open())
  {
    cout << "Invalid output file." << endl;
    return(1);
  }

  // Create reader object to store input optimization algorithm.
  cInpMap CtrlMap;
  CtrlMap.SetFeedback(Feedback);
  cOptAlgReadEntry *algread = new cOptAlgReadEntry;
  algread->inpmap           = &CtrlMap;
  CtrlMap.Insert("OPTIMIZATION.ALGORITHM",algread);

  // Read the input data and close the input file.
  if (Feedback)
    cout << endl << "\tReading the input data .................." << endl;
  ReadNeutralFile(CtrlMap,in);         // CtrlMap is defined in input.cpp
  CopyFileData(in,out);
  in.close( );

  //  Read problem input file.
  cInpMap ProbMap;
  cOptAlgorithm *alg  = algread->alg;
  cProblem      *prob = alg->GetProblem( );
  if (!prob)
  {
    cout << "An optimization problem must be defined!" << endl;
    return(1);
  }
  prob->LoadReadFunc(ProbMap);
  string pfext = cProblemFactory :: GetProbFileExt(prob);

  if (!pfext.empty( ))
  {
    aux = fname + pfext;
    in.open(aux.c_str());
    cout << "fname " << aux << endl; 
    if (!pfext.empty( ) && out.is_open( ))
    {
      ReadNeutralFile(ProbMap,in);         
      CopyFileData(in,out);
    }
    in.close( );
  }
 
  // Process the optimization and print the output data.

  if (Feedback)
    cout << endl << "\tProcessing the Optimization Problem ..........." << endl;
  alg->SetOutStream(out);
  alg->SetFeedback(Feedback);
  alg->Init( );               //Method defined in OptAlg 
  alg->Solver( );             //Virtual method defined in OptAlg. This method calls Evaluate( )]
  alg->PostProcessing( );     //Defined in OptAlg

  // Write end mark and close the output file.

  // Destroy the model.
  delete alg;

  // Print the total elapsed time.

  time_t end = time(NULL);
  double diff = difftime(end,start);
  cputime = clock( ) - cputime;
  if (Feedback)
  {
    cout << "Seed: " << Utl :: GetSeed( ) << endl;
    cout << endl << "\t\tCPU time: " << (double)(cputime/CLOCKS_PER_SEC) << " (s)" << endl;
    cout << endl << "\t\tWall clock: " << diff << " (s)" << endl;
  }
  out << "\n%SEED\n" << seed << "\n";
  out << "\n%CPU.TIME" << endl;
  out << (double)(cputime/CLOCKS_PER_SEC) << endl;

  out << "\n%WALL.CLOCK.TIME" << endl;
  out << diff << endl;
  out  << "\n%END";

  // Finalize MPI

#ifdef _MPI_
  MPI::Finalize( );
#endif

  return(0);
}

// ======================================================= End of file =====
