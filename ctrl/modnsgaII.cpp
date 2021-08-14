// -------------------------------------------------------------------------
// modNSGAII.cpp - implementation of class cmodNSGAII.
// -------------------------------------------------------------------------
// Copyright (c) 2013 LMCV/UFC
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
// Created:      30-Aug-2017    Marina Alves Maia
// -------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#ifdef _OMP_
#include "omp.h"
#endif

#ifdef _MPI_
#include "mpi.h"
#endif

#include "modnsgaII.h"
#include "problem.h"
#include "sel.h"
#include "group.h"
#include "individual.h"
#include "penalty.h"
#include "utl.h"
#include "input.h"
#include "gblvar.h"
#include "gbldef.h"

// -------------------------------------------------------------------------
// Class NSGAII-m:
//

// -------------------------------------------------------------------------
// Set read functions labels:

// =============================== LoadReadFunc ============================

void cmodNSGAII :: LoadReadFunc(cInpMap &im)
{
  // Call parent class load functions.
  cOptAlgorithm :: LoadReadFunc(im);

  // Register read functions.

  im.Insert("EXPLOREBOUNDS.RATE", makeReadObj(cmodNSGAII,ReadExploreBounds));
}


// ========================= ReadExploreBounds =========================

void cmodNSGAII :: ReadExploreBounds(istream &in)
{
  if (!(in >> ExpBoundsRate))
  {
    cout << "Error in the input of the explore boundaries rate." << endl;
    exit(0);
  }
}

// ============================== cmodNSGAII ==============================

cmodNSGAII :: cmodNSGAII(void) : cOptAlgorithm( )
{
  Type = modNSGAII;
  Prob = 0;
  ExpBoundsRate = 0;
}

// ============================= ~cmodNSGAII ==============================

cmodNSGAII :: ~cmodNSGAII(void)
{
}

// =============================== Solver ==================================

void cmodNSGAII :: Solver(void)
{
  // Solve the problem as many times as specified by the user.

  for (int opt = 0; opt < OptNum; opt++)
  {
    // Track number of individual evaluations.

    int EvalNum = 0;

    // Randomize rates.

    RandomRates( );

    // Create the population, mating pool and parent array.

    int newpopsize = 2*PopSize;
    cPopulation pop(PopSize, SolType, Prob);
    cPopulation son(PopSize, SolType, Prob);
    cPopulation parent(PopSize, SolType, Prob);
    cPopulation newpop(newpopsize, SolType, Prob);          // parent and son population merged together

    // Generate the initial population.

    #pragma omp parallel for num_threads(omp_maxthread)

    for (int i = 0; i < PopSize; i++)
    {
        // Fill genotypes with random values.

        pop[i]->Init( );

        // Evaluate the objective functions and constraints of the initial population

        pop[i]->Evaluate( );

     #pragma omp critical

        EvalNum++;
    }

     // Evaluate the normalization of all constraints for each individual.

      if (Pen)
        Pen->EvalPenObjFunc(&pop, TolViol);

      // Assign the rank number for each individual

       Rank(&pop);

       // Divide population in fronts according to the rank solution

       int numfronts;

       Fronts(pop, numfronts);

       CrowdDistance(pop, numfronts);

       // Select parents for crossover.

       Sel->Select(1.0, &pop, &parent, 0, &SolRank, &CrowdDistanceValue);

       // Crossover.

       Crossover(parent, son);

       Mutation(son);

       ExploreBounds(0, son, round(ExpBoundsRate*son.GetSize( )));

       // Evaluate the objective function and constraints of offsprings.

       // Main Loop starts here

       if (Feedback) cout << "Optimization: " << opt + 1 << endl;

       for (int gen = 0; gen < MaxGen; gen++)
       {
           if (((gen+1)%5 == 0 && Feedback) || gen == 0) cout << "Generation: " << gen + 1 << endl;

           for (int i = 0; i < PopSize; i++)
           {
             son[i]->Evaluate( );
             EvalNum++;
           }

           if (Pen)
             Pen->EvalPenObjFunc(&son, TolViol);

           // Merge parent and son into a new population

           Merge(pop, son, newpop);

           // Evaluate the normalization of all constraints for each individual.

           for (int p = 0; p < 2*PopSize; p++)
           {
               newpop[p]->Evaluate( );
           }

           if (Pen)
               Pen->EvalPenObjFunc(&newpop, TolViol);

           // Assign the rank number for each individual

             Rank(&newpop);

           // Sort population in Fronts using rank

            Fronts(newpop, numfronts);

            // Fill next generation population

            FillNextGen(pop, newpop);

           for (int p = 0; p < PopSize; p++)
           {
               pop[p]->Evaluate( );
           }

           // Evaluate the normalization of all constraints for each individual.

           if (Pen)
             Pen->EvalPenObjFunc(&pop, TolViol);

           // Rank

           Rank(&pop);

           // Fronts

           int allfronts;

           Fronts(pop, allfronts);

           CrowdDistance(pop, allfronts);

           // Select parents for crossover.

           Sel->Select(1.0, &pop, &parent, 0, &SolRank, &CrowdDistanceValue);

           // Crossover.

           Crossover(parent, son);

           Mutation(son);

           ExploreBounds(gen, son, round(ExpBoundsRate*son.GetSize()));

       }
       PrintPostVar(MaxGen,opt,EvalNum,&pop, &SolRank);
  }
}

// -------------------------------------------------------------------------
// Protected methods:
//

// ============================ Crossover ==================================

void cmodNSGAII :: Crossover(cPopulation &parent, cPopulation &son)
{
  // Perform the crossover operation.

  #pragma omp parallel for num_threads(omp_maxthread)
  for (int i = 0; i < son.GetSize( ); i = i+2)
  {
    double r = Utl::RandDec( );

    son[i  ]->Crossover(CrossType, r, parent[i],parent[i+1]);
    son[i+1]->Crossover(CrossType, r, parent[i+1],parent[i]);
  }
}

// ============================ ExploreBounds ==================================

void cmodNSGAII :: ExploreBounds(int gen, cPopulation &pop, int nbounds)
{

  int popsize = pop.GetSize( );

 // Replace nbounds individuals with their lower bound material variables

  #pragma omp parallel for num_threads(omp_maxthread)

  for (int i = 0; i < nbounds; i++)
  {
    double coord = Utl::RandDec( );
    int value = round((popsize-1)*coord);

  //  cout << "Modified solution " << value << endl;

    pop[value]->ExploreBounds( );

    pop[value]->Evaluate( );
  }

 // Replace nbounds individuals with their upper bound material variables

  for (int i = 0; i < nbounds; i++)
  {
    double coord = Utl::RandDec( );
    int value = round((popsize-1)*coord);

    // cout << "Modified solution " << value << endl;

    pop[value]->ExploreUpperBounds( );

    pop[value]->Evaluate( );
  }
}

// ============================== Mutation =================================

void cmodNSGAII :: Mutation(cPopulation &son)
{
  // Perform the mutation operation.
  if (MutProb)
  {
      for (int i = 0; i < son.GetSize( ); i++)
      {
          son[i]->Mutate(MutProb);
      }
   }
}

// ============================ checkdominance ==================================

int cmodNSGAII :: checkdominance(cOptSolution *a, cOptSolution *b)
{
    // 1 means a dominates b
    // 2 means b dominates a
    // 3 means there is no dominance relatioship between the solutions

    if (fabs(a->GetNormConst() - 0.0) < TolViol && fabs(b->GetNormConst() - 0.0) < TolViol)                        // Feasible individual
    {

          if ((fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) < TolViol) && a->GetObjFunc(1) < b->GetObjFunc(1)
                  && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) > TolViol)
          {              return 1;          }
          else if (a->GetObjFunc(0) < b->GetObjFunc(0) && (fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) < TolViol)
                   && fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) > TolViol)
          {              return 1;          }
          else if (a->GetObjFunc(0) < b->GetObjFunc(0) && a->GetObjFunc(1) < b->GetObjFunc(1)
                  && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) > TolViol && fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) > TolViol)
          {              return 1;          }
          else if ((fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) < TolViol) && a->GetObjFunc(1) > b->GetObjFunc(1) &&
                   fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) > TolViol)
          {              return 2;          }
          else if (a->GetObjFunc(0) > b->GetObjFunc(0) && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) < TolViol
                   && fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) > TolViol)
          {              return 2;          }
          else if (a->GetObjFunc(0) > b->GetObjFunc(0) && a->GetObjFunc(1) > b->GetObjFunc(1)
                  && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) > TolViol && fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) > TolViol)
          {              return 2;          }
          else if ((fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) < TolViol) && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) < TolViol)
          {              return 3;          }
          else if (a->GetObjFunc(0) < b->GetObjFunc(0) && a->GetObjFunc(1) > b->GetObjFunc(1)
                   && fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) > TolViol && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) > TolViol)
          {             return 3;          }
          else if (a->GetObjFunc(0) > b->GetObjFunc(0) && a->GetObjFunc(1) < b->GetObjFunc(1)
                   && fabs(a->GetObjFunc(0) - b->GetObjFunc(0)) > TolViol && fabs(a->GetObjFunc(1) - b->GetObjFunc(1)) > TolViol)
          {             return 3;          }
    }
    else if (fabs(a->GetNormConst() - 0.0) > TolViol && fabs(b->GetNormConst() - 0.0) < TolViol)               // Unfeasible individual
    {              return 2;    }
    else if (fabs(a->GetNormConst() - 0.0) < TolViol && fabs(b->GetNormConst() - 0.0) > TolViol)               // Feasible individual
    {              return 1;    }
    else if (fabs(a->GetNormConst() - 0.0) > TolViol && fabs(b->GetNormConst() - 0.0) > TolViol)
    {
        if (a->GetNormConst() < b->GetNormConst() && fabs(a->GetNormConst() - b->GetNormConst()) > TolViol)    // Unfeasible individual
        {              return 1;          }
        else if (a->GetNormConst() > b->GetNormConst() && fabs(a->GetNormConst() - b->GetNormConst()) > TolViol)
        {              return 2;          }
        else if (fabs(a->GetNormConst() - b->GetNormConst()) < TolViol)
        {             return 3;           }
    }
}

// ============================== Rank =====================================

void cmodNSGAII :: Rank(cGroup *pop)
{
    int popsize = pop->GetSize();

    SolRank.resize(popsize);

    // Evaluate rank number

    for (int i = 0; i < popsize; i++)
    {
        int ranktemp = 0;

     //   cout << "\nIndividual " << i << endl;

        for (int j = 0; j < popsize; j++)
        {
            if (i != j)
            {
                int aux = 0;

                aux = checkdominance(pop->GetSol(i), pop->GetSol(j));

                if (aux == 1)
                {
                    ranktemp = ranktemp + 0;
                }
                else if (aux == 2)
                {
                    ranktemp = ranktemp + 1;
            //        cout << "Solution " << j << " dominates it" << endl;
                }
                else if (aux == 3)
                {
                    ranktemp = ranktemp + 0;
                }
            }
        }

        SolRank[i] = ranktemp;
    //    cout << "Rank of individual " << i << ": " << SolRank[i] << "\n" << endl;

     }
}

// ============================= Fronts ====================================

void cmodNSGAII :: Fronts(cPopulation &population, int &numfronts)
{
    int populationsize = population.GetSize();

    int contadorrepetidos = 0;
    int uniqueranks = 0;

    // Evaluates the number of fronts

    int *RanksUnicos = new int [populationsize];
    int *RanksRepetidos = new int [populationsize];

    for (int i = 0; i < populationsize; i++)
    {
        int rep = 0;
        for (int j = 0; j < populationsize; j++)
        {
            if (SolRank[i] == SolRank[j] && i != j)
            {
                rep = rep + 1;
            }
            else rep = rep;
        }

        if (rep == 0)
        {
            RanksUnicos[uniqueranks] = SolRank[i];
            uniqueranks = uniqueranks + 1;
        }
        else if (rep != 0)
        {
            if (contadorrepetidos == 0)
            {
                RanksRepetidos[0] = SolRank[i];
                contadorrepetidos += 1;
            }

            int aux = 0;

            for (int k = 0; k < contadorrepetidos; k++)
            {
              if (SolRank[i] == RanksRepetidos[k])
              {
                    aux = aux + 1;
              }
              else aux = aux;
            }

            if (aux == 0 && contadorrepetidos != 0)
            {
                RanksRepetidos[contadorrepetidos] = SolRank[i];
                contadorrepetidos = contadorrepetidos + 1;
            }

        }
    }

    numfronts = uniqueranks + contadorrepetidos;

    int *RankGeral = new int [numfronts];

   // cout << "\nRanks Repetidos" << endl;
    for (int i = 0; i < contadorrepetidos; i ++)
    {
        RankGeral[i] = RanksRepetidos[i];
        // cout << RankGeral[i] << endl;
    }

  //  cout << "\nRanks Unicos" << endl;
    for (int i = 0; i < uniqueranks; i++)
    {
        RankGeral[contadorrepetidos+i] = RanksUnicos[i];
     //   cout << RankGeral[contadorrepetidos+i] << endl;
    }

 //   cout << "\nVetor de ranks por inteiro nao ordenado" << endl;
/*
    for (int i = 0; i < numfronts; i++)
    {
        cout << RankGeral[i] << endl;
    }*/

    // Ordenar RankGeral em ordem crescente

    for (int i = 0; i < numfronts; i++)
    {
        for (int j = i+1; j < numfronts; j++)
        {
            if (RankGeral[i] > RankGeral[j])
            {
                int temp = RankGeral[i];
                RankGeral[i] = RankGeral[j];
                RankGeral[j] = temp;
            }
        }
    }

    // Creates a struct for each front

    FrontClass = new sFrontClass[numfronts];

    int cont = 0;

    // Evaluates the size of each front

    for (int j = 0; j < numfronts; j++)
    {
        int cont = 0;
        for (int i = 0; i < populationsize; i++)
        {
            if (SolRank[i] == RankGeral[j])
            {
                cont += 1;
            }
        }
        FrontClass[j].FrontSize = cont;
    }

    for (int i = 0; i < numfronts; i++)
    {
        FrontClass[i].Front.resize(FrontClass[i].FrontSize);
    }

    // Filling in the Fronts

    for (int j = 0; j < numfronts; j++)
    {
        int cont = 0;
        for (int i = 0; i < populationsize; i++)
        {
            if (SolRank[i] == RankGeral[j])
            {
                FrontClass[j].Front[cont] = i;
                cont += 1;
            }
        }
    }

    // Debugging purposes
/*
   for (int i = 0; i < numfronts; i++)
    {
        cout << "\nMembros da frente " << i << endl;

        int tempfrontsize = FrontClass[i].FrontSize;

        for (int j = 0; j < tempfrontsize; j++)
        {
            cout << FrontClass[i].Front[j] << endl;
        }
    }
*/
}


// ============================= FillNextGen ===============================

void cmodNSGAII :: FillNextGen(cPopulation &pop, cPopulation &newpop)
{
    int popsize = pop.GetSize();

    int maxsize = 0;
    int contador = 0;
    int contfrentes = 0;

    maxsize = maxsize + FrontClass[contfrentes].FrontSize;

    while (maxsize <= popsize)
    {
        int tempfrontsize = FrontClass[contfrentes].FrontSize;
        for (int j = 0; j < tempfrontsize; j++)
        {pop[contador]->Copy(newpop[FrontClass[contfrentes].Front[j]]);
             contador += 1;
        }
        contfrentes += 1;
        maxsize = maxsize + FrontClass[contfrentes].FrontSize;
    }

  //  cout << "\nCurrent size of population " << contador << endl;

    int frontid = contfrentes;           // id of the front which surpassess the size of the next generation's population

   // cout << "\nFront " << frontid << " exceed max size of population of " << popsize << endl;

    if (contador < popsize)               // only access the algorithm if the next generation's population has exceeded its size
        {
             CrowdDistanceCut(newpop, frontid);

            // Sort crowding distance by descending order

            int frontsize = FrontClass[frontid].FrontSize;

            int *SortedCrowdDist = new int [frontsize];

            for (int l = 0; l < frontsize; l++)
            {
                SortedCrowdDist[l] = FrontClass[frontid].Front[l];
            }

            for (int i = 0; i < frontsize; i++)
            {
                for (int j = i+1; j < frontsize; j++)
                {
                    if (FrontClass[frontid].CrowdDistance[i] < FrontClass[frontid].CrowdDistance[j])
                    {
                        int temp = SortedCrowdDist[i];
                        SortedCrowdDist[i] = SortedCrowdDist[j];
                        SortedCrowdDist[j] = temp;
                    }
                }
            }

       /*   cout << "Ordered by CD" << endl;

         // Imprimir o vetor de solucoes CD ordenado

          cout << "\nImprimir o vetor de solucoes CD ordenado " << endl;

          for (int i = 0; i < frontsize; i++)
          {
              cout << SortedCrowdDist[i] << endl;
          }*/

          // Fill next generation with population sorted by CD

      //  cout << "\nTamanho da pop antes de add solucoes com CD calculado: " << contador << endl;

          int popsizerestante = popsize - contador;

          for (int i = 0; i < popsizerestante; i++)
          {
              pop[contador]->Copy(newpop[SortedCrowdDist[i]]);
              contador += 1;
          }

        //  cout << "Tamanho da pop apos a adicao das solucoes pelo CD: " << contador << endl;

          // Delete

          delete [] SortedCrowdDist;
        }
 }

// ============================ CrowdDistance =============================

void cmodNSGAII :: CrowdDistance(cPopulation &population, int numfronts)
{

  for (int frontid = 0; frontid < numfronts; frontid++)
  {
   int frontsize = FrontClass[frontid].FrontSize;

   int *FrontSol = new int [frontsize];

   // Passando o endereco das solucoes para um vetor auxiliar

    for (int i = 0; i < frontsize; i++)
    {FrontSol[i] = FrontClass[frontid].Front[i];}

   // cout << "\n======Frente " << frontid << " possui as solucoes: " << endl;
  /* if (frontsize > 2)
   {
       cout << "========================================" << endl;
   }

    for (int i = 0; i < frontsize; i++)
    {
        cout << FrontSol[i] << endl;
    }*/

    int *SortedSolObjFun1 = new int [frontsize];
    int *SortedSolObjFun2 = new int [frontsize];

    for (int i = 0; i < frontsize; i++)
    {
        SortedSolObjFun1[i] = FrontSol[i];
        SortedSolObjFun2[i] = FrontSol[i];
    }

    // Sort objective function 1 by ascending order

    for (int i = 0; i < frontsize; i++)
    {
        for (int j = i+1; j < frontsize; j++)
        {
            if (population.GetSol(SortedSolObjFun1[i])->GetObjFunc(0) > population.GetSol(SortedSolObjFun1[j])->GetObjFunc(0))
            {
                int temp = SortedSolObjFun1[i];
                SortedSolObjFun1[i] = SortedSolObjFun1[j];
                SortedSolObjFun1[j] = temp;
            }
         }
    }

    // Print
    /* cout << "\nVetor de funcoes objetivo 1 ordenado " << endl;
    for (int i = 0; i < frontsize; i++)
    {
        cout << "Solucao " << SortedSolObjFun1[i] << " tem fobj 1 = " << population.GetSol(SortedSolObjFun1[i])->GetObjFunc(0) << endl;
    }*/

    // Sort objective function 2 by ascending order


    for (int i = 0; i < frontsize; i++)
    {
        for (int j = i+1; j < frontsize; j++)
        {
            if (population.GetSol(SortedSolObjFun2[i])->GetObjFunc(1) > population.GetSol(SortedSolObjFun2[j])->GetObjFunc(1))
            {
                int temp = SortedSolObjFun2[i];
                SortedSolObjFun2[i] = SortedSolObjFun2[j];
                SortedSolObjFun2[j] = temp;
            }
         }
    }

    // Print

  /* cout << "\nVetor de funcoes objetivo 2 ordenado " << endl;
    for (int i = 0; i < frontsize; i++)
    {
        cout << "Solucao " << SortedSolObjFun2[i] << " tem fobj 2 = " << population.GetSol(SortedSolObjFun2[i])->GetObjFunc(1) << endl;
    }*/

    // Maximum and minimum values for objective functions

    double minimumfobj1 = population.GetSol(SortedSolObjFun1[0])->GetObjFunc(0);
    double maximumfobj1 = population.GetSol(SortedSolObjFun1[frontsize-1])->GetObjFunc(0);

    double minimumfobj2 = population.GetSol(SortedSolObjFun2[0])->GetObjFunc(1);
    double maximumfobj2 = population.GetSol(SortedSolObjFun2[frontsize-1])->GetObjFunc(1);

 /*   cout << "\nValores minimos e maximos das funcoes objetivo" << endl;
    cout << "Fobj 1 minimo: " << minimumfobj1 << endl;
    cout << "Fobj 1 maximo: " << maximumfobj1 << endl;
    cout << "Fobj 2 minimo: " << minimumfobj2 << endl;
    cout << "Fobj 2 maximo: " << maximumfobj2 << endl;*/

    double DifMaxmin1 = maximumfobj1 - minimumfobj1;
    double DifMaxmin2 = maximumfobj2 - minimumfobj2;

  /*  cout << "\n Diferenca" << endl;

    cout << "Dif maxima fobj1: " << DifMaxmin1 << endl;
    cout << "Dif maxima fobj2: " << DifMaxmin2 << endl;*/

    // Receber a distancia de multidao

    FrontClass[frontid].CrowdDistance.resize(frontsize);

    // Evaluate the CD

    if (frontsize == 1)
    {
        FrontClass[frontid].CrowdDistance[0] = 1.0e6;
    }
    else if (frontsize == 2)
    {
        FrontClass[frontid].CrowdDistance[0] = 1.0e6;
        FrontClass[frontid].CrowdDistance[1] = 1.0e6;
    }
    else if (frontsize > 2)
    {
        // Obj 1

        for (int j = 0; j < frontsize; j++)
        {
            for (int i = 0; i < frontsize; i++)
            {
                if (SortedSolObjFun1[j] == FrontClass[frontid].Front[i])
                {
                    if (j == 0 || j == frontsize-1)
                    {
                        FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                    }
                    else
                    {
                             if (DifMaxmin1 == 0)
                            {
                                FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                            }
                            else if (fabs(DifMaxmin1) != 0)
                            {
                                FrontClass[frontid].CrowdDistance[i] = (population.GetSol(SortedSolObjFun1[j+1])->GetObjFunc(0) -
                                        population.GetSol(SortedSolObjFun1[j-1])->GetObjFunc(0))/fabs(DifMaxmin1);
                            }
                        }
                    }
                }
            }


        // Obj 2

        for (int j = 0; j < frontsize; j++)
        {
            for (int i = 0; i < frontsize; i++)
            {
                if (SortedSolObjFun2[j] == FrontClass[frontid].Front[i])
                {
                    if (j == 0 || j == frontsize-1)
                    {
                        FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                    }
                    else
                    {
                            if (DifMaxmin2 == 0)
                            {
                                FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                            }
                            else if (fabs(DifMaxmin2) != 0)
                            {
                                FrontClass[frontid].CrowdDistance[i] = FrontClass[frontid].CrowdDistance[i] + (population.GetSol(SortedSolObjFun2[j+1])->GetObjFunc(1) -
                                        population.GetSol(SortedSolObjFun2[j-1])->GetObjFunc(1))/fabs(DifMaxmin2);
                            }
                        }
                    }
                }
            }
    }


  // Print teste

/* cout << "\nValor de crowd distance das solucoes" << endl;
  for (int i = 0; i < frontsize; i++)
  {
      cout << "Solucao " << FrontClass[frontid].Front[i] << " possui CD = " << FrontClass[frontid].CrowdDistance[i] << endl;
  }*/


  // Delete
  delete [] FrontSol;
  delete [] SortedSolObjFun1;
  delete [] SortedSolObjFun2;
  }

  // Passando os CDs da struct para o vetor

  int populationsize = population.GetSize();
  CrowdDistanceValue.resize(populationsize);

  for (int cont = 0; cont < populationsize; cont++)
  {
  for (int i = 0; i < numfronts; i++)
   {
       int tempfrontsize = FrontClass[i].FrontSize;

       for (int j = 0; j < tempfrontsize; j++)
       {
           if (FrontClass[i].Front[j] == cont)
           {
               CrowdDistanceValue[cont] = FrontClass[i].CrowdDistance[j];
           }
       }
   }
  }

  // Print para verificacao

 /* cout << "Verificando o vetor" << endl;

  for (int i = 0; i < populationsize; i++)
  {
      cout << "Solucao " << i << " possui CD = " << CrowdDistanceValue[i] << " e rank = " << SolRank[i] << endl;
  }*/

}

// ============================ CrowdDistanceCut =============================

void cmodNSGAII :: CrowdDistanceCut(cPopulation &population, int frontidd)
{

    for (int frontid = frontidd; frontid < frontidd+1; frontid++)
    {
     int frontsize = FrontClass[frontid].FrontSize;

     int *FrontSol = new int [frontsize];

     // Passando o endereco das solucoes para um vetor auxiliar

      for (int i = 0; i < frontsize; i++)
      {FrontSol[i] = FrontClass[frontid].Front[i];}

     // cout << "\n==Frente " << frontid << " possui as solucoes tal: " << endl;
    /* if (frontsize > 2)
     {
         cout << "========================================" << endl;
     }

      for (int i = 0; i < frontsize; i++)
      {
          cout << FrontSol[i] << endl;
      }*/

      int *SortedSolObjFun1 = new int [frontsize];
      int *SortedSolObjFun2 = new int [frontsize];

      for (int i = 0; i < frontsize; i++)
      {
          SortedSolObjFun1[i] = FrontSol[i];
          SortedSolObjFun2[i] = FrontSol[i];
      }

      // Sort objective function 1 by ascending order

      for (int i = 0; i < frontsize; i++)
      {
          for (int j = i+1; j < frontsize; j++)
          {
              if (population.GetSol(SortedSolObjFun1[i])->GetObjFunc(0) > population.GetSol(SortedSolObjFun1[j])->GetObjFunc(0))
              {
                  int temp = SortedSolObjFun1[i];
                  SortedSolObjFun1[i] = SortedSolObjFun1[j];
                  SortedSolObjFun1[j] = temp;
              }
           }
      }



      // Print
      /* cout << "\nVetor de funcoes objetivo 1 ordenado " << endl;
      for (int i = 0; i < frontsize; i++)
      {
          cout << "Solucao " << SortedSolObjFun1[i] << " tem fobj 1 = " << population.GetSol(SortedSolObjFun1[i])->GetObjFunc(0) << endl;
      }*/

      // Sort objective function 2 by ascending order


      for (int i = 0; i < frontsize; i++)
      {
          for (int j = i+1; j < frontsize; j++)
          {
              if (population.GetSol(SortedSolObjFun2[i])->GetObjFunc(1) > population.GetSol(SortedSolObjFun2[j])->GetObjFunc(1))
              {
                  int temp = SortedSolObjFun2[i];
                  SortedSolObjFun2[i] = SortedSolObjFun2[j];
                  SortedSolObjFun2[j] = temp;
              }
           }
      }

      // Print

    /* cout << "\nVetor de funcoes objetivo 2 ordenado " << endl;
      for (int i = 0; i < frontsize; i++)
      {
          cout << "Solucao " << SortedSolObjFun2[i] << " tem fobj 2 = " << population.GetSol(SortedSolObjFun2[i])->GetObjFunc(1) << endl;
      }*/

      // Maximum and minimum values for objective functions

      double minimumfobj1 = population.GetSol(SortedSolObjFun1[0])->GetObjFunc(0);
      double maximumfobj1 = population.GetSol(SortedSolObjFun1[frontsize-1])->GetObjFunc(0);

      double minimumfobj2 = population.GetSol(SortedSolObjFun2[0])->GetObjFunc(1);
      double maximumfobj2 = population.GetSol(SortedSolObjFun2[frontsize-1])->GetObjFunc(1);

   /*   cout << "\nValores minimos e maximos das funcoes objetivo" << endl;
      cout << "Fobj 1 minimo: " << minimumfobj1 << endl;
      cout << "Fobj 1 maximo: " << maximumfobj1 << endl;
      cout << "Fobj 2 minimo: " << minimumfobj2 << endl;
      cout << "Fobj 2 maximo: " << maximumfobj2 << endl;*/

      double DifMaxmin1 = maximumfobj1 - minimumfobj1;
      double DifMaxmin2 = maximumfobj2 - minimumfobj2;

    /*  cout << "\n Diferenca" << endl;

      cout << "Dif maxima fobj1: " << DifMaxmin1 << endl;
      cout << "Dif maxima fobj2: " << DifMaxmin2 << endl;*/

      // Receber a distancia de multidao

      FrontClass[frontid].CrowdDistance.resize(frontsize);

      // Evaluate the CD

      if (frontsize == 1)
      {
          FrontClass[frontid].CrowdDistance[0] = 1.0e6;
      }
      else if (frontsize == 2)
      {
          FrontClass[frontid].CrowdDistance[0] = 1.0e6;
          FrontClass[frontid].CrowdDistance[1] = 1.0e6;
      }
      else if (frontsize > 2)
      {
          // Obj 1

          for (int j = 0; j < frontsize; j++)
          {
              for (int i = 0; i < frontsize; i++)
              {
                  if (SortedSolObjFun1[j] == FrontClass[frontid].Front[i])
                  {
                      if (j == 0 || j == frontsize-1)
                      {
                          FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                      }
                      else
                      {
                               if (DifMaxmin1 == 0)
                              {
                                  FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                              }
                              else if (fabs(DifMaxmin1) != 0)
                              {
                                  FrontClass[frontid].CrowdDistance[i] = (population.GetSol(SortedSolObjFun1[j+1])->GetObjFunc(0) -
                                          population.GetSol(SortedSolObjFun1[j-1])->GetObjFunc(0))/fabs(DifMaxmin1);
                              }
                          }
                      }
                  }
              }


          // Obj 2

          for (int j = 0; j < frontsize; j++)
          {
              for (int i = 0; i < frontsize; i++)
              {
                  if (SortedSolObjFun2[j] == FrontClass[frontid].Front[i])
                  {
                      if (j == 0 || j == frontsize-1)
                      {
                          FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                      }
                      else
                      {
                              if (DifMaxmin2 == 0)
                              {
                                  FrontClass[frontid].CrowdDistance[i] = 1.0e6;
                              }
                              else if (fabs(DifMaxmin2) != 0)
                              {
                                  FrontClass[frontid].CrowdDistance[i] = FrontClass[frontid].CrowdDistance[i] + (population.GetSol(SortedSolObjFun2[j+1])->GetObjFunc(1) -
                                          population.GetSol(SortedSolObjFun2[j-1])->GetObjFunc(1))/fabs(DifMaxmin2);
                              }
                          }
                      }
                  }
              }
      }


    // Print teste

  /* cout << "\nValor de crowd distance das solucoes" << endl;
    for (int i = 0; i < frontsize; i++)
    {
        cout << "Solucao " << FrontClass[frontid].Front[i] << " possui CD = " << FrontClass[frontid].CrowdDistance[i] << endl;
    }*/


    // Delete
    delete [] FrontSol;
    delete [] SortedSolObjFun1;
    delete [] SortedSolObjFun2;
    }


    // Print para verificacao

   /* cout << "Verificando o vetor" << endl;

    for (int i = 0; i < populationsize; i++)
    {
        cout << "Solucao " << i << " possui CD = " << CrowdDistanceValue[i] << " e rank = " << SolRank[i] << endl;
    }*/
}

// ============================= Migration =================================

void cmodNSGAII :: Migration(cPopulation &pop)
{
#ifdef _MPI_

  // Get deme size.

  int popsize = pop.GetSize( );

  // Get size and rank.

  int rank = MPI::COMM_WORLD.Get_rank( );
  int size = MPI::COMM_WORLD.Get_size( );

  // Check population size constraint.

  if ((size-1)*MigrationNum >= pop.GetSize( ))
  {
    cout << "Number of received individuals in one deme is greater than the    original population size." << endl;
    exit(0);
  }

  // Send and receive individuals.

  for (int i = 0; i < MigrationNum; i++)
    for (int deme = 0; deme < size; deme++)
    {
      if (rank == deme)
        pop[i]->Send( );
      else
        pop[popsize-1-i]->Receive(deme);
    }

#endif
}

// ================================ Merge ==================================

void cmodNSGAII :: Merge(cPopulation &pop, cPopulation &son, cPopulation &newpop)
{
  // Get necessary data.

  int popsize = pop.GetSize( );
  int sonsize = son.GetSize( );

  #pragma omp parallel for num_threads(omp_maxthread)

 // cout << "Populacao de pais e filhos unidos de tamanho " << popsize + sonsize << "\n " << endl;

  for (int i = 0; i < popsize; i++)
  {
    newpop[i]->Copy(pop[i]);
  }

  for (int j = 0; j < sonsize; j++)
  {
    newpop[j+popsize]->Copy(son[j]);
  }
}

// ============================= RandomRates ==============================

void cmodNSGAII :: RandomRates(void)
{
    if (!MutProb && MaxMut)
      MutProb = Utl::RandDouble(MinMut, MaxMut);
}

// ======================================================= End of file =====
