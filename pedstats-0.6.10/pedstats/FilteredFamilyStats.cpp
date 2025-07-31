////////////////////////////////////////////////////////////////////// 
// pedstats/FilteredFamilyStats.cpp 
// (c) 2000-2007 Goncalo Abecasis (c) 2002-2007 Jan Wigginton
// 
// This file is distributed as part of the PEDSTATS source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile PEDSTATS.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#include "FilteredFamilyStats.h"
#include "FilteredPedigree.h"
#include "MathConstant.h"

// Written by Jan Wigginton

FilteredFamilyStats::FilteredFamilyStats(FilteredPedigree & p) : ped(p)
   {
   }

void FilteredFamilyStats::Initialize()
   {
   minFemales = minMales = minFounders = minNonFounders = 0;
   minMembers = minGenerations = 0;
   maxFemales = maxMales = maxFounders = maxNonFounders = 0;
   maxMembers = maxGenerations = 0;

   femaleCount = maleCount = founderCount = nonFounderCount = personCount = 0;
   averageGenerations = averageFemales = averageMales = 0;
   averageFounders = averageNonFounders = averageMembers = 0;

   familyMembers.Dimension(familiesKept);
   familyMembers.Zero();

   familyFounders.Dimension(familiesKept);
   familyFounders.Zero();

   familyNonFounders.Dimension(familiesKept);
   familyNonFounders.Zero();

   familyGenerations.Dimension(familiesKept);
   familyGenerations.Zero();

   familyFemales.Dimension(familiesKept);
   familyFemales.Zero();

   familyMales.Dimension(familiesKept);
   familyMales.Zero();
   }

void FilteredFamilyStats::AccumulateStatistics()
   {
   IntArray remaining_members(ped.familyCount);
   remaining_members.Zero();

   for (int f = 0; f < ped.familyCount; f++)
      for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         remaining_members[f] += (!ped.filtered[i]);

   filteredFamilies.Dimension(ped.familyCount);
   filteredFamilies.Zero();

   familiesKept = 0;
   for (int f = 0; f < ped.familyCount; f++)
      {
      filteredFamilies[f]  = (remaining_members[f] == 0);
      if (!filteredFamilies[f]) familiesKept++;
      }

   Initialize();
   AccumulateFamilyCounts();
   AccumulateOverallStatistics();
   }

void FilteredFamilyStats::AccumulateFamilyCounts()
   {
   if (familiesKept == 0) return;

   int kept_index = 0;
   for (int f = 0; f < ped.familyCount; f++)
      {
      if (filteredFamilies[f]) continue;

      familyIndices.Push(f);

      for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         {
         familyGenerations[kept_index] = max(Generation(ped.persons[i]) + 1,
                                   familyGenerations[kept_index]);
         if(ped.filtered[i]) continue;

         ped[i].isFounder() ? familyFounders[kept_index]++ : familyNonFounders[kept_index]++;

         familyFemales[kept_index] += (ped[i].sex == SEX_FEMALE);
         familyMales[kept_index]   += (ped[i].sex == SEX_MALE);
         familyMembers[kept_index] ++;
         }

      kept_index ++;
      }
   }

void FilteredFamilyStats::AccumulateOverallStatistics()
   {
   if (familiesKept == 0) return;

   maxFemales = familyFemales.Max();
   minFemales = familyFemales.Min();
   femaleCount = familyFemales.Sum();

   maxMales = familyMales.Max();
   minMales = familyMales.Min();
   maleCount = familyMales.Sum();

   maxFounders = familyFounders.Max();
   minFounders = familyFounders.Min();
   founderCount = familyFounders.Sum();

   maxNonFounders = familyNonFounders.Max();
   minNonFounders = familyNonFounders.Min();
   nonFounderCount = familyNonFounders.Sum();

   maxMembers = familyMembers.Max();
   minMembers = familyMembers.Min();
   personCount = familyMembers.Sum();

   maxGenerations = familyGenerations.Max();
   minGenerations = familyGenerations.Min();

   double div  = max(familiesKept, 1);
   averageGenerations = familyGenerations.Sum() / div;
   averageFemales = femaleCount / div;
   averageMales = maleCount / div;
   averageFounders = founderCount / div;
   averageNonFounders = nonFounderCount / div;
   averageMembers = personCount / div;
   }

int FilteredFamilyStats::Generation(Person * p)
   {
   if (p->isFounder())
      return 0;

   int father = Generation(p->father);
   int mother = Generation(p->mother);

   return (mother > father ? mother : father) + 1;
   }



 
 
 
 
 
 
