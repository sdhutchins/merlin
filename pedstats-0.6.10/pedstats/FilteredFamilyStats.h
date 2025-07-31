////////////////////////////////////////////////////////////////////// 
// pedstats/FilteredFamilyStats.h 
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
 
#ifndef __FILTEREDFAMSTATS_H__
#define __FILTEREDFAMSTATS_H__

#include "IntArray.h"
#include "PedigreePerson.h"

class FilteredPedigree;

class FilteredFamilyStats
   {
   public:

      //Indicator array for families from original pedigree with all members filtered out
      IntArray filteredFamilies;

      //Indices in original pedigree of families retained after filtering
      IntArray familyIndices;

      // Number of members, generations per remaining family after filtering
      IntArray familyMembers, familyGenerations;
      // Number of non-founders, founders per remaining family after filtering
      IntArray familyNonFounders, familyFounders;
      // Number of females, males per remaining family after filtering
      IntArray familyFemales, familyMales;

      // Overall statistics after filtering
      int maxMembers, minMembers;
      int maxFounders, minFounders, founderCount;
      int maxNonFounders, minNonFounders, nonFounderCount;
      int maxFemales, minFemales, femaleCount;
      int maxMales, minMales, maleCount;
      int minGenerations, maxGenerations;

      double averageFounders, averageNonFounders;
      double averageFemales, averageMales;
      double averageGenerations, averageMembers;

      // Number of families, individuals after filtering
      int familiesKept, personCount;

       FilteredFamilyStats(FilteredPedigree & p);
      ~FilteredFamilyStats() {};

      void AccumulateStatistics();
      int  Generation(Person * p);

   private:

      FilteredPedigree & ped;

      void Initialize();

      void AccumulateFamilyCounts();
      void AccumulateOverallStatistics();

   };

#endif
 
 
 
 
 
 
