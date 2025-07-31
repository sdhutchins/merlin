////////////////////////////////////////////////////////////////////// 
// pedstats/PedigreePairs.h 
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
 
// PedigreePairs.h
/* Written by Jan Wigginton */

#ifndef __PEDIGREEPAIRS_H__
#define __PEDIGREEPAIRS_H__

#include "SimplePairs.h"
#include "FilteredPedigree.h"

class PedigreePairs
   {
   public:

      FilteredPedigree & ped;

      SimplePairList sibs, halfSibs;
      SimplePairList cousins, avuncular;
      SimplePairList grandparents, parents;
      SimplePairList others;

      int sibCt, halfSibCt;
      int cousinCt, avuncularCt;
      int parentCt, grandParentCt;
      int otherCt;

      int totalCt;

      bool doOthers;
     // static double markerPctCutoff;

      // constructors
      PedigreePairs(FilteredPedigree & p);

      // utility
      void Build();
      void UpdateCount();
      void PrintCounts(FILE * file, const char * label = NULL);

      // split functions
      void SplitOnAffection(int affection, PedigreePairs & unaffected,
         PedigreePairs & discordant, PedigreePairs & affected);
      void SplitOnGender(PedigreePairs & females, PedigreePairs & males,
         PedigreePairs & opposite);
      void SplitOnGenderOrder(PedigreePairs & mf_pairs, PedigreePairs & fm_pairs);

      //filter functions
      void FilterOnCovariate(int covariate, PedigreePairs & filtered);
      void FilterOnTrait(int trait, PedigreePairs & filtered);
      void FilterOnDiagnosed(int affection, PedigreePairs & filtered);
      void FilterOnGenotyped(int marker, PedigreePairs & filtered);

      //count functions
      void CountAffectionTypes(int var, IntArray & unaffected, IntArray & discordant,
         IntArray & affected, bool count_unordered = true);
      void GetCorrelations(int var, int is_covariate, Vector & corr, IntArray & counts,
         bool do_unordered = true);
   };

#endif





 
 
 
 
