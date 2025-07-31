////////////////////////////////////////////////////////////////////// 
// pedstats/SimplePairs.h 
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
 
/* Written by Jan Wigginton */

#ifndef __SIMPLEPAIR_H__
#define __SIMPLEPAIR_H__

#include "PedigreePerson.h"
#include "FilteredPedigree.h"
#include "Constant.h"


#define TRAIT_TOL 1e-8

class SimplePair
   {
   public:

      // i1 is ancestor, i2 is descendant when determined by PedigreePairs::Build
      int i1, i2;

      void operator=(SimplePair & rhs) { i1 = rhs.i1; i2 = rhs.i2; }

      SimplePair(){};
      SimplePair(int i_1, int i_2) { i1 = i_1 ; i2 = i_2; }
   };


class SimplePairList
   {
   public:

      SimplePairList(FilteredPedigree & p);
      SimplePairList(FilteredPedigree & p, const char *x_tag, const char *y_tag);
      ~SimplePairList();

      FilteredPedigree   & ped;
      SimplePair * list;

      int          count;
      String xTag, yTag;
      String gender;

      void Clear();
      void Append(int i, int j);
      void Append(SimplePair & rhs);
      SimplePair & operator[](int i) { return list[i]; }

      // split and count functions
      void SplitOnAffection(int affection, SimplePairList & unaffected,
        SimplePairList & discordant, SimplePairList & affected);
      void SplitOnGender(SimplePairList & females, SimplePairList & males,
        SimplePairList & opposite);
      void SplitOnGenderOrder(SimplePairList & mf_pairs, SimplePairList & fm_pairs);
      void CountAffectionTypes(int var, int & unaffected, int  & discordant, int & affected);

      // pairwise trait values
      void CheckSameTraitValue(int var, const char * relationship_label, int digits,
            bool is_covariate, bool twins_only = true);
      void CheckValueGap(int var, double gap, const char * relationship_label, int digits,
            bool is_covariate, bool min_gap);

      // filter functions
      void FilterOnCovariate(int covariate, SimplePairList & filtered);
      void FilterOnTrait(int trait, SimplePairList & filtered);
      void FilterOnDiagnosed(int affection, SimplePairList & filtered);
      void FilterOnGenotyped(int marker, SimplePairList & filtered);
      void FilterOnBothFemale(SimplePairList & females);
      void FilterOnBothMale(SimplePairList & males);

      void   OrderPairsOnGender(SimplePairList & new_list);
      double Correlation(int var, bool is_covariate, bool mirrored_pairs,
         int & pairs_typed, Vector * x_values = NULL, Vector * y_values = NULL);

      // Builds a list of all sib pairs in pedigree
      void   BuildSibList();
   private:

      int size;

      void Grow();

      bool RetrievePairValues(int var, Person * p1, Person * p2, double & val_1, double & val_2, bool is_covariate);
      void UpdateTags(SimplePairList & pairs, const char * g_tag = NULL);
      double GetTraitTol(int var, bool is_covariate);
      void Initialize(const char * x_tag, const char * y_tag);
   };

bool isOrderedParent(Person * p1, Person * p2);
bool isParent(Person * p1, Person * p2);
bool isOrderedGrandParentOffspring(Person * p1, Person * p2);
bool isGrandParentOffspring(Person * p1, Person * p2);
bool isOrderedAvuncular(Person * p1, Person * p2);
bool isAvuncular(Person * p1, Person * p2);
bool isSib(Person * p1, Person * p2);
bool isHalfSib(Person * p1, Person * p2);
bool isCousin(Person * p1, Person * p2);
bool isOtherRelative(Person * p1, Person * p2);

bool FilterOnTrait(int var, Person & p1, Person & p2);
bool FilterOnCovariate(int var, Person & p1, Person & p2);
bool FilterOnDiagnosed(int affection, Person & p1, Person & p2);
bool FilterOnGenotype(int marker, Person & p1, Person & p2);
bool FilterOnBothMale(Person * p1, Person * p2);
bool FilterOnBothFemale(Person * p1, Person * p2);
bool FilterOnOppositeSex(Person * p1, Person * p2);

#endif
 
 
 
 
 
 
 
 
