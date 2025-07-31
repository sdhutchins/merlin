////////////////////////////////////////////////////////////////////// 
// pedstats/Pedstats.h 
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
 
#ifndef __PEDSTATS_H__
#define __PEDSTATS_H__

#include "Constant.h"
#include "MathStats.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathConstant.h"
#include "Parameters.h"
#include "PedigreePairs.h"
#include "Error.h"
#include "IBD.h"
#include  "PDF.h"
#include "QuickIndex.h"
#include "Manners.h"
#include "FilteredPedigree.h"
#include "PedstatsQuality.h"
#include "PedstatsHWE.h"
#include "PedstatsAgeCheck.h"
#include "FilteredFamilyStats.h"
#include "GenotypeLists.h"

#include <stdio.h>

#define EXACT_MIN_ALLELES 10
#define HET_WINDOW_SIZE 10

void   ReadParameters(int argc, char * argv[]);
void   FamilyStatistics(FilteredFamilyStats & stats);
void   CheckFamilyConnections();
void   RewritePedigree();
void   RunInheritanceChecks();
void   FileReport(bool wrote_pdf);
void   CheckForLogicalSummaries(FilteredFamilyStats & stats);

void   TraitStatistics(FILE * file, FilteredPedigree & ped, SimplePairList & pairs,
         const char * add = NULL);
void   CovarStatistics(FILE * file, FilteredPedigree & ped, SimplePairList & pairs,
         const char * add = NULL);
void   GetTraitStatistics(FILE * file, FilteredPedigree & ped, SimplePairList & pairs,
         int sex_filter, const char * label, bool is_covariate = false);

void   AffectionStatistics(FILE * file, FilteredPedigree & ped, const char * add = NULL);
void   GetAffectionStatistics(FILE * file, FilteredPedigree & ped, int sex_filter,
         const char * label);

void   MarkerSummaries(PDF * pdf);
void   MarkerStatistics(FILE * file, FilteredPedigree & ped, const char * add = NULL);
void   GetMarkerStatistics(FILE * file, FilteredPedigree & ped, int sex_filter,
         const char * label);

// Subroutines for Hardy-Weinberg testing
void   HWStatistics(PDF * pdf = NULL);

// Other optional checks
void   IBDStatistics();
void   ByFamilyStatistics();
void   ByAffectionStatistics();
void   GetByAffectionStatistics(int a, FILE * file, const char * label, int affection_target);

// Summaries of data quality
void QualityStatistics();

// Age checks
void   CheckAges(SimplePairList  & sib_pairs);

// Pair statistics
void   PairStatistics(FILE * file, FilteredPedigree & ped, PedigreePairs & pairs,
         const char * add = NULL, bool traits_only = false);
void   GetPairStatistics(FILE * file, PedigreePairs & pairs,
         const char * label, bool traits_only);
void   PairTraitStatistics(FILE * file, PedigreePairs & pairs, bool is_covariate);
void   GetPairTraitStatistics(FILE * file, PedigreePairs & pairs, bool is_covariate,
         StringArray * corr_txt, IntArray * cts, bool do_unordered = true);
void   PairAffectionStatistics(FILE * file, PedigreePairs & pairs);
void   PrintDetailHeader(FILE * file, const char * additional = NULL);

// Miscellaneous/utility functions
void   LoadAncestors(Person * p1, IntArray & ancestors);
int    Generation(Person * p);
void   MakeIndex(IntArray & index, IntArray & counts);

// global filtering functions
bool   FilterOnCovariate(int var, Person & p1, Person & p2);
bool   FilterOnTrait(int var, Person & p1, Person & p2);
bool   FilterOnGenotype(int marker, Person & p1, Person & p2);

bool   isParent(Person * p1, Person * p2);
bool   isOrderedParent(Person * p1, Person * p2);
bool   isGrandParentOffspring(Person * p1, Person * p2);
bool   isOrderedGrandParentOffspring(Person * p1, Person * p2);
bool   isAvuncular(Person * p1, Person * p2);
bool   isOrderedAvuncular(Person * p1, Person * p2);
bool   isOtherRelative(Person * p1, Person * p2);

bool   isCousin(Person * p1, Person * p2);
bool   isSib(Person * p1, Person * p2);
bool   isHalfSib(Person * p1, Person * p2);

bool   FilterOnBothMale(Person * p1, Person * p2);
bool   FilterOnBothFemale(Person * p1, Person * p2);
bool   FilterOnOppositeSex(Person * p1, Person * p2);

void   ZeroFounderParents(Pedigree & ped);

#endif








 
 
 
 
 
