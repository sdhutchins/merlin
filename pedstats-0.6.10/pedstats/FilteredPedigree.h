////////////////////////////////////////////////////////////////////// 
// pedstats/FilteredPedigree.h 
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
 
#ifndef __FILTEREDPEDIGREE__
#define __FILTEREDPEDIGREE__

#include "Pedigree.h"

// Written by Jan Wigginton

class FilteredPedigree : public Pedigree
   {
   public:

      FilteredPedigree();
      virtual ~FilteredPedigree() {};

      // Number of individuals in pedigree after filtering
      int remainingCt;

      // True if individuals are to be excluded from pedigree
      bool filterOn;

      // Set elements to true (1) for individuals to be excluded from consideration.
      IntArray filtered;

      // Affection to filter on
      static String affectedFor;

      // Minimum phenotypes, covariates, and genotyped markers for persons retained in pedigree
      static int minPhenos, minCovariates, minGenos;

      // Minimum, maximum family size after filtering families on size
      static int minFamSize, maxFamSize;

      // Labels describing parameters for filter
      String textFilterLabel, pdfFilterLabel;

      // Initializes overall filter
      void InitializeFilter();
      
      // Initializes filter for current pedigree based on settings for an existing pedigree
      void InitializeFilter(FilteredPedigree & ped);

      // Determines which individuals are filtered
      void ApplyFilters(bool silent = false);

      // Builds a new pedigree consisting of just family members of family[id]
      void ExtractFamily(int id, FilteredPedigree & single_fam_ped);

      // Builds a new pedigree consisting of just full entries for pedigree members
      // known to be affected for affections[a], blank elements for all others.
      void ExtractOnAffection(int a, FilteredPedigree & new_ped, int target_status);

      // Merges additional filter with current filtering rules for pedigree
      void ApplyFilter(IntArray & filter);

      virtual void Sort();

   private:

      bool logFiltering;
    
      FILE * filterLog;

      StringArray filterLabels;
      int targetAffection;
      int filtersApplied;

      void RunAffectionFilter(bool silent = false);
      void RunGenotypeFilter(bool silent = false);
      void RunPhenotypeFilter(bool is_covariate, bool silent = false);

      int  BuildAffectionFilter(IntArray & affection_filter);
      int  BuildGenotypeFilter(IntArray & genotype_filter);
      int  BuildPhenotypeFilter(IntArray & phenotype_filter, bool is_covariate);

      void BuildFilterLabels();
      void CheckFilterBounds();

   };

#endif



 
 
 
 
 
 
