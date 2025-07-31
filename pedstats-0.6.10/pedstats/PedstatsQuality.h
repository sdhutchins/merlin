////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsQuality.h 
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
 
// PedstatsQuality.h
// Written by Jan Wigginton

#ifndef __PEDSTATSQUALITY_H__
#define __PEDSTATSQUALITY_H__

#include "StringArray.h"
#include "MathMatrix.h"
#include "FilteredPedigree.h"
#include "IntArray.h"

class PedstatsQuality
   {
   public:

      PedstatsQuality(FilteredPedigree & p) ;
      ~PedstatsQuality();

      void Summarize(int window_size, int sex_filter = SEX_UNKNOWN);

   private:

      //Alias for pedigree to summarize
      FilteredPedigree & ped;

      //Size of print window
      int windowSize;

      //Number of markers with some filtered people genotyped
      int nHetMarkers;

      //Average heterozygosity among markers with some filtered people typed
      //Average proportion of filtered individuals genotyped among all genotyped markers
      double avgHet, avgGenoProp;

      //Top/bottom genotype proportions by marker corresponding names, number of genotypes
      Vector      maxGenoProps, minGenoProps;
      StringArray maxGenoPropNames, minGenoPropNames;
      IntArray    nTypedMinHet, nTypedMaxHet;

      //Top/bottom heterozygosity by marker, corresponding names, number of genotypes
      Vector      maxHets, minHets;
      StringArray maxHetNames, minHetNames;
      IntArray    nTypedMinGenoProp, nTypedMaxGenoProp;

      // Total number of genotypes, possible heterozygote genotypes for all filtered individuals, all markers
      int totalGenotypes, totalPossibleHets;

      void GetAllStats(Vector & het_stats, Vector & geno_stats, IntArray & n_typed, IntArray & n_het_typed,
           int sex_filter = SEX_UNKNOWN);

      void GetExtremeStats(Vector & stats, IntArray & counts, bool use_het,
            int window_size);
      void GetExtremeIndices(const Vector & stats, IntArray & max_indices,
            IntArray & min_indices, int window_size);

      void SetValues(Vector & values, StringArray & names, IntArray & counts,
            const       IntArray & indices, const Vector & stats, const IntArray & stat_cts,
            int window_size);

      void Dimension(int window_size);
      void DimensionArrays(Vector & array, StringArray & names, IntArray & counts);

      StringArray & GrabNameArray(bool use_het, bool use_max);
      Vector      & GrabValueVector(bool use_het, bool use_max);
      IntArray    & GrabCountArray(bool use_het, bool use_max);

      void  PrintExtremeValues(bool use_het, int window_size, int sex_filter);
   };

#endif
 
 
 
 
 
