////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsQuality.cpp 
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
 
// Written by Jan Wigginton

#include "PedstatsQuality.h"
#include "QuickIndex.h"
#include "MathConstant.h"

PedstatsQuality::PedstatsQuality(FilteredPedigree & p) : ped(p)
   {
   //Dimension(window_size);
   }

PedstatsQuality::~PedstatsQuality()
   {
   }

void PedstatsQuality::Dimension(int window_size)
   {
   windowSize = window_size;

   DimensionArrays(maxHets, maxHetNames, nTypedMaxHet);
   DimensionArrays(minHets, minHetNames, nTypedMinHet);

   DimensionArrays(maxGenoProps, maxGenoPropNames, nTypedMaxGenoProp);
   DimensionArrays(minGenoProps, minGenoPropNames, nTypedMinGenoProp);

   avgHet = avgGenoProp = _NAN_;
   }

void PedstatsQuality::DimensionArrays(Vector & array, StringArray & names, IntArray & counts)
   {
   array.Clear();
   array.Dimension(windowSize);
   array.Set(_NAN_);

   names.Clear();
   names.Dimension(windowSize);

   counts.Clear();
   counts.Dimension(windowSize);
   counts.Zero();
   }

void PedstatsQuality::Summarize(int window_size, int sex_filter)
   {
   if (Person::markerCount == 0) return;

   if (sex_filter == SEX_UNKNOWN || sex_filter == -99)
      {
      String fill('=', 12);
      printf("DATA QUALITY\n%s\n", (const char *) fill);
      }

   Dimension(window_size);

   String added;
   switch (sex_filter)
     {
     case SEX_MALE : added = "\nAMONG MALES: \n"; break;
     case SEX_FEMALE : added = "\nAMONG FEMALES: \n"; break;
     case SEX_UNKNOWN : added = "\nALL DATA: \n"; break;
     default: added = ""; sex_filter = SEX_UNKNOWN; break;
     }

   printf("%s", (const char *) added);

   Vector het_stats, geno_stats;
   IntArray n_typed, n_het_typed;

   GetAllStats(het_stats, geno_stats, n_typed, n_het_typed, sex_filter);

   int window_size_1 = min(windowSize, geno_stats.SafeCount());
   int window_size_2 = min(windowSize, het_stats.SafeCount());

   GetExtremeStats(geno_stats, n_typed, false, window_size_1);
   GetExtremeStats(het_stats, n_het_typed, true, window_size_2);

   PrintExtremeValues(false, window_size_1, sex_filter);
   printf("\n");
   PrintExtremeValues(true, window_size_2, sex_filter);
   printf("\n");
   }

void PedstatsQuality::GetAllStats(Vector & het_stats, Vector & geno_stats,
 IntArray & n_typed, IntArray & n_het_typed, int sex_filter)
   {
   het_stats.Dimension(Person::markerCount);
   geno_stats.Dimension(Person::markerCount);
   n_typed.Dimension(Person::markerCount);
   n_het_typed.Dimension(Person::markerCount);

   het_stats.Set(_NAN_);
   geno_stats.Set(_NAN_);
   n_typed.Zero();
   n_het_typed.Zero();

   int n_genos, n_het_genos;
   int n_possible_genotypes = 0, n_het_total = 0;

   for (int i = 0; i < ped.count; i++)
      {
      if (ped.filtered[i]) continue;
      if (sex_filter != SEX_UNKNOWN && ped[i].sex != sex_filter) continue;
      n_possible_genotypes ++;
      }

   if (n_possible_genotypes == 0) return;

   int n_possible_hets;
   for (int m = 0;  m < Person::markerCount; m++)
      {
      n_possible_hets = n_genos = n_het_genos = 0;

      for (int i = 0; i < ped.count; i++)
         {
         if (ped.filtered[i]) continue;
         if (sex_filter != SEX_UNKNOWN && ped[i].sex != sex_filter) continue;

         if (ped[i].isGenotyped(m) == false) continue;
         
         n_genos++;
     
         if (ped.chromosomeX && ped[i].sex == SEX_MALE) continue;

         n_possible_hets ++;
         n_het_genos += ped[i].markers[m].isHeterozygous();
         }

      n_typed[m] = n_genos;
      n_het_typed[m] = n_possible_hets;

      n_het_total += n_het_genos;

      geno_stats[m] = (double) n_genos / (double) n_possible_genotypes;

      if (n_possible_hets > 0)
         het_stats[m] = (double) n_het_genos / (double) n_possible_hets;
      }

   totalGenotypes = n_typed.Sum();
   totalPossibleHets = n_het_typed.Sum();

   nHetMarkers = het_stats.SafeCount();

   avgGenoProp = (double) n_typed.Sum() / (double)(Person::markerCount * n_possible_genotypes);
   
   if (totalPossibleHets > 0)
      avgHet = (double) n_het_total / (double) totalPossibleHets;
   }

void PedstatsQuality::GetExtremeStats(Vector & stats, IntArray & counts, bool use_het,
 int window_size)
   {
   if (window_size == 0) return;

   IntArray min_indices(window_size), max_indices(window_size);

   GetExtremeIndices(stats, max_indices, min_indices, window_size);

   Vector & min_stats = GrabValueVector(use_het, false);
   StringArray & min_names = GrabNameArray(use_het, false);
   IntArray & min_counts = GrabCountArray(use_het, false);

   Vector & max_stats = GrabValueVector(use_het, true);
   StringArray & max_names = GrabNameArray(use_het, true);
   IntArray & max_counts = GrabCountArray(use_het, true);

   SetValues(min_stats, min_names, min_counts, min_indices, stats, counts, window_size);
   SetValues(max_stats, max_names, max_counts, max_indices, stats, counts, window_size);
   }

void PedstatsQuality::GetExtremeIndices(const Vector & stats, IntArray & max_indices,
 IntArray & min_indices, int window_size)
   {
   QuickIndex index(stats);

   min_indices.Clear();
   max_indices.Clear();

   int j = 0, ct = 0, n_stats = stats.Length();
   while (j < n_stats)
      {
      if (stats[index[j]] != _NAN_)
         {
         min_indices.Push(index[j]);
         ct ++;
         }
      j++;
      if (ct >= window_size) break;
      }

   j = ct = 0;
   int last = stats.Length() - 1;
   while (j < n_stats)
      {
      if (stats[index[last - j]] != _NAN_)
         {
         max_indices.Push(index[last - j]);
         ct ++;
         }
      j++;

      if (ct >= window_size) break;
      }
   }

StringArray & PedstatsQuality::GrabNameArray(bool use_het, bool use_max)
   {
   if (use_het)
      return use_max ? maxHetNames : minHetNames;

   return use_max ? maxGenoPropNames : minGenoPropNames;
   }

Vector & PedstatsQuality::GrabValueVector(bool use_het, bool use_max)
   {
   if (use_het)
      return use_max ? maxHets: minHets;

   return use_max ? maxGenoProps : minGenoProps;
   }

IntArray & PedstatsQuality::GrabCountArray(bool use_het, bool use_max)
   {
   if (use_het)
      return use_max ? nTypedMaxHet : nTypedMinHet;

   return use_max ? nTypedMaxGenoProp : nTypedMinGenoProp;
   }

void PedstatsQuality::SetValues(Vector & extreme_vals, StringArray & extreme_names,
 IntArray & extreme_cts, const IntArray & indices, const Vector & stats,
 const IntArray & stat_cts, int window_size)
   {
   for (int w = 0; w < window_size; w++)
      {
      int index = indices[w];
      extreme_vals[w] = stats[index];
      extreme_names[w] =  ped.markerNames[index];
      extreme_cts[w] = stat_cts[index];
      }
   }

void PedstatsQuality::PrintExtremeValues(bool use_het, int window_size,
 int sex_filter)
   {
   printf("\n%s %s BY MARKER\n", "HIGHEST AND LOWEST",
            use_het ? "HETEROZYGOSITIES" : "GENOTYPING RATES");

   if (window_size == 0)
      {
      printf("No data\n");
      return;
      }

   StringArray & max_names  = GrabNameArray(use_het, true);
   Vector      & max_values = GrabValueVector(use_het, true);
   IntArray    & max_counts = GrabCountArray(use_het, true);

   StringArray & min_names  = GrabNameArray(use_het, false);
   Vector      & min_values = GrabValueVector(use_het, false);
   IntArray    & min_counts = GrabCountArray(use_het, false);

   String geno_tag = use_het && ped.chromosomeX ? "N_GENO*" : "N_GENO";

   String header;
   header.printf("\n%15s %6s %7s %8s | %13s %6s %7s %8s\n",
            "MARKER", "RANK", use_het ? "HET" : "PROP",  (const char *) geno_tag,
            "MARKER", "RANK", use_het ? "HET" : "PROP", (const char *) geno_tag);
   String fill('-', header.Length() - 2);
   printf("%s\n%s\n", (const char *) header, (const char *) fill);

   int top = use_het ? nHetMarkers : Person::markerCount;

   for (int i = 0; i < window_size; i++)
      {
      printf("%15s %6d %6.1f%% %8d | %13s %6d %6.1f%% %8d\n",
        (const char *) max_names[i], i + 1,   max_values[i] * 100.0, max_counts[i],
        (const char *) min_names[i], top - i, min_values[i] * 100.0, min_counts[i]);
      }

   double average = use_het ? avgHet : avgGenoProp;
   int total_genotypes = use_het ? totalPossibleHets : totalGenotypes;
   printf("\n%15s %5d %6.1f%% %8d\n", "Totals", top, average * 100.0, total_genotypes);

   if (ped.chromosomeX && use_het)
      printf("\n* NOTE:  Heterozygosity calculations and genotype counts exclude males\n");
   if (ped.filterOn)
      printf("\n**  NOTE:  Sample restricted to %s", (const char *) ped.textFilterLabel);
   }








 
 
 
 
 
