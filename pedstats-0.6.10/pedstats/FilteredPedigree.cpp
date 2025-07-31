////////////////////////////////////////////////////////////////////// 
// pedstats/FilteredPedigree.cpp 
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
 
#include "FilteredPedigree.h"
#include "Constant.h"
#include "Error.h"
#include "Sort.h"

#define N_POSSIBLE_FILTERS 4

int    FilteredPedigree::minGenos = 0;
int    FilteredPedigree::minPhenos = 0;
int    FilteredPedigree::minCovariates = 0;
int    FilteredPedigree::minFamSize = 0;
int    FilteredPedigree::maxFamSize = 99999;
String FilteredPedigree::affectedFor = "";

FilteredPedigree::FilteredPedigree()
   {
   targetAffection = -1;
   filterOn = false;
   filtersApplied = 0;
   filterLog = NULL;
   }

void FilteredPedigree::InitializeFilter()
   {
   filtered.Dimension(count);
   filtered.Zero();
   filtersApplied = 0;
   }

void FilteredPedigree::InitializeFilter(FilteredPedigree & ped)
   {
   filterOn = ped.filterOn;
   InitializeFilter();
   }

void FilteredPedigree::ExtractFamily(int id, FilteredPedigree & single_fam_ped)
   {
   Pedigree::ExtractFamily(id, (Pedigree &) single_fam_ped);

   single_fam_ped.InitializeFilter(*this);
   single_fam_ped.ApplyFilters(true);
   }

void FilteredPedigree::ExtractOnAffection(int a, FilteredPedigree & new_ped, int target_status)
   {
   Pedigree::ExtractOnAffection(a, (Pedigree &) new_ped, target_status);
   new_ped.InitializeFilter(*this);
   new_ped.ApplyFilters(true);
   }

// Merges additional filter with current filtering rules for pedigree
void FilteredPedigree::ApplyFilter(IntArray & filter)
   {
   if (filter.Length() != count || filtered.Length() != count)
      error("FilteredPedigree:Size of pedigree filter doesn't match number of persons in pedigree");

   for (int i = 0; i < count; i++)
      filtered[i] |= filter[i];
   }

void FilteredPedigree::Sort()
   {
   Pedigree::Sort();
   //InitializeFilter();
   }

void FilteredPedigree::ApplyFilters(bool silent)
   {
   logFiltering = !silent;

   CheckFilterBounds();
   InitializeFilter();

   filterOn = (FilteredPedigree::minPhenos > 0 || FilteredPedigree::minCovariates > 0
              || FilteredPedigree::minGenos > 0 || FilteredPedigree::targetAffection != -1);

   remainingCt = count;
   filtersApplied = 0;

   if (filterOn == false) return;

   if (logFiltering)
      {
      filterLog = fopen("pedstats.filterlog", "wt");
      if (filterLog == NULL)
         error("FilteredPedigree::Unable to open filter log %s\n", "pedstats.filterlog");
      }

   if (targetAffection != -1)
      RunAffectionFilter(silent);

   if (minGenos > 0)
      RunGenotypeFilter(silent);

   if (minPhenos > 0)
      RunPhenotypeFilter(false, silent);

   if (minCovariates > 0)
      RunPhenotypeFilter(true, silent);

   BuildFilterLabels();

   if (filterLog)
      {
      fprintf(filterLog, "\n\n\n");
      fclose(filterLog);
      }
   }

void FilteredPedigree::CheckFilterBounds()
   {
   if (minGenos > Person::markerCount)
      error("Minimum number of genotypes exceeds number of markers in this dataset");

   if (minCovariates > Person::covariateCount)
      error("Minimum number of measured covariates exceeds number of covariates in this dataset");

   if (minPhenos > Person::traitCount)
      error("Minimum number of trait phenotypes exceeds number of traits in this dataset");

   if (affectedFor != "")
      {
      targetAffection = 0;
      while (targetAffection < Person::affectionCount && Person::affectionNames[targetAffection] != affectedFor)
         targetAffection ++;

      if (targetAffection == Person::affectionCount)
         error("Unable to filter pedigree on affection variable %s. This \n"
               "variable is not present in this data set. Affection variable\n"
               "names are case sensitive and must correspond exactly to that\n"
               "listed in the dat file\n");
      }

   filterLabels.Dimension(N_POSSIBLE_FILTERS);
   }

void FilteredPedigree::RunAffectionFilter(bool silent)
   {
   IntArray affection_filter(count);
   int to_filter = BuildAffectionFilter(affection_filter);

   filterLabels[filtersApplied].printf("affected by %s", (const char *) affectedFor);

   if (!silent)
      {
      printf("\nFiltering sample to remove individuals not diagnosed or unaffected by \n%s...\n",
       (const char *) affectedFor);
      printf("%d of %d individuals removed from pedigree sample\n", to_filter, remainingCt);
      }

   ApplyFilter(affection_filter);
   filtersApplied ++;
   }

void FilteredPedigree::RunGenotypeFilter(bool silent)
   {
   IntArray genotype_filter(count);
   int to_filter = BuildGenotypeFilter(genotype_filter);

   filterLabels[filtersApplied].printf("genotyped on at least %d marker%s",
                                           minGenos, minGenos != 1 ? "s" : "");

   if (!silent)
      {
      printf("\nFiltering sample to remove individuals genotyped on less than %d "
             "marker%s ...\n", minGenos, minGenos != 1 ? "s" : "");
      printf("%d of %d individuals removed from pedigree sample\n", to_filter, remainingCt);
      }

   ApplyFilter(genotype_filter);
   remainingCt -= to_filter;
   filtersApplied ++;
   }

void FilteredPedigree::RunPhenotypeFilter(bool is_covariate, bool silent)
   {
   IntArray phenotype_filter(count);
   int to_filter = BuildPhenotypeFilter(phenotype_filter, is_covariate);

   int cutoff = is_covariate ? minCovariates : minPhenos;

   filterLabels[filtersApplied].printf("%s for at least %d %s%s",
                   is_covariate ? "controlled" : "phenotyped",
                   cutoff, is_covariate ? "covariate" : "trait", cutoff != 1 ? "s" : "");
   if (!silent)
      {
      printf("\nFiltering sample to remove individuals %s on less than %d trait%s ...\n",
             is_covariate ? "controlled" : "phenotyped",
             cutoff, cutoff != 1 ? "s" : "");
      printf("%d of %d individuals removed from pedigree sample\n", to_filter, remainingCt);
      }

   ApplyFilter(phenotype_filter);
   remainingCt -= to_filter;
   filtersApplied ++;
   }

int FilteredPedigree::BuildGenotypeFilter(IntArray & genotype_filter)
   {
   int to_filter = 0;
   genotype_filter.Zero();

   if (logFiltering)
      {
      String fill('=', 44);
      fprintf(filterLog, "\n\nGENOTYPE FILTER: ");
      fprintf(filterLog, "\n\n%15s %15s %12s\n\n%s\n", "FAMILY", "PERSON", "N_GENOTYPES", (const char *) fill);
      }

   bool filter_person;
   Person * p;

   for (int i = 0; i < count; i++)
      {
      p = persons[i];

      filter_person = (p->ngeno < minGenos);
      genotype_filter[i] = filter_person;

      to_filter += (filtered[i] == 0 && filter_person == 1);

      if (logFiltering && filter_person)
         fprintf(filterLog, "%15s %15s %12d\n", (const char *) p->famid, (const char *) p->pid, p->ngeno);
      }

   return to_filter;
   }

int FilteredPedigree::BuildPhenotypeFilter(IntArray & phenotype_filter, bool is_covariate)
   {
   phenotype_filter.Zero();

   int n_phenos, to_filter = 0;
   int min_values = is_covariate ? minCovariates : minPhenos;
   int upper_index = is_covariate ? Person::covariateCount : Person::traitCount;

   bool filter_person;
   Person * p;

   if (logFiltering)
      {
      String fill('=', 44);
      fprintf(filterLog, "\n\n%s FILTER: ", is_covariate ? "COVARIATE" : "TRAIT");
      fprintf(filterLog, "\n\n%15s %15s %12s\n\n%s\n", "FAMILY", "PERSON",
       is_covariate ? "N_COVARIATES" : "N_TRAITS", (const char *) fill);
      }

   for (int i = 0; i < count; i++)
      {
      n_phenos = 0;
      p = persons[i];

      for (int j = 0; j < upper_index; j++)
         n_phenos += is_covariate ? p->isControlled(j) : p->isPhenotyped(j);

      phenotype_filter[i] = (n_phenos < min_values);
      filter_person       = (n_phenos < min_values);

      to_filter += (filtered[i] == 0 && filter_person);

      if (logFiltering && filter_person)
         fprintf(filterLog, "%15s %15s %12d\n", (const char *) p->famid, (const char *) p->pid, n_phenos);
      }

   return to_filter;
   }

int FilteredPedigree::BuildAffectionFilter(IntArray & affection_filter)
   {
   int to_filter = 0;
   affection_filter.Zero();

   bool filter_person;
   Person * p;

   if (logFiltering)
      {
      String fill('=', 44);
      fprintf(filterLog, "\n\nAFFECTION FILTER: ");
      fprintf(filterLog, "\n\n%15s %15s %12s\n\n%s\n", "FAMILY", "PERSON", "AFFECTION",
       (const char *) fill);
      }

   for (int i = 0; i < count; i++)
      {
      p = persons[i];

      filter_person       = (p->affections[targetAffection] != 2);
      affection_filter[i] = (p->affections[targetAffection] != 2);

      to_filter += (filtered[i] == 0 && filter_person);
      if (logFiltering && filter_person)
         fprintf(filterLog, "%15s %15s %10d\n", (const char *) p->famid, (const char *) p->pid,
          p->affections[targetAffection]);
      }
   return to_filter;
   }

void FilteredPedigree::BuildFilterLabels()
   {
   StringArray switches(3);
   for (int i = 0; i < 3; i++)
      switches[i] = "";

   if (filtersApplied > 1)
      {
      switches[filtersApplied - 2] = " and ";
      if (filtersApplied > 2)
         {
         switches[filtersApplied - 3] = ", ";
         if (filtersApplied > 3)
            switches[filtersApplied - 4] = ", ";
         }
      }

   textFilterLabel.printf("individuals %s%s\n%s%s%s%s\n%s%s%s",
                         (const char *) filterLabels[0],
                         (const char *) switches[0], "  ",
                         (const char *) filterLabels[1], (const char *) switches[1],
                         (const char *) filterLabels[2],
                         switches[2] == " and " ? " " : "  ", (const char *) switches[2],
                         (const char *) filterLabels[3]);

   pdfFilterLabel.printf("* NOTE:  Sample restricted to individuals %s%s%s%s%s%s%s",
                         (const char *) filterLabels[0], (const char *) switches[0],
                         (const char *) filterLabels[1], (const char *) switches[1],
                         (const char *) filterLabels[2], (const char *) switches[2],
                         (const char *) filterLabels[3]);
   }




 
 
 
 
