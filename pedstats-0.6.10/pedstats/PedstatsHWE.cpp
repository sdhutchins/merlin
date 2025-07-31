////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsHWE.cpp 
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
 
#include "PedstatsHWE.h"
#include "MathStats.h"
#include "MathConstant.h"
#include "Kinship.h"
#include "PDFhistogram.h"

#include <math.h>

bool    PedstatsHWE::showAll = false;
double  PedstatsHWE::significanceCutoff = 0.05;
FILE  * PedstatsHWE::outputFile = stdout;

PedstatsHWE::PedstatsHWE(FilteredPedigree & p) : ped(p)
   {
   InitializeStatistics();
   setType = stAll;
   filesOpened = false;
   }

PedstatsHWE::~PedstatsHWE()
   {
   if (filesOpened)
      fclose(outputFile);
   }

void PedstatsHWE::OpenFiles(const char * filename)
   {
   outputFile = fopen(filename, "wt");
   if (outputFile == NULL)
      error("PedstatsHWE::Unable to open file %s for output", filename);
   filesOpened = true;
   }

void PedstatsHWE::InitializeStatistics()
   {
   markerTests = attemptedTests = 0;
   failedTests = pooledTests = failedPooledTests = 0;
   failedTests_01 = failedTests_05 = 0;

   sparseNames.Clear();
   monomorphicNames.Clear();
   genotypedProportions.Clear();
   genotypedProportions.Dimension(ped.count);
   genotypedProportions.Zero();
   }

void PedstatsHWE::InitializeArrays()
   {
   observed.Dimension(alleleCount, alleleCount);
   observed.Zero();

   alleleFreqs.Dimension(alleleCount);
   alleleFreqs.Zero();

   pooledAlleles.Dimension(alleleCount);
   pooledAlleles.Zero();
   }

bool PedstatsHWE::TestMarker(int m)
   {
   currMarker = m;
   poolAllele = -1;

   hetProbs.Zero();
   alleleCount = ped.CountAlleles(m);
   testPerformed = false;
   chisq = _NAN_;

   attemptedTests ++;

   if (alleleCount <= 1)
      {
      monomorphicNames.Push(Person::markerNames[m]);
      return false;
      }

   InitializeArrays();
   DeterminePooledAlleleInfo(m);

   if (totalGenotypes <= 3)
      {
      sparseNames.Push(Person::markerNames[m]);
      return false;
      }

   if (realAlleles <= 1)
      {
      monomorphicNames.Push(Person::markerNames[m]);
      return false;
      }

   chisq = DetermineGenotypeCounts();

   useExact = realAlleles == 2 && totalGenotypes < 5000;
   pValue = useExact ? RunExactTest() :
                    realAlleles == 1.0 ? 1.0 : chidist(chisq, realAlleles * (realAlleles - 1) * 0.5);

   markerTests ++;

   bool fail = pValue < significanceCutoff;
   pooledTests += (pooledAlleles.Sum() > 1);
   failedTests += (fail);
   failedTests_01 += pValue < 0.01;
   failedTests_05 += pValue < 0.05;
   failedPooledTests += (fail && pooledAlleles.Sum() > 1);
   testPerformed = true;

   return fail;
   }

int PedstatsHWE::DetermineExpectedHeterozygotes()
   {
   double expected_hets = 0;
   int n_alleles = alleleFreqs.Length();

   for (int i = 0; i < n_alleles; i++)
      for (int j = i + 1; j < n_alleles; j++)
         expected_hets += (totalGenotypes * 2.0 * alleleFreqs[i] * alleleFreqs[j]);

   return (int) expected_hets;
   }

void PedstatsHWE::DeterminePooledAlleleInfo(int m)
   {
   smallestAllele = largestAllele = -1;
   alleleCount = alleleFreqs.Length();

   DetermineAlleleFreqs(false);
   if (totalGenotypes < 3) return;
   DetermineAlleleRange();
   PoolAlleles(alleleFreqs);
   DetermineAlleleFreqs(true);

   expectedHets = DetermineExpectedHeterozygotes();

   if (totalGenotypes > 3 && realAlleles > 1)
      DetermineAlleleNames();
   }

void PedstatsHWE::PoolAlleles(Vector & freqs)
   {
   // Find the minimum useful allele frequency
   double min_freq = sqrt(3.0 / (double) totalGenotypes);

   // List all rare alleles and find their total frequency
   double pooled_sum = 0.0;

   // Also find the rarest allele with frequency >= min_freq
   int mini = -1, commonCount = 0, maxi = -1;

   for (int i = 0; i < alleleCount; i++)
      {
      if (freqs[i] == 0.0)
         continue;

      if (freqs[i] < min_freq)
         {
         pooled_sum += freqs[i];
         pooledAlleles[i] = true;
         if (maxi == -1 || freqs[i] > freqs[maxi])
            maxi = i;
         }
      else
         {
         if (mini == -1 || freqs[i] < freqs[mini])
            mini = i;
         commonCount ++;
         }
      }

   // If all rare alleles combined are still too rare...
   if (commonCount == 0)
      {
      pooled_sum -= freqs[maxi];
      pooledAlleles[maxi] = false;
      }
   else if (pooled_sum > 0.0 && pooled_sum < min_freq && commonCount > 1)
      {
      pooled_sum += freqs[mini];
      pooledAlleles[mini] = true;
      }

   // Glob all rare alleles into 1
   poolAllele = alleleCount - 1;
   if (pooled_sum > 0.0)
      while (pooledAlleles[poolAllele] == false)
         poolAllele--;
   }

void PedstatsHWE::DetermineAlleleFreqs(bool final_count)
   {
   totalGenotypes = 0;
   alleleFreqs.Zero();
   observed.Zero();

   for (int i = 0; i < ped.count; i++)
      {
      if (personsSelected[i] == 0 || !(ped[i].isGenotyped(currMarker)))
         continue;

      int one = ped[i].markers[currMarker].one - 1;
      int two = ped[i].markers[currMarker].two - 1;

      if (final_count)
         {
         if (pooledAlleles[one]) one = poolAllele;
         if (pooledAlleles[two]) two = poolAllele;

        (one > two) ? observed[one][two] ++ : observed[two][one] ++;
        }

      alleleFreqs[one] ++;
      alleleFreqs[two] ++;
      totalGenotypes ++;
      }

   if (totalGenotypes <= 3)
      return;

   // Work out allele frequencies
   alleleFreqs.Multiply(1.0 / (2.0 * totalGenotypes));

   if (final_count)
      {
      realAlleles = 0;
      for (int i = 0; i < alleleCount; i++)
         realAlleles += (alleleFreqs[i] > 0);
      }
   }

void PedstatsHWE::DetermineAlleleRange()
   {
   for (int i = 0; i < alleleCount; i++)
      if (alleleFreqs[i] > 0)
         {
         if (smallestAllele == -1) smallestAllele = i+1;
         largestAllele = i + 1;
         }
   }

void PedstatsHWE::DetermineAlleleNames()
   {
   alleleNames.Dimension(realAlleles);

   pooledNames = "";
   int pool_ct = pooledAlleles.Sum();

   MarkerInfo * marker_info = ped.markerInfo[currMarker];

   if (pool_ct > 1)
      for (int i = 0; i < alleleCount; i++)
         if (pooledAlleles[i])
            {
            pooledNames += marker_info->alleleLabels[i + 1];
            pooledNames += " ";
            }

   int j = 0;
   for (int i = 0; i < alleleCount; i++)
      {
      if (alleleFreqs[i] != 0.0)
         {
         if (i == poolAllele && pool_ct > 1)
            alleleNames[j] = "* P";
         else
            alleleNames[j].printf("%s", (const char *) marker_info->alleleLabels[i+1]);
         j++;
         }
      }
   }

double PedstatsHWE::DetermineGenotypeCounts()
   {
   realExpected.Dimension(realAlleles, realAlleles);
   realObserved.Dimension(realAlleles, realAlleles);

   realObserved.Zero();
   realExpected.Zero();

   double expected, chi = 0.0;
   int allele_count = alleleFreqs.Length();

   int l = 0, k = 0;
   for (int i = 0; i < allele_count; i++)
      if (alleleFreqs[i] > 0)
         {
         l = 0;
         for (int j = 0; j <= i; j++)
            if (alleleFreqs[j] > 0)
               {
               expected = alleleFreqs[i] * alleleFreqs[j] * totalGenotypes * (i == j ? 1 : 2);
               chi += square(expected - observed[i][j]) / expected;

               realExpected[k][l] = expected;
               realObserved[k][l] = observed[i][j];
               l ++;
               }
         k ++;
         }
   return chi;
   }

double PedstatsHWE::RunExactTest()
   {
   IntArray allele(2);
   allele.Zero();

   int k = 0;
   for (int j = 0; j < alleleCount; j++)
      if (alleleFreqs[j] != 0)
         {
         if (k > 1)
             error("HWE::Expected two alleles, current allele is %d%d",
                   k + 1, k == 2 ? "rd" : "th");

         allele[k] = j;
         k++;
         }

   double n_rare_alleles = 2.0 * totalGenotypes;
   n_rare_alleles *= alleleFreqs[allele[0]] < alleleFreqs[allele[1]] ?
                                                  alleleFreqs[allele[0]] :
                                                  alleleFreqs[allele[1]] ;
   n_rare_alleles += 1e-6;

   if (hetProbs.Length() < n_rare_alleles + 1)
      hetProbs.Dimension((int) n_rare_alleles * 2);

   return ExactHWE(int(n_rare_alleles), int(realObserved[1][0]));
   }

double PedstatsHWE::ExactHWE(int rare , int heterozygotes)
   {
   if (rare > totalGenotypes)
      rare = 2 * totalGenotypes - rare;

   if (heterozygotes > rare)
      error("PedstatsHWE: %d heterozygotes but only %d minor alleles");

   if ((rare & 1) ^(heterozygotes & 1))
      error("PedstatsHWE: %s number of minor alleles, but %s number of heterozygotes",
           (rare & 1) ? "odd" : "even",
           (heterozygotes & 1) ? "odd" : "even");

   hetProbs.Dimension(rare + 1);
   hetProbs.Zero();

   // start at midpoint
   int mid = rare * (2 * totalGenotypes - rare) / (2 * totalGenotypes);

   // check to ensure that midpoint and rare alleles have same parity
   if ((rare & 1) ^ (mid & 1))
      mid++;

   int het = mid;
   int hom_r = (rare - mid) / 2;
   int hom_c = totalGenotypes - het - hom_r;

   hetProbs[mid] = 1.0;
   double sum = hetProbs[mid];

   for (het = mid; het > 1; het -= 2)
      {
      hetProbs[het - 2] = hetProbs[het] * het * (het - 1.0)
                            /(4.0 * (hom_r + 1.0) * (hom_c + 1.0));
      sum += hetProbs[het - 2];

      // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
      hom_r++;
      hom_c++;
      }

   het = mid;
   hom_r = (rare - mid) / 2;
   hom_c = totalGenotypes - het - hom_r;

   for (het = mid; het <= rare - 2; het += 2)
      {
      hetProbs[het + 2] = hetProbs[het] * 4.0 * hom_r * hom_c
                            /((het + 2.0) * (het + 1.0));
      sum += hetProbs[het + 2];

      // add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
      hom_r--;
      hom_c--;
      }

   hetProbs /= sum;

   double p_rank = 0.0;
   for (int het = 0; het <= rare; het++)
      {
      if (hetProbs[het] > hetProbs[heterozygotes])
         continue;
      p_rank += hetProbs[het];
      }

   return min(1.0, p_rank);
   }

void PedstatsHWE::TestAllMarkers(PedstatsHWESetType t, PDF * pdf)
   {
   setType = t;
   InitializeStatistics();
   SelectTestSample();
   PrintTableHeader();

   for (int m = 0; m < Person::markerCount; m++)
      {
      bool fail = TestMarker(m);
      if (!testPerformed || !fail && !showAll) continue;

      if (pdf != NULL) GraphTestResult(*pdf);
      PrintTestResult();
      }

   PrintTableFooter();
   if (setType == stSelected)
      PrintSampleSelection();
   }

void PedstatsHWE::PrintSampleSelection()
   {
   if (personsSelected.Sum() <= 0) return;

   String select_filename = "pedstats.hweselection";
   FILE * select_file = fopen((const char *) select_filename, "wt");
   if (select_file == NULL)
      error("PedstatsHWE : Unable to open file %s", (const char *) select_filename);

   String fill('=', 47 + ped.filterOn * 2);
   fprintf(select_file, "\nINDIVIDUALS SELECTED FOR HARDY-WEINBERG TESTING%s\n%s\n\n", ped.filterOn ? " *" : "", 
       (const char *) fill);

   fprintf(select_file, "%15s %15s %15s %15s\n", "FAMILY", "PERSON", "STATUS", "GENOTYPE PROP");

   for (int i = 0; i < ped.count; i++)
      if (personsSelected[i])
         fprintf(select_file, "%15s %15s %15s %15.5f\n", (const char *) ped[i].famid,
            (const char *) ped[i].pid, ped[i].isFounder() ? "FOUNDER" : "NON-FOUNDER",
            genotypedProportions[i]);


   if (ped.filterOn)
      fprintf(select_file, "\n** NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);
   else 
      fprintf(select_file, "\n\n");

   if (select_file)
      fclose(select_file);
   }

void PedstatsHWE::SelectTestSample()
   {
   personsSelected.Clear();
   personsSelected.Dimension(ped.count);
   for (int i = 0; i < ped.count; i++)
      personsSelected[i] = (ped.filtered[i] == 0 && (!ped.chromosomeX || ped[i].sex != SEX_MALE));

   if (setType == stFounders)
      for (int i = 0; i < ped.count; i++)
         personsSelected[i] &= ped[i].isFounder();
   else if (setType == stSelected)
      SelectIndependentSet();
   
   ConstructSampleLabels();
   }

void PedstatsHWE::ConstructSampleLabels()
   {
   textSetTypeLabel = ped.chromosomeX ? "Among All Females" : "Among All Individuals";
   pdfSetTypeLabel =  textSetTypeLabel;

   if (setType == stFounders)
      {
      textSetTypeLabel = ped.chromosomeX ? "Among Founder Females" : "Among Founders";
      pdfSetTypeLabel = ped.chromosomeX ? "Among Founder Females" : "Among Founders";
      }
   else if (setType == stSelected)
      {
      textSetTypeLabel.printf("Using %d Unrelated %s", personsSelected.Sum(), 
                               ped.chromosomeX ? "Females" : "Individuals");
      pdfSetTypeLabel.printf("Among %d Unrelated %s", personsSelected.Sum() , 
                              ped.chromosomeX ? "Females" : "Individuals");
      }

   if (ped.filterOn)
      {
      textSetTypeLabel += " *";
      pdfSetTypeLabel += " *";
      }
   }
void PedstatsHWE::SelectIndependentSet()
   {
   //Determine which markers will be tested
   IntArray tested_markers(ped.markerCount);
   for (int m = 0; m < ped.markerCount; m++)
      tested_markers[m] = IsTested(m);

   if (tested_markers.Sum() == 0) return;

   // Determine per-family genotyping thresholds for inclusion in sampling set
   DetermineGenotypedProportions(tested_markers);
   DetermineTypedCutoffs();

   // Limit set to individuals with minimum proportion of tested markers typed
   LimitToGenotypedPersons(tested_markers);

   for (int fam = 0; fam < ped.familyCount; fam++)
      if (CheckFounderTyping(fam))
         SelectFounders(fam);
      else
         SelectIndependentMembers(fam, tested_markers);
   }

void PedstatsHWE::DetermineGenotypedProportions(IntArray & tested_markers)
   {
   int n_possible = tested_markers.Sum();

   genotypedProportions.Dimension(ped.count);
   genotypedProportions.Zero();

   int n_genotypes;
   for (int i = 0; i < ped.count; i++)
      {
      n_genotypes = 0;
      for (int m = 0; m < ped.markerCount; m++)
         n_genotypes += (tested_markers[m] && ped[i].isGenotyped(m));

      if (n_possible > 0)
         genotypedProportions[i] = (double) n_genotypes / (double) n_possible;
      }
   }

void PedstatsHWE::DetermineTypedCutoffs()
   {
   typedCutoffs.Dimension(ped.familyCount);
   typedCutoffs.Zero();

   double max_geno_prop;
   for (int fam = 0; fam < ped.familyCount; fam++)
      {
      max_geno_prop = 0.0;
      for (int i = ped.families[fam]->first; i <= ped.families[fam]->last; i++)
         {
         if (personsSelected[i] == 0) continue;

         max_geno_prop = max(max_geno_prop, genotypedProportions[i]);
         }

      typedCutoffs[fam] = max_geno_prop > 0.0 ? 0.9 * max_geno_prop : 1.0;
      }
   }

void PedstatsHWE::LimitToGenotypedPersons(const IntArray & tested_markers)
   {
   if (tested_markers.Sum() == 0) return;

   for (int fam = 0; fam < ped.familyCount; fam++)
      for (int i = ped.families[fam]->first; i <= ped.families[fam]->last; i++)
         {
         if (personsSelected[i] == 0) continue;

         bool eliminate = (genotypedProportions[i] <= typedCutoffs[fam]);
         personsSelected[i] &= (eliminate == false);
         }
   }

bool PedstatsHWE::CheckFounderTyping(int fam)
   {
   Family * family = ped.families[fam];

   //check if each family founder is typed and has not been filtered
   for (int i = family->first; i <= family->last; i++)
      if (ped[i].isFounder() && personsSelected[i] == 0)
         return false;

   return true;
   }

void PedstatsHWE::SelectFounders(int fam)
   {
   Family * family = ped.families[fam];

   for (int i = family->first; i <= family->last; i++)
      personsSelected[i] &= ped[i].isFounder();
   }

void PedstatsHWE::SelectIndependentMembers(int fam, const IntArray & markers_tested)
   {
   Family * family = ped.families[fam];
   Kinship kin;
   kin.Setup(*family);

   // unselected individuals have already been eliminated
   IntArray eliminated(family->count);
   for (int i = 0; i < family->count; i++)
      eliminated[i] = (personsSelected[i + family->first] == 0);

   while (eliminated.Sum() < family->count)
      {
      int pick = family->first, n_relatives, min_relatives = 99999;

      for (int i = family->first; i <= family->last; i++)
         {
         if (eliminated[i - family->first]) continue;

         n_relatives = 0;
         for (int j = family->first; j <= family->last; j++)
            {
            if (i == j || eliminated[j - family->first]) continue;
            n_relatives += (kin(ped[i],ped[j]) > 0);
            }

         if (n_relatives < min_relatives)
            {
            min_relatives = n_relatives;
            pick = i;
            }
         }

      personsSelected[pick] &= 1;

      for (int i = family->first; i <= family->last; i++)
         {
         if (eliminated[i - family->first]) continue;
         bool eliminate = kin(ped[pick], ped[i]) > 0;
         eliminated[i - family->first] |= eliminate;
         if (i != pick && eliminate)
            personsSelected[i] = 0;
         }
      }
   }

bool PedstatsHWE::IsTested(int m)
   {
   // monomorphic markers
   int allele_count = ped.CountAlleles(m);
   if (allele_count <= 1)
      return false;

   // mostly untyped markers
   int i = -1, total_genotypes = 0;
   int upper = ped.count - 1;
   while (total_genotypes < 3 && i < upper)
      {
      i ++;
      if (personsSelected[i] == 0) continue;
      total_genotypes += ped[i].isGenotyped(m);
      }
   if (total_genotypes < 3)
      return false;

   // markers that will be monomorphic after pooling
   total_genotypes = 0;
   Vector allele_freqs(allele_count);
   allele_freqs.Zero();
   for (int i = 0; i < ped.count; i++)
      {
      if (personsSelected[i] == 0 || !ped[i].isGenotyped(m)) continue;

      int one = ped[i].markers[m].one - 1;
      int two = ped[i].markers[m].two - 1;

      allele_freqs[one] ++;
      allele_freqs[two] ++;
      total_genotypes ++;
      }

   allele_freqs.Multiply(1.0 / (2.0 * total_genotypes));
   double min_freq = sqrt(3.0 / (double) total_genotypes);

   int real_alleles = 0;
   for (int i = 0; i < allele_count; i++)
      real_alleles += (allele_freqs[i] >= min_freq);

   if (real_alleles == 1)
      {
      double sum = 0.0;
      for (int i = 0; i < allele_count; i++)
         if (allele_freqs[i] >= min_freq)
            sum += allele_freqs[i];
      if (sum < min_freq)
         return false;
      }

   return true;
   }

void PedstatsHWE::PrintTableHeader()
   {
   String fill('=', textSetTypeLabel.Length() + 21);

   fprintf(outputFile, "\n\n\nHARDY-WEINBERG CHECK %s\n%s\n",
           (const char *) textSetTypeLabel.AsUpper(), (const char *) fill);
   }

void PedstatsHWE::PrintTestResult()
   {
   //int allele_count = pooledAlleles.Count();
   int real_alleles = alleleNames.Length();
   int rare_count = pooledAlleles.Sum();

   String allele_ranges, allele_counts, hom_counts, p_value_printed;

   if (rare_count > 1)
      allele_counts.printf("%2d, %2d pooled", real_alleles, rare_count);
   else
      allele_counts.printf("%2d", real_alleles);

   int allelect_col_width = 13;
   int range_col_width = 7;
   int hom_col_width = ped.count > 9999 ? 17 : 15;
   int spacer = hom_col_width == 17 ? 1 : 2;

   if (showAll && markerTests == 1 || !showAll && failedTests == 1)
      fprintf(outputFile, "%15s %*s %6s %6s%*s%*s %*s %7s\n",
                "", hom_col_width, "N_HOM", "N_HET", "E_HET", spacer, "",
                allelect_col_width, "N_ALLELES", range_col_width, "ALLELES", "P-VALUE");

   if (pValue > 0.0001)
      p_value_printed.printf("%7.4f", pValue);
   else
      p_value_printed.printf("%7.*e", pValue < 9.9e-99 ? 0 : 1, pValue);

   int total_hom = 0;
   for (int i = 0; i < real_alleles; i++)
      total_hom += int(realObserved[i][i]);

   if (useExact)
      {
      String name_0 = (alleleNames[0] == "* P" ? "*" : (const char *) alleleNames[0]);
      String name_1 = (alleleNames[1] == "* P" ? "*" : (const char *) alleleNames[1]);

      bool reverse = (realObserved[0][0] > realObserved[1][1]);

      int rare_hom = reverse ? int(realObserved[1][1]) : int(realObserved[0][0]);

      hom_counts.printf("%d, %d rare", total_hom, rare_hom);

      allele_ranges.printf("%2s/%-2s", (const char *) (reverse ? name_1 : name_0),
         (const char *) (reverse ? name_0 : name_1));

      fprintf(outputFile, "%15.15s %*s %6d %6d%*s%*s %*s %7s %s\n",
          (const char *) Person::markerNames[currMarker],
          hom_col_width, (const char *) hom_counts.Trim(), int(realObserved[1][0]), expectedHets,
          spacer, "", allelect_col_width, (const char *) allele_counts,
          range_col_width, (const char *) allele_ranges.Trim(),
          (const char *) p_value_printed.Trim(), "E");
      }
   else
      {
      allele_ranges.printf("%2d-%-2d", smallestAllele, largestAllele);

      fprintf(outputFile, "%15s %*d %6d %6d%*s%*s %*s %7s %s\n",
             (const char *) Person::markerNames[currMarker],
             hom_col_width, total_hom, totalGenotypes - total_hom, expectedHets, spacer, "",
             allelect_col_width, (const char *) allele_counts, range_col_width,
             (const char *) allele_ranges.Trim(), (const char *) p_value_printed.Trim(), "A");
      }
   }

void PedstatsHWE::PrintTableFooter()
   {
   if (markerTests == 0)
      {
      fprintf(outputFile, "No markers tested.\n\n");
      return;
      }
   
   if (markerTests == 1)
      fprintf(outputFile, failedTests == 0 ? "All tested markers okay.\n\n" : "\n");
   else
      fprintf(outputFile, failedTests == 0 ? "All %d tested markers okay.\n\n" : "\n", markerTests);


   if ((!showAll && failedPooledTests > 0) || (showAll && pooledTests > 0))
      fprintf(outputFile, "*  Denotes pooled alleles\n");

   if (ped.filterOn)
      fprintf(outputFile, "\n** NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);

   fprintf(outputFile, "\n%15s %15s %15s %15s %15s\n%15s %15d %15d %15d %15d\n", "", "Attempted",
           "Performed", "Failed [0.05]", "Failed [0.01]", "Total Tests", attemptedTests, markerTests,
           failedTests_05, failedTests_01);

   int sparse_ct = sparseNames.Length();
   if (sparse_ct > 0)
      {
      fprintf(outputFile, "The following markers were not tested for Hardy-Weinberg "
              "due to a low \nnumber of genotyped individuals:\n");

      for (int i = 0, line = 80; i < sparse_ct; i++)
         {
         if (line + sparseNames[i].Length() + 1 > 79)
            fprintf(outputFile, "\n "), line = 3;

         fprintf(outputFile, "%s ", (const char *) sparseNames[i]);
         line += sparseNames[i].Length() + 1;
         }
      }
   fprintf(outputFile, "\n\n");

   int exact_sparse = monomorphicNames.Length();
   if (exact_sparse > 0)
      {
      fprintf(outputFile, "The following markers could not be tested for Hardy-Weinberg "
              "\nbecause they were monomorphic.\n");

      for (int i = 0, line = 80; i < exact_sparse; i++)
         {
         if (line + monomorphicNames[i].Length() + 1 > 79)
            fprintf(outputFile, "\n "), line = 3;

         fprintf(outputFile, "%s ", (const char *) monomorphicNames[i]);
         line += monomorphicNames[i].Length() + 1;
         }
      }
   fprintf(outputFile, "\n\n");
   }

void PedstatsHWE::GraphTestResult(PDF & pdf)
   {
   pdf.page.OpenPage();

   PDFGrid hw_graph;
   hw_graph.printUpper = false;
   hw_graph.title.printf("Observed");

   hw_graph.Dimension(realObserved.rows, realObserved.cols);
   hw_graph.xAxis.label = "Allele 1";
   hw_graph.yAxis.label = "Allele 2";

   for (int i = alleleNames.Length() - 1; i >= 0; i--)
      {
      hw_graph.xAxis.SetTickLabel(i, alleleNames[i]);
      hw_graph.yAxis.SetTickLabel(i, alleleNames[i]);
      }

   double ratio, ratio_g, ratio_b, ratio_r;
   double expected_max = realExpected.SafeMax();

   for (int i = 0; i < realObserved.rows; i++)
      for (int j = 0 ; j <= i; j++)
         {
         hw_graph.SetCellText(i, j, realObserved[i][j], 1);
         ratio = min(1.1, realObserved[i][j] / expected_max);
         ratio *= 2.5;

         ratio_b = min(1.0, 3.3 - ratio);
         ratio_b = max(0.0, ratio_b);

         ratio_g = max(0.0, 2.0 - ratio);
         ratio_g = min(ratio_g , 1.0);

         ratio_r = max(0.0, 1.0 - ratio) ;
         ratio_r = min(1.0, ratio_r);

         hw_graph.SetCellColor(i, j, ratio_r , ratio_g , ratio_b);
         }

   hw_graph.DrawInGrid(pdf, 0, 0, 1, 3, 0.35);

   useExact ? GraphExact(pdf) : GraphAsymptotic(pdf, hw_graph, expected_max);

   String title, subtitle;
   title.printf("%sHardy-Weinberg Test %s for %s", useExact ? "Exact " : "",
                (const char *) pdfSetTypeLabel, (const char *) Person::markerNames[currMarker]);

   String p_print;
   p_print.printf(pValue > 0.0001 ? "%.4f" : "%.1e", pValue);

   if (!useExact)
      subtitle.printf("( Chi-squared: %5.4lf, p = %s%s )", chisq, (const char *) p_print,
                    ped.filterOn ? "**" : "" );
   else
      subtitle.printf("( Observed heterozygotes: %d, p = %s%s )", int(realObserved[1][0]),
        (const char *) p_print,  ped.filterOn ? "**" : "" );

   MakeSignificanceKey(pdf, pValue);

   pdf.page.AddHeader((const char *) title, (const char *) subtitle);
   if (ped.filterOn)
      pdf.page.AddFootNote(0.05, (const char *) ("*" + ped.pdfFilterLabel), 0.08);

   if (poolAllele != -1)
      {
      String rare_note;
      rare_note.printf("*  NOTE:  P represents pool of alleles %s", (const char *) pooledNames);
      pdf.page.AddFootNote(0.05, (const char *) rare_note, 0.10);
      }

   pdf.page.ClosePage();
   }

void PedstatsHWE::GraphAsymptotic(PDF & pdf, PDFGrid & hw_graph, double expected_max)
   {
   double ratio, ratio_g, ratio_b, ratio_r;
   hw_graph.title.printf("Expected");

   for (int i = 0; i < realExpected.rows; i++)
      for (int j = 0 ; j <= i; j++)
         {
         hw_graph.SetCellText(i, j, realExpected[i][j], 1);
         ratio = min(1.1, realExpected[i][j]/expected_max);
         ratio *= 2.5;

         ratio_b = min(1.0, 3.3 - ratio);
         ratio_b = max(0.0, ratio_b);

         ratio_g  = max(0.0, 2.0 - ratio);
         ratio_g = min(1.0, ratio_g);

         ratio_r = max(0.0, 1.0 - ratio);
         ratio_r = min(1.0, ratio_r);

         hw_graph.SetCellColor(i, j, ratio_r ,  ratio_g , ratio_b);
         }

   hw_graph.DrawInGrid(pdf, 0, 1, 1, 3, 0.35);

   hw_graph.title.printf("Residuals");

   double resid, diff;
   for (int i = 0; i < realExpected.rows; i++)
      for (int j = 0 ; j <= i; j++)
         {
         diff = realObserved[i][j] - realExpected[i][j];
         resid = diff / sqrt(realExpected[i][j]);
         ratio = fabs(resid);

         ratio_r = 1.0;
         ratio_g = ratio < 1.0 ? 1.0 : 2.5 - ratio;
         ratio_b = max(0.0, 1.0 - ratio);

         hw_graph.SetCellText(i, j, diff, 1);
         hw_graph.SetCellColor(i, j, ratio_r, ratio_g, ratio_b);
         }

   hw_graph.DrawInGrid(pdf, 0, 2, 1, 3, 0.35);
   }

void PedstatsHWE::GraphExact(PDF & pdf)
   {
   int observed_hets = int(realObserved[1][0]);
   int min_hom = realObserved[1][1] > realObserved[0][0] ? 0 : 1;
   int rare_copies = int(observed_hets + 2 * realObserved[min_hom][min_hom] + 1e-6);
   int n_poss_hets = hetProbs.Length();

   PDFHistogram prob_graph;
   prob_graph.title.printf("Heterozygote probability distribution");
   prob_graph.subTitle.printf("( %d rare allele copies )", rare_copies);

   int max = n_poss_hets - 1;
   while (hetProbs[max] <= 1e-6)
      max--;

   int min = 0;
   while (hetProbs[min] <= 1e-6)
      min++;

   if (max <= min)
      {
      max = n_poss_hets - 1;
      min = 0;
      }

   Matrix graph_probs(3, max - min + 1);
   graph_probs.Zero();

   int series;
   for (int i = min; i <= max; i++)
      {
      graph_probs[0][i - min] = i * 1.0;
      series = hetProbs[i] <= hetProbs[observed_hets] ? 1 : 2;
      graph_probs[series][i - min]= hetProbs[i];
      }

   prob_graph.stacked = true;
   prob_graph.SetCategoricalDataCounts(graph_probs);
   prob_graph.SetDataValueTag(0, observed_hets, "*");
   prob_graph.useLegend = false;

   prob_graph.xAxis.minMin = min - 3;
   prob_graph.xAxis.maxMax = max + 3;

   prob_graph.yAxis.label = "Probability";
   prob_graph.xAxis.label = "Number of Heterozygotes";
   prob_graph.SetSeriesColor(0, 1.0, 0.0, 0.0);
   prob_graph.SetSeriesColor(1, 0.0, 0.0, 1.0);

   double x0 = 0.99 * pdf.page.GetWidth() * (0.33);
   double y0 = 0.9 * pdf.page.GetHeight() * (0.35);
   double x1 = 0.99 * pdf.page.GetWidth() * (1.0 );
   double y1 = 0.9 * pdf.page.GetHeight() * ((0.0 + 1.0)/1.0) - 0.35;

   prob_graph.DrawInBox(pdf, x0, y0, x1, y1);
   }

void PedstatsHWE::MakeSignificanceKey(PDF & pdf, double p_value)
   {
   double page_width = pdf.page.GetWidth();
   double page_height = pdf.page.GetHeight();

   double box_width = 0.25 * page_width;
   double box_height = 0.15 * page_height;

   double x0 = 0.39 * page_width;
   double y0 = 0.15 * page_height;
   double y1 = y0 + 0.15 * page_height;

   double pt_1 = x0 + 0.2 * box_width;
   double pt_2 = x0 + 0.5 * box_width;
   double pt_3 = x0 + 0.8 * box_width;

   pdf.page.SetFontSize(8);
   pdf.page.hTextAlignment = taCenter;
   pdf.page.WriteText(x0 + 0.1 * box_width, y1 - 0.15 * box_height, "Threshold");
   pdf.page.WriteText(x0 + 0.1 * box_width, y1 - 0.2 * box_height, "-----------------");
   pdf.page.WriteText(pt_1, y1 - 0.4 * box_height, "p < 0.05");
   pdf.page.WriteText(pt_2, y1 - 0.4 * box_height, "p < 0.01");
   pdf.page.WriteText(pt_3, y1 - 0.4 * box_height, "p < 0.001");

   double pt_y3 = y1 - 0.6 * box_height;
   double s_w = 0.1 * box_width;
   double s_h = 0.1 * box_height;

   p_value > 0.001 ? DrawGreenChecks(pdf, pt_3, pt_y3, s_w, s_h) : DrawRedCross(pdf, pt_3, pt_y3, s_w, s_h);
   p_value > 0.01  ? DrawGreenChecks(pdf, pt_2, pt_y3, s_w, s_h) : DrawRedCross(pdf, pt_2, pt_y3, s_w, s_h);
   p_value > 0.05  ? DrawGreenChecks(pdf, pt_1, pt_y3, s_w, s_h) : DrawRedCross(pdf, pt_1, pt_y3, s_w, s_h);
   }

void PedstatsHWE::DrawGreenChecks(PDF & pdf, double x, double y, double width, double height)
   {
   pdf.page.SetLineColor(0.0, 1.0, 0.0);
   double line_width = 4.0;

   pdf.page.DrawCheck(x, y, width, height, 1.0, line_width);
   pdf.page.DrawCheck(x , y - 0.1 * height, width, height, 0.5, line_width);
   }

void PedstatsHWE::DrawRedCross(PDF & pdf, double x, double y, double width, double height)
   {
   pdf.page.SetLineColor(1.0, 0.0, 0.0);
   pdf.page.SetLineWidth(4.0);
   pdf.page.DrawLine(x - width * 0.5, y - height * 0.5, x + width * 0.5, y + height * 0.5);
   pdf.page.DrawLine(x - width * 0.5, y + height * 0.5, x + width * 0.5, y - height * 0.5);
   }



 
 
 
 
