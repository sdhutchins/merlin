////////////////////////////////////////////////////////////////////// 
// pedstats/Pedstats.cpp 
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
 
#include "Pedstats.h"
#include "PedstatsPDF.h"

String pedName = "pedstats.ped";
String dataName = "pedstats.dat";
String ibdName = "pedstats.ibd";
String pdfName = "pedstats.pdf";
String markerOutputName = "pedstats.markerinfo";

FilteredPedigree ped;
IBDTable ibd;
PedigreePairs all_pairs(ped);
FilteredFamilyStats family_stats(ped);
PedstatsHWE hwe(ped);
PedstatsQuality qual(ped);

bool hardy = false, hardyAll = false, hardyFounders = false, hardyUnrelated = false;
bool markerTables = false;

bool zeroFounderParents  = false;
bool verboseGroups = false;

bool nocheck = false;
bool pairs = false;
bool byGender = false;
bool byFamily = false;
bool byAffection = false;
bool trimPedigree = false;
bool rewritePedigree = false;

bool writePDF = false;
bool writeFamilyPDF = false;
bool writeTraitPDF = false;
bool writeMarkerPDF = false;
bool writeAffPDF = false;

int trait_filter = 0;
int covariate_filter = 0;
int affection_filter = 0;

int markersReported = 50;

int main(int argc, char * argv[])
   {
   #ifndef PEDVERSION
      printf("Pedigree Statistics\n(c) 1999-2006 Goncalo Abecasis, 2002-2006 Jan Wigginton\n");
   #else
      printf("Pedigree Statistics - " PEDVERSION "\n(c) 1999-2006 Goncalo Abecasis, 2002-2006 Jan Wigginton\n");
   #endif

    ReadParameters(argc, argv);

   ped.Prepare(dataName);
   ped.Load(pedName);
   ped.LumpAlleles(0.0,false);
   
   if (trimPedigree) ped.Trim();
   ped.ApplyFilters();

   ibd.Load(ped, ibdName);
   all_pairs.Build();

   if (PedstatsAgeCheck::ageLabel != "") CheckAges(all_pairs.sibs);

   family_stats.AccumulateStatistics();

   CheckForLogicalSummaries(family_stats);

   FamilyStatistics(family_stats);
   TraitStatistics(stdout, ped, all_pairs.sibs);
   CovarStatistics(stdout, ped, all_pairs.sibs);
   AffectionStatistics(stdout, ped);

   if (pairs)    PairStatistics(stdout, ped, all_pairs);
   if (!nocheck) RunInheritanceChecks();

   PDF * pdf = NULL;

   if (writeTraitPDF || writeFamilyPDF || writeMarkerPDF || writeAffPDF)
      {
      pdf = new PDF();
      pdf->OpenFile((const char *) pdfName);
      pdf->page.SetSize(psLetterR);

      if (writeFamilyPDF) GraphFamilyStatistics(*pdf, family_stats);
      if (writeTraitPDF)  GraphTraitStatistics(*pdf);
      if (writeTraitPDF)  GraphCovariateStatistics(*pdf);
      if (writeAffPDF)    GraphAffectionStatistics(*pdf);
      }

   MarkerSummaries(pdf);
   if (byFamily) ByFamilyStatistics();
   IBDStatistics();

   FileReport(pdf != NULL);

   if (pdf != NULL)
      delete pdf;

    return 0;
   };

BEGIN_LONG_PARAMETERS(longParameters)
   LONG_PARAMETER_GROUP("Pedigree File")
      LONG_PARAMETER("ignoreMendelianErrors", &nocheck)
      LONG_PARAMETER("chromosomeX", &Pedigree::chromosomeX)
      LONG_PARAMETER("trim", &trimPedigree)
   LONG_PARAMETER_GROUP("Hardy-Weinberg")
      LONG_PARAMETER("hardyWeinberg", &hardy)
      LONG_PARAMETER("showAll", &PedstatsHWE::showAll)
      LONG_DOUBLEPARAMETER("cutoff", &PedstatsHWE::significanceCutoff)
   LONG_PARAMETER_GROUP("HW Sample")
      LONG_PARAMETER("checkFounders", &hardyFounders)
      LONG_PARAMETER("checkAll", &hardyAll)
      LONG_PARAMETER("checkUnrelated", &hardyUnrelated)
   LONG_PARAMETER_GROUP("Output")
      LONG_PARAMETER("pairs", &pairs)
      LONG_PARAMETER("rewritePedigree", &rewritePedigree)
      LONG_PARAMETER("markerTables", &markerTables)
      LONG_PARAMETER("verbose", &verboseGroups)
   LONG_PARAMETER_GROUP("Grouping")
      LONG_PARAMETER("bySex", &byGender)
      LONG_PARAMETER("byFamily", &byFamily)
   LONG_PARAMETER_GROUP("Age Checking")
      LONG_STRINGPARAMETER("age", &PedstatsAgeCheck::ageLabel)
      LONG_STRINGPARAMETER("birth", &PedstatsAgeCheck::birthLabel)
   LONG_PARAMETER_GROUP("Generations")
      LONG_DOUBLEPARAMETER("minGap", &PedstatsAgeCheck::minGap)
      LONG_DOUBLEPARAMETER("maxGap", &PedstatsAgeCheck::maxGap)
      LONG_DOUBLEPARAMETER("sibGap", &PedstatsAgeCheck::sibGap)
   LONG_PARAMETER_GROUP("PDF Options")
      LONG_PARAMETER("pdf", &writePDF)
      LONG_PARAMETER("familyPDF", &writeFamilyPDF)
      LONG_PARAMETER("traitPDF", &writeTraitPDF)
      LONG_PARAMETER("affPDF", &writeAffPDF)
      LONG_PARAMETER("markerPDF", &writeMarkerPDF)
   LONG_PARAMETER_GROUP("Filter")
      LONG_INTPARAMETER("minGenos", &FilteredPedigree::minGenos)
      LONG_INTPARAMETER("minPhenos", &FilteredPedigree::minPhenos)
      LONG_INTPARAMETER("minCovariates", &FilteredPedigree::minCovariates)
      LONG_STRINGPARAMETER("affectedFor", &FilteredPedigree::affectedFor)
   BEGIN_LEGACY_PARAMETERS()
       LONG_DOUBLEPARAMETER("minGeneration", &PedstatsAgeCheck::minGap)
       LONG_DOUBLEPARAMETER("maxGeneration", &PedstatsAgeCheck::maxGap)
       LONG_PARAMETER("zeroFounders", &zeroFounderParents)
   END_LONG_PARAMETERS();

void ReadParameters(int argc, char * argv[])
   {
   ParameterList pl;

   pl.Add(new StringParameter('p', "Pedigree File", pedName));
   pl.Add(new StringParameter('d', "Data File", dataName));
   pl.Add(new StringParameter('i', "IBD File", ibdName));
   pl.Add(new StringParameter('a', "Adobe PDF File", pdfName));
   pl.Add(new StringParameter('x', "Missing Value Code", Pedigree::missing));
   pl.Add(new LongParameters("Additional Options", longParameters));

   pl.Read(argc, argv);

   if (hwe.showAll)
      hwe.significanceCutoff = 1.0;
   
   hardyUnrelated |= hardy;
   hardyAll |= hardy;

   pl.Status();
   }

void RunInheritanceChecks()
   {
   // First to the regular check, which will stop if a problem is found
   ped.InheritanceCheck();

   bool pass = true;

   // Then do the more thorough check for any non-trivial pedigrees
   for (int m = 0; m < ped.markerCount; m++)
      for (int i = 0; i < ped.familyCount; i++)
         if ((ped.families[i]->count > ped.families[i]->founders))
            pass &= GenotypeList::EliminateGenotypes(ped, ped.families[i], m);

   //If the second level checks were not successful, then stop
   if (!pass)
      error("Mendelian inheritance errors detected\n");
   }

void CheckForLogicalSummaries(FilteredFamilyStats & family_stats)
   {
   writeTraitPDF |= writePDF;
   writeFamilyPDF |= writePDF;
   writeAffPDF |= writePDF;
   writeMarkerPDF |= (writePDF && ped.markerCount <= markersReported);

   writeAffPDF &= (Person::affectionCount > 0);
   writeMarkerPDF &= (ped.markerCount > 0);
   writeTraitPDF &= (Person::traitCount > 0 || Person::covariateCount > 0);
   hardyAll |= hardy;
   hardyUnrelated |= hardy;

   if (family_stats.nonFounderCount == 0 && hardy)
      hardy = hardyFounders = false;

   markerTables |= (ped.markerCount > markersReported);
   }

void FamilyStatistics(FilteredFamilyStats & stats)
   {
   if (stats.familiesKept == 0)
      {
      printf("WARNING - Pedigree is empty %s", ped.filterOn ? "after filtering" : "");
      return;
      }

   IntArray counts(stats.maxMembers + 1);
   IntArray gens(stats.maxGenerations + 1);
   IntArray female_cts(stats.maxFemales + 1), male_cts(stats.maxMales + 1);

   counts.Zero();
   gens.Zero();
   female_cts.Zero();
   male_cts.Zero();

   for (int f = 0; f < stats.familiesKept; f++)
      {
      counts[stats.familyMembers[f]]++;
      gens[stats.familyGenerations[f]]++;
      female_cts[stats.familyFemales[f]]++;
      male_cts[stats.familyMales[f]]++;
      }

   IntArray idxcount(stats.maxMembers), idxgen(stats.maxGenerations + 1);
   IntArray idx_females(stats.maxFemales + 1), idx_males(stats.maxMales + 1);

   MakeIndex(idxcount, counts);
   MakeIndex(idxgen, gens);
   MakeIndex(idx_females, female_cts);
   MakeIndex(idx_males, male_cts);

   double fdiv = max((double) stats.familiesKept, 1.0) / (double) 100.0;

   printf("\nPEDIGREE STRUCTURE %s\n"
          "==================%s\n"
          "   Individuals: %d \n"
          "      Founders: %d founders, %d nonfounders\n"
          "        Gender: %d females, %d males\n",
          ped.filterOn ?  "(AFTER FILTERING)" : "",
          ped.filterOn ?  "=================" : "",
          stats.founderCount + stats.nonFounderCount, stats.founderCount, stats.nonFounderCount,
          stats.femaleCount, stats.maleCount);

   printf("      Families: %d\n\n"
          "  Family Sizes\n"
          "       Average: %.2f (%d to %d)\n"
          "  Distribution: %d (%.1f%%), %d (%.1f%%) and %d (%.1f%%)\n\n"
          "  Generations%s\n"
          "       Average: %.2f (%d to %d)\n"
          "  Distribution: %d (%.1f%%), %d (%.1f%%) and %d (%.1f%%)\n\n",
          stats.familiesKept,
          stats.averageMembers, stats.minMembers, stats.maxMembers,
          idxcount[0], (double) counts[idxcount[0]] / fdiv,
          idxcount[1], (double) counts[idxcount[1]] / fdiv,
          idxcount[2], (double) counts[idxcount[2]] / fdiv,
          ped.filterOn ? "*" : "",
          stats.averageGenerations, stats.minGenerations, stats.maxGenerations,
          idxgen[0], (double) gens[idxgen[0]] / fdiv,
          idxgen[1], (double) gens[idxgen[1]] / fdiv,
          idxgen[2], (double) gens[idxgen[2]] / fdiv);

   if (byGender)
      printf( " Females Per Family\n"
          "       Average: %.2f (%d to %d)\n"
          "  Distribution: %d (%.1f%%), %d (%.1f%%) and %d (%.1f%%)\n\n"
          "  Males Per Family\n"
          "       Average: %.2f (%d to %d)\n"
          "  Distribution: %d (%.1f%%), %d (%.1f%%) and %d (%.1f%%)\n\n",
          stats.averageFemales, stats.minFemales, stats.maxFemales,
          idx_females[0], (double) female_cts[idx_females[0]] / fdiv,
          idx_females[1], (double) female_cts[idx_females[1]] / fdiv,
          idx_females[2], (double) female_cts[idx_females[2]] / fdiv,
          stats.averageMales, stats.minMales, stats.maxMales,
          idx_males[0], (double) male_cts[idx_males[0]] / fdiv,
          idx_males[1], (double) male_cts[idx_males[1]] / fdiv,
          idx_males[2], (double) male_cts[idx_males[2]] / fdiv);

   if (ped.filterOn)
      printf("\n* NOTE: Generation structure reported for filtered data may include\n"
              "  filtered individuals\n\n\n");

   CheckFamilyConnections();
   if (rewritePedigree)
      RewritePedigree();
   }

void CheckFamilyConnections()
   {
   printf("Checking family connectedness ... \n");
   int connected = true;

   IntArray group_membership;

   for (int f = 0, groups; f < ped.familyCount; f++)
      {
      group_membership.Clear();

      if ( (groups = ped.families[f]->ConnectedGroups(&group_membership)) != 1 )
         {
         printf(" Family %s (%d individuals) consists of %d disjoint groups\n",
          (const char *) ped.families[f]->famid, ped.families[f]->count, groups),
         connected = false;
         if (verboseGroups)
            {
            for (int g = 1; g <= groups; g++)
               {
               printf("\n\tGroup %d for family %s includes\n", g, (const char *) ped.families[f]->famid);

               for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
                  if (group_membership[ped[i].traverse] == g)
                     printf("\t\tPerson %s\n", (const char *) ped[i].pid);
               }
            printf("\n");
            }
         }
      }

   printf(connected ? "   All individuals in each family are connected.\n\n" :
                      "  \nUse the --rewritePedigree option to break families"
                      " into connected groups\n");
   }

void RewritePedigree()
   {
   printf("Rewriting data file as [pedstats.dat]... ");
   ped.WriteDataFile("pedstats.dat");

   printf("Done!\n"
             "Rewriting pedigree file as [pedstats.ped]... ");

   FILE * pedfile = fopen("pedstats.ped", "wt");

   if (zeroFounderParents)
      ZeroFounderParents(ped);

   IntArray groups;
   for (int f = 0; f < ped.familyCount; f++)
      {
      ped.families[f]->ConnectedGroups(&groups);

      for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         {
         String orig_id = ped[i].pid;

         ped.WritePerson(pedfile, i, (const char *) (ped[i].famid + "_" + groups[ped[i].traverse]),
                         orig_id);
         }
      }

   printf("Done!\n");
   fclose(pedfile);
   }

void TraitStatistics(FILE * file, FilteredPedigree & ped, SimplePairList & pairs,
 const char * added_label)
   {
   if (Person::traitCount == 0) return;

   String extra = added_label != NULL ? added_label : "";
   if (ped.filterOn)
      extra += " *";

   String fill('=', extra.Length() + 29);
   fprintf(file, "\n\nQUANTITATIVE TRAIT STATISTICS%s\n%s\n",
                 (const char *) extra, (const char *) fill);

   GetTraitStatistics(file, ped, pairs, SEX_UNKNOWN, "ALL DATA:", false);

   if (byGender)
      {
      GetTraitStatistics(file, ped, pairs, SEX_FEMALE, "FEMALES:", false);
      GetTraitStatistics(file, ped, pairs, SEX_MALE, "MALES:", false);
      }

   if (ped.filterOn)
      fprintf(file, "\n* NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);
   }

void CovarStatistics(FILE * file, FilteredPedigree & ped, SimplePairList & pairs, const char * added_label)
   {
   if (Person::covariateCount == 0) return;

   String extra = added_label != NULL ? added_label : "";
   if (ped.filterOn)
      extra += " *";
   String fill('=', extra.Length() + 20);

   fprintf(file, "\n\nCOVARIATE STATISTICS%s\n%s\n",(const char *) extra, (const char *) fill);

   GetTraitStatistics(file, ped, pairs, SEX_UNKNOWN, "ALL DATA:", true);

   if (byGender)
      {
      GetTraitStatistics(file, ped, pairs, SEX_FEMALE, "FEMALES:", true);
      GetTraitStatistics(file, ped, pairs, SEX_MALE, "MALES:", true);
      }

   if (ped.filterOn)
      fprintf(file, "\n* NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);
   }

void GetTraitStatistics(FILE * file, FilteredPedigree & ped, SimplePairList & sib_pairs,
 int sex_filter, const char * label, bool is_covariate)
   {
   if (byGender)
      fprintf(file, "\n%-15s\n", label);

   fprintf(file, "\n%11s %19s %8s %8s %8s %8s %8s\n",
           "", "[All Phenotypes]", "Min", "Max", "Mean", "Var", "SibCorr");

   int count = 0, possible = 0, fcount = 0, fpossible = 0;
   int filter_limit = is_covariate ? Person::covariateCount : Person::traitCount;

   StringArray founder_stats(filter_limit);

   for (int t = 0; t < filter_limit; t++)
      {
      trait_filter = t;
      int tcount = 0, tpossible = 0, tfcount = 0, tfpossible = 0;

      double mean = 0.0, var = 0.0;
      double min_value = 1e300, max_value = -1e300;

      double mean_f = 0.0, var_f = 0.0;
      double min_fvalue = 1e300, max_fvalue = -1e300;

      bool typed;
      double val;
      for (int i = 0; i < ped.count; i++)
         {
         if (ped.filtered[i]) continue;

         if (sex_filter != SEX_UNKNOWN && ped[i].sex != sex_filter) continue;

         tpossible ++;
         if (ped[i].isFounder()) tfpossible ++;

         typed = is_covariate ? ped[i].isControlled(t) : ped[i].isPhenotyped(t);
         if (!typed ) continue;

         val = is_covariate ? ped[i].covariates[t] : ped[i].traits[t];

         min_value = min(min_value, val);
         max_value = max(max_value, val);
         mean += val;
         var += val * val;

         if (ped[i].isFounder())
            {
            tfcount ++;
            min_fvalue = min(min_fvalue, val);
            max_fvalue = max(max_fvalue, val);
            mean_f += val;
            var_f += val * val;
            }
         tcount ++;
         }

      possible += tpossible;
      fpossible += tfpossible;
      count += tcount;
      fcount += tfcount;

      var = (tcount <= 1) ? 0.0 : (double) (tcount * var - mean * mean)
                          /(double)((tcount) * (tcount -1.0));
      mean /= max(tcount, 1);

      var_f = (tfcount <= 1) ? 0.0 : (double) (tfcount * var_f - mean_f * mean_f)
                             /(double)((tfcount) * (tfcount - 1.0));
      mean_f /= max(tfcount, 1);

      int num_pairs = 0;
      double r_sib;

      if (sex_filter != SEX_UNKNOWN)
         {
         SimplePairList filtered_sibs(ped);
         sex_filter == SEX_MALE ?
                       sib_pairs.FilterOnBothMale(filtered_sibs) :
                       sib_pairs.FilterOnBothFemale(filtered_sibs);

         r_sib = filtered_sibs.Correlation(t, is_covariate, true, num_pairs);
         }
      else
         r_sib = sib_pairs.Correlation(t, is_covariate, true, num_pairs);

      String trait_name;
      trait_name = is_covariate ? Person::covariateNames[t] : Person ::traitNames[t];

      if (num_pairs > 1)
         fprintf(file, "%15.15s %8d %5.1f%% %8.*f %8.*f %8.*f %8.*f %8.3f\n",
                 (const char * ) trait_name, tcount,
                 (double) tcount / (double)  max(tpossible, 1) * 100,
                 fabs(min_value) < 1e4 ? 3 : 0, min_value,
                 fabs(max_value) < 1e4 ? 3 : 0, max_value,
                 fabs(mean) < 1e4 ? 3 : 0, mean,
                 var < 1e4 ? 3 : 0, var, r_sib);
      else if (tcount > 0)
         fprintf(file, "%15.15s %8d %5.1f%% %8.*f %8.*f %8.*f %8.*f %8s\n",
                 (const char * ) trait_name, tcount,
                 (double) tcount / (double)  max(tpossible, 1) * 100,
                  fabs(min_value) < 1e4 ? 3 : 0, min_value,
                  fabs(max_value) < 1e4 ? 3 : 0, max_value,
                  fabs(mean) < 1e4 ? 3 : 0, mean,
                  var < 1e4 ? 3 : 0, var, "-");
      else
         fprintf(file, "%15.15s %8d %5.1f%% %8s %8s %8s %8s %8s\n",
                 (const char * ) trait_name, tcount,
                 (double) tcount / (double)  max(tpossible, 1) * 100,
                 "-", "-", "-", "-", "-");

      if (tfcount > 0)
         founder_stats[t].printf("%15.15s %8d %5.1f%% %8.*f %8.*f %8.*f %8.*f %8s",
            (const char * ) trait_name, tfcount,
            (double) tfcount / (double)  max(tfpossible, 1) * 100,
            fabs(min_fvalue) < 1e4 ? 3 : 0, min_fvalue,
            fabs(max_fvalue) < 1e4 ? 3 : 0, max_fvalue,
            fabs(mean_f) < 1e4 ? 3 : 0, mean_f,
            var_f < 1e4 ? 3 : 0, var_f, "-");
      else
         founder_stats[t].printf("%15.15s %8d %5.1f%% %8s %8s %8s %8s %8s",
            (const char * ) trait_name, tfcount,
            (double) tfcount / (double)  max(tfpossible, 1) * 100,
            "-", "-", "-", "-", "-");
      }

   fprintf(file, "%15.15s %8d %5.1f%%\n\n", "Total",
           count, (double) count / (double) max(possible, 1) * 100);

   fprintf(file, "\n%10s %20s %8s %8s %8s %8s %8s\n",
          "", "[Founders Only]", "Min", "Max", "Mean", "Var", "SibCorr");

   for (int i = 0; i < filter_limit; i++)
      fprintf(file, "%60s\n", (const char *) founder_stats[i]);

   fprintf(file, "%15.15s %8d %5.1f%%\n\n", "Total",
           fcount, (double) fcount / (double) max(fpossible, 1) * 100);
   }

void AffectionStatistics(FILE * file, FilteredPedigree & ped, const char * added_label)
   {
   if (Person::affectionCount == 0) return;

   String extra = added_label != NULL ? added_label : "";
   if (ped.filterOn)
      extra += " *";

   String fill('=', extra.Length() + 20);
   fprintf(file, "\n\nAFFECTION STATISTICS%s\n%s\n",
          (const char *) extra, (const char *) fill);

   GetAffectionStatistics(file, ped, SEX_UNKNOWN, "ALL DATA:");

   if (byGender)
      {
      GetAffectionStatistics(file, ped, SEX_FEMALE, "FEMALES:");
      GetAffectionStatistics(file, ped, SEX_MALE, "MALES:");
      }

   if (ped.filterOn)
      fprintf(file, "\n* NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);
   }

void GetAffectionStatistics(FILE * file, FilteredPedigree & ped, int sex_filter, const char * label)
   {
   if (byGender) fprintf(file, "\n%-15s", label);
   fprintf(file, "\n%15s %15s %15s %10s\n", "", "[Diagnostics]", "[Founders]", "Prevalence");

   int count = 0, possible = 0, fcount = 0, fpossible = 0;

   for (int a = 0; a < Person::affectionCount; a++)
      {
      int tcount = 0, tpossible = 0, tfcount = 0, tfpossible = 0;
      int affecteds = 0;

      for (int i = 0; i < ped.count; i++)
         {
         if (ped.filtered[i]) continue;

         if (sex_filter != SEX_UNKNOWN && ped[i].sex != sex_filter) continue;
         tpossible ++;
         tfpossible += ped[i].isFounder();

         if (!ped[i].isDiagnosed(a)) continue;

         affecteds += ped[i].affections[a] == 2;
         tfcount += ped[i].isFounder();
         tcount ++;
         }

      possible += tpossible;
      fpossible += tfpossible;
      count += tcount;
      fcount += tfcount;

      fprintf(file, "%15.15s %8d %5.1f%% %8d %5.1f%% %9.1f%%\n",
              (const char *) Person::affectionNames[a],
              tcount, (double) tcount / (double)  max(tpossible, 1) * 100,
              tfcount, (double) tfcount / (double) max(tfpossible, 1) * 100,
              (double) affecteds / (double) max(tcount, 1) * 100 );
      }

   fprintf(file, "%15s %8d %5.1f%% %8d %5.1f%%\n\n", "Total",
           count, (double) count / (double) max(possible, 1) * 100,
           fcount, (double) fcount / (double) max(fpossible, 1) * 100 );
   }

void MarkerStatistics(FILE * file, FilteredPedigree & ped, const char * added_label)
   {
   if (Person::markerCount == 0) return;

   String extra = added_label != NULL ? added_label : "";
   if (ped.filterOn)
      extra += " *";
   String fill('=', extra.Length() + 26);
   fprintf(file, "\n\nMARKER GENOTYPE STATISTICS%s\n%s\n",(const char *) extra, (const char *) fill);

   GetMarkerStatistics(file, ped, SEX_UNKNOWN, "ALL DATA:");

   if (byGender)
      {
      GetMarkerStatistics(file, ped, SEX_FEMALE, "FEMALES:");
      GetMarkerStatistics(file, ped, SEX_MALE, "MALES:");
      }

   if (ped.filterOn)
      fprintf(file, "\n* NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);

   fprintf(file, "Total markers: %d\n", Person::markerCount);
   }

void GetMarkerStatistics(FILE * file, FilteredPedigree & ped, int sex_filter, const char * label)
   {
   if (byGender) fprintf(file, "%-15s", label);

   fprintf(file, "\n%15s %15s %15s %10s\n","", "[Genotypes]", "[Founders]", "Hetero");

   int genotypes = 0, possible = 0, fgenotypes = 0, fpossible = 0,
       hetero = 0, maxHetero = 0;

   for (int m = 0; m < Person::markerCount; m++)
      {
      int mGenotypes = 0, mPossible = 0, mfGenotypes = 0, mfPossible = 0,
          mHetero = 0, mMaxHetero = 0;
      int m0, m1;
      bool is_founder;

      for (int i = 0; i < ped.count; i++)
         {
         if (ped.filtered[i]) continue;
         if (sex_filter != SEX_UNKNOWN && ped[i].sex != sex_filter) continue;

         is_founder = ped[i].isFounder();

         mPossible ++;
         mfPossible += is_founder;

         m0 = ped[i].markers[m].one;
         m1 = ped[i].markers[m].two;

         if (m0 == 0 || m1 == 0) continue;

         mGenotypes ++;
         mfGenotypes += is_founder;

         if (ped.chromosomeX && ped[i].sex == SEX_MALE) continue;
         if (m0 != m1)
            mHetero ++;
         mMaxHetero ++;
         }

      possible += mPossible;
      fpossible += mfPossible;
      fgenotypes += mfGenotypes;
      genotypes += mGenotypes;
      hetero += mHetero;
      maxHetero += mMaxHetero;

      fprintf(file, "%15s %8d %5.1f%% %8d %5.1f%% %9.1f%%\n",
              (const char *) Person::markerNames[m],
              mGenotypes, (double) mGenotypes / (double) max(mPossible, 1) * 100,
               mfGenotypes, (double)mfGenotypes / (double)max(mfPossible, 1) * 100,
              (double) mHetero / (double) max(mMaxHetero, 1) * 100);
      }

   fprintf(file, "%15s %8d %5.1f%% %8d %5.1f%% %9.1f%%\n%s\n", "Total",
           genotypes, (double) genotypes / (double) max(possible, 1) * 100,
           fgenotypes, (double) fgenotypes / (double) max(fpossible, 1) * 100,
           (double) hetero / (double) max(maxHetero, 1) * 100,
           ped.chromosomeX ? "\nNOTE: Heterozygosity calculation excludes males\n" : "");
   }

void IBDStatistics()
   {
   if (Person::markerCount == 0) return;
   if (ibd.isEmpty()) return;

   String fill('=', 33);
   printf("\n\nIBD FAMILY INFORMATION (FROM FILE)\n%s\n\n", (const char *) fill);

   printf("%17s%10s %10s %16s %16s\n",
          " ", "[Singletons]", "[Trios]", "[2 generations]", "[3+ generations]");

   for (int m = 0; m < Person::markerCount; m++)
      {
      int ibdTwoGen = 0, ibdPossTwoGen = 0;
      int ibdTwoPlusGen = 0, ibdPossTwoPlusGen = 0;
      int ibdTrios = 0, ibdPossTrios = 0;
      int ibdSingletons = 0, ibdPossSingletons = 0;

      for (int f = 0; f < ped.familyCount; f++)
         {
         bool have_family = ibd.HaveFamily(m, ped.families[f]) && !ibd.isEmpty();

         if (ped.families[f]->count == 3 && ped.families[f]->founders == 2)
            {
            if (have_family) ibdTrios ++;
            ibdPossTrios ++;
            }
         else if (ped.families[f]->count == ped.families[f]->founders)
            {
            if (have_family) ibdSingletons ++;
            ibdPossSingletons ++;
            }

         int generations = 0;
         for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
            generations = max(family_stats.Generation(ped.persons[i]) + 1, generations);

         if (generations == 2)
            {
            if (have_family) ibdTwoGen ++;
            ibdPossTwoGen ++;
            }
         else if (generations > 2)
            {
            if (have_family) ibdTwoPlusGen ++;
            ibdPossTwoPlusGen ++;
            }
         }
      printf("%15.15s %6d/%-6d %6d/%-6d %7d/%-8d %7d/%-8d\n",
             (const char *) Person::markerNames[m],
             ibdSingletons, ibdPossSingletons, ibdTrios, ibdPossTrios,
             ibdTwoGen, ibdPossTwoGen, ibdTwoPlusGen, ibdPossTwoPlusGen);
     }

   printf("\n\nNOTE: All counts are number of families per category that have some\n"
          "  marker IBD information available.\n\n\n");

   if (ped.filterOn)
      printf("\n* NOTE:  Generation structure reported includes filtered individuals\n\n");
   }

void MakeIndex(IntArray & index, IntArray & counts)
   {
   index.Dimension(counts.Length());
   index.SetSequence(0, 1);

   // Insertion sort should be OK, for N < 20
   for (int i, j = 1; j < counts.Length(); j++)
      {
      int tmp = counts[j];  // index[j] = j
      for (i = j - 1; (i >= 0) && (counts[index[i]] < tmp); i--)
         index[i + 1] = index[i];
      index[i + 1] = j;
      }

   // Make sure we have at least three items
   for (int i = counts.Length() + 1; i < 4; i++)
      {
      counts.Dimension(i);
      index.Dimension(i);

      counts[i - 1] = 0;
      index[i - 1] = i - 1;
      }
   }

void PairStatistics(FILE * file, FilteredPedigree & ped, PedigreePairs & pairs,
 const char * added_label, bool traits_only)
   {
   String extra = added_label != NULL ? added_label : "";
   if (ped.filterOn)
      extra += " *";

   String fill('=', extra.Length() + 15);
   fprintf(file, "\n\nPAIR STATISTICS%s\n%s\n", (const char *) extra, (const char *) fill);

   GetPairStatistics(file, pairs, "ALL DATA:", traits_only);

   PedigreePairs female_prs(ped), male_prs(ped), opp_prs(ped);
   pairs.SplitOnGender(female_prs, male_prs, opp_prs);

   if (byGender)
      {
      GetPairStatistics(file, female_prs, "SAME SEX PAIRS (FEMALE):", traits_only);
      GetPairStatistics(file, male_prs, "SAME SEX PAIRS (MALE):", traits_only);
      GetPairStatistics(file, opp_prs, "OPPOSITE SEX PAIRS:", traits_only);
      }

   if (ped.filterOn)
      fprintf(file, "\n* NOTE:  Sample restricted to %s\n", (const char *) ped.textFilterLabel);
   }

void GetPairStatistics(FILE * file, PedigreePairs & pairs, const char * label, bool traits_only)
   {
   if (byGender)
      fprintf(file, "\n%-15s %s\n", label, pairs.totalCt > 0 ? "" : "NONE" );
   else if (pairs.totalCt <= 0)
      fprintf(file, "\n\n%s\n", "NO RELATIVE PAIRS");

   if (pairs.totalCt <= 0)
      return;

   if (!traits_only)
      pairs.PrintCounts(file);

   PairTraitStatistics(file, pairs, false);
   PairTraitStatistics(file, pairs, true);

   if (!traits_only)
      PairAffectionStatistics(file, pairs);
   }

void PairAffectionStatistics(FILE * file, PedigreePairs & pairs)
   {
   if (Person::affectionCount == 0) return;

   fprintf(file, "\n\nPair Counts by Affection Status:\n");

   fprintf(file, "%15s %8s %8s %8s %12s %12s %9s", " ",
          "Sib", "HalfSib","Cousin", "ParentChild", "Grandparent", "Avuncular");

   for (affection_filter = 0; affection_filter < Person::affectionCount; affection_filter++)
      {
      IntArray unaff, disc, aff;

      pairs.CountAffectionTypes(affection_filter, unaff, disc, aff);

      fprintf(file, "\n%15.15s\n", (const char *) ped.affectionNames[affection_filter]);

      fprintf(file, "%2s[%11s] %8d %8d %8d %12d %12d %9d\n", " ", "Unaffected",
              unaff[0], unaff[1], unaff[2], unaff[3], unaff[4], unaff[5]);

      fprintf(file, "%2s[%11s] %8d %8d %8d %12d %12d %9d\n", " ", "Discordant",
              disc[0], disc[1], disc[2], disc[3], disc[4], disc[5]);

      fprintf(file, "%2s[%11s] %8d %8d %8d %12d %12d %9d\n", " ", "Affected",
              aff[0], aff[1], aff[2], aff[3], aff[4], aff[5]);
      }

   if (pairs.sibs.gender == "Opposite")
      {
      bool header_printed = false;

      for (affection_filter = 0; affection_filter < ped.affectionCount; affection_filter++)
         {
         PedigreePairs mf_pairs(ped), fm_pairs(ped);

         if ( pairs.avuncularCt + pairs.grandParentCt + pairs.parentCt <= 0 ) continue;

         pairs.SplitOnGenderOrder(mf_pairs, fm_pairs);

         IntArray mf_unaff, mf_disc, mf_aff;
         IntArray fm_unaff, fm_disc, fm_aff;

         mf_pairs.CountAffectionTypes(affection_filter, mf_unaff, mf_disc, mf_aff, false);
         fm_pairs.CountAffectionTypes(affection_filter, fm_unaff, fm_disc, fm_aff, false);

         int total_diagnosed = 0;
         for (int i = 3; i < 6; i++)
             total_diagnosed += mf_unaff[i] + mf_disc[i] + mf_aff[i]
                             + fm_unaff[i] + fm_disc[i] + fm_aff[i];
         if (total_diagnosed == 0) continue;

         if (header_printed == false) PrintDetailHeader(file);
         header_printed = true;

         fprintf(file, "%15.15s\n", (const char *) Person::affectionNames[affection_filter]);

         fprintf(file, "%2s[%11s] %8d %8d %8d %8d %8d %8d\n", " ", "Unaffected",
           mf_unaff[3], fm_unaff[3], mf_unaff[4], fm_unaff[4], mf_unaff[5], fm_unaff[5]);

         fprintf(file, "%2s[%11s] %8d %8d %8d %8d %8d %8d\n", " ", "Discordant",
           mf_disc[3], fm_disc[3], mf_disc[4], fm_disc[4], mf_disc[5], fm_disc[5]);

         fprintf(file, "%2s[%11s] %8d %8d %8d %8d %8d %8d\n", " ", "Affected",
           mf_aff[3], fm_aff[3], mf_aff[4], fm_aff[4], mf_aff[5], fm_aff[5]);
         }
      }
   fprintf(file, "\n");
   }

void PrintDetailHeader(FILE * file, const char * additional)
   {
   fprintf(file, "\nOrdered Pair Detail %s:\n", additional ? additional : "");

   fprintf(file, "%15s %17s %17s %17s\n", " ",
           "ParentChild", "Grandparent", "Avuncular");

   fprintf(file, "%15s %8s %8s %8s %8s %8s %8s\n"," ",
           "F/D", "M/S", "GF/GD", "GM/GS", "U/N", "A/N");
   }

void PairTraitStatistics(FILE * file, PedigreePairs & pairs, bool is_covariate)
   {
   int filter_limit = (is_covariate ? Person::covariateCount : Person::traitCount);
   if (filter_limit == 0) return;

   IntArray * cts = new IntArray[filter_limit];
   StringArray * corr_txt = new StringArray[filter_limit];

   GetPairTraitStatistics(file, pairs, is_covariate, corr_txt, cts);

   String trait_name;

   fprintf(file, "\nPair Correlations for Each %s:\n", is_covariate ? "Covariate" : "Trait");
   fprintf(file, "%15s %8s %8s %8s %12s %12s %9s\n", " " ,
           "Sib", "HalfSib", "Cousin", "ParentChild", "Grandparent", "Avuncular");

   for(int i = 0; i < filter_limit; i++)
      {
      trait_name = is_covariate ? Person::covariateNames[i] : Person::traitNames[i] ;

      fprintf(file, "%15s %8s %8s %8s %12s %12s %9s\n",(const char *) trait_name,
             (const char *) corr_txt[i][0], (const char *) corr_txt[i][1],
             (const char *) corr_txt[i][2], (const char *) corr_txt[i][3],
             (const char *) corr_txt[i][4], (const char *) corr_txt[i][5]);
      }

   fprintf(file, "\nPair Counts for Each %s:\n", is_covariate ? "Covariate" : "Trait");
   fprintf(file, "%15s %8s %8s %8s %12s %12s %9s\n", " ",
           "Sib", "HalfSib", "Cousin", "ParentChild", "Grandparent", "Avuncular");

   for( int i = 0; i < filter_limit; i++)
      {
      trait_name = is_covariate ? Person::covariateNames[i] : Person::traitNames[i] ;

      fprintf(file, "%15s %8d %8d %8d %12d %12d %9d\n", (const char *) trait_name,
              cts[i][0], cts[i][1], cts[i][2], cts[i][3], cts[i][4], cts[i][5]);
      }

   if (pairs.sibs.gender == "Opposite" && (pairs.avuncularCt + pairs.grandParentCt + pairs.parentCt > 0))
      {
      PedigreePairs mf_pairs(ped), fm_pairs(ped);
      IntArray * mf_cts = new IntArray[filter_limit];
      IntArray * fm_cts = new IntArray[filter_limit];
      StringArray * mf_corr = new StringArray[filter_limit];
      StringArray * fm_corr = new StringArray[filter_limit];

      pairs.SplitOnGenderOrder(mf_pairs, fm_pairs);

      GetPairTraitStatistics(file, fm_pairs, is_covariate, fm_corr, fm_cts, false);
      GetPairTraitStatistics(file, mf_pairs, is_covariate, mf_corr, mf_cts, false);

      bool need_header = true;
      for( int i = 0; i < filter_limit; i++)
         {
         if (fm_cts[i][3] + fm_cts[i][4] + fm_cts[i][5]
             + mf_cts[i][3] + mf_cts[i][4] + mf_cts[i][5] <= 0) continue;

         if (need_header) PrintDetailHeader(file, "- Correlations");
         need_header = false;

         String trait_name = is_covariate ? Person::covariateNames[i] : Person::traitNames[i];
         fprintf(file, "%15s %8s %8s %8s %8s %8s %8s\n", (const char *) trait_name,
                 (const char *) mf_corr[i][3], (const char *) fm_corr[i][3],
                 (const char *) mf_corr[i][4], (const char *) fm_corr[i][4],
                 (const char *) mf_corr[i][5], (const char *) fm_corr[i][5]);
         }

      need_header = true;
      for( int i = 0; i < filter_limit; i++)
         {
         if (fm_cts[i][3] + fm_cts[i][4] + fm_cts[i][5]
             + mf_cts[i][3] + mf_cts[i][4] + mf_cts[i][5] <= 0) continue;

         if (need_header) PrintDetailHeader(file, "- Counts");
         need_header = false;

         String trait_name = is_covariate ? Person::covariateNames[i] : Person::traitNames[i];
         fprintf(file, "%15s %8d %8d %8d %8d %8d %8d\n", (const char *) trait_name,
                 mf_cts[i][3], fm_cts[i][3], mf_cts[i][4], fm_cts[i][4], mf_cts[i][5], fm_cts[i][5]);
         }

      if (mf_cts)
         delete [] mf_cts;
      if (fm_cts)
         delete [] fm_cts;
      if (mf_corr)
         delete [] mf_corr;
      if (fm_corr)
         delete [] fm_corr;
      }
   delete [] cts;
   delete [] corr_txt;
   }

void GetPairTraitStatistics(FILE * file, PedigreePairs & pairs, bool is_covariate,
 StringArray * corr_txt, IntArray * cts, bool do_unordered)
   {
   int filter_limit = (is_covariate ? Person::covariateCount : Person::traitCount);

   for (int i = 0; i < filter_limit; i++)
      {
      cts[i].Dimension(6);
      cts[i].Zero();
      corr_txt[i].Dimension(6);
      }

   for (int filter = 0; filter < filter_limit; filter++)
      {
      Vector correlations(6);
      pairs.GetCorrelations(filter, is_covariate, correlations, cts[filter], do_unordered);

      for (int j = 0; j < 6; j++)
         if (cts[filter][j] > 1)
            corr_txt[filter][j].printf("%7.4lf", correlations[j]);
         else
            corr_txt[filter][j].printf("%7s", "-");
      }
   }

void ByFamilyStatistics()
   {
   FILE * covar_file = NULL;
   FILE * trait_file = NULL;
   FILE * affect_file = NULL;
   FILE * marker_file = NULL;
   FILE * pair_file = NULL;

   if (Person::covariateCount > 0)
      {
      covar_file = fopen("pedstats.cov.fam", "wt");
      if (covar_file == NULL)
         error("Error opening file [%s]\n", "pedstats.cov.fam");

      fprintf(covar_file, "\nCOVARIATE STATISTICS BY FAMILY FOR PEDIGREE FILE: %s\n",
                 (const char * ) pedName);
      }

   if (Person::traitCount > 0)
      {
      trait_file = fopen("pedstats.trait.fam", "wt");
      if (trait_file == NULL)
         error("Error opening file [%s]\n", "pedstats.trait.fam");

      fprintf(trait_file, "\nTRAIT STATISTICS BY FAMILY FOR PEDIGREE FILE: %s\n",
                 (const char * ) pedName);
      }

   if (Person::affectionCount > 0)
      {
      affect_file = fopen("pedstats.aff.fam", "wt");
      if (affect_file == NULL)
         error("Error opening file [%s]\n", "pedstats.aff.fam");

      fprintf(affect_file, "\nAFFECTION STATISTICS BY FAMILY FOR PEDIGREE FILE: %s\n",
             (const char * ) pedName);
      }

   if (Person::markerCount > 0)
      {
      marker_file = fopen("pedstats.mark.fam", "wt");
      if (marker_file == NULL)
         error("Error opening file [%s]\n", "pedstats.mark.fam");

      fprintf(marker_file, "\nMARKER STATISTICS BY FAMILY FOR PEDIGREE FILE: %s\n",
             (const char * ) pedName);
      }

   if (pairs)
      {
      pair_file = fopen("pedstats.pair.fam", "wt");
      if (pair_file == NULL)
         error("Error opening file [%s]\n", "pedstats.pair.fam");

      fprintf(pair_file, "\nPAIR STATISTICS BY FAMILY FOR PEDIGREE FILE: %s\n",
             (const char * ) pedName);
      }

   for (int f = 0; f < ped.familyCount; f++)
      {
      FilteredPedigree curr_ped;
      curr_ped.pd = ped.pd;

      // this subroutine also builds up the filter
      ped.ExtractFamily(f, curr_ped);

      PedigreePairs curr_pairs(curr_ped);
      curr_pairs.Build();

      String added_label;
      added_label.printf(" FOR FAMILY %s", (const char *) ped.families[f]->famid);

      TraitStatistics(trait_file, curr_ped, curr_pairs.sibs, (const char *) added_label);
      CovarStatistics(covar_file, curr_ped, curr_pairs.sibs, (const char *) added_label);
      AffectionStatistics(affect_file, curr_ped, (const char *) added_label);
      MarkerStatistics(marker_file, curr_ped, (const char *) added_label);

      if (pairs)
         PairStatistics(pair_file, curr_ped, curr_pairs, (const char *) added_label);
      }

   if (trait_file)  fclose(trait_file);
   if (covar_file)  fclose(covar_file);
   if (affect_file) fclose(affect_file);
   if (marker_file) fclose(marker_file);
   if (pair_file)   fclose(pair_file);
   }

void MarkerSummaries(PDF * pdf)
   {
   if (ped.markerCount == 0) return;

   if (markerTables)
      {
      if (ped.markerCount > markersReported)
         printf("Switching to summary output mode because there are more than %d markers.\n"
             "See file pedstats.markerinfo for detailed marker information%s\n\n",
              markersReported, hardy || hardyAll || hardyFounders || hardyUnrelated ?
              " and Hardy-Weinberg\ntest results." : ".");

      QualityStatistics();
      fflush(stdout);
      hwe.OpenFiles((const char *) markerOutputName);
      }

   MarkerStatistics(hwe.outputFile, ped);

   if (writeMarkerPDF)
      {
      printf("\nWriting graphical output...\n");
      GraphMarkerStatistics(*pdf);
      }

   bool default_hwe_output  = !hwe.showAll || ped.markerCount <= markersReported || !writePDF || writeMarkerPDF;

   HWStatistics(default_hwe_output ? pdf : NULL);
   }

void QualityStatistics()
   {
   qual.Summarize(HET_WINDOW_SIZE, byGender ? SEX_UNKNOWN : -99);

   if (byGender)
      {
      qual.Summarize(HET_WINDOW_SIZE, SEX_FEMALE);
      qual.Summarize(HET_WINDOW_SIZE, SEX_MALE);
      }
   }

void HWStatistics(PDF * pdf)
   {
   if (hardyAll)
      hwe.TestAllMarkers(stAll, pdf);
   if (hardyFounders)
      hwe.TestAllMarkers(stFounders, pdf);
   if (hardyUnrelated)
      hwe.TestAllMarkers(stSelected, pdf);
   }

void CheckAges(SimplePairList & sib_pairs)
   {
   PedstatsAgeCheck age_check(ped, &sib_pairs);
   age_check.RunChecks();
   }

void FileReport(bool wrote_pdf)
   {
   if (wrote_pdf)
      {
      printf("Graphical output written to %s\n", (const char *) pdfName);
      if (writeMarkerPDF == false && ped.markerCount > markersReported)
         printf("\n\tNOTE: Graphical marker summaries have not been written\n"
               "\tfor this data set because it contains more than %d markers.\n"
               "\tIf you'd like to Pedstats to generate summaries for all\n"
               "\t%d markers, use the --markerPDF option.\n",
                markersReported, ped.markerCount);

      bool default_hwe_output  = !hwe.showAll || ped.markerCount <= markersReported || !writePDF || writeMarkerPDF;
      if (default_hwe_output == false)
         printf("\n\tNOTE: Graphical summaries for Hardy-Weinberg tests have not been written\n"
                "\tfor this data set because more than %d summaries have been requested. \n"
                "\tIf you'd like Pedstats to graphically summarize all %d markers, use the\n"
                "\t--markerPDF option in addition to the --showAll option.\n",
                markersReported, ped.markerCount);
      }

   if (markerTables)
      {
      String hwe_description = "";
      if (hardyFounders || hardyAll || hardyUnrelated)
         {
         String output;
         output.printf("with p-value less than %4.3f ", hwe.significanceCutoff);
         if (hwe.showAll)
            output = "of all markers ";
         hwe_description.printf(" and results for hwe tests \n%s", (const char *) output, hwe.significanceCutoff);
         }

      printf("\nDetailed marker summaries for all %d markers %swritten to file %s\n", ped.markerCount, (const char *) hwe_description,
         (const char *) markerOutputName);
      }

   if (hardyUnrelated)
      printf("\nIndividuals selected for unrelated sample Hardy-Weinberg test written to\nfile %s\n",
             "pedstats.hweselection");

   if (ped.filterOn)
      printf("\nFiltered individuals written to file %s\n", "pedstats.filterlog");

   if (byFamily)
      {
      if (Person::traitCount > 0)
         printf("\nTrait statistics by family written to file %s\n", "pedstats.trait.fam");
      if (Person::covariateCount > 0)
         printf("\nCovariate statistics by family written to file %s\n", "pedstats.cov.fam");
      if (Person::affectionCount > 0)
         printf("\nAffection statistics by family written to file %s\n", "pedstats.aff.fam");
      if (Person::markerCount > 0)
         printf("\nMarker statistics by family written to file %s\n", "pedstats.mark.fam");
      if (pairs && (Person::traitCount > 0 || Person::covariateCount > 0 || Person::markerCount > 0))
         printf("\nPair statistics by family written to file %s\n", "pedstats.pair.fam");
      }

   printf("\nIf you find this program useful in your work, please cite:\n\n"
          "\tWigginton, JE and Abecasis, GR (2005) PEDSTATS: descriptive\n"
          "\tstatistics, graphics and quality assessment for gene mapping data.\n"
          "\tBioinformatics. 21(16):3445-3447\n\n");

   printf("\n\n\n\n");
   }

void ZeroFounderParents(Pedigree & ped)
   {
   for (int f = 0; f < ped.familyCount; f++)
      {
      Family * fam = ped.families[f];

      for (int i = fam->first; i <= fam->last; i++)
         {
         ped[i].AssessStatus();
         if (ped[i].isFounder())
            {
            ped[i].fatid = "0";
            ped[i].motid = "0";
            }
         }
      }
   }

 
