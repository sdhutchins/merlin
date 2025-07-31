////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsPDF.cpp 
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
 
//Written by Jan Wigginton

#include "Pedstats.h"
#include "PedstatsPDF.h"
#include "QuickIndex.h"

// externally defined functions and variables
//
#define PED_FILTER_CUTOFF 10

extern String pedName;
extern FilteredPedigree ped;
extern PedigreePairs all_pairs;
extern int affection_filter, trait_filter, covariate_filter;
extern bool pairs;
extern bool byGender;

int    pr_graph = 0, pr_gender_graph = 0;

void GraphFamilyStatistics(PDF & pdf, FilteredFamilyStats & stats)
   {
   pdf.page.OpenPage();

   // move lower chart row up slightly if we need to add a footnote to the page
   double offset = ped.filterOn ? 0.04 : 0.0;

   // settings for family member charts
   PDFHistogram family_chart;
   family_chart.useLegend = false;
   family_chart.yAxis.useDiscreteValues = true;

   family_chart.yAxis.label = "Number of pedigrees";
   family_chart.Dimension(1, stats.familiesKept);
   family_chart.xAxis.minMin = -0.5;

   // overall
   family_chart.SetDataValues(0, stats.familyMembers);
   family_chart.SetSeriesColor(0, 1.0, 0.0, 0.0);
   family_chart.SetCategoricalBars(true, stats.minMembers, stats.maxMembers);

   family_chart.title = "Individuals";
   family_chart.xAxis.label = "Number of family members";

   family_chart.DrawInUpperLeft(pdf);

   // founders
   family_chart.SetDataValues(0, stats.familyFounders);
   family_chart.SetCategoricalBars(true, stats.minFounders, stats.maxFounders);
   family_chart.SetSeriesColor(0, 0.0, 1.0, 1.0);

   family_chart.title = "Founders";
   family_chart.xAxis.label = "Number of founders";

   family_chart.DrawInLowerRight(pdf, offset);

   // non-founders
   family_chart.SetDataValues(0, stats.familyNonFounders);
   family_chart.SetCategoricalBars(true, stats.minNonFounders, stats.maxNonFounders);
   family_chart.SetSeriesColor(0, 0.0, 1.0, 0.0);

   family_chart.title = "Non-founders";
   family_chart.xAxis.label = "Number of non-founders";

   family_chart.DrawInLowerLeft(pdf, offset);

   // number of generations
   family_chart.SetDataValues(0, stats.familyGenerations);
   family_chart.SetCategoricalBars(true, stats.minGenerations, stats.maxGenerations);
   family_chart.SetSeriesColor(0, 0.0, 0.0, 1.0);

   family_chart.title = "Generations";
   if (ped.filterOn)
      family_chart.subTitle.printf("( May include filtered individuals )");
   family_chart.xAxis.label = "Number of generations";

   family_chart.DrawInUpperRight(pdf);

   // make page header
   String title, sub_title;

   title.printf("Summary of Family Structure");
   sub_title.printf("( %d individual%s in %d famil%s%s )",
                     stats.personCount,  stats.personCount != 1 ? "s" : "",
                     stats.familiesKept, stats.familiesKept != 1 ? "ies" : "y",
                     ped.filterOn ? "*" : "");

   pdf.page.AddHeader((const char *) title, (const char *) sub_title);
   if (ped.filterOn)
     pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

   pdf.page.ClosePage();
   }

void GraphTraitStatistics(PDF & pdf)
   {
   if (Person::traitCount == 0) return;

   for (trait_filter = 0; trait_filter < Person::traitCount; trait_filter++)
      {
      pr_graph = pr_gender_graph = 0;
      GraphTrait(pdf, trait_filter, false);
      if (pairs || byGender)
         GraphVarForAllPairs(pdf, trait_filter, false);
      }
   }

void GraphTrait(PDF & pdf, int t, bool is_covariate)
   {
   Vector allTraitValues(ped.count);
   Vector founderTraitValues(ped.count), nonFounderTraitValues(ped.count);
   Vector femaleTraitValues(ped.count), maleTraitValues(ped.count);

   allTraitValues.Clear();
   founderTraitValues.Clear();
   nonFounderTraitValues.Clear();
   femaleTraitValues.Clear();
   maleTraitValues.Clear();

   double trait_value;

   int ct_total = 0, ct_founders = 0, ct_nonfounders = 0;
   int ct_females = 0, ct_males = 0;

   bool typed;
   for (int i = 0; i < ped.count; i++)
      {
      if (ped.filtered[i]) continue;

      typed = is_covariate ? ped[i].isControlled(t) : ped[i].isPhenotyped(t);
      if (!typed ) continue;

      ct_total ++;
      trait_value = is_covariate ? ped[i].covariates[t] : ped[i].traits[t];

      allTraitValues.Push(trait_value);

      if (ped[i].isFounder())
         {
         ct_founders ++;
         founderTraitValues.Push(trait_value);
         }
      else
         {
         nonFounderTraitValues.Push(trait_value);
         ct_nonfounders ++;
         }

      if (ped[i].sex == SEX_MALE)
         {
         maleTraitValues.Push(trait_value);
         ct_males ++;
         }
      else if (ped[i].sex == SEX_FEMALE)
         {
         femaleTraitValues.Push(trait_value);
         ct_females ++;
         }
      }

   pdf.page.OpenPage();

   allTraitValues.Dimension(ct_total);
   founderTraitValues.Dimension(ct_founders);
   nonFounderTraitValues.Dimension(ct_nonfounders);
   maleTraitValues.Dimension(ct_males);
   femaleTraitValues.Dimension(ct_females);

   PDFHistogram trait_value_chart;
   trait_value_chart.xAxis.label =  is_covariate ? "Covariate value" : "Trait value";
   trait_value_chart.yAxis.label = "Individuals";
   trait_value_chart.stacked = false;

   GraphVarOverall(pdf, trait_value_chart, allTraitValues);
   GraphVarBySex(pdf, trait_value_chart, femaleTraitValues, maleTraitValues, true);
   GraphVarByFounder(pdf, trait_value_chart, founderTraitValues, nonFounderTraitValues, true);
   GraphVarByPairs(pdf, is_covariate, all_pairs.sibs, true, true);

   String title;
   if (is_covariate)
      title.printf("Distribution of Values for %-15.15s",(const char *)  ped.covariateNames[t]);
   else
      title.printf("Distribution of Values for Trait %s",(const char *)  ped.traitNames[t]);

   String sub_title;
   if (ct_total > 0)
      sub_title.printf( "( Mean: %.2f Variance: %.2f %s )", allTraitValues.Average(),
                       (ct_total <= 1) ? 0.0 : allTraitValues.Var(),
                       ped.filterOn ? "*" : "");
   else if (ped.filterOn)
      title += "*";

   pdf.page.AddHeader((const char *) title, (const char *) sub_title);

   if (ped.filterOn)
      pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

   pdf.page.ClosePage();
   }

void GraphCovariateStatistics(PDF & pdf)
   {
   if (Person::covariateCount == 0) return;

   for (covariate_filter = 0; covariate_filter < Person::covariateCount; covariate_filter++)
      {
      pr_graph = 0;
      pr_gender_graph = 0;
      GraphTrait(pdf, covariate_filter, true);
      if (pairs || byGender)
         GraphVarForAllPairs(pdf, covariate_filter, true);
      }
   }

void GraphVarOverall(PDF & pdf, PDFHistogram & chart, Vector & values, MarkerInfo * marker_info)
   {
   bool use_categories = (marker_info != NULL);

   if (use_categories)
      {
      int n_alleles = values.Length() - 1;
      chart.Dimension(2, n_alleles);

      Matrix info_cts(2, n_alleles);

      int total = int(values.Sum());
      for (int i = 0; i < n_alleles; i++)
         {
         info_cts[0][i] = i + 1;
         info_cts[1][i] = total != 0 ?  (double) values[i + 1] / (double) total : 0.0 ;
         }

      chart.SetCategoricalDataCounts(info_cts);
      }
   else
      {
      chart.Dimension(1, values.dim);
      chart.SetDataValues(0, values);
      }

   chart.SetSeriesColor(0, 0.0, 0.0, 1.0);
   if (!chart.xAxis.useDiscreteValues) chart.SetBarSpacing(0.0);
   chart.narrowBars = true;
   chart.useLegend = false;

   chart.title = "Overall distribution";

   chart.DrawInUpperLeft(pdf);
   chart.narrowBars = false;
   }

void GraphVarBySex(PDF & pdf, PDFHistogram & chart, Vector & v_female, Vector & v_male,
   bool do_subtitle, MarkerInfo * marker_info)
   {
   bool use_categorical = (marker_info != NULL);

   int max_dim = max(v_female.dim, v_male.dim);

   int tot_female = int(v_female.Sum());
   int tot_male   = int(v_male.Sum());

   if (use_categorical)
      {
      int n_alleles = v_female.Length() - 1;
      chart.Dimension(3, n_alleles);

      Matrix data_counts(3, n_alleles);
      for (int i = 0; i < n_alleles; i++)
         {
         data_counts[0][i] = i + 1;
         data_counts[1][i] = (double) v_female[i + 1] / (double) tot_female;
         data_counts[2][i] = (double)   v_male[i + 1] / (double) tot_male;
         }

      chart.SetCategoricalDataCounts(data_counts);
      }
   else
      {
      chart.Dimension(2, max_dim);
      chart.SetDataValues(0, v_female);
      chart.SetDataValues(1, v_male);
      }

   chart.SetSeriesColor(0, 0.0, 0.0, 1.0);
   chart.SetSeriesColor(1, 0.0, 1.0, 0.0);
   if (!chart.xAxis.useDiscreteValues) chart.SetBarSpacing(0.1);
   chart.useLegend = true;

   chart.title = "Distribution by sex ";
   if (do_subtitle)
      {
      int male_ct = v_male.SafeCount();
      int female_ct = v_female.SafeCount();

      chart.subTitle.printf("( %d male%s, %d female%s )",
                              male_ct,   (male_ct != 1) ?  "s" : "",
                              female_ct, (female_ct != 1) ? "s": "" );
      }

   chart.SetSeriesLabel(0, "Females");
   chart.SetSeriesLabel(1, "Males");

   chart.DrawInUpperRight(pdf);
   chart.subTitle = "";
   }

void GraphVarByFounder(PDF & pdf, PDFHistogram & chart, Vector & v_founder, Vector & v_nonfounder,
                        bool do_subtitle, MarkerInfo * marker_info)
   {
   bool use_categorical = (marker_info != NULL);

   if  (use_categorical)
      {
      int n_alleles = v_founder.Length() - 1;
      chart.Dimension(3, n_alleles);

      Matrix data_counts(3, n_alleles);
      int tot_founder = int(v_founder.Sum());
      int tot_nonfounder = int(v_nonfounder.Sum());

      for (int i = 0; i < n_alleles; i++)
         {
         data_counts[0][i] = i + 1;
         data_counts[1][i] = (double) v_founder[i + 1] / (double) tot_founder;
         data_counts[2][i] = (double) v_nonfounder[i + 1] / (double) tot_nonfounder;
         }
      chart.SetCategoricalDataCounts(data_counts);
      }
   else
      {
      int max_dim = max(v_founder.dim, v_nonfounder.dim);
      chart.Dimension(2, max_dim);

      chart.SetDataValues(0, v_founder);
      chart.SetDataValues(1, v_nonfounder);
      }

   if (!chart.xAxis.useDiscreteValues) chart.SetBarSpacing(0.1);

   chart.title = "Distribution by founder status";

   if (do_subtitle)
      {
      int founder_ct = v_founder.SafeCount();
      int nonfounder_ct = v_nonfounder.SafeCount();

      chart.subTitle.printf("( %d founder%s, %d non-founder%s )",
                     founder_ct, founder_ct != 1 ?  "s" : "",
                     nonfounder_ct, nonfounder_ct != 1 ? "s": "");
      }

   chart.SetSeriesLabel(0, "Founder");
   chart.SetSeriesLabel(1, "Non-Founder");

   chart.SetSeriesColor(0, 0.0, 0.0, 1.0);
   chart.SetSeriesColor(1, 0.0, 1.0, 0.0);

   chart.DrawInLowerRight(pdf, ped.filterOn ? 0.04 : 0.0);
   chart.subTitle = "";
   }

void GraphVarByPairs(PDF & pdf, bool is_covariate, SimplePairList & pair_list,
 bool pairs_page, bool is_mirrored, const char * alt_title)
   {
   PDFLineChart chart;

   chart.yAxis.frameData = chart.xAxis.frameData = true;
   chart.symmetricAxes = true;
   chart.title = alt_title;

   if (chart.title != "")
      {
      if (pairs_page && chart.title != "Avuncular plot")
         chart.title += " pairs";
      }
   else
      {
      String lc_tag = pair_list.yTag;
      lc_tag.ToLower();
      chart.title.printf("%s-%s plot", (const char *) pair_list.xTag, (const char *) lc_tag);
      }

   chart.xAxis.label.printf("%s %s value",( const char *) pair_list.xTag,
                           is_covariate ? "covariate" : "trait");

   if (is_mirrored)
      chart.yAxis.label.printf("Co-%s %s value", (const char *) pair_list.yTag,
                               is_covariate ? "covariate" : "trait");
   else
      chart.yAxis.label.printf("%s %s value", (const char *) pair_list.yTag.Capitalize(),
                               is_covariate ? "covariate" : "trait");

   int actual_pairs = 0;
   double r_pair = GetPairData(pair_list, chart, actual_pairs, is_covariate, is_mirrored);

   String temp;
   if (actual_pairs > 0)
      {
      temp.printf("( %d phenotyped pair%s ", actual_pairs, actual_pairs != 1 ? "s" : "");
      String add = " )";
      if (r_pair != _NAN_ && actual_pairs > 1 )
         add.printf(": r = %5.2lf )", r_pair);
      temp += add;
      }
   else
      temp.printf("( No phenotyped pairs )");
   chart.subTitle += temp;

   chart.ShowLine(0, false);
   chart.ShowMarker(0,true);

   if (actual_pairs <= 4)
      chart.SetMarkerRadius(0, 5.0);
   else if (actual_pairs < 20)
      chart.SetMarkerRadius(0, 4.0);
   else if (actual_pairs < 40)
      chart.SetMarkerRadius(0, 2.0);

   chart.useLegend = false;

   pairs_page ? DrawSixth(pdf, &chart, pr_graph) :  DrawQuadrant(pdf, &chart, pr_gender_graph);
   }

double GetPairData(SimplePairList & pair_list, PDFLineChart & chart, int & pairs,
   bool is_covariate, bool is_mirrored)
   {
   Vector x_vals(0);
   Matrix y_vals(1, 0);

   int var = is_covariate ? covariate_filter : trait_filter;
   double r_pair = pair_list.Correlation(var, is_covariate, is_mirrored, pairs, &x_vals, &y_vals[0]);

   int data_points = is_mirrored ? 2 * pairs : pairs;
   x_vals.Dimension(data_points);
   y_vals.Dimension(1, data_points);
   chart.Dimension(1, data_points);

   chart.SetDataValues(x_vals,  y_vals);

   return r_pair;
   }

void GraphMarkerStatistics(PDF & pdf)
   {
   if (Person::markerCount == 0) return;

   for (int m = 0; m < Person::markerCount; m++)
      GraphMarker(pdf, m);
   }

void GraphMarker(PDF & pdf, int m)
   {
   int mGenotypes = 0, mPossible = 0;
   int lastAllele = 0, firstAllele = 100000;
   int ct_males = 0, ct_females = 0;
   int ct_founders = 0, ct_nonfounders = 0;
   int hetCt = 0, homCt = 0;

   MarkerInfo * marker_info = ped.markerInfo[m];

   int n_alleles = marker_info->alleleLabels.Length() - 1;
   bool use_labels = false;
   for (int i = 1; i <= n_alleles; i++)
      use_labels |= (atoi(marker_info->alleleLabels[i])  == 0);

   if (use_labels == false)
      marker_info = NULL;

   int poss_alleles = use_labels ? n_alleles + 1 : 2 * ped.count;

   Vector overallAlleles(poss_alleles);
   Vector founderAlleles(poss_alleles), nonFounderAlleles(poss_alleles);
   Vector maleAlleles(poss_alleles), femaleAlleles(poss_alleles);

   if (use_labels == false)
      {
      overallAlleles.Clear();
      founderAlleles.Clear();
      nonFounderAlleles.Clear();
      femaleAlleles.Clear();
      maleAlleles.Clear();
      }
   else
      {
      overallAlleles.Zero();
      founderAlleles.Zero();
      nonFounderAlleles.Zero();
      femaleAlleles.Zero();
      maleAlleles.Zero();
      }

   int m0, m1;
   for (int i = 0; i < ped.count; i++)
      {
      if (ped.filtered[i]) continue;

      mPossible  ++;

      m0 = ped[i].markers[m].one;
      m1 = ped[i].markers[m].two;

      if (m0 == 0 || m1 == 0) continue;


      if (use_labels)
         {
         overallAlleles[m0] ++;
         overallAlleles[m1] ++;
         }
      else
         {
         overallAlleles.Push(m0);
         overallAlleles.Push(m1);
         }

      mGenotypes ++;


      if (ped[i].isFounder())
         {
         if (use_labels)
            {
            founderAlleles[m0]++;
            founderAlleles[m1]++;
            }
         else
            {
            founderAlleles.Push(m0);
            founderAlleles.Push(m1);
            }

         ct_founders ++;
         }
      else
         {
         if (use_labels)
            {
            nonFounderAlleles[m0] ++;
            nonFounderAlleles[m1] ++;
            }
         else
            {
            nonFounderAlleles.Push(m0);
            nonFounderAlleles.Push(m1);
            }

         ct_nonfounders ++;
         }

      if (ped[i].sex == SEX_MALE)
         {
         if (use_labels)
            {
            maleAlleles[m0]++;
            maleAlleles[m1]++;
            }
         else
            {
            maleAlleles.Push(m0);
            maleAlleles.Push(m1);
            }
         ct_males ++;
         }
      else if (ped[i].sex == SEX_FEMALE)
         {
         if (use_labels)
            {
            femaleAlleles[m0]++;
            femaleAlleles[m1]++;
            }
         else
            {
            femaleAlleles.Push(m0);
            femaleAlleles.Push(m1);
            }

         ct_females ++;
         }

      firstAllele = min(min(m0, m1), firstAllele);
      lastAllele  = max(max(m0, m1), lastAllele);

      if (ped.chromosomeX && ped[i].sex == SEX_MALE) continue;

      m0 != m1 ? hetCt ++ : homCt ++;
      }

   //printf("Number of observed alleles = %d\n", observedAlleles.Length());

   if (use_labels == false)
      {
      overallAlleles.Dimension( 2 * mGenotypes);
      founderAlleles.Dimension( 2 * ct_founders);
      nonFounderAlleles.Dimension( 2 * ct_nonfounders);
      maleAlleles.Dimension( 2 * ct_males);
      femaleAlleles.Dimension( 2 * ct_females);
      }

   pdf.page.OpenPage();

   PDFHistogram marker_chart;

   // Overall distribution
   marker_chart.stacked = false;
   marker_chart.yAxis.maxMax = 1.1;
   //marker_chart.xAxis.SetStringTickLabels(true);

   if (use_labels)
      {
      marker_chart.SetUseSeriesPercentages(false);
      marker_chart.SetCategoricalBars(true, 1.0, n_alleles);
      marker_chart.xAxis.SetStringTickLabels(true);

      for (int i = 0; i < n_alleles; i++)
         marker_chart.xAxis.SetTickLabel(i, (const char *) marker_info->alleleLabels[i + 1]);
      marker_chart.xAxis.useTicks = true;

      marker_chart.xAxis.SetStep(1.0);
      marker_chart.xAxis.SetMin(0.5);
      marker_chart.xAxis.SetMax(n_alleles + 0.5);
      }
   else
      {
      marker_chart.SetUseSeriesPercentages(true);
      marker_chart.SetCategoricalBars(true, firstAllele, lastAllele);
      }

   marker_chart.xAxis.label = use_labels ? "Allele label" : "Allele number";
   marker_chart.yAxis.label = "Allele frequency";

   GraphVarOverall(pdf, marker_chart, overallAlleles, marker_info);

   // By sex distribution
   marker_chart.subTitle.printf("( %d male%s, %d female%s )",
                                 ct_males, (ct_males != 1) ? "s" : "",
                                 ct_females, (ct_females != 1) ? "s" : "");
   GraphVarBySex(pdf, marker_chart, femaleAlleles, maleAlleles, false, marker_info);

   // By founder status distribution
   marker_chart.subTitle.printf("( %d founder%s, %d non-founder%s )",
                                   ct_founders, ct_founders != 1 ? "s" : "",
                                   ct_nonfounders, ct_nonfounders != 1 ? "s" : "");

   GraphVarByFounder(pdf, marker_chart, founderAlleles, nonFounderAlleles, false, marker_info);

   // Informativeness chart
   marker_chart.Dimension(2, max(hetCt, homCt));

   Matrix info_cts(3, 2);

   double div = max((double) hetCt + homCt, 1.0);

   info_cts[0][0] = 1;
   info_cts[0][1] = 2;
   info_cts[1][0] = (double) homCt/ div;
   info_cts[1][1] = 0.0;
   info_cts[2][0] = 0.0;
   info_cts[2][1] = (double) hetCt/ div;

   marker_chart.SetUseSeriesPercentages(false);
   marker_chart.SetCategoricalDataCounts(info_cts);

   marker_chart.SetCategoricalBars(true, 1.0, 2.0);
   marker_chart.stacked = true;

   marker_chart.xAxis.SetStringTickLabels(true);
   marker_chart.xAxis.SetTickLabel(0, "Homozygote");
   marker_chart.xAxis.SetTickLabel(1, "Heterozygote");
   marker_chart.xAxis.useTicks = true;

   marker_chart.xAxis.SetMin(0.5);
   marker_chart.xAxis.SetMax(2.5);

   marker_chart.SetSeriesColor(0, 1.0, 0.0, 0.0);
   marker_chart.SetSeriesColor(1, 0.0, 1.0, 0.0);

   marker_chart.title = "Informativeness";
   marker_chart.useLegend = false;

   if (ped.chromosomeX) marker_chart.subTitle = "( X-linked marker: Among females )";
   marker_chart.yAxis.label = "Proportion of population";
   marker_chart.xAxis.label = "";

   marker_chart.DrawInLowerLeft(pdf, ped.filterOn ? 0.04 : 0.0);
   marker_chart.subTitle = "";

   String title, sub_title;
   title.printf("Marker Allele Frequencies for: %-15.15s", (const char *) Person::markerNames[m]);
   sub_title.printf("( %d of %d %s genotyped%s)", mGenotypes, mPossible,
                     mPossible != 1 ? "individuals" : "individual",
                     ped.filterOn ? "*" : " ");

   pdf.page.AddHeader((const char *) title, (const char *) sub_title);
   if (ped.filterOn)
     pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

   pdf.page.ClosePage();
   }

void GraphAffectionStatistics(PDF & pdf)
   {
   for (affection_filter = 0; affection_filter < Person::affectionCount; affection_filter++)
      {
      pr_graph = pr_gender_graph = 0;

      GraphAffection(pdf, affection_filter);
      if (pairs || byGender)
         GraphAffectionForAllPairs(pdf, affection_filter);
      }
   }

void GraphAffection(PDF & pdf, int a)
   {
   if (Person::affectionCount == 0) return;

   Vector nonFounderAff(ped.count), founderAff(ped.count);
   Vector maleAff(ped.count), femaleAff(ped.count);
   Vector overallAff(ped.count);

   overallAff.Clear();
   nonFounderAff.Clear();
   founderAff.Clear();
   maleAff.Clear();
   femaleAff.Clear();

   int ct_total = 0, ct_possible = 0;
   int ct_founders = 0, ct_nonfounders = 0;
   int ct_males = 0, ct_females = 0;

   for (int i = 0; i < ped.count; i++)
      {
      if (ped.filtered[i]) continue;

      ct_possible ++;
      if (!ped[i].isDiagnosed(a)) continue;
      ct_total ++;

      double aff_stat = ped[i].affections[a];

      overallAff.Push(aff_stat);

      if (ped[i].isFounder())
         {
         ct_founders ++;
         founderAff.Push(aff_stat);
         }
      else
         {
         ct_nonfounders ++;
         nonFounderAff.Push(aff_stat);
         }

      if (ped[i].sex == SEX_MALE)
         {
         ct_males ++;
         maleAff.Push(aff_stat);
         }
      else if (ped[i].sex == SEX_FEMALE)
         {
         ct_females ++;
         femaleAff.Push(aff_stat);
         }
      }

   overallAff.Dimension(ct_total);
   founderAff.Dimension(ct_founders);
   nonFounderAff.Dimension(ct_nonfounders);
   femaleAff.Dimension(ct_females);
   maleAff.Dimension(ct_males);

   pdf.page.OpenPage();

   PDFHistogram affection_chart;

   affection_chart.stacked = false;
   affection_chart.SetCategoricalBars(true, 1.0, 2.0);
   affection_chart.yAxis.useDiscreteValues = true;
   affection_chart.yAxis.SetMaxDigits(0);

   affection_chart.xAxis.SetStringTickLabels(true);

   affection_chart.xAxis.SetTickLabel(2, "");
   affection_chart.xAxis.SetTickLabel(1, "Affected");
   affection_chart.xAxis.SetTickLabel(0, "Unaffected");

   affection_chart.yAxis.label = "Individuals";

   GraphVarOverall(pdf, affection_chart, overallAff);
   GraphVarBySex(pdf, affection_chart, femaleAff, maleAff, true);
   GraphVarByFounder(pdf, affection_chart, founderAff, nonFounderAff, true);
   GraphPairAffection(pdf, affection_chart, all_pairs.sibs, true, a);

   String sub_title, title;

   title.printf("Affection Status Statistics for: %-15.15s",
                  (const char *) Person::affectionNames[a]);

   sub_title.printf("( %d of %d individual%s diagnosed%s)", ct_total,
                     ct_possible, ct_possible != 1 ? "s" : "",
                     ped.filterOn ? "*" : " ");

   pdf.page.AddHeader((const char *) title, (const char *) sub_title);
   if (ped.filterOn)
      pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

   pdf.page.ClosePage();
   }

void GraphPairAffection(PDF & pdf, PDFHistogram & chart, SimplePairList & pair_list,
   bool pairs_page, int index, const char * alt_title)
   {
   int num_pairs = GetPairAffectionData(pair_list, chart, index);

   chart.stacked = true;
   chart.xAxis.SetMin(0.0);
   chart.SetCategoricalBars(true, 1.0, 3.0);
   chart.xAxis.SetStringTickLabels(true);
   chart.yAxis.useDiscreteValues = true;
   chart.yAxis.SetMaxDigits(0);

   if (pr_graph > 0) chart.xAxis.useAlternatingLabels = true;

   chart.useLegend = false;

   chart.xAxis.SetTickLabel(0, "Unaffected");
   chart.xAxis.SetTickLabel(1, "Discordant");
   chart.xAxis.SetTickLabel(2, "Affected");
   chart.xAxis.label = "";

   chart.title = alt_title;
   if (chart.title == "")
      {
      if (pair_list.xTag != "Uncle/Aunt")
         {
         if (pair_list.xTag.AsLower() != pair_list.yTag.AsLower())
            chart.title.printf("%s-%s pairs", (const char * )pair_list.xTag, (const char *) pair_list.yTag);
         else
            chart.title.printf("%s pairs", (const char*) pair_list.xTag);
         }
      else
         chart.title = "Avuncular pairs";
      }

   chart.yAxis.label = chart.title;

   if (num_pairs > 0)
      chart.subTitle.printf("( %d diagnosed %s )", num_pairs, num_pairs != 1 ? "pairs" : "pair");
   else
      chart.subTitle.printf("( No diagnosed pairs )");

   chart.SetSeriesLabel(0, "Unaffected");
   chart.SetSeriesLabel(1, "Discordant");
   chart.SetSeriesLabel(2, "Affected");

   chart.SetSeriesColor(0, 0.0, 1.0, 0.0);
   chart.SetSeriesColor(1, 1.0, 1.0, 0.0);
   chart.SetSeriesColor(2, 1.0, 0.0, 0.0);

   pairs_page ? DrawSixth(pdf, &chart, pr_graph) : DrawQuadrant(pdf, &chart, pr_gender_graph);

   chart.xAxis.useTicks = true;
   chart.xAxis.SetNumericTickLabels(true);
   }

int GetPairAffectionData(SimplePairList & pair_list, PDFHistogram & chart, int index)
   {
   int unaffected = 0, discordant = 0, affected = 0;
   pair_list.CountAffectionTypes(index, unaffected, discordant, affected);

   int max_dim = max(unaffected, max(discordant, affected));
   chart.Dimension(3, max_dim);

   Matrix affection_cts(4, 3);

   affection_cts[0][0] = 1.0;
   affection_cts[0][1] = 2.0;
   affection_cts[0][2] = 3.0;
   affection_cts[1][0] = unaffected;
   affection_cts[1][1] = 0.0;
   affection_cts[1][2] = 0.0;
   affection_cts[2][0] = 0.0;
   affection_cts[2][1] = discordant;
   affection_cts[2][2] = 0.0;
   affection_cts[3][0] = 0.0;
   affection_cts[3][1] = 0.0;
   affection_cts[3][2] = affected;

   chart.SetCategoricalDataCounts(affection_cts);

   return (unaffected + discordant + affected);
   }

void GraphAffectionForAllPairs(PDF & pdf, int filter_var)
   {
   PDFHistogram chart;

   PedigreePairs filtered_pairs(ped);
   all_pairs.FilterOnDiagnosed(filter_var, filtered_pairs);

   if (pairs)
      {
      pdf.page.OpenPage();

      GraphPairAffection(pdf, chart, filtered_pairs.grandparents, true, filter_var);
      GraphPairAffection(pdf, chart, filtered_pairs.parents, true, filter_var);
      GraphPairAffection(pdf, chart, filtered_pairs.avuncular, true, filter_var);

      GraphPairAffection(pdf, chart, filtered_pairs.sibs, true, filter_var);
      GraphPairAffection(pdf, chart, filtered_pairs.halfSibs, true, filter_var);
      GraphPairAffection(pdf, chart, filtered_pairs.cousins, true, filter_var);

      String title;
      title.printf("Relative Pair Distributions for %s%s",
                (const char *) Person::affectionNames[filter_var],
                ped.filterOn ? "*" : "");

      pdf.page.AddHeader((const char *) title, "");
      if (ped.filterOn)
         pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

      pdf.page.ClosePage();
      }

   if (byGender && filtered_pairs.totalCt > PED_FILTER_CUTOFF)
      {
      if (filtered_pairs.sibCt > PED_FILTER_CUTOFF)
         DoAffectionPairsByGender(pdf, chart, filtered_pairs.sibs, true);
      if (filtered_pairs.halfSibCt > PED_FILTER_CUTOFF)
         DoAffectionPairsByGender(pdf, chart, filtered_pairs.halfSibs, true);
      if (filtered_pairs.cousinCt > PED_FILTER_CUTOFF)
         DoAffectionPairsByGender(pdf, chart, filtered_pairs.cousins, true);

      if (filtered_pairs.grandParentCt > PED_FILTER_CUTOFF)
         DoAffectionPairsByGender(pdf, chart, filtered_pairs.grandparents, false);
      if (filtered_pairs.parentCt > PED_FILTER_CUTOFF)
         DoAffectionPairsByGender(pdf, chart, filtered_pairs.parents, false);
      if (filtered_pairs.avuncularCt > PED_FILTER_CUTOFF)
         DoAffectionPairsByGender(pdf, chart, filtered_pairs.avuncular, false);
      }
   }

void GraphVarForAllPairs(PDF & pdf, int filter_var, bool is_covariate)
   {
   PedigreePairs filtered_pairs(ped);
   if (is_covariate)
      all_pairs.FilterOnCovariate(filter_var, filtered_pairs);
   else
      all_pairs.FilterOnTrait(filter_var, filtered_pairs);

   if (pairs)
      {
      pdf.page.OpenPage();

      GraphVarByPairs(pdf, is_covariate, filtered_pairs.grandparents, true, false);
      GraphVarByPairs(pdf, is_covariate, filtered_pairs.parents, true, false);
      GraphVarByPairs(pdf, is_covariate, filtered_pairs.avuncular, true, false, "Avuncular plot");

      GraphVarByPairs(pdf, is_covariate, filtered_pairs.sibs, true, true);
      GraphVarByPairs(pdf, is_covariate, filtered_pairs.halfSibs, true, true);
      GraphVarByPairs(pdf, is_covariate, filtered_pairs.cousins, true, true);

      String title;
      title.printf("Relative Pair Plots for %s%s", is_covariate ?
                (const char *) ped.covariateNames[filter_var] :
                (const char *) ped.traitNames[filter_var],
                ped.filterOn ? "*" : "");

      pdf.page.AddHeader((const char *) title, "");
      if (ped.filterOn)
         pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

      pdf.page.ClosePage();
      }

   if (byGender && filtered_pairs.totalCt > PED_FILTER_CUTOFF)
      {
      if (filtered_pairs.sibCt > PED_FILTER_CUTOFF)
         DoVarPairsByGender(pdf, is_covariate, filtered_pairs.sibs, true);
      if (filtered_pairs.halfSibCt > PED_FILTER_CUTOFF)
         DoVarPairsByGender(pdf, is_covariate, filtered_pairs.halfSibs, true);
      if (filtered_pairs.cousinCt > PED_FILTER_CUTOFF)
         DoVarPairsByGender(pdf, is_covariate, filtered_pairs.cousins, true);

      if (filtered_pairs.grandParentCt > PED_FILTER_CUTOFF)
         DoVarPairsByGender(pdf, is_covariate, filtered_pairs.grandparents, false);
      if (filtered_pairs.parentCt > PED_FILTER_CUTOFF)
         DoVarPairsByGender(pdf, is_covariate, filtered_pairs.parents, false);
      if (filtered_pairs.avuncularCt > PED_FILTER_CUTOFF)
         DoVarPairsByGender(pdf, is_covariate, filtered_pairs.avuncular, false);
      }
   }

void DoVarPairsByGender( PDF & pdf, bool is_covariate, SimplePairList & pair_list,
 bool is_mirrored)
   {
   pdf.page.OpenPage();
   String alt_title;

   pr_gender_graph = 0;

   SimplePairList female_pairs(ped), male_pairs(ped), opp_pairs(ped);
   pair_list.SplitOnGender(female_pairs, male_pairs, opp_pairs);

   GetGenderPlotTitle(male_pairs, alt_title);
   GraphVarByPairs(pdf, is_covariate, male_pairs, false, is_mirrored, (const char *) alt_title);

   GetGenderPlotTitle(female_pairs, alt_title);
   GraphVarByPairs(pdf, is_covariate, female_pairs, false, is_mirrored, (const char *) alt_title);

   if (is_mirrored)
      {
      GetGenderPlotTitle(opp_pairs, alt_title);
      String final_title = "Mirrored " + alt_title;
      GraphVarByPairs(pdf, is_covariate, opp_pairs, false, true, (const char *) final_title);

      SimplePairList ordered_opp_pairs(ped);
      opp_pairs.OrderPairsOnGender(ordered_opp_pairs);
      final_title = "Ordered " + alt_title;
      GraphVarByPairs(pdf, is_covariate, ordered_opp_pairs, false, false, (const char *) final_title);
      }
   else
      {
      SimplePairList fm_pairs(ped), mf_pairs(ped);
      opp_pairs.SplitOnGenderOrder(mf_pairs, fm_pairs);

      GetGenderPlotTitle(mf_pairs, alt_title);
      GraphVarByPairs(pdf, is_covariate, mf_pairs, false,  is_mirrored, (const char *) alt_title);

      GetGenderPlotTitle(fm_pairs, alt_title);
      GraphVarByPairs(pdf, is_covariate, fm_pairs, false, is_mirrored, (const char *) alt_title);
      }

   String title;
   if (pair_list.xTag != "Avuncular" && pair_list.xTag != "Uncle/Aunt")
      title.printf("%s-%s Plots by Gender for %s%s",
                   (const char *) pair_list.xTag, (const char *) pair_list.yTag.Capitalize(),
                   is_covariate ?
                   (const char *) ped.covariateNames[covariate_filter] :
                   (const char *) ped.traitNames[trait_filter],
                   ped.filterOn ? "*" : "");
   else
      title.printf("%s Plots by Gender for %s%s", "Avuncular",
                   is_covariate ?
                   (const char *) ped.covariateNames[covariate_filter] :
                   (const char *) ped.traitNames[trait_filter],
                   ped.filterOn ? "*" : "");

   pdf.page.AddHeader((const char *) title, "");
   if (ped.filterOn)
      pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);

   pdf.page.ClosePage();
   }

void DoAffectionPairsByGender(PDF & pdf, PDFHistogram & chart, SimplePairList & pair_list,
bool is_mirrored)
   {
   pdf.page.OpenPage();

   // page title
   String title;
   if (is_mirrored)
      title.printf("%s Pair Distributions by Gender for %s", (const char * ) pair_list.xTag,
               (const char *) Person::affectionNames[affection_filter] );
   else if (pair_list.xTag != "Uncle/Aunt")
      title.printf("%s-%s Pair Distributions by Gender for %s", (const char * ) pair_list.xTag,
               (const char *) pair_list.yTag.Capitalize(),
               (const char *) Person::affectionNames[affection_filter]);
   else
      title.printf("%s Pair Distributions by Gender for %s", "Avuncular",
               (const char *) Person::affectionNames[affection_filter]);

   if (ped.filterOn)
      {
      title += "*";
      pdf.page.AddFootNote(0.07, (const char *) ped.pdfFilterLabel);
      }
   pdf.page.AddHeader((const char *) title, "");


   // page body
   String alt_title;
   pr_gender_graph = 0;

   SimplePairList female_pairs(ped), male_pairs(ped), opp_pairs(ped);
   pair_list.SplitOnGender(female_pairs, male_pairs, opp_pairs);

   GetGenderPlotTitle(male_pairs, alt_title);
   GraphPairAffection(pdf, chart, male_pairs, false, affection_filter, alt_title);

   GetGenderPlotTitle(female_pairs, alt_title);
   GraphPairAffection(pdf, chart, female_pairs, false, affection_filter, alt_title);

   if (is_mirrored)
      {
      GetGenderPlotTitle(opp_pairs, alt_title);
      GraphPairAffection(pdf, chart, opp_pairs, false, affection_filter, alt_title);
      }
   else
      {
      SimplePairList mf_pairs(ped), fm_pairs(ped);
      opp_pairs.SplitOnGenderOrder(mf_pairs, fm_pairs);

      GetGenderPlotTitle(mf_pairs, alt_title);
      GraphPairAffection(pdf, chart, mf_pairs, false, affection_filter, alt_title);

      GetGenderPlotTitle(fm_pairs, alt_title);
      GraphPairAffection(pdf, chart, fm_pairs, false, affection_filter, alt_title);
      }

   pdf.page.ClosePage();
   }

void GetGenderPlotTitle(SimplePairList & pairs, String & alt_title)
   {
   String xtag = pairs.xTag;
   String ytag;

   if (pairs.xTag == "Cousin")
      {
      if (pairs.gender == "Opposite")
         alt_title =  "Male-female";
      else
         alt_title.printf("%s-%s", (const char * ) pairs.gender, (const char *) pairs.gender);

      return;
      }

   bool avuncular = (pairs.xTag == "Uncle/Aunt" || pairs.xTag == "Avuncular");
   bool parents = (pairs.xTag == "Parent");
   bool sibs = (pairs.xTag == "Sib");
   bool half_sibs = (pairs.xTag == "Half sib");

   if (sibs || half_sibs)
      {
      if (pairs.gender != "Female")
         xtag = sibs ? "Brother" : "Half brother";
      else
         xtag = sibs ? "Sister" : "Half sister";

      if (pairs.gender != "Male")
         ytag = sibs ? "Sister" : "Half sister";
      else
         ytag = xtag;
      }
   else
      {
      if (pairs.gender == "Male" || pairs.gender == "Opposite" && !(pr_gender_graph%2))
         xtag = avuncular ? "Uncle" : (parents ? "Father" : "Grandfather");
      else
         xtag = avuncular ? "Aunt" : (parents ? "Mother" : "Grandmother");

      if (pairs.gender == "Female" || pairs.gender == "Opposite" && !(pr_gender_graph%2))
         ytag = avuncular ? "Niece" : (parents ? "Daughter" : "Granddaughter" );
      else
         ytag = avuncular ? "Nephew" : (parents ? "Son" : "Grandson");

      pairs.xTag = xtag;
      pairs.yTag = ytag;
      }

   alt_title.printf("%s-%s pairs", (const char *) xtag, (const char * ) ytag);
   }

void DrawSixth(PDF & pdf, PDFChartBasics * chart, int & pr_graph)
   {
   switch (pr_graph)
      {
      case 0: chart->DrawInLowerLeft(pdf, ped.filterOn ? 0.04 : 0.0); break;
      case 1: chart->DrawInGrid(pdf, 0, 0, 2, 3); break;
      case 2: chart->DrawInGrid(pdf, 0, 1, 2, 3); break;
      case 3: chart->DrawInGrid(pdf, 0, 2, 2, 3); break;
      case 4: chart->DrawInGrid(pdf, 1, 0, 2, 3); break;
      case 5: chart->DrawInGrid(pdf, 1, 1, 2, 3); break;
      case 6: chart->DrawInGrid(pdf, 1, 2, 2, 3); break;
      }
   pr_graph ++;
   }

void DrawQuadrant(PDF & pdf, PDFChartBasics * chart, int & pr_gender_graph)
   {
   switch (pr_gender_graph%4)
      {
      case 0: chart->DrawInUpperLeft(pdf); break;
      case 1: chart->DrawInUpperRight(pdf); break;
      case 2: chart->DrawInLowerLeft(pdf, ped.filterOn ? 0.04 : 0.0); break;
      case 3: chart->DrawInLowerRight(pdf, ped.filterOn ? 0.04 : 0.0); break;
      }
   pr_gender_graph ++;
   }


 
 

 
 
