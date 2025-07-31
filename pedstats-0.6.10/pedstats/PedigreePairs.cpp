////////////////////////////////////////////////////////////////////// 
// pedstats/PedigreePairs.cpp 
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

#include "PedigreePairs.h"
#include "Error.h"

PedigreePairs::PedigreePairs(FilteredPedigree  & p): ped(p), sibs(p, "Sib", "sib"),
 halfSibs(p, "Half sib", "half sib"), cousins(p, "Cousin", "cousin"),
 avuncular(p, "Uncle/Aunt", "nephew/niece"),
 grandparents(p, "Grandparent", "grandchild"),
 parents(p, "Parent", "child"),
 others(p, "Other", "other")
   {
   doOthers = false;
   halfSibCt = sibCt = cousinCt = 0;
   avuncularCt = grandParentCt = parentCt = 0;
   otherCt = totalCt = 0;
   }

// This builds the set of 6 pair lists -> note for pair types with intrinsic ordering,
// the ancestor index is listed as the first element.
void PedigreePairs::Build()
   {
   for (int f = 0; f < ped.familyCount; f++)
      for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         {
         if (ped.filtered[i]) continue;
         for (int j = i+1; j <= ped.families[f]->last; j++)
            {
            if (ped.filtered[j]) continue;

            if (isSib(&ped[i], &ped[j]))
               sibs.Append(i, j);
            else if (isOrderedParent(&ped[i], &ped[j]))
               parents.Append(i, j);   // i is parent of j
            else if (isOrderedParent(&ped[j], &ped[i]))
               parents.Append(j, i);  // j is parent of i
            else if (isHalfSib(&ped[i], &ped[j]))
               halfSibs.Append(i, j);
            else if (isOrderedAvuncular(&ped[i], &ped[j]))
               avuncular.Append(i, j);
            else if (isOrderedAvuncular(&ped[j], &ped[i]))
               avuncular.Append(j, i);
            else if (isCousin(&ped[i], &ped[j]))
               cousins.Append(i, j);
            else if (isOrderedGrandParentOffspring(&ped[i], &ped[j]))
               grandparents.Append(i, j);
            else if (isOrderedGrandParentOffspring(&ped[j], &ped[i]))
               grandparents.Append(j, i);
            else if (isOtherRelative(&ped[j], &ped[i]))
               others.Append(i, j);
            }
         }
   UpdateCount();
   }

void PedigreePairs::UpdateCount()
   {
   sibCt = sibs.count;
   halfSibCt = halfSibs.count;
   cousinCt = cousins.count;
   parentCt = parents.count;
   grandParentCt = grandparents.count;
   avuncularCt = avuncular.count;

   otherCt = others.count;

   totalCt = sibCt + halfSibCt + cousinCt + parentCt + grandParentCt + avuncularCt;

   if (doOthers)
     totalCt += otherCt;
   }

void PedigreePairs::SplitOnGender(PedigreePairs & female_pairs, PedigreePairs & male_pairs,
PedigreePairs & opp_pairs)
   {
   sibs.SplitOnGender(female_pairs.sibs, male_pairs.sibs, opp_pairs.sibs);
   halfSibs.SplitOnGender(female_pairs.halfSibs, male_pairs.halfSibs, opp_pairs.halfSibs);
   cousins.SplitOnGender(female_pairs.cousins, male_pairs.cousins, opp_pairs.cousins);
   avuncular.SplitOnGender(female_pairs.avuncular, male_pairs.avuncular, opp_pairs.avuncular);
   grandparents.SplitOnGender(female_pairs.grandparents, male_pairs.grandparents, opp_pairs.grandparents);
   parents.SplitOnGender(female_pairs.parents, male_pairs.parents, opp_pairs.parents);

   female_pairs.UpdateCount();
   male_pairs.UpdateCount();
   opp_pairs.UpdateCount();
   }

void PedigreePairs::SplitOnGenderOrder(PedigreePairs & mf_pairs, PedigreePairs & fm_pairs)
   {
   avuncular.SplitOnGenderOrder(mf_pairs.avuncular, fm_pairs.avuncular);
   grandparents.SplitOnGenderOrder(mf_pairs.grandparents, fm_pairs.grandparents);
   parents.SplitOnGenderOrder(mf_pairs.parents, fm_pairs.parents);

   mf_pairs.UpdateCount();
   fm_pairs.UpdateCount();
   }

void PedigreePairs::SplitOnAffection(int affection, PedigreePairs & unaffected,
PedigreePairs & discordant, PedigreePairs & affected)
   {
   sibs.SplitOnAffection(affection, unaffected.sibs, discordant.sibs, affected.sibs);
   halfSibs.SplitOnAffection(affection, unaffected.halfSibs, discordant.halfSibs, affected.halfSibs);
   cousins.SplitOnAffection(affection, unaffected.cousins, discordant.cousins, affected.cousins);
   avuncular.SplitOnAffection(affection, unaffected.avuncular, discordant.avuncular, affected.avuncular);
   grandparents.SplitOnAffection(affection, unaffected.grandparents, discordant.grandparents, affected.grandparents);
   parents.SplitOnAffection(affection, unaffected.parents, discordant.parents, affected.parents);

   unaffected.UpdateCount();
   discordant.UpdateCount();
   affected.UpdateCount();
   }

void PedigreePairs::FilterOnCovariate(int covariate, PedigreePairs & filtered)
   {
   sibs.FilterOnCovariate(covariate, filtered.sibs);
   halfSibs.FilterOnCovariate(covariate, filtered.halfSibs);
   cousins.FilterOnCovariate(covariate, filtered.cousins);
   parents.FilterOnCovariate(covariate, filtered.parents);
   grandparents.FilterOnCovariate(covariate, filtered.grandparents);
   avuncular.FilterOnCovariate(covariate, filtered.avuncular);

   filtered.UpdateCount();
   }

void PedigreePairs::FilterOnTrait(int trait, PedigreePairs & filtered)
   {
   sibs.FilterOnTrait(trait, filtered.sibs);
   halfSibs.FilterOnTrait(trait, filtered.halfSibs);
   cousins.FilterOnTrait(trait, filtered.cousins);
   parents.FilterOnTrait(trait, filtered.parents);
   grandparents.FilterOnTrait(trait, filtered.grandparents);
   avuncular.FilterOnTrait(trait, filtered.avuncular);

   filtered.UpdateCount();
   }

void PedigreePairs::FilterOnDiagnosed(int affection, PedigreePairs & filtered)
   {
   sibs.FilterOnDiagnosed(affection, filtered.sibs);
   halfSibs.FilterOnDiagnosed(affection, filtered.halfSibs);
   cousins.FilterOnDiagnosed(affection, filtered.cousins);
   parents.FilterOnDiagnosed(affection, filtered.parents);
   grandparents.FilterOnDiagnosed(affection, filtered.grandparents);
   avuncular.FilterOnDiagnosed(affection, filtered.avuncular);

   filtered.UpdateCount();
   }

void PedigreePairs::FilterOnGenotyped(int marker, PedigreePairs & filtered)
   {
   sibs.FilterOnGenotyped(marker, filtered.sibs);
   halfSibs.FilterOnGenotyped(marker, filtered.halfSibs);
   cousins.FilterOnGenotyped(marker, filtered.cousins);
   parents.FilterOnGenotyped(marker, filtered.parents);
   grandparents.FilterOnGenotyped(marker, filtered.grandparents);
   avuncular.FilterOnGenotyped(marker, filtered.avuncular);
   others.FilterOnGenotyped(marker, filtered.others);

   filtered.UpdateCount();
   }

void PedigreePairs::PrintCounts(FILE * file, const char * label)
   {
   if (label)
     fprintf(file, "%s\n", (const char *) label);

   if (totalCt == 0)
      fprintf(file, " No relative pairs\n\n");
   else
      {
      fprintf(file, "\nRelative Pair Counts:\n" );

      if (sibCt)
         fprintf(file, "%25s: %8d %s\n", "Sib-pairs", sibCt, sibCt == 1 ? "pair" : "pairs");
      if (halfSibCt)
         fprintf(file, "%25s: %8d %s\n", "Half-Sibs", halfSibCt, halfSibCt == 1 ? "pair" : "pairs");
      if (cousinCt)
         fprintf(file, "%25s: %8d %s\n", "Cousins", cousinCt, cousinCt == 1 ? "pair" : "pairs");

      if (sibs.gender != "Opposite")
         {
         if (parentCt)
            fprintf(file, "%25s: %8d %s\n", "Parent-Child", parentCt,
                    parentCt == 1 ? "pair" : "pairs");
         if (grandParentCt)
            fprintf(file, "%25s: %8d %s\n", "Grandparent-Grandchild",
                    grandParentCt, grandParentCt == 1 ? "pair" : "pairs");
         if (avuncularCt)
            fprintf(file, "%25s: %8d %s\n", "Avuncular", avuncularCt,
                    avuncularCt == 1 ? "pair" : "pairs");
         }
      else
         {
         PedigreePairs fm_pairs(ped), mf_pairs(ped);
         SplitOnGenderOrder(mf_pairs, fm_pairs);

         if (parentCt)
            fprintf(file, "%25s: %8d %s (%d F/D, %d M/S)\n", "Parent-Child",
               parentCt, parentCt == 1 ? "pair" : "pairs",
               mf_pairs.parentCt, fm_pairs.parentCt);

         if (grandParentCt)
            fprintf(file, "%25s: %8d %s (%d GF/GD, %d GM/GS)\n", "Grandparent-Grandchild",
               grandParentCt, grandParentCt == 1 ? "pair" : "pairs",
               mf_pairs.grandParentCt, fm_pairs.grandParentCt);

         if (avuncularCt)
            fprintf(file, "%25s: %8d %s (%d U/N, %d A/N)\n", "Avuncular",
               avuncularCt, avuncularCt == 1 ? "pair" : "pairs",
               mf_pairs.avuncularCt, fm_pairs.avuncularCt);
         }
      }
   }

void PedigreePairs::GetCorrelations(int var, int is_covariate, Vector & corr,
IntArray & counts, bool do_unordered)
   {
   corr.Clear();
   corr.Dimension(6);

   counts.Dimension(6);
   counts.Zero();

   if (do_unordered)
      {
      corr[0] = sibs.Correlation(var, is_covariate, true, counts[0]);
      corr[1] = halfSibs.Correlation(var, is_covariate, true, counts[1]);
      corr[2] = cousins.Correlation(var, is_covariate, true, counts[2]);
      }

   corr[3] = parents.Correlation(var, is_covariate, false, counts[3]);
   corr[4] = grandparents.Correlation(var, is_covariate, false, counts[4]);
   corr[5] = avuncular.Correlation(var, is_covariate, false, counts[5]);
   }

void PedigreePairs::CountAffectionTypes(int var, IntArray & unaffected,
IntArray & discordant, IntArray & affected, bool count_unordered)
   {
   unaffected.Dimension(6);
   discordant.Dimension(6);
   affected.Dimension(6);

   unaffected.Zero();
   discordant.Zero();
   affected.Zero();

   if (count_unordered)
     {
     sibs.CountAffectionTypes(var, unaffected[0], discordant[0], affected[0]);
     halfSibs.CountAffectionTypes(var, unaffected[1], discordant[1], affected[1]);
     cousins.CountAffectionTypes(var, unaffected[2], discordant[2], affected[2]);
     }

   parents.CountAffectionTypes(var, unaffected[3], discordant[3], affected[3]);
   grandparents.CountAffectionTypes(var, unaffected[4], discordant[4], affected[4]);
   avuncular.CountAffectionTypes(var, unaffected[5], discordant[5], affected[5]);
   }



 
 
 
 
 
 
 
 
