////////////////////////////////////////////////////////////////////// 
// pedstats/SimplePairs.cpp 
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

#include "SimplePairs.h"
#include "FilteredPedigree.h"
#include "PedigreePairs.h"
#include "MathConstant.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

SimplePairList::SimplePairList(FilteredPedigree & p, const char *x_tag, const char *y_tag): ped(p)
   {
   Initialize(x_tag, y_tag);
   }

SimplePairList::SimplePairList(FilteredPedigree & p): ped(p)
   {
   Initialize("", "");
   }

void SimplePairList::Initialize(const char * x_tag, const char * y_tag)
   {
   count = 0;
   size = 1024;
   list = new SimplePair[size];
   xTag = x_tag;
   yTag = y_tag;
   gender = "None";
   }

SimplePairList::~SimplePairList()
   {
   if (list)
      delete [] list;
   list = NULL;
   }

void SimplePairList::Append(SimplePair & pair)
   {
   if (count == size)
      Grow();
   list[count++] = pair;
   }

void SimplePairList::Append(int i, int j)
   {
   SimplePair pair(i, j);
   Append(pair);
   }

void SimplePairList::Grow()
   {
   int newSize = size * 2;
   SimplePair * newList = new SimplePair [newSize];

   memcpy(newList, list, size * sizeof(SimplePair));

   delete [] list;

   list = newList;
   size = newSize;
   }

void SimplePairList::Clear()
   {
   count = 0;
   }

void SimplePairList::BuildSibList()
   {
   for (int f = 0; f < ped.familyCount; f++)
      for (int i = ped.families[f]->first; i <= ped.families[f]->last; i++)
         {
         if (ped.filtered[i]) continue;
         for (int j = i+1; j <= ped.families[f]->last; j++)
            {
            if (ped.filtered[j]) continue;

            if (isSib(&ped[i], &ped[j]))
               Append(i, j);
            }
         }
   }

void SimplePairList::SplitOnAffection(int a, SimplePairList & unaffected,
SimplePairList & discordant, SimplePairList & affected)
   {
   unaffected.Clear();
   discordant.Clear();
   affected.Clear();

   int aff_1, aff_2;

   for (int i = 0; i < count; i++)
      {
      aff_1 = ped[list[i].i1].affections[a];
      aff_2 = ped[list[i].i2].affections[a];

      if (aff_1 == 1 && aff_2 == 1)
         unaffected.Append(list[i]);
      else if ((aff_1 == 1 && aff_2 == 2) || (aff_1 == 2 && aff_2 == 1))
         discordant.Append(list[i]);
      else if (aff_1 == 2 && aff_2 == 2)
         affected.Append(list[i]);
      }

   UpdateTags(unaffected);
   UpdateTags(affected);
   UpdateTags(discordant);
   }

void SimplePairList::CountAffectionTypes(int var, int & unaffected,
int & discordant, int & affected)
   {
   int aff_1, aff_2;

   for (int i = 0; i < count; i++)
      {
      aff_1 = ped[list[i].i1].affections[var];
      aff_2 = ped[list[i].i2].affections[var];

      if (aff_1 == 1 && aff_2 == 1)
         unaffected++;
      else if ((aff_1 == 1 && aff_2 == 2) || (aff_1 == 2 && aff_2 == 1))
         discordant++;
      else if (aff_1 == 2 && aff_2 == 2)
         affected++;
      }
   }

void SimplePairList::SplitOnGender(SimplePairList & females, SimplePairList & males,
SimplePairList & opposites)
   {
   males.Clear();
   females.Clear();
   opposites.Clear();

   int g_1, g_2;

   for (int i = 0; i < count; i++)
      {
      g_1 = ped[list[i].i1].sex;
      g_2 = ped[list[i].i2].sex;

      if (g_1 == SEX_MALE && g_2 == SEX_MALE)
         males.Append(list[i]);
      else if ((g_1 == SEX_FEMALE && g_2 == SEX_MALE))
         opposites.Append(list[i]);
      else if (g_1 == SEX_MALE && g_2 == SEX_FEMALE)
         opposites.Append(list[i]);
      else if (g_1 == SEX_FEMALE && g_2 == SEX_FEMALE)
         females.Append(list[i]);
      }

   UpdateTags(males, "Male");
   UpdateTags(females, "Female");
   UpdateTags(opposites, "Opposite");
   }

void SimplePairList::UpdateTags(SimplePairList & pairs, const char * g_tag)
  {
  pairs.xTag = xTag;
  pairs.yTag = yTag;
  pairs.gender = g_tag != NULL ? g_tag : (const char * ) gender;
  }

// NOTE: This function assumes that first person in pair is ancestor, second is descendant
void SimplePairList::SplitOnGenderOrder(SimplePairList & mf_pairs, SimplePairList & fm_pairs)
   {
   mf_pairs.Clear();
   fm_pairs.Clear();

   for (int i = 0; i < count; i++)
      {
      if (ped[list[i].i1].sex == SEX_MALE && ped[list[i].i2].sex == SEX_FEMALE)
         mf_pairs.Append(list[i]);
      else if (ped[list[i].i1].sex == SEX_FEMALE && ped[list[i].i2].sex == SEX_MALE)
         fm_pairs.Append(list[i]);
      }

   UpdateTags(fm_pairs, "Opposite");
   UpdateTags(mf_pairs, "Opposite");
   }

void SimplePairList::OrderPairsOnGender(SimplePairList & new_list)
  {
  new_list.Clear();
  for (int i = 0; i < count; i++)
     new_list.Append(list[i]);

  for (int i = 0; i < count; i++)
     {
     if (ped[new_list[i].i1].sex == SEX_FEMALE && ped[new_list[i].i2].sex  == SEX_MALE)
        {
        int swap = new_list[i].i1;
        new_list[i].i1 = new_list[i].i2;
        new_list[i].i2 = swap;
        }
     }

  String tag;
  tag  = "Male " + xTag;
  new_list.xTag = tag;

  tag = "Female " + yTag;
  new_list.yTag = tag;
  }

void SimplePairList::FilterOnCovariate(int var, SimplePairList & filtered)
   {
   filtered.Clear();

   for (int i = 0; i < count; i++)
     if (::FilterOnCovariate(var, ped[list[i].i1], ped[list[i].i2]))
        filtered.Append(list[i]);

   UpdateTags(filtered);
   }

void SimplePairList::FilterOnTrait(int var, SimplePairList & filtered)
   {
   filtered.Clear();

   for (int i = 0; i < count; i++)
      if (::FilterOnTrait(var, ped[list[i].i1], ped[list[i].i2]))
        filtered.Append(list[i]);

   UpdateTags(filtered);
   }

void SimplePairList::FilterOnDiagnosed(int var, SimplePairList & filtered)
  {
  filtered.Clear();

  for (int i = 0; i < count; i++)
     if (::FilterOnDiagnosed(var, ped[list[i].i1], ped[list[i].i2]))
         filtered.Append(list[i]);

  UpdateTags(filtered);
  }

void SimplePairList::FilterOnGenotyped(int marker, SimplePairList & filtered)
   {
   filtered.Clear();

   for (int i = 0; i < count; i++)
      if (::FilterOnGenotype(marker, ped[list[i].i1], ped[list[i].i2]))
         filtered.Append(list[i]);

   UpdateTags(filtered);
   }

double SimplePairList::Correlation(int var, bool is_covariate, bool mirrored_pairs,
int & pairs_typed, Vector * x_values, Vector * y_values)
   {
   double s_xx = 0.0, s_yy = 0.0, s_xy = 0.0;
   double x_bar = 0.0, y_bar = 0.0;
   double x_val = 0.0, y_val = 0.0;

   Person * p1, * p2;

   pairs_typed = 0;
   int x_pts = 0, y_pts = 0;

   if (x_values != NULL)
      (*x_values).Dimension(2 * count);

   if (y_values != NULL)
      (*y_values).Dimension(2 * count);

   bool both_phenotyped;
   for (int j = 0; j < count; j++)
      {
      p1 = &ped[list[j].i1];
      p2 = &ped[list[j].i2];

      both_phenotyped = is_covariate ? ::FilterOnCovariate(var, *p1, *p2)
                                     : ::FilterOnTrait(var, *p1, *p2);
      if (!both_phenotyped) continue;

      x_val = is_covariate ? p1->covariates[var] : p1->traits[var];
      y_val = is_covariate ? p2->covariates[var] : p2->traits[var];

      s_xx += (x_val * x_val);
      s_yy += (y_val * y_val);
      s_xy += (x_val * y_val);

      x_bar += x_val;
      y_bar += y_val;

      pairs_typed++;

      if (x_values)
         {
         (*x_values)[x_pts++] = x_val;
         if (mirrored_pairs)
            (*x_values)[x_pts++] = y_val;
         }

      if (y_values)
         {
         (*y_values)[y_pts++]  =  y_val;
         if (mirrored_pairs)
            (*y_values)[y_pts++] = x_val;
         }

      if (mirrored_pairs)
         {
         s_xx += (y_val * y_val);
         s_yy += (x_val * x_val);
         s_xy += (x_val * y_val);

         x_bar += y_val;
         y_bar += x_val;
         pairs_typed++;
         }
      }

   if (x_values != NULL)
      (*x_values).Dimension(x_pts);
   if (y_values != NULL)
      (*y_values).Dimension(y_pts);

   double r_pair = _NAN_;
   if (pairs_typed > 1)
      {
      x_bar /= pairs_typed;
      y_bar /= pairs_typed;
      s_xx -= (pairs_typed * x_bar * x_bar);
      s_yy -= (pairs_typed * y_bar * y_bar);
      s_xy -= (pairs_typed * x_bar * y_bar);

      double denom = s_xx * s_yy;
      if (denom > 0.0) r_pair =  s_xy / sqrt(denom);
      }

   if (mirrored_pairs)
      pairs_typed /= 2;

   return r_pair;
   }

void SimplePairList::FilterOnBothMale(SimplePairList & males)
   {
   males.Clear();

   for (int i = 0; i < count; i++)
      if (::FilterOnBothMale(&ped[list[i].i1], &ped[list[i].i2]))
         males.Append(list[i]);

   UpdateTags(males, "Males");
   }

void SimplePairList::FilterOnBothFemale(SimplePairList & females)
   {
   females.Clear();

   for (int i = 0; i < count; i++)
      if (::FilterOnBothFemale(&ped[list[i].i1], &ped[list[i].i2]))
         females.Append(list[i]);

   UpdateTags(females, "Females");
   }

/////////////////////////////////////////////////////////////////////////////////
//VALUE CHECKING FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////

void SimplePairList::CheckSameTraitValue(int var, const char * relationship_label,
  int digits, bool is_covariate, bool twins_only)
   {
   Person  * p1, * p2;

   double  val_1, val_2;
   double  trait_tol = GetTraitTol(var, is_covariate);

   String  var_label = is_covariate ? (const char *) Person::covariateNames[var]
                                    : (const char *) Person::traitNames[var];

   for (int i = 0; i < count; i++)
      {
      p1 = &ped[list[i].i1];
      p2 = &ped[list[i].i2];

      if (twins_only && !p1->isTwin(*p2)) continue;

      if (!RetrievePairValues(var, p1, p2, val_1, val_2, is_covariate)) continue;

      if ( fabs(val_1 - val_2) > trait_tol);

      printf("\n\nIn family %s, %s values for %s \n  %s (%.*f) and %s (%.*f) are different.\n",
                (const char *) p1->famid, (const char *) var_label, relationship_label,
                (const char *) p1->pid, digits, val_1,
                (const char *) p2->pid, digits, val_2);
      }
   }

void SimplePairList::CheckValueGap(int var, double gap, const char * relationship_label,
 int digits, bool is_covariate, bool min_gap)
   {
   Person  * p1, * p2;

   double  val_1, val_2, diff;

   String var_label = is_covariate ? (const char *) Person::covariateNames[var]
                                   : (const char *) Person::traitNames[var];

   for (int i = 0; i <count; i++)
      {
      p1 = &ped[list[i].i1];
      p2 = &ped[list[i].i2];

      if (!RetrievePairValues(var, p1, p2, val_1, val_2, is_covariate)) continue;

      diff = fabs(val_1 - val_2);

      if ( min_gap && diff >= gap ) continue;
      if (!min_gap && diff <= gap ) continue;

      printf("\n\nIn family %s, %s values for %s \n  %s (%.*f) and %s (%.*f) are %s\n  than might be expected.\n",
                (const char *) p1->famid, (const char *) var_label, relationship_label,
                (const char *) p1->pid, digits, val_1, (const char *) p2->pid, digits, val_2,
                min_gap ? "closer" : "further apart");
      }
   }

double SimplePairList::GetTraitTol(int var, bool is_covariate)
   {
   Person  * p1, * p2;

   double  val_1, val_2;
   double  min_val = 1e300, max_val = -1e300;
   double  pr_min, pr_max;

   bool first_pr = true;

   for (int i = 0; i < count; i++)
      {
      p1 = &ped[list[i].i1];
      p2 = &ped[list[i].i2];

      if (!RetrievePairValues(var, p1, p2, val_1, val_2, is_covariate)) continue;

      if (val_1 < val_2)
         {
         pr_min = val_1;
         pr_max = val_2;
         }
      else
         {
         pr_min = val_2;
         pr_max = val_1;
         }

      min_val = first_pr ? pr_min : min(min_val, pr_min);
      max_val = first_pr ? pr_max : max(max_val, pr_max);

      first_pr = false;
      }

   double trait_tol = !first_pr && (max_val  - min_val) < 1.0 ? TRAIT_TOL * (max_val - min_val) : TRAIT_TOL;

   return trait_tol;
   }

bool SimplePairList::RetrievePairValues(int var, Person * p1, Person * p2,
  double & val_1, double & val_2, bool is_covariate)
   {
   bool both_phenotyped = is_covariate ?  ::FilterOnCovariate(var, *p1, *p2)
                                       :  ::FilterOnTrait(var, *p1, *p2);
   if (!both_phenotyped) return false;

   val_1 = is_covariate ? p1->covariates[var] : p1->traits[var];
   val_2 = is_covariate ? p2->covariates[var] : p2->traits[var];

   return true;
   }


//GLOBAL FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

// returns true if p2 is child of p1
bool isOrderedParent(Person * p1, Person * p2)
   {
   return (!p2->isFounder() && ( p2->father == p1 || p2->mother == p1));
   }

bool isParent(Person * p1, Person * p2)
   {
   Person * ancestor = (p1->traverse < p2->traverse) ? p1 : p2;
   Person * descendant = (p1->traverse < p2->traverse) ? p2 : p1;

   return isOrderedParent(ancestor, descendant);
   }

// returns true if p2  is grandchild of p1
bool isOrderedGrandParentOffspring(Person * p1, Person * p2)
   {
   if (p2->isFounder()) return false;

   if (!p2->father->isFounder() && (p2->father->father == p1 || p2->father->mother == p1))
      return true;

   if (!p2->mother->isFounder() && (p2->mother->father == p1 || p2->mother->mother == p1))
      return true;

   return false;
   }

bool isGrandParentOffspring(Person * p1, Person * p2)
   {
   Person * ancestor = (p1->traverse < p2->traverse) ? p1 : p2;
   Person * descendant = (p1->traverse < p2->traverse) ? p2 : p1;

   return isOrderedGrandParentOffspring(ancestor, descendant);
   }

// returns true if p2 is nephew/niece of p1
bool isOrderedAvuncular(Person * p1, Person * p2)
   {
   if (!p2->isFounder() && !isParent(p1, p2))
      if (p1->isSib(*p2->father) || p1->isSib(*p2->mother))
         return true;

   return false;
   }

bool isAvuncular(Person * p1, Person * p2)
   {
   return (isOrderedAvuncular(p1, p2) || isOrderedAvuncular(p2, p1));
   }

bool isSib(Person * p1, Person * p2)
   {
   if (p1 == p2) return false;

   return p1->isSib(*p2);
   }

bool isHalfSib(Person * p1, Person * p2)
   {
   if (p1 == p2) return false;

   return p1->isHalfSib(*p2);
   }

bool isCousin(Person * p1, Person * p2)
   {
   if (p1->isFounder() || p2->isFounder() || p1->isSib(*p2) || p1->isHalfSib(*p2))
      return false;

   if (!p1->mother->isFounder())
      if (!p2->mother->isFounder() && p1->mother->isSib(*p2->mother) ||
          !p2->father->isFounder() && p1->mother->isSib(*p2->father) )
         return true;

   if (!p1->father->isFounder())
      if (!p2->mother->isFounder() && p1->father->isSib(*p2->mother) ||
          !p2->father->isFounder() && p1->father->isSib(*p2->father) )
         return true;

   return false;
   }

bool FilterOnTrait(int var, Person & p1, Person & p2)
   {
   return ((p1.traits[var] != _NAN_) && (p2.traits[var] != _NAN_));
   }

bool FilterOnCovariate(int var, Person & p1, Person & p2)
   {
   return ((p1.covariates[var] != _NAN_) && (p2.covariates[var] != _NAN_));
   }

bool FilterOnDiagnosed(int affection, Person & p1, Person & p2)
   {
   return (p1.isDiagnosed(affection) && p2.isDiagnosed(affection));
   }

bool FilterOnGenotype(int marker, Person & p1, Person & p2)
   {
   return (p1.isGenotyped(marker) && p2.isGenotyped(marker));
   }

bool FilterOnBothMale(Person * p1, Person * p2)
   {
   return (( p1->sex == SEX_MALE) && (p2->sex == SEX_MALE));
   }

bool FilterOnBothFemale(Person * p1, Person * p2)
   {
   return ((p1->sex == SEX_FEMALE) && (p2->sex == SEX_FEMALE));
   }

bool FilterOnOppositeSex(Person * p1, Person * p2)
   {
   return  ((p1->sex == SEX_MALE && p2->sex == SEX_FEMALE) ||
            (p1->sex == SEX_FEMALE && p2->sex == SEX_MALE) ) ;
   }

void LoadAncestors(Person * p1, IntArray & ancestors)
   {
   ancestors.Push(p1->serial);
   if (p1->isFounder())
      return;

   if (p1->mother)
      LoadAncestors(p1->mother, ancestors);
   if (p1->father)
      LoadAncestors(p1->father, ancestors);
   }

bool isOtherRelative(Person * p1, Person * p2)
   {
   IntArray ancestors;
   bool     related = false;

   if (p1->isFounder() && p2->isFounder())
      return false;

   LoadAncestors(p1, ancestors);
   LoadAncestors(p2, ancestors);

   ancestors.Sort();

   int count = ancestors.Length();
   for (int i = 1; i < count; i++)
      {
      if (ancestors[i-1] == ancestors[i])
         {
         related = true;
         break;
         }
      }
   return related;
   }


 
 
 
 
 
 
 
 
