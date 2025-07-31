////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsAgeCheck.cpp 
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
 
#include "PedstatsAgeCheck.h"
#include "MathConstant.h"
#include <math.h>

String PedstatsAgeCheck::ageLabel = "";
String PedstatsAgeCheck::birthLabel = "";

double PedstatsAgeCheck::minGap = 13.0;
double PedstatsAgeCheck::maxGap = 70.0;
double PedstatsAgeCheck::sibGap = 30.0;

PedstatsAgeCheck::PedstatsAgeCheck(FilteredPedigree & p, SimplePairList * sib_pairs) : ped(p)
   {
   currFamily = NULL;
   offspring = NULL;
   keepSibList = (sib_pairs != NULL);

   sibPairs = sib_pairs != NULL ? sib_pairs : new SimplePairList(p);

   if (sib_pairs == NULL)
      sibPairs->BuildSibList();
   }

PedstatsAgeCheck::~PedstatsAgeCheck()
   {
   if (offspring)
      delete [] offspring;

   if (!keepSibList && sibPairs)
      delete sibPairs;
   }

void PedstatsAgeCheck::RunChecks()
   {
   if (maxGap < minGap)
      error("PedstatsAgeCheck::Minimum generation gap (%lf) exceeds maximum generation gap (%lf)", minGap, maxGap);

   if (birthLabel != "")
      {
      compLabel = birthLabel;
      birthYear = true;
      RunCheck();
      }
   if (ageLabel != "")
      {
      compLabel = ageLabel;
      birthYear = false;
      RunCheck();
      }
   }

void PedstatsAgeCheck::RunCheck()
   {
   bool found_name = false;

   for (int j = 0; j < Person::traitCount; j++)
      {
      if (Person::traitNames[j] != compLabel) continue;
      found_name = true;
      covariateCheck = false;
      varIndex = j;
      }

   for (int j = 0; j < Person::covariateCount; j++)
      {
      if (Person::covariateNames[j] != compLabel) continue;
      found_name = true;
      covariateCheck = true;
      varIndex = j;
      }

   if (!found_name)
      error("An age ordering check was requested for the variable %s, but a variable\n"
            "with this name was not found in this data set. Age labels are case sensitive\n"
            "and must correspond exactly to that listed in your dat file", (const char *) compLabel);

   printf("\n");
   CheckAgeValues(false);
   if (errorCount > 0) 
     printf("\n\n");
   CheckAgeValues(true);
   if (errorCount > 0)
      printf("\n\n");
   CheckSibValues();
   }

void PedstatsAgeCheck::CheckAgeValues(bool max_check)
   {
   errorCount = 0;
   compDir =  birthYear ^ max_check ? -1 : 1;
   compGap =  max_check ? maxGap : minGap;
   compFact = birthYear ? -1 : 1;

   printf("\nChecking for gaps %s %.1f among relative pairs for %s %s ... ",
          max_check ? "greater than" : "less than", compGap,
          covariateCheck ? "covariate" : "trait", (const char *) compLabel);

   for (int f = 0; f < ped.familyCount; f++)
      {
      currFamily = ped.families[f];
      InitializeFamilyCheck();
      RunFamilyCheck(max_check);
      }
   }

void PedstatsAgeCheck::InitializeFamilyCheck()
   {
   if (offspring)
     delete [] offspring;

   ageBoundary.Dimension(currFamily->count);
   ageBoundary.Zero();

   reasonCount.Clear();
   extremeChild.Clear();

   reason.Clear();
   reason.Dimension(currFamily->count);

   reasonCount.Dimension(currFamily->count);
   reasonCount.Zero();

   extremeChild.Dimension(currFamily->count);
   extremeChild.Set(-1);

   offspring = new IntArray[currFamily->count];
   }

void PedstatsAgeCheck::RunFamilyCheck(bool max_check)
   {
   for (int i = currFamily->count - 1; i >= 0; i--)
      {
      Person * person = ped.persons[currFamily->path[i]];
      int father = person->father != NULL ? (person->father)->traverse : -1;
      int mother = person->mother != NULL ? (person->mother)->traverse : -1;

      bool is_known = covariateCheck ? person->isControlled(varIndex) : person->isPhenotyped(varIndex);
      double person_age = covariateCheck ? person->covariates[varIndex] : person->traits[varIndex];
      if (is_known)
         {
         if (fabs(ageBoundary[i]) < 1e-5)
            InitializeLeaf(person);
         else if (compDir * person_age < compDir * ageBoundary[i])
            ReportAgeError(person, max_check);
         }

      if (mother != -1)
         UpdateAge(person, mother);

      if (father != -1)
         UpdateAge(person, father);
      }
   }

void PedstatsAgeCheck::InitializeLeaf(Person * p)
   {
   ageBoundary[p->traverse] = covariateCheck ? p->covariates[varIndex] : p->traits[varIndex];
   reasonCount[p->traverse] ++;
   }

void PedstatsAgeCheck::ReportAgeError(Person * p, bool max_check)
   {
   int person_index = p->traverse;
   PrintAgeError(p, max_check);

   bool first = true;
   int offspring_ct = offspring[person_index].Count();

   for (int c = 0; c < offspring_ct; c++)
      {
      if (offspring[person_index][c] == extremeChild[person_index]) continue;

      IntArray grandchildren;
      CheckChild(p, offspring[person_index][c], first, false, &grandchildren);

      int grandchild_ct = grandchildren.Count();
      for (int g = 0; g < grandchild_ct; g++)
         CheckChild(p, grandchildren[g], first, true);
      }

   reasonCount[p->traverse] ++;
   }

void PedstatsAgeCheck::PrintAgeError(Person * p, bool max_check)
   {
   int person_index = p->traverse;

   printf("\n\n\n\nIn family %s, individual %s ", (const char *) currFamily->famid,
                                              (const char *) p->pid);
   if (birthYear)
      printf("was born in %.0f.\n"
       "  However, %s should have been born %s %.0f since",
       covariateCheck ? p->covariates[varIndex] : p->traits[varIndex], (const char *) p->pid, max_check ? "no sooner than" : "no later than",
       ageBoundary[person_index]);
   else
      printf("has age %.1f.\n"
       "  However, %s should have age %s %.1f, since",
          covariateCheck ? p->covariates[varIndex] : p->traits[varIndex] , (const char *) p->pid, max_check ? "no more than" : "at least",
          ageBoundary[person_index]);

   printf("\n  %s", (const char *) reason[person_index]);
   
   errorCount ++;
   }

void PedstatsAgeCheck::CheckChild(Person * p, int child_index, bool & first_listed, bool is_grandchild, 
 IntArray * grandchildren)
   {
   double   child_age;
   Person * child = ped.persons[currFamily ->path[child_index]];

   if (!RetrieveChildInformation(child, child_age, &offspring[child_index], grandchildren))
       return;

   CheckChildValue(p, child, child_age, first_listed, is_grandchild);
   }

bool PedstatsAgeCheck::RetrieveChildInformation(Person * child, double & child_age, IntArray * child_offspring,
 IntArray * grandchildren)
   {
   bool child_known = covariateCheck ? child->isControlled(varIndex) : child->isPhenotyped(varIndex);

   if (child_known)
       child_age = covariateCheck ? child->covariates[varIndex] : child->traits[varIndex];
   else if (grandchildren)
      grandchildren->Append(*child_offspring);

   return child_known;
   }

void PedstatsAgeCheck::CheckChildValue(Person * p, Person * child, double child_age,
  bool & first_listed, bool is_grandchild)
   {
   double comp_gap = is_grandchild ? compGap * 2.0 : compGap;

   int child_index = child->traverse;
   double person_age = covariateCheck ? p->covariates[varIndex] : p->traits[varIndex];

   if (compDir * person_age < compDir * (child_age + compFact * comp_gap))
      {
      if (first_listed)
         {
         printf("\n\n  Additionally,\n");
         first_listed = false;
         }

      String printed_age = "UKNOWN", printed_age_child = "UNKNOWN";

      printed_age.printf("%.*f", birthYear ? 0 : 1, person_age);
      printed_age_child.printf("%.*f", birthYear ? 0 : 1, child_age);

      printf("\n  %s (%s) is the %s%s of %s (%s)\n",
             (const char * ) p->pid,(const char *) printed_age,
             is_grandchild ? "grand" : "",
             p->sex == SEX_MALE ? "father" : "mother",
             (const char *) child->pid, (const char *) printed_age_child);

      errorCount ++;
      if (reasonCount[child_index] > 1)
         printf("  %s\n", (const char *) reason[child_index]);
      }
   }

void PedstatsAgeCheck::UpdateAge(Person * p, int parent_index)
   {
   int person_index = p->traverse;

   bool person_known = covariateCheck ? p->isControlled(varIndex) : p->isPhenotyped(varIndex);

    
   Person * ancestor = ped.persons[currFamily->path[parent_index]];

   offspring[parent_index].Push(person_index);

   if (!person_known && compGap == maxGap && ageBoundary[person_index] < 1e-5 ) return;
   
   if (compDir * ageBoundary[parent_index] < compDir * (ageBoundary[person_index] + compFact * compGap)
       || fabs(ageBoundary[parent_index]) < 1e-5)
      {
      String printed_age = "UKNOWN", printed_age_parent = "UNKNOWN";

      double person_age = covariateCheck ? p->covariates[varIndex] : p->traits[varIndex];
      double parent_age = covariateCheck ? ancestor->covariates[varIndex] : ancestor->traits[varIndex];

      bool parent_known = covariateCheck ? ancestor->isControlled(varIndex) : ancestor->isPhenotyped(varIndex);

      if (person_known)
         printed_age.printf("%.*f", birthYear ? 0 : 1, person_age);

      if (parent_known)
         printed_age_parent.printf("%.*f", birthYear ? 0 : 1, parent_age);

      reason[parent_index].printf("%s (%s) is the %s of %s (%s).",
                           (const char *) ancestor->pid, (const char *) printed_age_parent,
                           ancestor->sex == SEX_MALE ? "father" : "mother",
                           (const char *) p->pid, (const char *) printed_age);

      if (reasonCount[person_index] > 1 || !person_known)
          {
          reason[parent_index] += "\n  ";
          reason[parent_index] += (const char *) reason[person_index];
          }

      reasonCount[parent_index] = reasonCount[person_index] + 1;
      extremeChild[parent_index] = person_index;
      ageBoundary[parent_index] = ageBoundary[person_index] + compFact * compGap;
      }
   }

void PedstatsAgeCheck::CheckSibValues()
   {
   printf("\nChecking for values that differ among twin pairs for %s %s", 
           covariateCheck ? "covariate" : "trait", (const char *) compLabel);
   sibPairs->CheckSameTraitValue(varIndex, "twins", birthYear ? 0 : 1, covariateCheck, true);

   printf("\nChecking for differences greater than %.*f in sibling values for %s %s", 
           birthYear ? 0 : 1, sibGap, covariateCheck ? "covariate" : "trait", (const char *) compLabel),
   sibPairs->CheckValueGap(varIndex, sibGap, "siblings", birthYear ? 0 : 1, covariateCheck, false);
   }
 
 
 
 
 
 
