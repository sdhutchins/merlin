////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsAgeCheck.h 
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
 
#ifndef __PEDSTATSAGECHECK_H__
#define __PEDSTATSAGECHECK_H__

#include "StringArray.h"
// PedstatsAgeCheck.h
// Written by Jan Wigginton

#include "MathMatrix.h"
#include "FilteredPedigree.h"
#include "SimplePairs.h"
#include "IntArray.h"

class PedstatsAgeCheck
   {
   public:

      static String  ageLabel, birthLabel;
      static double  minGap, maxGap, sibGap;

      void RunChecks();

      PedstatsAgeCheck(FilteredPedigree & p, SimplePairList * pairs = NULL);
      ~PedstatsAgeCheck();

   private:

      FilteredPedigree & ped;
      SimplePairList * sibPairs;

      Family * currFamily;

      bool   birthYear;
      bool   covariateCheck;
      int    varIndex;

      int      compDir, compFact;
      double   compGap;
      String   compLabel;

      bool keepSibList;
      int  errorCount; 

      // Minimum (maximum) allowable values for an individual, given family information
      Vector ageBoundary;

      //Messages used to describe error conditions
      StringArray reason;
   
      //Entry i has lines in error message corresponding to ith individual in pedigreee
      IntArray reasonCount;

      //Indices for offspring of each individual
      IntArray * offspring;

      // Index of (youngest or oldest ?) child for each individual
      IntArray extremeChild;

      void RunCheck();
      void RunFamilyCheck(bool max_check);
      void InitializeFamilyCheck();
      void CheckAgeValues(bool max_check);
      void InitializeLeaf(Person * p);

      void ReportAgeError(Person * p, bool max_check);
      void PrintAgeError(Person * p, bool max_check);
      void CheckChild(Person * p, int child,  bool & first, bool is_grand_child, IntArray * grandchildren = NULL);
      void UpdateAge(Person * p, int parent_index);

      void CheckSibValues();
      void CheckChildValue(Person * p, Person * child, double child_age, bool & first_listed,  bool is_grandchild);
      bool RetrieveChildInformation(Person * child, double & child_known, IntArray * child_offspring,
            IntArray * grandchildren);

   };


#endif


 
 
 
 
 
 
 
