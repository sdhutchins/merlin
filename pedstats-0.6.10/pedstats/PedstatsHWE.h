////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsHWE.h 
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
 
#ifndef __PEDSTATSHWE_H__
#define __PEDSTATSHWE_H__

#include "StringArray.h"
#include "MathMatrix.h"
#include "FilteredPedigree.h"
#include "IntArray.h"
#include "PDF.h"
#include "PDFgrid.h"

enum PedstatsHWESetType    {stAll = 0, stFounders = 1, stSelected = 2};

class PedstatsHWE
   {
   public:

      PedstatsHWE(FilteredPedigree & p);
      ~PedstatsHWE();

      // Cutoff for reporting test results
      static double significanceCutoff;

      // File used for output (default setting is stdout)
      static FILE * outputFile;

      //Indicator to display all test results (default is to display only significant)
      static bool showAll;

      // Index of current marker with respect to Pedigree globals
      int currMarker;

      // Number of alleles at current marker
      int alleleCount;

      // Number of alleles after pooling
      int realAlleles;

      // Number of allele that has pooled alleles added to counts
      int poolAllele;

      // Names of pooled alleles
      String pooledNames;

      // Names for all alleles
      StringArray alleleNames;

      // Allele frequencies after pooling
      Vector alleleFreqs;

      // Indices in allele freq array corresponding to rare alleles
      IntArray pooledAlleles;

      // True when current marker is tested with exact test
      bool useExact;

      // False when current marker was not tested due to error conditions such
      // as low allele counts, only one allele, etc
      bool testPerformed;

      // P-value, chi-square from test of current marker
      double pValue, chisq;

      // Expected heterozygotes for current marker
      int expectedHets;

      // Smallest, largest allele for current marker
      int smallestAllele, largestAllele;

      // Number of failed tests, number of tests done with some alleles pooled for the test
      // Number of failed tests with some alleles pooled
      int failedTests, pooledTests, failedPooledTests;

      // Number of failed tests at 0.01 significance cutoff
      int failedTests_01, failedTests_05;

      // Number of markers tested, number of attempted tests
      int markerTests, attemptedTests;

      // Total genotypes
      int totalGenotypes;

      // Vector of heterozygote probabilities for the exact test
      Vector hetProbs;

      // Observed and expected genotype counts with respect to pooled alleles
      Matrix realExpected, realObserved;

      // Names of markers that are monomorphic or have not been tested due to sparse cells
      StringArray sparseNames, monomorphicNames;

      // Redirect output to a text file
      void OpenFiles(const char * filename);

      bool TestMarker(int m);
      void TestAllMarkers(PedstatsHWESetType t, PDF * pdf = NULL);

      // Indices in pedigree of individuals selected for this test
      IntArray personsSelected;

   private:

      // Sample type for test (all individuals, founders, independent sample)
      PedstatsHWESetType setType;

      //Family-specific cutoff for proportion of tested markers that must be typed in order
      //to considered for inclusion in selected sample -- indexed by family
      Vector typedCutoffs;

      //Proportion of markers genotyped in each individual in the sample
      Vector genotypedProportions;

      // label for sample type
      String textSetTypeLabel, pdfSetTypeLabel;

      // Alias for pedigree
      FilteredPedigree & ped;

      // Raw observed genotype counts
      Matrix observed;

      // True when outputFile is opened and text redirected from stdout
      bool filesOpened;

      double  RunExactTest();
      double  ExactHWE(int rare , int heterozygotes);

      int    DetermineExpectedHeterozygotes();
      double DetermineGenotypeCounts();

      void   DetermineAlleleNames();
      void   DetermineAlleleRange();
      void   DetermineAlleleFreqs(bool final_count);
      void   DeterminePooledAlleleInfo(int m);

      void PoolAlleles(Vector & freqs);
      void InitializeArrays();
      void InitializeStatistics();

      void SelectTestSample();
      bool CheckFounderTyping(int fam);
      void SelectIndependentSet();
      void SelectFounders(int fam);
      void ConstructSampleLabels();

      void SelectIndependentMembers(int fam, const IntArray & markers_tested);
      void LimitToGenotypedPersons(const IntArray & tested_markers);
      void DetermineGenotypedProportions(IntArray & tested_markers);
      void DetermineTypedCutoffs();

      bool IsTested(int m);

      void GraphTestResult(PDF & pdf);
      void GraphAsymptotic(PDF & pdf, PDFGrid & hw_graph, double expected_max);
      void GraphExact(PDF & pdf);

      void MakeSignificanceKey(PDF & pdf, double p_value);
      void DrawGreenChecks(PDF & pdf, double x, double y, double width, double height);
      void DrawRedCross(PDF & pdf, double x, double y, double width, double height);

      void PrintTestResult();
      void PrintTableHeader();
      void PrintTableFooter();
      void PrintSampleSelection();

   };

#endif


 
 
 
 
 
