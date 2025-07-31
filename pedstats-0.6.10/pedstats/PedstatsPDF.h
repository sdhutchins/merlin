////////////////////////////////////////////////////////////////////// 
// pedstats/PedstatsPDF.h 
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
 
#ifndef __PEDSTATS_PDF_H__
#define __PEDSTATS_PDF_H__

#include "PDFhistogram.h"
#include "PDFlinechart.h"
#include "PDFgrid.h"
#include "PedigreePairs.h"
#include "PedstatsHWE.h"
#include "FilteredFamilyStats.h"

// Code for main graph pages
//
void GraphAffectionStatistics(PDF &pdf);
void GraphAffection(PDF & pdf, int a);

void GraphFamilyStatistics(PDF & pdf, FilteredFamilyStats & stats);

void GraphMarkerStatistics(PDF & pdf);
void GraphMarker(PDF & pdf, int m);

void GraphCovariateStatistics(PDF & pdf);
void GraphTraitStatistics(PDF & pdf);
void GraphTrait(PDF & pdf, int t, bool is_covariate);

// Code for single graphs
//
void GraphVarOverall(PDF & pdf, PDFHistogram & chart, Vector & values, MarkerInfo * marker_info = NULL);
void GraphVarBySex(PDF & pdf, PDFHistogram & chart, Vector & v_female,
     Vector & v_male, bool do_subtitle, MarkerInfo * marker_info = NULL);
void GraphVarByFounder(PDF & pdf, PDFHistogram & chart, Vector & v_founder,
     Vector & v_nonfounder, bool do_subtitle, MarkerInfo * marker_info = NULL);

// Code for single (pair) graphs
//
void   GraphPairAffection(PDF & pdf, PDFHistogram & chart, SimplePairList & pair_list,
       bool pairs_page, int index, const char * alt_title = "");
int    GetPairAffectionData(SimplePairList & pair_list, PDFHistogram & chart, int index);

void   GraphVarByPairs(PDF & pdf, bool is_covariate, SimplePairList & pairs, bool prs_page,
       bool is_mirrored, const char * alt_title = "");
double GetPairData(SimplePairList & pair_list, PDFLineChart & chart,
       int & sib_pairs, bool is_covariate = false, bool is_mirrored = true);

// Code for all pair pages (--pairs)
//
void GraphVarForAllPairs(PDF & pdf, int filter_var, bool is_covariate);
void GraphAffectionForAllPairs(PDF & pdf, int filter_var);

// Code for gender-specific pair pages (--bySex)
//
void DoAffectionPairsByGender(PDF & pdf, PDFHistogram & chart, SimplePairList & pair_list,
     bool is_mirrored);
void DoVarPairsByGender(PDF & pdf, bool is_covariate, SimplePairList & pair_list,
     bool is_mirrored);
void GetGenderPlotTitle(SimplePairList & pairs, String & alt_title);

// Code for data quality graphing

void GraphQualityHistogram(PDF & pdf, Vector & props, const char *label, bool graph_hets);
void GraphDataQuality(PDF & pdf, Vector & het_props, Vector & geno_props, IntArray & unfiltered);
void GraphQualityByIndividual(PDF & pdf, Vector & props, IntArray & unfiltered, const char * label, bool graph_hets);

// Auxilliary graphing functions
//
void DrawSixth(PDF & pdf, PDFChartBasics * chart, int & pr_graph);
void DrawQuadrant(PDF & pdf, PDFChartBasics * chart, int & pr_gender_graph);


#endif



 
 
 
 
 
 
 
 
 
