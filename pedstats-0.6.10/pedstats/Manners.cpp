////////////////////////////////////////////////////////////////////// 
// pedstats/Manners.cpp 
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
 
#include "Manners.h"

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

void SetupCrashHandlers()
   {
   signal(SIGINT, (signal_handler) UserBreak);
   signal(SIGSEGV, (signal_handler) OutOfMemory);
   signal(SIGABRT, (signal_handler) OutOfMemory);
   }

void OutOfMemory(int)
   {
   printf("\n\nPEDSTATS HAS CRASHED\n\n"
          "The operating system has decided to terminate Pedstats,\n"
          "probably due to a memory access problem.\n\n"
          "To help improve this program, please report bugs by\n"
          "e-mailing wiggie@umich.edu\n\n");

   exit(EXIT_FAILURE);
   }

void UserBreak(int)
   {
   printf("\n\nPEDSTATS STOPPED BY USER\n\n");

   exit(EXIT_FAILURE);
   }

 
 
 
 
 
 
