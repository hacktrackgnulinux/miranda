/* matrix.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdio.h>
#include "utils.h"
#include "miranda.h"


/* Set the contents of a matrix to zeros */

void
clear_matrix (double **m1, int i1, int j1, int i2, int j2)
{
  int i = 0;
  int j = 0;

  for (i = i1; i <= i2; i++)
    {
      for (j = j1; j <= j2; j++)
	{
	  m1[i][j] = 0;
	}

    }

}

/* Print out the contents of a double matrix */

int
dump_matrix (int len1, int len2, double **matrix)
{

  int i, j = 0;

  for (i = 0; i <= len1; i++)
    {
      for (j = 0; j <= len2; j++)
	{
	  printf ("%2.2lf ", matrix[i][j]);
	}
      printf ("\n");
    }

  return (1);
}

/* Print out the contents of an integer matrix */

int
dump_matrix2 (int len1, int len2, int **matrix)
{

  int i, j = 0;

  for (i = 0; i <= len1; i++)
    {
      for (j = 0; j <= len2; j++)
	{
	  printf ("%d ", matrix[i][j]);
	}
      printf ("\n");
    }

  return (1);
}
/* miranda.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/* 										  */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details				  */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*									          */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:				  */
/* miranda@cbio.mskcc.org (reaches both).					  */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"

int
main (int argc, char *argv[])
{

  char filename1[200];
  char filename2[200];
  char fileout[200];
  FILE *fp1 = 0;
  FILE *fp2 = 0;
  FILE *fpout = stdout;

  /* Set Default Parameter Values */

  scale = 2.0;			/* The 5' miRNA scaling parameter             */
  nomodel = 0;			/* Strict alignment model on/off              */
  gap_open = -8;		/* Gap-open Penalty                           */
  gap_extend = -2;		/* Gap-extend Penalty                         */
  score_threshold = 50;		/* SW Score Threshold for reporting hits      */
  energy_threshold = -20;	/* Energy Threshold (DG) for reporting hits   */
  verbosity = 1;		/* Verbose mode on/off                        */
  outfile = 0;			/* Dump to file on/off                        */
  truncated = 0;		/* Truncate sequences on/off                  */
  do_shuffle = 0;		/* Generate statistics using seq shuffling    */
  no_energy = 0;		/* Turn off Vienna Energy Calcs - FASTER      */
  average = 0;			/* Some statistics for shuffled searches      */
  stdev = 0;
  z_threshold = 5.0;		/* Z-Score threshold >=           */
  shuffle_window = 10;		/* Size of shuffling window       */
  total_shuffles = 100;		/* Total number of shuffles       */
  uniform = 0;			/* Uniform Shuffling mode on/off  */
  total_hits = 0;		/* Generic counter for alignments */

  /* Command-line parsing begins here */
  parse_command_line (argc, argv, &filename1, &filename2, &fileout);

  /* Now check our input and output files can be accessed / created */

  if ((fp1 = fopen (filename1, "r")) == NULL)
    {
      fprintf (stderr, "Error: Cannot open file %s\n", filename1);
      exit (1);
    }

  if ((fp2 = fopen (filename2, "r")) == NULL)
    {
      fprintf (stderr, "Error: Cannot open file %s\n", filename2);
      exit (1);
    }


  if ((outfile) && ((fpout = fopen (fileout, "w")) == NULL))
    {
      fprintf (stderr, "Error: Cannot create output file %s\n", fileout);
      exit (1);
    }

  if (verbosity)
    {
      print_parameters (filename1, filename2, fpout);
    }

  /* Everything looks good.... Start the Scan! */

  find_targets (fp1, fp2, fpout, filename2);
  exit (0);
}
/* output.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"

  /* Command-line parsing begins here */

int
parse_command_line (int argc, char *argv[], char *filename1, char *filename2,
		    char *fileout)
{

  int i = 0;
  char *endptr;

  for (i = 0; i < argc; i++)
    {

      if ((!strcmp (argv[i], "--version")) || (!strcmp (argv[i], "-v"))
	  || (!strcmp (argv[i], "--license"))
	  || (!strcmp (argv[i], "-license")))
	{
	  print_banner (stdout);
	  print_license (stdout);
	  exit (0);
	}

      if ((!strcmp (argv[i], "--help")) || (!strcmp (argv[i], "-h"))
	  || (!strcmp (argv[i], "--h")) || (!strcmp (argv[i], "-help")) || (!strcmp (argv[i], "-usage")))
	{
	  print_options ();
	  exit (0);
	}
    }

  if (argc > 2)
    {

      /* This should contain a microRNA FASTA Sequence (query) */
      strcpy (filename1, argv[1]);

      /* This should contain UTR FASTA Sequence(s) (reference) */
      strcpy (filename2, argv[2]);

      for (i = 3; i < argc; i++)
	{
	  if (!strcmp (argv[i], "-s") && (argc > i + 1))
	    {
	      total_shuffles = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-go") && (argc > i + 1))
	    {
	      gap_open = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-ge") && (argc > i + 1))
	    {
	      gap_extend = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-scale") && (argc > i + 1))
	    {
	      scale = strtod (argv[i + 1], &endptr);
	    }

	  if (!strcmp (argv[i], "-shuffle"))
	    {
	      do_shuffle = 1;
	    }

	  if (!strcmp (argv[i], "-noenergy"))
	    {
	      no_energy = 1;
	    }

	  if (!strcmp (argv[i], "-loose"))
	    {
	      nomodel = 1;
	    }

	  if (!strcmp (argv[i], "-w") && (argc > i + 1))
	    {
	      shuffle_window = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-out") && (argc > i + 1))
	    {
	      strcpy (fileout, argv[i + 1]);
	      outfile = 1;
	    }


	  if (!strcmp (argv[i], "-en") && (argc > i + 1))
	    {
	      energy_threshold = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-sc") && (argc > i + 1))
	    {
	      score_threshold = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-z") && (argc > i + 1))
	    {
	      z_threshold = strtod (argv[i + 1], &endptr);
	    }


	  if (!strcmp (argv[i], "-trim") && (argc > i + 1))
	    {
	      truncated = atoi (argv[i + 1]);
	    }

	  if (!strcmp (argv[i], "-uniform"))
	    {
	      uniform = 1;
	    }

	  if (!strcmp (argv[i], "-quiet"))
	    {
	      verbosity = 0;
	    }
	}

      if (!outfile)
	{
	  /* Print the GPL Friendly Banner */
	  print_banner (stdout);
	  print_small_license (stdout);
	}

  } else
    {

      /* No input, so print banner AND usage, then quit */
      print_banner (stdout);
      print_small_license (stdout);
      print_usage (stdout);
      exit (0);
    }

  return (1);
}

void
print_license (FILE * fpout)
{
  fprintf
    (fpout,
     "   This program is free software; you can redistribute it and/or modify\n");
  fprintf (fpout,
	   "   it under the terms of the GNU General Public License as published by\n");
  fprintf (fpout,
	   "   the Free Software Foundation; either version 2 of the License, or (at\n");
  fprintf (fpout, "   your option) any later version.\n");
  fprintf (fpout, "\n");
  fprintf
    (fpout,
     "   This program is distributed in the hope that it will be useful,\n");
  fprintf (fpout,
	   "   but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
  fprintf (fpout,
	   "   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU\n");
  fprintf (fpout, "   General Public License for more details.\n");
  fprintf (fpout, "\n");
  fprintf
    (fpout,
     "   You should have received a copy of the GNU General Public License\n");
  fprintf (fpout,
	   "   along with this program; if not, write to the Free Software\n");
  fprintf (fpout,
	   "   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307\n");
  fprintf (fpout, "   USA\n\n");
}

void
print_small_license (FILE * fpout)
{
  fprintf (fpout, "   %s comes with ABSOLUTELY NO WARRANTY;\n", PACKAGE);
  fprintf
    (fpout,
     "   This is free software, and you are welcome to redistribute it\n");
  fprintf (fpout,
	   "   under certain conditions; type `miranda --license' for details.\n\n");
}

/* Print out the banner, version and GPL information */

void
print_banner (FILE * fpout)
{
  fprintf (fpout, "\n\n");
  fprintf
    (fpout,
     "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
  fprintf (fpout, "%s v%s    microRNA Target Scanning Algorithm\n", PACKAGE,
	   VERSION);
  fprintf (fpout,
	   "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
  fprintf (fpout, "(c) 2003 Memorial Sloan-Kettering Cancer Center, New York\n");
  fprintf (fpout, "\nAuthors: Anton Enright, Bino John, Chris Sander and Debora Marks\n");
  fprintf (fpout, "(mirnatargets@cbio.mskcc.org - reaches all authors)\n");
  fprintf (fpout, "\nSoftware written by: Anton Enright\n");
  fprintf (fpout, "Distributed for anyone to use under the GNU Public License (GPL),\n");
  fprintf (fpout, "See the files \'COPYING\' and \'LICENSE\' for details\n");
  fprintf (fpout, "\n");
  fprintf (fpout, "If you use this software please cite:\n");
  fprintf (fpout,
	   "Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;\n");
  fprintf (fpout, "(2003) Genome Biology; 5(1):R1.\n");
  fprintf (fpout, "\n");
}

/* When no input is given print out program usage */

void
print_usage ()
{ 

  printf ("miRanda is an miRNA target scanner which aims to predict mRNA\n");
  printf ("targets for microRNAs using dynamic-programming alignment and\n");
  printf ("thermodynamics.\n\n");
  printf ("Usage:\tmiranda query.fasta reference.fasta\n");
  printf ("\nWhere:\n\t\'query\' is a FASTA file with a microRNA query\n");
  printf
    ("\t\'reference\' is a FASTA file containing the sequence(s)\n\tto be scanned.\n\n");
}


/* Routine to print out hit alignments and information from the hit_struct */

void
print_options ()
{

  char *inttobool[2][3];
  char *inttoboolr[2][3];
  strcpy ((char *) inttobool[0], "off");
  strcpy ((char *) inttobool[1], "on");
  strcpy ((char *) inttoboolr[1], "off");
  strcpy ((char *) inttoboolr[0], "on");


  print_banner (stdout);
  print_small_license (stdout);
  print_usage ();
  printf ("OPTIONS\n\n");
  printf (" --help -h\tDisplay this message\n");
  printf (" --version -v\tDisplay version information\n");
  printf (" --license\tDisplay license information\n");
  printf ("\nCore algorithm parameters:\n");
  printf (" -sc S\t\tSet score threshold to S\t\t[DEFAULT: %3.1lf]\n",
	  score_threshold);
  printf
    (" -en -E\t\tSet energy threshold to -E kcal/mol\t[DEFAULT: %3.1lf]\n",
     energy_threshold);
  printf (" -scale Z\tSet scaling parameter to Z\t\t[DEFAULT: %3.1lf]\n",
	  scale);
  printf (" -loose\t\tRemove strict duplex heuristics\t\t[DEFAULT: %s]\n",
	  (char *) inttobool[nomodel]);
  printf ("\nAlignment parameters:\n");
  printf (" -go -X\t\tSet gap-open penalty to -X\t\t[DEFAULT: %3.1lf]\n",
	  gap_extend);
  printf (" -ge -X\t\tSet gap-extend penalty to -X\t\t[DEFAULT: %3.1lf]\n",
	  gap_open);
  printf ("\nGeneral Options:\n");
  printf (" -out file\tOutput results to file\t\t\t[DEFAULT: %s]\n",
	  (char *) inttobool[outfile]);
  printf (" -quiet\t\tDo not output alignments\t\t[DEFAULT: %s]\n",
	  (char *) inttoboolr[verbosity]);
  printf (" -trim T\tTrim reference sequences to T nt\t[DEFAULT: %s]\n",
	  (char *) inttobool[truncated]);
  printf (" -noenergy\tDo not perform thermodynamics\t\t[DEFAULT: %s]\n",
	  (char *) inttobool[no_energy]);
  printf ("\nGenerating statistics from sequence shuffling:\n");
  printf
    (" -shuffle\tGenerate statistics using seq shuffling\t[DEFAULT: %s]\n\t\tNote: This is much slower than a normal scan\n",
     (char *) inttobool[do_shuffle]);
  printf (" -s\t\tTotal number of shuffles to perform\t[DEFAULT: %d]\n",
	  total_shuffles);
  printf (" -w\t\tShuffle window size\t\t\t[DEFAULT: %d]\n", shuffle_window);
  printf (" -uniform\tUniform shuffle instead of windowed\t[DEFAULT: %s]\n",
	  (char *) inttobool[uniform]);
  printf (" -z Z\t\tZ-Score threshold\t\t\t[DEFAULT: %3.1lf]\n", z_threshold);
  printf ("\n\n");
  printf ("This software will be further developed under the open source model,\n");
  printf ("coordinated by Anton Enright and Chris Sander (miranda@cbio.mskcc.org).\n");
  printf ("\nPlease send bug reports to: miranda@cbio.mskcc.org.\n\n");
}



void
print_parameters (char *filename1, char *filename2, FILE * fpout)
{

  if (outfile)
    {
      print_banner (fpout);
      print_small_license (fpout);
    }
  /* Display current parameter settings */
  fprintf (fpout, "Current Settings:\n");
  fprintf
    (fpout,
     "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
  fprintf (fpout, "Query Filename:\t\t%s\n", filename1);
  fprintf (fpout, "Reference Filename:\t%s\n", filename2);
  fprintf (fpout, "Gap Open Penalty:\t%lf\nGap Extend:\t\t%lf\n", gap_open,
	   gap_extend);

  fprintf (fpout, "Score Threshold\t\t%lf\n", score_threshold);
  fprintf (fpout, "Energy Threshold\t%lf kcal/mol\n", energy_threshold);

  if (do_shuffle)
    {
      fprintf (fpout, "Z-Score Threshold\t%lf\n", z_threshold);
      fprintf (fpout, "\n");
      fprintf (fpout, "Shuffling Turned on:\n");
      fprintf (fpout, "Shuffles:\t\t%d\n", total_shuffles);
      if (!uniform)
	{
	  fprintf (fpout, "Window Size:\t\t%d\n", shuffle_window);
      } else
	{
	  fprintf (fpout, "Uniform Shuffle\t%d\n", uniform);
	}
    }

  fprintf (fpout, "Scaling Parameter:\t%lf\n", scale);
  fprintf
    (fpout,
     "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");

}


void
printhit (char *query, char *reference, hit_struct * hit, char *sequence1,
	  char *sequence2, int direction, double z_score, double energy,
	  FILE * fpout)
{

  double similarity = 0;
  double identity = 0;
  int alignment_length = 0;
  int i = 0;

  alignment_length = strlen (hit->alignment[0]);

  for (i = 0; i < alignment_length; i++)
    {

      if (hit->alignment[1][i] == '|')
	{
	  similarity++;
	  identity++;
	}
      if (hit->alignment[1][i] == ':')
	{
	  similarity++;
	}
    }

  similarity = (similarity / (double) alignment_length) * 100;
  identity = (identity / (double) alignment_length) * 100;

  if (direction == FORWARD)
    {
      fprintf
	(fpout,
	 "\n   Forward:\tScore: %lf  Q:%d to %d  R:%d to %d Align Len (%d) (%3.2lf%%) (%3.2lf%%)\n\n",
	 hit->score, hit->query_start + 1, hit->query_end + 1,
	 hit->ref_start + 1, hit->ref_end + 1, alignment_length, identity,
	 similarity);
    }
  if (direction == REVERSE)
    {
      fprintf (fpout, "\n   Reverse:\tScore: %lf  Q:%d to %d  R:%d to %d\n\n",
	       hit->score, hit->query_end + 1, hit->query_start + 1,
	       hit->ref_start + 1, hit->ref_end + 1);
    }

  revstring (hit->alignment[0]);
  revstring (hit->alignment[1]);
  revstring (hit->alignment[2]);

  fprintf (fpout,
	   "   Query:    3' %s%s%s 5'\n                %s%s%s\n   Ref:      5' %s%s%s 3'\n\n",
	   hit->rest[0], hit->alignment[0], hit->rest[3], hit->rest[2],
	   hit->alignment[1], hit->rest[5], hit->rest[1], hit->alignment[2],
	   hit->rest[4]);

  if (do_shuffle)
    {
      fprintf (fpout, "   Z-Score: %2.3lf\n", z_score);
  } else
    {
      z_score = 0;
    }
  if (!no_energy)
    {
      fprintf (fpout, "   Energy:  %lf kCal/Mol\n", energy);
  } else
    {
      energy = 0;
    }

  fprintf (fpout, "\nScores for this hit:\n");
  fprintf (fpout,
	   ">%s\t%s\t%2.2lf\t%2.2lf\t%2.2lf\t%d %d\t%d %d\t%d\t%3.2lf%%\t%3.2lf%%\n\n",
	   query, reference, hit->score, energy, z_score,
	   hit->query_start + 1, hit->query_end + 1, hit->ref_start + 1,
	   hit->ref_end + 1, alignment_length, identity, similarity);


}
/* scan.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "miranda.h"


/* Load Sequences and Set-up the Alignment Run */

int
find_targets (FILE * fp1, FILE * fp2, FILE * fpout, char *filename)
{

  /* The three key alignment matrices */
  double **matrix1;		/* Core Scoring Matrix        */
  int **matrix2;		/* Traceback Matrix           */
  int **matrix3;		/* Sub-optimal path heirarchy */
  int processed = 0;
  int i = 0;
  int seqlen1;
  int seqlen2;
  hit_struct hit;
  hit_struct *hit_ptr;
  score_struct *scores;
  final_score finalscore;
  final_score *fscore_ptr;

  /* Sequence Information, IDs and Descriptions */
  char *query;			/* The Query Sequence (miRNA)             */
  char *query_des;
  char *query_id;
  char *reference;		/* The Reference Sequence (UTR)           */
  char *reference2;		/* A Second Copy of the Ref for Shuffling */
  char *reference_des;
  char *reference_id;
  char *rev_query;		/* Another copy of the query sequence     */

  /* Scoring Information */

  double end_score;
  double maximum = 0;
  double *dist;

  /* File IO            */
  long stream_pos_q = 0;
  long current_pos_q = 0;
  long stream_pos_r = 0;
  long current_pos_r = 0;

  hit_ptr = &hit;
  fscore_ptr = &finalscore;

  /* Memory Allocation for Sequences */
  query = (char *) calloc (100000, sizeof (char));
  query_des = (char *) calloc (1000, sizeof (char));
  query_id = (char *) calloc (1000, sizeof (char));
  reference = (char *) calloc (100000, sizeof (char));
  reference_des = (char *) calloc (1000, sizeof (char));
  reference2 = (char *) calloc (100000, sizeof (char));
  reference_id = (char *) calloc (1000, sizeof (char));
  rev_query = (char *) calloc (100000, sizeof (char));

  /* Array to store distribution of shuffled alignments */
  dist = (double *) calloc (total_shuffles, sizeof (double));

  /* Prepare the generic base lookup array */
  initialize_bases ();

  /* Read the query sequence(s) (microRNA(s)) from a FASTA file */
  while ((current_pos_q =
	  readinseq (stream_pos_q, fp1, query, query_des, query_id)))
    {

      if (verbosity)
	{
	  fprintf (fpout, "Read Sequence:%s %s(%d nt)\n", query_id, query_des,
		   (int) strlen (query));
	}


      /* We are doing alignments like this:
       * 
       *           microRNA
       *   3'-<<<<<<<<<<<<<<<<<<<-5'
       *        |||o|||||  ||||| 
       *   5'->>>>>>>>>>>>>>>>>>>-3'
       *      Reference Sequence 
       *
       *
       * Hence we should reverse one of the two sequences 
       */

      /* Reverse the query (microRNA) sequence */
      strcpy (rev_query, query);
      revstring (query);

      /* Loop over all reference sequences in FASTA file */
      /* Do full scan for each                           */

      fclose (fp2);
      if ((fp2 = fopen (filename, "r")) == NULL)
	{
	  fprintf (stderr, "Error: Cannot open file %s\n", filename);
	  exit (1);
	}
      stream_pos_r = 0;

      while ((current_pos_r =
	      readinseq (stream_pos_r, fp2, reference, reference_des,
			 reference_id)))
	{
	  /* Keep track of the number of sequences scanned so far */
	  processed++;

	  if (verbosity)
	    {
	      fprintf (fpout, "Read Sequence:%s %s(%d nt)\n", reference_id,
		       reference_des, (int) strlen (reference));
	    }

	  if (truncated)
	    {
	      reference[truncated] = '\0';
	    }

	  /* Get sequence lengths for query and reference */
	  seqlen1 = strlen (query);
	  seqlen2 = strlen (reference);
	  strcpy (reference2, reference);

	  /* Initialize the hit / alignment constructs for this sequence */
	  hit.alignment[0] =
	    (char *) calloc (seqlen1 + seqlen2, sizeof (char));
	  hit.alignment[1] =
	    (char *) calloc (seqlen1 + seqlen2, sizeof (char));
	  hit.alignment[2] =
	    (char *) calloc (seqlen1 + seqlen2, sizeof (char));
	  hit.rest[0] = (char *) calloc (30, sizeof (char));
	  hit.rest[1] = (char *) calloc (30, sizeof (char));
	  hit.rest[2] = (char *) calloc (30, sizeof (char));
	  hit.rest[3] = (char *) calloc (30, sizeof (char));
	  hit.rest[4] = (char *) calloc (30, sizeof (char));
	  hit.rest[5] = (char *) calloc (30, sizeof (char));


	  /* Structure for sub-optimal score list */
	  scores =
	    (score_struct *) calloc (seqlen1 * seqlen2,
				     sizeof (score_struct));

	  /* Initialize the three alignment matrices */
	  matrix1 = calloc ((seqlen1 + 1), sizeof (double *));
	  matrix2 = calloc ((seqlen1 + 1), sizeof (int *));
	  matrix3 = calloc ((seqlen1 + 1), sizeof (int *));

	  for (i = 0; i < seqlen1 + 1; i++)
	    {
	      matrix1[i] = calloc ((seqlen2 + 1), sizeof (double));
	      matrix2[i] = calloc ((seqlen2 + 1), sizeof (int));
	      matrix3[i] = calloc ((seqlen2 + 1), sizeof (int));
	      matrix1[i][0] = matrix2[i][0] = matrix3[i][0] = 0;
	    }

	  for (i = 0; i < seqlen2 + 1; i++)
	    {
	      matrix1[0][i] = matrix2[0][i] = matrix3[0][i] = 0;
	    }


	  if (verbosity && do_shuffle)
	    {
	      fprintf
		(fpout,
		 "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	      fprintf (fpout,
		       "Generating Alignment Distribution of Shuffled Sequences\n");
	      fprintf (fpout,
		       "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	    }

	  if (uniform)
	    {
	      shuffle_window = seqlen2;
	    }


	  if (do_shuffle)
	    {
	      irand (0);
	      for (i = 0; i < total_shuffles; i++)
		{
		  end_score = 0;
		  shuffle (reference2, seqlen2, shuffle_window);
		  dist[i] =
		    build_matrix_quick (matrix1, matrix2, query, reference2,
					seqlen1, seqlen2);
		}

	      for (i = 0; i <= total_shuffles; i++)
		{
		  average += (dist[i]);
		  if (dist[i] > maximum)
		    {
		      maximum = dist[i];
		    }
		}
	      average = average / (double) total_shuffles;

	      for (i = 0; i <= total_shuffles; i++)
		{
		  stdev += ((dist[i] - average) * (dist[i] - average));
		}
	      stdev = stdev / (double) (total_shuffles - 1);
	      stdev = sqrt (stdev);
	    }





	  if (verbosity)
	    {
	      if (do_shuffle)
		{
		  fprintf (fpout, "done\t");
		  fprintf (fpout,
			   "Average: %3.2lf\tSt. Dev: %3.2lf\tMax: %3.2lf\n",
			   average, stdev, maximum);
		}

	      fprintf
		(fpout,
		 "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	      fprintf (fpout, "Performing Scan: %s vs %s\n", query_id,
		       reference_id);
	      fprintf (fpout,
		       "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	    }

	  end_score =
	    do_alignment (matrix1, matrix2, matrix3, query, reference, scores,
			  hit_ptr, seqlen1, seqlen2, 1, fscore_ptr, FORWARD,
			  query_id, reference_id, fpout);

	  if (verbosity)
	    {
	      fprintf (fpout, "Score for this Scan:\n");
	    }

	  if (end_score > 0.0)
	    {
	      fprintf
		(fpout,
		 "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions\n");

	      if (!no_energy)
		{
		  fprintf
		    (fpout,
		     ">>%s\t%s\t%2.2lf\t-%2.2lf\t%2.2lf\t%2.2lf\t%d\t%d\t%d\t%s\n",
		     query_id, reference_id, finalscore.total_score,
		     end_score, finalscore.max_score, finalscore.max_hit,
		     processed, seqlen1, seqlen2, finalscore.positional);
	      } else
		{
		  fprintf
		    (fpout,
		     ">>%s\t%s\t%2.2lf\t0.0\t%2.2lf\t0.0\t%d\t%d\t%d\t%s\n",
		     query_id, reference_id, finalscore.total_score,
		     finalscore.max_score, processed, seqlen1, seqlen2,
		     finalscore.positional);

		}
	      fflush (fpout);
	  } else
	    {
	      fprintf (fpout, "No Hits Found above Threshold\n");
	    }


	  if (verbosity)
	    {
	      fprintf (fpout, "Complete\n\n");
	    }

	  fflush (fpout);
	  stream_pos_r = current_pos_r;

	  for (i = 0; i < seqlen1 + 1; i++)
	    {
	      free (matrix1[i]);
	      free (matrix2[i]);
	      free (matrix3[i]);
	    }

	  free (matrix1);
	  free (matrix2);
	  free (matrix3);
	  free (hit.alignment[0]);
	  free (hit.alignment[1]);
	  free (hit.alignment[2]);
	  free (hit.rest[0]);
	  free (hit.rest[1]);
	  free (hit.rest[2]);
	  free (hit.rest[3]);
	  free (hit.rest[4]);
	  free (hit.rest[5]);
	  free (scores);

	}
      stream_pos_q = current_pos_q;
      stream_pos_r = 0;
      current_pos_r = 0;

    }

  fprintf (fpout, "Run Complete\n");
  fflush (fpout);
  if (outfile)
    {
      fclose (fpout);
    }


  return (1);
}

double
do_alignment (double **m1, int **m2, int **m3, char *query, char *reference,
	      score_struct * scores, hit_struct * hit, int seqlen1,
	      int seqlen2, int verbose, final_score * finalscore,
	      int direction, char *query_id, char *reference_id, FILE * fpout)
{

  int i = 0;
  int j = 0;
  int z = 0;
  double energy = 0;
  double scan_score = 0;
  double z_score = 0;
  double identity = 0;
  int hit_cluster[100];
  int valid_hits = 0;
  int good_call = 0;
  int cmin = 0;
  int cmax = 0;
  int diff = 0;
  char strop1[200];
  char strop2[200];
  int count1 = 0;
  int count2 = 0;
  int count3 = 0;
  int count4 = 0;
  int fail = 0;
  int mypos = 0;


  total_hits = 0;
  finalscore->no_hits = 0;
  finalscore->max_hit = 0;
  finalscore->max_score = 0;
  finalscore->scan_score = 0;
  finalscore->total_score = 0;
  finalscore->positional[0] = '\0';


  build_matrix (m1, m2, m3, query, reference, seqlen1, seqlen2, scores);
  qsort (scores, total_hits, sizeof (score_struct), cmpscores);

  for (i = 0; i < total_hits; i++)
    {
      good_call = 1;
      clear_hit (hit, seqlen1, seqlen2);
      if (scores[i].score > score_threshold)
	{
	  traceback (m1, m2, query, reference, scores[i].i, scores[i].j, 1,
		     hit, 0);

	  /* traceback (m1, m2, query, reference, seqlen1, scores[i].j, 1,
	     hit, 0); */

	  if (hit->query_start >= 1)
	    {
	      for (j = 0; j <= hit->query_start - 1; j++)
		{
		  diff = hit->query_start - j;
		  hit->rest[0][j] = query[j];
		  hit->rest[1][j] = reference[hit->ref_start - diff];
		  hit->rest[2][j] = ' ';
		}
	    }

	  if ((hit->query_end) < seqlen1)
	    {
	      for (j = hit->query_end; j < seqlen1; j++)
		{
		  diff = j - hit->query_end;
		  hit->rest[3][j - hit->query_end] = query[j];
		  hit->rest[4][j - hit->query_end] =
		    reference[hit->ref_end + diff];
		  hit->rest[5][j - hit->query_end] = ' ';
		}
	    }

	  strop1[0] = '\0';
	  strop2[0] = '\0';
	  sprintf (strop1, "%s%s%s", hit->rest[3], hit->alignment[0],
		   hit->rest[0]);
	  sprintf (strop2, "%s%s%s", hit->rest[5], hit->alignment[1],
		   hit->rest[2]);

	  mypos = 0;
	  fail = 0;
	  count1 = count2 = count3 = count4 = 0;

	  if (!nomodel)
	    {
	      for (j = 0; j < strlen (strop1); j++)
		{
		  if (strop1[j] != '-')
		    {
		      mypos++;
		    }
		  if ((mypos >= 1) && (mypos <= 3))
		    {
		      if (strop2[j] != ' ')
			{
			  count1++;
			}
		    }
		  if ((mypos >= 2) && (mypos <= 11))
		    {
		      if (strop2[j] != ' ')
			{
			  count2++;
			}
		    }
		  if ((mypos >= 8) && (mypos <= (seqlen1 - 5)))
		    {
		      if (strop2[j] == ' ')
			{
			  count3++;
			}
		    }
		  if ((mypos >= (seqlen1 - 4)) && (mypos <= (seqlen1)))
		    {
		      if (strop2[j] != ' ')
			{
			  count4++;
			}
		    }
		}

	      if (!nomodel)
		{
		  if (count1 < 1)
		    {
		      fail = 1;
		    }
		  if (count2 < 5)
		    {
		      fail = 1;
		    }
		  if (count3 < 1)
		    {
		      fail = 1;
		    }
		  if (count4 < 2)
		    {
		      fail = 1;
		    }

		}
	    }

	  identity=0;
	  for (j = 0; j < strlen (hit->alignment[0]); j++)
	    {
	      if (hit->alignment[1][j] == '|')
		{
		  identity++;
		}
	    }

	  identity = (identity / strlen (hit->alignment[0])) * 100;

	  if ((identity >= 40) && (!fail))
	    {

	      if (do_shuffle){
	          z_score = (hit->score - average) / stdev;
	        } else {
		  z_score = 1000000;
	        }

	      if (!no_energy)
		{
		  energy = get_energy (hit);
	      } else
		{
		  energy = -1000000;
		}


	      if ((energy < energy_threshold) && (z_score >= z_threshold))
		{

		  if (valid_hits == 0)
		    {
		      hit_cluster[valid_hits] = hit->ref_start;
		      valid_hits++;
		      cmax = hit->ref_start;
		      cmin = hit->ref_start;
		  } else
		    {


		      for (z = 0; z < valid_hits; z++)
			{
			  if (hit_cluster[z] > cmax)
			    {
			      cmax = hit_cluster[z];
			    }

			  if (hit_cluster[z] < cmin)
			    {
			      cmin = hit_cluster[z];
			    }

			  if ((hit->ref_start >= (hit_cluster[z] - 25))
			      && (hit->ref_start <= hit_cluster[z]))
			    {
			      good_call = 0;
			    }
			  if ((hit->ref_start <= (hit_cluster[z] + 25))
			      && (hit->ref_start >= hit_cluster[z]))
			    {
			      good_call = 0;
			    }
			}

		      if (hit->ref_start > (cmin + 1000))
			{
			  good_call = 0;
			}
		      if (hit->ref_start < (cmin - 1000))
			{
			  good_call = 0;
			}

		    }

		  if (good_call)
		    {

		      hit_cluster[valid_hits] = hit->ref_start;
		      valid_hits++;

		      scan_score += (energy * -1);
		      finalscore->no_hits++;
		      sprintf (finalscore->positional, "%s %d",
			       finalscore->positional, hit->ref_start);


		      if (energy < finalscore->max_hit)
			{
			  finalscore->max_hit = energy;
			}


		      finalscore->total_score += hit->score;
		      if (hit->score > finalscore->max_score)
			{
			  finalscore->max_score = hit->score;
			}
		      /*fprintf(fpout,"PASS: %d %d %d %d %s\n",count1,count2,count3,count4,strop2); */
		      printhit (query_id, reference_id, hit, query, reference,
				direction, z_score, energy, fpout);
		    }
		}
	    }
	}
      scores[i].score = scores[i].path = scores[i].i = scores[i].j = 0;
    }
  return (scan_score);
}
/* seqio.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"

long
readinseq (long pos, FILE * fp, char *sequence, char *description, char *name)
{

  char c;
  int i;
  int j;
  int flag = 0;

  while ((c = fgetc (fp)) != EOF)
    {

      if (c == '>')
	{
	  /* ungetc(c, fp); */
	  i = 0;
	  j = 0;
	  while ((c = fgetc (fp)) != '\n')
	    {
	      if (c == '\t')
		{
		  c = ' ';
		}

	      if (c == ' ')
		{
		  flag = 1;
		}

	      if (!flag)
		{
		  name[i] = c;
		  i++;
	      } else
		{
		  description[j] = c;
		  j++;
		}

	    }
	  name[i] = '\0';
	  description[j] = '\0';

	}
      i = 0;
      while (((c = fgetc (fp)) != EOF) && (c != '>'))
	{
	  if ((c != '\n') && (c != ' ') && (c != '\r'))
	    {

	      if (((c >= 'a') && (c <= 'z')))
		{
		  c += 'A' - 'a';
		}
	      if (((c >= 'A') && (c <= 'Z')) || (c == '*'))
		{
		  sequence[i] = c;
		  i++;
	      } else
		{
		  sequence[i] = 'X';
		  i++;
		  fprintf (stderr,
			   "Error: Sequence file contains non-standard or lowercase characters %c\n",
			   c);
		}
	    }
	}
      sequence[i] = '\0';
      ungetc (c, fp);
      return (1);
    }

  return (0);
}
/* statistics.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "utils.h"
#include "miranda.h"

void
shuffle (char *sequence, int seqlen, int wsize)
{
  int k = 0;
  char tmp;
  char *top;
  int i, j = 0;
  int mm = seqlen % wsize;

  for (k = 0; k < (seqlen - wsize) + 1; k += wsize)
    {
      top = &sequence[k];
      for (i = wsize; i > 0; i--)
	{
	  j = nrand (i);
	  tmp = top[j];
	  top[j] = top[i - 1];
	  top[i - 1] = tmp;
	}
    }
  top = &sequence[seqlen - mm];
  for (i = mm; i > 0; i--)
    {
      j = nrand (i);
      tmp = top[j];
      top[j] = top[i - 1];
      top[i - 1] = tmp;
    }


}

void
irand (int n)
{				/* initialize random number generator */
  if (n == 0)
    {
      n = time (NULL);
      n = n % 16381;
      if ((n % 2) == 0)
	n++;
    }
  srand48 (n);
}


int
nrand (int n)
{				/* returns a random number between 0 and n-1
				 * where n < 64K) */
  int rn;
  rn = lrand48 ();
  rn = rn >> 16;
  rn = (rn % n);
  return rn;
}

int
getfreq (char *sequence, int seqlen, double *frequency)
{

  int i = 0;
  for (i = 0; i < 256; i++)
    {
      frequency[i] = 0;
    }

  for (i = 0; i < seqlen; i++)
    {
      frequency[toupper (sequence[i])]++;
    }

  for (i = 0; i < 256; i++)
    {
      frequency[i] = frequency[i] / seqlen;
    }

  return (1);
}
/* swat.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"
#include "scmatrix.h"

void
initialize_bases ()
{
  bases['C'] = 0;
  bases['G'] = 1;
  bases['A'] = 2;
  bases['T'] = 3;
  bases['U'] = 4;
  bases['X'] = 5;
  bases['N'] = 5;
  bases['c'] = 0;
  bases['g'] = 1;
  bases['a'] = 2;
  bases['t'] = 3;
  bases['u'] = 4;
  bases['x'] = 5;
  bases['n'] = 5;

}

double
build_matrix_quick (double **m1, int **m2, char *sequence1,
		    char *sequence2, int seqlen1, int seqlen2)
{

  double penalty1 = 0;
  double penalty2 = 0;
  int i, j = 0;
  double maxi = 0;

  for (i = 1; i <= seqlen1; i++)
    {
      for (j = 1; j <= seqlen2; j++)
	{

	  if (m2[i][j - 1] == LEFT)
	    {
	      penalty1 = gap_extend;
	  } else
	    {
	      penalty1 = gap_open;
	    }

	  if (m2[i - 1][j] == UP)
	    {
	      penalty2 = gap_extend;
	  } else
	    {
	      penalty2 = gap_open;
	    }

	  m1[i][j] =
	    max (m1[i - 1][j - 1] +
		 (score (sequence1[i - 1], sequence2[j - 1])),
		 m1[i][j - 1] + (penalty1), m1[i - 1][j] + (penalty2));
	  m2[i][j] = CURR;

	  if (m1[i][j] < 0)
	    {
	      m1[i][j] = 0;
	    }

	  if (m1[i][j] > maxi)
	    {
	      maxi = m1[i][j];
	    }

	}
    }
  return (maxi);
}


void
build_matrix (double **m1, int **m2, int **m3, char *sequence1,
	      char *sequence2, int seqlen1, int seqlen2,
	      score_struct * scores)
{

  double penalty1 = 0;
  double penalty2 = 0;
  int i, j = 0;
  int path = 0;

  for (i = 1; i <= seqlen1; i++)
    {
      for (j = 1; j <= seqlen2; j++)
	{

	  if (m2[i][j - 1] == LEFT)
	    {
	      penalty1 = gap_extend;
	  } else
	    {
	      penalty1 = gap_open;
	    }

	  if (m2[i - 1][j] == UP)
	    {
	      penalty2 = gap_extend;
	  } else
	    {
	      penalty2 = gap_open;
	    }

	  if (i < (strlen (sequence1) - 10))
	    {
	      m1[i][j] =
		max (m1[i - 1][j - 1] +
		     (score (sequence1[i - 1], sequence2[j - 1])),
		     m1[i][j - 1] + (penalty1), m1[i - 1][j] + (penalty2));
	      m2[i][j] = CURR;
	      m3[i][j] = 0;
	  } else
	    {
	      m1[i][j] =
		max (m1[i - 1][j - 1] +
		     (scale * score (sequence1[i - 1], sequence2[j - 1])),
		     m1[i][j - 1] + (scale * penalty1),
		     m1[i - 1][j] + (scale * penalty2));
	      m2[i][j] = CURR;
	      m3[i][j] = 0;
	    }

	  if (m1[i][j] < 0)
	    {
	      m1[i][j] = 0;
	  } else
	    {
	      if (m2[i][j] == DIAG)
		{
		  if (m3[i - 1][j - 1] <= 0)
		    {
		      total_hits++;
		      m3[i][j] = path = total_hits;
		  } else
		    {
		      m3[i][j] = path = m3[i - 1][j - 1];
		    }
	      } else if (m2[i][j] == LEFT)
		{
		  if (m3[i][j - 1] <= 0)
		    {
		      total_hits++;
		      m3[i][j] = path = total_hits;
		  } else
		    {
		      m3[i][j] = path = m3[i][j - 1];

		    }
	      } else
		{
		  if (m3[i - 1][j] <= 0)
		    {
		      total_hits++;
		      m3[i][j] = path = total_hits;
		  } else
		    {
		      m3[i][j] = path = m3[i - 1][j];
		    }
		}

	      if (m1[i][j] >= scores[path].score)
		{
		  scores[path].score = m1[i][j];
		  scores[path].path = path;
		  scores[path].i = i;
		  scores[path].j = j;
		}
	    }
	}
    }
}

void
traceback (double **m1, int **m2, char *sequence1, char *sequence2, int i,
	   int j, int remove, hit_struct * hit_ptr, int length)
{

  double score = 0;

  if (length == 0)
    {
      hit_ptr->query_end = i;
      hit_ptr->ref_end = j;
    }
  if (m1[i][j] > 0)
    /*if (i>0 && j >0) */
    {

      if (m2[i][j] == DIAG)
	{
	  length++;
	  traceback (m1, m2, sequence1, sequence2, i - 1, j - 1, 1, hit_ptr,
		     length);

	  score =
	    (double)
	    match[bases[(int) sequence1[i - 1]]][bases
						 [(int) (sequence2[j - 1])]];
	  if (i >= (strlen (sequence1) - 10))
	    {
	      score = score * scale;
	    }

	  hit_ptr->score += score;
	  hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
	  hit_ptr->alignment[2][length - 1] = sequence2[j - 1];

	  if (score > 2)
	    {
	      hit_ptr->alignment[1][length - 1] = '|';
	  } else if ((score >= 0) && (score <= 4))
	    {
	      hit_ptr->alignment[1][length - 1] = ':';
	  } else
	    {
	      hit_ptr->alignment[1][length - 1] = ' ';
	    }
	  if (remove)
	    {
	      m1[i][j] = 0;
	    }
      } else if (m2[i][j] == UP)
	{

	  if (m2[i - 1][j] == UP)
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_extend;
	      } else
		{
		  hit_ptr->score += (gap_extend * scale);
		}
	  } else
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_open;
	      } else
		{
		  hit_ptr->score += (gap_open * scale);
		}
	    }


	  length++;
	  traceback (m1, m2, sequence1, sequence2, i - 1, j, 1, hit_ptr,
		     length);
	  hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
	  hit_ptr->alignment[1][length - 1] = ' ';
	  hit_ptr->alignment[2][length - 1] = '-';

	  if (remove)
	    {
	      m1[i][j] = 0;
	    }
      } else
	{

	  if (m2[i][j - 1] == LEFT)
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_extend;
	      } else
		{
		  hit_ptr->score += (gap_extend * scale);
		}
	  } else
	    {
	      if (i < (strlen (sequence1) - 10))
		{
		  hit_ptr->score += gap_open;
	      } else
		{
		  hit_ptr->score += (gap_open * scale);
		}
	    }

	  length++;
	  traceback (m1, m2, sequence1, sequence2, i, j - 1, 1, hit_ptr,
		     length);
	  hit_ptr->alignment[0][length - 1] = '-';
	  hit_ptr->alignment[1][length - 1] = ' ';
	  hit_ptr->alignment[2][length - 1] = sequence2[j - 1];

	  if (remove)
	    {
	      m1[i][j] = 0;
	    }
	}
  } else
    {
      hit_ptr->query_start = i;
      hit_ptr->ref_start = j;
      hit_ptr->alignment[0][length] = '\0';
      hit_ptr->alignment[1][length] = '\0';
      hit_ptr->alignment[2][length] = '\0';
    }
}

int
build_sub_matrix (int **matrix)
{

  int i = 0;
  int j = 0;

  for (i = 0; i < 6; i++)
    {
      for (j = 0; j < 6; j++)
	{
	  matrix[toupper (baselist[i])][toupper (baselist[j])] =
	    (int) match[i][j];
	}
    }

  return (1);
}


double
score (char nt1, char nt2)
{
  return (double) match[bases[(int) nt1]][bases[(int) nt2]];
}
/* thermo.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * This file includes headers and links to the RNAlib library of 
 * Ivo Hofackers Vienna Package
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "inverse.h"
#include "RNAstruct.h"
#include "treedist.h"
#include "stringdist.h"
#include "profiledist.h"

double
get_energy (hit_struct * hit)
{
  double energy = 0;
  int i = 0;
  int j = 0;
  char foldsequence[5000];
  revstring (hit->alignment[0]);
  revstring (hit->rest[1]);
  revstring (hit->rest[4]);
  foldsequence[0] = '\0';

  for (i = 0; i < strlen (hit->rest[0]); i++)
    {
      foldsequence[j] = hit->rest[0][i];
      j++;
    }

  for (i = 0; i < strlen (hit->alignment[0]); i++)
    {
      if (hit->alignment[0][i] != '-')
	{
	  foldsequence[j] = hit->alignment[0][i];
	  j++;
	}
    }

  for (i = 0; i < strlen (hit->rest[3]); i++)
    {
      foldsequence[j] = hit->rest[3][i];
      j++;
    }

  for (i = 0; i < 7; i++)
    {
      foldsequence[j] = 'X';
      j++;
    }

  for (i = 0; i < strlen (hit->rest[4]); i++)
    {
      foldsequence[j] = hit->rest[4][i];
      j++;
    }


  for (i = 0; i < strlen (hit->alignment[2]); i++)
    {
      if (hit->alignment[2][i] != '-')
	{
	  foldsequence[j] = hit->alignment[2][i];
	  j++;
	}
    }

  for (i = 0; i < strlen (hit->rest[1]); i++)
    {
      foldsequence[j] = hit->rest[1][i];
      j++;
    }

  foldsequence[j] = '\0';

  /* printf("FOLD: %s\n",foldsequence); */
  energy = vfold (foldsequence);
  revstring (hit->alignment[0]);
  revstring (hit->rest[1]);
  revstring (hit->rest[4]);
  return (energy);

}

double
vfold (char *sequence)
{
  void *struct1;
  double e1;
  struct1 = (char *) space (sizeof (char) * (strlen (sequence) + 1));
  temperature = 30;
  initialize_fold (strlen (sequence));
  e1 = fold (sequence, struct1);
  /*PS_rna_plot (sequence, struct1, "out.ps"); */
  free_arrays ();
  free (struct1);
  return (e1);
}
/* utils.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <config.h>
#include <stdio.h>
#include <string.h>
#include "miranda.h"
#include "utils.h"

int
cmpscores (const void *p1, const void *p2)
{

  score_struct *s1;
  score_struct *s2;

  s1 = (score_struct *) p1;
  s2 = (score_struct *) p2;

  if (s1->score < s2->score)
    {
      return (1);
  } else if (s1->score > s2->score)
    {
      return (-1);
  } else
    {
      return (0);
    }

}

void
clear_hit (hit_struct * hit, int seqlen1, int seqlen2)
{

  hit->score = 0;
  hit->query_start = 0;
  hit->query_end = 0;
  hit->ref_start = 0;
  hit->ref_end = 0;

  memset (hit->alignment[0], '\0', seqlen1 + seqlen2);
  memset (hit->alignment[1], '\0', seqlen1 + seqlen2);
  memset (hit->alignment[2], '\0', seqlen1 + seqlen2);

  memset (hit->rest[0], '\0', 30);
  memset (hit->rest[1], '\0', 30);
  memset (hit->rest[2], '\0', 30);
  memset (hit->rest[3], '\0', 30);
  memset (hit->rest[4], '\0', 30);
  memset (hit->rest[5], '\0', 30);

}

double
max (double a, double b, double c)
{

  if ((a >= b) && (a >= c))
    {
      CURR = DIAG;
      return (a);
  } else if (b >= c)
    {
      CURR = LEFT;
      return (b);
  } else
    {
      CURR = UP;
      return (c);
    }
}

void
revstring (char s[])
{
  int c, i, j;

  for (i = 0, j = strlen (s) - 1; i < j; i++, j--)
    {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
    }
}

int
getbase (int c)
{

  switch ((int) c)
    {
    case 'C':
      return (0);
    case 'G':
      return (1);
    case 'A':
      return (2);
    case 'T':
      return (3);
    case 'U':
      return (4);
    case 'X':
      return (5);
    default:
      return (-1);
    }

}
