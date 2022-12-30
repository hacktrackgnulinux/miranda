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
