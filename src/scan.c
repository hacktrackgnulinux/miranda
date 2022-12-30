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
#include <ctype.h>
#include "utils.h"
#include "miranda.h"

/* Load Sequences and Set-up the Alignment Run */

int
find_targets (FILE * fp1, FILE * fp2, FILE * fpout, char *filename)
{

  /* The three key alignment matrices */
  int **matrix1;		/* Best score of all three states (nt-nt,nt-gap,gap-nt  */
  int ***matrix2;		/* Traceback Matrix           */
  int **matrix3;		/* best score for state nt-nt */
  int **matrix4;		/* best score for state gap-nt */
  int **matrix5;		/* best score for state nt-gap */
  int **nt_nt_score;            /* nt nt match matrix for easy lookup */

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

      /* Reverse the wquery (microRNA) sequence */
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
	  seqlen1_global=seqlen1;
          seqlen2 = strlen (reference);

          length_3p_for_weighting=seqlen1-length_5p_for_weighting;
          overlap_cutoff=(seqlen1);

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
	  matrix1 = calloc ((seqlen1 + 1), sizeof (int *));
	  matrix2 = (int ***) calloc (4, sizeof (int **));
	  matrix3 = calloc ((seqlen1 + 1), sizeof (int *));
	  matrix4 = calloc ((seqlen1 + 1), sizeof (int *));
	  matrix5 = calloc ((seqlen1 + 1), sizeof (int *));
          nt_nt_score=calloc ((seqlen1 + 1), sizeof (int *));
	  
	for (i=0; i <4;i++){
	matrix2[i]= (int **) calloc ((seqlen1 + 1), sizeof (int *));
	}


	  for (i = 0; i < seqlen1 + 1; i++)
	    {
	      matrix1[i] = calloc ((seqlen2 + 1), sizeof (int));

	      matrix2[0][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
	      matrix2[1][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
	      matrix2[2][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
	      matrix2[3][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));

	      matrix3[i] = calloc ((seqlen2 + 1), sizeof (int));
	      matrix4[i] = calloc ((seqlen2 + 1), sizeof (int));
	      matrix5[i] = calloc ((seqlen2 + 1), sizeof (int));
	      nt_nt_score[i] = calloc ((seqlen2 + 1), sizeof (int));
	      matrix1[i][0] =  matrix3[i][0] = nt_nt_score[i][0]=0;

               matrix2[0][i][0] =matrix2[1][i][0]=matrix2[2][i][0]=matrix2[3][i][0]=0;
	       matrix4[i][0]=0;matrix5[i][0]=0;

/*	       matrix4[i][0]=-10000;matrix5[i][0]=gap_open+(i*gap_extend); */
            }

	  for (i = 0; i < seqlen2 + 1; i++)
	    {
	      matrix1[0][i] = matrix3[0][i] = nt_nt_score[0][i]=0;
	      matrix4[0][i] =0;matrix5[0][i] =0;
/*	      matrix4[0][i] =gap_open+(i*gap_extend);matrix5[0][i] =-10000; */
               matrix2[0][0][i] =matrix2[1][0][i]=matrix2[2][0][i]=matrix2[3][0][i]=0;
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
          	   /*	  dist[i] =
		    build_matrix_quick (matrix1, matrix2, query, reference2,
					seqlen1, seqlen2);

                  this functionality is disabled as a quick fix  */
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
	    do_alignment (matrix1, matrix2, matrix3, matrix4, matrix5, nt_nt_score,query, reference, scores,
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
	      free (matrix2[0][i]);
	      free (matrix2[1][i]);
	      free (matrix2[2][i]);
	      free (matrix2[3][i]);
	      free (matrix3[i]);
	      free (matrix4[i]);
	      free (matrix5[i]);
	      free (nt_nt_score[i]);
	    }
            free (matrix2[0]);
            free (matrix2[1]);
            free (matrix2[2]);
            free (matrix2[3]);

	  free (matrix1);
	  free (matrix2);
	  free (matrix3);
	  free (matrix4);
	  free (matrix5);
	  free (nt_nt_score);
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
do_alignment (int **best, int ***track, int **a_nt_nt, int **b_gap_nt, int **c_nt_gap, int **nt_nt_score, 
char *query, char *reference,score_struct * scores, hit_struct * hit, int seqlen1,int seqlen2, int verbose, 
final_score * finalscore,int direction, char *query_id, char *reference_id, FILE * fpout)
{

  int i = 0;
  int j = 0;
  int z = 0;
  int utr_offset3p;
  int utr_offset5p;
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
  int hit_score;
  int tmp_integer;
int *good_ones_starts_j, *good_ones_ends_j,good_ones_count;
good_ones_count=-1;

good_ones_starts_j=(int *) calloc (seqlen2, sizeof (int));
good_ones_ends_j=(int *) calloc (seqlen2, sizeof (int));

  total_hits = 0;
  finalscore->no_hits = 0;
  finalscore->max_hit = 0;
  finalscore->max_score = 0;
  finalscore->scan_score = 0;
  finalscore->total_score = 0;
  finalscore->positional[0] = '\0';

  get_nt_nt_seq_scores (nt_nt_score,query,reference,seqlen1,seqlen2);

  build_matrix (best, track, a_nt_nt, b_gap_nt,c_nt_gap,nt_nt_score, query, 
reference, seqlen1, seqlen2, scores);

  qsort (scores, total_hits, sizeof (score_struct), cmpscores);

  for (i = 0; i <=total_hits; i++)
    {
      utr_offset3p=0;
      utr_offset5p=0;
      good_call = 1;
      clear_hit (hit, seqlen1, seqlen2);
      hit_score=scores[i].score;
      if (hit_score > score_threshold)
	{
/*
printf ("score i j path ->> %d %d %d %d\n",scores[i].score,scores[i].i,scores[i].j,scores[i].path); 
continue;*/

	  traceback (best, track, query, reference, scores[i].i, scores[i].j, 1,
		     hit, hit_score);

	good_call=testfor_overlap (good_ones_starts_j,good_ones_ends_j,&good_ones_count,
		  hit->ref_start,hit->ref_end);
                  if (good_call==1){
                   good_ones_starts_j[good_ones_count]=hit->ref_start;
                   good_ones_ends_j[good_ones_count]=hit->ref_end;
                  }                  
	  if (hit->query_start >= 1)
	    {
	      for (j = 0; j <= hit->query_start - 1; j++)
		{
		  diff = hit->query_start - j;
		  hit->rest[0][j] = tolower(query[j]);
		  tmp_integer=hit->ref_start - diff;
                 if (tmp_integer >=0){
		  hit->rest[1][j] = tolower(reference[tmp_integer]);
                  utr_offset3p++;
                 }else{

                 hit->rest[1][j] = '-';
                 }

		  hit->rest[2][j] = ' ';
		}
	    } /* if hit->query_start >= 1) */

	  if ((hit->query_end) < seqlen1)
	    {
	      for (j = hit->query_end; j < seqlen1; j++)
		{
		  diff = j - hit->query_end;
		  hit->rest[3][j - hit->query_end] = tolower(query[j]);
		 tmp_integer=hit->ref_end + diff;
                 if (tmp_integer < seqlen2){
		  hit->rest[4][j - hit->query_end] =
		    tolower(reference[tmp_integer]);
                    utr_offset5p++;
                 } else{
                  hit->rest[4][j - hit->query_end] = '-';
                 } /* if (tmp_integer < */
		  hit->rest[5][j - hit->query_end] = ' ';
		}
	    }/* if hit->query_end <seqlen1) */

/* Adjusting for offset due to local alignment in next two lines */

          hit->ref_end+=(utr_offset5p-1);
          hit->ref_start-=utr_offset3p;

	  strop1[0] = '\0';
	  strop2[0] = '\0';
	  sprintf (strop1, "%s%s%s", hit->rest[3], hit->alignment[0],
		   hit->rest[0]);
	  sprintf (strop2, "%s%s%s", hit->rest[5], hit->alignment[1],
		   hit->rest[2]);

/*
printf ("   QueTE:    3' %s5'\n                %s\n   RTT:      5' %s 3'\n\n",
           hit->alignment[0], 
           hit->alignment[1], 
           hit->alignment[2]);
*/

/*printf("strop1 %s:\n",strop1);
printf("strop2 %s:\n",strop2);
*/
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
		    }/* if  mypos >= 1) && (mypos <= 3)) */

		  if ((mypos >= 2) && (mypos <= 11))
		    {
		      if (strop2[j] != ' ')
			{
			  count2++;
			}
		    }  /* if ((mypos >= 2) && (mypos <= 11) */

		  if ((mypos >= 8) && (mypos <= (seqlen1 - 5)))
		    {
		      if (strop2[j] == ' ')
			{
			  count3++;
			}
		    } /* ((mypos >= 8) && (mypos <= (seqlen1 - 5 */

		  if ((mypos >= (seqlen1 - 4)) && (mypos <= (seqlen1)))
		    {
		      if (strop2[j] != ' ')
			{
			  count4++;
			}
		    }
		}/* for j = 0; j < strlen (str */

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

		}  /* if !nomodel) 1 */
	    } /* if !nomodel) 2 */

	  identity=0;
	  for (j = 0; j < strlen (hit->alignment[0]); j++)
	    {
	      if (hit->alignment[1][j] == '|')
		{
		  identity++;
		}
	    }/* for (j = 0; j < strlen (hit->alignment[0]); j++) */

	  identity = (identity / strlen (hit->alignment[0])) * 100;

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

			} /* for (z = 0; z < valid_hits; z++) */


		    } /* if (valid_hits == 0) */


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
		      printhit (query_id, reference_id, hit, query, reference,
				direction, z_score, energy, fpout);
                     
		    } /* if (good_call */
		} /* if ((energy < energy_threshold) && (z_score >= z_threshold)) */
	} /* if (scores[i].score > score_threshold) */
      scores[i].score = scores[i].path = scores[i].i = scores[i].j = 0;
    } /* for i =0 to total hits */

free (good_ones_starts_j);
free (good_ones_ends_j);

  return (scan_score);
}
