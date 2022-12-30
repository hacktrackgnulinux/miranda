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

/*----------------------------------------------------------------*/
void get_nt_nt_seq_scores (int **nt_nt_scores,char *sequence1, char *sequence2, int seqlen1, int seqlen2) {
int i,j;

for (i = 1; i <= seqlen1; i++)  {

      for (j = 1; j <= seqlen2; j++){

	if (i <= length_3p_for_weighting  || i == seqlen1){
		nt_nt_scores[i][j]=score (sequence1[i - 1], sequence2[j - 1]);
	} else {
	        nt_nt_scores[i][j]=(scale * score5p (sequence1[i - 1], sequence2[j - 1]));
        } /* if (i <= length_3p_for_weighting */

     }/* for j */
}/* for i */
} /* sub routine get nt_nt_seq_scores*/

/*----------------------------------------------------------------*/

void
build_matrix ( int **best, int ***track, int **a_nt_nt,int **b_gap_nt,int 
**c_nt_gap,int **nt_nt_score, char *sequence1,char *sequence2, int seqlen1, int seqlen2, score_struct * scores)

{
int i,j,max,open,ext,good_count,tmp_a,tmp_b,tmp_c;

long int tmp;

good_count=-1;


  for (i = 1; i <= seqlen1; i++) {
      for (j = 1; j <= seqlen2; j++) {

if (i <= length_3p_for_weighting  || i == seqlen1){ open=gap_open;ext=gap_extend;
}else{open=scale*gap_open;ext=scale*gap_extend;}

	max=max_finder_fourstates( a_nt_nt[i-1][j-1],b_gap_nt[i-1][j-1],c_nt_gap[i-1][j-1]);

	a_nt_nt[i][j]=max + nt_nt_score[i][j];
        track[1][i][j]=CURR;
        if (a_nt_nt[i][j] <=0){a_nt_nt[i][j]=0;track[1][i][j]=0;}

        tmp_a=(a_nt_nt[i][j-1]+open); tmp_b=(b_gap_nt[i][j-1] +ext);
        if (tmp_b > tmp_a){
        b_gap_nt[i][j]=tmp_b;track[2][i][j]=2;
        }else{
	b_gap_nt[i][j]=tmp_a;track[2][i][j]=1;
        }

	tmp_a=(a_nt_nt[i-1][j]+open); tmp_c=(c_nt_gap[i-1][j] +ext);
        if (tmp_c > tmp_a){
         c_nt_gap[i][j]=tmp_c;track[3][i][j]=3;
        }else{
         c_nt_gap[i][j]=tmp_a;track[3][i][j]=1;    
        }
       
        best[i][j]=max_finder_and_track_threestates (a_nt_nt[i][j],b_gap_nt[i][j],c_nt_gap[i][j]);       
        track[0][i][j]=CURR;

        if (best[i][j] > score_threshold){
	          good_count++;
                 scores[good_count].score = best[i][j];
                  scores[good_count].path = good_count;
                  scores[good_count].i = i;
                  scores[good_count].j = j;
        }/* best[i][j] > score_threshold ) */
	    } /* for j */
	} /* for i */

        qsort (scores, (good_count+1), sizeof (score_struct), cmpscores);
        nullify_j_overlaps (good_count,scores);
        qsort (scores, (good_count+1), sizeof (score_struct), cmpscores);

/*
for (index=0;index <=total_hits;index++){
printf ("stuct %d %d %d %d\n",scores[index].score,scores[index].i,scores[index].path,scores[index].j);
}

printf ("\n");
printf ("a_nt_nt->\n");
print2dmatrix(a_nt_nt,seqlen1,seqlen2);

printf ("b_gap_nt->\n");
print2dmatrix(b_gap_nt,seqlen1,seqlen2);

printf ("c_nt_gap->\n");
print2dmatrix(c_nt_gap,seqlen1,seqlen2);

printf ("best->\n");
print2dmatrix(best,seqlen1,seqlen2);
printf ("track->\n");
print2dmatrix(track[0],seqlen1,seqlen2);
*/


} /* sub rooutine build_matrix */
/*----------------------------------------------------------------*/

/*----------------------------------------------------------------*/

void
traceback (int **best, int ***track, char *sequence1, char *sequence2, int i,
	   int j, int remove, hit_struct * hit_ptr, int hit_score) {

      int length,track_array;
      length = 0;
      hit_ptr->query_end = i;
      hit_ptr->ref_end = j;
      hit_ptr->score=hit_score;

track_array=track[0][i][j];

while (best[i][j] >0){
if (track_array ==1){
	  length++;
          hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
          hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
	  hit_ptr->alignment[1][length - 1]=ali_rep[bases[(int) sequence1[i - 1]]][bases[(int) (sequence2[j - 1])]];
	track_array=track[track_array][i][j];
          i--;j--;

} else if (track_array ==2){
	  length++;
	  hit_ptr->alignment[0][length - 1] = '-';
          hit_ptr->alignment[2][length - 1] = sequence2[j - 1];
          hit_ptr->alignment[1][length - 1]=' ';
	track_array=track[track_array][i][j];
          j--;
} else if (track_array ==3){
	  length++;
          hit_ptr->alignment[0][length - 1] = sequence1[i - 1];
          hit_ptr->alignment[2][length - 1] = '-';
          hit_ptr->alignment[1][length - 1]=' ';
	track_array=track[track_array][i][j];
          i--;
}else {
      hit_ptr->alignment[0][length] = '\0';
      hit_ptr->alignment[1][length] = '\0';
      hit_ptr->alignment[2][length] = '\0';
     break;
}

}/*while */

      hit_ptr->query_start = i;
      hit_ptr->ref_start = j;


} /* traceback function */
/*----------------------------------------------------------------*/
int testfor_overlap (int *good_ones_starts_j,int *good_ones_ends_j,int *good_ones_count,
                 int test_start,int test_end){
int index,good_call;
 int ref_start,ref_end,min_end,max_start,overlap;
good_call=1;

if (*good_ones_count <0){good_call=1;*good_ones_count=*good_ones_count+1;return(good_call);}

for (index=0;index <=  *good_ones_count;index++){
ref_start=good_ones_starts_j[index];
ref_end=good_ones_ends_j[index];

min_end= ref_end < test_end ? ref_end : test_end;
max_start=ref_start > test_start ? ref_start : test_start;

overlap=(min_end-max_start);
if (overlap >5){
	good_call=0;
break;
}/* if overlap */
}/* for index */
if (good_call==1){*good_ones_count=*good_ones_count+1;}

/*
printf ("ref_s ref_e test_s test_e min_end max_start good_call good_ones_count %d %d %d %d %d %d %d %d\n",
ref_start,ref_end,test_start,test_end,min_end,max_start,good_call,*good_ones_count);
*/

return(good_call);
} /* subroutine  nullify_overlap */
/*----------------------------------------------------------------*/



void nullify_j_overlaps (int total_elements,score_struct * scores){

int index;
int outer_j;
int inner_index;

for (index=0;index <=total_elements;index++){
int score_inner=scores[index].score;

if (score_inner ==0){continue;}
total_hits++;
outer_j=scores[index].j;

	for (inner_index=total_elements;inner_index >index;inner_index--){
             int inner_j=scores[inner_index].j;
             int diff;diff=abs(inner_j-outer_j);
             if (diff <=overlap_cutoff) {
	      scores[inner_index].score=0;
	     }
        }/* for (inner_index=total_elements;inner_index >index;inner_index--) */
}/* index total_elements */
} /* void nullify_j_overlaps */
/*----------------------------------------------------------------*/



double
build_matrix_quick (int **m1, int **m2, char *sequence1,
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


int score (char nt1, char nt2)
{
  return (int) match[bases[(int) nt1]][bases[(int) nt2]];
}

int score5p (char nt1, char nt2)
{
  return (int) match5p[bases[(int) nt1]][bases[(int) nt2]];
}

