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

/* ------------------------------------------------*/

int max_finder_fourstates (a, b, c)
{

int max;

  if ((a >= b) && (a >= c))
    {
      max=a;
      CURR=1; 
  } else if (b >= c)
    {
      max=b;
      CURR=2;
  } else
    {
      max =c;
      CURR=3;
    }

return(max);
}/* max finder sub routine*/
/* ------------------------------------------------*/

/* ------------------------------------------------*/

int max_finder_threestates (a, b)
{
int max;
max = a >=b ? a:b;
/*if (max <0){max =0;}*/
return(max);
}/* max finder sub routine*/
/* ------------------------------------------------*/


int max_finder_and_track_threestates (a, b, c)
{

int max;

  if ((a >= b) && (a >= c))
    {
      max=a; 
      CURR=1;
  } else if (b >= c)
    {
      max=b;
      CURR=2;
  } else
    {
      max =c;
      CURR=3;
    }

if (max <=0){max =0;CURR=0;}
return(max);
}/* max finder sub routine*/
/* ------------------------------------------------*/



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

void print2dmatrix (int **m,int seqlen1,int seqlen2) {

int i,j;

  for (i = 0; i <= seqlen1; i++) {
      for (j = 0; j <= seqlen2; j++) {
       printf ("%+4d ",m[i][j]);
     } /* for j */
	printf ("\n");
} /* for i */

}


void print1dmatrix (int *m,int seqlen1) {

int i;

  for (i = 0; i <= seqlen1; i++) {
       printf ("%+4d ",m[i]);
  } /* for i */
	printf ("\n");
}


