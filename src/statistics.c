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
