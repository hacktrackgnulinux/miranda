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
