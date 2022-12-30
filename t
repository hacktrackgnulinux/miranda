  score_threshold = 50;		/* SW Score Threshold for reporting hits      */
	      score_threshold = atoi (argv[i + 1]);
  printf (" -sc S\t\tSet score threshold to S\t\t[DEFAULT: %3.1lf]\n",
	  score_threshold);
  fprintf (fpout, "Score Threshold\t\t%lf\n", score_threshold);
	  char *sequence2, int direction, double z_score, double energy,
	 hit->score, hit->query_start + 1, hit->query_end + 1,
	       hit->score, hit->query_end + 1, hit->query_start + 1,
      fprintf (fpout, "   Z-Score: %2.3lf\n", z_score);
      z_score = 0;
	   query, reference, hit->score, energy, z_score,
  score_struct *scores;
  final_score finalscore;
  final_score *fscore_ptr;
  double end_score;
  fscore_ptr = &finalscore;
	  /* Structure for sub-optimal score list */
	  scores =
	    (score_struct *) calloc (seqlen1 * seqlen2,
				     sizeof (score_struct));
		  end_score = 0;
	  end_score =
	    do_alignment (matrix1, matrix2, matrix3, query, reference, scores,
			  hit_ptr, seqlen1, seqlen2, 1, fscore_ptr, FORWARD,
	  if (end_score > 0.0)
		     query_id, reference_id, finalscore.total_score,
		     end_score, finalscore.max_score, finalscore.max_hit,
		     processed, seqlen1, seqlen2, finalscore.positional);
		     query_id, reference_id, finalscore.total_score,
		     finalscore.max_score, processed, seqlen1, seqlen2,
		     finalscore.positional);
	  free (scores);
	      score_struct * scores, hit_struct * hit, int seqlen1,
	      int seqlen2, int verbose, final_score * finalscore,
  double scan_score = 0;
  double z_score = 0;
  finalscore->no_hits = 0;
  finalscore->max_hit = 0;
  finalscore->max_score = 0;
  finalscore->scan_score = 0;
  finalscore->total_score = 0;
  finalscore->positional[0] = '\0';
  build_matrix (m1, m2, m3, query, reference, seqlen1, seqlen2, scores);
  qsort (scores, total_hits, sizeof (score_struct), cmpscores);
      if (scores[i].score > score_threshold)
	  traceback (m1, m2, query, reference, scores[i].i, scores[i].j, 1,
	  /* traceback (m1, m2, query, reference, seqlen1, scores[i].j, 1,
	          z_score = (hit->score - average) / stdev;
		  z_score = 1000000;
	      if ((energy < energy_threshold) && (z_score >= z_threshold))
		      scan_score += (energy * -1);
		      finalscore->no_hits++;
		      sprintf (finalscore->positional, "%s %d",
			       finalscore->positional, hit->ref_start);
		      if (energy < finalscore->max_hit)
			  finalscore->max_hit = energy;
		      finalscore->total_score += hit->score;
		      if (hit->score > finalscore->max_score)
			  finalscore->max_score = hit->score;
				direction, z_score, energy, fpout);
      scores[i].score = scores[i].path = scores[i].i = scores[i].j = 0;
  return (scan_score);
		 (score (sequence1[i - 1], sequence2[j - 1])),
	      score_struct * scores)
		     (score (sequence1[i - 1], sequence2[j - 1])),
		     (scale * score (sequence1[i - 1], sequence2[j - 1])),
	      if (m1[i][j] >= scores[path].score)
		  scores[path].score = m1[i][j];
		  scores[path].path = path;
		  scores[path].i = i;
		  scores[path].j = j;
  double score = 0;
	  score =
	      score = score * scale;
	  hit_ptr->score += score;
	  if (score > 2)
	  } else if ((score >= 0) && (score <= 4))
		  hit_ptr->score += gap_extend;
		  hit_ptr->score += (gap_extend * scale);
		  hit_ptr->score += gap_open;
		  hit_ptr->score += (gap_open * scale);
		  hit_ptr->score += gap_extend;
		  hit_ptr->score += (gap_extend * scale);
		  hit_ptr->score += gap_open;
		  hit_ptr->score += (gap_open * scale);
score (char nt1, char nt2)
cmpscores (const void *p1, const void *p2)
  score_struct *s1;
  score_struct *s2;
  s1 = (score_struct *) p1;
  s2 = (score_struct *) p2;
  if (s1->score < s2->score)
  } else if (s1->score > s2->score)
  hit->score = 0;
