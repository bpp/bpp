/*
    Copyright (C) 2016-2018 Tomas Flouri, Bruce Rannala and Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp.h"

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

/* used for creating FigTree.tre */
typedef struct nodepinfo_s 
{
  /* 95% HPD CI */
  double lo;
  double hi;

  /* mean age */
  double age;

  /* mean theta */
  double theta;
} nodepinfo_t; 

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static long get_long(const char * line, long * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%ld%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static long get_double(const char * line, double * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%lf%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  return *x > *y;
}

static double eff_ict(double * y, long n, double mean, double stdev)
{
  /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
     sequence method */

  long i,j;
  double tint = 1;
  double rho, rho0 = 0;

  /* TODO: ADDED NOW */
  double * x = (double *)xmalloc((size_t)n * sizeof(double));
  for (i = 0; i < n; ++i)
    x[i] = (y[i]-mean)/stdev;


  if (stdev/(fabs(mean)+1) < 1E-9)
  {
   tint = n;
  }
  else
  {
    for (i = 1; i < n-10; ++i)
    {
      rho = 0;
      for (j = 0; j < n - i; ++j)
        rho += x[j]*x[i+j];

      rho /= (n-1);

      if (i > 10 && rho+rho0 < 0)
        break;

      tint += rho*2;
      rho0 = rho;
    }
  }

  free(x);

  return tint;
}

static void hpd_interval(double * x,
                         long n,
                         double * ltail,
                         double * rtail,
                         double alpha)
{
  long lrow = (long)(n*alpha/2);
  long urow = (long)(n*(1-alpha/2));
  long diffrow = urow - lrow;
  long l,r;

  long left = lrow;
  

  double w = x[urow] - x[lrow];

  *ltail = x[lrow];
  *rtail = x[urow];

  if (n <= 2) return;

  for (l=0,r=l+diffrow; r < n; l++,r++)
  {
    if (x[r] - x[l] < w)
    {
      left = l;
      w = x[r] - x[l];
    }
  }

  *ltail = x[left];
  *rtail = x[left + diffrow];
}

static char * cb_attributes(const snode_t * node)
{
  char * s = NULL;

  /* inner node */
  if (node->left)
  {
    nodepinfo_t * info = (nodepinfo_t *)(node->data);
    if (node->parent)
    {
      nodepinfo_t * pinfo = (nodepinfo_t *)(node->parent->data);

      if (opt_est_theta)
        xasprintf(&s,"[&height_95%%_HPD={%.8f, %.8f}, theta=%.7f]: %f",
                  info->lo, info->hi, info->theta, pinfo->age - info->age);
      else
        xasprintf(&s,"[&height_95%%_HPD={%.8f, %.8f}]: %f",
                  info->lo, info->hi, pinfo->age - info->age);
    }
    else
    {
      if (opt_est_theta)
        xasprintf(&s,"[&height_95%%_HPD={%.8f, %.8f}, theta=%.7f]",
                  info->lo, info->hi, info->theta);
      else
        xasprintf(&s,"[&height_95%%_HPD={%.8f, %.8f}]", info->lo, info->hi);
    }
  }
  else
  {
    nodepinfo_t * info = (nodepinfo_t *)(node->data);
    nodepinfo_t * pinfo = (nodepinfo_t *)(node->parent->data);
    xasprintf(&s,"%s: %f", node->label, pinfo->age - info->age);
  }

  return s;
}

static void write_figtree(stree_t * stree,
                          double * mean,
                          double * hpd025,
                          double * hpd975)
{
  long i;
  FILE * fp_tree = NULL;
  
  fp_tree = xopen("FigTree.tre","w");

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    stree->nodes[i]->data = (void *)xmalloc(sizeof(nodepinfo_t));

  long theta_count = 0;
  if (opt_est_theta)
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    {
      nodepinfo_t * info = (nodepinfo_t *)(stree->nodes[i]->data);
      if (stree->nodes[i]->theta >= 0)
        info->theta = mean[theta_count++];
      else
        info->theta = -1;
    }


  for (i = 0; i < stree->tip_count; ++i)
  {
    nodepinfo_t * info = (nodepinfo_t *)(stree->nodes[i]->data);
    info->lo = info->hi = info->age = 0;
  }

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->tau > 0)
    {
      nodepinfo_t * info = (nodepinfo_t *)(stree->nodes[i]->data);
      info->lo  = hpd025[theta_count + i - stree->tip_count];
      info->hi  = hpd975[theta_count + i - stree->tip_count];
      info->age = mean[theta_count + i - stree->tip_count];
    }
    else
      assert(0);

  char * newick = stree_export_newick(stree->root,cb_attributes);

  fprintf(fp_tree,
          "#NEXUS\n"
          "BEGIN TREES;\n"
          "\tUTREE 1 = %s\n"
          "END;\n\n"
          "[Species tree with tau as branch lengths and theta as labels, for FigTree.\n"
          "In FigTree, choose 95%%HPD for Node Bars and label for Node Labels]\n",
          newick);


  free(newick);

  fclose(fp_tree);

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    free(stree->nodes[i]->data);

  return;
}

void allfixed_summary(FILE * fp_out, stree_t * stree)
{
  long i,j,count;
  long sample_num;
  long rc = 0;
  FILE * fp;

  /* TODO: pretty-fy output */

  fp = xopen(opt_mcmcfile,"r");

  /* skip line containing header */
  getnextline(fp);
  assert(strlen(line) > 4);
  char * header = xstrdup(line);

  /* compute number of columns in the file */
  long col_count = 0;

  /* compute number of theta parameters */
  if (opt_est_theta)
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      if (stree->nodes[i]->theta >= 0)
        col_count++;

  /* compute number of tau parameters */
  for (i = 0; i < stree->inner_count; ++i)
    if (stree->nodes[stree->tip_count+i]->tau)
      col_count++;

  if (opt_est_locusrate && opt_print_locusrate)
    col_count += opt_locus_count;

  if (opt_est_heredity && opt_print_hscalars)
    col_count += opt_locus_count;

  /* add one more for log-L if usedata is on */
  if (opt_usedata)
    col_count++;


  /* allocate storage matrix */
  double ** matrix = (double **)xmalloc((size_t)col_count * sizeof(double *));
  for (i = 0; i < col_count; ++i)
    matrix[i] = (double *)xmalloc((size_t)opt_samples * sizeof(double));

  double * mean = (double *)xmalloc((size_t)col_count * sizeof(double));

  double * hpd025 = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * hpd975 = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * tint = (double *)xmalloc((size_t)col_count * sizeof(double));
  double * stdev = (double *)xmalloc((size_t)col_count * sizeof(double));

  long line_count = 0;
  long bad_count = 0;
  long lineno = 0;
  long prevbad = 0;

  /* read data line by line and store in matrix */
  while (getnextline(fp))
  {
    double x;
    char * p = line;

    ++lineno;

    /* skip sample number */
    count = get_long(p,&sample_num);
    if (!count) goto l_unwind;

    p += count;

    /* read remaining elements of current row */

    for (i = 0; i < col_count; ++i)
    {
      count = get_double(p,&x);
      if (!count)
      {
        if (prevbad || (line_count == 0))
        {
          if (line_count == 0)
            fprintf(stderr,
                    "ERROR: First record has mismatching number of columns\n");
          else
            fprintf(stderr,
                    "ERROR: Found two consecutive records with mismatching "
                    "number of columns (lines %ld and %ld)\n", lineno-1,lineno);
          goto l_unwind;
        }
        else
        {
          fprintf(stderr,
                  "WARNING: Found and ignored record with mismatching number "
                  "of columns (line %ld)\n", lineno);
          prevbad = 1;
          assert(line_count > 0);
          bad_count++;
          break;
        }
      }

      p += count;

      matrix[i][line_count] = x;
    }
    if (i == col_count)
    {
      line_count++;
      prevbad = 0;
    }
  }
  assert(line_count > 0);
  if (bad_count)
    fprintf(stderr, "Skipped a total of %ld erroneous records...\n", bad_count);

  fprintf(stdout, "          %s\n", header+4);
  fprintf(fp_out, "          %s\n", header+4);
  free(header);


  /* compute means */
  fprintf(stdout, "mean    ");
  fprintf(fp_out, "mean    ");
  for (i = 0; i < col_count; ++i)
  {
    double sum = 0;
    for (j = 0; j < opt_samples; ++j)
      sum += matrix[i][j];

    mean[i] = sum/opt_samples;
    fprintf(stdout, "  %f", mean[i]);
    fprintf(fp_out, "  %f", mean[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* compute standard deviation */
  for (i = 0; i < col_count; ++i)
  {
    double sd = 0;
    for (j = 0; j < opt_samples; ++j)
      sd += (matrix[i][j]-mean[i]) * (matrix[i][j]-mean[i]);

    stdev[i] = sqrt(sd/(opt_samples-1));
  }
  
  /* compute tint */
  for (i = 0; i < col_count; ++i)
    tint[i] = eff_ict(matrix[i],opt_samples,mean[i],stdev[i]);

  /* compute and print medians */
  fprintf(stdout, "median  ");
  fprintf(fp_out, "median  ");
  long median_line = opt_samples / 2;

  for (i = 0; i < col_count; ++i)
  {
    qsort(matrix[i], opt_samples, sizeof(double), cb_cmp_double);

    double median = matrix[i][median_line];
    if ((opt_samples & 1) == 0)
    {
      median += matrix[i][median_line-1];
      median /= 2;
    }

    fprintf(stdout, "  %f", median);
    fprintf(fp_out, "  %f", median);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print standard deviation */
  fprintf(stdout, "S.D     ");
  fprintf(fp_out, "S.D     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", stdev[i]);
    fprintf(fp_out, "  %f", stdev[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print minimum values */
  fprintf(stdout, "min     ");
  fprintf(fp_out, "min     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", matrix[i][0]);
    fprintf(fp_out, "  %f", matrix[i][0]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print maximum values */
  fprintf(stdout, "max     ");
  fprintf(fp_out, "max     ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", matrix[i][opt_samples-1]);
    fprintf(fp_out, "  %f", matrix[i][opt_samples-1]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print line at 2.5% of matrix */
  fprintf(stdout, "2.5%%    ");
  fprintf(fp_out, "2.5%%    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", matrix[i][(long)(opt_samples*.025)]);
    fprintf(fp_out, "  %f", matrix[i][(long)(opt_samples*.025)]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print line at 97.5% of matrix */
  fprintf(stdout, "97.5%%   ");
  fprintf(fp_out, "97.5%%   ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", matrix[i][(long)(opt_samples*.975)]);
    fprintf(fp_out, "  %f", matrix[i][(long)(opt_samples*.975)]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* compute and print HPD 2.5% and 97.5% */
  for (i = 0; i < col_count; ++i)
    hpd_interval(matrix[i],opt_samples,hpd025+i,hpd975+i,0.05);

  /* print 2.5% HPD */
  fprintf(stdout, "2.5%%HPD ");
  fprintf(fp_out, "2.5%%HPD ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", hpd025[i]);
    fprintf(fp_out, "  %f", hpd025[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print 97.5% HPD */
  fprintf(stdout, "97.5%%HPD");
  fprintf(fp_out, "97.5%%HPD");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", hpd975[i]);
    fprintf(fp_out, "  %f", hpd975[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print ESS */
  fprintf(stdout, "ESS*    ");
  fprintf(fp_out, "ESS*    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", opt_samples/tint[i]);
    fprintf(fp_out, "  %f", opt_samples/tint[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
    
  /* print Eff */
  fprintf(stdout, "Eff*    ");
  fprintf(fp_out, "Eff*    ");
  for (i = 0; i < col_count; ++i)
  {
    fprintf(stdout, "  %f", 1/tint[i]);
    fprintf(fp_out, "  %f", 1/tint[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* success */
  rc = 1;

l_unwind:
  for (i = 0; i < col_count; ++i)
    free(matrix[i]);
  free(matrix);

  if (rc)
  {
    /* write figtree file */
    write_figtree(stree,mean,hpd025,hpd975);
    fprintf(stdout, "FigTree tree is in FigTree.tre\n");
  }

  free(mean);
  free(hpd025);
  free(hpd975);
  free(tint);
  free(stdev);

  fclose(fp);

  if (!rc)
    fatal("Error while reading/summarizing %s", opt_mcmcfile);

}

