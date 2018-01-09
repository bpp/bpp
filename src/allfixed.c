/*
    Copyright (C) 2016-2017 Tomas Flouri and Ziheng Yang

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

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "bpp.h"

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

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

void allfixed_summary(stree_t * stree)
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
  printf("          %s\n", line+4);

  /* compute number of columns in the file */
  long col_count = 0;

  /* compute number of theta parameters */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->theta >= 0)
      col_count++;

  /* compute number of tau parameters */
  for (i = 0; i < stree->inner_count; ++i)
    if (stree->nodes[stree->tip_count+i]->tau)
      col_count++;

  /* add one more for log-L */
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

  /* read data line by line and store in matrix */
  while (getnextline(fp))
  {
    double x;
    char * p = line;

    /* skip sample number */
    count = get_long(p,&sample_num);
    if (!count) goto l_unwind;

    p += count;

    /* read remaining elements of current row */

    for (i = 0; i < col_count; ++i)
    {
      count = get_double(p,&x);
      if (!count) goto l_unwind;

      p += count;

      matrix[i][line_count] = x;
    }

    line_count++;
  }

  /* compute means */
  printf("mean    ");
  for (i = 0; i < col_count; ++i)
  {
    double sum = 0;
    for (j = 0; j < opt_samples; ++j)
      sum += matrix[i][j];

    mean[i] = sum/opt_samples;
    printf("  %f", mean[i]);
  }
  printf("\n");

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
  printf("median  ");
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

    printf("  %f", median);
  }
  printf("\n");

  /* print standard deviation */
  printf("S.D     ");
  for (i = 0; i < col_count; ++i)
  {
    printf("  %f", stdev[i]);
  }
  printf("\n");

  /* print minimum values */
  printf("min     ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][0]);
  printf("\n");

  /* print maximum values */
  printf("max     ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][opt_samples-1]);
  printf("\n");

  /* print line at 2.5% of matrix */
  printf("2.5%%    ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][(long)(opt_samples*.025)]);
  printf("\n");

  /* print line at 97.5% of matrix */
  printf("97.5%%   ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", matrix[i][(long)(opt_samples*.975)]);
  printf("\n");

  /* compute and print HPD 2.5% and 97.5% */
  for (i = 0; i < col_count; ++i)
    hpd_interval(matrix[i],opt_samples,hpd025+i,hpd975+i,0.05);

  /* print 2.5% HPD */
  printf("2.5%%HPD ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", hpd025[i]);
  printf("\n");

  /* print 97.5% HPD */
  printf("97.5%%HPD");
  for (i = 0; i < col_count; ++i)
    printf("  %f", hpd975[i]);
  printf("\n");

  /* print ESS */
  printf("ESS*    ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", opt_samples/tint[i]);
  printf("\n");
    
  /* print Eff */
  printf("Eff*    ");
  for (i = 0; i < col_count; ++i)
    printf("  %f", 1/tint[i]);
  printf("\n");

  /* success */
  rc = 1;

l_unwind:
  for (i = 0; i < col_count; ++i)
    free(matrix[i]);
  free(matrix);

  free(mean);
  free(hpd025);
  free(hpd975);
  free(tint);
  free(stdev);

  fclose(fp);

  if (!rc)
    fatal("Error while reading/summarizing %s", opt_mcmcfile);

  return;
}

