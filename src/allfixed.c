/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

  ret = ws + end - start;
  free(s);
  return ret;
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

  ret = ws + end - start;
  free(s);
  return ret;
}

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  if ( *x > *y) return 1;
  if ( *x < *y) return -1;
  return 0;
}

#if 1
static double eff_ict(double * y, long n, double mean, double stdev)
{
  /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
     sequence method */

  long i,j;
  double tint = 1;
  double rho, rho0 = 0;
  long maxlag = 2000;
  long minNr = 10;

  double * x = (double *)xmalloc((size_t)n * sizeof(double));
  for (i = 0; i < n; ++i)
    x[i] = (y[i]-mean)/stdev;

  if (stdev/(fabs(mean)+1) < 1E-9)
  {
   tint = n;
  }
  else
  {
    for (i = 1; i < MIN(maxlag,n-minNr); ++i)
    {
      rho = 0;
      for (j = 0; j < n - i; ++j)
        rho += x[j]*x[i+j];

      rho /= (n-i);

      if (i > minNr && rho+rho0 < 0)
        break;

      tint += rho*2;
      rho0 = rho;
    }
  }

  free(x);

  return tint;
}
#else
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
#endif

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
        xasprintf(&s,"%s[&height_95%%_HPD={%.8f, %.8f}, theta=%.7f]: %f", node->label, info->lo, info->hi, info->theta, pinfo->age - info->age);
      else
        xasprintf(&s,"%s[&height_95%%_HPD={%.8f, %.8f}]: %f", node->label, info->lo, info->hi, pinfo->age - info->age);
    }
    else
    {
      if (opt_est_theta)
        xasprintf(&s,"%s[&height_95%%_HPD={%.8f, %.8f}, theta=%.7f]", node->label, info->lo, info->hi, info->theta);
      else
        xasprintf(&s,"%s[&height_95%%_HPD={%.8f, %.8f}]", node->label, info->lo, info->hi);
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

static char * cb_msci_thetasandtaus(const snode_t * node)
{
  char * s = NULL;
  char * stheta = NULL;

  if (node->theta > 0)
    xasprintf(&stheta, "#%f", node->theta);
  else
    stheta = xstrdup("");

  if (!node->left && !node->right)
  {
    /* tip or hybrid mirror node */
    if (node->hybrid)
    {
      /* hybrid mirror node */
      assert(node_is_mirror(node));

      if (node_is_hybridization(node))
      {
        /* hybridization event */
        if (node->hphi >= 0)
          xasprintf(&s,
                    "%s[&phi=%f,tau-parent=%s]:%f %s",
                    node->label,
                    node->hphi,
                    node->htau ? "yes" : "no",
                    node->tau,
                    stheta);
        else
          xasprintf(&s,
                    "%s[tau-parent=%s]:%f %s",
                    node->label,
                    node->htau ? "yes" : "no",
                    node->tau,
                    stheta);
      }
      else
      {
        /* bidirectional introgression */
        if (node->hphi >= 0)
          xasprintf(&s,
                    "%s[&phi=%f]:%f %s",
                    node->label,
                    node->hphi,
                    node->tau,
                    stheta);
        else
          xasprintf(&s, "%s", node->label);
      }
    }
    else
      xasprintf(&s, "%s %s", node->label, stheta);
  }
  else
  {
    /* inner or hybrid node */

    if (node->hybrid)
    {
      /* hybrid non-mirror node */
      assert(!node_is_mirror(node));
      if (node_is_hybridization(node))
      {
        /* hybridization event */
        assert(node->left && !node->right);

        if (node->hphi >= 0)
            xasprintf(&s,
                      "%s[&phi=%f,tau-parent=%s]:%f %s",
                      node->label ? node->label : "",
                      node->hphi,
                      node->htau ? "yes" : "no",
                      node->tau,
                      stheta);
          else
            xasprintf(&s,
                      "%s[tau-parent=%s]:%f %s",
                      node->label ? node->label : "",
                      node->htau ? "yes" : "no",
                      node->tau,
                      stheta);
      }
      else
      {
        /* bidirectional introgression */
        assert(node->left && node->right);
        assert(node->right->hybrid && !node->right->left && !node->right->right);
        if (node->hphi >= 0)
        {
          xasprintf(&s,
                    "%s[&phi=%f]:%f %s",
                    node->label,
                    node->hphi,
                    node->tau,
                    stheta);
        }
        else
        {
            xasprintf(&s,
                      "%s:%f %s",
                      node->label,
                      node->tau,
                      stheta);
        }
      }
    }
    else
    {
      xasprintf(&s,
                "%s:%f %s",
                node->label ? node->label : "",
                node->tau,
                stheta);
    }
  }
  if (stheta)
    free(stheta);
  return s;
}

/* Ziheng 2020-10-2 
*  (i) Space for stree->nodes[i]->data is allocated for for snodes_total, including mirror nodes.
*      Where is the space freed?
*  (ii) I copied mean values for tau and theta into stree->nodes, so the old values are overwritten.
*       Since this is after the mcmc sample is summarized, it is at the end of the run, 
*       so stree probably does not have useful tau and theta values.
*/

long opt_msci_faketree_binarize;

static void write_figtree(FILE * fp_out,
                          stree_t * stree,
                          double * mean,
                          double * hpd025,
                          double * hpd975)
{
  long i,j, theta_count = 0, tau_count = 0;
  int maxlen = 0;
  FILE * fp_tree = NULL;
  unsigned int snodes_total = stree->tip_count + stree->inner_count + stree->hybrid_count;

  fp_tree = xopen(opt_msci ? "FakeTree.tre" : "FigTree.tre", "w");

  for (i = 0; i < snodes_total; ++i)
    stree->nodes[i]->data = (void *)xmalloc(sizeof(nodepinfo_t));

  /* copy theta */
  if (opt_est_theta)
    for (i = 0; i < snodes_total; ++i)
    {
      nodepinfo_t* info = (nodepinfo_t*)(stree->nodes[i]->data);
      if (stree->nodes[i]->theta >= 0)
      {
        info->theta = mean[theta_count++];
        /* Ziheng 2020-10-2 */
        stree->nodes[i]->theta = info->theta;
      }
      else
        info->theta = -1;
    }
  for (i = 0; i < stree->tip_count; ++i)
  {
    nodepinfo_t * info = (nodepinfo_t *)(stree->nodes[i]->data);
    info->lo = info->hi = info->age = 0;
  }

  /* copy tau */
  for (i = stree->tip_count; i < snodes_total; ++i)
    if (stree->nodes[i]->tau > 0)
    {
      if (opt_msci && i >= stree->tip_count + stree->inner_count)  /* hybrid mirror node */
      {
        stree->nodes[i]->tau = stree->nodes[i]->hybrid->tau;
        nodepinfo_t* info = (nodepinfo_t*)(stree->nodes[i]->data);
        nodepinfo_t* tmp = (nodepinfo_t*)(stree->nodes[i]->hybrid->data);
        info->age = tmp->age;
      }
      else 
      {
        nodepinfo_t* info = (nodepinfo_t*)(stree->nodes[i]->data);
        info->lo = hpd025[theta_count - stree->tip_count + i];
        info->hi = hpd975[theta_count - stree->tip_count + i];
        info->age = mean[theta_count - stree->tip_count + i];
        stree->nodes[i]->tau = info->age;
        tau_count++;
      }
    }
    else
      assert(0);

  /* copy phi for msci model */
  for (i = 0; i < stree->hybrid_count; i++)
  {
    snode_t * mnode = stree->nodes[stree->tip_count + stree->inner_count + i];
    
    #if 0
    /* old code before the introduction of has_phi */
    if (!node_is_bidirection(mnode))
    {
      if (mnode->hybrid->htau == 0 && mnode->htau == 1)
        mnode = mnode->hybrid;
    }

    mnode->hphi = mean[theta_count + tau_count + i];
    mnode->hybrid->hphi = 1 - mnode->hphi;
    #else

    /* new correct code */
    snode_t * pnode = mnode->has_phi ? mnode : mnode->hybrid;
    pnode->hphi = mean[theta_count + tau_count + i];
    pnode->hybrid->hphi = 1 - pnode->hphi;
    #endif
  }

  /*** Ziheng 2020-10-2 ***/
  for (maxlen = 0, i = 0; i < snodes_total; ++i)
  {
    if (maxlen < (int)strlen(stree->nodes[i]->label))
      maxlen = strlen(stree->nodes[i]->label);
  }
  fprintf(fp_out, "List of nodes, taus and thetas:\n");
  fprintf(fp_out, "Node (+1)       Tau      Theta    Label\n");
  for (i = 0; i < snodes_total; ++i) {
    fprintf(fp_out,
            "%-9ld %9.6f  %9.6f    %*s ",
            i,
            stree->nodes[i]->tau,
            stree->nodes[i]->theta,
            maxlen,
            stree->nodes[i]->label ? stree->nodes[i]->label : "-");
    fprintf(fp_out, "[");
    for (j = 0; j < stree->tip_count; ++j)
      if (stree->pptable[j][i])
        fprintf(fp_out, " %s", stree->nodes[j]->label);
    fprintf(fp_out, " ]\n");
  }
  char * newick;
  if(!opt_msci)
    newick = stree_export_newick(stree->root, cb_attributes);
  else
  {
    opt_msci_faketree_binarize = 0;
    fprintf(stdout, "\nSpecies tree network:\n");
    fprintf(fp_out, "\nSpecies tree network:\n");
    newick = msci_export_newick(stree->root, NULL);
    fprintf(stdout, "%s\n", newick);
    fprintf(fp_out, "%s\n", newick);
    free(newick);

    /* print network with attributes */
    fprintf(stdout,
            "\nSpecies tree network with attributes (thetas and tau 95%% HPD), "
           "and branch lengths (tau difference):\n");
    fprintf(fp_out,
            "\nSpecies tree network with attributes (thetas and tau 95%% HPD), "
           "and branch lengths (tau difference):\n");
    newick = msci_export_newick(stree->root, cb_attributes);
    fprintf(stdout, "%s\n", newick);
    fprintf(fp_out, "%s\n", newick);
    free(newick);

    /* print network suitable for running simulations */
    fprintf(stdout, "\nSpecies tree network with taus and thetas:\n");
    fprintf(fp_out, "\nSpecies tree network with taus and thetas:\n");
    newick = msci_export_newick(stree->root, cb_msci_thetasandtaus);
    fprintf(stdout, "%s\n", newick);
    fprintf(fp_out, "%s\n", newick);
    free(newick);

    opt_msci_faketree_binarize = 1;
    newick = msci_export_newick(stree->root, cb_attributes);
    opt_msci_faketree_binarize = 0;
  }

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

  for (i = 0; i < snodes_total; ++i)
    free(stree->nodes[i]->data);

  return;
}


void allfixed_summary(FILE * fp_out, stree_t * stree)
{
  long i, j, count;
  long sample_num;
  long rc = 0;
  FILE * fp;
  unsigned int snodes_total = stree->tip_count + stree->inner_count;
  
  if (opt_msci)
    snodes_total += stree->hybrid_count;

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
    for (i = 0; i < snodes_total; ++i)
      if (stree->nodes[i]->theta >= 0)
        col_count++;

  /* compute number of tau parameters */
  for (i = 0; i < stree->inner_count; ++i)
    if (stree->nodes[stree->tip_count+i]->tau)
      col_count++;

  if (opt_migration && !opt_est_geneflow)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        if (opt_mig_bitmatrix[i][j])
          col_count++;
  }

  /* compute number of phi parameters */
  if (opt_msci)
    col_count += stree->hybrid_count;

  /* column for mubar */
  if (opt_est_mubar && opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    ++col_count;

  /* add one more for log-L if usedata is on */
  if (opt_usedata)
    col_count++;

  /* column for vbar (or nu for Gamma-Dirichlet) */
  if (opt_clock != BPP_CLOCK_GLOBAL)
    ++col_count;


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
                    "ERROR: First record has mismatching number of columns (expected %ld)\n",col_count);
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

  /* BDI label-switching routine to resolve unidentifiability issues */
  if (opt_msci && opt_usedata)
  {
    long bidir_count = 0;
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      if (stree->nodes[i]->hybrid &&
          node_is_bidirection(stree->nodes[i]) &&
          stree->nodes[i]->prop_tau)
      {
        bidir_count++;
      }

    if (bidir_count && opt_linkedtheta != BPP_LINKEDTHETA_MSCI)
      lswitch(stree, header, matrix, col_count, fp_out);
  }

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

  if (rc && stree->tip_count > 1)
  {
    /* write figtree file */
    write_figtree(fp_out, stree, mean, hpd025, hpd975);
    if (!opt_msci)
      fprintf(stdout, "\nFigTree tree is in FigTree.tre\n");
    else 
      fprintf(stdout, "\nFigTree tree is in FakeTree.tre\n");
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
