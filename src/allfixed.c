/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#define DIST_GAMMA      0
#define DIST_INVGAMMA   1
#define DIST_BETA       2

#define A1B1_BINS 1000
#define A1B1_BINS_M 60
#define A1B1_TAIL 0.05

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

  ret = ws + end - start;
  free(s);
  return ret;
}

static long get_doubleordash(const char * line, double * value)
{
  int len=0;
  long ret = 0;
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

  if (!strcmp(start,"-"))
  {
    ret = ws+end-start;
    free(s);
    return ret;
  }

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
double eff_ict(double * y, long n, double mean, double stdev, double * rho1)
{
  /* This calculates Efficiency or Tint using Geyer's (1992) initial positive
     sequence method */

  long i,j;
  double tint = -1;
  double rho, rho0 = 1;
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

      if (i == 1)
        *rho1 = rho;

      if (i % 2 == 1)
      {
        if (i > minNr && rho+rho0 < 0)
          break;
        else
          tint += (rho0 + rho) * 2;
      }
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
  size_t i;
  char * tmplabel = NULL;
  char * s = NULL;
  double theta;

  /* inner node */
  if (node->left)
  {
    nodepinfo_t * info = (nodepinfo_t *)(node->data);
    if (node->linked_theta)
    {
      nodepinfo_t * linfo = (nodepinfo_t *)(node->linked_theta->data);
      theta = linfo->theta;
    }
    else
      theta = info->theta;
      
    /* replace commas with pipes in label */
    tmplabel = xstrdup(node->label);
    for (i = 0; i < strlen(tmplabel); ++i)
      if (tmplabel[i] == ',')
        tmplabel[i] = '|';

    if (node->parent)
    {
      nodepinfo_t * pinfo = (nodepinfo_t *)(node->parent->data);
      if (opt_est_theta)
        xasprintf(&s,
                  "%s[&height_95%%_HPD={%.8f, %.8f}, theta=%.7f]: %f",
                  tmplabel, info->lo, info->hi, theta, pinfo->age-info->age);
      else
        xasprintf(&s,
                  "%s[&height_95%%_HPD={%.8f, %.8f}]: %f",
                  tmplabel, info->lo, info->hi, pinfo->age-info->age);
    }
    else
    {
      if (opt_est_theta)
        xasprintf(&s,
                  "%s[&height_95%%_HPD={%.8f, %.8f}, theta=%.7f]",
                  tmplabel, info->lo, info->hi, theta);
      else
        xasprintf(&s,
                  "%s[&height_95%%_HPD={%.8f, %.8f}]",
                  tmplabel, info->lo, info->hi);
    }
    free(tmplabel);
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
  char * tmplabel = NULL;
  size_t i;

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
      
      tmplabel = NULL;
      if (node->label)
      {
        tmplabel = xstrdup(node->label);
        for (i = 0; i < strlen(tmplabel); ++i)
          if (tmplabel[i] == ',')
            tmplabel[i] = '|';
      }

      xasprintf(&s,
                "%s:%f %s",
                tmplabel ? tmplabel : "",
                node->tau,
                stheta);

      if (tmplabel)
        free(tmplabel);
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

  char * treefile = NULL;
  xasprintf(&treefile,
            opt_msci ? "%s.FakeTree.tre" : "%s.FigTree.tre",
            opt_jobname);

  fp_tree = xopen(treefile, "w");
  free(treefile);

  for (i = 0; i < snodes_total; ++i)
    stree->nodes[i]->data = (void *)xcalloc(sizeof(nodepinfo_t),1);

  /* copy theta */
  if (opt_est_theta)
  {
    for (i = 0; i < snodes_total; ++i)
    {
      nodepinfo_t* info = (nodepinfo_t*)(stree->nodes[i]->data);
      if (stree->nodes[i]->theta >= 0 && stree->nodes[i]->linked_theta == NULL)
      {
        info->theta = mean[theta_count++];
        stree->nodes[i]->theta = info->theta;
      }
    }
    for (i = 0; i < snodes_total; ++i)
      if (stree->nodes[i]->theta >= 0 && stree->nodes[i]->linked_theta)
        stree->nodes[i]->theta = stree->nodes[i]->linked_theta->theta;
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

  //for (i = 0; i < snodes_total; ++i)
  //  free(stree->nodes[i]->data);

  return;
}

static long tokens_count(const char * hdr)
{
  size_t tokenlen;
  char * p = xstrdup(hdr);
  char * s = p;
  long tokens = 0;

  while (*s)
  {
    /* skip spaces */
    tokenlen = strspn(s," \t\r\n");
    s += tokenlen;

    tokenlen = strcspn(s," \t\r\n");
    if (tokenlen)
      ++tokens;
    s += tokenlen;
  }

  free(p);
  return tokens;
}

static char ** header_explode(const char * hdr, long * count)
{
  long tcount = 0;
  size_t tokenlen;
  char * p = xstrdup(hdr);
  char * s = p;
  char ** tokens;

  tcount = tokens_count(hdr);
  assert(tcount>1);
  --tcount;

  tokens = (char **)xmalloc((size_t)tcount*sizeof(char *));

  assert(hdr[0] == 'G');
  assert(hdr[1] == 'e');
  assert(hdr[2] == 'n');
  assert(hdr[3] == '\t');

  s += 3;

  tcount = 0;
  while (*s)
  {
    /* skip spaces */
    tokenlen = strspn(s," \t\r\n");
    s += tokenlen;

    tokenlen = strcspn(s," \t\r\n");

    tokens[tcount++] = xstrndup(s,tokenlen);

    s += tokenlen;
  }

  free(p);

  *count = tcount;
  return tokens;
}

static void header_tokens_shorten(char ** tokens, long keep_labels, long count)
{
  long i;

  if (keep_labels)
  {
    /* keep original labels instead of mapped numbers */
    for (i = 0; i < count; ++i)
    {
      char * start = tokens[i];
      char * p = strchr(start,':');
      if (p)
      {
        char * q = strchr(p+1,':');
        assert(q);
        assert(*(q+1));
        assert(p != q);
        long len = strlen(q+1);
        memmove(p+1,q+1,len);
        p[len+1] = 0;
        char * tmp = xstrdup(start);
        free(tokens[i]);
        tokens[i] = tmp;
      }
    }
    return;
  }

  for (i = 0; i < count; ++i)
  {
    char * start = tokens[i];
    char * p = strchr(start,':');
    if (p)
    {
      p = strchr(p+1,':');
      if (p)
      {
        char * tmp = xstrndup(tokens[i],p-start);
        free(tokens[i]);
        tokens[i] = tmp;
      }
    }
  }
}

static double xfloor1(double x)
{
  double r = floor(x);
  if (r < 1) r = 1;
  return r;
}

static char * center(const char * s, int space)
{
  char * r;
  int len = (int)strlen(s);
  int left,right;

  left  = (space-len)/2;
  right = space-len-left;

  xasprintf(&r, "%*s%s%*s", left, "", s, right, "");

  return r;
}

static long find_theta_index(char ** tokens, long cols, const char * wstr)
{
  long i,k,m;

  const char * theta_prefix = "theta:";
  const char * tok;
  char * wtmp;
  char * p;
  size_t tp_len = strlen(theta_prefix);
  size_t toklen;

  /* get the wstr number into m */
  wtmp = xstrdup(wstr);
  p = wtmp + strlen(wtmp);
  assert(*p == 0);
  --p;
  assert(*p == '1');
  --p;
  assert(*p == 'a' || *p == 'b');
  --p;
  assert(*p == '_');

  while (*(--p) != '>')
  {
    assert(*p >= '0' && *p <= '9');
  }
  p++;
  p = xstrdup(p);
  m = atol(p);

  free(wtmp);
  free(p);


  for (i = 0; i < cols; ++i)
  {
    /* get token and length */
    tok = tokens[i];
    toklen = strlen(tok);

    /* check if its a theta */
    if (!(toklen > tp_len && !strncmp(tok, theta_prefix, tp_len)))
      continue;

    /* skip b1's */
    assert(tok[toklen-1] == '1');
    assert(tok[toklen-3] == '_');
    if (tok[toklen-2] == 'b')
      continue;
    else
      assert(tok[toklen-2] == 'a');

    /* get the theta string number and convert it to long */
    char * x = xstrndup(tok+tp_len, toklen-tp_len);
    k = atol(x);
    free(x);

    if (k == m) break;
  }
  if (i == cols)
    return -1;

  return i;
}
static void update_col_sizes(int * label_size,
                             char * pname,
                             int * pname_size,
                             int prec,
                             double mean,
                             double stdev,
                             double et025,
                             double et975,
                             double hpd025,
                             double hpd975,
                             double c,
                             double effu,
                             double effy)
{
  int digits;

  /* mean */
  digits = (int)floor(log10(xfloor1(mean))+1);
  digits += prec+1;
  label_size[0] = MAX(label_size[0],digits);

  /* standard deviation */
  digits = (int)floor(log10(xfloor1(stdev))+1);
  digits += prec+1;
  label_size[1] = MAX(label_size[1],digits);

  /* 2.5% */
  digits = (int)floor(log10(xfloor1(et025))+1);
  digits += prec+1;
  label_size[2] = MAX(label_size[2],digits);

  /* 97.5% */
  digits = (int)floor(log10(xfloor1(et975))+1);
  digits += prec+1;
  label_size[3] = MAX(label_size[3],digits);

  /* 2.5% HPD */
  digits = (int)floor(log10(xfloor1(hpd025))+1);
  digits += prec+1;
  label_size[4] = MAX(label_size[4],digits);

  /* 97.5% HPD */
  digits = (int)floor(log10(xfloor1(hpd975))+1);
  digits += prec+1;
  label_size[5] = MAX(label_size[5],digits);

  /* Effu */
  digits = (int)floor(log10(xfloor1(effu))+1);
  digits += prec+1;
  label_size[6] = MAX(label_size[6],digits);

  /* Effy */
  digits = (int)floor(log10(xfloor1(effy))+1);
  digits += prec+1;
  label_size[7] = MAX(label_size[7],digits);

  /* c */
  digits = (int)floor(log10(xfloor1(c))+1);
  digits += prec+1;
  label_size[8] = MAX(label_size[8],digits);

  /* get column name */
  *pname_size = MAX(*pname_size,strlen(pname));
}

static void print_a1b1_summary_line(FILE * fp_out,
                                    int prec,
                                    int pname_size,
                                    char * pname,
                                    int * label_size,
                                    double mean,
                                    double stdev,
                                    double et025,
                                    double et975,
                                    double hpd025,
                                    double hpd975,
                                    double c,
                                    double effu,
                                    double effy)
{
  char * s;

  /* print line */
  fprintf(stdout, "%-*s", pname_size,pname);
  fprintf(fp_out, "%-*s", pname_size,pname);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* mean */
  xasprintf(&s, "%.*f", prec, mean);
  fprintf(stdout, "%*s", label_size[0], s);
  fprintf(fp_out, "%*s", label_size[0], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* S.D */
  xasprintf(&s, "%.*f", prec, stdev);
  fprintf(stdout, "%*s", label_size[1], s);
  fprintf(fp_out, "%*s", label_size[1], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* 2.5% */
  xasprintf(&s, "%.*f", prec, et025);
  fprintf(stdout, "%*s", label_size[2], s);
  fprintf(fp_out, "%*s", label_size[2], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* 97.5% */
  xasprintf(&s, "%.*f", prec, et975);
  fprintf(stdout, "%*s", label_size[3], s);
  fprintf(fp_out, "%*s", label_size[3], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* 2.5% HPD */
  xasprintf(&s, "%.*f", prec, hpd025);
  fprintf(stdout, "%*s", label_size[4], s);
  fprintf(fp_out, "%*s", label_size[4], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* 97.5% HPD */
  xasprintf(&s, "%.*f", prec, hpd975);
  fprintf(stdout, "%*s", label_size[5], s);
  fprintf(fp_out, "%*s", label_size[5], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* Effu */
  xasprintf(&s, "%.*f", prec, effu);
  fprintf(stdout, "%*s", label_size[6], s);
  fprintf(fp_out, "%*s", label_size[6], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* Effy */
  xasprintf(&s, "%.*f", prec, effy);
  fprintf(stdout, "%*s", label_size[7], s);
  fprintf(fp_out, "%*s", label_size[7], s);
  free(s);

  fprintf(stdout, "  ");
  fprintf(fp_out, "  ");

  /* c */
  xasprintf(&s, "%.*f", prec, c);
  fprintf(stdout, "%*s", label_size[8], s);
  fprintf(fp_out, "%*s", label_size[8], s);
  free(s);

  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

}

static void summarize_a1b1(stree_t * stree, FILE * fp_out)
{
  int prec;
  int * label_size = NULL;
  unsigned int snodes_total;
  long i,j,k = 0;
  long count;
  long cols = 0;
  long line_count = 0;
  long bad_count = 0;
  long na_count = 0;
  long lineno = 0;
  long prevbad = 0;
  long sample_num;
  long param_count = 0;
  long * a1_index = NULL;
  long * b1_index = NULL;
  long * a2_index = NULL;
  long * b2_index = NULL;
  long * prdist1 = NULL;
  long * prdist2 = NULL;
  FILE * fp;
  char ** pname_mcmc = NULL;
  char ** pname = NULL;

  char * label[] =
  {
    "mean",       /* 0 */
    "S.D",        /* 1 */
    "2.5%",       /* 2 */
    "97.5%",      /* 3 */
    "2.5%HPD",    /* 4 */
    "97.5%HPD",   /* 5 */
    "Effu",       /* 6 */
    "Effy",       /* 7 */
    "c"           /* 8 */
  };
  snodes_total = stree->tip_count + stree->inner_count + stree->hybrid_count;
  fp = xopen(opt_a1b1file,"r");

  /* skip line containing header */
  getnextline(fp);
  assert(strlen(line) > 4);
  char * header = xstrdup(line);

  long token_count;
  char ** tokens = header_explode(header,&token_count);

  char ** numeric_tokens = (char **)xmalloc((size_t)token_count*sizeof(char *));
  for (i = 0; i < token_count; ++i)
    numeric_tokens[i] = xstrdup(tokens[i]);
  header_tokens_shorten(tokens,opt_keep_labels,token_count);
  /* the following is required for finding theta indices associated with Ws */
  header_tokens_shorten(numeric_tokens,0,token_count);

  /* compute number of theta parameters */
  for (i = 0; i < snodes_total; ++i)
    if (stree->nodes[i]->theta >= 0 && stree->nodes[i]->linked_theta == NULL)
      cols += 2;

#if  1
  if (opt_migration && !opt_est_geneflow)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        if (opt_mig_bitmatrix[i][j])
          cols += 2;
  }
#endif
  /* compute number of phi parameters */
  if (opt_msci)
    cols += 2*stree->hybrid_count;

  /* number of columns including converted Ms (from Ws) */
  /* allocate storage matrix */
  double ** matrix = (double **)xmalloc((size_t)cols * sizeof(double *));
  for (i = 0; i < cols; ++i)
    matrix[i] = (double *)xmalloc((size_t)opt_samples * sizeof(double));

  /* at most cols+cols/2 parameters to process (cols/2 is max # of Ms) */
  double * wmean = (double *)xmalloc((size_t)cols * sizeof(double));
  double * wstdev= (double *)xmalloc((size_t)cols * sizeof(double));
  double * mean = (double *)xmalloc((size_t)cols * sizeof(double));
  double * et025 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * et975 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * hpd025 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * hpd975 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * tint = (double *)xmalloc((size_t)cols * sizeof(double));
  double * rho1 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * c = (double *)xmalloc((size_t)cols * sizeof(double));
  double * stdev = (double *)xmalloc((size_t)cols * sizeof(double));
  double * effu = (double *)xmalloc((size_t)cols * sizeof(double));
  double * effy = (double *)xmalloc((size_t)cols * sizeof(double));

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

    for (i = 0; i < cols; ++i)
    {
      x = -1;
      count = get_doubleordash(p,&x);
      if (!count)
      {
        if (prevbad || (line_count == 0))
        {
          printf("i = %ld\n", i);
          if (line_count == 0)
            fprintf(stderr,
                    "ERROR: First record has mismatching number of columns (expected %ld)\n",cols);
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

      /* skip if missing data */
      if (x == -1)
        ++na_count;

      matrix[i][line_count] = x;
    }
    if (i == cols)
    {
      line_count++;
      prevbad = 0;
    }
  }
  assert(line_count > 0);
  if (bad_count)
    fprintf(stderr, "Skipped a total of %ld erroneous records...\n", bad_count);
  if (na_count)
    fprintf(stderr, "Found %ld missing data entries...\n", na_count);

  long label_count = sizeof(label) / sizeof(label[0]);
  label_size = (int *)xcalloc((size_t)label_count,sizeof(int));
  int pname_size = 5;  /* param */
  prec = 6;

  for (i = 0; i < label_count; ++i)
    label_size[i] = MAX(prec, (int)strlen(label[i]));

  a1_index = (long *)xmalloc((size_t)cols*sizeof(long));
  b1_index = (long *)xmalloc((size_t)cols*sizeof(long));
  a2_index = (long *)xmalloc((size_t)cols*sizeof(long));
  b2_index = (long *)xmalloc((size_t)cols*sizeof(long));
  prdist1  = (long *)xmalloc((size_t)cols*sizeof(long));
  prdist2  = (long *)xmalloc((size_t)cols*sizeof(long));

  for (i = 0; i < cols; ++i)
    a2_index[i] = b2_index[i] = -1;
  
  const char * theta_prefix = "theta:";
  const char * phi_prefix = "phi:";
  const char * w_prefix = "W:";
  size_t tp_len = strlen(theta_prefix);
  size_t pp_len = strlen(phi_prefix);
  size_t wp_len = strlen(w_prefix);

  pname      = (char **)xmalloc((size_t)cols*sizeof(char *));
  pname_mcmc = (char **)xmalloc((size_t)cols*sizeof(char *));
  long dist_theta = DIST_INVGAMMA;
  if (opt_theta_prior == BPP_THETA_PRIOR_GAMMA)
  {
    if (opt_theta_prop == BPP_THETA_PROP_MG_GAMMA)
      dist_theta = DIST_GAMMA;
    else
      assert(opt_theta_prop == BPP_THETA_PROP_MG_INVG);
  }
  else
    assert(opt_theta_prior == BPP_THETA_PRIOR_INVGAMMA);

  long dist_phi = DIST_BETA;
  long dist_w = DIST_GAMMA;

  for (i = 0; i < cols; ++i)
  {
    char * tok = tokens[i];
    long trim_suffix = 0;
    size_t toklen = strlen(tok);

    if (toklen > tp_len && !strncmp(tok, theta_prefix, tp_len))
    {
      /* theta */
      a1_index[param_count] = i;
      b1_index[param_count] = ++i;
      prdist1[param_count] = dist_theta;
    }
    else if (toklen > pp_len && !strncmp(tok, phi_prefix, pp_len))
    {
      /* phi */
      a1_index[param_count] = i;
      b1_index[param_count] = ++i;
      prdist1[param_count] = dist_phi;
    }
    else if (toklen > wp_len && !strncmp(tok, w_prefix, wp_len))
    {
      /* W */
      a1_index[param_count] = i;
      b1_index[param_count] = ++i;

      a2_index[param_count] = find_theta_index(numeric_tokens, cols, numeric_tokens[i]);
      b2_index[param_count] = a2_index[param_count]+1;

      prdist1[param_count] = dist_w;
      prdist2[param_count] = dist_theta;
    }
    else
    {
      assert(0);
    }

    if (!trim_suffix)
    {
      assert(toklen > 3);
      assert(tok[toklen-1] == '1');
      assert(tok[toklen-2] == 'a');
      assert(tok[toklen-3] == '_');
      pname_mcmc[param_count++] = xstrndup(tok,toklen-3);
    }
    else
    {
      pname_mcmc[param_count++] = xstrdup(tok);
    }
  }

  k=0;
  assert(param_count);
  for (i = 0; i < param_count; ++i)
  {
    double * a1_list = matrix[a1_index[i]];
    double * b1_list = matrix[b1_index[i]];
    double * a2_list = NULL;
    double * b2_list = NULL;

    conditional_to_marginal(a1_list,
                            b1_list,
                            opt_samples,
                            A1B1_BINS,
                            prdist1[i],
                            A1B1_TAIL,
                            mean+k,
                            stdev+k,
                            et025+k,
                            et975+k,
                            hpd025+k,
                            hpd975+k,
                            c+k,
                            effu+k,
                            effy+k);
    pname[k] = xstrdup(pname_mcmc[i]); 
    
    update_col_sizes(label_size,
                     pname[k],
                     &pname_size,
                     prec,
                     mean[k],
                     stdev[k],
                     et025[k],
                     et975[k],
                     hpd025[k],
                     hpd975[k],
                     c[k],
                     effu[k],
                     effy[k]);
    ++k;



    if (a2_index[i] >= 0)
    {
      a2_list = matrix[a2_index[i]];
      b2_list = matrix[b2_index[i]];

      conditional_to_marginal_M(a2_list,
                                b2_list,
                                a1_list,
                                b1_list,
                                opt_samples,
                                A1B1_BINS_M,
                                prdist2[i],
                                A1B1_TAIL,
                                wmean+k,
                                wstdev+k,
                                mean+k,
                                stdev+k,
                                et025+k,
                                et975+k,
                                hpd025+k,
                                hpd975+k,
                                c+k,
                                effu+k,
                                effy+k);
      pname[k] = xstrdup(pname_mcmc[i]); 
      assert(pname[k][0] == 'W');
      pname[k][0] = 'M';
      
      update_col_sizes(label_size,
                       pname[k],
                       &pname_size,
                       prec,
                       mean[k],
                       stdev[k],
                       et025[k],
                       et975[k],
                       hpd025[k],
                       hpd975[k],
                       c[k],
                       effu[k],
                       effy[k]);
      ++k;
    }
  }

  /* print centered header */
  int linesize = pname_size;
  char * s = center("param", pname_size);
  fprintf(stdout, "%s",s);
  fprintf(fp_out, "%s",s);
  free(s);

  for (i = 0; i < label_count; ++i)
  {
    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");
    s = center(label[i],label_size[i]);
    fprintf(stdout, "%s", s);
    fprintf(fp_out, "%s", s);
    free(s);

    linesize += 2 + label_size[i];
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
  for (i = 0; i < linesize; ++i)
  {
    fprintf(stdout, "-");
    fprintf(fp_out, "-");
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  for (i = 0; i < k; ++i)
    print_a1b1_summary_line(fp_out,
                            prec,
                            pname_size,
                            pname[i],
                            label_size,
                            mean[i],
                            stdev[i],
                            et025[i],
                            et975[i],
                            hpd025[i],
                            hpd975[i],
                            c[i],
                            effu[i],
                            effy[i]);
    

  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

l_unwind:
  free(label_size);
  free(a1_index);
  free(b1_index);
  free(a2_index);
  free(b2_index);
  free(header);
  free(prdist1);
  free(prdist2);

  for (i = 0; i < param_count; ++i)
    free(pname_mcmc[i]);
  free(pname_mcmc);

  for (i = 0; i < k; ++i)
    free(pname[i]);
  free(pname);

  for (i = 0; i < token_count; ++i)
  {
    free(tokens[i]);
    free(numeric_tokens[i]);
  }
  free(tokens);
  free(numeric_tokens);

  for (i = 0; i < cols; ++i)
    free(matrix[i]);
  free(matrix);

  free(wmean);
  free(wstdev);
  free(mean);
  free(stdev);
  free(et025);
  free(et975);
  free(hpd025);
  free(hpd975);
  free(tint);
  free(rho1);
  free(effu);
  free(effy);
  free(c);

  fclose(fp);
}

void allfixed_summary(FILE * fp_out, stree_t * stree)
{
  int prec;
  long i, j, count;
  long sample_num;
  long rc = 0;
  FILE * fp;
  unsigned int snodes_total = stree->tip_count + stree->inner_count;
  
  if (opt_msci)
    snodes_total += stree->hybrid_count;

  char * hdr[] = {"Node-Index", "Node-Type", "Node-Label"};
  char * ntypes[] = {"Tip", "Root", "Hybrid", "Inner"};

  /* TODO: pretty-fy output */
  int index_digits = MAX((int)strlen(hdr[0]),(int)(floor(log10(snodes_total+1)+1)));

  /* print I T L */
  fprintf(stdout, "%*s", index_digits, hdr[0]);
  fprintf(fp_out, "%*s", index_digits, hdr[0]);
  fprintf(stdout, "  %s  %s\n", hdr[1], hdr[2]);
  fprintf(fp_out, "  %s  %s\n", hdr[1], hdr[2]);


  long linewidth = index_digits+4+(long)strlen(hdr[1])+(long)strlen(hdr[2]);
  for (j = 0; j < linewidth; ++j)
  {
    fprintf(stdout,"-");
    fprintf(fp_out,"-");
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
    

  for (i = 0; i < snodes_total; ++i)
  {
    char * ntype;
    if (i < stree->tip_count)
      ntype = ntypes[0];
    else if (stree->nodes[i] == stree->root)
      ntype = ntypes[1];
    else if (stree->nodes[i]->hybrid)
      ntype = ntypes[2];
    else
      ntype = ntypes[3];

    if (strchr(stree->nodes[i]->label,','))
    {
      fprintf(stdout, "%-*ld  %-*s  MRCA( %s )\n",
              index_digits,i+1,
              (int)strlen(hdr[1]),ntype,stree->nodes[i]->label);
      fprintf(fp_out, "%-*ld  %-*s  MRCA( %s )\n",
              index_digits,i+1,
              (int)strlen(hdr[1]),ntype,stree->nodes[i]->label);
    }
    else
    {
      if (opt_msci && stree->nodes[i]->hybrid)
      {
        fprintf(stdout, "%-*ld  %-*s  %s (parental node index %d)\n",
                index_digits,i+1,
                (int)strlen(hdr[1]),ntype,stree->nodes[i]->label,
                stree->nodes[i]->parent->node_index+1);
        fprintf(fp_out, "%-*ld  %-*s  %s (parental node index: %d)\n",
                index_digits,i+1,(
                int)strlen(hdr[1]),ntype,stree->nodes[i]->label,
                stree->nodes[i]->parent->node_index+1);
      }
      else
      {
        fprintf(stdout, "%-*ld  %-*s  %s\n",
                index_digits,i+1,
                (int)strlen(hdr[1]),ntype,stree->nodes[i]->label);
        fprintf(fp_out, "%-*ld  %-*s  %s\n",
                index_digits,i+1,
                (int)strlen(hdr[1]),ntype,stree->nodes[i]->label);
      }
    }
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  fp = xopen(opt_mcmcfile,"r");
  /* skip line containing header */
  getnextline(fp);
  assert(strlen(line) > 4);
  char * header = xstrdup(line);

  long token_count;
  char ** tokens = header_explode(header,&token_count);
  header_tokens_shorten(tokens,opt_keep_labels,token_count);

  /* compute number of columns in the file */
  long col_count = 0;

  /* compute number of theta parameters */
  if (opt_est_theta)
    for (i = 0; i < snodes_total; ++i)
      if (stree->nodes[i]->theta >= 0 && stree->nodes[i]->linked_theta == NULL)
        col_count++;

  /* compute number of tau parameters */
  for (i = 0; i < stree->inner_count; ++i) {
    if (stree->nodes[stree->tip_count+i]->tau) {
      col_count++;
      if (opt_est_locusrate == MUTRATE_ONLY && opt_datefile) 
     	 col_count++;
    }
    
  }

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
  if ((opt_est_mubar && opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL) || opt_datefile )
    ++col_count;

  /* add one more for log-L if usedata is on */
  if (opt_usedata)
    col_count++;

  /* column for vbar (or nu for Gamma-Dirichlet) */
  if (opt_clock != BPP_CLOCK_GLOBAL)
    ++col_count;

  /* number of columns including converted Ms (from Ws) */
  long cols = opt_migration ?
                           col_count+opt_migration_count : col_count;
  /* allocate storage matrix */
  double ** matrix = (double **)xmalloc((size_t)cols * sizeof(double *));
  for (i = 0; i < cols; ++i)
    matrix[i] = (double *)xmalloc((size_t)opt_samples * sizeof(double));

  double * mean = (double *)xmalloc((size_t)cols * sizeof(double));

  double * hpd025 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * hpd975 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * tint = (double *)xmalloc((size_t)cols * sizeof(double));
  double * rho1 = (double *)xmalloc((size_t)cols * sizeof(double));
  double * stdev = (double *)xmalloc((size_t)cols * sizeof(double));

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

  char * label[] =
  {
    "mean",       /*  0 */
    "median",     /*  1 */
    "S.D",        /*  2 */
    "min",        /*  3 */
    "max",        /*  4 */
    "2.5%",       /*  5 */
    "97.5%",      /*  6 */
    "2.5%HPD",    /*  7 */
    "97.5%HPD",   /*  8 */
    "ESS*",       /*  9 */
    "Eff*",        /* 10 */
    "rho1"        /* 11 */
  };

  long label_count = sizeof(label) / sizeof(label[0]);
  int * label_size = (int *)xcalloc((size_t)label_count,sizeof(int));
  int pname_size = 5;  /* param */
  prec = 6;
  if (opt_datefile || opt_print_locus)
    prec = 12;

  for (i = 0; i < label_count; ++i)
    label_size[i] = MAX(prec, (int)strlen(label[i]));


  /* now created additional columns for M */
  if (opt_migration)
  {
    char ** new_tokens = NULL;
    new_tokens = (char **)xmalloc((size_t)cols * sizeof(char *));

    /* copy tokens */
    for (i = 0; i < col_count; ++i)
      new_tokens[i] = tokens[i];
    free(tokens);
    tokens = new_tokens;

    /* move lnL to the end */
    long start = col_count;
    if (!strcmp(tokens[col_count-1],"lnL"))
    {
      tokens[cols-1] = tokens[col_count-1];
      for (i = 0; i < opt_samples; ++i)
        matrix[cols-1][i] = matrix[col_count-1][i];
      start = col_count-1;
    }

    /* get column indices of Ws and corresponding thetas */
    long * w_indices;
    long * theta_indices;
    w_indices = (long *)xmalloc((size_t)opt_migration_count*sizeof(long));
    theta_indices = (long *)xmalloc((size_t)opt_migration_count*sizeof(long));

    long w_count = 0;
    long theta_count = 0;
    for (i = 0; i < col_count; ++i)
    {
      if (strlen(tokens[i]) > 2 && tokens[i][0] == 'W' && tokens[i][1] == ':')
      {
        w_indices[w_count++] = i;
        char * tgt_index = strchr(tokens[i],'>');
        assert(tgt_index);
        assert(tgt_index[1]);

        char * nodestr;
        xasprintf(&nodestr, "theta:%s", tgt_index+1);
        for (j = 0; j < col_count; ++j)
        {
          if (!strcmp(tokens[j],nodestr))
          {
            theta_indices[theta_count++] = j;
            break;
          }
        }
        assert(j < col_count);
        free(nodestr);
      }
    }
    assert(w_count == opt_migration_count);
    assert(theta_count == w_count);

    /* go through Ws and calculate Ms */
    for (i = 0; i < opt_migration_count; ++i)
    {
      for (j = 0; j < opt_samples; ++j)
      {
        matrix[start+i][j] = matrix[w_indices[i]][j] *
                             matrix[theta_indices[i]][j]/4.;
      }
      tokens[start+i] = xstrdup(tokens[w_indices[i]]);
      tokens[start+i][0] = 'M';
    }
    free(w_indices);
    free(theta_indices);

    /* update number of columns and number of tokens (headers) */
    col_count = cols;
    token_count += opt_migration_count;
  }


  for (i = 0; i < col_count; ++i)
  {
    long median_line = opt_samples / 2;
    int digits;

    /* mean */
    double sum = 0;
    for (j = 0; j < opt_samples; ++j)
      sum += matrix[i][j];
    mean[i] = sum/opt_samples;

    digits = (int)floor(log10(xfloor1(mean[i]))+1);
    digits += prec+1;
    label_size[0] = MAX(label_size[0],digits);

    /* standard deviation */
    double sd = 0;
    for (j = 0; j < opt_samples; ++j)
      sd += (matrix[i][j]-mean[i]) * (matrix[i][j]-mean[i]);
    stdev[i] = sqrt(sd/(opt_samples-1));

    digits = (int)floor(log10(xfloor1(stdev[i]))+1);
    digits += prec+1;
    label_size[2] = MAX(label_size[2],digits);

    /* tint */
    double _rho1 = 0;
    tint[i] = eff_ict(matrix[i],opt_samples,mean[i],stdev[i], &_rho1);
    rho1[i] = _rho1;

    /* qsort */
    qsort(matrix[i], opt_samples, sizeof(double), cb_cmp_double);

    /* median */
    double median = matrix[i][median_line];
    if ((opt_samples & 1) == 0)
    {
      median += matrix[i][median_line-1];
      median /= 2;
    }
    digits = (int)floor(log10(xfloor1(median))+1);
    digits += prec+1;
    label_size[1] = MAX(label_size[1],digits);

    /* min */
    digits = (int)floor(log10(xfloor1(matrix[i][0]))+1);
    digits += prec+1;
    label_size[3] = MAX(label_size[3],digits);

    /* max */
    digits = (int)floor(log10(xfloor1(matrix[i][opt_samples-1]))+1);
    digits += prec+1;
    label_size[4] = MAX(label_size[4],digits);

    /* 2.5% */
    digits = (int)floor(log10(xfloor1(matrix[i][(long)(opt_samples*.025)]))+1);
    digits += prec+1;
    label_size[5] = MAX(label_size[5],digits);

    /* 97.5% */
    digits = (int)floor(log10(xfloor1(matrix[i][(long)(opt_samples*.975)]))+1);
    digits += prec+1;
    label_size[6] = MAX(label_size[6],digits);

    /* fill 95% HPD arrays */
    hpd_interval(matrix[i],opt_samples,hpd025+i,hpd975+i,0.05);

    /* 2.5% HPD */
    digits = (int)floor(log10(xfloor1(hpd025[i]))+1);
    digits += prec+1;
    label_size[7] = MAX(label_size[7],digits);

    /* 97.5% HPD */
    digits = (int)floor(log10(xfloor1(hpd975[i]))+1);
    digits += prec+1;
    label_size[8] = MAX(label_size[8],digits);

    /* ESS */
    digits = (int)floor(log10(xfloor1(opt_samples/tint[i]))+1);
    digits += prec+1;
    label_size[9] = MAX(label_size[9],digits);

    /* Eff */
    digits = (int)floor(log10(xfloor1(1/tint[i]))+1);
    digits += prec+1;
    label_size[10] = MAX(label_size[10],digits);

    /* rho1 */
    digits = (int)floor(log10(xfloor1(rho1[i]))+1);
    digits += prec+1;
    label_size[11] = MAX(label_size[11],digits);

    /* get column name */
    pname_size = MAX(pname_size,strlen(tokens[i]));
  }

  /* print centered header */
  int linesize = pname_size;
  char * s = center("param", pname_size);
  fprintf(stdout, "%s",s);
  fprintf(fp_out, "%s",s);
  free(s);


  for (i = 0; i < label_count; ++i)
  {
    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");
    s = center(label[i],label_size[i]);
    fprintf(stdout, "%s", s);
    fprintf(fp_out, "%s", s);
    free(s);

    linesize += 2 + label_size[i];
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
  for (i = 0; i < linesize; ++i)
  {
    fprintf(stdout, "-");
    fprintf(fp_out, "-");
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  /* print each row */
  for (i = 0; i < col_count; ++i)
  {
    if (!strcmp(tokens[i],"lnL") && i == col_count-1)
    {
      fprintf(stdout, "\n");
      fprintf(fp_out, "\n");
    }

    fprintf(stdout, "%-*s", pname_size,tokens[i]);
    fprintf(fp_out, "%-*s", pname_size,tokens[i]);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* mean */
    xasprintf(&s, "%.*f", prec, mean[i]);
    fprintf(stdout, "%*s", label_size[0], s);
    fprintf(fp_out, "%*s", label_size[0], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* median */
    long median_line = opt_samples / 2;
    double median = matrix[i][median_line];
    if ((opt_samples & 1) == 0)
    {
      median += matrix[i][median_line-1];
      median /= 2;
    }
    xasprintf(&s, "%.*f", prec, median);
    fprintf(stdout, "%*s", label_size[1], s);
    fprintf(fp_out, "%*s", label_size[1], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* S.D */
    xasprintf(&s, "%.*f", prec, stdev[i]);
    fprintf(stdout, "%*s", label_size[2], s);
    fprintf(fp_out, "%*s", label_size[2], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* min */
    xasprintf(&s, "%.*f", prec, matrix[i][0]);
    fprintf(stdout, "%*s", label_size[3], s);
    fprintf(fp_out, "%*s", label_size[3], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* max */
    xasprintf(&s, "%.*f", prec, matrix[i][opt_samples-1]);
    fprintf(stdout, "%*s", label_size[4], s);
    fprintf(fp_out, "%*s", label_size[4], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* 2.5% */
    xasprintf(&s, "%.*f", prec, matrix[i][(long)(opt_samples*.025)]);
    fprintf(stdout, "%*s", label_size[5], s);
    fprintf(fp_out, "%*s", label_size[5], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* 97.5% */
    xasprintf(&s, "%.*f", prec, matrix[i][(long)(opt_samples*.975)]);
    fprintf(stdout, "%*s", label_size[6], s);
    fprintf(fp_out, "%*s", label_size[6], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* 2.5% HPD */
    xasprintf(&s, "%.*f", prec, hpd025[i]);
    fprintf(stdout, "%*s", label_size[7], s);
    fprintf(fp_out, "%*s", label_size[7], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* 97.5% HPD */
    xasprintf(&s, "%.*f", prec, hpd975[i]);
    fprintf(stdout, "%*s", label_size[8], s);
    fprintf(fp_out, "%*s", label_size[8], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* ESS */
    xasprintf(&s, "%.*f", prec, opt_samples/tint[i]);
    fprintf(stdout, "%*s", label_size[9], s);
    fprintf(fp_out, "%*s", label_size[9], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* Eff */
    xasprintf(&s, "%.*f", prec, 1/tint[i]);
    fprintf(stdout, "%*s", label_size[10], s);
    fprintf(fp_out, "%*s", label_size[10], s);
    free(s);

    fprintf(stdout, "  ");
    fprintf(fp_out, "  ");

    /* rho1 */
    xasprintf(&s, "%.*f", prec, rho1[i]);
    fprintf(stdout, "%*s", label_size[11], s);
    fprintf(fp_out, "%*s", label_size[11], s);
    free(s);

    fprintf(stdout, "\n");
    fprintf(fp_out, "\n");

  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  free(label_size);

  /* success */
  rc = 1;

l_unwind:
  for (i = 0; i < token_count; ++i)
    free(tokens[i]);
  free(tokens);
  free(header);

  for (i = 0; i < col_count; ++i)
    free(matrix[i]);
  free(matrix);

  if (rc && stree->tip_count > 1)
  {
    /* write figtree file */
    write_figtree(fp_out, stree, mean, hpd025, hpd975);
    if (!opt_msci)
      fprintf(stdout, "\nFigTree tree is in %s.FigTree.tre\n", opt_jobname);
    else 
      fprintf(stdout, "\nFigTree tree is in %s.FakeTree.tre\n", opt_jobname);
  }

  free(mean);
  free(hpd025);
  free(hpd975);
  free(tint);
  free(rho1);
  free(stdev);
  fclose(fp);

  if (!rc)
    fatal("Error while reading/summarizing %s", opt_mcmcfile);

  if (opt_a1b1file)
  {
    fprintf(stdout, "\nSummarizing parameter estimates using file %s ...\n\n", opt_a1b1file);
    fprintf(fp_out, "\nSummarizing parameter estimates using file %s ...\n\n", opt_a1b1file);
    summarize_a1b1(stree,fp_out);
  }
}
