/*
    Copyright (C) 2016-2024 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static long dmodels_count = 0;       /* number of species delimitation models */
static long dmodel_curindex = 0;
static long * hist = NULL;
static double * dprior = NULL;
static char ** dmodels = NULL;

static double * prior_A11 = NULL;

static snode_t ** trav = NULL;
static int trav_size = 0;

const char * const * bs_last_addr = NULL;

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static double facn;

static int longint_len(long x)
{
  return x ? (int)floor(log10(labs(x)))+1 : 1;
}

double * getpriorA11()
{
  return prior_A11;
}

static double tree_counter(long tip_count, long rooted)
{
  long i;
  double trees = 1;
  rooted = !!rooted;

  for (i = 4; i <= tip_count+rooted; ++i)
    trees *= 2*i-5;

  return trees;
}

static double lh_counter(long tip_count)
{
  long i;
  double lh = 1;

  for (i = 3; i <= tip_count; ++i)
    lh *= i*(i-1) / 2;

  return lh;
}

static double factorial(long n)
{
  long i;
  double f = 1;

  for (i = 2; i <= n; ++i)
    f *= i;
    
  return f;
}

static void print_pinfo(double * prior, long * a, long n, long k)
{
  long i;
  double nd;

  /* count number of delimitations */
  for (i = 0, nd = facn; i < k; ++i)
    if (a[i] > 1) nd /= factorial(a[i]);

  long start = 0;
  long nsame = 1;

  for (i = 1; i < k; ++i)
  {
    if (a[i] != a[start])
    {
      if (nsame > 1) nd /= factorial(nsame);
      start = i; nsame = 1;
    }
    else ++nsame;
  }
  if (nsame > 1) nd /= factorial(nsame);

  /* calculate number of rooted tress with k tips */
  double tree_count = tree_counter(k,BPP_TRUE);

  /* calculate guide trees */
  double guide_count = 1;
  for (i = 0; i < k; ++i)
    if (a[i] > 2) guide_count *= tree_counter(a[i],BPP_TRUE);

  /* calculate weight for labeled histories if prior defined */
  double wlh = 1;
  if (k > 3 && (opt_delimit_prior == BPP_SPECIES_PRIOR_LH ||
                opt_delimit_prior == BPP_SPECIES_PRIOR_SLH))
    wlh = lh_counter(k) / tree_count;


  printf("%*s  nD = %3.0f  Trees = %4.0f  Guide = %4.0f  pro=%4.0f\n",
         (int)(2*n+5-k*2),
         "",
         nd,
         tree_count,
         guide_count,
         nd*tree_count*guide_count*wlh);

  prior[k-1] += nd*tree_count*guide_count*wlh;
}

static void print_partition(long * a, long k, int spc, double * prior, long n)
{
  long i;

  printf("[%*ld ] ", spc+1,k);
  for (i = 0; i < k; ++i)
    printf(" %ld", a[i]);

  print_pinfo(prior,a,n,k);

}

/* Jerome Kelleher's optimal algorithm for generating ascending partitions.
   See also Table 1 column "Number of Delimitations" in Yang & Rannala (2014
   MBE 12: 3125-3135)
   Number of delimitations should be equal to n!/\prod_k #permutations

   Tomas: Check also the "Permutation of multisets" as described here:
   https://en.wikipedia.org/wiki/Permutation
*/
void partition_fast(long n)
{
  int spc = longint_len(n);
  long i,l;
  long k = 1;
  long x;
  long y = n-1;
  long count = 0;
  long * a;

  printf("Counting delimitations...\n");
  if (!prior_A11)
    prior_A11 = (double *)xcalloc((size_t)(n+1),sizeof(double));
  
  /* initialize facn to factorial of n */
  facn = factorial(n);
  
  a = (long *)xmalloc((size_t)(n+1) * sizeof(long));

  for (i = 0; i <= n; ++i)
    a[i] = i;

  while (k)
  {
    x = a[k-1] + 1;
    k -= 1;

    while (2*x <= y)
    {
      a[k] = x;
      y -= x;
      k += 1;
    }
    l = k+1;
    while (x <= y)
    {
      a[k] = x;
      a[l] = y;
      count++;
      print_partition(a,k+2,spc,prior_A11,n);
      x += 1;
      y -= 1;
    }
    a[k] = x+y;
    y = x + y - 1;
    count++;
    print_partition(a,k+1,spc,prior_A11,n);
  }

  double sum = 0;
  for (i = 0; i < n; ++i)
    sum += prior_A11[i];
  for (i = 0; i < n; ++i)
    prior_A11[i] /= sum;

  free(a);
}

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  memcpy(temp,line,line_size*sizeof(char));
  free(line);
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

static char * cb_serialize_support(const snode_t * node)
{
  char * s = NULL;

  /* inner node */
  if (node->left)
    xasprintf(&s, "#%f", node->support);
  else
    xasprintf(&s, "%s", node->label);

  return s;
}

void delimit_summary(FILE * fp_out, stree_t * stree)
{
  long i,j,np;
  long line_count = 0;
  char model[2048];
  FILE * fp;

  /* open MCMC file for reading */
  fp = xopen(opt_mcmcfile,"r");

  double * posterior = (double *)xcalloc(dmodels_count,sizeof(double));

  assert(stree->inner_count < 2048);

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    stree->nodes[i]->support = 0;

  /* skip first header line */
  getnextline(fp);

  while (getnextline(fp))
  {
    sscanf(line,"%ld\t%ld\t %s", &i, &np, model);

    i = delimit_getindexfromstring(model);
    posterior[i]++;

    line_count++;
  }
  assert(line_count);

  for (i = 0; i < dmodels_count; ++i)
    posterior[i] /= line_count;

  for (i = 0; i < dmodels_count; ++i)
  {
    for (j = 0; j < stree->inner_count; ++j)
      if (dmodels[i][j] == '1')
        stree->nodes[stree->tip_count+j]->support += posterior[i];
  }

  fprintf(stdout, "Summarizing the species-delimitation sample in file %s\n\n",
          opt_mcmcfile);
  fprintf(fp_out, "Summarizing the species-delimitation sample in file %s\n\n",
          opt_mcmcfile);

  fprintf(stdout, "Number of species-delimitation models = %ld\n\n",
          dmodels_count);
  fprintf(fp_out, "Number of species-delimitation models = %ld\n\n",
          dmodels_count);

  fprintf(stdout, "     model    prior    posterior\n");
  fprintf(fp_out, "     model    prior    posterior\n");
  for (i = 0; i < dmodels_count; ++i)
  {
    fprintf(stdout,
            "%4ld %s   %f   %f\n",
            i+1,
            dmodels[i],
            dprior[i],
            posterior[i]);
    fprintf(fp_out,
            "%4ld %s   %f   %f\n",
            i+1,
            dmodels[i],
            dprior[i],
            posterior[i]);
  }

  fprintf(stdout, "\nOrder of ancestral nodes:\n");
  fprintf(fp_out, "\nOrder of ancestral nodes:\n");
  for (j = 0; j < stree->inner_count; ++j)
  {
    fprintf(stdout, "  %s\n", stree->nodes[stree->tip_count+j]->label);
    fprintf(fp_out, "  %s\n", stree->nodes[stree->tip_count+j]->label);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  
  /* print guide tree */
  char * newick = stree_export_newick(stree->root, cb_serialize_support);
  fprintf(stdout,
          "Guide tree with posterior probability for presence of nodes:\n");
  fprintf(fp_out,
          "Guide tree with posterior probability for presence of nodes:\n");
  fprintf(stdout, "%s;\n", newick);
  fprintf(fp_out, "%s;\n", newick);

  free(newick);
  free(posterior);
  fclose(fp);
}

long delimitation_getparam_count()
{
  return dmodels_count;
}

char * delimitation_getparam_string()
{
  return dmodels[dmodel_curindex];
}

long delimitation_getcurindex()
{
  return dmodel_curindex;
}


int cb_strcmp(const void * s1, const void * s2)
{
  const char * key = s1;
  const char * const * arg = s2;

  bs_last_addr = arg;

  return strcmp(key, *arg);
}

long delimit_getindexfromstring(char * model)
{
  bsearch(model,
          dmodels,
          dmodels_count,
          sizeof(char *),
          (int(*)(const void *, const void*))cb_strcmp);
  
  return (char **)bs_last_addr - dmodels;
}

long delimit_getindex(stree_t * stree)
{
  unsigned int i;
  char * model = (char *)xmalloc((stree->inner_count+1) * sizeof(char));

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    model[i - stree->tip_count] = stree->nodes[i]->tau > 0 ? '1' : '0';
  model[stree->tip_count-1] = 0;
  
  
  bsearch(model,
          dmodels,
          dmodels_count,
          sizeof(char *),
          (int(*)(const void *, const void*))cb_strcmp);
    
  free(model);

  return (char **)bs_last_addr - dmodels;
}

void delimit_setindex(long index)
{
  dmodel_curindex = index;
}

void delimit_resetpriors()
{
  long i;
  
  for (i = 0; i < dmodels_count; ++i)
    dprior[i] = 0;
}

static long delimitations_count_recursive(snode_t * node)
{
  long x,y;

  if (!node->left)
    return 1;

  x = delimitations_count_recursive(node->left);
  y = delimitations_count_recursive(node->right);

  return x*y+1;
}

long delimitations_count(stree_t * stree)
{
  return delimitations_count_recursive(stree->root);
}

static void print()
{
  int i;
  const long thread_index = 0;

  char * delimitation = dmodels[dmodels_count++]; 

  for (i=0; i<trav_size; ++i)
    delimitation[i] = (trav[i]->mark[thread_index] & FLAG_MISC) ? '1' : '0';
  delimitation[trav_size] = 0;
}

void delimitation_set(stree_t * stree, long index)
{
  long i;
  char * delimitation = dmodels[index]; 

  /* TODO: 
   
   1. This might be different when we change the species tree 
   2. Change to bitvectors
  
  */
           
  for (i = 0; i < stree->inner_count; ++i)
    trav[i]->tau = (delimitation[i] == '1') ? 1 : 0;

  dmodel_curindex = index;
}

static void explore(snode_t ** start, snode_t ** end)
{
  const long thread_index = 0;

  while (end != start)
  {
    snode_t * curnode = *end;

    if (curnode->parent->mark[thread_index] & FLAG_MISC)
    {
      curnode->mark[thread_index] |= FLAG_MISC;   /* set flag */
      print();

      explore(end, trav+trav_size-1);

      curnode->mark[thread_index] &= ~FLAG_MISC;  /* unset flag */
    }

    end--;
  }
}

static void preorder_recursive(snode_t * node,
                               unsigned int * index,
                               snode_t ** outbuffer)
{
  if (!node->left)
    return;

  outbuffer[*index] = node;

  *index = *index + 1;

  preorder_recursive(node->left,  index, outbuffer);
  preorder_recursive(node->right, index, outbuffer);
}

long delimitations_init(stree_t * stree)
{
  unsigned int index = 0;
  long i;
  const long thread_index = 0;

  dmodels_count = delimitations_count(stree);
  printf("Total species delimitations: %ld\n", dmodels_count);

  trav = (snode_t **)xmalloc(stree->inner_count * sizeof(snode_t *));

  preorder_recursive(stree->root, &index, trav);
  assert(index == stree->inner_count);
  trav_size = index;

  dmodels = (char **)xmalloc(dmodels_count * sizeof(char *));
  for (i = 0; i < dmodels_count; ++i)
    dmodels[i] = (char *)xmalloc((trav_size+1)*sizeof(char));
  hist = (long *)xmalloc(dmodels_count * sizeof(long));

  dprior = (double *)xmalloc(dmodels_count * sizeof(double));

  for (i = 0; i < stree->inner_count; ++i)
    trav[i]->mark[thread_index] &= ~FLAG_MISC;    /* unset flag */

  /* print current delimitation model, i.e. all zeros 000..0 */
  dmodels_count = 0;
  print();

  /* print next model with only one population at the root, i.e. 1000...0 */
  stree->root->mark[thread_index] |= FLAG_MISC;   /* set flag */
  print();

  /* recursively explore all other possibilities */
  explore(trav,trav+trav_size-1);

  stree->root->mark[thread_index] &= ~FLAG_MISC;  /* unset flag */

  #if 0
  for (i = 0; i < dmodels_count; ++i)
    printf("%s\n", dmodels[i]);
  #endif

  /* setup priors */
  for (i = 0; i < dmodels_count; ++i)
  {
    long n;
    double p;

    delimitation_set(stree,i);
    n = hist[i] = histories(stree);

    /* 
       For analysis 10 we use either a uniform prior over all rooted trees
       or a dirichlet prior
       TODO: For analysis 11 we need to add another prior, check
       InitializeDelimitationModel() in old bpp
    */

    switch (opt_delimit_prior)
    {
      case BPP_SPECIES_PRIOR_LH:
        p = 1.0;
        break;

      case BPP_SPECIES_PRIOR_UNIFORM:
        p = 1.0 / n;
        break;

      default:
        fatal("Unknown species delimitation prior");
    }

    p = (p < 1e-300) ? -500 : log(p);
    dprior[i] = exp(p) * n;
  }

  double sum = 0;
  for (i = 0; i < dmodels_count; ++i)
    sum += dprior[i];

  double norm = 1.0 / sum;
  for (i = 0; i < dmodels_count; ++i)
    dprior[i] *= norm;
  
  for (i = 0; i < dmodels_count; ++i)
    printf("delimitation model %3ld: %s  LH: %ld   prior  %f\n",
           i, dmodels[i], hist[i],dprior[i]);

  return dmodels_count;
}

void delimitations_fini()
{
  long i;

  if (opt_method == METHOD_10)
  {
    free(trav);
    for (i = 0; i < dmodels_count; ++i)
      free(dmodels[i]);
    free(dmodels);
    free(hist);
    free(dprior);
  }

  if (opt_method == METHOD_11)
  {
    free(prior_A11);
  }
}

static double binomial(double n, int k, double * scale)
{
   double c = 1, i, large = 1e99;

   *scale=0;
   if((int)k!=k) 
      fatal("k is not a whole number in Binomial.");
   if (k == 0) return(1);
   if (n>0 && (k<0 || k>n)) return (0);

   if(n>0 && (int)n==n) k = MIN(k, (int)n - k);
   for (i = 1; i <= k; i++) {
      c *= (n - k + i) / i;
      if (c > large) {
         *scale += log(c); c = 1;
      } 
   }
   return(c);

}

static int fill_lr(snode_t * node, int * lr)
{
  unsigned int i = node->node_index;
  int rc = 0;

  if (lr[i] != -1)
    return rc;

  if (!node->left || node->tau == 0) 
  {
    lr[i] = 0;
    rc = 1;
  }
  else if (lr[node->left->node_index] != -1 && lr[node->right->node_index] != -1)
  {
    lr[i] = lr[node->left->node_index] + lr[node->right->node_index] + 1;
    rc = 1;
  }

  return rc;
}

static int histories_recursive(snode_t * snode, int * lr)
{
  int rc = 0;

  if (!snode->left)
    return rc;

  rc |= histories_recursive(snode->left,lr);
  rc |= histories_recursive(snode->right,lr);

  if (!snode->tau)
    return rc;

  rc |= fill_lr(snode->left,lr);
  rc |= fill_lr(snode->right,lr);

  return rc;
}

double lnprior_species_model(stree_t * stree)
{
  double p;

  switch (opt_delimit_prior)
  {
    case BPP_SPECIES_PRIOR_SLH:
    case BPP_SPECIES_PRIOR_LH:
      p = 1.0;
      break;

    case BPP_SPECIES_PRIOR_SUNIFORM:
    case BPP_SPECIES_PRIOR_UNIFORM:
      p = 1.0 / histories(stree);
      break;

    default:
      fatal("Unknown species delimitation prior");
  }

  /* TODO: A11 */
  if (opt_method == METHOD_11 && opt_delimit_prior >= 2)
  {
    long i;
    long tau_count = 0;

    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      if (stree->nodes[i]->tau > 0)
        ++tau_count;

    p /= prior_A11[tau_count];
  }

  p = (p < 1e-300) ? -500 : log(p);
  
  return p;
}


long histories(stree_t * stree)
{
  unsigned int i,j,k;
  double n = 1;
  double y;

  int * lr = (int *)xmalloc((size_t)(stree->tip_count + stree->inner_count) * 
                            sizeof(int));

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    lr[i] = -1;

  for (i = 0; i < stree->tip_count; ++i)
  {
    if (!histories_recursive(stree->root,lr))
      break;
  }

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->tau == 0) continue;

    j = stree->nodes[i]->left->node_index;
    k = stree->nodes[i]->right->node_index;

    if (lr[j] == -1 || lr[k] == -1)
      fatal("Internal error when computing delimitation histories");

    if (lr[j] && lr[k])
    {
      n *= binomial(lr[j] + lr[k], lr[j], &y);
      if (y)
        fatal("y not expected");
    }

  }
  free(lr);
  return (long)n;
}
