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

static long dmodels_count = 0;
static long dmodel_curindex = 0;
static long * hist = NULL;
static double * dprior = NULL;
static char ** dmodels = NULL;

static snode_t ** trav = NULL;
static int trav_size = 0;

const char * const * bs_last_addr = NULL;

long delimitation_getparam_count()
{
  return dmodels_count;
}

char * delimitation_getparam_string()
{
  return dmodels[dmodel_curindex];
}


int cb_strcmp(const void * s1, const void * s2)
{
  const char * key = s1;
  const char * const * arg = s2;

  bs_last_addr = arg;

  return strcmp(key, *arg);
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

  char * delimitation = dmodels[dmodels_count++]; 

  for (i=0; i<trav_size; ++i)
    delimitation[i] = trav[i]->mark ? '1' : '0';
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
  while (end != start)
  {
    snode_t * curnode = *end;

    if (curnode->parent->mark)
    {
      curnode->mark = 1;
      print();

      explore(end, trav+trav_size-1);

      curnode->mark = 0;
    }

    end--;
  }
}

static void preorder_recursive(snode_t * node,
                               unsigned int * index,
                               snode_t ** buffer)
{
  if (!node->left)
    return;

  buffer[*index] = node;

  *index = *index + 1;

  preorder_recursive(node->left,  index, buffer);
  preorder_recursive(node->right, index, buffer);
}

long delimitations_init(stree_t * stree)
{
  unsigned int index = 0;
  long i;

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
    trav[i]->mark = 0;

  /* print current delimitation model, i.e. all zeros 000..0 */
  dmodels_count = 0;
  print();

  /* print next model with only one population at the root, i.e. 1000...0 */
  stree->root->mark = 1;
  print();

  /* recursively explore all other possibilities */
  explore(trav,trav+trav_size-1);

  stree->root->mark = 0;

  for (i = 0; i < dmodels_count; ++i)
    printf("%s\n", dmodels[i]);


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
      case BPP_DELIMIT_PRIOR_DIRICHLET:
        p = 1.0;
        break;

      case BPP_DELIMIT_PRIOR_UNIFORM:
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
    printf("delimitation model %3ld: %s  prior  %f\n",
           i, dmodels[i], dprior[i]);

  return dmodels_count;
}

void delimitations_fini()
{
  long i;

  free(trav);
  for (i = 0; i < dmodels_count; ++i)
    free(dmodels[i]);
  free(dmodels);
  free(hist);
  free(dprior);
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

long histories(stree_t * stree)
{
  unsigned int i;

  int * lr = (int *)xmalloc((size_t)(stree->tip_count + stree->inner_count) * 
                            sizeof(int));

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    lr[i] = -1;

  for (i = 0; i < stree->tip_count; ++i)
  {
    if (!histories_recursive(stree->root,lr))
      break;
  }

  double n = 1;
  double y;
//  printf("\n");
  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    unsigned int j,k;

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
