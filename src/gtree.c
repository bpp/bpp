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

static int replace_item(pop_t ** list, pop_t * from, pop_t * to, int count)
{
  int i;

  for (i = 0; i < count; ++i)
  {
    if (list[i] == from)
    {
      list[i] = to;
      break;
    }
  }

  return i != count;
}

static int remove_item(pop_t ** list, pop_t * item, int count)
{
  int i;

  for (i = 0; i < count; ++i)
  {
    if (list[i] == item)
    {
      if (i < count-1)
        list[i] = list[count-1];
      break;
    }
  }
  return (i != count);
}

static int cb_cmp_spectime(const void * a, const void * b)
{
  
  pop_t * const * x = a;
  pop_t * const * y = b;

  printf("sorting %f %f\n", (*x)->tau, (*y)->tau);

  if ((*x)->tau - (*y)->tau > 0) return 1;
  return -1;
}

void gtree_simulate(rtree_t * stree,
                    pop_t ** poplist,
                    msa_t ** msalist,
                    int msa_id)
{
  unsigned int i,j,k;
  unsigned int population_count;
  unsigned int epoch_count;
  unsigned int lineage_count = 0;
  double t, tmax, sum;
  double * ci;
  pop_t ** epoch;
  pop_t ** population;

  int n = msalist[msa_id]->count;

  /* get a list of inner nodes (epochs) and sort them in ascending order of
     speciation times */
  #if 1
  printf("Getting %d epochs\n", stree->inner_count);
  #endif
  epoch = (pop_t **)xmalloc(stree->inner_count*sizeof(pop_t *));
  memcpy(epoch,
         poplist + stree->tip_count,
         stree->inner_count * sizeof(pop_t *));

  #if 1
  for (i = 0; i < stree->inner_count; ++i)
    printf("printing %f\n", epoch[i]->tau);
  #endif

  qsort(&(epoch[0]), stree->inner_count, sizeof(pop_t *), cb_cmp_spectime);
  epoch_count = stree->inner_count;

  #if 1
  printf("Epochs:\n");
  for (i = 0; i < epoch_count; ++i)
  {
    printf("\t%s: tau %f theta %f\n", epoch[i]->label, epoch[i]->tau, epoch[i]->theta);
  }
  #endif

  /* current active populations are the extant species */
  population = (pop_t **)xmalloc(stree->tip_count*sizeof(pop_t *));
  memcpy(population, poplist, stree->tip_count*sizeof(pop_t *));
  population_count = stree->tip_count;

  /* start at present time */
  t = 0;

  /* count total number of lineages for this locus */
  for (i = 0; i < stree->tip_count; ++i)
  {
    lineage_count += population[i]->seq_count[msa_id];
  }

  /* allocate space for storing coalescent rates for each population */
  ci = (double *)xmalloc(population_count * sizeof(double));

  /* current epoch index */
  unsigned int e = 0;

  rnode_t ** nodelist = (rnode_t **)xcalloc(msalist[msa_id]->count,
                                            sizeof(rnode_t *));
  for (i = 0; i < msalist[msa_id]->count; ++i)
    nodelist[i] = (rnode_t *)xcalloc(1,sizeof(rnode_t));

  for (; ; --population_count)
  {
      
    while (1)
    {
      /* calculate poisson rates: ci[j] is coalescent rate for population j */
      for (j=0, sum=0; j < population_count; ++j)
      {
        k = population[j]->seq_count[msa_id];

        if (k >= 2)
        {
          ci[j] = k*(k-1)/population[j]->theta;
          sum += ci[j];
          printf("ci[%d] = %f (k: %d theta: %f)\n", j, ci[j], k, population[j]->theta);
        }
      }

      tmax = epoch[e]->tau;
        

      printf("sum: %f (e = %d)\n", sum, e);
      t += legacy_rndexp(1/sum);
      printf("t: %f\n", t);

      /* if the generated time is larger than the current epoch, and we have at
         least two populations, then break and merge the lineages of the two
         populations in the current epoch */
      if (t > tmax && population_count != 1) break;

      /* TODO: if (C+M < 1e-300) {} */

      /* we will coalesce two lineages from a population chosen at random
         with the rates as probabilities. First choose the population */
      double r = legacy_rndu()*sum;
      double tmp = 0;
      for (j = 0; j < population_count; ++j)
      {
        tmp += ci[j];
        if (r < tmp) break;
      }

      assert(j < population_count);

      /* now choose two lineages from population j */
      k = population[j]->seq_count[msa_id] * (population[j]->seq_count[msa_id]-1) * legacy_rndu();

      int k1,k2;

      k1 = k / (population[j]->seq_count[msa_id]-1);
      k2 = k % (population[j]->seq_count[msa_id]-1);
      if (k2 >= k1)
        k2++;
      else
        SWAP(k1,k2);

      printf("DEBUG: Merging %d and %d  (k=%d)\n", k1,k2,k);

      if (population[j]->node_index < stree->tip_count)
      {
        if (population[j]->seq_indices[msa_id][k1] != -1)
          printf("\t\t %d = %s\n", k1,msalist[msa_id]->label[population[j]->seq_indices[msa_id][k1]]);
        if (population[j]->seq_indices[msa_id][k2] != -1)
          printf("\t\t %d = %s\n", k2, msalist[msa_id]->label[population[j]->seq_indices[msa_id][k2]]);
        population[j]->seq_indices[msa_id][k1] = -1;
      }
      /* In pop i, replace k1 by new node, remove k2 */
      population[j]->seq_count[msa_id]--;
      if (population[j]->node_index < stree->tip_count) 
      {
        if (k2 != population[j]->seq_count[msa_id])
          population[j]->seq_indices[msa_id][k2] = population[j]->seq_indices[msa_id][population[j]->seq_count[msa_id]];
        }
      /*
      if (k2 != population[j]->seq_count[msa_id])
        do something
      */


      if (--n == 1) break;
    }

    t = tmax;

    if (population_count == 1 || n == 1) break;

    /* place current epoch in the list of populations and remove
       its two children */
    //if (!remove_item(population,epoch[i]->left,population_count))

    
    /* get left and right descendant populations of epoch[e] */
    pop_t * lpop= poplist[stree->nodes[epoch[e]->node_index]->left->pop_index];
    pop_t * rpop= poplist[stree->nodes[epoch[e]->node_index]->right->pop_index];

    if (!replace_item(population,lpop,epoch[e],population_count))
      fatal("Internal error during gene tree construction");
    if (!remove_item(population,rpop,population_count))
      fatal("Internal error during gene tree construction");

    /* add up the lineages of the two descendant population into epoch[e] */
    epoch[e]->seq_count[msa_id] = lpop->seq_count[msa_id]+rpop->seq_count[msa_id];

    if (e != epoch_count-1)
      ++e;
  }

  free(population);
  free(ci);
  free(epoch);
  free(nodelist);
}
