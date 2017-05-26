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

#define PROP_THRESHOLD 10

static int cb_cmp_nodelabel(void * a, void * b)
{
  snode_t * node = (snode_t *)a;
  char * label = (char * )b;

  return (!strcmp(node->label,label));
}

static int cb_cmp_pairpair(void * a, void * b)
{
  pair_t * p1 = (pair_t *)a;
  pair_t * p2 = (pair_t *)b;

  return (!strcmp(p1->label,p2->label));
}

static int cb_cmp_pairlabel(void * a, void * b)
{
  pair_t * pair = (pair_t *)a;
  char * label  = (char *)b;

  return (!strcmp(pair->label,label));
}

#if 0
/* species node data callback deallocator */
void cb_stree_dealloc(void * data)
{
  if (!data) return;

  stree_data_t * sdata = (stree_data_t *)data;

  if (sdata->seq_count)
    free(sdata->seq_count);

  if (sdata->label)
    free(sdata->label);

  free(data);
}
#endif

/* hashtable for indexing species tree labels */
static hashtable_t * species_hash(stree_t * tree)
{
  long i;

  /* using a hash table check whether there are duplicate nodes */
  hashtable_t * ht = hashtable_create(tree->tip_count);

  for (i = 0; i < tree->tip_count; ++i)
  {
    /* attempt to place the pair in the hash table and die with an
       error if a pair with the same label already exists */
    if (!hashtable_insert(ht,
                          (void *)(tree->nodes[i]),
                          hash_fnv(tree->nodes[i]->label),
                          hashtable_ptrcmp))
    {
      /* this should never happen because duplicate taxa were
         already checked for during tree parsing */
      fatal("Duplicate taxon or node label (%s) in file %s",
            tree->nodes[i]->label,
            opt_streefile);
    }
  }

  return ht;
}

hashtable_t * maplist_hash(list_t * maplist, hashtable_t * sht)
{
  if (!maplist) return NULL;

  list_item_t * li = maplist->head;

  hashtable_t * ht = hashtable_create(maplist->count);

  while (li)
  {
    mapping_t * mapping = (mapping_t *)(li->data);

    snode_t * node = hashtable_find(sht,
                                    (void *)(mapping->species),
                                    hash_fnv(mapping->species),
                                    cb_cmp_nodelabel);
    if (!node)
      fatal("Cannot find node with population label %s", mapping->species);

    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));

    pair->label = xstrdup(mapping->individual);
    pair->data = (void *)node;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(pair->label),
                          cb_cmp_pairpair))
      fatal("Duplicate mapping (%s -> %s) found in file %s",
            mapping->individual,
            mapping->species,
            opt_mapfile);

    printf("Mapped %s -> %s\n", pair->label, node->label);
    li = li->next;
  }
  
  return ht;
}

void populations_create(stree_t * stree, int msa_count)
{
  unsigned int i;

  #if 0
  /* create one species structure for every node */
  for (i = 0; i < stree->inner_count + stree->tip_count; ++i)
  {
    poplist[i] = (pop_t *)xcalloc(1, sizeof(pop_t));
    stree->nodes[i]->pop_index = i;
    poplist[i]->node_index = i;
  }
  #endif

  #if 0
  /* populations get the same labels as species tree tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    poplist[i]->label = xstrdup(nodes[i]->label);
    printf("tip node species label: %s\n", poplist[i]->label);
  }
  #endif

  /* labels at inner nodes are concatenated labels of their children */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    /* get the species structures of the two children */

    snode_t * node = stree->nodes[i];
    snode_t * lnode = stree->nodes[i]->left;
    snode_t * rnode = stree->nodes[i]->right;

    /* allocate necessary memory for label */
    if (node->label)
      free(node->label);
    node->label = (char *)xmalloc(strlen(lnode->label)+
                                  strlen(rnode->label)+1);

    /* concatenate */
    node->label[0] = 0;
    strcat(node->label,lnode->label);
    strcat(node->label,rnode->label);

    printf("inner node species label: %s\n", node->label);
  }

  /* create perloci sequence counters for tip and inner nodes */
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    
    stree->nodes[i]->seq_count = (int *)xcalloc(msa_count,sizeof(int));
    stree->nodes[i]->seq_indices = (int **)xcalloc(msa_count,sizeof(int *));
  }
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

static void stree_populate(stree_t * stree, msa_t ** msalist, list_t * maplist, int count)
{
  int i,j;
  snode_t * node;

  /* create a list of populations and associate them with species tree nodes */
  populations_create(stree, count);

  /* create one hash table of species and one for sequence->species mappings */
  hashtable_t * sht = species_hash(stree);
  hashtable_t * mht = maplist_hash(maplist,sht);

  /* go through the alignments and match each sequence with the corresponding
     species */
  for (i = 0; i < count; ++i)
  {
    msa_t * msa = msalist[i];

    /* go through the sequences of the current locus and match each sequence
       with the corresponding population using its species tag */


    for (j = 0; j < msa->count; ++j)
    {
      /* first get the species tag */
      char * label = msa->label[j];
      label = strchr(label, '^');
      if (!label)
        fatal("Cannot find species tag on sequence %s of locus %d",
              msa->label[j], i);

      /* skip the '^' mark */
      label++;
      if (!(*label))
        fatal("Sequence %s of locus %d contains no label",
              msa->label[j], i);

      pair_t * pair;
      pair = hashtable_find(mht,
                            (void *)label,
                            hash_fnv(label),
                            cb_cmp_pairlabel);
      if (!pair)
        fatal("Cannot find");

      node = (snode_t *)(pair->data);

#if 1 
      printf("Matched %s -> %s\n", label, node->label);
#endif

      /* increment sequence count for loci i in corresponding population */
      node->seq_count[i]++;
    }
    for (j = 0; j < stree->tip_count; ++j)
    {
      snode_t * node = stree->nodes[j];

      node->seq_indices[i] = (int *)xcalloc(node->seq_count[i],sizeof(int));
    }
  }

  /* repeat the exact thing, but now fill sequence indices instead */

  for (i = 0; i < count; ++i) 
  {
    msa_t * msa = msalist[i];
    
    int * counter = (int *)xcalloc(stree->tip_count, sizeof(int));
    for (j = 0; j < msa->count; ++j)
    {
      /* first get the species tag */
      char * label = msa->label[j];
      label = strchr(label, '^');
      if (!label)
        fatal("Cannot find species tag on sequence %s of locus %d",
              msa->label[j], i);

      /* skip the '^' mark */
      label++;
      if (!(*label))
        fatal("Sequence %s of locus %d contains no label",
              msa->label[j], i);

      pair_t * pair;
      pair = hashtable_find(mht,
                            (void *)label,
                            hash_fnv(label),
                            cb_cmp_pairlabel);
      if (!pair)
        fatal("Cannot find");

      node = (snode_t *)(pair->data);

      /* increment sequence count for loci i in corresponding population */
      int * k = &(counter[node->node_index]);
      node->seq_indices[i][*k] = j;
      (*k)++;
    }
    free(counter);
  }

#ifdef DEBUG
  printf("Population per-loci sequence count\n");
  for (i = 0; i < stree->tip_count; ++i)
  {
    node = stree->nodes[i];

    printf("%s\n", node->label);
    for (j = 0; j < count; ++j)
    {
      printf("  %d -> %d\n", j, node->seq_count[j]);
    }
  }
#endif


  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);
}

static void stree_init_tau_recursive(snode_t * node,
                                     double prop)
{
  /* end recursion if node is a tip */
  if (!node->left) return;

  /* get species record associate with species tree node */
  double tau_parent = node->parent->tau;

  node->tau = tau_parent * (prop + (1 - prop - 0.02)*legacy_rndu());
  printf("tau: %f\n", node->tau);

  stree_init_tau_recursive(node->left,prop);
  stree_init_tau_recursive(node->right,prop);
}

static void stree_init_tau(stree_t * stree)
{
  /* Initialize speciation times for each extinct species */ 

  double prop = (stree->root->leaves > PROP_THRESHOLD) ? 0.9 : 0.5;

  /* set the speciation time for root */
  stree->root->tau = opt_tau_beta/(opt_tau_alpha-1) * (0.9 + 0.2*legacy_rndu());
  printf("tau root: %f\n", stree->root->tau);

  /* recursively set the speciation time for the remaining inner nodes */
  stree_init_tau_recursive(stree->root->left,prop);
  stree_init_tau_recursive(stree->root->right,prop);
}

static void stree_init_theta(stree_t * stree, int msa_count)
{
  /* initialize population sizes for extinct populations and populations
     with more than one lineage at some locus */

  int j;
  unsigned int i,k = 0;

  /* go through tip nodes and setup thetas only for those that have
     two sequences in some loci */
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode_t * node = stree->nodes[i];

    for (j = 0; j < msa_count; ++j)
      if (node->seq_count[j] >= 2)
        break;

    /* if no loci exists with two or more sequences of such species then move
       to the next tip node */
    if (j == msa_count) continue;

    /* otherwise set theta around the mean of the inverse gamma prior */
    node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                        (0.9 + 0.2 * legacy_rndu());
    ++k;
    printf("theta tip: %f\n", node->theta);
  }

  /* go through inner nodes and setup thetas */
  /* TODO: Note, that we compute the theta for the root first to be in line
     with the original bpp code */
  stree->root->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                       (0.9 + 0.2 * legacy_rndu());
  ++k;
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count-1; ++i)
  {
    snode_t * node = stree->nodes[i];

    node->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                  (0.9 + 0.2 * legacy_rndu());
    ++k;
    printf("theta: %f\n", node->theta);
  }

  printf("Updated theta: %d\n", k);
}

void stree_init(stree_t * stree, msa_t ** msa, list_t * maplist, int msa_count)
{
  /* create populations on each node of the species tree, and compute and store
     the per-locus number of sequences for each species */
  stree_populate(stree, msa, maplist, msa_count);

  /* Initialize tip count proportion variable */

  /* Initialize population sizes */
  stree_init_theta(stree, msa_count);

  /* Initialize speciation times and create extinct species groups */
  stree_init_tau(stree);
}

