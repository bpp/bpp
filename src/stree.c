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
  rnode_t * node = (rnode_t *)a;
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
static hashtable_t * species_hash(rtree_t * tree)
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

    rnode_t * node = hashtable_find(sht,
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

static pop_t ** populations_create(rtree_t * stree, int msa_count)
{
  unsigned int i;

  rnode_t ** nodes = stree->nodes;

  pop_t ** poplist = (pop_t **)xmalloc((stree->tip_count + stree->inner_count) *
                                       sizeof(pop_t *));


  /* create one species structure for every node */
  for (i = 0; i < stree->inner_count + stree->tip_count; ++i)
  {
    poplist[i] = (pop_t *)xcalloc(1, sizeof(pop_t));
    stree->nodes[i]->pop_index = i;
    poplist[i]->node_index = i;
  }

  /* populations get the same labels as species tree tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    poplist[i]->label = xstrdup(nodes[i]->label);
    printf("tip node species label: %s\n", poplist[i]->label);
  }

  /* labels at inner nodes are concatenated labels of their children */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    /* get the species structures of the two children */
    pop_t * lpop = poplist[nodes[i]->left->pop_index];
    pop_t * rpop = poplist[nodes[i]->right->pop_index];

    /* allocate necessary memory for label */
    poplist[i]->label = (char *)xmalloc(strlen(lpop->label)+
                                        strlen(rpop->label)+1);

    /* concatenate */
    poplist[i]->label[0] = 0;
    strcat(poplist[i]->label,lpop->label);
    strcat(poplist[i]->label,rpop->label);

    printf("inner node species label: %s\n", poplist[i]->label);
  }

  /* create perloci sequence counters for tip and inner nodes */
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    poplist[i]->seq_count = (int *)xcalloc(msa_count,sizeof(int));
    poplist[i]->seq_indices = (int **)xcalloc(msa_count,sizeof(int *));
  }


  return poplist;
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

static pop_t ** stree_populate(rtree_t * stree, msa_t ** msalist, list_t * maplist, int count)
{
  int i,j,k;
  rnode_t * node;
  pop_t ** poplist;

  /* create a list of populations and associate them with species tree nodes */
  poplist = populations_create(stree, count);

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

      node = (rnode_t *)(pair->data);

#ifdef DEBUG
      printf("Matched %s -> %s\n", label, node->label);
#endif

      /* increment sequence count for loci i in corresponding population */
      pop_t * pop = poplist[node->pop_index];
      pop->seq_count[i]++;
    }
    for (j = 0; j < stree->tip_count; ++j)
    {
      pop_t * pop = poplist[j];
      pop->seq_indices[i] = (int *)xcalloc(pop->seq_count[i],sizeof(int));
    }
  }

  /* repeat the exact thing, but now fill sequence indices instead */

  for (i = 0; i < count; ++i) 
  {
    msa_t * msa = msalist[i];
    
    k = 0;
    int * counter = (int *)calloc(stree->tip_count, sizeof(int));
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

      node = (rnode_t *)(pair->data);

#ifdef DEBUG
      printf("Matched %s -> %s\n", label, node->label);
#endif

      /* increment sequence count for loci i in corresponding population */
      pop_t * pop = poplist[node->pop_index];
      pop->seq_indices[i][counter[node->pop_index]++] = j;
    }
    free(counter);
  }

#ifdef DEBUG
  printf("Population per-loci sequence count\n");
  for (i = 0; i < stree->tip_count; ++i)
  {
    printf("%s\n", poplist[i]->label);
    for (j = 0; j < count; ++j)
    {
      printf("  %d -> %d\n", j, poplist[i]->seq_count[j]);
    }
  }
#endif


  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);

  return poplist;
}

static void stree_init_tau_recursive(rnode_t * node,
                                     pop_t ** poplist,
                                     double prop)
{
  /* end recursion if node is a tip */
  if (!node->left) return;

  /* get species record associate with species tree node */
  pop_t * pop = poplist[node->pop_index];
  double tau_parent = poplist[node->parent->pop_index]->tau;

  pop->tau = tau_parent * (prop + (1 - prop - 0.02)*legacy_rndu());
  printf("tau: %f\n", pop->tau);

  stree_init_tau_recursive(node->left,poplist,prop);
  stree_init_tau_recursive(node->right,poplist,prop);
}

static void stree_init_tau(rtree_t * stree, pop_t ** poplist)
{
  /* Initialize speciation times for each extinct species */ 

  double prop = (stree->root->leaves > PROP_THRESHOLD) ? 0.9 : 0.5;

  /* get species record associate with species tree root*/
  pop_t * pop = poplist[stree->root->pop_index];

  /* set the speciation time for root */
  pop->tau = opt_tau_beta / (opt_tau_alpha - 1) * (0.9 + 0.2*legacy_rndu());
  printf("tau root: %f\n", pop->tau);

  /* recursively set the speciation time for the remaining inner nodes */
  stree_init_tau_recursive(stree->root->left,poplist,prop);
  stree_init_tau_recursive(stree->root->right,poplist,prop);
}

static void stree_init_theta(rtree_t * stree, pop_t ** poplist, int msa_count)
{
  /* initialize population sizes for extinct populations and populations
     with more than one lineage at some locus */

  int j;
  unsigned int i,k = 0;

  /* go through tip nodes and setup thetas only for those that have
     two sequences in some loci */
  for (i = 0; i < stree->tip_count; ++i)
  {
    for (j = 0; j < msa_count; ++j)
      if (poplist[i]->seq_count[j] >= 2)
        break;

    /* if no loci exists with two or more sequences of such species then move
       to the next tip node */
    if (j == msa_count) continue;

    /* otherwise set theta around the mean of the inverse gamma prior */
    poplist[i]->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                        (0.9 + 0.2 * legacy_rndu());
    ++k;
    printf("theta tip: %f\n", poplist[i]->theta);
  }

  /* go through inner nodes and setup thetas */
  /* TODO: Note, that we compute the theta for the root first to be in line
     with the original bpp code */
  pop_t * pop = poplist[stree->root->pop_index];
  pop->theta = opt_theta_beta / (opt_theta_alpha - 1) *
               (0.9 + 0.2 * legacy_rndu());
  ++k;
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count-1; ++i)
  {
    poplist[i]->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                        (0.9 + 0.2 * legacy_rndu());
    ++k;
    printf("theta: %f\n", poplist[i]->theta);
  }

  printf("Updated theta: %d\n", k);
}

pop_t ** stree_init(rtree_t * stree, msa_t ** msa, list_t * maplist, int msa_count)
{
  /* create populations on each node of the species tree, and compute and store
     the per-locus number of sequences for each species */
  pop_t ** poplist = stree_populate(stree, msa, maplist, msa_count);

  /* Initialize tip count proportion variable */

  /* Initialize population sizes */
  stree_init_theta(stree, poplist, msa_count);

  /* Initialize speciation times and create extinct species groups */
  stree_init_tau(stree,poplist);

  return poplist;
}

