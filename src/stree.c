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

static void populations_create(rtree_t * stree, int aln_count)
{
  unsigned int i;

  stree_data_t * sdata;
  rnode_t ** nodes = stree->nodes;

  /* create one species structure for every node */
  for (i = 0; i < stree->inner_count + stree->tip_count; ++i)
  {
    sdata = (stree_data_t *)xcalloc(1,sizeof(stree_data_t));
    nodes[i]->data = (void *)sdata;
  }

  /* populations get the same labels as species tree tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    sdata = (stree_data_t *)(nodes[i]->data);
    sdata->label = xstrdup(nodes[i]->label);
    printf("tip node species label: %s\n", sdata->label);
  }

  /* labels at inner nodes are concatenated labels of their children */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    sdata = (stree_data_t *)(nodes[i]->data);

    /* get the species structures of the two children */
    stree_data_t * ldata = (stree_data_t *)(nodes[i]->left->data);
    stree_data_t * rdata = (stree_data_t *)(nodes[i]->right->data);

    /* allocate necessary memory for label */
    sdata->label = (char *)xmalloc(strlen(ldata->label)+
                                   strlen(rdata->label)+1);

    /* concatenate */
    sdata->label[0] = 0;
    strcat(sdata->label,ldata->label);
    strcat(sdata->label,rdata->label);

    printf("inner node species label: %s\n", sdata->label);
  }

  /* create perloci sequence counters for tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    sdata = (stree_data_t *)(nodes[i]->data);

    sdata->seq_count = (int *)xmalloc(aln_count * sizeof(int));
    memset(sdata->seq_count, 0, aln_count * sizeof(int));
  }
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

static void stree_populate(rtree_t * stree, msa_t ** msalist, list_t * maplist, int count)
{
  int i,j;
  rnode_t * node;

  /* create population structures */
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

      node = (rnode_t *)(pair->data);

#ifdef DEBUG
      printf("Matched %s -> %s\n", label, node->label);
#endif

      /* increment sequence count for loci i in corresponding population */
      stree_data_t * data = (stree_data_t *)(node->data);
      data->seq_count[i]++;
    }
  }

#ifdef DEBUG
  printf("Population per-loci sequence count\n");
  for (i = 0; i < stree->tip_count; ++i)
  {
    stree_data_t * data = (stree_data_t *)(stree->nodes[i]->data);
    printf("%s\n", data->label);
    for (j = 0; j < count; ++j)
    {
      printf("  %d -> %d\n", j, data->seq_count[j]);
    }
  }
#endif


  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);
}

static void stree_init_tau_recursive(rnode_t * node, double prop)
{
  /* end recursion if node is a tip */
  if (!node->left) return;

  /* get species record associate with species tree node */
  stree_data_t * sdata = (stree_data_t *)(node->data);
  stree_data_t * parent_sdata = (stree_data_t *)(node->parent->data);

  double tau_parent = parent_sdata->tau;

  sdata->tau = tau_parent * (prop + (1 - prop - 0.02)*legacy_rndu());
  printf("tau: %f\n", sdata->tau);

  stree_init_tau_recursive(node->left,prop);
  stree_init_tau_recursive(node->right,prop);
}

static void stree_init_tau(rtree_t * stree)
{
  /* Initialize speciation times for each extinct species */ 

  double prop = (stree->root->leaves > PROP_THRESHOLD) ? 0.9 : 0.5;

  /* get species record associate with species tree root*/
  stree_data_t * sdata = (stree_data_t *)(stree->root->data);

  /* set the speciation time for root */
  sdata->tau = opt_tau_beta / (opt_tau_alpha - 1) * (0.9 + 0.2*legacy_rndu());
  printf("tau root: %f\n", sdata->tau);

  /* recursively set the speciation time for the remaining inner nodes */
  stree_init_tau_recursive(stree->root->left,prop);
  stree_init_tau_recursive(stree->root->right,prop);
}

static void stree_init_theta(rtree_t * stree, int msa_count)
{
  /* initialize population sizes for extinct populations and populations
     with more than one lineage at some locus */

  int j;
  unsigned int i,k = 0;

  /* go through tip nodes and setup thetas only for those that have
     two sequences in some loci */
  for (i = 0; i < stree->tip_count; ++i)
  {
    stree_data_t * sdata = (stree_data_t *)(stree->nodes[i]->data);

    for (j = 0; j < msa_count; ++j)
      if (sdata->seq_count[j] >= 2)
        break;

    /* if no loci exists with two or more sequences of such species then move
       to the next tip node */
    if (j == msa_count) continue;

    /* otherwise set theta around the mean of the inverse gamma prior */
    sdata->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                   (0.9 + 0.2 * legacy_rndu());
    ++k;
    printf("theta tip: %f\n", sdata->theta);
  }

  /* go through inner nodes and setup thetas */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    stree_data_t * sdata = (stree_data_t *)(stree->nodes[i]->data);
    sdata->theta = opt_theta_beta / (opt_theta_alpha - 1) *
                   (0.9 + 0.2 * legacy_rndu());
    ++k;
    printf("theta: %f\n", sdata->theta);
  }

  printf("Updated theta: %d\n", k);
}

void stree_init(rtree_t * stree, msa_t ** msa, list_t * maplist, int msa_count)
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

