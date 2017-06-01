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

/* association of gene nodes to species populations. offset[i] is the 'nodes'
   index at which the gene tree nodes for lineages inside species i are
   available. count dictates the number of available lineages */

typedef struct pop_s
{
  snode_t * snode;

  unsigned int seq_count;
  int * seq_indices;
  gnode_t ** nodes;

} pop_t;

static hashtable_t * sht;
static hashtable_t * mht;

static void gtree_traverse_postorder(gnode_t * node,
                                     int (*cbtrav)(gnode_t *),
                                     unsigned int * index,
                                     gnode_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;

  gtree_traverse_postorder(node->left, cbtrav, index, outbuffer);
  gtree_traverse_postorder(node->right, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

static void gtree_traverse_preorder(gnode_t * node,
                                    int (*cbtrav)(gnode_t *),
                                    unsigned int * index,
                                    gnode_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;

  outbuffer[*index] = node;
  *index = *index + 1;

  gtree_traverse_preorder(node->left, cbtrav, index, outbuffer);
  gtree_traverse_preorder(node->right, cbtrav, index, outbuffer);

}

int gtree_traverse(gnode_t * root,
                   int traversal,
                   int (*cbtrav)(gnode_t *),
                   gnode_t ** outbuffer,
                   unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return BPP_FAILURE;

  /* we will traverse an unrooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  if (traversal == TREE_TRAVERSE_POSTORDER)
    gtree_traverse_postorder(root, cbtrav, trav_size, outbuffer);
  else if (traversal == TREE_TRAVERSE_PREORDER)
    gtree_traverse_preorder(root, cbtrav, trav_size, outbuffer);
  else
    fatal("Invalid traversal value.");

  return BPP_SUCCESS;
}
static void dealloc_data(gnode_t * node,
                         void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

void gtree_destroy(gtree_t * tree,
                   void (*cb_destroy)(void *))
{
  unsigned int i;
  gnode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];
    dealloc_data(node,cb_destroy);

    if (node->label)
      free(node->label);

    free(node);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static char * export_newick_recursive(const gnode_t * root,
                                      char * (*cb_serialize)(const gnode_t *))
{
  char * newick;
  int size_alloced;
  assert(root != NULL);

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = strlen(newick);
    }
    else
    {
      size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = export_newick_recursive(root->left,cb_serialize);
    if (subtree1 == NULL)
    {
      return NULL;
    }
    char * subtree2 = export_newick_recursive(root->right,cb_serialize);
    if (subtree2 == NULL)
    {
      free(subtree1);
      return NULL;
    }

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s:%f",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("Memory allocation during newick export failed.");

  return newick;
}

char * gtree_export_newick(const gnode_t * root,
                           char * (*cb_serialize)(const gnode_t *))
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = strlen(newick);
    }
    else
    {
      size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = export_newick_recursive(root->left,cb_serialize);
    if (subtree1 == NULL)
      fatal("Unable to allocate enough memory.");

    char * subtree2 = export_newick_recursive(root->right,cb_serialize);
    if (subtree2 == NULL)
      fatal("Unable to allocate enough memory.");

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s:%f;",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("Memory allocation during newick export failed");

  return newick;
}


static void fill_nodes_recursive(gnode_t * node,
                                 gnode_t ** array,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!node->left)
  {
    array[*tip_index] = node;
    *tip_index = *tip_index + 1;
    return;
  }

  fill_nodes_recursive(node->left,  array, tip_index, inner_index);
  fill_nodes_recursive(node->right, array, tip_index, inner_index);

  array[*inner_index] = node;
  *inner_index = *inner_index + 1;
}

static unsigned int gtree_count_tips(gnode_t * root)
{
  unsigned int count = 0;

  if (root->left)
    count += gtree_count_tips(root->left);
  if (root->right)
    count += gtree_count_tips(root->right);

  if (!root->left && !root->right)
    return 1;

  return count;
}

static gtree_t * gtree_wraptree(gnode_t * root,
                                unsigned int tip_count)
{
  unsigned int i;

  gtree_t * tree = (gtree_t *)xmalloc(sizeof(gtree_t));

  if (tip_count < 2 && tip_count != 0)
    fatal("Invalid number of tips in input tree (%u).", tip_count);

  if (tip_count == 0)
  {
    /* if tip counts is set to 0 then recursively count the number of tips */
    tip_count = gtree_count_tips(root);
    if (tip_count < 2)
    {
      fatal("Input tree contains no inner nodes.");
    }
  }

  tree->nodes = (gnode_t **)xmalloc((2*tip_count-1)*sizeof(gnode_t *));
  
  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;

  fill_nodes_recursive(root->left, tree->nodes, &tip_index, &inner_index);
  fill_nodes_recursive(root->right,tree->nodes, &tip_index, &inner_index);
  tree->nodes[inner_index] = root;

  tree->tip_count = tip_count;
  tree->edge_count = 2*tip_count-2;
  tree->inner_count = tip_count-1;
  tree->root = root;

  for (i = 0; i < 2*tip_count-1; ++i)
    tree->nodes[i]->node_index = i;

  return tree;
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

static int cb_cmp_pairlabel(void * a, void * b)
{
  pair_t * pair = (pair_t *)a;
  char * label = (char *)b;

  return (!strcmp(pair->label,label));
}

static void fill_pop(pop_t * pop, stree_t * stree, msa_t * msa, int msa_id)
{
  unsigned int i,j;

    /* go through the sequences of the current locus and match each sequence
       with the corresponding population using its species tag */
  for (j = 0; j < stree->tip_count; ++j)
    pop[j].snode = stree->nodes[j];

  for (j = 0; j < (unsigned int)(msa->count); ++j)
  {
    /* first get the species tag */
    char * label = msa->label[j];
    label = strchr(label, '^');
    if (!label)
      fatal("Cannot find species tag on sequence %s of locus %d",
            msa->label[j], msa_id);

    /* skip the '^' mark */
    label++;
    if (!(*label))
      fatal("Sequence %s of locus %d contains no label",
            msa->label[j], msa_id);

    pair_t * pair;
    pair = hashtable_find(mht,
                          (void *)label,
                          hash_fnv(label),
                          cb_cmp_pairlabel);
    if (!pair)
      fatal("Cannot find species mapping for sequence %s of locus %d",
            label, msa_id);

    snode_t * node = (snode_t *)(pair->data);

    i = node->node_index;
    assert(node == pop[i].snode);

    /* increment sequence count for current loci in population i */
    pop[i].seq_count++;
  }

  /* now repeat that for obtaining the sequence indices */
  int * counter = (int *)xcalloc(stree->tip_count, sizeof(int));
  for (i = 0; i < stree->tip_count; ++i)
  {
    pop[i].seq_indices = (int *)xcalloc(pop[i].seq_count,sizeof(int));
    pop[i].nodes = (gnode_t **)xcalloc(pop[i].seq_count,sizeof(gnode_t *));
  }
  
  for (j = 0; j < (unsigned int)(msa->count); ++j)
  {
    /* first get the species tag */
    char * label = msa->label[j];
    label = strchr(label, '^');
    if (!label)
      fatal("Cannot find species tag on sequence %s of locus %d",
            msa->label[j], msa_id);

    /* skip the '^' mark */
    label++;
    if (!(*label))
      fatal("Sequence %s of locus %d contains no label",
            msa->label[j], msa_id);

    pair_t * pair;
    pair = hashtable_find(mht,
                          (void *)label,
                          hash_fnv(label),
                          cb_cmp_pairlabel);
    if (!pair)
      fatal("Cannot find species mapping for sequence %s of locus %d",
            label, msa_id);

    snode_t * node = (snode_t *)(pair->data);
    i = node->node_index;
    assert(node == pop[i].snode);

    /* increment sequence count for loci i in corresponding population */
    int * k = counter+i;
    pop[i].seq_indices[*k] = j;
    (*k)++;
  }

  free(counter);
}

static void replace(pop_t * pop, int count, snode_t * epoch)
{
  int i,j;

  /* delete the two children of epoch from the pop list, by replacing the child
     with the smaller index with epoch, and deleting the second child. It also
     merges the sequences (lineages) of the two children into the new population
  */

  /* replace left descendant of current epoch with epoch */
  for (i = 0; i < count; ++i)
  {
    if (pop[i].snode == epoch->left)
      break;
  }
  assert(i != count);

  /* delete right descendant of epoch */
  for (j = 0; j < count; ++j)
  {
    /* replace indices and nodes */
    if (pop[j].snode == epoch->right)
      break;
  }
  assert( j != count);

  /* order them such that i is the smaller index */
  if (j < i)
    SWAP(i,j);

  /* allocate indices and nodes arrays for the new population */
  int * indices = (int *)xmalloc((pop[i].seq_count+pop[j].seq_count) *
                                 sizeof(int));
  gnode_t ** nodes = (gnode_t **)xmalloc((pop[i].seq_count+pop[j].seq_count) *
                                         sizeof(gnode_t *));

  /* fill indices and nodes by merging the information from its two children
     populations */
  memcpy(indices,pop[i].seq_indices,pop[i].seq_count*sizeof(int));
  memcpy(indices+pop[i].seq_count,
         pop[j].seq_indices,
         pop[j].seq_count*sizeof(int));
  memcpy(nodes,pop[i].nodes,pop[i].seq_count*sizeof(gnode_t *));
  memcpy(nodes+pop[i].seq_count,
         pop[j].nodes,
         pop[j].seq_count*sizeof(gnode_t *));
  pop[i].seq_count += pop[j].seq_count;

  /* deallocate indices and nodes array of the two children populations */
  free(pop[i].seq_indices);
  free(pop[i].nodes);
  free(pop[j].seq_indices);
  free(pop[j].nodes);

  /* replace population i with new population */
  pop[i].snode = epoch;
  pop[i].seq_indices = indices;
  pop[i].nodes = nodes;

  /* if population j was not the last one, replace the last population in the
     list with j */
  if (j < count-1)
    memcpy(pop+j,pop+count-1,sizeof(pop_t));
}

static int cb_cmp_spectime(const void * a, const void * b)
{
  
  snode_t * const * x = a;
  snode_t * const * y = b;

  if ((*x)->tau - (*y)->tau > 0) return 1;
  return -1;
}

static gtree_t * gtree_simulate(stree_t * stree, msa_t * msa, int msa_id)
{
  int lineage_count = 0;
  unsigned int i,j,k;
  unsigned int epoch_count;
  double t, tmax, sum;
  double * ci;
  pop_t * pop;
  snode_t ** epoch;
  gnode_t * inner = NULL;

  /* get a list of inner nodes (epochs) */
  epoch = (snode_t  **)xmalloc(stree->inner_count*sizeof(snode_t *));
  memcpy(epoch,
         stree->nodes + stree->tip_count,
         stree->inner_count * sizeof(snode_t *));

  /* sort epochs in ascending order of speciation times */
  qsort(&(epoch[0]), stree->inner_count, sizeof(snode_t *), cb_cmp_spectime);
  epoch_count = stree->inner_count;

  /* create one hash table of species and one for sequence->species mappings */
  pop = (pop_t *)xcalloc(stree->tip_count, sizeof(pop_t));
  fill_pop(pop,stree,msa,msa_id);

  /* start at present time */
  t = 0;

  lineage_count = msa->count;

  /* allocate space for storing coalescent rates for each population */
  ci = (double *)xmalloc(stree->tip_count * sizeof(double));

  /* current epoch index */
  unsigned int e = 0;

  /* create a list of tip nodes for the target gene tree */
  gnode_t ** gtips = (gnode_t **)xcalloc(msa->count,
                                         sizeof(gnode_t *));
  for (i = 0; i < (unsigned int)(msa->count); ++i)
    gtips[i] = (gnode_t *)xcalloc(1,sizeof(gnode_t));

  /* fill each population with one gene tip node for each lineage */
  for (i = 0, j=0; i < stree->tip_count; ++i)
  {
    memcpy(pop[i].nodes,gtips+j,pop[i].seq_count*sizeof(gnode_t *));
    j += pop[i].seq_count;
  }

  /* set the clv index for each gene tip node. The index must be equal to the
     number of the sequence it represents in the current msa, and will be used
     later for setting up and computing the CLVs */
  for (i = 0; i < stree->tip_count; ++i)
  {
    for (j = 0; j < pop[i].seq_count; ++j)
    {
      int index = pop[i].seq_indices[j];
      pop[i].nodes[j]->label = xstrdup(msa->label[index]);
      pop[i].nodes[j]->clv_index = index;
    }
  }
  
  /* initialize counter for clv indices for inner nodes. They will have CLV
     indices in the range  [tip_count,2*tip_count-1] */
  unsigned int clv_index = msa->count;

  /* initialize the number of currently available populations for choosing
     to coalesce lineages in */
  unsigned int pop_count = stree->tip_count;

  /* TODO: The below loop is written to match exactly the loop in the original
     BPP. It is possible to implement a simpler, more understandable routine but
     it might change the structure of the randomly generated gene trees */

  /* loop until we are left with only 1 lineage in one ancestral population */
  for (; ; --pop_count)
  {
      
    while (1)
    {
      /* calculate poisson rates: ci[j] is coalescent rate for population j */
      for (j=0, sum=0; j < pop_count; ++j)
      {
        k = pop[j].seq_count;

        if (k >= 2)
        {
          ci[j] = k*(k-1)/pop[j].snode->theta;
          sum += ci[j];
        }
      }

      /* set max waiting time for this epoch */
      tmax = epoch[e]->tau;

      /* generate random waiting time from exponential distribution */
      t += legacy_rndexp(1/sum);

      /* if the generated time is larger than the current epoch, and we are not
         yet at the root of the species tree, then break and, subsequently, 
         merge the lineages of the two populations into the current epoch */
      if (t > tmax && pop_count != 1) break;

#ifdef DEBUG_GTREE_SIMULATE
      fprintf(stdout, "[Debug]: Coalescent waiting time: %f\n", t);
#endif

      /* TODO: Implement migration routine */

      /* select an available population at random using the poisson rates as
         weights */
      double r = legacy_rndu()*sum;
      double tmp = 0;
      for (j = 0; j < pop_count; ++j)
      {
        tmp += ci[j];
        if (r < tmp) break;
      }

      assert(j < pop_count);

      /* now choose two lineages from selected population j in exactly the same
         way as the original BPP */
      k = pop[j].seq_count * (pop[j].seq_count-1) * legacy_rndu();

      unsigned int k1 = k / (pop[j].seq_count-1);
      unsigned int k2 = k % (pop[j].seq_count-1);

      if (k2 >= k1)
        k2++;
      else
        SWAP(k1,k2);

#ifdef DEBUG_GTREE_SIMULATE
      fprintf(stdout, "[Debug]: Coalesce (%d,%d) into %d\n",
              pop[j].nodes[k1]->clv_index,
              pop[j].nodes[k2]->clv_index,
              clv_index);
#endif

      pop[j].seq_count--;
      
      /* allocate and fill new inner node as the parent of the gene tree nodes
         representing lineages k1 and k2 */
      inner = (gnode_t *)xcalloc(1,sizeof(gnode_t));
      inner->parent = NULL;
      inner->left  = pop[j].nodes[k1];
      inner->right = pop[j].nodes[k2];
      inner->clv_index = clv_index++;
      inner->left->parent = inner;
      inner->right->parent = inner;
      inner->time = t;

      /* in pop j, replace k1 by the new node, and remove k2 (replace it with
         last node in list) */
      pop[j].seq_indices[k1] = -1;   /* TODO HERE clv vector */
      pop[j].nodes[k1] = inner;

      if (k2 != pop[j].seq_count)
      {
        pop[j].seq_indices[k2] = pop[j].seq_indices[pop[j].seq_count];
        pop[j].nodes[k2] = pop[j].nodes[pop[j].seq_count];
      }
      
      /* break if this was the last available lineage in the current locus */
      if (--lineage_count == 1) break;
    }

    t = tmax;

    if (pop_count == 1 || lineage_count == 1) break;

    /* place current epoch in the list of populations, remove its two children
       and add up the lineages of the two descendants */
    replace(pop,pop_count,epoch[e]);
    
    if (e != epoch_count-1)
      ++e;
  }

  for(i = 0; i < pop_count; ++i)
  {
    free(pop[i].seq_indices);
    free(pop[i].nodes);
  }

  /* wrap the generated tree structure (made up of linked nodes) into gtree_t */
  gtree_t * gtree = gtree_wraptree(inner, (unsigned int)(msa->count));

#ifdef DEBUG_GTREE_SIMULATE
  /* create a newick string from constructed gene tree */
  for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
    if (gtree->nodes[i] != gtree->root)
      gtree->nodes[i]->length = gtree->nodes[i]->parent->time -
                                gtree->nodes[i]->time;

  fprintf(stdout, "[Debug]: Printing newick string for current gene tree\n");
  char * newick = gtree_export_newick(gtree->root,NULL);
  fprintf(stdout,"%s\n", newick);
  free(newick);
#endif

  /* cleanup */
  free(gtips);
  free(pop);
  free(ci);
  free(epoch);

  return gtree;
}

gtree_t ** gtree_init(stree_t * stree,
                      msa_t ** msalist,
                      list_t * maplist,
                      int msa_count)
{
  int i;
  gtree_t ** gtree;

  assert(msa_count > 0);

  gtree = (gtree_t **)xmalloc(msa_count*sizeof(gtree_t *));

  /* create mapping hash tables */
  sht = species_hash(stree);
  mht = maplist_hash(maplist,sht);

  /* generate random starting gene trees for each alignment */
  for (i = 0; i < msa_count; ++i)
    gtree[i] = gtree_simulate(stree, msalist[i],i);

  /* destroy the hash tables */
  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);

  return gtree;
}
