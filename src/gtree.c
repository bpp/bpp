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

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))


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

static double * sortbuffer = NULL;

static gnode_t *** travbuffer = NULL;

#if 0

/* 
   This functions are purely for debugging. I use them to update all partials
   (even if not necessary) and to reset the 'leaves' property of gene tree
   nodes, in order to ensure that computations are correct
*/

static void all_partials_recursive(gnode_t * node,
                                   unsigned int * trav_size,
                                   gnode_t ** outbuffer)
{
  if (!node->left)
    return;

  all_partials_recursive(node->left,  trav_size, outbuffer);
  all_partials_recursive(node->right, trav_size, outbuffer);

  outbuffer[*trav_size] = node;
  *trav_size = *trav_size + 1;
}

static void gtree_all_partials(gnode_t * root,
                               gnode_t ** travbuffer,
                               unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return;

  all_partials_recursive(root, trav_size, travbuffer);
}

static void gtree_reset_leaves_recursive(gnode_t * node)
{
  if (!node->left) return;
    node->leaves = 1;
    

  gtree_reset_leaves_recursive(node->left);
  gtree_reset_leaves_recursive(node->right);
  
  node->leaves = node->left->leaves + node->right->leaves;
}

static void gtree_reset_leaves(gnode_t * root)
{
  if (!root->left) return;

  gtree_reset_leaves_recursive(root);
}
#endif

static int return_partials_recursive(gnode_t * node,
                                     unsigned int * trav_size,
                                     gnode_t ** outbuffer)
{
  int mark = 0;

  if (!node->left)
    return node->mark & FLAG_PARTIAL_UPDATE;   /* TODO - SHOULD THIS BE 0 ? */

  mark |= return_partials_recursive(node->left,  trav_size, outbuffer);
  mark |= return_partials_recursive(node->right, trav_size, outbuffer);

  if ((node->mark & FLAG_PARTIAL_UPDATE) || mark)
  {
    outbuffer[*trav_size] = node;
    *trav_size = *trav_size + 1;
  }

  return (node->mark & FLAG_PARTIAL_UPDATE) | mark;
}

gnode_t ** gtree_return_partials(gnode_t * root,
                                 unsigned int msa_index,
                                 unsigned int * trav_size)
{
  gnode_t ** trav = travbuffer[msa_index];

  *trav_size = 0;
  if (!root->left) return NULL;

  return_partials_recursive(root, trav_size, trav);

  return trav;
}


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
      size_alloced = xasprintf(&newick, "%s:%f", root->label, root->length);
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
      size_alloced = xasprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
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
      size_alloced = xasprintf(&newick, "%s:%f", root->label, root->length);
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
      size_alloced = xasprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
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

static void fill_nodes_recursive(gnode_t * node, gnode_t ** array)
{
  array[node->clv_index] = node;

  if (!node->left)
    return;

  fill_nodes_recursive(node->left,array);
  fill_nodes_recursive(node->right,array);
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
  
  fill_nodes_recursive(root, tree->nodes);

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

static void fill_seqin_counts_recursive(snode_t * node, int msa_index)
{
  if (!node->left)
    return;

  fill_seqin_counts_recursive(node->left,msa_index);
  fill_seqin_counts_recursive(node->right,msa_index);

  snode_t * lnode = node->left;
  snode_t * rnode = node->right;

  node->seqin_count[msa_index] = lnode->seqin_count[msa_index] +
                                 rnode->seqin_count[msa_index] -
                                 lnode->event_count[msa_index] -
                                 rnode->event_count[msa_index];
}

void fill_seqin_counts(stree_t * stree, int msa_index)
{
  fill_seqin_counts_recursive(stree->root,msa_index);
}

static int cb_trav_full(snode_t * x)
{
  if (!x->left)
    return 0;

  return 1;
}

static gtree_t * gtree_simulate(stree_t * stree, msa_t * msa, int msa_index)
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
  stree_traverse(stree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_trav_full,
                 epoch,
                 &epoch_count);

  assert(epoch_count == stree->inner_count);

  snode_t ** sortptr = epoch;

  for (j = 0, i = 0; i < stree->inner_count; ++i)
  {
    if (epoch[i]->tau == 0)
    {
      if (i != j)
        SWAP(epoch[i], epoch[j]);
      
      ++j;
    }
  }
  sortptr += j;
  qsort(&(sortptr[0]), stree->inner_count-j, sizeof(snode_t *), cb_cmp_spectime);

  epoch_count = stree->inner_count;

  if (opt_debug)
  {
    for (i = 0; i < stree->inner_count; ++i)
      printf("[Debug]: Epoch %d (%s) time - %f\n",
             i, epoch[i]->label, epoch[i]->tau);
  }

  /* create one hash table of species and one for sequence->species mappings */
  pop = (pop_t *)xcalloc(stree->tip_count, sizeof(pop_t));
  fill_pop(pop,stree,msa,msa_index);

  for (i = 0; i < stree->tip_count; ++i)
    stree->nodes[i]->seqin_count[msa_index] = pop[i].seq_count;
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
  {
    gtips[i] = (gnode_t *)xcalloc(1,sizeof(gnode_t));
    gtips[i]->pmatrix_index = i;
    gtips[i]->scaler_index = PLL_SCALE_BUFFER_NONE;
    gtips[i]->leaves = 1;
  }

  /* fill each population with one gene tip node for each lineage */
  for (i = 0, j=0; i < stree->tip_count; ++i)
  {
    memcpy(pop[i].nodes,gtips+j,pop[i].seq_count*sizeof(gnode_t *));
    for (k = 0; k < pop[i].seq_count; ++k)
      pop[i].nodes[k]->pop = pop[i].snode;

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
    /* set max waiting time for this epoch */
    if (pop_count == 1)
      tmax = -1;
    else
      tmax = epoch[e]->tau;
      
    while (1)
    {
      if (!tmax) break;

      /* calculate poisson rates: ci[j] is coalescent rate for population j */
      for (j=0, sum=0; j < pop_count; ++j)
      {
        k = pop[j].seq_count;

        if (k >= 2)
        {
          ci[j] = k*(k-1)/pop[j].snode->theta;
          sum += ci[j];
        }
        else
          ci[j] = 0;
      }

      if (sum < 1e-300)
        break;

      /* generate random waiting time from exponential distribution */
      t += legacy_rndexp(1/sum);

      /* if the generated time is larger than the current epoch, and we are not
         yet at the root of the species tree, then break and, subsequently, 
         merge the lineages of the two populations into the current epoch */
      if (t > tmax && pop_count != 1) break;

      if (opt_debug)
        fprintf(stdout, "[Debug]: Coalescent waiting time: %f\n", t);

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

      if (opt_debug)
      {
        fprintf(stdout, "[Debug]: Coalesce (%3d,%3d) into %3d (age: %f)\n",
                pop[j].nodes[k1]->clv_index,
                pop[j].nodes[k2]->clv_index,
                clv_index,
                t);
      }

      pop[j].seq_count--;
      
      /* allocate and fill new inner node as the parent of the gene tree nodes
         representing lineages k1 and k2 */
      inner = (gnode_t *)xcalloc(1,sizeof(gnode_t));
      inner->parent = NULL;
      inner->left  = pop[j].nodes[k1];
      inner->right = pop[j].nodes[k2];
      inner->clv_index = clv_index;
      inner->scaler_index = PLL_SCALE_BUFFER_NONE;
      inner->pmatrix_index = clv_index;
      inner->left->parent = inner;
      inner->right->parent = inner;
      inner->time = t;
      inner->pop = pop[j].snode;
      inner->leaves = inner->left->leaves + inner->right->leaves;
      clv_index++;

      pop[j].snode->event_count[msa_index]++;
      dlist_item_t * dlitem = dlist_append(pop[j].snode->event[msa_index],inner);
      inner->event = dlitem;

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
    {
      ++e;
    }
  }

  for(i = 0; i < pop_count; ++i)
  {
    free(pop[i].seq_indices);
    free(pop[i].nodes);
  }

  /* wrap the generated tree structure (made up of linked nodes) into gtree_t */
  gtree_t * gtree = gtree_wraptree(inner, (unsigned int)(msa->count));

  /* create a newick string from constructed gene tree and print it on screen */
  if (opt_debug)
  {
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
      if (gtree->nodes[i] != gtree->root)
        gtree->nodes[i]->length = gtree->nodes[i]->parent->time -
                                  gtree->nodes[i]->time;

    fprintf(stdout, "[Debug]: Printing newick string for current gene tree\n");
    char * newick = gtree_export_newick(gtree->root,NULL);
    fprintf(stdout,"%s\n", newick);
    free(newick);
  }

  /* cleanup */
  free(gtips);
  free(pop);
  free(ci);
  free(epoch);


  /* Update the number of lineages coming into each ancestral population as the
     sum of lineages coming into its two children populations minus the
     coalescent events that have occured */
  fill_seqin_counts(stree,msa_index);

  if (opt_debug)
  {
    fprintf(stdout,
            "[Debug] # coalescent events, lineages coming in locus %d:\n", 
            msa_index);
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      fprintf(stdout, "  %-*s : %-3d %-3d (age: %f)\n", 
              (int)(strlen(stree->root->label)),
              stree->nodes[i]->label,
              stree->nodes[i]->event_count[msa_index],
              stree->nodes[i]->seqin_count[msa_index],
              stree->nodes[i]->tau);
  }

  return gtree;
}

static void reset_gene_leaves_count_recursive(snode_t * node, unsigned int locus_count)
{
  unsigned int j;

  if (!node->left)
    return;

  reset_gene_leaves_count_recursive(node->left,locus_count);
  reset_gene_leaves_count_recursive(node->right,locus_count);

  for (j = 0; j < locus_count; ++j)
    node->gene_leaves[j] = node->left->gene_leaves[j] +
                           node->right->gene_leaves[j];

}
void reset_gene_leaves_count(stree_t * stree)
{
  unsigned int i,j;

  /* gene leaves is the same as sequences coming in for tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
    for (j = 0; j < stree->locus_count; ++j)
      stree->nodes[i]->gene_leaves[j] = stree->nodes[i]->seqin_count[j];

  reset_gene_leaves_count_recursive(stree->root, stree->locus_count);

  for (j = 0; j < stree->locus_count; ++j)
    stree->root->gene_leaves[j] = stree->root->left->gene_leaves[j] +
                                  stree->root->right->gene_leaves[j];

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
  printf("Generating gene trees....");
  for (i = 0; i < msa_count; ++i)
    gtree[i] = gtree_simulate(stree, msalist[i],i);
  printf(" Done\n");

  /* destroy the hash tables */
  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);

  /* allocate sort buffer */
  int max_count = 0;
  for (i = 0; i < msa_count; ++i)
    if (msalist[i]->count > max_count)
      max_count = msalist[i]->count;

  /* alloate buffer for sorting coalescent times plus two for the beginning
     and end of epoch */
  sortbuffer = (double *)xmalloc((max_count+2) * sizeof(double));
  
  /* allocate traversal buffers */
  travbuffer = (gnode_t ***)xmalloc(msa_count * sizeof(gnode_t **));
  for (i = 0; i < msa_count; ++i)
    travbuffer[i] = (gnode_t **)xmalloc(gtree[i]->inner_count *
                                        sizeof(gnode_t *));

  /* reset number of gene leaves associated with each species tree subtree */
  reset_gene_leaves_count(stree);

  return gtree;
}

void gtree_update_branch_lengths(gtree_t ** gtree_list, int count)
{
  unsigned int i,j;

  for (j = 0; j < (unsigned int)count; ++j)
  {
    gtree_t * gtree = gtree_list[j];
    for (i = 0; i < gtree->tip_count+gtree->inner_count; ++i)
      if (gtree->nodes[i] != gtree->root)
        gtree->nodes[i]->length = gtree->nodes[i]->parent->time -
                                  gtree->nodes[i]->time;
  }
}

static int cb_cmp_double_asc(const void * a, const void * b)
{
  const double * x = a;
  const double * y = b;

  if (*x > *y) return 1;

  return -1;
}

double gtree_update_logprob_contrib(snode_t * snode, int msa_index)
{
    unsigned int j,k,n;
    double logpr = 0;
    double heredity = 1;
    double T2h = 0;
    dlist_item_t * event;

    sortbuffer[0] = snode->tau;
    j = 1;
    for (event = snode->event[msa_index]->head; event; event = event->next)
    {
      gnode_t * gnode = (gnode_t *)(event->data);
      sortbuffer[j++] = gnode->time;
    }
    if (snode->parent)
      sortbuffer[j++] = snode->parent->tau;


    /* TODO: Probably split the following qsort case into two:
       in case snode->parent then sort j-2 elements, otherwise
       j-1 elements.
    */

    /* if there was at least one coalescent event, sort */
    if (j > 1)
      qsort(sortbuffer+1,j-1,sizeof(double),cb_cmp_double_asc);

    #if 0
    printf("Population: %s tau: %f theta: %f events: %d seqin_count: %d\n",
           snode->label, snode->tau, snode->theta,
           snode->event_count[msa_index], snode->seqin_count[msa_index]);

    if (snode->parent)
      n = j-1;
    else
      n = j;

    for (k = 1; k < n; ++k)
    {
      printf("\t t = %f\n", sortbuffer[k]);
    }
    if (snode->parent)
      printf("\t tau = %f\n", sortbuffer[k]);
    printf("\n");
    #endif

    /* skip the last step in case the last value of n was supposed to be 1 */
    if ((unsigned int)(snode->seqin_count[msa_index]) == j-1) --j;
    for (k=1,n=snode->seqin_count[msa_index]; k < j; ++k, --n)
    {
      T2h += n*(n-1)*(sortbuffer[k] - sortbuffer[k-1])/heredity;
    }

    if (snode->event_count[msa_index])
      logpr += snode->event_count[msa_index] * log(2.0/snode->theta);

    if (T2h)
      logpr -= T2h/snode->theta;

    /* TODO: Be careful about which functions update the logpr contribution 
       and which do not */
    snode->old_logpr_contrib[msa_index] = snode->logpr_contrib[msa_index];
    snode->logpr_contrib[msa_index] = logpr;

    return logpr;
}

double gtree_logprob(stree_t * stree, int msa_index)
{
  unsigned int i;

  double logpr = 0;

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    logpr += gtree_update_logprob_contrib(stree->nodes[i],msa_index);

  return logpr;
}

double reflect(double t, double minage, double maxage)
{
  int side = 0;
  double n,excess = 0;
  const double EPSILON = 1e-100;

  if (maxage-minage < EPSILON)
    fatal("Internal error when proposing gene tree noge age");


  if (t < minage)
  {
    excess = minage - t;
    side = 0;
  }
  else if (t > maxage)
  {
    excess = t - maxage;
    side = 1;
  }

  if (excess)
  {
    double diff = maxage - minage;

    n = floor(excess / diff);

    if (fmod(n,2.0) > 0.1)
      side = !side;

    excess -= n*diff;

    t = side ? maxage-excess : minage+excess;
  }

  return t;
}

void unlink_event(gnode_t * node, int msa_index)
{
  /* first re-link the event before the current node with the one after */
  if (node->event->prev)
    node->event->prev->next = node->event->next;
  else
    node->pop->event[msa_index]->head = node->event->next;

  /* now re-link the event after the current node with the one before */
  if (node->event->next)
    node->event->next->prev = node->event->prev;
  else
    node->pop->event[msa_index]->tail = node->event->prev;
}

static long propose_ages(locus_t * locus, gtree_t * gtree, stree_t * stree, int msa_index)
{
  unsigned int i,k,j;
  long accepted = 0;
  double tnew,minage,maxage,oldage;
  double logpr;
  double logl;
  snode_t * pop;
  snode_t * oldpop;

  /* TODO: Instead of traversing the gene tree nodes this way, traverse the
     coalescent events for each population in the species tree instead. This
     will reduce the amount of required quick-sorts.
  */

  for (i = gtree->tip_count; i < gtree->inner_count+gtree->tip_count; ++i)
  {
    gnode_t * node = gtree->nodes[i];

    /* find minimum children age */
    minage = MAX(node->left->time,node->right->time);

    if (node->left->pop != node->right->pop)
    {
      snode_t * lpop;
      snode_t * rpop;

      /* find most recent ancestral population to the two child populations */
      /* TODO: Speed this up by using a lookup table */
      for (lpop = node->left->pop; lpop; lpop = lpop->parent)
      {
        for (rpop = node->right->pop; rpop; rpop = rpop->parent)
          if (rpop == lpop) break;
        if (rpop == lpop) break;
      }

      assert(rpop == lpop && rpop != NULL);

      minage = MAX(minage,lpop->tau);
    }

    /* compute max age. TODO: 999 is placed for compatibility with old bpp */
    maxage = node->parent ? node->parent->time : 999;

    assert(maxage > minage);

    tnew = node->time + opt_finetune_gtage * legacy_rnd_symmetrical();
    tnew = reflect(tnew, minage, maxage);

    /* find the first ancestral pop with age higher than the proposed tnew */
    /* TODO: Improvement: probably this can start from lpop/rpop (LCA of
        populations of two daughter nodes) */
    for (pop = node->left->pop; pop->parent; pop = pop->parent)
      if (pop->parent->tau > tnew)
        break;

    /* save old age and old ancestral population */
    oldage = node->time;
    oldpop = node->pop;

    node->time = tnew;

    /* If we need to change the population of node, we have to also update the
       coalescent events list of the current and new population */
    if (node->pop != pop)
    {
      /* remove current gene node from the list of coalescent events of its old
         population */

      unlink_event(node,msa_index);

      /* decrease the number of coalescent events for the current population */
      node->pop->event_count[msa_index]--;
        
      /* change population for the current gene tree node */
      node->pop = pop;

      /* now add the coalescent event to the new population, at the end */
      dlist_item_append(node->pop->event[msa_index],node->event);

      node->pop->event_count[msa_index]++;

      /* increase or decrease the number of incoming lineages to all populations in the path
      from old population to the new population, depending on the case  */
      if (tnew > oldage)
      {
        /* increase the number of incoming lineages to all populations in the
           path from old population (excluding) to the new population */
        for (pop = oldpop; pop != node->pop; pop = pop->parent)
          pop->parent->seqin_count[msa_index]++;
      }
      else
      {
        /* decrease the number of incoming lineages to all populations in the
           path from new population (excluding) to the old population */
        for (pop = node->pop; pop != oldpop; pop = pop->parent)
          pop->parent->seqin_count[msa_index]--;
      }
    }

    /* quick recomputation  of logpr */
    /* IMPORTANT NOTE:
       In the normal runs, DEBUG_LOGPROB is not defined, and *only* the gene
       tree log probability is only updated by re-computing the contributions
       for the coalescent events in the modified populations (partial 
       computation). If DEBUG_LOGPROB is defined, then the whole gene tree
       probability is recomputed. This is useful for debugging and checking that
       the partial computation is functioning properly */
    #ifndef DEBUG_LOGPROB
    logpr = gtree->logpr;
    if (oldpop == node->pop)
    {
      logpr -= node->pop->logpr_contrib[msa_index];
      logpr += gtree_update_logprob_contrib(node->pop,msa_index);
    }
    else
    {
      snode_t * start;
      snode_t * end;

      if (tnew > oldage)
      {
        start = oldpop;
        end   = node->pop->parent;
      }
      else
      {
        start = node->pop;
        end   = oldpop->parent;
      }
      for (pop = start; pop != end; pop = pop->parent)
      {
        logpr -= pop->logpr_contrib[msa_index];
        logpr += gtree_update_logprob_contrib(pop,msa_index);
      }
    }
    #else
      logpr = gtree_logprob(stree,msa_index);
      /* assertion just to remind us that this is not what we should put in
         the official version */
      assert(0);
    #endif

    /* now update branch lengths and prob matrices */
    gnode_t * temp;
    k = 0;
    travbuffer[msa_index][k++] = node->left;
    travbuffer[msa_index][k++] = node->right;
    if (node->parent)
      travbuffer[msa_index][k++] = node;
    locus_update_matrices_jc69(locus,travbuffer[msa_index],k);
      

    /* fill traversal buffer with root-path starting from current node */
    for (k=0, temp = node; temp; temp = temp->parent)
    {
      travbuffer[msa_index][k++] = temp;

      /* swap clv index to compute partials in a new location. This is useful
         when the proposal gets rejected, as we only have swap clv indices */
      temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
    }

    /* update partials */
    locus_update_partials(locus,travbuffer[msa_index],k);
    
    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

    /* acceptance ratio */
    double acceptance = logpr - gtree->logpr + logl - gtree->logl;

    if (opt_debug)
      fprintf(stdout, "[Debug] (age) lnacceptance = %f\n", acceptance);

    if (acceptance >= 0 || legacy_rndu() < exp(acceptance))
    {
      /* accepted */
      accepted++;

      /* update new log-likelihood and gene tree log probability */
      gtree->logpr = logpr;
      gtree->logl = logl;
    }
    else
    {
      /* rejected */

      /* need to reset clv indices to point to the old clv buffer */
      for (j = 0; j < k; ++j)
      {
        temp = travbuffer[msa_index][j];
        temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
      }
      
      /* now reset branch lengths and pmatrices */
      node->time = oldage;
      k = 0;
      travbuffer[msa_index][k++] = node->left;
      travbuffer[msa_index][k++] = node->right;
      if (node->parent)
        travbuffer[msa_index][k++] = node;
      locus_update_matrices_jc69(locus,travbuffer[msa_index],k);

      /* reset to old population, and reset gene tree log probability
         contributes for each modified species tree node */
      if (node->pop == oldpop)
      {
        node->pop->logpr_contrib[msa_index] = node->pop->old_logpr_contrib[msa_index];
      }
      else
      {
        /* remove current gene node from the list of coalescent events of its old
           population */

        unlink_event(node,msa_index);

        /* decrease the number of coalescent events for the current population */
        node->pop->event_count[msa_index]--;
          
        /* change population for the current gene tree node */
        SWAP(node->pop,oldpop);

        /* now add the coalescent event back to the old population, at the end */
        dlist_item_append(node->pop->event[msa_index],node->event);

        node->pop->event_count[msa_index]++;

        /* increase or decrease the number of incoming lineages to all
           populations in the path from old population to the new population,
           depending on the case  */
        if (oldage >= tnew)
        {
          /* increase the number of incoming lineages to all populations in the
             path from old population (excluding) to the new population */
          for (pop=oldpop; pop != node->pop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]++;
        }
        else
        {
          /* decrease the number of incoming lineages to all populations in the
             path from new population (excluding) to the old population */
          for (pop = node->pop; pop != oldpop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]--;
        }

        /* now restore the old log gene tree probability contribution for each
           affected species tree node */
        
        snode_t * start = (tnew > oldage) ? node->pop :  oldpop;
        snode_t * end   = (tnew > oldage) ? oldpop->parent : node->pop->parent;

        for (pop = start; pop != end; pop = pop->parent)
          pop->logpr_contrib[msa_index] = pop->old_logpr_contrib[msa_index];
      }
    }
  }
  return accepted;
}

double gtree_propose_ages(locus_t ** locus, gtree_t ** gtree, stree_t * stree)
{
  unsigned int i;
  long proposal_count = 0;
  long accepted = 0;

  for (i = 0; i < stree->locus_count; ++i)
  {
    /* TODO: Fix this to account mcmc.moveinnode in original bpp */
    proposal_count += gtree[i]->inner_count;
    accepted += propose_ages(locus[i],gtree[i],stree,i);
  }

  if (!accepted)
    return 0;

  return ((double)accepted/proposal_count);

}


static int perform_spr(gtree_t * gtree, gnode_t * curnode, gnode_t * target)
{
  gnode_t * sibling;
  gnode_t * father;
  gnode_t * oldroot = gtree->root;
  int ret = 0;

  sibling = (curnode->parent->left == curnode) ? 
                curnode->parent->right : curnode->parent->left;

  father = curnode->parent;

  assert(father != target);
  assert(target != sibling);

  /*             

                        /\                                          /\
                       /  \                                        /  \
                      /    \                                      /    \
                     /      \                           sibling  *      \
                    /        \                                  / \      \
            father *          \                                           \
                  / \          *                                           *
         curnode *   \        / \                                         / \
                / \   \          \                                           \
               /   \   * sibling  \                                           * father
                      / \          *  target                                 / \
                                  / \                              curnode  *   \
                                                                           / \   * target
                                                                          /   \ / \

  */

  /* if father is root then sibling becomes the new root */
  if (father == gtree->root)
  {
    sibling->leaves = gtree->root->leaves;   /* NEW */
    gtree->root = sibling;
    sibling->parent = NULL;
    ret = 2;
  }
  else
  {
    /* unlink parent from grandfather and link sibling to grandfather */
    if (father->parent->left == father)
      father->parent->left = sibling;
    else
      father->parent->right = sibling;
    sibling->parent = father->parent;

    /* update number of leaves all nodes from father's parent and up */
    gnode_t * temp;
    for (temp = father->parent; temp; temp = temp->parent)
      temp->leaves = temp->left->leaves +
                     temp->right->leaves;
  }

  /* regraft */
  father->parent = target->parent;
  target->parent = father;
  father->left = target;
  father->right = curnode;
  father->leaves = father->left->leaves + father->right->leaves;
  if (target == gtree->root)
  {
    father->leaves = gtree->root->leaves;
    gtree->root = father;
    ret = 1;
  }
  else
  {
    if (father->parent->left == target)
      father->parent->left = father;
    else
      father->parent->right = father;

    /* update number of leaves all nodes from father's parent and up */
    gnode_t * temp;
    for (temp = father->parent; temp; temp = temp->parent)
      temp->leaves = temp->left->leaves +
                     temp->right->leaves;
  }


  if (gtree->root != oldroot)
  {
    SWAP(gtree->root->pop,oldroot->pop);
    SWAP(gtree->root->time,oldroot->time);
    
    gtree->root->left->parent = oldroot;
    gtree->root->right->parent = oldroot;

    oldroot->left->parent = gtree->root;
    oldroot->right->parent = gtree->root;

    SWAP(gtree->root->left,oldroot->left);
    SWAP(gtree->root->right,oldroot->right);

    if (oldroot->left == oldroot)
      oldroot->left = gtree->root;
    if (oldroot->right == oldroot)
      oldroot->right = gtree->root;

    gnode_t * temp = oldroot->parent;

    if (temp->left == oldroot)
      temp->left = gtree->root;
    if (temp->right== oldroot)
      temp->right= gtree->root;

    if ((gtree->root->parent = temp) == gtree->root)
      gtree->root->parent = oldroot;
    oldroot->parent = NULL;


    /* TODO: ERROR */
    //SWAP(gtree->root->pmatrix_index, oldroot->pmatrix_index);
    //SWAP(gtree->root->scaler_index, oldroot->scaler_index);
    //SWAP(gtree->root->clv_valid, oldroot->clv_valid);
    //SWAP(gtree->root->clv_index, oldroot->clv_index);
    SWAP(gtree->root->leaves, oldroot->leaves);
    SWAP(gtree->root->event, oldroot->event);
    SWAP(gtree->root->event->data, oldroot->event->data);
    //SWAP(gtree->root->mark, oldroot->mark);

    gtree->root = oldroot;

    return ret;
  }

  return 0;
}

void gtree_fini(int msa_count)
{
  int i;

  /* free all module memory allocations */

  free(sortbuffer);
  for (i = 0; i < msa_count; ++i)
    free(travbuffer[i]);
  free(travbuffer);
}

static long propose_spr(locus_t * locus,
                        gtree_t * gtree,
                        stree_t * stree,
                        int msa_index)
{
  unsigned int i,j,k,m,n;
  unsigned int source_count, target_count;
  long accepted = 0;
  gnode_t * curnode;
  gnode_t * sibling;
  gnode_t * father;
  gnode_t * p;
  double minage,maxage,tnew;
  snode_t * pop;


  /*          

                                 *
                                / \
                       father  *   \
                              / \
                    curnode  /   \
                            *     *  sibling
                           / \   / \
                          /   \



  */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    curnode = gtree->nodes[i];
    if (curnode == gtree->root) continue;

    sibling = (curnode->parent->left == curnode) ? 
                curnode->parent->right : curnode->parent->left;
    father  = curnode->parent;
    
    assert(curnode->parent);

    /* find youngest population with subtree lineages more than current node subtree lineages */
    for (pop = curnode->pop; pop->gene_leaves[msa_index] <= curnode->leaves; pop = pop->parent)
      if (!pop->parent)
        break;

    /* TODO: Set age limits. 999 is set for backwards compatibility with bpp */
    minage = MAX(curnode->time, pop->tau);
    maxage = 999;

    tnew = father->time + opt_finetune_gtspr*legacy_rnd_symmetrical();
    tnew = reflect(tnew,minage,maxage);

    for (pop = curnode->pop; pop->parent; pop = pop->parent)
      if (pop->parent->tau > tnew) break;

    snode_t * pop_target = pop;

    /* identify target branches on which we can attach the pruned tree */
    /* TODO: We process the root node first to keep backwards compatibility with
       old bpp results */

    n = pop_target->node_index;
    target_count = 0;
    if (tnew >= gtree->root->time)
    {
      travbuffer[msa_index][target_count++] = gtree->root;
    }
    else
    {
      for (j = 0; j < gtree->tip_count + gtree->inner_count; ++j)
      {
        p = gtree->nodes[j];
        m = p->pop->node_index;
        if (p != curnode && p != gtree->root && p->time <= tnew &&
            p->parent->time > tnew && stree->pptable[m][n])
          travbuffer[msa_index][target_count++] = (p == father) ? sibling : p;
      }
    }

    source_count = 1;
    if (father != gtree->root)
    {
      n = father->pop->node_index;
      for (j = 0; j < gtree->tip_count + gtree->inner_count; ++j)
      {
        p = gtree->nodes[j];
        m = p->pop->node_index;
        if (p != curnode && p != gtree->root && p != sibling && p != father &&
            p->time <= father->time && p->parent->time > father->time &&
            stree->pptable[m][n])
          source_count++;
      }
    }

    assert(target_count);
    assert(source_count);
    

    /* randomly select a target node */
    gnode_t * target = travbuffer[msa_index][(int)(target_count*legacy_rndu())];

    /* regraft subtree */

    snode_t * oldpop = father->pop;
    double oldage = father->time;
    father->time = tnew;
    if (father->pop != pop_target)
    {
      /* TODO: update coalescent events */

      /* remove current gene node from the list of coalescent events of its old
         population */

      unlink_event(father,msa_index);

      /* decrease the number of coalescent events for the current population */
      father->pop->event_count[msa_index]--;
        
      /* change population for the current gene tree node */
      father->pop = pop_target;

      /* now add the coalescent event to the new population, at the end */
      dlist_item_append(father->pop->event[msa_index],father->event);

      father->pop->event_count[msa_index]++;

      /* increase or decrease the number of incoming lineages to all populations
         in the path from old population to the new population, depending on the
         case  */
      if (tnew > oldage)
      {
        /* increase the number of incoming lineages to all populations in the
           path from old population (excluding) to the new population */
        for (pop = oldpop; pop != father->pop; pop = pop->parent)
          pop->parent->seqin_count[msa_index]++;
      }
      else
      {
        /* decrease the number of incoming lineages to all populations in the
           path from new population (excluding) to the old population */
        for (pop = father->pop; pop != oldpop; pop = pop->parent)
          pop->parent->seqin_count[msa_index]--;
      }
    }
    
    int spr_required = (target != sibling && target != father);

    /* if the following holds we need to change tree topology */
    int root_changed = 0;
    if (spr_required)
      root_changed = perform_spr(gtree,curnode,target);

    if (root_changed)
      father = curnode->parent;

    /* recompute logpr */
    double logpr = gtree->logpr;
    if (oldpop == father->pop)
    {
      logpr -= father->pop->logpr_contrib[msa_index];
      logpr += gtree_update_logprob_contrib(father->pop,msa_index);
    }
    else
    {
      snode_t * start;
      snode_t * end;

      if (tnew > oldage)
      {
        start = oldpop;
        end   = father->pop->parent;
      }
      else
      {
        start = father->pop;
        end   = oldpop->parent;
      }
      for (pop = start; pop != end; pop = pop->parent)
      {
        logpr -= pop->logpr_contrib[msa_index];
        logpr += gtree_update_logprob_contrib(pop,msa_index);
      }
    }

    k = 0;
    travbuffer[msa_index][k++] = father->left;
    travbuffer[msa_index][k++] = father->right;
    if (father->parent)
      travbuffer[msa_index][k++] = father;
    if (spr_required)
      travbuffer[msa_index][k++] = sibling;
    locus_update_matrices_jc69(locus,travbuffer[msa_index],k);

    /* locate all nodes  whose CLV need to be updated */
    k = 0;
    if (!spr_required)
    {
      /* fill traversal buffer with root-path starting from father */
      gnode_t * temp;
      for (k=0, temp = father; temp; temp = temp->parent)
      {
        travbuffer[msa_index][k++] = temp;

        /* swap clv index to compute partials in a new location. This is useful
           when the proposal gets rejected, as we only have swap clv indices */
        temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
      }
    }
    else
    {
      /* if an SPR was done, we have two root-paths; one starting from sibling's
         parent and one from father */

      gnode_t * temp;

      /* mark the root-path starting from father */
      for (temp = father; temp; temp = temp->parent)
        temp->mark |= FLAG_MISC;    /* set FLAG_MISC */

      /* fill traversal buffer with nodes on the root-path starting from
         sibling's parent and stop when the lowest common ancestor is found */

      /* TODO: Should this also check for temp != NULL, in case father was root ? */
      for (temp=sibling->parent; !(temp->mark & FLAG_MISC); temp=temp->parent)
      {
        travbuffer[msa_index][k++] = temp;

        /* swap clv index to compute partials in a new location. This is useful
           when the proposal gets rejected, as we only have swap clv indices */
        temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
      }

      /* now fill the remaining traversal buffer with the root-path starting
         from father and reset markings */
      for (temp = father; temp; temp = temp->parent)
      {
        temp->mark &= ~FLAG_MISC;   /* unset FLAG_MISC */
        travbuffer[msa_index][k++] = temp;

        /* swap clv index to compute partials in a new location. This is useful
           when the proposal gets rejected, as we only have swap clv indices */
        temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
      }
    }

    /* update partials */
    locus_update_partials(locus,travbuffer[msa_index],k);

    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    double logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

    /* acceptance ratio */
    double acceptance = log((double)target_count / source_count) + logpr - gtree->logpr + logl - gtree->logl;

    if (opt_debug)
      printf("[Debug] (spr) lnacceptance = %f\n", acceptance);

    if (acceptance >= 0 || legacy_rndu() < exp(acceptance))
    {
      gtree->logpr = logpr;
      gtree->logl = logl;

      accepted++;
    }
    else
    {
      /* rejected */

      /* need to reset clv indices to point to the old clv buffer */
      for (j = 0; j < k; ++j)
      {
        gnode_t * temp = travbuffer[msa_index][j];
        temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
      }
      
      /* now reset branch lengths and pmatrices */

      /* TODO: Now this is correct, but be very careful about the exact steps of
         placing nodes in the traversal buffer. First place the fathers two
         children, then do the SPR to restore topology, then add father's parent
         and then if an SPR was required, add the sibling */

      father->time = oldage;
      k = 0;

      travbuffer[msa_index][k++] = father->left;
      travbuffer[msa_index][k++] = father->right;  /* target or sibling */

      if (spr_required)
      {
        /* if root_changed == 2 it means that the old sibling is now the root,
           and because of the old bpp compatibility code, the node the
           represents the root is *not* the node that used to be the old
           sibling. Therefore, we must pass the root as the target of the SPR
           (and not the node pointing to the old sibling). In the other cases,
           the target is just the old sibling */
        if (root_changed == 2)
          root_changed = perform_spr(gtree,curnode,gtree->root);
        else
          root_changed = perform_spr(gtree,curnode,sibling);
      }

      if (root_changed)
        father = curnode->parent;

      if (father->parent)
        travbuffer[msa_index][k++] = father;

      if (spr_required)
      {
        sibling = (curnode->parent->left == curnode) ? 
                    curnode->parent->right : curnode->parent->left;
        travbuffer[msa_index][k++] = sibling;
      }

      locus_update_matrices_jc69(locus,travbuffer[msa_index],k);

      if (father->pop == oldpop)
      {
        father->pop->logpr_contrib[msa_index] = father->pop->old_logpr_contrib[msa_index];
      }
      else
      {
        /* remove current gene node from the list of coalescent events of its old
           population */

        unlink_event(father,msa_index);

        /* decrease the number of coalescent events for the current population */
        father->pop->event_count[msa_index]--;
          
        /* change population for the current gene tree node */
        SWAP(father->pop,oldpop);

        /* now add the coalescent event back to the old population, at the end */
        dlist_item_append(father->pop->event[msa_index],father->event);

        father->pop->event_count[msa_index]++;

        /* increase or decrease the number of incoming lineages to all populations in the path
        from old population to the new population, depending on the case  */
        if (oldage >= tnew)
        {
          /* TODO : Now this is covered here!!! */

          /* increase the number of incoming lineages to all populations in the
             path from old population (excluding) to the new population */
          for (pop=oldpop; pop != father->pop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]++;
        }
        else
        {
          /* decrease the number of incoming lineages to all populations in the
             path from new population (excluding) to the old population */
          for (pop = father->pop; pop != oldpop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]--;
        }

        /* now restore the old log gene tree probability contribution for each
           affected species tree node */
        
        snode_t * start = (tnew > oldage) ? father->pop :  oldpop;
        snode_t * end   = (tnew > oldage) ? oldpop->parent : father->pop->parent;

        for (pop = start; pop != end; pop = pop->parent)
          pop->logpr_contrib[msa_index] = pop->old_logpr_contrib[msa_index];
      }
    }
  }
  return accepted;
}

double gtree_propose_spr(locus_t ** locus, gtree_t ** gtree, stree_t * stree)
{
  unsigned int i;
  long proposal_count = 0;
  long accepted = 0;

  for (i = 0; i < stree->locus_count; ++i)
  {
    /* TODO: Fix this to account mcmc.moveinnode in original bpp */
    proposal_count += gtree[i]->edge_count;
    accepted += propose_spr(locus[i],gtree[i],stree,i);
  }

  if (!accepted)
    return 0;

  return ((double)accepted/proposal_count);
}
