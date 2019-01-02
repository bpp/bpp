/*
    Copyright (C) 2016-2019 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))
#define SWAP_PMAT_INDEX(e,i) (i) = (((e)+(i))%((e)<<1))
#define SWAP_SCALER_INDEX(n,i) (((n)+((i)-1))%(2*(n)-2))

#define GET_HINDEX(t,p) (((node_is_mirror((p)) ? \
                          (p)->node_index : (p)->hybrid->node_index)) - \
                        ((t)->tip_count+(t)->inner_count))

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

gnode_t * dbg_msci_y = NULL;
gnode_t * dbg_msci_a = NULL;
gnode_t * dbg_msci_s = NULL;
gnode_t * dbg_msci_t = NULL;

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

void gtree_reset_leaves(gnode_t * node)
{
  if (!node->left)
  {
    node->leaves = 1;
    return;
  }

  gtree_reset_leaves(node->left);
  gtree_reset_leaves(node->right);
  
  node->leaves = node->left->leaves + node->right->leaves;
}

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
    if (node->hpath)
      free(node->hpath);

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
  long size_alloced;
  assert(root != NULL);

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = (long)strlen(newick);
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
  long size_alloced;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = (long)strlen(newick);
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
    assert(tip_count);
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

  if (stree->tip_count == 1)
  {
    pop[0].seq_count = msa->count;
    pop[0].seq_indices = (int *)xcalloc((size_t)(msa->count),sizeof(int));
    pop[0].nodes = (gnode_t **)xcalloc((size_t)(msa->count),sizeof(gnode_t *));
    for (i = 0; i < msa->count; ++i)
      pop[0].seq_indices[i] = i;
    return;
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

    /* increment sequence count for current loci in population i */
    pop[i].seq_count++;
  }

  /* now repeat that for obtaining the sequence indices */
  int * counter = (int *)xcalloc(stree->tip_count, sizeof(int));
  for (i = 0; i < stree->tip_count; ++i)
  {
    size_t alloc_size = pop[i].seq_count;
    if (opt_migration)
      alloc_size = msa->count;

    pop[i].seq_indices = (int *)xcalloc(alloc_size,sizeof(int));
    pop[i].nodes = (gnode_t **)xcalloc(alloc_size,sizeof(gnode_t *));
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

static void network_bd_distribute_lineages(pop_t * pop, long mindex, long hindex)
{
  size_t new_seqcount;
  int * new_indices;
  gnode_t ** new_nodes;

  if (pop[mindex].seq_count)
  {
    /* allocate new arrays */
    new_seqcount = pop[hindex].seq_count + pop[mindex].seq_count;
    new_indices = (int *)xmalloc(new_seqcount * sizeof(int));
    new_nodes = (gnode_t **)xmalloc(new_seqcount * sizeof(gnode_t *));

    /* copy indices from hybrid and mirror node to new arrays */
    memcpy(new_indices,
           pop[hindex].seq_indices,
           pop[hindex].seq_count * sizeof(int));
    memcpy(new_indices+pop[hindex].seq_count,
           pop[mindex].seq_indices,
           pop[mindex].seq_count * sizeof(int));
    
    /* copy nodes from hybrid and mirror node to new arrays */
    memcpy(new_nodes,
           pop[hindex].nodes,
           pop[hindex].seq_count * sizeof(gnode_t *));
    memcpy(new_nodes+pop[hindex].seq_count,
           pop[mindex].nodes,
           pop[mindex].seq_count * sizeof(gnode_t *));
    
    /* free old arrays at hybrid node and set the new ones */
    free(pop[hindex].seq_indices);
    free(pop[hindex].nodes);
    pop[hindex].seq_indices = new_indices;
    pop[hindex].nodes       = new_nodes;
    pop[hindex].seq_count   = new_seqcount;
  }
}

static void replace_hybrid(stree_t * stree, pop_t * pop, unsigned int * count, snode_t * epoch)
{
  long i,j;

  /* this function does the same thing as replace(...) (see comments there)
     but adapted to hybridization/bidirectional introgression nodes on networks.
     The main difference is that such nodes do not necessarily have two children */
  assert(opt_network && epoch->hybrid);

  /* check if its the mirrored node (i.e. no children) */
  if (node_is_mirror(epoch))
  {
    /* first find the hybridization node in the pop list */
    for (j = 0; j < *count; ++j)
      if (pop[j].snode == epoch->hybrid)
        break;
    assert(j != *count);

    pop_t * hpop = pop+j;  /* population corresponding to hybrid. node */

    /* samples gene tree nodes (lineages) and distribute them to pop[i].snode
       and pop[i].snode->hybrid (ie left and right populations of hybridization)
       according to probability gamma */
    int * temp = (int *)xmalloc((size_t)(hpop->seq_count) * sizeof(int));


    /* hpop.seq_count is the number of lineages coming from the child of the
       hybridization node H. We now need to partition these lineages into two
       populations (H_left and H_right) by sampling.  */

    /* Sample lineages: Fill temp with 0 if lineage goes to epoch (mirror node)
       or 1 if epoch->hybrid (hybridization) */
    for (j = 0; j < hpop->seq_count; ++j)
      temp[j] = (legacy_rndu() <= hpop->snode->hgamma) ? 1 : 0;

    /* count how many lineages ended up in hybridization node */
    unsigned int hpop_seqcount = 0;
    for (j = 0; j < hpop->seq_count; ++j)
      hpop_seqcount += temp[j];

    gnode_t ** hnodes;
    gnode_t ** mnodes;
    int * hindices;
    int * mindices;
    unsigned int hnodes_count = 0;
    unsigned int mnodes_count = 0;

    /* allocate arrays holding lineages (gene tree nodes) for the two pops */
    hnodes = (gnode_t **)xmalloc((size_t)(hpop_seqcount)*sizeof(gnode_t *));
    mnodes = (gnode_t **)xmalloc((size_t)(hpop->seq_count - hpop_seqcount) *
                                 sizeof(gnode_t *));

    /* allocate arrays holding sequence indices for the two pops */
    hindices = (int *)xmalloc((size_t)(hpop_seqcount)*sizeof(int *));
    mindices = (int *)xmalloc((size_t)(hpop->seq_count-hpop_seqcount)*sizeof(int));

    /* fill hnodes and mnodes */
    unsigned int hindex = GET_HINDEX(stree,epoch);
    assert(hindex >= 0 && hindex < stree->hybrid_count);

    /* make sure epoch is a mirror node */
    if (node_is_bidirection(epoch))
    {
      assert(!epoch->left && !epoch->right);
    }
    else
    {
      assert(node_is_hybridization(epoch));
      assert(!epoch->left && !epoch->right);
    }

    for (j = 0; j < hpop->seq_count; ++j)
      if (temp[j])
      {
        hnodes[hnodes_count] = hpop->nodes[j];
        hindices[hnodes_count++] = hpop->seq_indices[j];

        hpop->nodes[j]->hpath[hindex] = BPP_HPATH_LEFT;
      }
      else
      {
        mnodes[mnodes_count] = hpop->nodes[j];
        mindices[mnodes_count++] = hpop->seq_indices[j];

        hpop->nodes[j]->hpath[hindex] = BPP_HPATH_RIGHT;
        if (node_is_bidirection(epoch))
        {
          unsigned int hindex2 = GET_HINDEX(stree,epoch->parent);
          hpop->nodes[j]->hpath[hindex2] = BPP_HPATH_LEFT;
        }
      }
    assert(hnodes_count == hpop_seqcount);

    /* update the population corresponding to the hybridization node */
    free(hpop->nodes);
    free(hpop->seq_indices);
    hpop->nodes = hnodes;
    hpop->seq_indices = hindices;
    hpop->seq_count = hnodes_count;

    /* update the population corresponding to the hybridization node */
    pop[*count].snode = epoch;
    pop[*count].nodes = mnodes;
    pop[*count].seq_indices = mindices;
    pop[*count].seq_count = mnodes_count;
    
    if (node_is_bidirection(epoch))
    {
      long h1_index;
      long h2_index;
      long m1_index;
      long m2_index;
      assert(node_is_mirror(epoch->parent->hybrid));

      /* first find the other hybridization mirror node in the pop list */
      for (j = 0; j < *count; ++j)
        if (pop[j].snode == epoch->parent->hybrid)
          break;

      /* other side of introgression found */
      if (j < *count)
      {
        m1_index = j;
        m2_index = *count;

        /* now find the two hybridization nodes */

        for (j = 0; j < *count; ++j)
          if (pop[j].snode == epoch->parent)
            break;
        assert(j < *count);
        h1_index = j;

        for (j = 0; j < *count; ++j)
          if (pop[j].snode == epoch->hybrid)
            break;
        assert(j < *count);
        h2_index = j;

        /* move lineages stored in mirror nodes to the opposite hybrid nodes, i.e.
           
           m1 -> h2 and m2 -> h1;
        */

        network_bd_distribute_lineages(pop, m1_index, h2_index);
        network_bd_distribute_lineages(pop, m2_index, h1_index);

        if (m1_index > m2_index) SWAP(m1_index,m2_index);
        if (h1_index > h2_index) SWAP(h1_index,h2_index);

        assert(h1_index < h2_index && h1_index < m1_index && h1_index < m2_index);
        assert(m1_index < m2_index);
        assert(m2_index == *count);

        free(pop[m1_index].seq_indices);
        free(pop[m1_index].nodes);
        free(pop[m2_index].seq_indices);
        free(pop[m2_index].nodes);
        pop[m1_index].seq_indices = pop[m2_index].seq_indices = NULL;
        pop[m1_index].nodes = pop[m2_index].nodes = NULL;
        pop[m1_index].seq_count = pop[m2_index].seq_count = 0;
        pop[m1_index].snode = pop[m2_index].snode = NULL;
        if (h2_index > m1_index)
        {
          pop[m1_index].snode = pop[h2_index].snode;
          pop[m1_index].nodes = pop[h2_index].nodes;
          pop[m1_index].seq_indices = pop[h2_index].seq_indices;
          pop[m1_index].seq_count = pop[h2_index].seq_count;
        }

        assert(*count > 1);
        *count = *count - 2;

      }
    }

    /* we increase by two because on the next iteration in gtree_simulate, count
       will be decreased and as such we will have +1 new items (this node) */
    *count = *count + 2; 
    free(temp);

    return;
  }

  /* find and replace left descendant of current epoch with epoch */
  assert(!node_is_mirror(epoch));
  if (node_is_bidirection(epoch))
    assert(epoch->left && epoch->right);
  else
    assert(epoch->left && !epoch->right);

  for (i = 0; i < *count; ++i)
  {
    if (pop[i].snode == epoch->left)
      break;
  }
  assert(i != *count);

  /* replace population i with new population */
  pop[i].snode = epoch;

  *count = *count + 1;
}

static void replace(pop_t * pop, int count, snode_t * epoch)
{
  int i,j;

  /* delete the two children of epoch from the pop list, by replacing the child
     with the smaller index with epoch, and deleting the second child. It also
     merges the sequences (lineages) of the two children into the new population
  */
  assert(!epoch->hybrid);

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
  pop[i].seq_indices = NULL;
  pop[i].nodes = NULL;

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

static void fill_hybrid_seqin_counts(stree_t * stree, gtree_t * gtree, int msa_index)
{
  unsigned int i,j;

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * mnode = stree->nodes[stree->tip_count+stree->inner_count+i];
    snode_t * hnode = mnode->hybrid;

    /* hybridization */
    if (node_is_hybridization(mnode))
    {
      assert(!node_is_mirror(hnode));
      for (j = 0; j < gtree->tip_count; ++j)
      {
        if (gtree->nodes[j]->hpath[i] == BPP_HPATH_LEFT)
          hnode->seqin_count[msa_index]++;
        else if (gtree->nodes[j]->hpath[i] == BPP_HPATH_RIGHT)
          mnode->seqin_count[msa_index]++;
      }

      for (j = gtree->tip_count; j < gtree->tip_count+gtree->inner_count; ++j)
      {
        /* TODO: The NONE/NONE check is probably unnecessary, change to assertion */
        gnode_t * x = gtree->nodes[j];
        if (x->left->hpath[i] == BPP_HPATH_NONE &&
            x->right->hpath[i] == BPP_HPATH_NONE)
        {
          if (x->hpath[i] == BPP_HPATH_LEFT)
            hnode->seqin_count[msa_index]++;
          else if (x->hpath[i] == BPP_HPATH_RIGHT)
            mnode->seqin_count[msa_index]++;
        }
      }
    }
    else     /* bidirection */
    {
      assert(node_is_bidirection(mnode));
      assert(!node_is_mirror(hnode));
      unsigned int hindex2 = GET_HINDEX(stree,mnode->parent);

      for (j = 0; j < gtree->tip_count; ++j)
      {
        if (gtree->nodes[j]->hpath[i] == BPP_HPATH_LEFT)
        {
          if (gtree->nodes[j]->hpath[hindex2] == BPP_HPATH_NONE)
            hnode->seqin_count[msa_index]++;
        }
        else if (gtree->nodes[j]->hpath[i] == BPP_HPATH_RIGHT)
        {
          mnode->seqin_count[msa_index]++;
          mnode->parent->seqin_count[msa_index]++;
        }
      }

      for (j = gtree->tip_count; j < gtree->tip_count+gtree->inner_count; ++j)
      {
        /* TODO: The NONE/NONE check is probably unnecessary, change to assertion */
        gnode_t * x = gtree->nodes[j];
        if (x->left->hpath[i] == BPP_HPATH_NONE &&
            x->right->hpath[i] == BPP_HPATH_NONE)
        {
          if (x->hpath[i] == BPP_HPATH_LEFT)
          {
            if (x->hpath[hindex2] == BPP_HPATH_NONE)
              hnode->seqin_count[msa_index]++;
          }
          else if (x->hpath[i] == BPP_HPATH_RIGHT)
          {
            mnode->seqin_count[msa_index]++;
            mnode->parent->seqin_count[msa_index]++;
          }
        }
      }
    } /* end of else bidirection */
  } /* end of for */
}

static void fill_seqin_counts_recursive(snode_t * node, int msa_index)
{
  if (!node->left)
    return;

  if (node->left)
    fill_seqin_counts_recursive(node->left,msa_index);
  if (node->right)
    fill_seqin_counts_recursive(node->right,msa_index);

  snode_t * lnode = node->left;
  snode_t * rnode = node->right;

  if (opt_network)
  {
    if (node->hybrid)
    {
      if (node_is_hybridization(node))
      {

        /* get child node in species tree. Note that it should always be the left
           child of the hybridization node */
        /* TODO: Simplify the below IFs, the current status is only for testing
           correctness */
        snode_t * child = NULL;
        if (node->left) 
        {
          assert(!child);
          child = node->left;
        }
        if (node->right) 
        {
          assert(!child);
          child = node->right;
        }
        if (node->hybrid->left)
        {
          assert(!child);
          child = node->hybrid->left;
        }
        if (node->hybrid->right)
        {
          assert(!child);
          child = node->hybrid->right;
        }

        assert((child->seqin_count[msa_index] - child->event_count[msa_index]) ==
               (node->seqin_count[msa_index] + node->hybrid->seqin_count[msa_index]));
      }
      else
      {
        assert(node_is_bidirection(node));

        if (node_is_mirror(node))
        {
        }
        else
        {
          assert(node->left && node->right);

          /* cannot do the checks here as the opposite node might have not yet been processed */
        }

      }
      return;
    }


    node->seqin_count[msa_index] = 0;
    if (lnode)
      node->seqin_count[msa_index] += lnode->seqin_count[msa_index] -
                                      lnode->event_count[msa_index];
    if (rnode)
      node->seqin_count[msa_index] += rnode->seqin_count[msa_index] -
                                      rnode->event_count[msa_index];
  }
  else
  {
    /* if no networks then this is valid */
    node->seqin_count[msa_index] = lnode->seqin_count[msa_index] +
                                   rnode->seqin_count[msa_index] -
                                   lnode->event_count[msa_index] -
                                   rnode->event_count[msa_index];
  }
}

void fill_seqin_counts(stree_t * stree, gtree_t * gtree, int msa_index)
{
  if (opt_network) 
    fill_hybrid_seqin_counts(stree, gtree, msa_index);

  fill_seqin_counts_recursive(stree->root,msa_index);

  /* TODO: This is only an optional check, delete */
  if (opt_network)
  {
    long i;
    for (i = 0; i < stree->hybrid_count; ++i)
    {
      snode_t * hnode;
      snode_t * mnode;

      mnode = stree->nodes[stree->tip_count+stree->inner_count+i];
      hnode = mnode->hybrid;

      if (node_is_bidirection(hnode))
      {
        assert(node_is_mirror(mnode));
        assert(!node_is_mirror(hnode));

        assert(mnode->parent && mnode->parent->hybrid);
        assert(mnode->seqin_count[msa_index] <= mnode->parent->seqin_count[msa_index]);

        assert(hnode->right && hnode->right->hybrid);
        assert(hnode->seqin_count[msa_index] >= hnode->right->seqin_count[msa_index]);

        assert(hnode->seqin_count[msa_index] + hnode->right->hybrid->seqin_count[msa_index] ==
               hnode->left->seqin_count[msa_index] - hnode->left->event_count[msa_index] +
               hnode->right->hybrid->left->seqin_count[msa_index] - hnode->right->hybrid->left->event_count[msa_index]);
      }
    }
  }
}

static int cb_trav_full(snode_t * x)
{
  if (!x->left)
    return 0;

  return 1;
}

static void epoch_reorder(snode_t ** epoch,
                          unsigned int epoch_count,
                          snode_t * hnode)
{
  unsigned int i;
  unsigned int hindex;

  /* first find the index of hnode in epoch */
  for (hindex = 0; hindex < epoch_count; ++hindex)
    if (epoch[hindex] == hnode)
      break;

  assert(hindex < epoch_count);

  if (node_is_bidirection(hnode))       /* bidirections */
  {
    /* we need to have mirror nodes right after their hybrid counterparts */
    for (i = 0; i < hindex; ++i)
    {
      if (epoch[i] == epoch[hindex]->hybrid)
      {
        SWAP(epoch[i],epoch[hindex]);
        break;
      }
    }
  }
  else          /* hybridization */
  {
    assert(!node_is_mirror(hnode));

    /* make sure parents are after hybridization event */
    if (!epoch[hindex]->parent->htau ||
        epoch[hindex]->parent->tau == epoch[hindex]->tau)
    {
      for (i = 0; i < hindex; ++i)
      {
        if (epoch[i] == epoch[hindex]->parent)
        {
          SWAP(epoch[i],epoch[hindex]);
          hindex = i;
          break;
        }
      }
    }

    if (!epoch[hindex]->hybrid->parent->htau ||
        epoch[hindex]->hybrid->parent->tau == epoch[hindex]->tau)
    {
      for (i = 0; i < hindex; ++i)
      {
        if (epoch[i] == epoch[hindex]->hybrid->parent)
        {
          SWAP(epoch[i],epoch[hindex]);
          hindex = i;
          break;
        }
      }
    }

    /* make sure its mirror node is directly after it */
    unsigned int mindex;
    for (mindex = 0; mindex < epoch_count; ++mindex)
      if (epoch[mindex] == hnode->hybrid)
        break;
    assert(mindex < epoch_count);
    assert(epoch[mindex]->tau == epoch[hindex]->tau);

    if (mindex < hindex)
    {
      assert(mindex == hindex-1);
      SWAP(epoch[hindex],epoch[mindex]);
    }
    else if (mindex > hindex+1)
    {
      assert(mindex - hindex <= 3);
      while (mindex != hindex+1)
      {
        SWAP(epoch[mindex],epoch[mindex-1]);
        mindex--;
      }
    }
  }
}

gtree_t * gtree_simulate(stree_t * stree, msa_t * msa, int msa_index)
{
  int lineage_count = 0;
  int scaler_index = 0;
  unsigned int i,j,k;
  unsigned int epoch_count;
  double t, tmax, csum, msum=0;
  double * ci;
  double * migrate = NULL;
  pop_t * pop;
  snode_t ** epoch;
  gnode_t * inner = NULL;

  if (opt_migration)
    migrate = (double *)xmalloc((size_t)stree->tip_count * sizeof(double));

  /* get a list of inner nodes (epochs) */
  epoch_count = stree->inner_count;
  if (opt_network)
    epoch_count = stree->inner_count + stree->hybrid_count;

  epoch = (snode_t **)xmalloc((size_t)epoch_count*sizeof(snode_t *));
  memcpy(epoch,
         stree->nodes + stree->tip_count,
         (epoch_count) * sizeof(snode_t *));

  /* sort epochs in ascending order of speciation times */
  if (opt_network)
  {
    memcpy(epoch,
           stree->nodes + stree->tip_count,
           (size_t)(epoch_count) * sizeof(snode_t *));
  }
  else
  {
    stree_traverse(stree->root,
                   TREE_TRAVERSE_POSTORDER,
                   cb_trav_full,
                   epoch,
                   &epoch_count);

    assert(epoch_count == stree->inner_count);
  }

  snode_t ** sortptr = epoch;

  /* move zeroes to beginning and store their number in j */
  for (j = 0, i = 0; i < epoch_count; ++i)
  {
    if (epoch[i]->tau == 0)
    {
      if (i != j)
        SWAP(epoch[i], epoch[j]);
      
      ++j;
    }
  }

  /* sort nonzero entries */
  sortptr += j;
  qsort(&(sortptr[0]), epoch_count-j, sizeof(snode_t *), cb_cmp_spectime);

  if (opt_network)
  {
    /* if we have a hybridization event, composed of nodes H,S,T, e.g.

                                  S   T  
                                   \ /
                                    * H
                                    |

       we need to make sure that in the case either S or T have the same tau as
       H, then they must appear *after* H in the 'epoch' array (see Ziheng's diagram)
    */

    for (i = 0; i < stree->inner_count; ++i)
    {
      snode_t * x = stree->nodes[stree->tip_count + i];
      if (x->hybrid)
        epoch_reorder(epoch,epoch_count,x);
    }
  }
  if (opt_debug)
  {
    for (i = 0; i < epoch_count; ++i)
      printf("[Debug]: Epoch %d (%s) time - %f\n",
             i, epoch[i]->label, epoch[i]->tau);
  }

  /* create one hash table of species and one for sequence->species mappings */
  pop = (pop_t *)xcalloc((size_t)(stree->tip_count+stree->hybrid_count),
                         sizeof(pop_t));
  fill_pop(pop,stree,msa,msa_index);

  for (i = 0; i < stree->tip_count; ++i)
    stree->nodes[i]->seqin_count[msa_index] = pop[i].seq_count;

  if (!opt_est_theta)
  {
    stree->notheta_logpr = 0;
    if (opt_est_heredity)
      stree->notheta_logpr += stree->notheta_hfactor;

    stree->notheta_old_logpr = 0;

    stree->notheta_sfactor = 0;
    for (i = 0; i < opt_locus_count; ++i)
    {
      long seqs = 0;
      for (j = 0; j < stree->tip_count; ++j)
        seqs += stree->nodes[j]->seqin_count[i];

      stree->notheta_logpr   += (seqs-1)*0.6931471805599453;
      stree->notheta_sfactor += (seqs-1)*0.6931471805599453;
    }
  }
  /* start at present time */
  t = 0;

  lineage_count = msa->count;

  /* allocate space for storing coalescent rates for each population */
  ci = (double *)xmalloc((size_t)(stree->tip_count + stree->hybrid_count) *
                         sizeof(double));

  /* current epoch index */
  unsigned int e = 0;

  /* create a list of tip nodes for the target gene tree */
  gnode_t ** gtips = (gnode_t **)xcalloc((size_t)(msa->count),
                                         sizeof(gnode_t *));
  for (i = 0; i < (unsigned int)(msa->count); ++i)
  {
    gtips[i] = (gnode_t *)xcalloc(1,sizeof(gnode_t));
    gtips[i]->pmatrix_index = i;
    gtips[i]->scaler_index = PLL_SCALE_BUFFER_NONE;
    gtips[i]->leaves = 1;

    if (opt_network)
    {
      /* TODO: Allocating hpath for tips is unnecessary, as they never pass
         through a hybridization event. However, if we don't allocate hpath
         we will need to check that a node is inner when accessing hpath */
      gtips[i]->hpath = (int *)xmalloc((size_t)(stree->hybrid_count) *
                                       sizeof(int));
      for (j = 0; j < stree->hybrid_count; ++j)
        gtips[i]->hpath[j] = BPP_HPATH_NONE;
    }
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
    if (pop_count == 1 && pop[0].snode == stree->root)
      tmax = -1;
    else
      tmax = epoch[e]->tau;
      
    while (1)
    {
      if (!tmax) break;

      /* calculate poisson rates: ci[j] is coalescent rate for population j */
      for (j=0, csum=0; j < pop_count; ++j)
      {
        k = pop[j].seq_count;

        if (k >= 2)
        {
          ci[j] = k*(k-1)/pop[j].snode->theta;
          csum += ci[j];
        }
        else
          ci[j] = 0;
      }

      if (opt_migration)
      {
        assert(!opt_network);
        msum = 0;
        long matrix_span = stree->tip_count + stree->inner_count;
        for (j = 0; j < pop_count; ++j)
        {
          for (k = 0, migrate[j] = 0; k < pop_count; ++k)
          {
            long mindexk = pop[k].snode->node_index;
            long mindexj = pop[j].snode->node_index;
            migrate[j] += pop[j].seq_count *
                          opt_migration_matrix[mindexk * matrix_span + mindexj] /
                          pop[j].snode->theta*4;
          }
          msum += migrate[j];
        }
      }

      if (csum+msum < 1e-300)
        break;

      /* generate random waiting time from exponential distribution */
      t += legacy_rndexp(1/(csum+msum));

      
      /* if the generated time is larger than the current epoch, and we are not
         yet at the root of the species tree, then break and, subsequently, 
         merge the lineages of the two populations into the current epoch */
      if (t > tmax && (pop_count != 1 || pop[0].snode != stree->root)) break;

      if (opt_debug)
        fprintf(stdout, "[Debug]: Coalescent waiting time: %f\n", t);

      /* TODO: Implement migration routine */

      /* select an available population at random using the poisson rates as
         weights */
      double r = legacy_rndu()*(csum+msum);
      if (r < csum)
      {
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

        if (opt_scaling)
          inner->scaler_index = scaler_index++;
        else
          inner->scaler_index = PLL_SCALE_BUFFER_NONE;

        inner->pmatrix_index = clv_index;
        inner->left->parent = inner;
        inner->right->parent = inner;
        inner->time = t;
        inner->pop = pop[j].snode;
        inner->leaves = inner->left->leaves + inner->right->leaves;
        clv_index++;

        if (opt_network)
        {
          inner->hpath = (int*)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));
          for (i = 0; i < stree->hybrid_count; ++i)
            inner->hpath[i] = BPP_HPATH_NONE;
        }

        pop[j].snode->event_count[msa_index]++;
        dlist_item_t * dlitem = dlist_append(pop[j].snode->event[msa_index],inner);
        inner->event = dlitem;
        if (!opt_est_theta)
          pop[j].snode->event_count_sum++;

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
      else
      {
        r -= csum;
        double tmp = 0;
        for (j = 0; j < pop_count; ++j)
        {
          tmp += migrate[j];
          if (r < tmp) break;
        }

        if (j == pop_count)
          fatal("The impossible just happened! Report gree_simulate()");
        tmp -= migrate[j];

        long mindexk = pop[0].snode->node_index;
        long mindexj = pop[j].snode->node_index;
        long matrix_span = stree->tip_count + stree->inner_count;
        for (k = 0; k < pop_count-1; ++k)
        {
          mindexk = pop[k].snode->node_index;
          tmp += pop[j].seq_count * 
                 opt_migration_matrix[mindexk * matrix_span + mindexj] /
                 pop[j].snode->theta*4;
          if (r < tmp)
            break;
        }
        mindexk = pop[k].snode->node_index;

        opt_migration_events[mindexk * matrix_span + mindexj] += 1;

        /* i is a migrant from population k to j */
        i = (long)(pop[j].seq_count * legacy_rndu());

        /* shift up lineages in pop j */
        pop[k].seq_indices[pop[k].seq_count] = pop[j].seq_indices[i];
        pop[k].nodes[pop[k].seq_count++] = pop[j].nodes[i];
        if (i != --pop[j].seq_count)
        {
          pop[j].seq_indices[i] = pop[j].seq_indices[pop[j].seq_count];
          pop[j].nodes[i] = pop[j].nodes[pop[j].seq_count];
        }
      }
    }

    t = tmax;

    if (lineage_count == 1 || (pop_count == 1 && pop[0].snode == stree->root)) break;

    /* place current epoch in the list of populations, remove its two children
       and add up the lineages of the two descendants */
    if (opt_network && epoch[e]->hybrid)
      replace_hybrid(stree,pop,&pop_count,epoch[e]);
    else
      replace(pop,pop_count,epoch[e]);
    
    if (e != epoch_count-1)
    {
      ++e;
    }
  }

  if (migrate)
    free(migrate);

  for(i = 0; i < pop_count; ++i)
  {
    free(pop[i].seq_indices);
    free(pop[i].nodes);
  }

  /* wrap the generated tree structure (made up of linked nodes) into gtree_t */
  gtree_t * gtree = gtree_wraptree(inner, (unsigned int)(msa->count));

  /* set path flags for gene tree root lineage if root coalesces before root
     population */
  if (opt_network)
  {
    snode_t * father = gtree->root->pop->parent;
    while (father)
    {
      if (father->hybrid)
      {
        long bidir = node_is_bidirection(father);

        unsigned int hindex = GET_HINDEX(stree, father);
        assert(hindex >= 0 && hindex < stree->hybrid_count);

        if (legacy_rndu() <= father->hgamma)
           gtree->root->hpath[hindex] = BPP_HPATH_LEFT;
        else
        {
          gtree->root->hpath[hindex] = BPP_HPATH_RIGHT;
          if (bidir)
          {
            assert(!node_is_mirror(father));
            unsigned int hindex2 = GET_HINDEX(stree,father->right);
            gtree->root->hpath[hindex2] = BPP_HPATH_LEFT;
            father = father->right->hybrid;
          }
          else  /* hybridization */
          {
            father = father->hybrid;
          }
        }        
      }
      father = father->parent;
    }
  }

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
  fill_seqin_counts(stree,gtree,msa_index);

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

/* TODO: The following is a very inefficient way of computing gene leave counts
   for hybridization nodes. For each hybridization node H the method does:
    1. Go through all gene tree tips of locus msa_index.
      a. For each tip visit its parent and repeat until we reach the root.
      b. If any of the visited nodes' hpath flag passes H then increase
         the gene leaves variable for node H  */
static void reset_hybrid_gene_leaves_count(stree_t * stree,
                                           gtree_t * gtree,
                                           int msa_index)
{
  unsigned int i,j;
  gnode_t * x;
  snode_t * mnode;
  snode_t * hnode;
  long bidir = 0;

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    mnode = stree->nodes[stree->tip_count+stree->inner_count+i];
    hnode = mnode->hybrid;
    assert(!node_is_mirror(hnode));

    hnode->gene_leaves[msa_index] = 0;
    mnode->gene_leaves[msa_index] = 0;
  }

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    mnode = stree->nodes[stree->tip_count+stree->inner_count+i];

    //assert(node_is_hybridization(mnode));
    bidir = 0;
    if (node_is_bidirection(mnode))
      bidir = 1;

    hnode = mnode->hybrid;
    assert(!node_is_mirror(hnode));


    for (j = 0; j < gtree->tip_count; ++j)
      for (x = gtree->nodes[j]; x; x = x->parent)
      {
        if (x->hpath[i] == BPP_HPATH_LEFT)
        {
          if (bidir)
          {
            unsigned int hindex2 = GET_HINDEX(stree,mnode->parent);
            if (x->hpath[hindex2] == BPP_HPATH_NONE)
              hnode->gene_leaves[msa_index]++;
          }
          else
            hnode->gene_leaves[msa_index]++;
          break;
        }
        if (x->hpath[i] == BPP_HPATH_RIGHT)
        {
          mnode->gene_leaves[msa_index]++;
          if (bidir)
          {
            assert(mnode->parent && mnode->parent->hybrid);
            mnode->parent->gene_leaves[msa_index]++;
          }
          break;
        }
      }
  }
}

static void reset_gene_leaves_count_recursive(snode_t * node, unsigned int locus_count)
{
  unsigned int j;

  if (!node->left && !node->right)
    return;

  if (node->left)
    reset_gene_leaves_count_recursive(node->left,locus_count);
  if (node->right)
    reset_gene_leaves_count_recursive(node->right,locus_count);

  if (opt_network && node->hybrid)
  {
    /* TODO: The below if block is only for testing correctness - it is not
       necessary and should be removed, but the 'return' statement must be
       kept */

    if (node_is_hybridization(node))
    {
      if (node_is_mirror(node))
        assert(!node->left && !node->right);
      else
        assert(node->left && !node->right);

      for (j = 0; j < locus_count; ++j)
        if (!node_is_mirror(node))
          assert((node->gene_leaves[j] + node->hybrid->gene_leaves[j]) ==
                 node->left->gene_leaves[j]);
    }
    else
    {
      /* bidirection */
      assert(node_is_bidirection(node));

      for (j = 0; j < locus_count; ++j)
        if (!node_is_mirror(node))
        {
          assert(node->left && node->right);

          /* TODO: we cannot perform this check as gene_leaves[j] from node->hybrid->left->hybrid
             may not be yet filled */
        }
    }

    return;
  }

  for (j = 0; j < locus_count; ++j)
    node->gene_leaves[j] = node->left->gene_leaves[j] +
                           node->right->gene_leaves[j];

}
void reset_gene_leaves_count(stree_t * stree, gtree_t ** gtree)
{
  unsigned int i,j;

  /* gtree is only necessary when opt_network */
  if (opt_network)
  {
    assert(gtree);
    for (i = 0; i < stree->locus_count; ++i)
      reset_hybrid_gene_leaves_count(stree,gtree[i],i);
  }

  /* gene leaves is the same as sequences coming in for tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
    for (j = 0; j < stree->locus_count; ++j)
      stree->nodes[i]->gene_leaves[j] = stree->nodes[i]->seqin_count[j];

  reset_gene_leaves_count_recursive(stree->root, stree->locus_count);

  /* TODO: IS THIS NECESSARY ? */
  if (stree->tip_count > 1)
  {
    for (j = 0; j < stree->locus_count; ++j)
      stree->root->gene_leaves[j] = stree->root->left->gene_leaves[j] +
                                    stree->root->right->gene_leaves[j];
  }
}

void gtree_alloc_internals(gtree_t ** gtree, long msa_count)
{
  long i;
  size_t minsize;

  /* allocate sort buffer */
  int max_count = 0;
  for (i = 0; i < msa_count; ++i)
    if (gtree[i]->tip_count > max_count)
      max_count = gtree[i]->tip_count;

  /* alloate buffer for sorting coalescent times plus two for the beginning
     and end of epoch */
  sortbuffer = (double *)xmalloc((size_t)(max_count+2) * sizeof(double));
  
  /* allocate traversal buffers */
  travbuffer = (gnode_t ***)xmalloc((size_t)msa_count * sizeof(gnode_t **));
  for (i = 0; i < msa_count; ++i)
  {
    /* in the gene tree SPR it is possible that this scenario happens:

                 *
                / \
       father  *   *
              / \ / \
             1  2 3  4

       father moving to 3 or 4, in whih case we need to update 4 CLVs
    */

    minsize = MAX(gtree[i]->inner_count,4);
    
    travbuffer[i] = (gnode_t **)xmalloc(minsize*sizeof(gnode_t *));
  }

}

void gtree_simulate_init(stree_t * stree, list_t * maplist)
{
  if (stree->tip_count == 1)
  {
    sht = NULL;
    mht = NULL;
  }
  else
  {
    sht = species_hash(stree);
    mht = maplist_hash(maplist,sht);
  }
}

void gtree_simulate_fini()
{
  if (sht)
    hashtable_destroy(sht,NULL);
  if (mht)
    hashtable_destroy(mht,cb_dealloc_pairlabel);
}
                        

gtree_t ** gtree_init(stree_t * stree,
                      msa_t ** msalist,
                      list_t * maplist,
                      int msa_count)
{
  int i;
  gtree_t ** gtree;

  assert(msa_count > 0);

  gtree = (gtree_t **)xmalloc((size_t)msa_count*sizeof(gtree_t *));

  /* create mapping hash tables */
  if (stree->tip_count == 1)
  {
    sht = NULL;
    mht = NULL;
  }
  else
  {
    sht = species_hash(stree);
    mht = maplist_hash(maplist,sht);
  }

  /* generate random starting gene trees for each alignment */
  printf("Generating gene trees....");
  for (i = 0; i < msa_count; ++i)
    gtree[i] = gtree_simulate(stree, msalist[i],i);
  printf(" Done\n");

  /* destroy the hash tables */
  if (stree->tip_count > 1)
  {
    hashtable_destroy(sht,NULL);
    hashtable_destroy(mht,cb_dealloc_pairlabel);
  }

  /* allocate static internal arrays sortbuffer and travbuffer */
  gtree_alloc_internals(gtree,msa_count);

  /* reset number of gene leaves associated with each species tree subtree */
  reset_gene_leaves_count(stree,gtree);

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

void logprob_revert_notheta(snode_t * snode, long msa_index)
{
  snode->t2h_sum -= snode->t2h[msa_index];
  snode->t2h[msa_index] = snode->old_t2h[msa_index];
  snode->t2h_sum += snode->t2h[msa_index];
  snode->notheta_logpr_contrib = snode->notheta_old_logpr_contrib;
}

double gtree_update_logprob_contrib(snode_t * snode,
                                    double heredity,
                                    long msa_index)
{
    unsigned int j,k,n;
    double logpr = 0;
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

    if (opt_network && snode->hybrid)
    {
      if (node_is_bidirection(snode) && !node_is_mirror(snode))
        logpr += (snode->seqin_count[msa_index] - snode->right->seqin_count[msa_index]) * log(snode->hgamma);
      else
        logpr += snode->seqin_count[msa_index] * log(snode->hgamma);
    }

    /* now distinguish between estimating theta and analytical computation */
    if (opt_est_theta)
    {
      if (snode->event_count[msa_index])
        logpr += snode->event_count[msa_index] * log(2.0/snode->theta);

      if (T2h)
        logpr -= T2h/snode->theta;

      /* TODO: Be careful about which functions update the logpr contribution 
         and which do not */
      snode->old_logpr_contrib[msa_index] = snode->logpr_contrib[msa_index];
      snode->logpr_contrib[msa_index] = logpr;
    }
    else
    {
      snode->old_t2h[msa_index] = snode->t2h[msa_index];

      snode->t2h[msa_index] = T2h;

      snode->t2h_sum -= snode->old_t2h[msa_index];
      snode->t2h_sum += snode->t2h[msa_index];
      

      if (snode->event_count_sum)
        logpr += opt_theta_alpha*log(opt_theta_beta) - lgamma(opt_theta_alpha) -
                 (opt_theta_alpha + snode->event_count_sum) *
                 log(opt_theta_beta + snode->t2h_sum) +
                 lgamma(opt_theta_alpha + snode->event_count_sum);
      else
        logpr -= opt_theta_alpha * log(1 + snode->t2h_sum / opt_theta_beta);

      /* TODO: this always updates the 'notheta_old_logpr_contrib'. Sometimes we
         do not want to this update because there could be multiple changes on
         the 'notheta_logpr_contrib' before deciding whether to accept or reject
         the proposal, e.g. by proposing new values on multiple loci. This leads
         to the problem that the 'notheta_old_logpr_contrib' is no longer the
         old value before any of the proposals started.  Currently this is fixed
         in the caller functions by storing the 'notheta_old_logpr_contrib' in
         some array allocated at the caller, but I should change this to only
         update the value through a flag passed to this function */
      snode->notheta_old_logpr_contrib = snode->notheta_logpr_contrib;
      snode->notheta_logpr_contrib = logpr;
    }

    return logpr;
}

double gtree_logprob(stree_t * stree, double heredity, long msa_index)
{
  unsigned int i;

  double logpr = 0;

  for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
    logpr += gtree_update_logprob_contrib(stree->nodes[i],heredity,msa_index);

  return logpr;
}

double reflect(double x, double a, double b)
{
  int side = 0;
  double n,excess = 0;
  const double EPSILON = 1e-200;

  /* returns a variable in range (a,b) by reflecting x back into the range */

  if (b-a < EPSILON)
    fatal("Internal error when proposing gene tree noge age");


  if (x < a)
  {
    excess = a - x;
    side = 0;
  }
  else if (x > b)
  {
    excess = x - b;
    side = 1;
  }

  if (excess)
  {
    double diff = b - a;

    n = floor(excess / diff);

    if (fmod(n,2.0) > 0.1)
      side = !side;

    excess -= n*diff;

    x = side ? b-excess : a+excess;
  }

   /* The following is to fix the problem of x landing on the boundary,
   due to small chances and rounding errors */
   while ((x - a < EPSILON) || (b - x < EPSILON))
      x = a + (b - a)*legacy_rndu();

  return x;
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

static void interchange_flags(stree_t * stree,
                              gnode_t * a,
                              gnode_t * b,
                              double tbnew,
                              double tbold,
                              snode_t * old_b_pop)
{
  snode_t * start;
  snode_t * end;
  snode_t * pop;

  if (old_b_pop == b->pop) return;

  if (tbnew < tbold)
  {
    /* old_b_pop *must* be an ancestor of b->pop */
    assert(stree->pptable[b->pop->node_index][old_b_pop->node_index]);

    start = b->pop;
    end = old_b_pop;

    pop = start;
    while (1)
    {
      pop = pop->parent;
      if (pop->hybrid)
      {
        assert(!node_is_mirror(pop));

        unsigned int hindex = GET_HINDEX(stree,pop);
        assert(hindex >= 0 && hindex < stree->hybrid_count);

        assert(a->hpath[hindex] != BPP_HPATH_NONE);
        assert(b->hpath[hindex] == BPP_HPATH_NONE);

        if (a->hpath[hindex] == BPP_HPATH_RIGHT)
          pop = pop->hybrid;

        b->hpath[hindex] = a->hpath[hindex];
        a->hpath[hindex] = BPP_HPATH_NONE;
      }
      if (pop == end) break;
    }
  }
  else
  {
    /* old_b_pop *must* be an ancestor of b->pop */
    assert(stree->pptable[old_b_pop->node_index][b->pop->node_index]);

    start = old_b_pop;
    end = b->pop;

    pop = start;
    while (1)
    {
      pop = pop->parent;
      if (pop->hybrid)
      {
        assert(!node_is_mirror(pop));

        unsigned int hindex = GET_HINDEX(stree,pop);
        assert(hindex >= 0 && hindex < stree->hybrid_count);

        assert(b->hpath[hindex] != BPP_HPATH_NONE);
        assert(a->hpath[hindex] == BPP_HPATH_NONE);

        if (b->hpath[hindex] == BPP_HPATH_RIGHT)
          pop = pop->hybrid;

        a->hpath[hindex] = b->hpath[hindex];
        b->hpath[hindex] = BPP_HPATH_NONE;
      }
      if (pop == end) break;
    }

  }
}


static void split_flags(stree_t * stree, gnode_t * a)
{
  long i;
  gnode_t * b     = a->parent;
  gnode_t * c     = b->parent;
  snode_t * start = b->pop;
  snode_t * end   = c ? c->pop : stree->root;
  snode_t * pop;

  /* splits the flags of a long branch a-c when b is placed in between:

     *  c              *  c
      \                 \
       \                 \
        \        ->       *  b
         \                 \
          \                 \
           * a               *  a

     IMPORTANT: The function expects that we already have edge a-b and b-c (b is
     already inserted)

     IMPORTANT: cleans flags for b beforehand
  */

  /* reset flags for edge b-c */
  for (i = 0; i < stree->hybrid_count; ++i)
    b->hpath[i] = BPP_HPATH_NONE;

  if (start == end) return;

  pop = start;
  while (1)
  {
    pop = pop->parent;
    if (pop->hybrid)
    {
      assert(!node_is_mirror(pop));

      unsigned hindex = GET_HINDEX(stree,pop);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      assert(a->hpath[hindex] != BPP_HPATH_NONE);

      if (a->hpath[hindex] == BPP_HPATH_RIGHT)
        pop = pop->hybrid;

      b->hpath[hindex] = a->hpath[hindex];
      a->hpath[hindex] = BPP_HPATH_NONE;
    }
    if (pop == end) break;
  }
}

static void join_flags(stree_t * stree, gnode_t * a, int * bpath, snode_t * old_b_pop)
{
  gnode_t * c     = a->parent;
  snode_t * start = old_b_pop;
  snode_t * end   = c ? c->pop : stree->root;
  snode_t * pop;

  /* joins the flags of two branches a-b and b-c when b is pruned:

     *  c              *  c
      \                 \
       \                 \
        * b      ->       \
         \                 \
          \                 \
           * a               *  a

     IMPORTANT: The function expects that we already have an edge a-c (b is
     already pruned off). The flags and populations of b from when b-c used
     to be an edge are given as bpath and old_b_pop.
  */

  if (start == end) return;

  assert(stree->pptable[start->node_index][end->node_index]);

  pop = start;
  while (1)
  {
    pop = pop->parent;
    if (pop->hybrid)
    {
      assert(!node_is_mirror(pop));

      unsigned int hindex = GET_HINDEX(stree,pop);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      assert(bpath[hindex] != BPP_HPATH_NONE);

      if (bpath[hindex] == BPP_HPATH_RIGHT)
        pop = pop->hybrid;

      a->hpath[hindex] = bpath[hindex];
    }
    if (pop == end) break;
  }
}

/* sample flags from pop(x) (exclusive) to pop(parent->x) (exclusive), ie
   sample flags for one branch with two ends fixed in the populations */
static double sample_hpath(stree_t * stree, gnode_t * x)
{
  long i;
  unsigned int hindex;
  long bidir = 0;
  snode_t * pop;
  snode_t * start = x->pop;
  snode_t * end   = x->parent ? x->parent->pop : stree->root;
  double contrib = 0;



  if (start == end) 
  {
    for (i = 0; i < stree->hybrid_count; ++i)
      x->hpath[i] = BPP_HPATH_NONE;
    return contrib;
  }

  /* array to keep track of which hybridization nodes were visited */
  int * visited = (int *)xcalloc((size_t)stree->hybrid_count,sizeof(int));

  assert(start->parent);

  /* sample flags for populations in between (excluding the two end points ) */
  //for (pop = start->parent; pop != end; pop = pop->parent)
  pop = start;
  while (1)
  {
    assert(pop && pop->parent);
    pop = pop->parent;
    if (pop->hybrid)
    {
      bidir = node_is_bidirection(pop);
      assert(!node_is_mirror(pop));

      hindex = GET_HINDEX(stree,pop);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      visited[hindex] = 1;

      if (stree->pptable[pop->node_index][end->node_index] &&
          stree->pptable[pop->hybrid->node_index][end->node_index])
      {
        if (legacy_rndu() <= pop->hgamma)
        {
          x->hpath[hindex] = BPP_HPATH_LEFT;
          contrib += log(pop->hgamma);
        }
        else
        {
          x->hpath[hindex] = BPP_HPATH_RIGHT;
          if (bidir)
          {
            /* do two steps (nodes) in one go */

            /* count contribution for first part of bidirection */
            contrib += log(pop->hybrid->hgamma);

            /* move to opposite node, set flag, and set visited */
            pop = pop->hybrid->parent;
            hindex = GET_HINDEX(stree,pop);
            x->hpath[hindex] = BPP_HPATH_LEFT;
            visited[hindex] = 1;
          }
          else
          {
            pop = pop->hybrid;
            contrib += log(pop->hgamma);
          }

        }
      }
      else
      {
        assert(stree->pptable[pop->node_index][end->node_index] ||
               stree->pptable[pop->hybrid->node_index][end->node_index]);

        if (stree->pptable[pop->node_index][end->node_index])
          x->hpath[hindex] = BPP_HPATH_LEFT;
        else
        {
          x->hpath[hindex] = BPP_HPATH_RIGHT;
          if (bidir)
          {
            /* do two steps (nodes) in one go */

            assert(pop->hybrid && pop->hybrid->parent);  /* TODO DEBUG */

            /* move to opposite node, set flag, and set visited */
            pop = pop->hybrid->parent;
            hindex = GET_HINDEX(stree,pop);
            x->hpath[hindex] = BPP_HPATH_LEFT;
            visited[hindex] = 1;
          }
          else
          {
            pop = pop->hybrid;
            assert(pop->parent);  /* TODO DEBUG */
          }
        }
      }
    }
    if (pop == end) break;
  }

  /* reset to BPP_HPATH_NONE for those hybridization nodes that were *NOT*
     visited, since the branch could have a flag for some of them from its
     previous state */
  for (i = 0; i < stree->hybrid_count; ++i)
    if (!visited[i])
      x->hpath[i] = BPP_HPATH_NONE;
  free(visited);

  return contrib;
}

static double sample_hpath_reverse(stree_t * stree, gnode_t * x, int * old_hpath)
{
  unsigned int hindex;
  long bidir = 0;
  snode_t * pop;
  snode_t * start = x->pop;
  snode_t * end   = x->parent ? x->parent->pop : stree->root;
  double contrib = 0;

  if (start == end) return contrib;

  assert(start->parent);

  /* sample flags for populations in between (excluding the two end points ) */
  for (pop = start->parent; pop != end; pop = pop->parent)
  {
    assert(pop);
    if (pop->hybrid)
    {
      bidir = node_is_bidirection(pop);
      assert(!node_is_mirror(pop));

      hindex = GET_HINDEX(stree,pop);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      if (stree->pptable[pop->node_index][end->node_index] &&
          stree->pptable[pop->hybrid->node_index][end->node_index])
      {
        if (old_hpath[hindex] == BPP_HPATH_LEFT)
        {
          contrib += log(pop->hgamma);
        }
        else
        {
          if (bidir)
          {
            contrib += log(pop->hybrid->hgamma);
            pop = pop->hybrid->parent;
          }
          else
          {
            pop = pop->hybrid;
            contrib += log(pop->hgamma);
          }
        }
      }
      else
      {
        assert(stree->pptable[pop->node_index][end->node_index] ||
               stree->pptable[pop->hybrid->node_index][end->node_index]);

        if (stree->pptable[pop->node_index][end->node_index])
        {
          //assert(old_hpath[hindex] == BPP_HPATH_LEFT);
        }
        else
        {
          //assert(old_hpath[hindex] == BPP_HPATH_RIGHT);
          if (bidir)
            pop = pop->hybrid->parent;
          else
            pop = pop->hybrid;
        }
      }
    }
    if (pop == end) break;
  }
  return contrib;

}

static void decrease_gene_leaves_count(stree_t * stree, gnode_t * x, int msa_index)
{
  unsigned int hindex;
  unsigned int leaves_count = x->leaves;
  snode_t * start = x->pop;
  snode_t * end = stree->root;

  if (start == end) return;
  
  while (start != end)
  {
    /* skip all ancestral gene nodes in population 'start' */
    while (x->parent && x->parent->pop == start)
    {
      x = x->parent;
    }

    /* move to ancestral population */
    start = start->parent;

    if (start->hybrid)
    {
      assert(!node_is_mirror(start));

      /* move according to path flags */
      hindex = GET_HINDEX(stree,start);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      /* find correct parent node according to hpath flag */
      assert(x->hpath[hindex] != BPP_HPATH_NONE);
      assert(start->left);
      if (x->hpath[hindex] == BPP_HPATH_RIGHT)
        start = start->hybrid;
    }
    start->gene_leaves[msa_index] -= leaves_count;
  }
}

static void increase_gene_leaves_count(stree_t * stree, gnode_t * x, int msa_index)
{
  unsigned int hindex;
  unsigned int leaves_count = x->leaves;
  snode_t * start = x->pop;
  snode_t * end = stree->root;

  if (start == end) return;

  while (start != end)
  {
    /* skip all ancestral gene nodes in population 'start' */
    while (x->parent && x->parent->pop == start)
    {
      x = x->parent;
    }

    /* move to ancestral population */
    start = start->parent;
    if (start->hybrid)
    {
      assert(!node_is_mirror(start));

      /* move according to path flags */
      hindex = GET_HINDEX(stree,start);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      /* find correct parent node according to hpath flag */
      assert(x->hpath[hindex] != BPP_HPATH_NONE);
      assert(start->left);
      if (x->hpath[hindex] == BPP_HPATH_RIGHT)
        start = start->hybrid;
    }

    start->gene_leaves[msa_index] += leaves_count;
  }
}

/* decreases seqin_count from pop(x) (exclusive) to pop(parent->x) (inclusive)
   unless pop(x) == pop(parent->x), for the species locus index */
static void decrease_seqin_count(stree_t * stree, gnode_t * x, int msa_index)
{
  unsigned int hindex;
  snode_t * start = x->pop;
  snode_t * end   = x->parent ? x->parent->pop : stree->root;

  if (start == end) return;

  while (start != end)
  {
    start = start->parent;

    if (start->hybrid)
    {
      assert(!node_is_mirror(start));

      /* 
        In the case of bidirectional introgression where the lineage goes right
        (i.e. to the node opposite), the algorithm will always visit the mirror
        node first, and immediately after the opposite hybrid node */
      
      /* move according to path flags */
      hindex = GET_HINDEX(stree,start);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      /* find correct parent node according to hpath flag */
      assert(x->hpath[hindex] != BPP_HPATH_NONE);
      assert(start->left);
      if (x->hpath[hindex] == BPP_HPATH_RIGHT)
        start = start->hybrid;
    }
    start->seqin_count[msa_index]--;
  }
}

/* increases seqin_count from pop(x) (exclusive) to pop(parent->x) (inclusive)
   unless pop(x) == pop(parent->x), for the species locus index */
static void increase_seqin_count(stree_t * stree, gnode_t * x, int msa_index)
{
  unsigned int hindex;
  snode_t * start = x->pop;
  snode_t * end   = x->parent ? x->parent->pop : stree->root;

  if (start == end) return;

  while (start != end)
  {
    start = start->parent;

    if (start->hybrid)
    {
      assert(!node_is_mirror(start));
      
      /* move according to path flags */
      hindex = GET_HINDEX(stree,start);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      assert(x->hpath[hindex] != BPP_HPATH_NONE);
      assert(start->left);
      if (x->hpath[hindex] == BPP_HPATH_RIGHT)
        start = start->hybrid;
    }

    start->seqin_count[msa_index]++;
  }
}

static snode_t * network_mrca_population(stree_t * stree, gnode_t * left, gnode_t * right)
{
  long i;
  snode_t * mrca = stree->root;

  /* use pop-pop table to find mrca of the two child populations */

  if (left->pop == right->pop) return left->pop;

  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
  {
    if (stree->pptable[left->pop->node_index][i] && 
        stree->pptable[right->pop->node_index][i] &&
        stree->nodes[i]->tau < mrca->tau)
    {
      mrca = stree->nodes[i];
    }
  }

  return mrca;
}

static long propose_ages(locus_t * locus, gtree_t * gtree, stree_t * stree, int msa_index)
{
  unsigned int i,k,j;
  unsigned int stree_total_nodes;
  long accepted = 0;
  double lnacceptance;
  double tnew,minage,maxage,oldage;
  double logpr;
  double logl;
  snode_t * pop;
  snode_t * oldpop;

  int * old_hpath_x = NULL;
  int * old_hpath_c1 = NULL;
  int * old_hpath_c2 = NULL;
  double hgamma_contrib = 0;
  double hgamma_contrib_reverse = 0;
  double hpop_contrib = 0;
  double hpop_contrib_reverse = 0;

  stree_total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  /* TODO: Instead of traversing the gene tree nodes this way, traverse the
     coalescent events for each population in the species tree instead. This
     will reduce the amount of required quick-sorts.
  */

  for (i = gtree->tip_count; i < gtree->inner_count+gtree->tip_count; ++i)
  {
    gnode_t * node = gtree->nodes[i];

    if (opt_network)
    {
      /* store sum of incoming lineages and coalescent events for detecting
         for which populations we need to recompute the MSC density */
      for (j = 0; j < stree_total_nodes; ++j)
      {
        stree->nodes[j]->hx = stree->nodes[j]->event_count[msa_index] +
                              stree->nodes[j]->seqin_count[msa_index];
      }

      /* TODO: The following loop corrects for bidirectional nodes */
      for (j = 0; j < stree->hybrid_count; ++j)
      {
        snode_t * snode = stree->nodes[stree->tip_count+stree->inner_count+j];

        if (node_is_bidirection(snode))
        {
          assert(snode->event_count[msa_index] == 0);
          snode->parent->hx -= snode->seqin_count[msa_index];
        }
      }
      hgamma_contrib = 0;
      hgamma_contrib_reverse = 0;
      hpop_contrib = 0;
      hpop_contrib_reverse = 0;
    }

    /* constraint min bound of proposed age by maximum age between children */
    minage = MAX(node->left->time,node->right->time);

    if (opt_network)
    {
      /* if the children are in different populations then further constraint
         minage by the tau of their most recent common ancestor population */
      snode_t * mrca = network_mrca_population(stree, node->left, node->right);

      minage = MAX(minage,mrca->tau);
    }
    else        /* not a network */
    {
      if (node->left->pop != node->right->pop)
      {
        snode_t * lpop = NULL;
        snode_t * rpop = NULL;

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
    }

    /* compute max age. TODO: 999 is placed for compatibility with old bpp */
    maxage = node->parent ? node->parent->time : 999;

    assert(maxage > minage);

    tnew = node->time + opt_finetune_gtage * legacy_rnd_symmetrical();
    tnew = reflect(tnew, minage, maxage);

    assert(tnew != 0);


    /* find the first ancestral pop with age higher than the proposed tnew */
    /* TODO: Improvement: probably this can start from lpop/rpop (LCA of
        populations of two daughter nodes) */
    if (opt_network)
    {
      /* allocate temporary storage */
      long cand_count = 0;
      snode_t ** candidates = (snode_t **)xmalloc((size_t)stree_total_nodes *
                                                  sizeof(snode_t *));
      snode_t * lpop = node->left->pop;
      snode_t * rpop = node->right->pop;

      /* find all feasible populations compatible with the new age */
      for (j = 0; j < stree_total_nodes; ++j)
      {
        snode_t * x = stree->nodes[j];
        if (stree->pptable[lpop->node_index][x->node_index] &&
            stree->pptable[rpop->node_index][x->node_index] &&
            (x->tau <= tnew) &&
            (!x->parent || x->parent->tau > tnew))
        {
          if (node->parent &&
              (!stree->pptable[x->node_index][node->parent->pop->node_index]))
            continue;

          candidates[cand_count++] = x;
        }
      }
      assert(cand_count > 0);
      
      /* randomly select one such population */
      j = (int)(cand_count*legacy_rndu());
      assert(j < cand_count);
      pop = candidates[j];

      /* compute the denominator for the hasting's correction */ 
      hpop_contrib = log(1.0 / cand_count);

      /* now compute the numberator for hasting's correction */
      cand_count = 0;
      for (j = 0; j < stree_total_nodes; ++j)
      {
        snode_t * x = stree->nodes[j];
        if (stree->pptable[lpop->node_index][x->node_index] &&
            stree->pptable[rpop->node_index][x->node_index] &&
            (x->tau <= node->time) &&
            (!x->parent || x->parent->tau > node->time))
        {
          if (node->parent &&
              (!stree->pptable[x->node_index][node->parent->pop->node_index]))
            continue;

          cand_count++;
        }
      }

      hpop_contrib_reverse = log(1.0 / cand_count);

      free(candidates);
    }
    else
    {   /* not a network */
      for (pop = node->left->pop; pop->parent; pop = pop->parent)
        if (pop->parent->tau > tnew)
          break;
    }

    if (opt_network)
    {
      /* Save old flags for the three nodes, to be used for rollback */
      old_hpath_x  = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));
      old_hpath_c1 = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));
      old_hpath_c2 = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));

      memcpy(old_hpath_x,  node->hpath, stree->hybrid_count*sizeof(int));
      memcpy(old_hpath_c1, node->left->hpath, stree->hybrid_count*sizeof(int));
      memcpy(old_hpath_c2, node->right->hpath, stree->hybrid_count*sizeof(int));

      /* Subtract seqin_count (nin) for the three branches */
      decrease_seqin_count(stree,node->left,msa_index);
      decrease_seqin_count(stree,node->right,msa_index);
      decrease_seqin_count(stree,node,msa_index);
      
      /* Subtract pop sizes */
      decrease_gene_leaves_count(stree,node->left,msa_index);
      decrease_gene_leaves_count(stree,node->right,msa_index);
    }

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
      if (!opt_est_theta)
        node->pop->event_count_sum--;
        
      /* change population for the current gene tree node */
      node->pop = pop;

      /* now add the coalescent event to the new population, at the end */
      dlist_item_append(node->pop->event[msa_index],node->event);

      node->pop->event_count[msa_index]++;
      if (!opt_est_theta)
        node->pop->event_count_sum++;

      /* increase or decrease the number of incoming lineages to all populations in the path
      from old population to the new population, depending on the case  */
      if (tnew > oldage)
      {
        /* increase the number of incoming lineages to all populations in the
           path from old population (excluding) to the new population */
        if (!opt_network)
        {
          for (pop = oldpop; pop != node->pop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]++;
        }
      }
      else  /* tnew < oldage */
      {
        /* decrease the number of incoming lineages to all populations in the
           path from new population (excluding) to the old population */
        if (!opt_network)
        {
          for (pop = node->pop; pop != oldpop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]--;
        }
      }
    }

    if (opt_network)
    {
      /* reset and sample new flags */
      hgamma_contrib += sample_hpath(stree,node);
      hgamma_contrib += sample_hpath(stree,node->left);
      hgamma_contrib += sample_hpath(stree,node->right);

      snode_t * newpop = node->pop;
      node->pop = oldpop;

      hgamma_contrib_reverse += sample_hpath_reverse(stree,node,old_hpath_x);
      hgamma_contrib_reverse += sample_hpath_reverse(stree,node->left,old_hpath_c1);
      hgamma_contrib_reverse += sample_hpath_reverse(stree,node->right,old_hpath_c2);

      node->pop = newpop;

      /* increase seqin_count (nin) for the three branches */
      increase_seqin_count(stree,node->left,msa_index);
      increase_seqin_count(stree,node->right,msa_index);
      increase_seqin_count(stree,node,msa_index);

      /* append pop sizes */
      increase_gene_leaves_count(stree,node->left,msa_index);
      increase_gene_leaves_count(stree,node->right,msa_index);

      node->pop->hx = -1;
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
    if (opt_est_theta)
      logpr = gtree->logpr;
    else
      logpr = stree->notheta_logpr;

    if (oldpop == node->pop)
    {
      /*TODO: BUG: separate to network and non-network case */
      if (!opt_network)
      {
        if (opt_est_theta)
          logpr -= node->pop->logpr_contrib[msa_index];
        else
          logpr -= node->pop->notheta_logpr_contrib;

        logpr += gtree_update_logprob_contrib(node->pop,
                                              locus->heredity[0],
                                              msa_index);
      }
      else
      {
        /* have a separate loop for hybrid mirror nodes to correct the problem
           of detecting populations that need MSC recomputation for
           bidirectional case */
        for (j = 0; j < stree->hybrid_count; ++j)
        {
          snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+j];

          /* correct for bidirectional introgression non-mirror nodes */
          if (node_is_bidirection(x) && x->parent->hx != -1)
              x->parent->hx += x->seqin_count[msa_index];

          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
            /* correct for bidirectional introgression non-mirror nodes */
            if (node_is_bidirection(x))
                x->parent->hx = -1;
          }
        }
        for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
        {
          snode_t * x = stree->nodes[j];
          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
          }
        }
      }
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

      if (opt_network)
      {
        /* have a separate loop for hybrid mirror nodes to correct the problem
           of detecting populations that need MSC recomputation for
           bidirectional case */
        for (j = 0; j < stree->hybrid_count; ++j)
        {
          snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+j];

          /* correct for bidirectional introgression non-mirror nodes */
          if (node_is_bidirection(x) && x->parent->hx != -1)
              x->parent->hx += x->seqin_count[msa_index];

          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
            /* correct for bidirectional introgression non-mirror nodes */
            if (node_is_bidirection(x))
                x->parent->hx = -1;
          }
        }
        for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        {
          snode_t * x = stree->nodes[j];
          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
          }
        }
      }
      else
      {
        for (pop = start; pop != end; pop = pop->parent)
        {
          if (opt_est_theta)
            logpr -= pop->logpr_contrib[msa_index];
          else
            logpr -= pop->notheta_logpr_contrib;

          logpr += gtree_update_logprob_contrib(pop,
                                                locus->heredity[0],
                                                msa_index);
        }
      }
    }
    #else
      assert(!opt_network);
      if (opt_est_theta)
        logpr = gtree_logprob(stree,locus->heredity[0],msa_index);
      else
        assert(0);

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
      if (opt_scaling)
        temp->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                               temp->scaler_index);
    }

    /* update partials */
    locus_update_partials(locus,travbuffer[msa_index],k);
    
    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

    if (opt_network)
    {
      lnacceptance = - (hgamma_contrib - hgamma_contrib_reverse);
      lnacceptance = lnacceptance - hpop_contrib + hpop_contrib_reverse;
    }
    else
      lnacceptance = 0;

    /* lnacceptance ratio */
    if (opt_est_theta)
      lnacceptance += logpr - gtree->logpr + logl - gtree->logl;
    else
      lnacceptance += logpr - stree->notheta_logpr + logl - gtree->logl;

    if (lnacceptance >= -1e-10 || legacy_rndu() < exp(lnacceptance))
    {
      /* accepted */
      accepted++;

      /* update new log-likelihood and gene tree log probability */
      if (opt_est_theta)
        gtree->logpr = logpr;
      else
        stree->notheta_logpr = logpr;

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
        if (opt_scaling)
          temp->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,temp->scaler_index);
      }
      
      /* now reset branch lengths and pmatrices */
      node->time = oldage;
      k = 0;
      travbuffer[msa_index][k++] = node->left;
      travbuffer[msa_index][k++] = node->right;
      if (node->parent)
        travbuffer[msa_index][k++] = node;
      locus_update_matrices_jc69(locus,travbuffer[msa_index],k);

      if (opt_network)
      {
        /* Subtract seqin_count (nin) for the three branches */
        decrease_seqin_count(stree,node->left,msa_index);
        decrease_seqin_count(stree,node->right,msa_index);
        decrease_seqin_count(stree,node,msa_index);
        
        /* Subtract pop sizes */
        decrease_gene_leaves_count(stree,node->left,msa_index);
        decrease_gene_leaves_count(stree,node->right,msa_index);
      }
          
      /* reset to old population, and reset gene tree log probability
         contributes for each modified species tree node */
      if (node->pop == oldpop)
      {
        if (opt_network)
        {
          memcpy(node->hpath, old_hpath_x,  stree->hybrid_count * sizeof(int));
          memcpy(node->left->hpath, old_hpath_c1, stree->hybrid_count * sizeof(int));
          memcpy(node->right->hpath, old_hpath_c2, stree->hybrid_count * sizeof(int));

          /* increase seqin_count (nin) for the three branches */
          increase_seqin_count(stree,node->left,msa_index);
          increase_seqin_count(stree,node->right,msa_index);
          increase_seqin_count(stree,node,msa_index);

          /* append pop sizes */
          increase_gene_leaves_count(stree,node->left,msa_index);
          increase_gene_leaves_count(stree,node->right,msa_index);
        }

        if (opt_network)
        {
          for (j = 0; j < stree_total_nodes; ++j)
          {
            snode_t * x = stree->nodes[j];
            if (x->hx)
            {
              if (opt_est_theta)
                x->logpr_contrib[msa_index] = x->old_logpr_contrib[msa_index];
              else
                logprob_revert_notheta(x,msa_index);
            }
          }
        }
        else
        {
          if (opt_est_theta)
            node->pop->logpr_contrib[msa_index] = node->pop->old_logpr_contrib[msa_index];
          else
            logprob_revert_notheta(node->pop,msa_index);
        }
      }
      else
      {
        /* remove current gene node from the list of coalescent events of its old
           population */

        unlink_event(node,msa_index);

        /* decrease the number of coalescent events for the current population */
        node->pop->event_count[msa_index]--;
        if (!opt_est_theta)
          node->pop->event_count_sum--;

        /* change population for the current gene tree node */
        SWAP(node->pop,oldpop);

        if (opt_network)
        {
          memcpy(node->hpath, old_hpath_x,  stree->hybrid_count * sizeof(int));
          memcpy(node->left->hpath, old_hpath_c1, stree->hybrid_count * sizeof(int));
          memcpy(node->right->hpath, old_hpath_c2, stree->hybrid_count * sizeof(int));

          /* increase seqin_count (nin) for the three branches */
          increase_seqin_count(stree,node->left,msa_index);
          increase_seqin_count(stree,node->right,msa_index);
          increase_seqin_count(stree,node,msa_index);

          /* append pop sizes */
          increase_gene_leaves_count(stree,node->left,msa_index);
          increase_gene_leaves_count(stree,node->right,msa_index);
        }

        /* now add the coalescent event back to the old population, at the end */
        dlist_item_append(node->pop->event[msa_index],node->event);

        node->pop->event_count[msa_index]++;
        if (!opt_est_theta)
          node->pop->event_count_sum++;

        /* increase or decrease the number of incoming lineages to all
           populations in the path from old population to the new population,
           depending on the case  */
        if (oldage >= tnew)
        {
          /* increase the number of incoming lineages to all populations in the
             path from old population (excluding) to the new population */
          if (!opt_network)
          {
            /* note that a few lines above we swapped the populations already
               ( SWAP(node->pop,oldpop) ) and so, oldpop is actually the new
               population, and node->pop is the old population, and so the for
               loop below is correct */
            for (pop=oldpop; pop != node->pop; pop = pop->parent)
              pop->parent->seqin_count[msa_index]++;
          }
        }
        else  /* tnew > oldage */
        {
          /* decrease the number of incoming lineages to all populations in the
             path from new population (excluding) to the old population */
          if (!opt_network)
          {
            for (pop = node->pop; pop != oldpop; pop = pop->parent)
              pop->parent->seqin_count[msa_index]--;
          }
        }

        /* now restore the old log gene tree probability contribution for each
           affected species tree node */
        
        snode_t * start = (tnew > oldage) ? node->pop :  oldpop;
        snode_t * end   = (tnew > oldage) ? oldpop->parent : node->pop->parent;

        if (opt_network)
        {
          for (j = 0; j < stree_total_nodes; ++j)
          {
            snode_t * x = stree->nodes[j];
            if (x->hx)
            {
              if (opt_est_theta)
                x->logpr_contrib[msa_index] = x->old_logpr_contrib[msa_index];
              else
                logprob_revert_notheta(x,msa_index);
            }
          }
        }
        else
        {
          for (pop = start; pop != end; pop = pop->parent)
          {
            if (opt_est_theta)
              pop->logpr_contrib[msa_index] = pop->old_logpr_contrib[msa_index];
            else
              logprob_revert_notheta(pop,msa_index);
          }
        }
      }
    }
    if (opt_network)
    {
      free(old_hpath_x);
      free(old_hpath_c1);
      free(old_hpath_c2);
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
  int ret = 0;
  gnode_t * sibling;
  gnode_t * father;
  gnode_t * oldroot = gtree->root;

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

    #if 1
    if (opt_network)
      SWAP(gtree->root->hpath, oldroot->hpath);
    #endif
    
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

static int branch_compat(stree_t * stree,
                         gnode_t * curnode,
                         gnode_t * target,
                         double tnew)
{
  unsigned int hindex;
  snode_t * target_pop;
  snode_t * end = target->parent ? target->parent->pop : stree->root;
  
  
  assert(target->time <= tnew);
  if (target->parent)
    assert(target->parent->time > tnew);

  target_pop = target->pop;
  while (target_pop != end)
  {
    snode_t * nextpop;

    if (!target_pop->parent) break;

    nextpop = target_pop->parent;
    if (nextpop->hybrid)
    {
      assert(!node_is_mirror(nextpop));

      hindex = GET_HINDEX(stree,nextpop);
      assert(hindex >= 0 && hindex < stree->hybrid_count);
      assert(target->hpath[hindex] != BPP_HPATH_NONE);

      if (target->hpath[hindex] == BPP_HPATH_RIGHT)
        nextpop = nextpop->hybrid;
    }

    if (nextpop->tau > tnew)
      break;

    target_pop = nextpop;
  }

  assert(target_pop);

  if (stree->pptable[curnode->pop->node_index][target_pop->node_index])
    return 1;

  return 0;
}

static gnode_t * network_fill_targets(stree_t * stree,
                                      gtree_t * gtree,
                                      int msa_index,
                                      gnode_t * curnode,
                                      gnode_t * father,
                                      gnode_t * sibling,
                                      double tnew,
                                      unsigned int * target_count,
                                      unsigned int * source_count,
                                      snode_t ** pop_target)
{
  unsigned int j,k,n,m;
  unsigned int ptarget_count = 0;
  unsigned int snodes_count;
  gnode_t * p;
  snode_t ** ptarget_list = NULL;

  *pop_target = NULL;

  /* make sure no species tree node is 'marked' */
  snodes_count = stree->tip_count + stree->inner_count + stree->hybrid_count;
  for (j = 0; j < snodes_count; ++j)
    assert(stree->nodes[j]->mark == 0);

  /* mark all species tree nodes nodes ancestor to curnode population, and with
     branches that include tnew. */
  for (j = 0; j < snodes_count; ++j)
  {
    snode_t * x = stree->nodes[j];

    if (stree->pptable[curnode->pop->node_index][x->node_index] &&
        //x->gene_leaves[msa_index] > curnode->leaves &&
        (x->tau <= tnew) && (x->parent && x->parent->tau > tnew))
    {
      x->mark = 1;
      ptarget_count++;
    }
  }
  if (stree->root->tau <= tnew)
  {
    stree->root->mark = 1;
    ptarget_count++;
  }

  /* fill a list with target populations */
  ptarget_list = (snode_t **)xmalloc((size_t)ptarget_count*sizeof(snode_t *));
  for (k = 0, j = 0; j < snodes_count; ++j)
  {
    if (stree->nodes[j]->mark)
    {
      stree->nodes[j]->mark = 0;
      ptarget_list[k++] = stree->nodes[j];
    }
  }
  assert(k == ptarget_count);

  /* identify target branches on which we can attach the pruned tree */
  *target_count = 0;
  if (tnew >= gtree->root->time)
  {
    travbuffer[msa_index][*target_count] = gtree->root;
    *target_count = *target_count + 1;
  }
  else
  {
    for (j = 0; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      p = gtree->nodes[j];
      m = p->pop->node_index;

      for (k = 0; k < ptarget_count; ++k)
      {
        n = ptarget_list[k]->node_index;

        if (p != curnode && p != gtree->root && p->time <= tnew &&
            p->parent->time > tnew && stree->pptable[m][n] &&
            branch_compat(stree,curnode,p,tnew))
        {
          travbuffer[msa_index][*target_count] = (p == father) ? sibling : p;
          *target_count = *target_count + 1;
          break;
        }
      }
    }
  }

  /* source count */
  *source_count = 1;
  if (father != gtree->root)
  {
      
    n = curnode->pop->node_index;
    for (j = 0; j < gtree->tip_count + gtree->inner_count; ++j)
    {
      p = gtree->nodes[j];
      m = p->parent ? p->parent->pop->node_index : stree->root->node_index;

      if (p != curnode && p != gtree->root && p != sibling && p != father &&
          p->time <= father->time && p->parent->time > father->time &&
          stree->pptable[n][m] && branch_compat(stree,curnode,p,father->time))
        *source_count = *source_count + 1;
    }
  }

  assert(*target_count);
  assert(*source_count);

  /* randomly select a target node */
  gnode_t * target = travbuffer[msa_index][(int)(*target_count*legacy_rndu())];
  assert(target);

  /* decide on the pop_target */
  snode_t * pop = target->pop;
  while (pop->parent && pop->parent->tau < tnew)
  {
    pop = pop->parent;
    if (pop->hybrid)
    {
      assert(!node_is_mirror(pop));

      unsigned int hindex = GET_HINDEX(stree,pop);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      int dbg_flag = target->hpath[hindex];
      if (target->hpath[hindex] == BPP_HPATH_NONE)
      {
        /* This case happens when the target node is the sibling of curnode
           (i.e. the other child of father) and tnew is older than father->time,
           and at the same time, the lineage from father to its parent passes
           through a hybridization event, and as such, target does not have that
           flag (it is in father).
        */

        assert(target == sibling);
        assert(father->time < tnew);
        assert(father->hpath[hindex] != BPP_HPATH_NONE);
        dbg_flag = father->hpath[hindex];
      }
      //assert(target->hpath[hindex] != BPP_HPATH_NONE);
      assert(dbg_flag != BPP_HPATH_NONE);
      //if (target->hpath[hindex] == BPP_HPATH_RIGHT)
      if (dbg_flag == BPP_HPATH_RIGHT)
        pop = pop->hybrid;
    }
  }
  assert(pop);
  *pop_target = pop;

  if (ptarget_list)
    free(ptarget_list);

  return target;
}

static long propose_spr(locus_t * locus,
                        gtree_t * gtree,
                        stree_t * stree,
                        int msa_index)
{
  unsigned int i,j,k,m,n;
  unsigned int source_count, target_count;
  unsigned int stree_total_nodes;
  long accepted = 0;
  gnode_t * curnode;
  gnode_t * sibling;
  gnode_t * father;
  gnode_t * p;
  double minage,maxage,tnew;
  double lnacceptance;
  double logpr;
  snode_t * pop;

  int * old_hpath_y = NULL;
  int * old_hpath_a = NULL;
  int * old_hpath_s = NULL;
  int * old_hpath_t = NULL;
  double hgamma_contrib = 0;
  double hgamma_contrib_reverse = 0;

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

  stree_total_nodes = stree->tip_count+stree->inner_count+stree->hybrid_count;

  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    curnode = gtree->nodes[i];
    if (curnode == gtree->root) continue;

    sibling = (curnode->parent->left == curnode) ? 
                curnode->parent->right : curnode->parent->left;
    father  = curnode->parent;
    
    assert(curnode->parent);

    if (opt_network)
    {
      /* store sum of incoming lineages and coalescent events for detecting
         for which populations we need to recompute the MSC density */
      for (j = 0; j < stree_total_nodes; ++j)
      {
        stree->nodes[j]->hx = stree->nodes[j]->event_count[msa_index] +
                              stree->nodes[j]->seqin_count[msa_index];
      }
      /* TODO: The following loop corrects for bidirectional nodes */
      for (j = 0; j < stree->hybrid_count; ++j)
      {
        snode_t * snode = stree->nodes[stree->tip_count+stree->inner_count+j];

        if (node_is_bidirection(snode))
        {
          assert(snode->event_count[msa_index] == 0);
          snode->parent->hx -= snode->seqin_count[msa_index];
        }
      }
      hgamma_contrib = 0;
      hgamma_contrib_reverse = 0;
    }

    /* find youngest population with subtree lineages more than current node subtree lineages */
    if (opt_network)
    {
      pop = stree->root;

      decrease_gene_leaves_count(stree,curnode,msa_index);
      curnode->pop->gene_leaves[msa_index] -= curnode->leaves;

      for (j = 0; j < stree_total_nodes; ++j)
      {
        snode_t * x = stree->nodes[j];

        if (stree->pptable[curnode->pop->node_index][x->node_index] &&
           (x->gene_leaves[msa_index] > 0) && (x->tau < pop->tau))
          pop = x;
      }
      curnode->pop->gene_leaves[msa_index] += curnode->leaves;
    }
    else
    {
      for (pop = curnode->pop;
           pop->gene_leaves[msa_index] <= curnode->leaves; 
           pop = pop->parent)
        if (!pop->parent)
          break;
    }

    /* TODO: Set age limits. 999 is set for backwards compatibility with bpp */
    minage = MAX(curnode->time, pop->tau);
    maxage = 999;

    tnew = father->time + opt_finetune_gtspr*legacy_rnd_symmetrical();
    tnew = reflect(tnew,minage,maxage);

    if (!opt_network)
    {
      for (pop = curnode->pop; pop->parent; pop = pop->parent)
        if (pop->parent->tau > tnew) break;
    }

    snode_t * pop_target = pop;

    /* identify target branches on which we can attach the pruned tree */
    /* TODO: We process the root node first to keep backwards compatibility with
       old bpp results */

    gnode_t * target;
    if (opt_network)
    {
      target = network_fill_targets(stree,
                                    gtree,
                                    msa_index,
                                    curnode,
                                    father,
                                    sibling,
                                    tnew,
                                    &target_count,
                                    &source_count,
                                    &pop_target);

      dbg_msci_t = target;
    }
    else
    {
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
      target = travbuffer[msa_index][(int)(target_count*legacy_rndu())];
    }

    if (opt_network)
    {
      /* Save old flags for the three nodes, to be used for rollback */
      old_hpath_y = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));
      old_hpath_a = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));
      old_hpath_s = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));
      old_hpath_t = (int *)xmalloc((size_t)(stree->hybrid_count)*sizeof(int));

      memcpy(old_hpath_y, father->hpath,  stree->hybrid_count*sizeof(int));
      memcpy(old_hpath_a, curnode->hpath, stree->hybrid_count*sizeof(int));
      memcpy(old_hpath_s, sibling->hpath, stree->hybrid_count*sizeof(int));
      memcpy(old_hpath_t, target->hpath,  stree->hybrid_count*sizeof(int));

      /* Subtract seqin_count (nin) for the three branches */
      decrease_seqin_count(stree,curnode,msa_index);
      
      dbg_msci_y = father;
      dbg_msci_s = sibling;
      dbg_msci_a = curnode;
    }

    /* regraft subtree */

    snode_t * oldpop = father->pop;
    double oldage = father->time;
    father->time = tnew;

    int spr_required = (target != sibling && target != father);
    gnode_t * dbg_old_father = father;

    if (father->pop != pop_target)
    {
      /* TODO: update coalescent events */

      /* remove current gene node from the list of coalescent events of its old
         population */

      unlink_event(father,msa_index);

      /* decrease the number of coalescent events for the current population */
      father->pop->event_count[msa_index]--;
      if (!opt_est_theta)
        father->pop->event_count_sum--;
        
      /* change population for the current gene tree node */
      father->pop = pop_target;

      /* now add the coalescent event to the new population, at the end */
      dlist_item_append(father->pop->event[msa_index],father->event);

      father->pop->event_count[msa_index]++;
      if (!opt_est_theta)
        father->pop->event_count_sum++;

      /* increase or decrease the number of incoming lineages to all populations
         in the path from old population to the new population, depending on the
         case  */
      if (tnew > oldage)
      {
        /* increase the number of incoming lineages to all populations in the
           path from old population (excluding) to the new population */
        if (!opt_network)
        {
          for (pop = oldpop; pop != father->pop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]++;
        }
      }
      else
      {
        /* decrease the number of incoming lineages to all populations in the
           path from new population (excluding) to the old population */
        if (!opt_network)
        {
          for (pop = father->pop; pop != oldpop; pop = pop->parent)
            pop->parent->seqin_count[msa_index]--;
        }
      }
    }
    
    /* if the following holds we need to change tree topology */
    int root_changed = 0;
    if (spr_required)
      root_changed = perform_spr(gtree,curnode,target);

    if (root_changed == 2)
    {
      assert(gtree->root == father);
      //assert(gtree->root == sibling);
      //assert(0);
    }

    gnode_t * original_father = father;

    if (root_changed)
      father = curnode->parent;

    /*******************************/
    if (root_changed == 1)
    {
       assert(father == target);
    }

    if (opt_network)
    {
      /* reset and sample new flags */
      hgamma_contrib += sample_hpath(stree,curnode);

      if (spr_required)
      {
        switch (root_changed)
        {
          case 0:
            join_flags(stree,sibling,old_hpath_y,oldpop);
            split_flags(stree,target);
            break;

          case 1:
            join_flags(stree,sibling,old_hpath_y,oldpop);
            split_flags(stree,original_father);
            break;

          case 2:
            join_flags(stree,gtree->root,old_hpath_y,oldpop);
            split_flags(stree,target);
            break;

          default:
            assert(0);
        }
      }
      else
      {
        /* redistribute flags between father and sibling */
        interchange_flags(stree,
                          sibling,
                          father,
                          tnew,
                          oldage,
                          oldpop);
      }

      snode_t * newpop = father->pop;
      father->pop = oldpop;

      hgamma_contrib_reverse += sample_hpath_reverse(stree,curnode,old_hpath_a);

      father->pop = newpop;

      /* increase seqin_count (nin) for the three branches */
      increase_seqin_count(stree, curnode, msa_index);

      /* append pop sizes */
      increase_gene_leaves_count(stree, curnode, msa_index);

      /* indicate we want to recompute the MSC density for this population */
      father->pop->hx = -1;
    }

    /* recompute logpr */
    if (opt_est_theta)
      logpr = gtree->logpr;
    else
      logpr = stree->notheta_logpr;

    if (oldpop == father->pop)
    {
      if (opt_network)
      {
        /* have a separate loop for hybrid mirror nodes to correct the problem
           of detecting populations that need MSC recomputation for
           bidirectional case */
        for (j = 0; j < stree->hybrid_count; ++j)
        {
          snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+j];

          /* correct for bidirectional introgression non-mirror nodes */
          if (node_is_bidirection(x) && x->parent->hx != -1)
              x->parent->hx += x->seqin_count[msa_index];

          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
            /* correct for bidirectional introgression non-mirror nodes */
            if (node_is_bidirection(x))
                x->parent->hx = -1;
          }
        }
        for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        {
          snode_t * x = stree->nodes[j];
          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,locus->heredity[0],msa_index);
          }
        }
      }
      else
      {
        if (opt_est_theta)
          logpr -= father->pop->logpr_contrib[msa_index];
        else
          logpr -= father->pop->notheta_logpr_contrib;

        logpr += gtree_update_logprob_contrib(father->pop,
                                              locus->heredity[0],
                                              msa_index);
      }
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

      /* TODO: The following code is identical to propose_ages(), with the
         exception that 'node' is 'father' here. Perhaps make it a function */
      if (opt_network)
      {
        /* have a separate loop for hybrid mirror nodes to correct the problem
           of detecting populations that need MSC recomputation for
           bidirectional case */
        for (j = 0; j < stree->hybrid_count; ++j)
        {
          snode_t * x = stree->nodes[stree->tip_count+stree->inner_count+j];

          /* correct for bidirectional introgression non-mirror nodes */
          if (node_is_bidirection(x) && x->parent->hx != -1)
              x->parent->hx += x->seqin_count[msa_index];

          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
            /* correct for bidirectional introgression non-mirror nodes */
            if (node_is_bidirection(x))
                x->parent->hx = -1;
          }
        }
        for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
        {
          snode_t * x = stree->nodes[j];
          if (x->seqin_count[msa_index] + x->event_count[msa_index] == x->hx)
          {
            x->hx = 0;  /* MSC density is not changed */
          }
          else
          {
            x->hx = 1;  /* MSC density has changed */
            if (opt_est_theta)
              logpr -= x->logpr_contrib[msa_index];
            else
              logpr -= x->notheta_logpr_contrib;

            logpr += gtree_update_logprob_contrib(x,
                                                  locus->heredity[0],
                                                  msa_index);
          }
        }
      }
      else
      {
        for (pop = start; pop != end; pop = pop->parent)
        {
          if (opt_est_theta)
            logpr -= pop->logpr_contrib[msa_index];
          else
            logpr -= pop->notheta_logpr_contrib;

          logpr += gtree_update_logprob_contrib(pop,
                                                locus->heredity[0],
                                                msa_index);
        }
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
        if (opt_scaling)
          temp->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,temp->scaler_index);
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
        if (opt_scaling)
          temp->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,temp->scaler_index);
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
        if (opt_scaling)
          temp->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,temp->scaler_index);
      }
    }

    /* update partials */
    locus_update_partials(locus,travbuffer[msa_index],k);

    /* compute log-likelihood */
    unsigned int param_indices[1] = {0};
    double logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

    if (opt_network)
      lnacceptance = hgamma_contrib_reverse - hgamma_contrib;
    else
      lnacceptance = 0;

    /* acceptance ratio */
    if (opt_est_theta)
      lnacceptance += log((double)target_count / source_count) +
                      logpr - gtree->logpr + logl - gtree->logl;
    else
    {
      lnacceptance += log((double)target_count / source_count) +
                      logpr - stree->notheta_logpr + logl - gtree->logl;
    }

    if (opt_debug)
      printf("[Debug] (spr) lnacceptance = %f\n", lnacceptance);

    if (lnacceptance >= -1e-10 || legacy_rndu() < exp(lnacceptance))
    {
      /* accepted */
      accepted++;

      if (opt_est_theta)
        gtree->logpr = logpr;
      else
        stree->notheta_logpr = logpr;

      gtree->logl = logl;
    }
    else
    {
      /* rejected */

      /* need to reset clv indices to point to the old clv buffer */
      for (j = 0; j < k; ++j)
      {
        gnode_t * temp = travbuffer[msa_index][j];
        temp->clv_index = SWAP_CLV_INDEX(gtree->tip_count,temp->clv_index);
        if (opt_scaling)
          temp->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,temp->scaler_index);
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

      if (opt_network)
      {
        assert(dbg_msci_a == curnode);
        
        /* Subtract seqin_count (nin) for the three branches */
        decrease_seqin_count(stree,curnode,msa_index);

        /* Subtract pop sizes */
        decrease_gene_leaves_count(stree,curnode,msa_index);
      }


      if (spr_required)
      {
        /* if root_changed == 2 it means that the old sibling is now the root,
           and because of the old bpp compatibility code, the node the
           represents the root is *not* the node that used to be the old
           sibling. Therefore, we must pass the root as the target of the SPR
           (and not the node pointing to the old sibling). In the other cases,
           the target is just the old sibling */
        if (root_changed == 2)
        {
          assert(gtree->root != sibling);
          assert(gtree->root == dbg_old_father);
          assert(sibling == curnode->parent);
          assert(father == sibling);
        }
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

      if (opt_network)
      {
        assert(dbg_msci_a == curnode);
        assert(dbg_msci_y == father);
        assert(dbg_msci_s == sibling);
        assert(dbg_msci_t == target);
      }

      if (father->pop == oldpop)
      {
        if (opt_network)
        {
          memcpy(father->hpath,  old_hpath_y, stree->hybrid_count*sizeof(int));
          memcpy(curnode->hpath, old_hpath_a, stree->hybrid_count*sizeof(int));
          memcpy(sibling->hpath, old_hpath_s, stree->hybrid_count*sizeof(int));
          memcpy(target->hpath,  old_hpath_t, stree->hybrid_count*sizeof(int));

          /* increase seqin_count (nin) for the three branches */
          increase_seqin_count(stree,curnode,msa_index);

          /* append pop sizes */
          increase_gene_leaves_count(stree,curnode,msa_index);
        }

        if (opt_network)
        {
          for (j = 0; j < stree_total_nodes; ++j)
          {
            snode_t * x = stree->nodes[j];
            if (x->hx)
            {
              if (opt_est_theta)
                x->logpr_contrib[msa_index] = x->old_logpr_contrib[msa_index];
              else
                logprob_revert_notheta(x,msa_index);
            }
          }
        }
        else
        {
          if (opt_est_theta)
            father->pop->logpr_contrib[msa_index] = father->pop->old_logpr_contrib[msa_index];
          else
          {
            /* TODO: The below code is the same as calling logprob_revert_notheta(father->pop,msa_index) */
            father->pop->t2h_sum -= father->pop->t2h[msa_index];
            father->pop->t2h[msa_index] = father->pop->old_t2h[msa_index];
            father->pop->t2h_sum += father->pop->t2h[msa_index];
            father->pop->notheta_logpr_contrib= father->pop->notheta_old_logpr_contrib;
          }
        }
      }
      else
      {
        /* remove current gene node from the list of coalescent events of its old
           population */

        unlink_event(father,msa_index);

        /* decrease the number of coalescent events for the current population */
        father->pop->event_count[msa_index]--;
        if (!opt_est_theta)
          father->pop->event_count_sum--;
          
        /* change population for the current gene tree node */
        SWAP(father->pop,oldpop);

        if (opt_network)
        {
          memcpy(father->hpath,  old_hpath_y, stree->hybrid_count*sizeof(int));
          memcpy(curnode->hpath, old_hpath_a, stree->hybrid_count*sizeof(int));
          memcpy(sibling->hpath, old_hpath_s, stree->hybrid_count*sizeof(int));
          memcpy(target->hpath,  old_hpath_t, stree->hybrid_count*sizeof(int));

          /* increase seqin_count (nin) for the three branches */
          increase_seqin_count(stree,curnode,msa_index);

          /* append pop sizes */
          increase_gene_leaves_count(stree,curnode,msa_index);
        }

        /* now add the coalescent event back to the old population, at the end */
        dlist_item_append(father->pop->event[msa_index],father->event);

        father->pop->event_count[msa_index]++;
        if (!opt_est_theta)
          father->pop->event_count_sum++;

        /* increase or decrease the number of incoming lineages to all populations in the path
        from old population to the new population, depending on the case  */
        if (oldage >= tnew)
        {
          /* TODO : Now this is covered here!!! */

          /* increase the number of incoming lineages to all populations in the
             path from old population (excluding) to the new population */
          if (!opt_network)
          {
            for (pop=oldpop; pop != father->pop; pop = pop->parent)
              pop->parent->seqin_count[msa_index]++;
          }
        }
        else
        {
          /* decrease the number of incoming lineages to all populations in the
             path from new population (excluding) to the old population */
          if (!opt_network)
          {
            for (pop = father->pop; pop != oldpop; pop = pop->parent)
              pop->parent->seqin_count[msa_index]--;
          }
        }

        /* now restore the old log gene tree probability contribution for each
           affected species tree node */
        
        snode_t * start = (tnew > oldage) ? father->pop :  oldpop;
        snode_t * end   = (tnew > oldage) ? oldpop->parent : father->pop->parent;

        if (opt_network)
        {
          for (j = 0; j < stree_total_nodes; ++j)
          {
            snode_t * x = stree->nodes[j];
            if (x->hx)
            {
              if (opt_est_theta)
                x->logpr_contrib[msa_index] = x->old_logpr_contrib[msa_index];
              else
                logprob_revert_notheta(x,msa_index);
            }
          }
        }
        else
        {
          for (pop = start; pop != end; pop = pop->parent)
          {
            if (opt_est_theta)
              pop->logpr_contrib[msa_index] = pop->old_logpr_contrib[msa_index];
            else
              logprob_revert_notheta(pop,msa_index);
          }
        }
      }

      if (opt_network)
      {
        assert(dbg_msci_y == father);
        assert(curnode->hpath[0] == old_hpath_a[0]);
        assert(father->hpath[0] == old_hpath_y[0]);
        assert(sibling->hpath[0] == old_hpath_s[0]);
        assert(target->hpath[0] == old_hpath_t[0]);
      }
    }

    /* delete rollback */
    if (opt_network)
    {
      free(old_hpath_y);
      free(old_hpath_a);
      free(old_hpath_s);
      free(old_hpath_t);
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

static long prop_locusrate(gtree_t ** gtree, stree_t * stree, locus_t ** locus)
{
  long i,j;
  long ref;
  long accepted = 0;
  double lnacceptance; 
  double new_locrate;
  double new_refrate;
  double old_locrate;
  double old_refrate;

  /* set reference locus as the one with the highest number of site patterns */
  for (i = 1, ref = 0; i < opt_locus_count; ++i)
    if (locus[i]->sites > locus[ref]->sites)
      ref = i;

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (i == ref) continue;

    gnode_t ** refnodes = gtree[ref]->nodes;
    for (j = 0; j < gtree[ref]->tip_count + gtree[ref]->inner_count; ++j)
      if (refnodes[j]->parent)
        SWAP_PMAT_INDEX(gtree[ref]->edge_count, refnodes[j]->pmatrix_index);

    gnode_t ** locnodes = gtree[i]->nodes;
    for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      if (locnodes[j]->parent)
        SWAP_PMAT_INDEX(gtree[i]->edge_count,locnodes[j]->pmatrix_index);

    old_locrate = locus[i]->mut_rates[0];
    old_refrate = locus[ref]->mut_rates[0];

    new_locrate = reflect(old_locrate+opt_finetune_locusrate*legacy_rnd_symmetrical(),
                          0, old_locrate + old_refrate);
    new_refrate = locus[ref]->mut_rates[0] - (new_locrate - old_locrate);

    locus[i]->mut_rates[0] = new_locrate;
    locus[ref]->mut_rates[0] = new_refrate;

    lnacceptance = (opt_locusrate_alpha - 1) *
                   log((new_locrate*new_refrate) / (old_locrate*old_refrate));

    /* update selected locus */
    locus_update_all_matrices_jc69(locus[i],gtree[i]);

    gnode_t ** gnodeptr = gtree[i]->nodes;
    for (j = gtree[i]->tip_count; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                              gnodeptr[j]->clv_index);
      if (opt_scaling)
        gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                      gnodeptr[j]->scaler_index);
    }

    locus_update_all_partials(locus[i],gtree[i]);

    /* update reference locus */
    locus_update_all_matrices_jc69(locus[ref],gtree[ref]);

    gnodeptr = gtree[ref]->nodes;
    for (j = gtree[ref]->tip_count; j < gtree[ref]->tip_count+gtree[ref]->inner_count; ++j)
    {
      gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[ref]->tip_count,
                                              gnodeptr[j]->clv_index);
      if (opt_scaling)
         gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[ref]->tip_count,
                                                      gnodeptr[j]->scaler_index);
    }
    locus_update_all_partials(locus[ref],gtree[ref]);

    unsigned int param_indices[1] = {0};

    double loc_logl = locus_root_loglikelihood(locus[i],
                                               gtree[i]->root,
                                               param_indices,
                                               NULL);
    double ref_logl = locus_root_loglikelihood(locus[ref],
                                               gtree[ref]->root,
                                               param_indices,
                                               NULL);

    lnacceptance += loc_logl - gtree[i]->logl + ref_logl - gtree[ref]->logl;

    if (opt_debug)
      fprintf(stdout, "[Debug] (locusrate) lnacceptance = %f\n", lnacceptance);

    if (lnacceptance >= -1e-10 || legacy_rndu() < exp(lnacceptance))
    {
      /* accept */
      accepted++;

      gtree[i]->logl = loc_logl;
      gtree[ref]->logl = ref_logl;
    }
    else
    {
      /* reject */
      locus[i]->mut_rates[0] = old_locrate;
      locus[ref]->mut_rates[0] = old_refrate;

      /* reset selected locus */
      gnodeptr = gtree[i]->nodes;
      for (j = gtree[i]->tip_count; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
      {
        gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[i]->tip_count,
                                                gnodeptr[j]->clv_index);
        if (opt_scaling)
          gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[i]->tip_count,
                                                        gnodeptr[j]->scaler_index);
      }

      /* reset reference locus */
      gnodeptr = gtree[ref]->nodes;
      for (j = gtree[ref]->tip_count; j < gtree[ref]->tip_count+gtree[ref]->inner_count; ++j)
      {
        gnodeptr[j]->clv_index = SWAP_CLV_INDEX(gtree[ref]->tip_count,
                                                gnodeptr[j]->clv_index);
        if (opt_scaling)
          gnodeptr[j]->scaler_index = SWAP_SCALER_INDEX(gtree[ref]->tip_count,
                                                        gnodeptr[j]->scaler_index);
      }
      
      for (j = 0; j < gtree[ref]->tip_count + gtree[ref]->inner_count; ++j)
        if (refnodes[j]->parent)
          SWAP_PMAT_INDEX(gtree[ref]->edge_count,refnodes[j]->pmatrix_index);
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        if (locnodes[j]->parent)
          SWAP_PMAT_INDEX(gtree[i]->edge_count,locnodes[j]->pmatrix_index);
    }
  }
  return accepted;
}

static long prop_heredity(gtree_t ** gtree, stree_t * stree, locus_t ** locus)
{
  long i,j;
  long accepted = 0;
  double hnew,hold;
  double lnacceptance;
  double logpr = 0;

  double hfactor = 0;

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (!opt_est_theta)
      logpr = stree->notheta_logpr;


    hold = locus[i]->heredity[0];
    hnew = hold + opt_finetune_locusrate*legacy_rnd_symmetrical();
    if (hnew < 0) hnew *= -1;
    
    locus[i]->heredity[0] = hnew;
    lnacceptance = (opt_heredity_alpha-1)*log(hnew/hold) -
                   opt_heredity_beta*(hnew-hold);

    if (opt_est_theta)
      logpr = gtree_logprob(stree,locus[i]->heredity[0],i);
    else
    {
      for (j = 0; j < stree->tip_count+stree->inner_count; ++j)
      {
        logpr -= stree->nodes[j]->notheta_logpr_contrib;
        logpr += gtree_update_logprob_contrib(stree->nodes[j],locus[i]->heredity[0],i);
      }
    }

    /* TODO: Perhaps we can avoid the check every 100-th term by using the log
       of heredity scaler from the beginning. E.g. if this loop is replaced by
       
       for (j = 0; j < opt_locus_count; ++j)
         logpr -= log(locus[j]->heredity[0]);

       then we only need to add and subtract the two corresponding heredity
       multipliers (the old and new)
    */
    if (!opt_est_theta)
    {
      hfactor = 0;
      double y = 1;
      for (j = 0; j < opt_locus_count; ++j)
      {
        y *= locus[j]->heredity[0];
        if ((j+1) % 100 == 0)
        {
          hfactor -= log(y);
          y = 1;
        }
      }
      hfactor -= log(y);
      logpr += hfactor - stree->notheta_hfactor;
    }

    if (opt_est_theta)
      lnacceptance += logpr - gtree[i]->logpr;
    else
      lnacceptance += logpr - stree->notheta_logpr;

    if (opt_debug)
      fprintf(stdout, "[Debug] (heredity) lnacceptance = %f\n", lnacceptance);
    if (lnacceptance >= -1e-10 || legacy_rndu() < exp(lnacceptance))
    {
      /* accepted */
      accepted++;

      if (opt_est_theta)
        gtree[i]->logpr = logpr;
      else
      {
        stree->notheta_logpr = logpr;
        stree->notheta_hfactor = hfactor;
      }
    }
    else
    {
      /* rejected */
      locus[i]->heredity[0] = hold;
      for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      {
        if (opt_est_theta)
          stree->nodes[j]->logpr_contrib[i] = stree->nodes[j]->old_logpr_contrib[i];
        else
          logprob_revert_notheta(stree->nodes[j],i);
      }
    }
  }
  return accepted;
}

double prop_locusrate_and_heredity(gtree_t ** gtree, stree_t * stree, locus_t ** locus)
{
  long accepted = 0;
  double divisor = 0;

  if (opt_est_locusrate)
    accepted = prop_locusrate(gtree,stree,locus);

  if (opt_est_heredity)
    accepted += prop_heredity(gtree,stree,locus);

  if (opt_est_locusrate)
    divisor = opt_locus_count-1;
  if (opt_est_heredity)
    divisor += opt_locus_count;

  return (accepted / divisor);
}

