/*
    Copyright (C) 2015 Tomas Flouri

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

static int indend_space = 4;

static void print_node_info(rtree_t * tree)
{
  printf (" %s", tree->label);
  printf (" %f", tree->length);
  printf("\n");
}

static void print_tree_recurse(rtree_t * tree,
                               int indend_level,
                               int * active_node_order)
{
  int i,j;

  if (!tree) return;

  for (i = 0; i < indend_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indend_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indend_space-1; ++j)
    printf ("-");
  if (tree->left || tree->right) printf("+");

  print_node_info(tree);

  if (active_node_order[indend_level-1] == 2)
    active_node_order[indend_level-1] = 0;

  active_node_order[indend_level] = 1;
  print_tree_recurse(tree->left,
                     indend_level+1,
                     active_node_order);
  active_node_order[indend_level] = 2;
  print_tree_recurse(tree->right,
                     indend_level+1,
                     active_node_order);

}

static int tree_indend_level(rtree_t * tree, int indend)
{
  if (!tree) return indend;

  int a = tree_indend_level(tree->left,  indend+1);
  int b = tree_indend_level(tree->right, indend+1);

  return (a > b ? a : b);
}

void rtree_show_ascii(rtree_t * tree)
{

  int indend_max = tree_indend_level(tree,0);

  int * active_node_order = (int *)malloc((size_t)(indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_node_info(tree);
  print_tree_recurse(tree->left,  1, active_node_order);
  active_node_order[0] = 2;
  print_tree_recurse(tree->right, 1, active_node_order);
  free(active_node_order);
}

static char * rtree_export_newick_recursive(rtree_t * root)
{
  char * newick;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (asprintf(&newick, "%s:%f", root->label, root->length) == -1)
      fatal("Unable to allocate enough memory.");
  }
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left);
    char * subtree2 = rtree_export_newick_recursive(root->right);

    if (asprintf(&newick, "(%s,%s)%s:%f", subtree1,
                                          subtree2,
                                          root->label,
                                          root->length) == -1)
      fatal("Unable to allocate enough memory.");

    free(subtree1);
    free(subtree2);
  }

  return newick;
}

char * rtree_export_newick(rtree_t * root)
{
  char * newick;

  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (asprintf(&newick, "%s:%f", root->label, root->length) == -1)
      fatal("Unable to allocate enough memory.");
  }
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left);
    char * subtree2 = rtree_export_newick_recursive(root->right);

    if (asprintf(&newick, "(%s,%s)%s:%f;", subtree1,
                                           subtree2,
                                           root->label,
                                           root->length) == -1)
      fatal("Unable to allocate enough memory.");

    free(subtree1);
    free(subtree2);
  }

  return newick;
}

static void rtree_traverse_recursive(rtree_t * node,
                                     int (*cbtrav)(rtree_t *),
                                     int * index,
                                     struct drand48_data * rstate,
                                     rtree_t ** outbuffer)
{
  double rand_double = 0;

  if (!node->left)
  {
    if (!cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
  {
    outbuffer[*index] = node;
    *index = *index + 1;
    return;
  }

  drand48_r(rstate, &rand_double);
  if (rand_double >= 0.5)
  {
    rtree_traverse_recursive(node->left, cbtrav, index, rstate, outbuffer);
    rtree_traverse_recursive(node->right, cbtrav, index, rstate, outbuffer);
  }
  else
  {
    rtree_traverse_recursive(node->right, cbtrav, index, rstate, outbuffer);
    rtree_traverse_recursive(node->left, cbtrav, index, rstate, outbuffer);
  }

}

int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   struct drand48_data * rstate,
                   rtree_t ** outbuffer)
{
  int index = 0;

  if (!root->left) return -1;

  /* we will traverse an rooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  rtree_traverse_recursive(root, cbtrav, &index, rstate, outbuffer);
  return index;
}

static void rtree_traverse_postorder_recursive(rtree_t * node,
                                               int (*cbtrav)(rtree_t *),
                                               int * index,
                                               rtree_t ** outbuffer)
{
  if (!node) return;

  rtree_traverse_postorder_recursive(node->left,  cbtrav, index, outbuffer);
  rtree_traverse_postorder_recursive(node->right, cbtrav, index, outbuffer);

  if (cbtrav(node))
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }
}


int rtree_traverse_postorder(rtree_t * root,
                             int (*cbtrav)(rtree_t *),
                             rtree_t ** outbuffer)
{
  int index = 0;

  if (!root->left) return -1;

  /* we will traverse an unrooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether to
     place the node in the list */

  rtree_traverse_postorder_recursive(root, cbtrav, &index, outbuffer);
  return index;
}

static int rtree_height_recursive(rtree_t * node)
{
  if (!node) return 1;

  int a = rtree_height_recursive(node->left);
  int b = rtree_height_recursive(node->right);

  return MAX(a,b)+1;
}


int rtree_height(rtree_t * root)
{
  return rtree_height_recursive(root);
}

static void rtree_query_tipnodes_recursive(rtree_t * node,
                                           rtree_t ** node_list,
                                           int * index)
{
  if (!node) return;

  if (!node->left)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  rtree_query_tipnodes_recursive(node->left,  node_list, index);
  rtree_query_tipnodes_recursive(node->right, node_list, index);
}

int rtree_query_tipnodes(rtree_t * root,
                         rtree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left)
  {
    node_list[index++] = root;
    return index;
  }

  rtree_query_tipnodes_recursive(root->left,  node_list, &index);
  rtree_query_tipnodes_recursive(root->right, node_list, &index);

  return index;
}

static void rtree_query_innernodes_recursive(rtree_t * root,
                                             rtree_t ** node_list,
                                             int * index)
{
  if (!root) return;
  if (!root->left) return;

  /* postorder traversal */

  rtree_query_innernodes_recursive(root->left,  node_list, index);
  rtree_query_innernodes_recursive(root->right, node_list, index);

  node_list[*index] = root;
  *index = *index + 1;
  return;
}

int rtree_query_innernodes(rtree_t * root,
                           rtree_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left) return 0;

  rtree_query_innernodes_recursive(root->left,  node_list, &index);
  rtree_query_innernodes_recursive(root->right, node_list, &index);

  node_list[index++] = root;

  return index;
}

void rtree_print_tips(rtree_t * node, FILE * out)
{
  if (node->left)  rtree_print_tips(node->left,out);
  if (node->right) rtree_print_tips(node->right,out);

  if (!node->left && !node->right)
    fprintf(out, "%s\n", node->label);
}


rtree_t * rtree_clone(rtree_t * node, rtree_t * parent)
{
  if (!node) return NULL;

  /* clone node */
  rtree_t * clone = (rtree_t *)xmalloc(sizeof(rtree_t));
  memcpy(clone,node,sizeof(rtree_t));
  clone->parent = parent;

  if (node->label)
    clone->label = xstrdup(node->label);

  /* clone the two subtrees */
  clone->left  = rtree_clone(node->left, clone);
  clone->right = rtree_clone(node->right, clone);

  return clone;
}

rtree_t ** rtree_tipstring_nodes(rtree_t * root,
                                 char * tipstring,
                                 unsigned int * tiplist_count)
{
  size_t i;
  unsigned int k;
  unsigned int commas_count = 0;

  char * taxon;
  unsigned long taxon_len;

  ENTRY * found = NULL;

  for (i = 0; i < strlen(tipstring); ++i)
    if (tipstring[i] == ',')
      commas_count++;

  rtree_t ** node_list = (rtree_t **)xmalloc((size_t)(root->leaves) *
                                             sizeof(rtree_t *));
  rtree_query_tipnodes(root, node_list);

  rtree_t ** out_node_list = (rtree_t **)xmalloc((size_t)(commas_count+1) *
                                                 sizeof(rtree_t *));

  /* create a hashtable of tip labels */
  hcreate(2 * (size_t)(root->leaves));

  for (i = 0; i < (unsigned int)(root->leaves); ++i)
  {
    ENTRY entry;
    entry.key  = node_list[i]->label;
    entry.data = node_list[i];
    hsearch(entry,ENTER);
  }

  char * s = tipstring;

  k = 0;
  while (*s)
  {
    /* get next tip */
    taxon_len = strcspn(s, ",");
    if (!taxon_len)
      fatal("Erroneous prune list format (double comma)/taxon missing");

    taxon = xstrndup(s, taxon_len);

    /* search tip in hash table */
    ENTRY query;
    query.key = taxon;
    found = NULL;
    found = hsearch(query,FIND);

    if (!found)
      fatal("Taxon %s in does not appear in the tree", taxon);

    /* store pointer in output list */
    out_node_list[k++] = (rtree_t *)(found->data);

    /* free tip label, and move to the beginning of next tip if available */
    free(taxon);
    s += taxon_len;
    if (*s == ',')
      s += 1;
  }

  /* kill the hash table */
  hdestroy();

  free(node_list);

  /* return number of tips in the list */
  *tiplist_count = commas_count + 1;

  /* return tip node list */
  return out_node_list;
}
