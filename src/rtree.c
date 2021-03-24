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

static int indent_space = 4;

char * cb_stree_print_none(const snode_t * node)
{
  char * s = NULL;

  if (node->left)
  {
    /* replacement to keep GCC warning from showing up when executing
       the following:
       
       xasprintf(&s, "");
    */

    s = (char *)xmalloc(sizeof(char));
    *s = 0;
  }
  else
  {
    xasprintf(&s, "%s", node->label);
  }

  return s;
}
char * cb_gtree_print_none(const gnode_t * node)
{
  char * s = NULL;

  if (node->left)
  {
    /* replacement to keep GCC warning from showing up when executing
       the following:
       
       xasprintf(&s, "");
    */

    s = (char *)xmalloc(sizeof(char));
    *s = 0;
  }
  else
  {
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

char * cb_stree_print_tau(const snode_t * node)
{
  char * s = NULL;

  if (node->left)
  {
    if (!node->tau)
    {
      /* replacement to keep GCC warning from showing up when executing
         the following:
         
         xasprintf(&s, "");
      */

      s = (char *)xmalloc(sizeof(char));
      *s = 0;
    }
    else
    {
      xasprintf(&s, ":%f", node->tau);
    }
  }
  else
  {
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

char * cb_stree_print_node_tau(const snode_t * node)
{
  char * s = NULL;

  if (node->left)
  {
    if (!node->tau)
    {
      /* replacement to keep GCC warning from showing up when executing
         the following:
         
         xasprintf(&s, "");
      */

      s = (char *)xmalloc(sizeof(char));
      *s = 0;
    }
    else
    {
      xasprintf(&s, ":%d:%f", node->node_index, node->tau);
    }
  }
  else
  {
    xasprintf(&s, "%s:%d", node->label, node->node_index);
  }

  return s;
}

char * cb_gtree_print_age(const gnode_t * node)
{
  char * s = NULL;

  if (node->left)
  {
    xasprintf(&s, ":%f", node->time);
  }
  else
  {
    xasprintf(&s, "%s", node->label);
  }

  return s;
}

char* cb_gtree_print_length(const gnode_t* node)
{
   char* s = NULL;

   if (node->left)
   {
      xasprintf(&s, ":%f", node->length);
   }
   else
   {
      xasprintf(&s, "%s:%f", node->label, node->length);
   }

   return s;
}


char * cb_gtree_print_node_age(const gnode_t * node)
{
  char * s = NULL;

  if (node->left)
  {
    xasprintf(&s, ":%d:%f", node->node_index, node->time);
  }
  else
  {
    xasprintf(&s, "%s:%d", node->label, node->node_index);
  }

  return s;
}

static void print_node_info(const snode_t * root, int options)
{
  if (options & RTREE_SHOW_LABEL)
    printf (" %s", root->label);
  if (options & RTREE_SHOW_BRANCH_LENGTH)
    printf (" %f", root->length);
  printf("\n");
}

static void print_tree_recurse(const snode_t * root,
                               int indent_level,
                               int * active_node_order,
                               int options)
{
  int i,j;

  if (!root) return;

  for (i = 0; i < indent_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indent_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indent_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indent_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indent_space-1; ++j)
    printf ("-");
  if (root->left || root->right) printf("+");

  print_node_info(root, options);

  if (active_node_order[indent_level-1] == 2)
    active_node_order[indent_level-1] = 0;

  active_node_order[indent_level] = 1;
  if (root->left)
    print_tree_recurse(root->left,
                       indent_level+1,
                       active_node_order,
                       options);
  active_node_order[indent_level] = 2;
  if (root->right)
    print_tree_recurse(root->right,
                       indent_level+1,
                       active_node_order,
                       options);

}

static unsigned int tree_indent_level(const snode_t * root, unsigned int indent)
{
  if (!root) return indent;

  unsigned int a = tree_indent_level(root->left,  indent+1);
  unsigned int b = tree_indent_level(root->right, indent+1);

  return (a > b ? a : b);
}

void stree_show_ascii(const snode_t * root, int options)
{

  unsigned int indent_max = tree_indent_level(root,0);

  int * active_node_order = (int *)xmalloc((indent_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  print_node_info(root, options);
  print_tree_recurse(root->left,  1, active_node_order, options);
  print_tree_recurse(root->right, 1, active_node_order, options);
  free(active_node_order);
}

static char * stree_export_newick_recursive(const snode_t * root,
                                            char * (*cb_serialize)(const snode_t *))
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
    char * subtree1 = stree_export_newick_recursive(root->left,cb_serialize);
    if (subtree1 == NULL)
    {
      return NULL;
    }
    char * subtree2 = stree_export_newick_recursive(root->right,cb_serialize);
    if (subtree2 == NULL)
    {
      free(subtree1);
      return NULL;
    }

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = xasprintf(&newick,
                               "(%s, %s)%s",
                               subtree1,
                               subtree2,
                               temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
                               "(%s, %s)%s:%f",
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

char * stree_export_newick(const snode_t * root, char * (*cb_serialize)(const snode_t *))
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
      size_alloced = xasprintf(&newick, "%s:%f", root->label, root->length);
  }
  else
  {
    char * subtree1 = stree_export_newick_recursive(root->left,cb_serialize);
    if (!subtree1)
      fatal("Unable to allocate enough memory.");

    char * subtree2 = stree_export_newick_recursive(root->right,cb_serialize);
    if (!subtree2)
      fatal("Unable to allocate enough memory.");

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = xasprintf(&newick, "(%s, %s)%s;", subtree1, subtree2, temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick, "(%s, %s)%s:%f;", subtree1, subtree2, root->label ? root->label : "", root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("memory allocation during newick export failed");

  return newick;
}


static void stree_traverse_postorder(snode_t * node,
                                     int (*cbtrav)(snode_t *),
                                     unsigned int * index,
                                     snode_t ** outbuffer)
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

  stree_traverse_postorder(node->left, cbtrav, index, outbuffer);
  stree_traverse_postorder(node->right, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

static void stree_traverse_preorder(snode_t * node,
                                    int (*cbtrav)(snode_t *),
                                    unsigned int * index,
                                    snode_t ** outbuffer)
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

  stree_traverse_preorder(node->left, cbtrav, index, outbuffer);
  stree_traverse_preorder(node->right, cbtrav, index, outbuffer);

}

int stree_traverse(snode_t * root,
                   int traversal,
                   int (*cbtrav)(snode_t *),
                   snode_t ** outbuffer,
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
    stree_traverse_postorder(root, cbtrav, trav_size, outbuffer);
  else if (traversal == TREE_TRAVERSE_PREORDER)
    stree_traverse_preorder(root, cbtrav, trav_size, outbuffer);
  else
    fatal("Invalid traversal value.");

  return BPP_SUCCESS;
}

#if 0
static void stree_query_tipnodes_recursive(pll_stree_t * node,
                                           pll_stree_t ** node_list,
                                           unsigned int * index)
{
  if (!node) return;

  if (!node->left)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  stree_query_tipnodes_recursive(node->left,  node_list, index);
  stree_query_tipnodes_recursive(node->right, node_list, index);
}

PLL_EXPORT unsigned int pll_stree_query_tipnodes(pll_stree_t * root,
                                                 pll_stree_t ** node_list)
{
  unsigned int index = 0;

  if (!root) return 0;
  if (!root->left)
  {
    node_list[index++] = root;
    return index;
  }

  stree_query_tipnodes_recursive(root->left,  node_list, &index);
  stree_query_tipnodes_recursive(root->right, node_list, &index);

  return index;
}

static void stree_query_innernodes_recursive(pll_stree_t * root,
                                             pll_stree_t ** node_list,
                                             unsigned int * index)
{
  if (!root) return;
  if (!root->left) return;

  /* postorder traversal */

  stree_query_innernodes_recursive(root->left,  node_list, index);
  stree_query_innernodes_recursive(root->right, node_list, index);

  node_list[*index] = root;
  *index = *index + 1;
  return;
}

PLL_EXPORT unsigned int pll_stree_query_innernodes(pll_stree_t * root,
                                        pll_stree_t ** node_list)
{
  unsigned int index = 0;

  if (!root) return 0;
  if (!root->left) return 0;

  stree_query_innernodes_recursive(root->left,  node_list, &index);
  stree_query_innernodes_recursive(root->right, node_list, &index);

  node_list[index++] = root;

  return index;
}
#endif


