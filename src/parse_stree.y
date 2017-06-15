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
%{
#include "bpp.h"

extern int stree_lex();
extern FILE * stree_in;
extern void stree_lex_destroy();
extern int stree_lineno;
extern int stree_colstart;
extern int stree_colend;

extern int stree_parse();
extern struct stree_buffer_state * stree__scan_string(const char * str);
extern void stree__delete_buffer(struct stree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

static void dealloc_data(snode_t * node,
                         void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

static void stree_graph_destroy(snode_t * root,
                                void (*cb_destroy)(void *))
{
  if (!root) return;

  stree_graph_destroy(root->left, cb_destroy);
  stree_graph_destroy(root->right, cb_destroy);

  dealloc_data(root, cb_destroy);

  free(root->label);
  free(root);
}

void stree_destroy(stree_t * tree,
                   void (*cb_destroy)(void *))
{
  unsigned int i,j;
  snode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];
    dealloc_data(node,cb_destroy);

    if (node->label)
      free(node->label);

    if (node->event)
    {
      for (j = 0; j < tree->locus_count; ++j)
        if (tree->nodes[i]->event[j])
        {
          dlist_clear(tree->nodes[i]->event[j],NULL);
          dlist_destroy(tree->nodes[i]->event[j]);
        }
      free(node->event);
    }
    
    if (node->event_count)
      free(node->event_count);

    if (node->seqin_count)
      free(node->seqin_count);

    if (node->logpr_contrib)
      free(node->logpr_contrib);

    if (node->old_logpr_contrib)
      free(node->old_logpr_contrib);

    if (node->gene_leaves)
      free(node->gene_leaves);

    free(node);
  }

  if (tree->pptable)
  {
    for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
      if (tree->pptable[i])
        free(tree->pptable[i]);
    free(tree->pptable);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static void stree_error(snode_t * node, const char * s)
{
  fatal("%s. (line %d column %d)", s, stree_lineno, stree_colstart);
}

%}

%union
{
  char * s;
  char * d;
  struct snode_s * tree;
}

%error-verbose
%parse-param {struct snode_s * tree}
%destructor { stree_graph_destroy($$,NULL); } subtree
%destructor { free($$); } STRING
%destructor { free($$); } NUMBER
%destructor { free($$); } label

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<tree> subtree
%start input
%%

input: OPAR subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  tree->left   = $2;
  tree->right  = $4;
  tree->label  = $6;
  tree->length = $7 ? atof($7) : 0;
  tree->parent = NULL;
  tree->leaves = $2->leaves + $4->leaves;
  free($7);

  tree->left->parent  = tree;
  tree->right->parent = tree;

};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$ = (snode_t *)calloc(1, sizeof(snode_t));
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  $$->leaves = $2->leaves + $4->leaves;
  free($7);

  $$->left->parent  = $$;
  $$->right->parent = $$;

}
       | label optional_length
{
  $$ = (snode_t *)calloc(1, sizeof(snode_t));
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->leaves = 1;
  $$->left   = NULL;
  $$->right  = NULL;
  tip_cnt++;
  free($2);
};


optional_label:  {$$ = NULL;} | label  {$$ = $1;};
optional_length: {$$ = NULL;} | COLON number {$$ = $2;};
label: STRING    {$$=$1;} | NUMBER {$$=$1;};
number: NUMBER   {$$=$1;};

%%

static void fill_nodes_recursive(snode_t * node,
                                 snode_t ** array,
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

static unsigned int stree_count_tips(snode_t * root)
{
  unsigned int count = 0;

  if (root->left)
    count += stree_count_tips(root->left);
  if (root->right)
    count += stree_count_tips(root->right);

  if (!root->left && !root->right)
    return 1;

  return count;
}

stree_t * stree_wraptree(snode_t * root,
                         unsigned int tip_count)
{
  unsigned int i;

  stree_t * tree = (stree_t *)xmalloc(sizeof(stree_t));

  if (tip_count < 2 && tip_count != 0)
    fatal("Invalid number of tips in input tree (%u).", tip_count);

  if (tip_count == 0)
  {
    /* if tip counts is set to 0 then recursively count the number of tips */
    tip_count = stree_count_tips(root);
    if (tip_count < 2)
    {
      fatal("Input tree contains no inner nodes.");
    }
  }

  tree->nodes = (snode_t **)xmalloc((2*tip_count-1)*sizeof(snode_t *));
  
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

stree_t * stree_parse_newick(const char * filename)
{
  stree_t * tree;

  struct snode_s * root;

  /* reset tip count */
  tip_cnt = 0;

  /* open input file */
  stree_in = fopen(filename, "r");
  if (!stree_in)
    fatal("Unable to open file (%s)", filename);

  /* create root node */
  root = (snode_t *)xcalloc(1, sizeof(snode_t));

  if (stree_parse(root))
  {
    stree_graph_destroy(root,NULL);
    root = NULL;
    fclose(stree_in);
    stree_lex_destroy();
    return NULL;
  }

  if (stree_in) fclose(stree_in);

  stree_lex_destroy();

  /* wrap tree */
  tree = stree_wraptree(root,tip_cnt);

  return tree;
}

stree_t * stree_parse_newick_string(const char * s)
{
  int rc;
  struct snode_s * root;
  stree_t * tree = NULL;

  /* reset tip count */
  tip_cnt = 0;

  root = (snode_t *)xcalloc(1, sizeof(snode_t));

  struct stree_buffer_state * buffer = stree__scan_string(s);
  rc = stree_parse(root);
  stree__delete_buffer(buffer);

  stree_lex_destroy();

  if (!rc)
  {
    tree = stree_wraptree(root,tip_cnt);
  }
  else
    free(root);

  return tree;
}
