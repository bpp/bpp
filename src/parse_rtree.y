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
%{
#include "bpp.h"

extern int rtree_lex();
extern FILE * rtree_in;
extern void rtree_lex_destroy();
extern int rtree_lineno;
extern int rtree_colstart;
extern int rtree_colend;

extern int rtree_parse();
extern struct rtree_buffer_state * rtree__scan_string(const char * str);
extern void rtree__delete_buffer(struct rtree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

static void dealloc_species_data(species_t * data)
{
  if (data)
  {
    if (data->label)
      free(data->label);

    if (data->perloci_seqcount)
      free(data->perloci_seqcount);

    free(data);
  }
}

static void dealloc_gtree_data(gtree_data_t * data)
{
  if (data)
    free(data);
}

void rtree_graph_destroy(rnode_t * root)
{
  if (!root) return;

  rtree_graph_destroy(root->left);
  rtree_graph_destroy(root->right);

  dealloc_species_data(root->species_data);
  dealloc_gtree_data(root->gtree_data);

  free(root->label);
  free(root);
}

void rtree_destroy(rtree_t * tree)
{
  unsigned int i;
  rnode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];

    if (node->label)
      free(node->label);

    free(node);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static void rtree_error(rnode_t * node, const char * s)
{
  fatal("%s. (line %d column %d)", s, rtree_lineno, rtree_colstart);
}

%}

%union
{
  char * s;
  char * d;
  struct rnode_s * tree;
}

%error-verbose
%parse-param {struct rnode_s * tree}
%destructor { rtree_graph_destroy($$); } subtree
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
  free($7);

  tree->left->parent  = tree;
  tree->right->parent = tree;

};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$ = (rnode_t *)calloc(1, sizeof(rnode_t));
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  free($7);

  $$->left->parent  = $$;
  $$->right->parent = $$;

}
       | label optional_length
{
  $$ = (rnode_t *)calloc(1, sizeof(rnode_t));
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
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

static void fill_nodes_recursive(rnode_t * node,
                                 rnode_t ** array,
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

static unsigned int rtree_count_tips(rnode_t * root)
{
  unsigned int count = 0;

  if (root->left)
    count += rtree_count_tips(root->left);
  if (root->right)
    count += rtree_count_tips(root->right);

  if (!root->left && !root->right)
    return 1;

  return count;
}

rtree_t * rtree_wraptree(rnode_t * root,
                         unsigned int tip_count)
{
  unsigned int i;

  rtree_t * tree = (rtree_t *)xmalloc(sizeof(rtree_t));

  if (tip_count < 2 && tip_count != 0)
    fatal("Invalid number of tips in input tree (%u).", tip_count);

  if (tip_count == 0)
  {
    /* if tip counts is set to 0 then recursively count the number of tips */
    tip_count = rtree_count_tips(root);
    if (tip_count < 2)
    {
      fatal("Input tree contains no inner nodes.");
    }
  }

  tree->nodes = (rnode_t **)xmalloc((2*tip_count-1)*sizeof(rnode_t *));
  
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
    tree->nodes[i]->index = i;

  return tree;
}

rtree_t * rtree_parse_newick(const char * filename)
{
  rtree_t * tree;

  struct rnode_s * root;

  /* reset tip count */
  tip_cnt = 0;

  /* open input file */
  rtree_in = fopen(filename, "r");
  if (!rtree_in)
    fatal("Unable to open file (%s)", filename);

  /* create root node */
  root = (rnode_t *)xcalloc(1, sizeof(rnode_t));

  if (rtree_parse(root))
  {
    rtree_graph_destroy(root);
    root = NULL;
    fclose(rtree_in);
    rtree_lex_destroy();
    return NULL;
  }

  if (rtree_in) fclose(rtree_in);

  rtree_lex_destroy();

  /* wrap tree */
  tree = rtree_wraptree(root,tip_cnt);

  return tree;
}

rtree_t * rtree_parse_newick_string(const char * s)
{
  int rc;
  struct rnode_s * root;
  rtree_t * tree = NULL;

  /* reset tip count */
  tip_cnt = 0;

  root = (rnode_t *)xcalloc(1, sizeof(rnode_t));

  struct rtree_buffer_state * buffer = rtree__scan_string(s);
  rc = rtree_parse(root);
  rtree__delete_buffer(buffer);

  rtree_lex_destroy();

  if (!rc)
  {
    tree = rtree_wraptree(root,tip_cnt);
  }
  else
    free(root);

  return tree;
}
