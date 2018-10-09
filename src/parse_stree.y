/*
    Copyright (C) 2016-2018 Tomas Flouri, Bruce Rannala and Ziheng Yang

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
static unsigned int hybrid_cnt = 0;
static unsigned int inner_cnt = 0;

static long get_double(const char * line, double * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%lf%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

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

  if (root->label)
    free(root->label);
  free(root);
}

void stree_destroy(stree_t * tree,
                   void (*cb_destroy)(void *))
{
  unsigned int i,j;
  snode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count + tree->hybrid_count; ++i)
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

    if (node->t2h)
      free(node->t2h);

    if (node->old_t2h)
      free(node->old_t2h);

    free(node);
  }

  /* safety check */
  assert(opt_network == !!tree->hybrid_count);

  if (tree->pptable)
  {
    for (i = 0; i < tree->tip_count + tree->inner_count + tree->hybrid_count; ++i)
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

static void parse_annotation(snode_t * snode, const char * annotation)
{
  double val;

  char * p = xstrdup(annotation);
  char * s = p;

  /* s now (possibly) has the format "token1,token2,...,tokenN"  where
     tokenX has the format option=value or &option=value */

  while (*s)
  {
    size_t comma = 0;
    /* find next comma and replace it with a 0 to create a string of the
       current token*/
    size_t tokenlen = strcspn(s,",");
    if (tokenlen)
    {
      if (s[tokenlen] == ',')
        comma = 1;
      s[tokenlen] = 0;
    }
    else
    {
      if (s == p)
        fatal("Error - annotation (%s) starts with a comma", annotation);
      else
        fatal("Consecutive comma symbols found in annotation (%s)", annotation);
    }
    
    /* find equal sign in current string */
    size_t optlen = strcspn(s,"=");
    if (!optlen)
      fatal("Erroneous format in annotation (%s) - token (%s) is missing a '=' "
            "sign and a value", annotation, s);

    if (strlen(s+optlen) == 0)
      fatal("Missing value for option (%s) in annotation (%s)", s, annotation);

    if (s[0] == '&')
    {
      ++s;
      --tokenlen;
      --optlen;
    }
    s[optlen] = 0;

    if (!strcasecmp(s,"gamma"))
    {
      if (!get_double(s+optlen+1,&val))
        fatal("Cannot parse value (%s) for token (%s) in annotation (%s) - "
              "value must be of type float", s+optlen, s, annotation);
      
      if (val < 0 || val > 1)
        fatal("Parameter gamma in annotations must be greater than 0 and "
              "smaller than 1");
      
      snode->hgamma = val;
    }
    else if (!strcasecmp(s,"tau-parent"))
    {
      /* temporarily set htau of parents node to the corresponding hybridization
         node. Later once the tree has been parsed and annotated, a separate
         loop will reset htau such that hybridization nodes get an htau=1 and
         their parents have their htau set correctly */
      if (!strcasecmp(s+optlen+1,"yes"))
        snode->htau = 1;
      else if (!strcasecmp(s+optlen+1, "no"))
        snode->htau = 0;
      else
        fatal("Invalid value (%s) for option (%s) in annotation (%s)",
              s+optlen+1, s, annotation);
    }
    else
      fatal("Invalid token (%s) in annotation (%s)", s, annotation);
    
    s += tokenlen+comma;
  }

  free(p);
  
}

static void annotate_bd_introgression(snode_t * xtip,
                                      snode_t * xinner,
                                      snode_t * y,
                                      snode_t * ylink)
{
  /* y is the parent of xtip, and ylink is the child of xinner

     W.l.o.g. the idea for annotations is the following:
       Annotation of xtip relates to the edge Y->X
       Annotation of xinner relates to R->X (R is parent of X)
       Annotation of y relates to P->Y (P is parent of Y)
       Annotation of ylink relates to X->Y

     Also, nodes xtip and ylink are to be deleted later on, therefore
     we do not annotate them

  */

  /* better safe than sorry */
  assert(!xtip->left && !xtip->right);          /* xtip is tip */
  assert(xinner->left || xinner->right);        /* xinner is inner */
  assert(y->left || y->right);                  /* y is inner */
  assert(ylink->left || ylink->right);          /* ylink is inner */
  
  if (xtip->data)
  {
    parse_annotation(xtip, (char *)(xtip->data));
    free(xtip->data);
    xtip->data = NULL;
  }
  if (xinner->data)
  {
    parse_annotation(xinner, (char *)(xinner->data));
    free(xinner->data);
    xinner->data = NULL;
  }
  if (y->data)
  {
    parse_annotation(y, (char *)(y->data));
    free(y->data);
    y->data = NULL;
  }
  if (ylink->data)
  {
    parse_annotation(ylink, (char *)(ylink->data));
    free(ylink->data);
    ylink->data = NULL;
  }

  if (xtip->htau == 1 || xinner->htau || y->htau || ylink->htau)
    fatal("Bidirectional introgression between (%s) and (%s) requires exactly "
          "one tau parameter shared by both nodes. Please remove *all* 'tau' "
          "annotations for these two nodes.",
          xtip->label, y->label);

  if (xtip->hgamma && xinner->hgamma)
  {
    if (xtip->hgamma + xinner->hgamma != 1)
      fatal("Gamma parameter annotations for bidirectional introgression event "
            "on node (%s) do not sum to 1",  xtip->label);
  }
  else if (xtip->hgamma)
  {
    /* gamma on edge Y->X */

    xinner->hybrid->hgamma = xtip->hgamma;
    xinner->hgamma = 1 - xtip->hgamma;

  }
  else if (xinner->hgamma)
    xinner->hybrid->hgamma = 1 - xinner->hgamma;
  else
    fatal("Missing gamma parameter annotation for bidirectional introgression "
          "event on node (%s)", xtip->label);

  if (y->hgamma && ylink->hgamma)
  {
    if (y->hgamma + ylink->hgamma != 1)
      fatal("Gamma parameter annotations for bidirectional introgression event "
            "on node (%s) do not sum to 1",  y->label);
  }
  else if (y->hgamma)
  {
    y->hybrid->hgamma = 1 - y->hgamma;
  }
  else if (ylink->hgamma)
  {
    y->hgamma = 1 - ylink->hgamma;
    y->hybrid->hgamma = ylink->hgamma;
  }
  else
    fatal("Missing gamma parameter annotation for bidirectional introgression "
          "event on node (%s)", y->label);
}
static void annotate_hybridization(snode_t * hinner, snode_t * htip)
{
  /* better safe than sorry */
  assert(hinner->parent);
  assert(hinner->hybrid->parent);

  if (hinner->parent == hinner->hybrid->parent)
  {
    /* parallel edges from parent - in this case parent must have a tau. For now
       we set to hinner/hinner->hybrid what we want the parent's htau to be */
    hinner->htau = hinner->hybrid->htau = 1;
  }

  if (hinner->data)
  {
    parse_annotation(hinner, (char *)(hinner->data));
    free(hinner->data);
    hinner->data = NULL;
  }

  if (htip->data)
  {
    parse_annotation(hinner->hybrid, (char *)(htip->data));
    free(htip->data);
    htip->data = NULL;
  }

  if (hinner->parent == hinner->hybrid->parent)
  {
    /* parallel edges from parent - in this case parent must have a tau.
       We check whether the htau for one of the two nodes has changed due to
       annotation. If it did, then we show an error */
    if (hinner->htau != 1 || hinner->hybrid->htau != 1)
      fatal("Hybridization event (%s) requires that the parent has a tau. "
            "Please remove any tau-parent annotation from the node",
            hinner->label);
   
  }

  if (hinner->hgamma && hinner->hybrid->hgamma)
  {
    if (hinner->hgamma + hinner->hybrid->hgamma != 1)
      fatal("Gamma parameter annotations for hybridization event (%s) do not "
            "sum to 1", hinner->label);
  }
  else if (hinner->hgamma)
  {
    hinner->hybrid->hgamma = 1 - hinner->hgamma;
    hinner->parent->hgamma = hinner->hgamma;
    hinner->hybrid->parent->hgamma = hinner->hybrid->hgamma;
  }
  else if (hinner->hybrid->hgamma)
  {
    hinner->hgamma = 1 - hinner->hybrid->hgamma;
    hinner->parent->hgamma = hinner->hgamma;
    hinner->hybrid->parent->hgamma = hinner->hybrid->hgamma;
  }
  else
  {
    fatal("Missing gamma parameter annotation for hybridization event (%s)",
          hinner->label);
  }

  hinner->parent->htau = hinner->htau;
  hinner->hybrid->parent->htau = hinner->hybrid->htau;
}

static void resolve_hybridization(stree_t * stree, long * dups)
{
  long i,j;
  long hybridization = 0;

  /* check for hybridization events 

          * R          ((A,(C)H)S,(H,B)T)R;
         / \
      S /   \          Method:
       *     \           1. Find a tip (htip) for which an inner node H with
      / \  H  \  T          identical label exists (and name it hinner) 
     /   *     *         2. Create hybrid nodes for hinner
    /    |     |\        3. Link the corresponding nodes to obtain the diagram
   /     |     | \          below        
   A     C     H  B

   Results in:
  
         R *
          / \
         /   \
     S  *     *  T
       / \   / \
      /   \ /   \
     /     * H   \
    /      |      \
   A       C       B  
  
  */

  for (i = 0; i < stree->tip_count; ++i)
  {
    if (!stree->nodes[i] || dups[i] == i) continue;

    /* if the duplicate of tip is a tip, ignore */
    if (dups[i] < stree->tip_count) continue;
    
    /* root cannot be duplicated */
    assert(stree->nodes[i]->parent);

    /* pointer to tip H node */
    snode_t * htip = stree->nodes[i];

    /* pointer to inner H node */
    snode_t * hinner = stree->nodes[dups[i]];

    /* ensure hinner is indeed an inner node and has a parent */ 
    if (!hinner->parent) continue;
    if (!hinner->left && !hinner->right) continue;

    hybrid_cnt++;
    hybridization = 1;

    hinner->hybrid = (snode_t *)xcalloc(1,sizeof(snode_t));
    hinner->hybrid->hybrid = hinner;
    assert(hinner->label);
    hinner->hybrid->label = xstrdup(hinner->label);

    if (htip->parent->left == htip)
      htip->parent->left = hinner->hybrid;
    else
      htip->parent->right = hinner->hybrid;
    hinner->hybrid->parent = htip->parent;
    
    /* check annotations */
    annotate_hybridization(hinner, htip);

    htip->parent = NULL;
    stree->nodes[htip->node_index] = NULL;
    stree_graph_destroy(htip,NULL);
    tip_cnt--;
  }

  /* rebuild stree->nodes array */
  if (hybridization)
  {
    snode_t ** tmp = (snode_t **)xcalloc((size_t)(tip_cnt+inner_cnt+hybrid_cnt),
                                         sizeof(snode_t *));
    for (i=0,j=0; i < stree->tip_count+stree->inner_count; ++i)
      if (stree->nodes[i])
        tmp[j++] = stree->nodes[i];
    
    assert(j == tip_cnt+inner_cnt);

    for (i=tip_cnt; i < tip_cnt+inner_cnt; ++i)
      if (tmp[i]->hybrid)
        tmp[j++] = tmp[i]->hybrid;

    assert (j == tip_cnt + inner_cnt + hybrid_cnt);

    free(stree->nodes);
    stree->nodes = tmp;
    for (i = 0; i < tip_cnt + inner_cnt + hybrid_cnt; ++i)
      stree->nodes[i]->node_index = i;

    stree->tip_count = tip_cnt;
    stree->inner_count = inner_cnt;
    stree->hybrid_count = hybrid_cnt;
    stree->edge_count = tip_cnt + inner_cnt + hybrid_cnt - 1;
  }
}

static void resolve_bd_introgression(stree_t * stree, long * dups)
{
  long i,j;
  snode_t * link1 = NULL;
  snode_t * link2 = NULL;

  /* check for bidirectional introgression 

          *           ((A,(B)Y)X,(X,B)Y)R;
         / \
      X /   \          Method:
       *     \           1. Find a tip (link1) for which an inner node with
      / \  Y  \  Y          identical label exists (in this case X) 
     /   *     *         2. Get the parent of tip node X (in this case node Y)
    /    |     |\        3. Check that inner node X has an unary child with the
   /     |     | \          same label as Y (name it link2)                           
   A     B     B  X      4. Create hybrid nodes for both X and Y                 
                         5. Link the corresponding nodes to obtain the diagram   
                            below
   
   Results in:

         R
        /\
       /  \
   X  *<-->*  Y
     /      \
    /        \
   A          B
              
  */

  /* check for bidirectional introgression */
  for (i = 0; i < stree->tip_count; ++i)
  {
    /* if no duplicate label for current tip node, move to next tip */
    if (!stree->nodes[i] || dups[i] == i) continue;

    /* if the duplicate of tip is a tip, ignore */
    if (dups[i] < stree->tip_count) continue;

    /* root cannot be duplicated */
    assert(stree->nodes[i]->parent);

    /* get pointers to the two inner nodes which are candidates for a
       bidirectional introgression */
    snode_t * x = stree->nodes[dups[i]];
    snode_t * y = stree->nodes[i]->parent;
    if (!y->label) continue;
    link1 = stree->nodes[i];

    assert(x->left || x->right);      /* confirm x is an inner node */
    assert(x->label);


    if (x->left && x->left->label && !strcmp(x->left->label,y->label))
      link2 = x->left;
    else if (x->right && x->right->label && !strcmp(x->right->label,y->label))
      link2 = x->right;

    if (link2)
    {
      /* link2 must be an unary node */
      assert((link2->left && !link2->right) || (!link2->left && link2->right));

      hybrid_cnt += 2;
      tip_cnt -= 2;
      inner_cnt--;

      /* we create an additional node for each introgression event */
      x->hybrid = (snode_t *)xcalloc(1,sizeof(snode_t));
      y->hybrid = (snode_t *)xcalloc(1,sizeof(snode_t));
      assert(x->label);
      assert(y->label);
      x->hybrid->label = xstrdup(x->label);
      y->hybrid->label = xstrdup(y->label);

      /* link hybridization nodes with their corresponding nodes */
      x->hybrid->hybrid = x;
      y->hybrid->hybrid = y;

      /* define parents of hybridization nodes */
      x->hybrid->parent = y->hybrid;
      y->hybrid->parent = x->hybrid;

      /* define one child (as left child node) for hybridization nodes */
      x->hybrid->left = y->hybrid;
      y->hybrid->left = x->hybrid;

      x->hybrid->leaves = y->leaves;
      y->hybrid->leaves = x->leaves;

      /* annotate */
      annotate_bd_introgression(link1, x, y, link2);

      /* now destroy redundant nodes/subtrees (i.e. link1 and link2) */
      if (link1->parent->left == link1)
      {
        link1->parent->left = link1->parent->right;
        link1->parent->right = NULL;
      }
      else
      {
        assert(link1 == link1->parent->right);
        link1->parent->right = NULL;
      }
      link1->parent = NULL;
      stree->nodes[link1->node_index] = NULL;
      stree_graph_destroy(link1,NULL);

      /* now the same for link2 */
      if (link2->parent->left == link2)
      {
        link2->parent->left = link2->parent->right;
        link2->parent->right = NULL;
      }
      else
      {
        assert(link2 == link2->parent->right);
        link2->parent->right = NULL;
      }
      link2->parent = NULL;
      stree->nodes[link2->node_index] = NULL;
      if (link2->left)
        stree->nodes[link2->left->node_index] = NULL;
      else
        stree->nodes[link2->right->node_index] = NULL;
      stree_graph_destroy(link2,NULL);
    }
  }

  /* rebuild stree->nodes array */
  if (hybrid_cnt)
  {
    snode_t ** tmp = (snode_t **)xcalloc((size_t)(tip_cnt+inner_cnt+hybrid_cnt),
                                         sizeof(snode_t *));
    for (i=0,j=0; i < stree->tip_count+stree->inner_count; ++i)
      if (stree->nodes[i])
        tmp[j++] = stree->nodes[i];
    
    assert(j == tip_cnt+inner_cnt);

    for (i=tip_cnt; i < tip_cnt+inner_cnt; ++i)
      if (tmp[i]->hybrid)
        tmp[j++] = tmp[i]->hybrid;

    assert (j == tip_cnt + inner_cnt + hybrid_cnt);

    free(stree->nodes);
    stree->nodes = tmp;
    for (i = 0; i < tip_cnt + inner_cnt + hybrid_cnt; ++i)
      stree->nodes[i]->node_index = i;

    stree->tip_count = tip_cnt;
    stree->inner_count = inner_cnt;
    stree->hybrid_count = hybrid_cnt;
    stree->edge_count = tip_cnt + inner_cnt + hybrid_cnt - 1;
  }
}

static void validate_stree(stree_t * stree)
{
  long i,j;

  /* this routine is called after hybridizations are resolved. It ensures that
     the speecies tree (now network):
     1. does not contain nodes with the same label
     2. No two unary nodes in relation parent-child (i.e. string of unaries)
  */

  long nodes_count = stree->tip_count + stree->inner_count; // + stree->hybrid_count;

  for (i = 0; i < nodes_count; ++i)
  {
    
    /* check if node is unary *and* has a parent */
    if (stree->nodes[i]->parent &&
        ((stree->nodes[i]->left && !stree->nodes[i]->right) ||
        (!stree->nodes[i]->left && stree->nodes[i]->right)))
    {
      /* if it's parent is also unary then we have two unary nodes in a row */
      if ((stree->nodes[i]->parent->left && !stree->nodes[i]->parent->right) ||
          (!stree->nodes[i]->parent->left && stree->nodes[i]->parent->right))
        fatal("Tree contains consecutive unary nodes");
    }

    if (!stree->nodes[i]->label) continue;
    /* 1. check for duplicate labels */
    for (j = i+1; j < nodes_count; ++j)
    {
      if (!stree->nodes[j]->label) continue;
      
      if (!strcmp(stree->nodes[i]->label,stree->nodes[j]->label))
        fatal("Tree has nodes with same label (%s)", stree->nodes[i]->label);
    }
  }
}

/* creates array dups, such that dups[i] is the index of a node with identical
   label to node stree->nodes[i]. If more than one duplicate exists for a given
   node, the procedure ends with an error. If no duplicate exists for a node,
   then dup[i] = i */
static long * create_duplicate_indices(stree_t * stree)
{
  long i,j;
  long freq = 0;

  long * dups = (long *)xmalloc((size_t)(stree->tip_count+stree->inner_count) *
                                sizeof(long));

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    dups[i] = i;

  /* check that a node does not appear more than twice */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (!stree->nodes[i]->label) continue;

    freq = 1;

    for (j = i+1; j < stree->tip_count + stree->inner_count; ++j)
    {
      if (!stree->nodes[j]->label) continue;
      if (!strcmp(stree->nodes[i]->label,stree->nodes[j]->label))
      {
        ++freq;
        dups[i] = j;
        dups[j] = i;
      }
    }

    if (freq > 2)
      fatal("Node with label %s appears more %ld times (allowed twice)",
            stree->nodes[i]->label, freq);
  }
  return dups;
}

static void resolve_network(stree_t * stree)
{
  long * dups;

  dups = create_duplicate_indices(stree);
  resolve_bd_introgression(stree,dups);
  free(dups);

  dups = create_duplicate_indices(stree);
  resolve_hybridization(stree,dups);
  free(dups);

  validate_stree(stree);
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
%token OBRA
%token CBRA
%token COMMA
%token COLON SEMICOLON
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label optional_annotation string_list
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

  inner_cnt++;
}
       | label SEMICOLON
{
  tree->left  = NULL;
  tree->right = NULL;
  tree->label = $1;
  tree->length = 0;
  tree->parent = NULL;
  tree->leaves = 1;

  tip_cnt++;
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

  inner_cnt++;
}
       | OPAR subtree CPAR label optional_annotation optional_length
{
  $$ = (snode_t *)calloc(1, sizeof(snode_t));
  $$->left   = $2;
  $$->right  = NULL;
  $$->label  = $4;
  $$->length = $6 ? atof($6) : 0;
  $$->leaves = $2->leaves;
  free($6);

  if ($5) 
    $$->data = (void *)($5);

  $$->left->parent = $$;
  inner_cnt++;
}
       | label optional_annotation optional_length
{
  $$ = (snode_t *)calloc(1, sizeof(snode_t));
  $$->label  = $1;
  $$->length = $3 ? atof($3) : 0;
  $$->leaves = 1;
  $$->left   = NULL;
  $$->right  = NULL;
  tip_cnt++;
  free($3);

  if ($2)
    $$->data = (void *)($2);
};

optional_annotation: {$$ = NULL;} | OBRA string_list CBRA {$$ = $2;};
string_list: STRING 
{
        $$ = $1;
}
       | STRING COMMA string_list
{
        $$ = (char *)xcalloc((strlen($1)+strlen($3)+2),sizeof(char));
        strcpy($$,$1);
        $$[strlen($1)] = ',';
        strcpy($$+strlen($1)+1,$3);
        free($1);
        free($3);
};
optional_label:  {$$ = NULL;} | label  {$$ = $1;};
optional_length: {$$ = NULL;} | COLON number {$$ = $2;};
label: STRING    {$$=$1;} | NUMBER {$$=$1;};
number: NUMBER   {$$=$1;};

%%

/* fill array in preorder */
static void fill_nodes_recursive(snode_t * node,
                                 snode_t ** array,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!node->left && !node->right)
  {
    array[*tip_index] = node;
    *tip_index = *tip_index + 1;
    return;
  }

  array[*inner_index] = node;
  *inner_index = *inner_index + 1;

  if (node->left)
    fill_nodes_recursive(node->left,  array, tip_index, inner_index);
  if (node->right)
    fill_nodes_recursive(node->right, array, tip_index, inner_index);
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

static void reorder(stree_t * stree)
{
  unsigned int i,j,k;
  unsigned int commas_count = 0;
  pair_t ** pairlist;

  /* compute number of commas in list of tips */
  for (i = 0; i < strlen(opt_reorder); ++i)
    if (opt_reorder[i] == ',')
      commas_count++;

  if (commas_count+1 != stree->tip_count)
    fatal("Labels (%d) specified in --reorder do not match species tree", commas_count);

  hashtable_t * ht = hashtable_create(stree->tip_count);

  pairlist = (pair_t **)xmalloc(stree->tip_count * sizeof(pair_t *));

  for (i = 0; i < stree->tip_count; ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = stree->nodes[i]->label;
    pair->data = (void *)(uintptr_t)i;
    pairlist[i] = pair;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(stree->nodes[i]->label),
                          cb_cmp_pairlabel))
      fatal("Duplicate taxon (%s)", stree->nodes[i]->label);
  }

  char * s = opt_reorder;

  k = 0;
  while (*s)
  {
    /* get next tip */
    size_t taxon_len = strcspn(s,",");
    if (!taxon_len)
      fatal("Erroneous format in --reorder (taxon missing)");

    char * taxon = xstrndup(s, taxon_len);

    /* search tip in hash table */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    cb_cmp_pairlabel);
    if (!query)
      fatal("Taxon %s does not appear in the tree", taxon);

    j = (unsigned int)(uintptr_t)(query->data);
    /* swap */
    if (j != k)
    {
      i = k;

      assert((unsigned int)(uintptr_t)(pairlist[i]->data) == i);
      assert(pairlist[j] == query);

      SWAP(stree->nodes[i],stree->nodes[j]);
      SWAP(pairlist[i],pairlist[j]);
      SWAP(pairlist[k]->data,pairlist[j]->data);
    }

    free(taxon);

    s += taxon_len;
    if (*s == ',')
      s += 1;

    ++k;
  }
  /* kill hashtable */
  hashtable_destroy(ht,free);

  free(pairlist);
}

stree_t * stree_wraptree(snode_t * root)
{
  unsigned int i;

  if (!root)
    fatal("No node found in tree");

  stree_t * stree = (stree_t *)xcalloc(1,sizeof(stree_t));
  
  unsigned int tip_count = stree_count_tips(root);
  assert(tip_count == tip_cnt);

  if (tip_count == 0)
    fatal("Invalid number of tips in input tree (%u).", tip_count);

  #if 0
  if (tip_count == 0)
  {
    /* if tip counts is set to 0 then recursively count the number of tips */
    /* TODO: Check difference with master version. Here stree_count_tips is
       always called at the beginning */
  }
  #endif

  stree->nodes = (snode_t **)xmalloc((tip_cnt+inner_cnt)*sizeof(snode_t *));
  
  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;

  /* fill tree->nodes in pre-order */
  fill_nodes_recursive(root, stree->nodes, &tip_index, &inner_index);

  stree->tip_count = tip_count;
  stree->edge_count = tip_count + inner_cnt + hybrid_cnt - 1;
  stree->inner_count = inner_cnt;
  stree->hybrid_count = hybrid_cnt;
  stree->root = root;
  stree->pptable = NULL;

  for (i = 0; i < tip_count + inner_cnt; ++i)
    stree->nodes[i]->node_index = i;

  resolve_network(stree);
  if (stree->hybrid_count)
    opt_network = 1;

  /* reorder tip nodes if specified */
  if (opt_reorder)
    reorder(stree);

  for (i = 0; i < tip_cnt + inner_cnt + hybrid_cnt; ++i)
    stree->nodes[i]->node_index = i;

  /* apply diploid information */
  if (opt_diploid)
  {
    if (opt_diploid_size != stree->tip_count)
      fatal("Number of 'diploid' assignments mismatch number of species");

    for (i = 0; i < stree->tip_count; ++i)
      stree->nodes[i]->diploid = opt_diploid[i];
  }

  if (opt_network)
  {
    /* for all nodes that are not hybridization events set htau to 1 */
    for (i = 0; i < stree->inner_count; ++i)
    {
      snode_t * snode = stree->nodes[stree->tip_count+i];

      if (!snode->hybrid) snode->htau = 1;
    }

    /* now reset the 'htau' values for parents of hybridization events, and set
       htau to hybridization nodes to 1 (inner nodes) and 0 (mirrored nodes) */
    for (i = 0; i < stree->hybrid_count; ++i)
    {
      snode_t * mnode = stree->nodes[stree->tip_count+stree->inner_count+i];
      snode_t * hnode = mnode->hybrid;
      assert(hnode);

      /* bidirection  */
      if (node_is_bidirection(hnode))
      {
        assert(hnode->node_index < stree->tip_count + stree->inner_count);

        /* if htau=0 for all four nodes involved in a bidirectional
           introgression then set htau=1 to inner node hnode. All other
           nodes keep htau=0 forever, and will share the tau of hnode */
        if (!hnode->htau && !mnode->htau && 
            !mnode->parent->htau && !mnode->parent->hybrid->htau)
          hnode->htau = 1;
      }
      else
      {
        hnode->parent->htau = hnode->htau;
        mnode->parent->htau = mnode->htau;

        hnode->htau = 1;
        mnode->htau = 0;

      }
    }
  }

  if (opt_network && stree->root->htau == 0)
    fatal("Error: species tree root requires a tau parameter [tau-parent=yes]");

  return stree;
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
  tree = stree_wraptree(root);

  return tree;
}

stree_t * stree_parse_newick_string(const char * s)
{
  int rc;
  struct snode_s * root;
  stree_t * tree = NULL;

  /* reset counts */
  tip_cnt = 0;
  hybrid_cnt = 0;
  inner_cnt = 0;

  root = (snode_t *)xcalloc(1, sizeof(snode_t));

  struct stree_buffer_state * buffer = stree__scan_string(s);
  rc = stree_parse(root);
  stree__delete_buffer(buffer);

  stree_lex_destroy();

  if (!rc)
  {
    tree = stree_wraptree(root);
  }
  else
    free(root);

  return tree;
}
