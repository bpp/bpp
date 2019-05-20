/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         stree_parse
#define yylex           stree_lex
#define yyerror         stree_error
#define yylval          stree_lval
#define yychar          stree_char
#define yydebug         stree_debug
#define yynerrs         stree_nerrs

/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 22 "parse_stree.y"

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

    if (node->hx)
      free(node->hx);

    if (opt_clock != BPP_CLOCK_GLOBAL)
      free(node->brate);

    if (node->mark)
      free(node->mark);

    free(node);
  }

  /* safety check */
  assert(opt_msci == !!tree->hybrid_count);

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

    if (!strcasecmp(s,"phi"))
    {
      if (!get_double(s+optlen+1,&val))
        fatal("Cannot parse value (%s) for token (%s) in annotation (%s) - "
              "value must be of type float", s+optlen, s, annotation);
      
      if (val < 0 || val > 1)
        fatal("Parameter phi in annotations must be greater than 0 and "
              "smaller than 1");
      
      snode->hphi = val;
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
  assert(xinner->left && xinner->right);        /* xinner is inner */
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

  if (xtip->hphi && xinner->hphi)
  {
    if (xtip->hphi + xinner->hphi != 1)
      fatal("Gamma parameter annotations for bidirectional introgression event "
            "on node (%s) do not sum to 1",  xtip->label);
  }
  else if (xtip->hphi)
  {
    /* phi on edge Y->X */

    xinner->hybrid->hphi = xtip->hphi;
    xinner->hphi = 1 - xtip->hphi;

  }
  else if (xinner->hphi)
    xinner->hybrid->hphi = 1 - xinner->hphi;
  else
    fatal("Missing phi parameter annotation for bidirectional introgression "
          "event on node (%s)", xtip->label);

  if (y->hphi && ylink->hphi)
  {
    if (y->hphi + ylink->hphi != 1)
      fatal("Gamma parameter annotations for bidirectional introgression event "
            "on node (%s) do not sum to 1",  y->label);
  }
  else if (y->hphi)
  {
    y->hybrid->hphi = 1 - y->hphi;
  }
  else if (ylink->hphi)
  {
    y->hphi = 1 - ylink->hphi;
    y->hybrid->hphi = ylink->hphi;
  }
  else
    fatal("Missing phi parameter annotation for bidirectional introgression "
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

  if (hinner->hphi && hinner->hybrid->hphi)
  {
    if (hinner->hphi + hinner->hybrid->hphi != 1)
      fatal("Gamma parameter annotations for hybridization event (%s) do not "
            "sum to 1", hinner->label);
  }
  else if (hinner->hphi)
  {
    hinner->hybrid->hphi = 1 - hinner->hphi;
    hinner->parent->hphi = hinner->hphi;
    hinner->hybrid->parent->hphi = hinner->hybrid->hphi;
  }
  else if (hinner->hybrid->hphi)
  {
    hinner->hphi = 1 - hinner->hybrid->hphi;
    hinner->parent->hphi = hinner->hphi;
    hinner->hybrid->parent->hphi = hinner->hybrid->hphi;
  }
  else
  {
    fatal("Missing phi parameter annotation for hybridization event (%s)",
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

    /* if simulation, copy branch length and theta */
    if (opt_simulate)
    {
      hinner->hybrid->theta  = htip->theta;
      hinner->hybrid->length = htip->length;
    }

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

          *  R        ((A,(B)Y)X,(X)Y)R;
         / \
      X /   \          Method:
       *     \           1. Find a tip (link1) for which an inner node with
      / \  Y  \  Y          identical label exists (in this case X) 
     /   *     *         2. Get the parent of tip node X (in this case node Y)
    /    |     |         3. Check that inner node X has an unary child with the
   /     |     |           same label as Y (name it link2)                           
   A     B     X         4. Create hybrid nodes for both X and Y                 
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

    /* if the duplicate of tip is a tip, fail as we checked this in
       create_duplicate_indices */
    assert(dups[i] >= stree->tip_count);

    /* root cannot be duplicated */
    assert(stree->nodes[i]->parent);

    /* get pointers to the two inner nodes which are candidates for a
       bidirectional introgression */
    snode_t * x = stree->nodes[dups[i]];        /* inner node x */
    snode_t * y = stree->nodes[i]->parent;      /* parent of tip node x */
    if (!y->label) continue;
    link1 = stree->nodes[i];                    /* tip node x */

    assert(x->left || x->right);      /* confirm x is an inner node */
    assert(x->label);


    if (x->left && x->left->label && !strcmp(x->left->label,y->label))
      link2 = x->left;
    else if (x->right && x->right->label && !strcmp(x->right->label,y->label))
      link2 = x->right;

    if (link2)
    {
      /* link2 must be an unary node */
      assert(link2->left && !link2->right);

      hybrid_cnt += 2;
      tip_cnt -= 1;
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

      /* define parents of (mirror) hybridization nodes */
      x->hybrid->parent = y;
      y->hybrid->parent = x;

      /* define one child (as left child node) for hybridization nodes */
      x->hybrid->left  = NULL;
      y->hybrid->left  = NULL;
      x->hybrid->right = NULL;
      y->hybrid->right = NULL;

      /* if simulation, copy branch length and theta */
      if (opt_simulate)
      {
        x->hybrid->tau   = link1->tau;
        x->hybrid->theta = link1->theta;
      }

      /* annotate */
      annotate_bd_introgression(link1, x, y, link2);

      /* now destroy redundant nodes/subtrees (i.e. link1 and link2) */
      assert(y == link1->parent);
      assert(y->left == link1 && y->right == NULL);
      y->left = NULL;
      link1->parent = NULL;
      stree->nodes[link1->node_index] = NULL;
      stree_graph_destroy(link1,NULL);

      /* link y with child of link2 */
      y->left = link2->left;
      y->left->parent = y;

      /* now destroy link2 */
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
      assert(link2->left && link2->right == NULL);
      link2->left = NULL;

      /* if simulation, copy branch length and theta */
      if (opt_simulate)
      {
        y->hybrid->tau   = link2->tau;
        y->hybrid->theta = link2->theta;
      }

      stree_graph_destroy(link2,NULL);

      x->right = y->hybrid;
      y->right = x->hybrid;

      y->leaves = y->left->leaves;
      x->leaves = x->left->leaves;
      x->hybrid->leaves = 0;
      y->hybrid->leaves = 0;


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
   node, the procedure ends with an error. If a duplicate tip label is found,
   the procedure ends with error. If no duplicate exists for a node, or a non-tip
   duplicate exists, then dup[i] = i. */
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
    else if (freq == 2)
    {
      /* if both are tips */
      if (i < stree->tip_count && dups[i] < stree->tip_count)
        fatal("Tip node with label %s appears more than once in the tree");
    }
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


/* Line 371 of yacc.c  */
#line 870 "parse_stree.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parse_stree.h".  */
#ifndef YY_STREE_PARSE_STREE_H_INCLUDED
# define YY_STREE_PARSE_STREE_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int stree_debug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     OPAR = 258,
     CPAR = 259,
     OBRA = 260,
     CBRA = 261,
     COMMA = 262,
     HASH = 263,
     COLON = 264,
     SEMICOLON = 265,
     STRING = 266,
     NUMBER = 267
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 387 of yacc.c  */
#line 818 "parse_stree.y"

  char * s;
  char * d;
  struct snode_s * tree;


/* Line 387 of yacc.c  */
#line 932 "parse_stree.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE stree_lval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int stree_parse (void *YYPARSE_PARAM);
#else
int stree_parse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int stree_parse (struct snode_s * tree);
#else
int stree_parse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_STREE_PARSE_STREE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 960 "parse_stree.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  9
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   39

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  13
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  10
/* YYNRULES -- Number of rules.  */
#define YYNRULES  19
/* YYNRULES -- Number of states.  */
#define YYNSTATES  45

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   267

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,    13,    16,    25,    33,    38,    39,    43,
      45,    49,    50,    52,    53,    56,    57,    60,    62,    64
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      14,     0,    -1,     3,    15,     7,    15,     4,    18,    19,
      20,    10,    -1,    21,    10,    -1,     3,    15,     7,    15,
       4,    18,    19,    20,    -1,     3,    15,     4,    21,    16,
      19,    20,    -1,    21,    16,    19,    20,    -1,    -1,     5,
      17,     6,    -1,    11,    -1,    11,     7,    17,    -1,    -1,
      21,    -1,    -1,     9,    22,    -1,    -1,     8,    22,    -1,
      11,    -1,    12,    -1,    12,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   846,   846,   863,   875,   892,   910,   927,   927,   928,
     932,   941,   941,   942,   942,   943,   943,   944,   944,   945
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "OPAR", "CPAR", "OBRA", "CBRA", "COMMA",
  "HASH", "COLON", "SEMICOLON", "STRING", "NUMBER", "$accept", "input",
  "subtree", "optional_annotation", "string_list", "optional_label",
  "optional_length", "optional_theta", "label", "number", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    13,    14,    14,    15,    15,    15,    16,    16,    17,
      17,    18,    18,    19,    19,    20,    20,    21,    21,    22
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     9,     2,     8,     7,     4,     0,     3,     1,
       3,     0,     1,     0,     2,     0,     2,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    17,    18,     0,     0,     0,     0,     7,     1,
       3,     0,     0,     0,    13,     0,     0,     0,     9,     0,
       0,    15,     7,     0,    11,     0,     8,    19,    14,     0,
       6,    13,    11,    13,    12,    10,    16,    15,    13,    15,
       5,    15,     0,     4,     2
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     4,     7,    14,    19,    33,    21,    30,     8,    28
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -29
static const yytype_int8 yypact[] =
{
       1,     5,   -29,   -29,    19,    -8,     5,     7,    16,   -29,
     -29,     2,     5,    12,    18,    14,     5,    24,    22,    25,
      21,    26,    16,    31,    14,    12,   -29,   -29,   -29,    21,
     -29,    18,    14,    18,   -29,   -29,   -29,    26,    18,    26,
     -29,    26,    20,   -29,   -29
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -29,   -29,    -5,    15,    11,     6,   -28,   -19,     0,    10
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
       5,    11,    10,    37,     1,    39,    15,    17,     6,    16,
      41,    23,     2,     3,    12,    22,     2,     3,    40,     9,
      42,    13,    43,    18,    34,     2,     3,    20,    24,    25,
      44,    26,    34,    27,    29,    32,    35,    31,    38,    36
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-29)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_uint8 yycheck[] =
{
       0,     6,    10,    31,     3,    33,     4,    12,     3,     7,
      38,    16,    11,    12,     7,    15,    11,    12,    37,     0,
      39,     5,    41,    11,    24,    11,    12,     9,     4,     7,
      10,     6,    32,    12,     8,     4,    25,    22,    32,    29
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    11,    12,    14,    21,     3,    15,    21,     0,
      10,    15,     7,     5,    16,     4,     7,    15,    11,    17,
       9,    19,    21,    15,     4,     7,     6,    12,    22,     8,
      20,    16,     4,    18,    21,    17,    22,    19,    18,    19,
      20,    19,    20,    20,    10
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (tree, YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value, tree); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, struct snode_s * tree)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, tree)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    struct snode_s * tree;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
  YYUSE (tree);
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
        break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, struct snode_s * tree)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, tree)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    struct snode_s * tree;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, tree);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, struct snode_s * tree)
#else
static void
yy_reduce_print (yyvsp, yyrule, tree)
    YYSTYPE *yyvsp;
    int yyrule;
    struct snode_s * tree;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       , tree);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, tree); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, struct snode_s * tree)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, tree)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    struct snode_s * tree;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (tree);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {
      case 11: /* STRING */
/* Line 1398 of yacc.c  */
#line 827 "parse_stree.y"
        { free(((*yyvaluep).s)); };
/* Line 1398 of yacc.c  */
#line 1890 "parse_stree.c"
        break;
      case 12: /* NUMBER */
/* Line 1398 of yacc.c  */
#line 828 "parse_stree.y"
        { free(((*yyvaluep).d)); };
/* Line 1398 of yacc.c  */
#line 1897 "parse_stree.c"
        break;
      case 15: /* subtree */
/* Line 1398 of yacc.c  */
#line 826 "parse_stree.y"
        { stree_graph_destroy(((*yyvaluep).tree),NULL); };
/* Line 1398 of yacc.c  */
#line 1904 "parse_stree.c"
        break;
      case 21: /* label */
/* Line 1398 of yacc.c  */
#line 829 "parse_stree.y"
        { free(((*yyvaluep).s)); };
/* Line 1398 of yacc.c  */
#line 1911 "parse_stree.c"
        break;

      default:
        break;
    }
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (struct snode_s * tree)
#else
int
yyparse (tree)
    struct snode_s * tree;
#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
/* Line 1792 of yacc.c  */
#line 847 "parse_stree.y"
    {
  tree->left   = (yyvsp[(2) - (9)].tree);
  tree->right  = (yyvsp[(4) - (9)].tree);
  tree->label  = (yyvsp[(6) - (9)].s);
  tree->length = (yyvsp[(7) - (9)].d) ? atof((yyvsp[(7) - (9)].d)) : 0;
  tree->parent = NULL;
  tree->leaves = (yyvsp[(2) - (9)].tree)->leaves + (yyvsp[(4) - (9)].tree)->leaves;
  tree->theta  = (yyvsp[(8) - (9)].d) ? atof((yyvsp[(8) - (9)].d)) : 0;
  free((yyvsp[(7) - (9)].d));
  if ((yyvsp[(8) - (9)].d)) free((yyvsp[(8) - (9)].d));

  tree->left->parent  = tree;
  tree->right->parent = tree;

  inner_cnt++;
}
    break;

  case 3:
/* Line 1792 of yacc.c  */
#line 864 "parse_stree.y"
    {
  tree->left  = NULL;
  tree->right = NULL;
  tree->label = (yyvsp[(1) - (2)].s);
  tree->length = 0;
  tree->parent = NULL;
  tree->leaves = 1;

  tip_cnt++;
}
    break;

  case 4:
/* Line 1792 of yacc.c  */
#line 876 "parse_stree.y"
    {
  (yyval.tree) = (snode_t *)calloc(1, sizeof(snode_t));
  (yyval.tree)->left   = (yyvsp[(2) - (8)].tree);
  (yyval.tree)->right  = (yyvsp[(4) - (8)].tree);
  (yyval.tree)->label  = (yyvsp[(6) - (8)].s);
  (yyval.tree)->length = (yyvsp[(7) - (8)].d) ? atof((yyvsp[(7) - (8)].d)) : 0;
  (yyval.tree)->leaves = (yyvsp[(2) - (8)].tree)->leaves + (yyvsp[(4) - (8)].tree)->leaves;
  (yyval.tree)->theta  = (yyvsp[(8) - (8)].d) ? atof((yyvsp[(8) - (8)].d)) : 0;
  free((yyvsp[(7) - (8)].d));
  if ((yyvsp[(8) - (8)].d)) free((yyvsp[(8) - (8)].d));

  (yyval.tree)->left->parent  = (yyval.tree);
  (yyval.tree)->right->parent = (yyval.tree);

  inner_cnt++;
}
    break;

  case 5:
/* Line 1792 of yacc.c  */
#line 893 "parse_stree.y"
    {
  (yyval.tree) = (snode_t *)calloc(1, sizeof(snode_t));
  (yyval.tree)->left   = (yyvsp[(2) - (7)].tree);
  (yyval.tree)->right  = NULL;
  (yyval.tree)->label  = (yyvsp[(4) - (7)].s);
  (yyval.tree)->length = (yyvsp[(6) - (7)].d) ? atof((yyvsp[(6) - (7)].d)) : 0;
  (yyval.tree)->leaves = (yyvsp[(2) - (7)].tree)->leaves;
  (yyval.tree)->theta  = (yyvsp[(7) - (7)].d) ? atof((yyvsp[(7) - (7)].d)) : 0;
  free((yyvsp[(6) - (7)].d));
  if ((yyvsp[(7) - (7)].d)) free((yyvsp[(7) - (7)].d));

  if ((yyvsp[(5) - (7)].s)) 
    (yyval.tree)->data = (void *)((yyvsp[(5) - (7)].s));

  (yyval.tree)->left->parent = (yyval.tree);
  inner_cnt++;
}
    break;

  case 6:
/* Line 1792 of yacc.c  */
#line 911 "parse_stree.y"
    {
  (yyval.tree) = (snode_t *)calloc(1, sizeof(snode_t));
  (yyval.tree)->label  = (yyvsp[(1) - (4)].s);
  (yyval.tree)->length = (yyvsp[(3) - (4)].d) ? atof((yyvsp[(3) - (4)].d)) : 0;
  (yyval.tree)->leaves = 1;
  (yyval.tree)->theta  = (yyvsp[(4) - (4)].d) ? atof((yyvsp[(4) - (4)].d)) : 0;
  (yyval.tree)->left   = NULL;
  (yyval.tree)->right  = NULL;
  tip_cnt++;
  free((yyvsp[(3) - (4)].d));
  if ((yyvsp[(4) - (4)].d)) free((yyvsp[(4) - (4)].d));

  if ((yyvsp[(2) - (4)].s))
    (yyval.tree)->data = (void *)((yyvsp[(2) - (4)].s));
}
    break;

  case 7:
/* Line 1792 of yacc.c  */
#line 927 "parse_stree.y"
    {(yyval.s) = NULL;}
    break;

  case 8:
/* Line 1792 of yacc.c  */
#line 927 "parse_stree.y"
    {(yyval.s) = (yyvsp[(2) - (3)].s);}
    break;

  case 9:
/* Line 1792 of yacc.c  */
#line 929 "parse_stree.y"
    {
        (yyval.s) = (yyvsp[(1) - (1)].s);
}
    break;

  case 10:
/* Line 1792 of yacc.c  */
#line 933 "parse_stree.y"
    {
        (yyval.s) = (char *)xcalloc((strlen((yyvsp[(1) - (3)].s))+strlen((yyvsp[(3) - (3)].s))+2),sizeof(char));
        strcpy((yyval.s),(yyvsp[(1) - (3)].s));
        (yyval.s)[strlen((yyvsp[(1) - (3)].s))] = ',';
        strcpy((yyval.s)+strlen((yyvsp[(1) - (3)].s))+1,(yyvsp[(3) - (3)].s));
        free((yyvsp[(1) - (3)].s));
        free((yyvsp[(3) - (3)].s));
}
    break;

  case 11:
/* Line 1792 of yacc.c  */
#line 941 "parse_stree.y"
    {(yyval.s) = NULL;}
    break;

  case 12:
/* Line 1792 of yacc.c  */
#line 941 "parse_stree.y"
    {(yyval.s) = (yyvsp[(1) - (1)].s);}
    break;

  case 13:
/* Line 1792 of yacc.c  */
#line 942 "parse_stree.y"
    {(yyval.d) = NULL;}
    break;

  case 14:
/* Line 1792 of yacc.c  */
#line 942 "parse_stree.y"
    {(yyval.d) = (yyvsp[(2) - (2)].d);}
    break;

  case 15:
/* Line 1792 of yacc.c  */
#line 943 "parse_stree.y"
    {(yyval.d) = NULL;}
    break;

  case 16:
/* Line 1792 of yacc.c  */
#line 943 "parse_stree.y"
    {(yyval.d) = (yyvsp[(2) - (2)].d);}
    break;

  case 17:
/* Line 1792 of yacc.c  */
#line 944 "parse_stree.y"
    {(yyval.s)=(yyvsp[(1) - (1)].s);}
    break;

  case 18:
/* Line 1792 of yacc.c  */
#line 944 "parse_stree.y"
    {(yyval.s)=(yyvsp[(1) - (1)].d);}
    break;

  case 19:
/* Line 1792 of yacc.c  */
#line 945 "parse_stree.y"
    {(yyval.d)=(yyvsp[(1) - (1)].d);}
    break;


/* Line 1792 of yacc.c  */
#line 2387 "parse_stree.c"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (tree, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (tree, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval, tree);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp, tree);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (tree, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, tree);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, tree);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2055 of yacc.c  */
#line 947 "parse_stree.y"


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
    opt_msci = 1;

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

  if (opt_msci)
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

  if (opt_msci && stree->root->htau == 0)
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
