/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

/* tokens returned by lexical analyzer */
#define TOKEN_NONE              0
#define TOKEN_ATTR              1
#define TOKEN_OPAR              2
#define TOKEN_CPAR              3
#define TOKEN_COMMA             4
#define TOKEN_COLON             5
#define TOKEN_SEMICOLON         6
#define TOKEN_HASH              7
#define TOKEN_STRING            8

const char * token_names[9] = {
  "Illegal", "TOKEN_ATTR", "TOKEN_OPAR", "TOKEN_CPAR", "TOKEN_COMMA",
  "TOKEN_COLON", "TOKEN_SEMICOLON", "TOKEN_HASH", "TOKEN_STRING"
};

static void node_destroy(node_t * root, void (*cb_data_destroy)(void *));

const unsigned int attrib_map[256] = 
 {
/* 0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F  */

   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 0 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 1 */    
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 2 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 3 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 4 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  /* 5 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 6 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  /* 7 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 8 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 9 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* A */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* B */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* C */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* D */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* E */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* F */
 };

const unsigned int string_map[256] = 
 {
/* 0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F  */

   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 0 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 1 */    
   0,  1,  0,  1,  1,  1,  1,  0,  0,  0,  1,  1,  0,  1,  1,  1,  /* 2 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1,  1,  /* 3 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 4 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  /* 5 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 6 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  /* 7 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 8 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 9 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* A */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* B */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* C */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* D */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* E */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* F */
 };

/* Syntax parsing table. Entry (i,j) indicates whether token j can follow
   token i (1) or not (0) */
const unsigned int syntax_table[9][9] =
 {
   /* NONE, ATTR, OPAR, CPAR, COMMA, COLON, SEMICOLON, HASH, STRING */
   {  0, 0, 1, 0, 0, 0, 0, 0, 1 },   /* NONE */
   {  0, 0, 0, 1, 1, 1, 1, 1, 0 },   /* ATTR */
   {  0, 0, 1, 0, 0, 0, 0, 0, 1 },   /* OPAR */
   {  0, 1, 0, 1, 1, 1, 1, 1, 1 },   /* CPAR */
   {  0, 0, 1, 0, 0, 0, 0, 0, 1 },   /* COMMA */
   {  0, 0, 0, 0, 0, 0, 0, 0, 1 },   /* COLON */
   {  0, 0, 0, 0, 0, 0, 0, 0, 0 },   /* SEMICOLON */
   {  0, 0, 0, 0, 0, 0, 0, 0, 1 },   /* HASH */
   {  0, 1, 0, 1, 1, 1, 1, 1, 0 },   /* STRING */
 };

typedef struct ltoken_s
{
  char * data;
  int type;
} ltoken_t;


static long parse_attr(char * s, list_t * token_list)
{
  long i = 0;
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));

  assert(*s == '[');
  ++s;

  /* match string */
  while (attrib_map[(int)(s[i])]) ++i;

  if (s[i] != ']')
    fatal("[ERROR] Cannot match closing bracket in attribute: %s", s-1);

  token->data = xstrndup(s-1,i+2);
  token->type = TOKEN_ATTR;

  list_append(token_list,(void *)token);

  /* return length of token */
  return i+2;
}

static void token_clear(void * tokenptr)
{
  ltoken_t * token = (ltoken_t *)tokenptr;
  if (token)
  {
    if (token->data)
      free(token->data);
    free(token);
  }
}

long parse_opar(char * s, list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup("(");
  token->type = TOKEN_OPAR;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_cpar(char * s, list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(")");
  token->type = TOKEN_CPAR;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_colon(char * s, list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(":");
  token->type = TOKEN_COLON;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_semicolon(char * s, list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(";");
  token->type = TOKEN_SEMICOLON;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_comma(char * s, list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(",");
  token->type = TOKEN_COMMA;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_theta(char * s, list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup("#");
  token->type = TOKEN_HASH;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_string(char * s, list_t * token_list)
{
  long i = 0;
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  
  /* match string */
  while (string_map[(int)(s[i])]) ++i;

  token->data = xstrndup(s,i);
  token->type = TOKEN_STRING;

  list_append(token_list,(void *)token);

  return i;
}

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
  if (root->attrib)
    free(root->attrib);
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

    if (opt_datefile && !opt_simulate) {
	if (node->node_index < tree->tip_count + opt_seqAncestral) {
		for (j = 0; j < opt_locus_count; j++) {
	
			if (node->epoch_count[j]) {
				free(node->date_count[j]);
				free(node->tip_date[j]);
			}	
		} 
		free(node->epoch_count);
		free(node->tip_date);
		free(node->date_count);
	
	}

    }

    if (opt_migration)
    {
      if (node->mig_target)
      {
        for (j = 0; j < tree->locus_count; ++j)
        {
          /* 
            IMPORTANT: Do not delete the dlist_item_t's with the below call:

            dlist_clear(tree->nodes[i]->mig_target[j],NULL); 

            as these elements are deleted in
            
              gtree_destroy(..)  -> miginfo_destroy(..)
          */

          dlist_destroy(tree->nodes[i]->mig_target[j]);
        }
        free(node->mig_target);
      }

      if (node->mig_source)
      {
        for (j = 0; j < tree->locus_count; ++j)
        {
          /*  IMPORTANT: Same applies here, see above

          dlist_clear(tree->nodes[i]->mig_source[j],NULL);

          */

          dlist_destroy(tree->nodes[i]->mig_source[j]);
        }
        free(node->mig_source);
      }
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

    if (opt_migration && node->migevent_count)
      free(node->migevent_count);

    if (opt_migration && node->migbuffer)
      free(node->migbuffer);

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

  if (tree->td)
    free(tree->td);

  if (tree->mi_tbuffer)
  {
    for (i = 0; i < opt_threads; ++i)
    {
      /* all dlist_items should be unlinked */
      if (tree->mi_tbuffer[i]) // && tree->mi_tbuffer[i]->me)
      {
        /* dlist_item_ts must already be unlinked, so we set count to 0 to avoid
           unlinking them again in miginfo_destroy. Locus number is set to -1 as
           it's irrelevant, the function will not enter the locus loop */
        tree->mi_tbuffer[i]->count = 0;
        miginfo_destroy(tree->mi_tbuffer[i], -1, MI_DLI_FREE);
      }
    }
    free(tree->mi_tbuffer);
  }
  if (tree->migcount_sum)
    free(tree->migcount_sum);

  free(tree->u_constraint);
  free(tree->l_constraint);
  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

int opt_precision = 6;
static char * ntree_export_newick_recursive(node_t * root, long print_bl)
{
  int i;
  char * newick;
  char * x;

  if (!root) return NULL;

  if (!root->children_count)
  {
    if (print_bl)
      xasprintf(&newick,
                "%s%s:%.*f",
                root->label,
                root->attr ? root->attr : "",
                opt_precision,
                root->length);
    else
      xasprintf(&newick,
                "%s%s",
                root->label,
                root->attr ? root->attr : "");
  }
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0],print_bl);
    xasprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i],print_bl);
      xasprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    if (print_bl)
      xasprintf(&x,
                "%s)%s%s:%.*f",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "", 
                opt_precision,
                root->length);
    else
      xasprintf(&x,
                "%s)%s%s",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "");
    free(newick);
    newick = x;
  }

  return newick;
}

char * ntree_export_newick(ntree_t * tree, long print_bl)
{
  int i;
  char * newick;
  char * x;

  node_t * root = tree->root;

  if (!root) return NULL;

  if (!root->children_count)
  {
    if (print_bl)
      xasprintf(&newick,
                "%s%s:%.*f",
                root->label,
                root->attr ? root->attr : "",
                opt_precision,
                root->length);
    else
      xasprintf(&newick,
                "%s%s",
                root->label,
                root->attr ? root->attr : "");
  }
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0],print_bl);
    xasprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i],print_bl);
      xasprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    if (print_bl)
      xasprintf(&x,
                "%s)%s%s:%.*f;",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "",
                opt_precision,
                root->length);
    else
      xasprintf(&x,
                "%s)%s%s",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "");
    free(newick);
    newick = x;
  }

  return newick;
}

list_t * parse_tree(char * s)
{
  int rc = 1;
  long count;
  long opar_count = 0;
  long cpar_count = 0;
  list_t * token_list = NULL;

  /* skip all white-space */
  size_t ws = strspn(s, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!s[ws] || s[ws] == '*' || s[ws] == '#')
  {
    rc = 0;
    goto l_unwind;
  }

  char * start = s+ws;

  token_list = (list_t *)xcalloc(1,sizeof(list_t));

  /* TODO: Find end by going backwards and set a zero char */

  while (*start)
  {
    size_t ws = strspn(start, " \t\r\n");
    start += ws;

    if (*start == '(')
    {
      if (opt_debug_parser)
        printf("Parsing OPAR...\n");
      count = parse_opar(start,token_list);
      start += count;
    }
    else if (*start == ')')
    {
      if (opt_debug_parser)
        printf("Parsing CPAR...\n");
      count = parse_cpar(start,token_list);
      start += count;
    }
    else if (*start == '[')
    {
      if (opt_debug_parser)
        printf("Parsing ATTR...\n");
      count = parse_attr(start,token_list);
      start += count;
    }
    else if (*start == ':')
    {
      if (opt_debug_parser)
        printf("Parsing COLON...\n");
      count = parse_colon(start,token_list);
      start += count;
    }
    else if (*start == ';')
    {
      if (opt_debug_parser)
        printf("Parsing SEMICOLON...\n");
      count = parse_semicolon(start,token_list);
      start += count;
    }
    else if (*start == ',')
    {
      if (opt_debug_parser)
        printf("Parsing COMMA...\n");
      count = parse_comma(start,token_list);
      start += count;
    }
    else if (*start == '#')
    {
      if (opt_debug_parser)
        printf("Parsing HASH...\n");
      count = parse_theta(start,token_list);
      start += count;
    }
    else if (string_map[(int)(*start)])
    {
      if (opt_debug_parser)
        printf("Parsing STRING...\n");
      /* TODO : eliminate starting with theta  */
      count = parse_string(start,token_list);
      start += count;
    }
    else
    {
      fatal("Illegal character");
    }

  }
  if (opt_debug_parser)
  {
    printf("Finished\nListing tokens:\n");

    list_item_t * li = token_list->head;
    while (li)
    {
      ltoken_t * token = (ltoken_t *)(li->data);

      printf("%s : %s\n", token_names[token->type], token->data);

      li = li->next;
    }
  }

  
  /* trivial checks */

  /* 1. No tokens means error */
  if (token_list->count == 0)
    fatal("No valid constructs found");


  /* 2. First token must be either OPAR or STRING */
  list_item_t * li = token_list->head;
  ltoken_t * token = (ltoken_t *)(li->data);
  if (token->type != TOKEN_OPAR && token->type != TOKEN_STRING)
    fatal("Newick string must start with '(' or a label");
    
  /* 3. If the first token is not an OPAR then no OPAR and CPAR must exist in
     the list */
  if (token->type != TOKEN_OPAR)
  {
    li = li->next;
    while (li)
    {
      token = (ltoken_t *)(li->data);
      if (token->type == TOKEN_OPAR || token->type == TOKEN_CPAR)
        break;
      li = li->next;
    }
    if (li)
      fatal("Illegal newick format specification - found opening/closing "
            "parenthesis but newick does not start with opening parenthesis");
  }

  /* 4. If first token is OPAR, then check that there is no OPAR when the
     number of open OPARs is 0 */
  li = token_list->head;
  token = (ltoken_t *)(li->data);
  if (token->type == TOKEN_OPAR)
  {
    opar_count++;
    long openpars = 1;
    li = li->next;
    while (li)
    {
      token = (ltoken_t *)(li->data);
      if (token->type == TOKEN_OPAR)
      {
        if (openpars == 0)
        {
          snprintf(bpp_errmsg, 200, "Invalid token type while parsing n-ary tree");
          rc = 0;
          goto l_unwind;
          #if 0
          fatal("Cannot open new parenthesis");
          #endif
        }
        ++openpars;
        ++opar_count;
      }
      else if (token->type == TOKEN_CPAR)
      {
        --openpars;
        ++cpar_count;
      }
      if (openpars < 0)
        fatal("Invalid newick format");

      li = li->next;
    }
  }

  if (opar_count != cpar_count)
  {
    snprintf(bpp_errmsg,
             200,
             "Mismatching number of opening (%ld) and closing (%ld) parentheses)",
             opar_count, cpar_count);
    rc = 0;
    goto l_unwind;
    #if 0
    fatal("Mismatching number of opening (%ld) and closing (%ld) parentheses)",
          opar_count, cpar_count);
    #endif
  }

l_unwind:
  if (!rc)
  {
    if (token_list)
    {
      list_clear(token_list,token_clear);
      free(token_list);
      token_list = NULL;
    }
  }
  return token_list;
}

static void parse_annotation(snode_t * snode, const char * annotation)
{
  double val;

  char * p = xstrdup(annotation);
  char * s = p;

  /* Please note annotation is in the format:

     [&attrib1=A1,&attrib2=A2,...,&attribN=AN]

     The ampersand (&) is optional.
     We thus, remove the opening/closing square brackets.
  */

  assert(*s == '[');
  *s = 0; s++;
  size_t tlen = strlen(s);
  assert(s[tlen-1] == ']');
  s[tlen-1] = 0;

  /* remove all spaces and tabs */
  char * dst = s;
  char * src = s;
  while (*src)
  {
    if (*src != ' ' && *src != '\t')
    {
      if (src != dst)
        *dst = *src;

      ++dst;
    }
    ++src;
  }
  *dst = 0;


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

    if (!strcasecmp(s,"defphi") || !strcasecmp(s,"&defphi"))
    {
      snode->has_phi = 1;
      s += tokenlen+comma;

      continue;
    }


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
                                      snode_t * yparent,
                                      snode_t * ychild)
{
  /* yparent is the parent of xtip, and ychild is the child of xinner

     W.l.o.g. the idea for annotations is the following:
       Annotation of xtip relates to the edge Y->X
       Annotation of xinner relates to R->X (R is parent of X)
       Annotation of yparent relates to P->Y (P is parent of Y)
       Annotation of ychild relates to X->Y

     Also, nodes xtip and ychild are to be deleted later on, therefore
     we do not annotate them

  */

  /* better safe than sorry */
  assert(!xtip->left && !xtip->right);          /* xtip is tip */
  assert(xinner->left && xinner->right);        /* xinner is inner */
  assert(yparent->left || yparent->right);      /* yparent is inner */
  assert(ychild->left || ychild->right);        /* ychild is inner */
  
  if (xtip->attrib)
  {
    parse_annotation(xtip, xtip->attrib);
    free(xtip->attrib);
    xtip->attrib = NULL;
  }
  if (xinner->attrib)
  {
    parse_annotation(xinner, xinner->attrib);
    free(xinner->attrib);
    xinner->attrib = NULL;
  }
  if (yparent->attrib)
  {
    parse_annotation(yparent, yparent->attrib);
    free(yparent->attrib);
    yparent->attrib = NULL;
  }
  if (ychild->attrib)
  {
    parse_annotation(ychild, ychild->attrib);
    free(ychild->attrib);
    ychild->attrib = NULL;
  }

  if (xtip->htau == 1 || xinner->htau || yparent->htau || ychild->htau)
    fatal("Bidirectional introgression between (%s) and (%s) requires exactly "
          "one tau parameter shared by both nodes. Please remove *all* 'tau' "
          "annotations for these two nodes.",
          xtip->label, yparent->label);

  if (xtip->hphi && xinner->hphi)
  {
    if (xtip->hphi + xinner->hphi != 1)
      fatal("Phi parameter annotations for bidirectional introgression event "
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
  {
    if (opt_simulate)
      fatal("Missing phi parameter annotation for bidirectional introgression "
            "event on node (%s)", xtip->label);
  }

  if (yparent->hphi && ychild->hphi)
  {
    if (yparent->hphi + ychild->hphi != 1)
      fatal("Phi parameter annotations for bidirectional introgression event "
            "on node (%s) do not sum to 1",  yparent->label);
  }
  else if (yparent->hphi)
  {
    yparent->hybrid->hphi = 1 - yparent->hphi;
  }
  else if (ychild->hphi)
  {
    yparent->hphi = 1 - ychild->hphi;
    yparent->hybrid->hphi = ychild->hphi;
  }
  else
  {
    if (opt_simulate)
      fatal("Missing phi parameter annotation for bidirectional introgression "
            "event on node (%s)", yparent->label);
  }
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

  if (hinner->attrib)
  {
    parse_annotation(hinner, hinner->attrib);
    free(hinner->attrib);
    hinner->attrib = NULL;
  }

  if (htip->attrib)
  {
    parse_annotation(hinner->hybrid, htip->attrib);
    free(htip->attrib);
    htip->attrib = NULL;
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
      fatal("Phi parameter annotations for hybridization event (%s) do not "
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
    if (opt_simulate)
      fatal("Missing phi parameter annotation for hybridization event (%s)",
            hinner->label);
  }
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

  unsigned int tip_cnt = stree->tip_count;
  unsigned int inner_cnt = stree->inner_count;
  unsigned int hybrid_cnt = stree->hybrid_count;

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
  snode_t * xtip = NULL;
  snode_t * ychild = NULL;

  /* check for bidirectional introgression 

  two possible notations:

  Diagram 1:  ((A,(B)Y)X,(X)Y)R;
  Diagram 2:  ((A,Y)X,(X,B)Y)R;


         *  R             *  R             
        / \              / \         Method:
     X /   \     or     /   \          1. Find a tip (xtip) X for which an inner
      *     \          /     \            node with same label exists (xinner)
     / \  Y  \ Y    X /       \  Y     2. Get the parent of tip node X (in this
    /   *     *      *         *          case node Y and call it yparent)
   /    |     |     / \       / \      3. Check that inner node X has a child
  /     |     |    /   \     /   \        with label Y which is either unary or 
  A     B     X    A    Y   X     B       tip (name it ychild)
                                       4. Create hybrid nodes for X and Y
                                       5. Link the corresponding nodes to obtain
                                          the diagram below
   Results in:

         R
        /\
       /  \
   X  *<-->*  Y
     /      \
    /        \
   A          B

  In the case of diagram two (four tips) we transform it to the first case by
  moving subtree B to become a child of tip node Y (i.e. tip Y becomes unary
  inner node). This change of 'status' of node Y leads to a decrease of one
  tip and increase of one inner node.


  */

  unsigned int tip_cnt = stree->tip_count;
  unsigned int inner_cnt = stree->inner_count;
  unsigned int hybrid_cnt = stree->hybrid_count;

  /* check for bidirectional introgression */
  for (i = 0; i < stree->tip_count; ++i)
  {
    /* if no duplicate label for current tip node, move to next tip */
    if (!stree->nodes[i] || stree->nodes[i]->left || dups[i] == i) continue;

    ychild = xtip = NULL;

    /* if the duplicate of tip is a tip, fail as we checked this in
       create_duplicate_indices */
    assert(dups[i] >= stree->tip_count);

    /* root cannot be duplicated */
    assert(stree->nodes[i]->parent);

    /* get pointers to the two inner nodes which are candidates for a
       bidirectional introgression */
    snode_t * xinner  = stree->nodes[dups[i]];        /* inner node x */
    snode_t * yparent = stree->nodes[i]->parent;      /* parent of tip node x */
    if (!yparent->label) continue;
    xtip = stree->nodes[i];                           /* tip node x */

    assert(xinner->left || xinner->right);      /* confirm x is an inner node */
    assert(xinner->label);


    if (xinner->left && xinner->left->label &&
        !strcmp(xinner->left->label,yparent->label))
      ychild = xinner->left;
    else if (xinner->right && xinner->right->label &&
             !strcmp(xinner->right->label,yparent->label))
      ychild = xinner->right;

    if (ychild)
    {
      if (!(ychild->left && !ychild->right))
      {
        /* if ychild is not an unary node it must be a tip */
        assert(!ychild->left && !ychild->right);
        
        /* also y must have two children */
        assert(yparent->left && yparent->right);

        /* swap subtrees */
        if (yparent->left == xtip)
        {
          ychild->left = yparent->right;
          yparent->right->parent = ychild;
          yparent->right = NULL;
        }
        else
        {
          assert(yparent->right == xtip);
          ychild->left = yparent->left;
          yparent->left->parent = ychild;
          yparent->left = yparent->right;
          yparent->right = NULL;
        }
        tip_cnt--;
        inner_cnt++;
      }
      /* ychild must be an unary node */
      assert(ychild->left && !ychild->right);

      hybrid_cnt += 2;
      tip_cnt -= 1;
      inner_cnt--;

      /* we create an additional node for each introgression event */
      xinner->hybrid = (snode_t *)xcalloc(1,sizeof(snode_t));
      yparent->hybrid = (snode_t *)xcalloc(1,sizeof(snode_t));
      assert(xinner->label);
      assert(yparent->label);
      xinner->hybrid->label = xstrdup(xinner->label);
      yparent->hybrid->label = xstrdup(yparent->label);

      /* link hybridization nodes with their corresponding nodes */
      xinner->hybrid->hybrid = xinner;
      yparent->hybrid->hybrid = yparent;

      /* define parents of (mirror) hybridization nodes */
      xinner->hybrid->parent = yparent;
      yparent->hybrid->parent = xinner;

      /* define one child (as left child node) for hybridization nodes */
      xinner->hybrid->left   = NULL;
      yparent->hybrid->left  = NULL;
      xinner->hybrid->right  = NULL;
      yparent->hybrid->right = NULL;

      /* if simulation, copy branch length and theta */
      if (opt_simulate)
      {
        xinner->hybrid->tau   = xtip->tau;
        xinner->hybrid->theta = xtip->theta;
      }

      /* annotate */
      annotate_bd_introgression(xtip, xinner, yparent, ychild);

      /* now destroy redundant nodes/subtrees (i.e. xtip and ychild) */
      assert(yparent == xtip->parent);
      assert(yparent->left == xtip && yparent->right == NULL);
      yparent->left = NULL;
      xtip->parent = NULL;
      stree->nodes[xtip->node_index] = NULL;
      stree_graph_destroy(xtip,NULL);

      /* link y with child of ychild */
      yparent->left = ychild->left;
      yparent->left->parent = yparent;

      /* now destroy ychild */
      if (ychild->parent->left == ychild)
      {
        ychild->parent->left = ychild->parent->right;
        ychild->parent->right = NULL;
      }
      else
      {
        assert(ychild == ychild->parent->right);
        ychild->parent->right = NULL;
      }
      ychild->parent = NULL;
      stree->nodes[ychild->node_index] = NULL;
      assert(ychild->left && ychild->right == NULL);
      ychild->left = NULL;

      /* if simulation, copy branch length and theta */
      if (opt_simulate)
      {
        yparent->hybrid->tau   = ychild->tau;
        yparent->hybrid->theta = ychild->theta;
      }

      stree_graph_destroy(ychild,NULL);

      xinner->right = yparent->hybrid;
      yparent->right = xinner->hybrid;

      yparent->leaves = yparent->left->leaves;
      xinner->leaves = xinner->left->leaves;
      xinner->hybrid->leaves = 0;
      yparent->hybrid->leaves = 0;
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
     the species tree (now network):
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

      /* NOTE: Consecutive unary nodes are still possible, see issue #146 on
         GitHub. I disabled this check for now */
      #if 0
      if ((stree->nodes[i]->parent->left && !stree->nodes[i]->parent->right) ||
          (!stree->nodes[i]->parent->left && stree->nodes[i]->parent->right))
        fatal("Tree contains consecutive unary nodes");
      #endif
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
      fatal("Node with label %s appears %ld times (allowed twice)",
            stree->nodes[i]->label, freq);
    else if (freq == 2)
    {
      /* if both are tips */
      if (i < stree->tip_count && dups[i] < stree->tip_count)
        fatal("Tip node with label %s appears more than once in the tree",
              stree->nodes[i]->label);
    }
  }
  return dups;
}

static void set_phi_values(stree_t * stree)
{
  unsigned int i;
  snode_t * snode;
  long thread_index = 0;

  for (i=0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];

    if (snode->hybrid)
    {
      if (fabs(snode->hphi + snode->hybrid->hphi - 1) > 1e-10)
      {
        if (node_is_bidirection(snode))
        {
          /* bidirection (model D) */

          /* we set phi for the vertical branch to U|(0.7,0.9). The snode is
             not the mirror node, hence it is the vertical branch */
          double a = 0.7; double b = 0.9;
          double r = (b-a)*legacy_rndu(thread_index) + a;
          
          snode->hphi = r;
          snode->hybrid->hphi = 1-r;
        }
        else
        {
          /* hybridization */

          /* for models A and C draw the value of phi from U(0,1). For model B
             set phi for the vertical branch to U(0.7,0.9) */
          if ((!snode->htau && !snode->hybrid->htau) ||
              (snode->htau && snode->hybrid->htau))
          {
            /* model A or C */
            snode->hphi = legacy_rndu(thread_index);
            snode->hybrid->hphi = 1 - snode->hphi;
          }
          else
          {
            /* model B */
            double a = 0.7; double b = 0.9;
            double r = (b-a)*legacy_rndu(thread_index) + a;
            if (snode->htau)
            {
              snode->hphi = r;
              snode->hybrid->hphi = 1 - r;
            }
            else
            {
              snode->hybrid->hphi = r;
              snode->hphi = 1 - r;
            }
          }
        }
      }
    }
  }
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

  if (opt_cfile && stree->hybrid_count)
    set_phi_values(stree);

  validate_stree(stree);
}

static void fill_node_lists_recursive(node_t * node,
                                      node_t ** tiplist,
                                      node_t ** innerlist,
                                      int * tipcount,
                                      int * innercount)
{
  int i;

  if (!node) return;

  /* tip node */
  if (!node->children_count)
  {
    tiplist[(*tipcount)++] = node;

    return;
  }

  /* inner node */
  for (i = 0; i < node->children_count; ++i)
  {
    fill_node_lists_recursive(node->children[i],
                              tiplist,
                              innerlist,
                              tipcount,
                              innercount);
  }
  innerlist[(*innercount)++] = node;
}

/* fill the the ntree_t->leaves and ntree_t->inner list of nodes by recursively
   traversing the tree graph starting with the root node */
static void fill_node_lists(node_t * node,
                            node_t ** tiplist,
                            node_t ** innerlist)
{
  int i;
  int tipcount = 0;
  int innercount = 0;

  if (!node) return;

  /* tip node */
  if (!node->children_count)
  {
    tiplist[tipcount++] = node;
    return;
  }

  /* inner nodes */
  for (i = 0; i < node->children_count; ++i)
  {
    fill_node_lists_recursive(node->children[i],
                              tiplist,
                              innerlist,
                              &tipcount,
                              &innercount);
  }
  innerlist[innercount++] = node;
}

ntree_t * ntree_wraptree(node_t * root, int tip_count, int inner_count)
{
  long i;

  ntree_t * tree = (ntree_t *)xcalloc(1,sizeof(ntree_t));

  tree->root = root;
  tree->tip_count = tip_count;
  tree->inner_count = inner_count;

  tree->leaves = (node_t **)xmalloc(tree->tip_count * sizeof(node_t *));
  tree->inner  = (node_t **)xmalloc(tree->inner_count * sizeof(node_t *));

  /* fill node lists in postorder traversal */
  fill_node_lists(tree->root, tree->leaves, tree->inner);

  for (i = 0; i < tree->tip_count; ++i)
    tree->leaves[i]->node_index = i;

  for (i = 0; i < tree->inner_count; ++i)
    tree->inner[i]->node_index = i;

  return tree;
}

static ntree_t * syntax_parse(list_t * token_list)
{
  ltoken_t * token;
  list_item_t * li;
  long prev_token_type = TOKEN_NONE;

  node_t * node = NULL;
  node_t * root = NULL;
  node_t * tmp  = NULL;
  node_t ** children;
  double value;

  int tip_count = 0;
  int inner_count = 0;

  li = token_list->head;

  while (li)
  {
    token = (ltoken_t *)(li->data);

    if (!syntax_table[prev_token_type][token->type])
    {
      #if 0
      printf("%ld %d\n", prev_token_type, token->type);
      #endif
      #if 0
      fatal("Invalid token type while parsing n-ary tree");
      #endif
      snprintf(bpp_errmsg, 200, "Invalid token type while parsing n-ary tree");
      if (root)
        node_destroy(root,NULL);
      return NULL;
    }

    switch (token->type)
    {
      case TOKEN_NONE:
        fatal("Unexpected TOKEN_NONE");
        break;

      case TOKEN_ATTR:
        node->attr = xstrdup(token->data);
        break;

      case TOKEN_OPAR:
        tmp = (node_t *)xcalloc(1,sizeof(node_t));
        tmp->parent = node;
        node = tmp;
        if (!root)
          root = node;
        inner_count++;
        break;

      case TOKEN_CPAR:
        assert(node);
        node = node->parent;

        assert(node);
        if (node->parent)
        {
          node->parent->children_count++;

          if (node->parent->children)
          {
            children = (node_t**)xmalloc((size_t)(node->parent->children_count)*
                                          sizeof(node_t *));
            memcpy(children,node->parent->children,
                   (node->parent->children_count-1)*sizeof(node_t *));
            children[node->parent->children_count-1] = node;
            free(node->parent->children);
            node->parent->children = children;
          }
          else
          {
            assert(node->parent->children_count == 1);
            node->parent->children = (node_t **)xmalloc(sizeof(node_t *));
            node->parent->children[0] = node;
          }
        }
        break;

      case TOKEN_COMMA:
        node = node->parent;
        break;

      case TOKEN_COLON:
        break;

      case TOKEN_SEMICOLON:
        break;

      case TOKEN_HASH:
        break;

      case TOKEN_STRING:
        switch (prev_token_type)
        {
          case TOKEN_NONE:
            assert(node == NULL);
            node = (node_t *)xcalloc(1,sizeof(node_t));
            node->parent = NULL;
            node->label = xstrdup(token->data);
            root = node;
            tip_count++;
            break;

          case TOKEN_OPAR:
            tmp = (node_t *)xcalloc(1,sizeof(node_t));
            tmp->parent = node;
            tmp->label = xstrdup(token->data);
            node = tmp;

            assert(node->parent);
            assert(node->parent->children_count == 0 && !node->parent->children);
            node->parent->children_count = 1;
            node->parent->children = (node_t **)xmalloc(sizeof(node_t *));
            node->parent->children[0] = node;
            tip_count++;
            break;

          case TOKEN_CPAR:
            node->label = xstrdup(token->data);
            break;

          case TOKEN_COMMA:
            tmp = (node_t *)xcalloc(1,sizeof(node_t));
            tmp->parent = node;
            tmp->label = xstrdup(token->data);
            node = tmp;

            node_t * parent = node->parent;
            assert(parent);
            assert(parent->children_count > 0 && parent->children);
            parent->children_count++;

            children = (node_t **)xmalloc((size_t)(parent->children_count) *
                                          sizeof(node_t *));
            memcpy(children,
                   parent->children,
                   (parent->children_count-1)*sizeof(node_t *));
            children[parent->children_count-1] = node;
            free(parent->children);
            parent->children = children;
            tip_count++;
            break;
            
          case TOKEN_COLON:
            if (!get_double(token->data,&value))
              fatal("ERROR: Expected floating point number (branch length) but "
                    "got: %s", token->data);
            if (value < 0)
              fatal("ERROR: Found negative branch length (%f)", value);

            node->length = value;
            break;

          case TOKEN_HASH:
            if (!get_double(token->data,&value))
              fatal("ERROR: Expected floating point number (theta) but got: %s",
                    token->data);
            if (value < 0)
              fatal("ERROR: Found negative theta (%f)", value);

            node->theta = value;
            break;

          default:
            fatal("Internal error when parsing syntax table");

        }
        break;

      default:
        fatal("Unknown token type");
    }

    prev_token_type = token->type;
    li = li->next;
  }

  return ntree_wraptree(root,tip_count,inner_count);
}

static void node_destroy(node_t * root, void (*cb_data_destroy)(void *))
{
  int i;

  if (!root) return;

  if (root->children)
  {
    for (i = 0; i < root->children_count; ++i)
      node_destroy(root->children[i],cb_data_destroy);

    free(root->children);
  }
  
  if (root->data && cb_data_destroy)
    cb_data_destroy(root->data);

  if (root->attr)
    free(root->attr);

  free(root->label);
  free(root);
}


void ntree_destroy(ntree_t * tree, void (*cb_data_destroy)(void *))
{
  if (!tree) return;

  node_destroy(tree->root,cb_data_destroy);

  if (tree->leaves)
    free(tree->leaves);

  if (tree->inner)
    free(tree->inner);

  free(tree);
}

int ntree_check_rbinary(ntree_t * tree)
{
  int i;
  /* checks whether tree is binary rooted */

  assert(!tree->root->parent);

  /* check all inner nodes that they have out-degree 2 */
  for (i = 0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i]->children_count != 2)
      return 0;
  }

  return 1;
}

int ntree_check_ubinary(ntree_t * tree)
{
  int i;
  /* checks whether tree is unrooted binary */

  assert(!tree->root->parent);

  for (i = 0; i < tree->inner_count; ++i)
  {
    if (!tree->inner[i]->parent)
    {
      /* root case */
      if (tree->inner[i]->children_count != 3) return 0;
    }
    else
    {
      /* other inner case */

      if (tree->inner[i]->children_count != 2) return 0;
    }
  }

  return 1;
}

int ntree_check_nary(ntree_t * tree)
{
  int i;
  /* checks whether tree is n-ary */

  assert(!tree->root->parent);

  for (i = 0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i]->children_count > 2) return 1;
  }

  return 0;
}

static snode_t * snode_from_node_recursive(node_t * node, snode_t * parent)
{
  snode_t * snode = (snode_t *)xcalloc(1,sizeof(snode_t));

  assert(node->children_count <= 2);

  if (node->children_count)
  {
    snode->left  = snode_from_node_recursive(node->children[0],snode);
    if (node->children_count > 1)
      snode->right = snode_from_node_recursive(node->children[1],snode);
    else
      snode->right = NULL;

  }
  else
  {
    snode->left = snode->right = NULL;
  }

  snode->parent = parent;
  snode->length = node->length;
  snode->theta  = node->theta;
  snode->tau    = node->tau;
  if (node->label)
    snode->label = xstrdup(node->label);

  if (node->attr)
    snode->attrib = xstrdup(node->attr);

  return snode;
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

  int hashSize = stree->tip_count;

  if (commas_count+1 != stree->tip_count) {
          if (!(opt_datefile && commas_count + 1 == stree->tip_count + stree->inner_count)) {
                fatal("Labels (%d) specified in --reorder do not match species tree (%d)", commas_count, stree->tip_count);
          } else {
                hashSize = stree->tip_count + stree->inner_count;
                opt_seqAncestral = stree->inner_count;

          }
  }


  hashtable_t * ht = hashtable_create(hashSize);

  pairlist = (pair_t **)xmalloc(hashSize * sizeof(pair_t *));

  for (i = 0; i < hashSize; ++i)
  {
    if (!stree->nodes[i]->label) {
	    fatal("Internal nodes are not all labeled. They must be labeled with sampling from ancestral populations.");
    }
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

static void fill_nodes_recursive(snode_t * snode,
                                 snode_t ** array,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!snode->left && !snode->right)
  {
    array[*tip_index] = snode;
    *tip_index = *tip_index + 1;
    return;
  }

  array[*inner_index] = snode;
  *inner_index = *inner_index + 1;

  if (snode->left)
    fill_nodes_recursive(snode->left,  array, tip_index, inner_index);
  if (snode->right)
    fill_nodes_recursive(snode->right, array, tip_index, inner_index);
}

stree_t * stree_from_ntree(ntree_t * ntree)
{
  /* assumes rooted binary ntree */
  unsigned int i;
  unsigned int tip_index = 0;
  unsigned int inner_index = ntree->tip_count;

  stree_t * stree = (stree_t *)xcalloc(1,sizeof(stree_t));

  /* trivial validation of tree */
  if (!ntree_check_rbinary(ntree))
  {
    if (ntree_check_ubinary(ntree))
    {
      fatal("BPP requires a rooted binary tree as input."
             "Given tree is unrooted binary");
    }
    else
    {
      /* check if tree is n-ary (n>2) */
      if (ntree_check_nary(ntree))
        fatal("BPP requires a rooted binary tree as input."
              "Given tree is multifurcating");
    }
  }

  stree->root =  snode_from_node_recursive(ntree->root,NULL);

  stree->tip_count = ntree->tip_count;
  stree->inner_count = ntree->inner_count;

  stree->edge_count = stree->tip_count+stree->inner_count-1;
  stree->hybrid_count = 0;
  stree->pptable = NULL;

  stree->nodes = (snode_t **)xmalloc((size_t)(stree->tip_count+stree->inner_count) *
                                     sizeof(snode_t *));

  fill_nodes_recursive(stree->root, stree->nodes, &tip_index, &inner_index);

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    stree->nodes[i]->node_index = i;

  resolve_network(stree);
  if (stree->hybrid_count)
    opt_msci = 1;

  /* reorder tip nodes if specified */
  if (opt_reorder)
    reorder(stree);

  for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
    stree->nodes[i]->node_index = i;

  /* apply diploid information */
  if (opt_diploid)
  {
    if (opt_diploid_size != stree->tip_count + opt_seqAncestral)
      fatal("Number of 'phase' assignments mismatch number of species");

    for (i = 0; i < stree->tip_count; ++i)
      stree->nodes[i]->diploid = opt_diploid[i];
  }

  if (opt_msci)
  {
    /* for all nodes that are not hybridization events set prop_tau to 1 */
    for (i = 0; i < stree->inner_count; ++i)
    {
      snode_t * snode = stree->nodes[stree->tip_count+i];

      if (!snode->hybrid) snode->prop_tau = 1;
    }

    /* now set the prop_tau on hybridization events and bidirection events.
       Hybridization nodes receive 1 (inner nodes) and 0 (mirored nodes).
       On bidirections only one of the two inner nodes receives prop_tau=1,
       the other node including the two mirror nodes receive 0 */
    for (i = 0; i < stree->hybrid_count; ++i)
    {
      snode_t * mnode = stree->nodes[stree->tip_count+stree->inner_count+i];
      snode_t * hnode = mnode->hybrid;
      assert(hnode);

      if (node_is_bidirection(hnode))
      {
        /* bidirection  */

        assert(hnode->node_index < stree->tip_count + stree->inner_count);

        assert((hnode->htau  && !mnode->htau &&  mnode->parent->htau && !mnode->parent->hybrid->htau) ||
               (!hnode->htau && !mnode->htau && !mnode->parent->htau && !mnode->parent->hybrid->htau));
        /* this means the bidirection was already processed from the other end node */
        if (hnode->htau) continue;

        /* if htau=0 for all four nodes involved in a bidirectional
           introgression then set htau=1 to inner node hnode. All other
           nodes keep htau=0 forever, and will share the tau of hnode */
        if (!hnode->htau && !mnode->htau && 
            !mnode->parent->htau && !mnode->parent->hybrid->htau)
        {
          /* Note: Changed the code to have both non-mirror nodes to htau=1.
             This is important for computing gene tree branch lengths with
             relaxed clock */
          hnode->htau = 1;
          mnode->parent->htau = 1;

          hnode->prop_tau = 1;
          mnode->prop_tau = 0;
          mnode->parent->prop_tau = 0;
          mnode->parent->hybrid->prop_tau = 0;
        }
        else
        {
          /* This should never occur */
          assert(0);
        }
      }
      else
      {
        /* hybridization */

        if (!hnode->htau)
          hnode->parent->prop_tau = 0;
        if (!mnode->htau)
          mnode->parent->prop_tau = 0;

        hnode->prop_tau = 1;
        mnode->prop_tau = 0;
      }
    }
  }
  else
  {
    /* not msci */
    for (i = 0; i < stree->inner_count; ++i)
      stree->nodes[stree->tip_count+i]->prop_tau = 1;
  }

  if (opt_msci && stree->root->prop_tau == 0)
    fatal("Error: species tree root requires a tau parameter [tau-parent=yes]");

  return stree;
}

static void network_assign_phis(stree_t * stree)
{
  long offset,i;

  offset = stree->tip_count+stree->inner_count;

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * mnode = stree->nodes[offset+i];
    snode_t * snode = mnode->hybrid;

    /* both edges have a phi */
    if (mnode->has_phi && snode->has_phi)
      fatal("[ERROR] Multiple phi definitions (defphi) found on node %s.\n"
            "[ERROR] Define the phi parameter only on one of the two branches.",
            snode->label);
    
    /* one of the edges has a phi */
    if (mnode->has_phi || snode->has_phi) continue;

    /* none of the edges has a phi. If just one of the edges is horizontal
       place phi on it, otherwise place it on the main edge */
    if (!mnode->htau && snode->htau)
      mnode->has_phi = 1;
    else
      snode->has_phi = 1;
  }
}

stree_t * bpp_parse_newick_string(const char * line)
{
  ntree_t * tree = NULL;
  stree_t * stree = NULL;
  list_t * token_list = NULL;

  #if 1
  /* old parser */
  char * s = xstrdup(line);

  if (!(token_list = parse_tree(s)))
    goto l_unwind;

  if (!(tree = syntax_parse(token_list)))
    goto l_unwind;

  stree = stree_from_ntree(tree);

  /* TODO: The PPTABLE is not yet allocated and thus MSci doesn't work */
  if (!opt_msci)
    stree_reset_leaves(stree);

  if (opt_msci)
  {
    /* postprocess to set phi parameters */
    network_assign_phis(stree);
  }

  #else
  /* old parser */

  stree_t * stree = stree_parse_newick_string(line);
  #endif

l_unwind:
  if (s)
    free(s);
  if (token_list)
  {
    list_clear(token_list,token_clear);
    free(token_list);
  }

  if (tree)
    ntree_destroy(tree,NULL);
  
  return stree;
}

ntree_t * bpp_parse_newick_string_ntree(const char * line)
{
  ntree_t * tree = NULL;
  list_t * token_list = NULL;

  char * s = xstrdup(line);

  if (!(token_list = parse_tree(s)))
    goto l_unwind;

  tree = syntax_parse(token_list);

l_unwind:
  if (token_list)
  {
    list_clear(token_list,token_clear);
    free(token_list);
  }
  free(s);

  return tree;
}
