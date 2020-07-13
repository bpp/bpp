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

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static void ntree_set_leaves_count_recursive(node_t * node)
{
  long i;

  if (!node->children_count)
  {
    node->leaves = 1;
    return;
  }

  node->leaves = 0;
  for (i = 0; i < node->children_count; ++i)
  {
    ntree_set_leaves_count_recursive(node->children[i]);
    node->leaves += node->children[i]->leaves;
  }
}

static void ntree_set_leaves_count(ntree_t * ntree)
{
  ntree_set_leaves_count_recursive(ntree->root);
}

static long is_emptyline(const char * line)
{
  size_t ws = strspn(line, " \t\r\n");
  if (!line[ws] || line[ws] == '*' || line[ws] == '#') return 1;
  return 0;
}

/* get string from current position pointed by line until a delimiter. Create
   new string with result, store it in value, and return number of characters
   read */
static long get_delstring(const char * line, const char * del, char ** value)
{
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

  /* skip all characters until a delimiter is found */
  char * end = start + strcspn(start, del);

  *end = 0;

  if (start==end)
  {
    free(s);
    return 0;
  }

  *value = xstrdup(start);

  free(s);
  return ws + end - start;
}

static long get_string(const char * line, char ** value)
{
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

  /* skip all characters except star and hash */
  while (*p && *p != '*' && *p != '#') p++;


  /* now go back, skipping all white-space until a character occurs */
  for (--p; *p == ' ' || *p == '\t' || *p == '\r' || *p =='\n'; --p);
  
  char * end = p+1;
  *end = 0;

  *value = xstrdup(start);
  free(s);

  return ws + end - start;
}

void constdefs_dealloc(void * data)
{
  if (data)
  {
    constdefs_t * defs = (constdefs_t *)data;

    if (defs->arg1)
      free(defs->arg1);
    if (defs->arg2)
      free(defs->arg2);

    free(defs);
  }
}

static long parse_define(const char * line, constdefs_t * defs)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tag = NULL;

  long count;

  /* get label of new group tag */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->arg1 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* get 'as' keyword */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;

  if (strcasecmp(tag,"as"))
    goto l_unwind;
  free(tag); tag = NULL;
  
  p += count;

  /* get newick string forming the group */
  count = get_string(p,&tag);
  if (!count) goto l_unwind;
  defs->arg2 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  ret = 1;

l_unwind:
  if (ret) ret = p-s;
  free(s);
  if (tag)
    free(tag);

  return ret;
}

static long parse_constraint(const char * line, constdefs_t * defs)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tree = NULL;

  long count;

  /* skip all whitespace and '=' if exists */
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
  if (*p == '=')
    ++p;

  /* get newick string forming the group */
  count = get_string(p,&tree);
  if (!count) goto l_unwind;
  defs->arg1 = xstrdup(tree);

  p += count;

  /* TODO: check tree notation using bpp_parse_newick_string_ntree(tree) */
  //assert(0);

  ret = 1;

l_unwind:
  if (ret) ret = p-s;
  free(s);
  if (tree)
    free(tree);

  return ret;
}

static constdefs_t * parse_constdefs(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  constdefs_t * defs = NULL;
  char * tag = NULL;

  defs = (constdefs_t *)xcalloc(1,sizeof(constdefs_t));

  long count;

  /* get entry type */
  count = get_delstring(p," \t\r\n", &tag);
  if (!count) goto l_unwind;
  
  if (!strcasecmp(tag,"constraint"))
  {
    defs->type = BPP_CONSTDEFS_CONSTRAINT;
  }
  else if (!strcasecmp(tag,"define"))
  {
    defs->type = BPP_CONSTDEFS_DEFINE;
  }
  else if (!strcasecmp(tag,"outgroup"))
  {
    defs->type = BPP_CONSTDEFS_OUTGROUP;
  }
  else
  {
    goto l_unwind;
  }
  free(tag); tag = NULL;

  p += count;

  if (defs->type == BPP_CONSTDEFS_CONSTRAINT ||
     defs->type == BPP_CONSTDEFS_OUTGROUP)
  {
    /* CONSTRAINT or OUTGROUP */

    count = parse_constraint(p, defs);
    if (!count) goto l_unwind;

    p += count;
  }
  else if (defs->type == BPP_CONSTDEFS_DEFINE)
  {
    /* DEFINE */

    count = parse_define(p, defs);
    if (!count) goto l_unwind;

    p += count;
  }

  if (is_emptyline(p)) ret = 1;

l_unwind:
  if (!ret)
  {
    constdefs_dealloc(defs);
    defs = NULL;
  }
  free(s);
  if (tag)
    free(tag);

  return defs;
}

list_t * parse_constfile(const char * constfile)
{
  long line_count = 0;
  long ret = 1;
  FILE * fp;
  list_t * list;
  constdefs_t * defs = NULL;

  fp = xopen(constfile,"r");

  list = (list_t *)xcalloc(1,sizeof(list_t));

  while (getnextline(fp))
  {
    ++line_count;

    if (is_emptyline(line)) continue;

    if (!(defs = parse_constdefs(line)))
    {
      ret = 0;
      fprintf(stderr, "Invalid entry in %s (line %ld)\n", constfile, line_count);
      goto l_unwind;
    }

    defs->lineno = line_count;

    /* insert item into list */
    list_append(list,(void *)defs);
  }

l_unwind:
  fclose(fp);
  if (ret == 0)
  {
    list_clear(list,constdefs_dealloc);
    free(list);
    list = NULL;
  }
  return list;
}

static hashtable_t * stree_hash(stree_t * stree)
{
  long i;
  hashtable_t * ht;

  ht = hashtable_create(stree->tip_count);

  /* index all tip labels in the hash table */
  for (i = 0; i < stree->tip_count; ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = stree->nodes[i]->label;
    pair->data = (void *)(uintptr_t)i;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(stree->nodes[i]->label),
                          cb_cmp_pairlabel))
      fatal("Duplicate taxon (%s)", stree->nodes[i]->label);
  }

  return ht;
}

static hashtable_t * ntree_hash(ntree_t * ntree)
{
  long i;
  hashtable_t * ht;

  ht = hashtable_create(ntree->tip_count);

  /* index all tip labels in the hash table */
  for (i = 0; i < ntree->tip_count; ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = ntree->leaves[i]->label;
    pair->data = (void *)(uintptr_t)i;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(ntree->leaves[i]->label),
                          cb_cmp_pairlabel))
      fatal("Duplicate taxon (%s)", ntree->leaves[i]->label);
  }

  return ht;
}

/* find LCA of tip nodes (labels) (n-tree version) */
static node_t * ntree_lca_nodes(ntree_t * ntree,
                                hashtable_t * ht,
                                char ** labels,
                                long label_count)
{
  long i,j;
  long ht_created = 0;
  node_t * lca = NULL;

  if (label_count > ntree->tip_count)
    return NULL;

  if (label_count == ntree->tip_count && label_count == 1)
  {
    if (!strcmp(ntree->leaves[0]->label,labels[0]))
      lca = ntree->leaves[0];
    return lca;
  }

  if (!ht)
  {
    ht = ntree_hash(ntree);
    ht_created = 1;
  }

  for (j = 0; j < label_count; ++j)
  {
    /* get next tip */
    char * taxon = xstrdup(labels[j]);
    
    /* remove prefix and suffix spaces */
    if (!strlen(taxon))
      fatal("Invalid taxon in list");
    
    /* find taxon in the hash list and obtain its index (i) */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    cb_cmp_pairlabel);
    if (!query)
    {
      free(taxon);
      goto l_unwind;
    }

    i = (unsigned int)(uintptr_t)(query->data);

    assert(i < ntree->tip_count);
    assert(!strcmp(ntree->leaves[i]->label, taxon));

    ntree->leaves[i]->mark = 1;

    free(taxon);
  }

  /* -*- Find LCA -*- */

  /* 1. Mark all root-paths of marked tip nodes */
  for (i = 0; i < ntree->tip_count; ++i)
  {
    if (ntree->leaves[i]->mark)
    {
      node_t * x = ntree->leaves[i]->parent;
      while (x)
      {
        x->mark = 1;
        x = x->parent;
      }
    }
  }

  /* 2. Descend from root and first first descendant (incl. root) that has both
         children marked */
  lca = ntree->root;
  while (1)
  {
    long mark_count = 0;
    node_t * next_child = NULL;
    assert(lca->children_count > 0);

    for (i = 0; i < lca->children_count; ++i)
      if (lca->children[i]->mark)
      {
        ++mark_count;
        next_child = lca->children[i];
      }

    if (mark_count > 1)
      break;

    assert(mark_count == 1);
    assert(next_child);

    lca = next_child;
  }

l_unwind:
  for (i = 0; i < ntree->tip_count; ++i)
    ntree->leaves[i]->mark = 0;
  for (i = 0; i < ntree->inner_count; ++i)
    ntree->inner[i]->mark = 0;

  if (ht_created)
  {
    hashtable_destroy(ht,free);
  }
  return lca;
}

/* find LCA of tip nodes (labels) (binary tree version) */
static snode_t * lca_nodes(stree_t * stree,
                           hashtable_t * ht,
                           char ** labels,
                           long label_count)
{
  long i,j;
  long ht_created = 0;

  if (!ht)
  {
    ht = stree_hash(stree);
    ht_created = 1;
  }

  for (j = 0; j < label_count; ++j)
  {
    /* get next tip */
    char * taxon = xstrdup(labels[j]);
    
    /* remove prefix and suffix spaces */
    if (!strlen(taxon))
      fatal("Invalid taxon in list");
    
    /* find taxon in the hash list and obtain its index (i) */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    cb_cmp_pairlabel);
    if (!query)
      fatal("Taxon %s does not appear in the species tree", taxon);

    i = (unsigned int)(uintptr_t)(query->data);

    assert(i < stree->tip_count);
    assert(!strcmp(stree->nodes[i]->label, taxon));

    stree->nodes[i]->mark[0] = 1;

    free(taxon);
  }

  /* -*- Find LCA -*- */

  /* 1. Mark all root-paths of marked tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    if (stree->nodes[i]->mark[0])
    {
      snode_t * x = stree->nodes[i]->parent;
      while (x)
      {
        x->mark[0] = 1;
        x = x->parent;
      }
    }
  }

  /* 2. Descend from root and first first descendant (incl. root) that has both
         children marked */
  snode_t * lca = stree->root;
  while (1)
  {
    assert(lca->left && lca->right);

    if (lca->left->mark[0] && lca->right->mark[0])
      break;

    assert(lca->left->mark[0] || lca->right->mark[0]);
    
    if (lca->left->mark[0])
      lca = lca->left;
    else
      lca = lca->right;
  }

  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
    stree->nodes[i]->mark[0] = 0;

  if (ht_created)
  {
    hashtable_destroy(ht,free);
  }
  return lca;
}

static void ntree_subtree_tiplabels_recursive(node_t * root,
                                              char ** labels,
                                              long * index)
{
  long i;

  if (!root->children_count)
  {
    labels[*index] = xstrdup(root->label);
    *index = *index + 1;

    return;
  }

  for (i = 0; i < root->children_count; ++i)
    ntree_subtree_tiplabels_recursive(root->children[i],labels,index);
}

/* get the labels of tips in the subtree rooted at node 'root' */
static char ** ntree_subtree_tiplabels(node_t * root)
{
  long index = 0;
  char ** labels = (char **)xmalloc((size_t)(root->leaves)*sizeof(char *));

  ntree_subtree_tiplabels_recursive(root, labels, &index);

  return labels;
}

/* we are given a binary species tree and an n-ary constraint tree. Check if
   the species tree has a subtree whose tips match exactly those of the
   constraint */
int is_subtree(stree_t * stree, ntree_t * ntree)
{
  long i,j;
  int ret = 0;
  assert(!opt_msci);
  assert(ntree->tip_count > 0);

  for (i = 0; i < ntree->inner_count; ++i)
  {
    char ** labels = ntree_subtree_tiplabels(ntree->inner[i]); 

    snode_t * lca = lca_nodes(stree,NULL,labels,ntree->inner[i]->leaves);
    assert(lca);

    for (j = 0; j < ntree->inner[i]->leaves; ++j)
      free(labels[j]);
    free(labels);

    if (lca->leaves != ntree->inner[i]->leaves)
      goto l_unwind;
  }

  ret = 1;
l_unwind:
  return ret;
}

static long tiplabel_exists(stree_t * stree, const char * label)
{
  long i;
  for (i = 0; i < stree->tip_count; ++i)
    if (!strcmp(stree->nodes[i]->label, label))
      break;

  return !(i == stree->tip_count);
}

/* given an n-ary constraint tree, a species tree and a list of definitions,
   find all tips in the constraint that do not appear in the species tree.
   Ensure that each such tip is actually a 'definition' that is an alias of
   another n-ary tree/constraint. Replace those tips with the corresponding
   trees */
static void ntree_replace_aliases(stree_t * stree,
                                  ntree_t ** ptrd,
                                  long lineno,
                                  char ** def_label,
                                  char ** def_string,
                                  long def_count)
{
  int needs_rewrap = 0;
  int new_tip_count = 0;
  int new_inner_count = 0;
  long i,j;
  ntree_t * d = *ptrd;

  new_tip_count = d->tip_count;
  new_inner_count = d->inner_count;

  /* check if all d tips exist in species tree. If not, check for previous
     definitions */
  for (i = 0; i < d->tip_count; ++i)
  {
    if (!tiplabel_exists(stree,d->leaves[i]->label))
    {
      for (j = 0; j < def_count; ++j)
        if (!strcmp(def_label[j],d->leaves[i]->label))
          break;
      if (j == def_count)
        fatal("Definition in %s (line %ld) contains undefined taxon (%s)",
               opt_constfile, lineno, d->leaves[i]->label);

      /* found definition. expand into a tree T */
      ntree_t * t = bpp_parse_newick_string_ntree(def_string[j]);
      assert(t);

      /* now we need to replace the leaf with the previous definition label to
         tree T */
      node_t * leaf = d->leaves[i];
      if (!leaf->parent)
      {
        assert(leaf == d->root);
        d->root = t->root;
      }
      else
      {
        for (j = 0; j < leaf->parent->children_count; ++j)
          if (leaf->parent->children[j] == leaf)
            break;
        assert(j < leaf->parent->children_count);

        t->root->parent = leaf->parent;
        leaf->parent->children[j] = t->root;
      }

      new_tip_count--;
      new_tip_count += t->tip_count;
      new_inner_count += t->inner_count;

      free(leaf->label);
      free(leaf);
      free(t->leaves);
      free(t->inner);
      free(t);
      needs_rewrap = 1;
    }
  }

  /* re-wrap tree structure, delete old and replace with new */
  if (needs_rewrap)
  {
    ntree_t * newtree = ntree_wraptree(d->root,new_tip_count,new_inner_count);
    free(d->leaves);
    free(d->inner);
    free(d);
    d = newtree;
    ntree_set_leaves_count(d);
    *ptrd = d;
  }
}

void constraint_mark_recursive(snode_t * node, long cvalue, long lineno)
{
  /* mark constraints in pre-order */

  if (node->constraint)
    return;

  node->constraint = cvalue;
  node->constraint_lineno = lineno;

  if (!node->left && !node->right)
    return;

  if (node->left)
    constraint_mark_recursive(node->left,cvalue,lineno);
  if (node->right)
    constraint_mark_recursive(node->right,cvalue,lineno);
}

void constraint_process_recursive(stree_t * t,
                                  node_t * cr,
                                  long * cvalue,
                                  long lineno)
{
  long i,j;

  if (cr->children_count == 0)
    return;

  for (i = 0; i < cr->children_count; ++i)
    constraint_process_recursive(t,cr->children[i],cvalue,lineno); 

    /* get labels of subtree */
  char ** labels = ntree_subtree_tiplabels(cr);

  snode_t * lca = lca_nodes(t,NULL,labels,cr->leaves);
  assert(lca);

  for (j = 0; j < cr->leaves; ++j)
    free(labels[j]);
  free(labels);

  /* ensure the direct descendants of lca have both either the same constraint
     or no constraint */
  if (lca->left && lca->right)
  {
    if (lca->left->constraint != lca->right->constraint)
      fatal("Conflicting constraints in file %s (lines %ld and %ld)",
            opt_constfile,
            lca->left->constraint ? 
              lca->left->constraint_lineno : lca->right->constraint_lineno,
            lineno);
  }

  /* mark all unmarked descendants in preorder traversal with a new unique value */
  *cvalue = *cvalue+1;

  if (lca->left)
    constraint_mark_recursive(lca->left,*cvalue,lineno);

  if (lca->right)
    constraint_mark_recursive(lca->right,*cvalue,lineno);
}

static int strip(char ** s)
{
  size_t ws;
  size_t len;
  char * x = *s;
  char * y = *s;

  ws = strspn(x," \t\r\n");
  x = x+ws;

  len = strlen(x);
  if (len < 1) return 0;

  y = x+len-1;
  while (*y == ' ' || *y == '\t' || *y == '\r' || *y == '\n')
    --y;
  ++y;
  *y = 0;

  y = xstrdup(x);
  free(*s);
  *s = y;

  return 1;
}

static char ** tokenize_csv(const char * csv, long * count)
{
  long comma_count = 0;
  long i;

  char * p = xstrdup(csv);
  char * s = p;

  /* count number of commas */
  for (i = 0; i < (long)strlen(s); ++i)
    if (s[i] == ',') ++comma_count;

  /* allocate space for tokenization */
  char ** tokens = (char **)xmalloc((size_t)(comma_count+1) * sizeof(char *));

  i = 0;
  while (*s)
  {
    size_t comma = 0;

    /* find next comma and replace it with a 0 to create a string of the
       current token */
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
        fatal("Comma separated list starts from a comma");
      else
        fatal("Consecutive comma symbols found");
    }

    char * taxon = xstrdup(s);
    if (!strip(&taxon))
      fatal("No taxon found");
    tokens[i++] = taxon;

    s += tokenlen+comma;
  }

  *count = i;

  free(p);

  return tokens;
}

static int constraints_merge(ntree_t ** t1ptr, ntree_t ** t2ptr)
{
  long i,j,k;
  char ** labels = NULL;

  ntree_t * t1 = *t1ptr;
  ntree_t * t2 = *t2ptr;

  /* We go through each internal node u of t2. Each such node u has n>1 children
     c1,c2,...,cn.

     We mark all rootpaths from leaves in L(u). We find the LCA. If the number of leaves on all marked
     children of the LCA equals to L(u) we are fine, otherwise error.
  */


  for (i = 0; i < t2->inner_count; ++i)
  {
    node_t * node = t2->inner[i];
    assert(node->leaves > 1);

    labels = ntree_subtree_tiplabels(node);
    node_t * lca = ntree_lca_nodes(t1,NULL,labels,node->leaves);
    assert(lca && lca->children_count > 1);
    
    /* mark the children of LCA lying on the rootpaths */
    for (j = 0; j < node->leaves; ++j)
    {
      for (k = 0; k < t1->tip_count; ++k)
        if (!strcmp(t1->leaves[k]->label,labels[j]))
          break;
      if (k == t1->tip_count)
        fatal("Internal error");

      node_t * temp = t1->leaves[k];

      /* mark children of LCA in rootpaths */
      while (temp)
      {
        if (temp->parent == lca)
        {
          temp->mark = 1;
          break;
        }

        temp = temp->parent;
      }
      assert(temp && temp->parent && temp->parent == lca);
    }
    for (j = 0; j < node->leaves; ++j)
      free(labels[j]);
    free(labels);

    /* now check compatibility. Leaves of marked children must be equal to the
       leaves of t2->inner[i] */
    long leaves = 0;
    long marked_children_count = 0;
    for (j = 0; j < lca->children_count; ++j)
    {
      if (lca->children[j]->mark)
      {
        leaves += lca->children[j]->leaves;
        marked_children_count++;
      }
    }
    if (leaves != node->leaves)
    {
      for (j = 0; j < lca->children_count; ++j)
        lca->children[j]->mark = 0;
      return 0;                         /* incompatible */
    }

    /* further resolve multifurcation */
    if (lca->children_count > marked_children_count)
    {
      node_t * new = (node_t *)xcalloc(1,sizeof(node_t));
      new->children_count = marked_children_count;
      new->children = (node_t **)xmalloc((size_t)marked_children_count *
                                          sizeof(node_t *));
      long dst = 0;
      long src = 0;
      for (j = 0; j < lca->children_count; ++j)
      {
        if (lca->children[j]->mark)
        {
          lca->children[j]->mark = 0;
          new->children[dst++] = lca->children[j];
          lca->children[j]->parent = new;
        }
        else
        {
          if (src != j)
            lca->children[src] = lca->children[j];

          src++;
        }
      }
      assert(src == (lca->children_count - marked_children_count));

      lca->children[src++] = new;
      new->parent = lca;
      lca->children_count = src;
      new->children_count = marked_children_count;

      free(t1->leaves); free(t1->inner);
      ntree_t * newtree = ntree_wraptree(t1->root,t1->tip_count,t1->inner_count+1);
      free(t1);
      t1 = newtree;
      ntree_set_leaves_count(t1);

      *t1ptr = t1;
    }
    else
    {
      /* clean marks */
      for (j = 0; j < lca->children_count; ++j)
        lca->children[j]->mark = 0;
    }

  }
  return 1;     /* compatible */
}

static int constraints_reduce(ntree_t ** t1ptr,
                              ntree_t ** t2ptr,
                              const char * c1,
                              const char * c2,
                              long line1,
                              long line2)
{
  long i,j;
  long common = 0;
  long diff = 0;

  ntree_t * t1 = *t1ptr;
  ntree_t * t2 = *t2ptr;

  /* compute union and difference */
  for (i = 0; i < t1->tip_count; ++i)
  {
    for (j = 0; j < t2->tip_count; ++j)
    {
      if (!strcmp(t1->leaves[i]->label,t2->leaves[j]->label))
        break;
    }
    if (j == t2->tip_count)
      ++diff;
    else
      ++common;
  }

  if (common == t1->tip_count)
  {
    /* possibly compatible */
    if (!constraints_merge(t2ptr,t1ptr))
      fatal("Conflicting constraints found:\n"
            "  Line %2ld: %s\n"
            "  Line %2ld: %s",
            line1,c1,line2,c2);
    return 0;
  }
  else if (common == t2->tip_count)
  {
    /* possibly compatible */
    if (!constraints_merge(t1ptr,t2ptr))
      fatal("Conflicting constraints found:\n"
            "  Line %2ld: %s\n"
            "  Line %2ld: %s",
            line1,c1,line2,c2);
    return 1;
  }
  else if (common && diff)
  {
    /* conflicting */
    return 2;
  }
  else if (!common && diff)
  {
    /* independent */
      return 3;
  }
  else
  {
    assert(0);
  }

  return 4;
}

static void complement_definition(stree_t * stree,
                                  ntree_t ** ptrd,
                                  constdefs_t * def)
{
  long i,j;
  long taxa_count;
  long space;
  size_t len;
  char * newdef;
  char * p;

  ntree_t * d = *ptrd;

  /* mark tips in the species tree */
  for (i = 0; i < d->tip_count; ++i)
  {
    for (j = 0; j < stree->tip_count; ++j)
      if (!strcmp(d->leaves[i]->label,stree->nodes[j]->label))
      {
        stree->nodes[j]->mark[0] = 1;
        break;
      }
    if (j == stree->tip_count)
      fatal("Invalid definition in file %s (line %ld)",
            opt_constfile, def->lineno);
  }

  space = taxa_count = 0;
  for (i = 0; i < stree->tip_count; ++i)
    if (!stree->nodes[i]->mark[0])
    {
      space += strlen(stree->nodes[i]->label);
      taxa_count++;
    }
  assert(taxa_count == stree->tip_count - d->tip_count);

  /* one comma less than number of taxa plus start/end parenthesis plus
     terminating zero plus semicolon */
  space += taxa_count+3;
  newdef = (char *)xcalloc((size_t)space,sizeof(char));

  p = newdef;
  taxa_count = 0;
  *p++ = '(';
  for (i = 0; i < stree->tip_count; ++i)
  {
    if (!stree->nodes[i]->mark[0])
    {
      if (taxa_count)
        *p++ = ',';
      
      ++taxa_count;
      len = strlen(stree->nodes[i]->label);
      memcpy(p,stree->nodes[i]->label,len);
      p += len;
    }
    else
      stree->nodes[i]->mark[0] = 0;
  }
  *p++ = ')';
  *p++ = ';';
  *p++ = '\0';

  /* delete old tree and definition */
  free(def->arg2);
  ntree_destroy(d,NULL);

  def->arg2 = newdef;
  d = bpp_parse_newick_string_ntree(def->arg2);
  ntree_set_leaves_count(d);
  *ptrd = d;
}

void definition_process(stree_t * stree,
                        ntree_t ** ptrd,
                        constdefs_t * def,
                        char ** def_label,
                        char ** def_string,
                        long def_count,
                        int cpl)
{
  /* check that all d tips exist in species tree, and replace aliases with
     previous definitions */
  ntree_replace_aliases(stree,ptrd,def->lineno,def_label,def_string,def_count);

  /* complement */
  if (cpl)
    complement_definition(stree,ptrd,def);

  ntree_t * d = *ptrd;

  /* check if constraint tree is a subtree of stree */
  if (!is_subtree(stree,d))
    fatal("Invalid definition in file %s (line %ld)", opt_constfile, def->lineno);

  def_label[def_count]  = xstrdup(def->arg1);
  def_string[def_count] = ntree_export_newick(d,0);
}

static void definitions_expand(list_t * constlist,
                               stree_t * stree,
                               char ** def_labels,
                               char ** def_expand,
                               FILE * fp_out)
{
  list_item_t * li;
  int cpl = 0;          /* complement */
  long def_count = 0;

  li = constlist->head;
  while (li)
  {
    constdefs_t * def = (constdefs_t *)(li->data);

    list_item_t * tmp = li->next;

    if (def->type == BPP_CONSTDEFS_DEFINE)
    {
      cpl = 0;
      ntree_t * t = NULL;

      if (strlen(def->arg2) > 3 && !strncasecmp(def->arg2, "not", 3))
      {
        t = bpp_parse_newick_string_ntree(def->arg2+3);
        cpl = 1;
      }
      else
        t = bpp_parse_newick_string_ntree(def->arg2);
      if (!t)
        fatal("Error while parsing definition (line %ld)", def->lineno);
      ntree_set_leaves_count(t);
      definition_process(stree,&t,def,def_labels, def_expand,def_count,cpl);
      
      /* expanded form 
         printf("    alias: %s\n", def_string[def_count]);
      */
      def_count++;
      ntree_destroy(t,NULL);

      fprintf(stdout,
              "Definition %2ld (line %ld): %s = %s\n",
              def_count, def->lineno, def->arg1, def->arg2);
      fprintf(fp_out,
              "Definition %2ld (line %ld): %s = %s\n",
              def_count, def->lineno, def->arg1, def->arg2);

      list_delitem(constlist,li,constdefs_dealloc);
    }

    li = tmp;
  }
}

static void constraints_validate(ntree_t * tree, long lineno)
{
  long i,j;

  /* check for duplicate tips */
  for (i = 0; i < tree->tip_count; ++i)
    for (j = i+1; j < tree->tip_count; ++j)
      if (!strcmp(tree->leaves[i]->label, tree->leaves[j]->label))
        fatal("Constraint on line %ld contains duplicate taxon label (%s)",
              lineno, tree->leaves[i]->label);
  
  /* check for unary nodes */
  for (i = 0; i < tree->inner_count; ++i)
    if (tree->inner[i]->children_count == 1)
      fatal("Constraint on line %ld contrains an unary node", lineno);
}

static void constraints_process(stree_t * stree,
                                list_t * constlist,
                                char ** def_labels,
                                char ** def_string,
                                long def_count,
                                FILE * fp_out)
{
  long i,j;
  long const_count = 0;
  long * lines;
  list_item_t * li;
  ntree_t ** trees;
  list_item_t ** liptr;

  /* find number of constraints */
  for (li = constlist->head; li; li = li->next)
  {
    constdefs_t * def = (constdefs_t *)(li->data);
    if (def->type == BPP_CONSTDEFS_CONSTRAINT)
      const_count++;
  }

  /* serialize list of constraints */
  trees = (ntree_t **)xmalloc((size_t)const_count * sizeof(ntree_t *));
  liptr = (list_item_t **)xmalloc((size_t)const_count * sizeof(list_item_t *));
  lines = (long *)xmalloc((size_t)const_count * sizeof(long));

  /* parse tree constraints, replace aliases and store list pointers */
  for (i=0,li = constlist->head; li; li = li->next)
  {
    constdefs_t * def = (constdefs_t *)(li->data);

    if (def->type == BPP_CONSTDEFS_CONSTRAINT)
    {
      trees[i] = bpp_parse_newick_string_ntree(def->arg1);
      liptr[i] = li;
      lines[i] = def->lineno;
      ntree_set_leaves_count(trees[i]);
      ntree_replace_aliases(stree,
                            trees+i,
                            def->lineno,
                            def_labels,
                            def_string,
                            def_count);
      ++i;
    }
  }

  /* check for duplicate tips or unary nodes */
  for (i = 0; i < const_count; ++i)
    constraints_validate(trees[i],lines[i]);

  fprintf(stdout, "List of constraints in file %s\n", opt_constfile);
  fprintf(fp_out, "List of constraints in file %s\n", opt_constfile);
  for (i = 0; i < const_count; ++i)
  {
    constdefs_t * def = (constdefs_t *)(liptr[i]->data);
    fprintf(stdout, "Constraint %2ld (line %2ld): %s\n",
           i+1, def->lineno, def->arg1);
    fprintf(fp_out, "Constraint %2ld (line %2ld): %s\n",
           i+1, def->lineno, def->arg1);

    if (!is_subtree(stree,trees[i]))
      fatal("Starting species tree contradicts constraint in file %s (line %ld)",
            opt_constfile, lines[i]);
  }

  /* pairwisely compare/reduce/merge constraints. Constraints removed due to
     merging are delete from the arrays and replaced with NULL, hence the
     break/continue in the loops */
  for (i = 0; i < const_count ;++i)
  {
    for (j = i+1; j < const_count; ++j)
    {
      if (!trees[i]) break;
      if (!trees[j]) continue;

      char * c1 = ntree_export_newick(trees[i],0);
      char * c2 = ntree_export_newick(trees[j],0);
      char * mc = NULL;
      int rc = constraints_reduce(trees+i,trees+j,c1,c2,lines[i],lines[j]);
      switch (rc)
      {
        case 0:
          mc = ntree_export_newick(trees[j],0);
          fprintf(stdout, "Merged constraints %s and %s -> %s\n", c1, c2, mc);
          fprintf(fp_out, "Merged constraints %s and %s -> %s\n", c1, c2, mc);
          free(mc);
          ntree_destroy(trees[i],NULL);
          mc = NULL; trees[i] = NULL;
          break;
        case 1:
          mc = ntree_export_newick(trees[i],0);
          fprintf(stdout, "Merged constraints %s and %s -> %s\n", c1, c2, mc);
          fprintf(fp_out, "Merged constraints %s and %s -> %s\n", c1, c2, mc);
          free(mc);
          ntree_destroy(trees[j],NULL);
          mc = NULL; trees[j] = NULL;
          break;
        case 2:
          #if 0
          printf("Constraints %s and %s are conflicting\n", c1, c2);
          #endif
          break;
        case 3:
          #if 0
          printf("Constraints %s and %s are independent\n", c1, c2);
          #endif
          break;
        default:
          fatal("Internal error - invalid code");
          break;
      }
      free(c1);
      free(c2);
    }
  }

  long cvalue = 0;

  /* apply constraints */
  for (i = 0; i < const_count; ++i)
  {
    if (!trees[i]) continue;

    char * stmp = ntree_export_newick(trees[i],0);
    fprintf(stdout, "Applying constraint: %s\n", stmp);
    fprintf(fp_out, "Applying constraint: %s\n", stmp);
    free(stmp);

    /* TODO: This check is probably no longer necessary */
    if (!is_subtree(stree,trees[i]))
      fatal("Starting species tree contradicts constraint in file %s (line %ld)",
            opt_constfile, lines[i]);

    constraint_process_recursive(stree,trees[i]->root,&cvalue,lines[i]);
  }

  for (i = 0; i < const_count; ++i)
    if (trees[i])
      ntree_destroy(trees[i],NULL);

  free(trees);
  free(liptr);
  free(lines);

  for (i = 0; i < def_count; ++i)
  {
    free(def_labels[i]);
    free(def_string[i]);
  }
  free(def_labels);
  free(def_string);
}

static void convert_outgroup_to_constraint(list_item_t * li,
                                           stree_t * stree,
                                           FILE * fp_out)
{
  long i,j;
  size_t len;
  char * p;

  constdefs_t * def = (constdefs_t *)(li->data);

  assert(def->type == BPP_CONSTDEFS_OUTGROUP);

  def->type = BPP_CONSTDEFS_CONSTRAINT;

  /* remove trailing whitespace and semicolon */
  len = strlen(def->arg1);
  if (len < 1)
    fatal("Invalid outgroup definition (line %ld)", def->lineno);
  p = def->arg1+len-1;
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
    --p;
  if (*p != ';') 
    ++p;
  *p = 0;
  
  long count = 0;
  char ** labels = tokenize_csv(def->arg1,&count);
  if (!labels || !count)
    fatal("No labels found in outgroup definition (line %ld)", def->lineno);

  /* check that all taxa in the array appear in the species tree */
  for (i = 0; i < count; ++i)
    if (!tiplabel_exists(stree,labels[i]))
      fatal("Invalid taxon %s in outgroup definition (line %ld)",
            labels[i], def->lineno);

  /* outgroup taxa must be at most two less than the tips in the tree */
  if (count >= stree->tip_count-1)
    fatal("At least two taxa must be part of the ingroup (line %ld)",
          def->lineno);

  for (i = 0; i < count; ++i)
  {
    for (j = 0; j < stree->tip_count; ++j)
      if (!strcmp(labels[i],stree->nodes[j]->label))
      {
        stree->nodes[j]->mark[0] = 1;
        break;
      }
    assert(j < stree->tip_count);
  }

  /* now construct a comma-separated sequence of the complementary taxa */

  /* first count the space we will need */
  long taxa_count = 0;
  long space = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    if (!stree->nodes[i]->mark[0])
    {
      space += strlen(stree->nodes[i]->label);
      taxa_count++;
    }
  } 
  /* one comma less than number of taxa plus start/end parenthesis plus
     terminating zero */
  space += taxa_count+2;
  char * constraint = (char *)xcalloc((size_t)space,sizeof(char));

  p = constraint;
  taxa_count = 0;
  *p++ = '(';
  for (i = 0; i < stree->tip_count; ++i)
  {
    if (!stree->nodes[i]->mark[0])
    {
      if (taxa_count)
        *p++ = ',';
      
      ++taxa_count;
      len = strlen(stree->nodes[i]->label);
      memcpy(p,stree->nodes[i]->label,len);
      p += len;
    }
  }
  *p++ = ')';
  *p++ = '\0';
  fprintf(stdout, "Outgroup: %s\n", def->arg1);
  fprintf(fp_out, "Outgroup: %s\n", def->arg1);
  fprintf(stdout, "Converting outgroup to constraint: %s\n", constraint);
  fprintf(fp_out, "Converting outgroup to constraint: %s\n", constraint);
  free(def->arg1);
  def->arg1 = constraint;

  /* reset marks */
  for (i = 0; i < stree->tip_count; ++i)
    stree->nodes[i]->mark[0] = 0;

  for (i = 0; i < count; ++i)
    free(labels[i]);
  free(labels);
}

void parse_and_set_constraints(stree_t * stree, FILE * fp_out)
{
  long def_count = 0;
  long og_count = 0;
  long const_count = 0;

  list_item_t * outgroup_li = NULL;

  list_t * constlist = parse_constfile(opt_constfile);
  if (!constlist)
    fatal("Invalid syntax found in file %s", opt_constfile);

  if (constlist->count && stree->tip_count < 2)
    fatal("Constraints require  species tree of more than 2 species");


  /* Go through all constraints, count them and check their syntax is correct */
  list_item_t * li = constlist->head;
  while (li)
  {
    constdefs_t * def = (constdefs_t *)(li->data);

    if (def->type == BPP_CONSTDEFS_DEFINE)
    {
      if (tiplabel_exists(stree,def->arg1))
        fatal("Definition %s in %s (line %ld) already exists as a taxon",
              def->arg1, opt_constfile, def->lineno);

      ntree_t * t = NULL;
      if (strlen(def->arg2) > 3 && !strncasecmp(def->arg2, "not", 3))
        t = bpp_parse_newick_string_ntree(def->arg2+3);
      else
        t = bpp_parse_newick_string_ntree(def->arg2);
      if (!t)
        fatal("Definition in %s (line %ld) is not a valid tree",
              opt_constfile, def->lineno);

      ntree_destroy(t,NULL);

      def_count++;

    }
    else if (def->type == BPP_CONSTDEFS_OUTGROUP)
    {
      outgroup_li = li;
      og_count++;
    }
    else if (def->type == BPP_CONSTDEFS_CONSTRAINT)
    {
      ntree_t * t = bpp_parse_newick_string_ntree(def->arg1);
      if (!t)
        fatal("Constraint in %s (line %ld) is not a valid tree", def->lineno);
      ntree_destroy(t,NULL);
      const_count++;
    }

    li = li->next;
  }

  opt_constraint_count = const_count + og_count;

  if (og_count > 1)
    fatal("Constraint file %s contains more than one outgroup definitions",
          opt_constfile);
  
  /* we want the outgroup to be processed last */
  if (outgroup_li)
    convert_outgroup_to_constraint(outgroup_li,stree,fp_out);

  char ** def_labels = (char **)xcalloc((size_t)def_count,sizeof(char *));
  char ** def_string = (char **)xcalloc((size_t)def_count,sizeof(char *));
  definitions_expand(constlist,stree,def_labels,def_string,fp_out);

  /* remove redundant constraints */
  constraints_process(stree,constlist,def_labels,def_string,def_count,fp_out);

  /* set constraints */
  list_clear(constlist,constdefs_dealloc);
  free(constlist);
}
