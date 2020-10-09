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

static int show_branches = 0;

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

void mscidefs_dealloc(void * data)
{
  if (data)
  {
    mscidefs_t * defs = (mscidefs_t *)data;

    if (defs->node1_1)
      free(defs->node1_1);
    if (defs->node1_2)
      free(defs->node1_2);
    if (defs->node2_1)
      free(defs->node2_1);
    if (defs->node2_2)
      free(defs->node2_2);
    if (defs->label1)
      free(defs->label1);
    if (defs->label2)
      free(defs->label2);

    free(defs);
  }
}

static long parse_tree(const char * line, mscidefs_t * defs)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tag = NULL;

  long count;

  count = get_string(p, &tag);
  if (!count) goto l_unwind;
  defs->node1_1 = xstrdup(tag);
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

static long parse_define(const char * line, mscidefs_t * defs)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tag = NULL;

  long count;

  /* get label of LCA */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->node1_1 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* get 'as' keyword */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;

  if (strcasecmp(tag,"as"))
    goto l_unwind;
  free(tag); tag = NULL;
  
  p += count;

  /* get all node labels forming the LCA */
  count = get_string(p,&tag);
  if (!count) goto l_unwind;
  defs->label1 = xstrdup(tag);
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

static long parse_hybrid(const char * line, mscidefs_t * defs)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tag = NULL;

  long count;

  /* denotes whether the hybridization consists of parallel edges */
  int parallel = 0;

  /* get node 1 end point 1 */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->node1_1 = xstrdup(tag);
  free(tag);

  p += count;

  /* get node 1 end point 2 */

  count = get_delstring(p, " \t\r\n,", &tag);
  if (!count) goto l_unwind;
  defs->node1_2 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* skip all whitespace */
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

  if (*p == ',')
  {
    ++p;

    /* get node 2 end point 1 */

    count = get_delstring(p, " \t\r\n", &tag);
    if (!count) goto l_unwind;
    defs->node2_1 = xstrdup(tag);
    free(tag); tag = NULL;

    p += count;

    /* get node 2 end point 2 */

    count = get_delstring(p, " \t\r\n", &tag);
    if (!count) goto l_unwind;
    defs->node2_2 = xstrdup(tag);
    free(tag); tag = NULL;

    p += count;
  }
  else
  {
    parallel = 1;
  }

  /* get 'as' keyword */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;

  if (strcasecmp(tag,"as"))
    goto l_unwind;
  free(tag); tag = NULL;
  
  p += count;

  /* get label for first node */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->label1 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* get label for second node */
  count = get_delstring(p," \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->label2 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  if (!parallel)
  {
    /* get 'tau' keyword */
    count = get_delstring(p, " \t\r\n=", &tag);
    if (!count) goto l_unwind;

    if (strcasecmp(tag,"tau"))
      goto l_unwind;
    free(tag); tag = NULL;
    
    p += count;

    /* skip all whitespace */
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

    if (*p != '=')
      goto l_unwind;
    ++p;

    /* get first tau value */
    count = get_delstring(p, " \t\r\n,", &tag);
    if (!count) goto l_unwind;
    if (strcasecmp(tag,"yes") && strcasecmp(tag,"no"))
      goto l_unwind;
    if (!strcasecmp(tag,"yes"))
      defs->has_tau1 = 1;
    else
      defs->has_tau1 = 0;
    free(tag); tag = NULL;

    p += count;
    
    /* skip all whitespace */
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

    if (*p != ',')
      goto l_unwind;
    ++p;

    /* get second tau value */
    count = get_delstring(p, " \t\r\n", &tag);
    if (!count) goto l_unwind;
    if (strcasecmp(tag,"yes") && strcasecmp(tag,"no"))
      goto l_unwind;
    if (!strcasecmp(tag,"yes"))
      defs->has_tau2 = 1;
    else
      defs->has_tau2 = 0;
    free(tag); tag = NULL;

    p += count;
  }

  /* read optional phi values */
  if (!is_emptyline(p))
  {
    /* get 'phi' keyword */
    count = get_delstring(p, " \t\r\n=", &tag);
    if (!count) goto l_unwind;

    if (strcasecmp(tag,"phi"))
      goto l_unwind;
    free(tag); tag = NULL;
    
    p += count;

    /* skip all whitespace */
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

    if (*p != '=')
      goto l_unwind;
    ++p;

    /* get phi value */
    double tmp = 0;
    count = get_double(p, &tmp);
    if (!count) goto l_unwind;
    if (tmp < 0)
      goto l_unwind;
    defs->phi1 = tmp;

    p += count;
  }

  ret = 1;

l_unwind:
  if (ret) ret = p-s;
  free(s);
  if (tag)
    free(tag);

  return ret;
}

static long parse_showbl(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tag = NULL;

  long count;

  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;

  p += count;

  if (strcasecmp(tag,"bl") || !is_emptyline(p))
    goto l_unwind;

  ret = 1;
  show_branches = 1;

l_unwind:
  if (ret) ret = p-s;
  free(s);
  if (tag)
    free(tag);

  return ret;
}

static long parse_bidir(const char * line, mscidefs_t * defs)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * tag = NULL;

  long count;

  /* get node 1 end point 1 */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->node1_1 = xstrdup(tag);
  free(tag);

  p += count;

  /* get node 1 end point 2 */

  count = get_delstring(p, " \t\r\n,", &tag);
  if (!count) goto l_unwind;
  defs->node1_2 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* skip all whitespace */
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

  if (*p != ',')
    goto l_unwind;
  ++p;

  /* get node 2 */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->node2_1 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* get node 2 end point 2 */

  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->node2_2 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* get 'as' keyword */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;

  if (strcasecmp(tag,"as"))
    goto l_unwind;
  free(tag); tag = NULL;
  
  p += count;

  /* get label for first node */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->label1 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  /* get label for second node */
  count = get_delstring(p, " \t\r\n", &tag);
  if (!count) goto l_unwind;
  defs->label2 = xstrdup(tag);
  free(tag); tag = NULL;

  p += count;

  if (!is_emptyline(p))
  {
    /* get 'phi' keyword */
    count = get_delstring(p, " \t\r\n=", &tag);
    if (!count || strcasecmp(tag,"phi")) goto l_unwind;
    free(tag); tag = NULL;
    
    p += count;

    /* skip all whitespace */
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

    if (*p != '=') goto l_unwind;
    ++p;

    /* get first phi value */
    count = get_delstring(p, " \t\r\n,", &tag);
    if (!count) goto l_unwind;

    p += count;

    double tmp = 0;
    count = get_double(tag, &tmp);
    if (!count || tmp < 0) goto l_unwind;
    defs->phi1 = tmp;

    /* skip all whitespace */
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

    if (*p != ',') goto l_unwind;
    ++p;

    /* get second phi value */
    count = get_double(p, &tmp);
    if (!count || tmp < 0) goto l_unwind;
    defs->phi2 = tmp;

    p += count;
  }

  ret = 1;

l_unwind:
  if (ret) ret = p-s;
  free(s);
  if (tag)
    free(tag);

  return ret;
}


static mscidefs_t * parse_mscidefs(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  mscidefs_t * defs = NULL;
  char * tag = NULL;

  defs = (mscidefs_t *)xcalloc(1,sizeof(mscidefs_t));

  defs->phi1 = defs->phi2 = -1;

  long count;

  /* get entry type */
  count = get_delstring(p," \t\r\n", &tag);
  if (!count) goto l_unwind;
  
  if (!strcasecmp(tag,"tree"))
  {
    defs->type = BPP_MSCIDEFS_TREE;
  }
  else if (!strcasecmp(tag,"define"))
  {
    defs->type = BPP_MSCIDEFS_DEFINE;
  }
  else if (!strcasecmp(tag,"hybridization"))
  {
    defs->type = BPP_MSCIDEFS_HYBRID;
  }
  else if (!strcasecmp(tag,"bidirection"))
  {
    defs->type = BPP_MSCIDEFS_BIDIR;
  }
  else if (!strcasecmp(tag,"show"))
  {
    defs->type = BPP_MSCIDEFS_SHOWBL;
  }
  else
  {
    goto l_unwind;
  }
  free(tag); tag = NULL;

  p += count;

  if (defs->type == BPP_MSCIDEFS_TREE)
  {
    /* TREE */

    count = parse_tree(p, defs);
    if (!count) goto l_unwind;

    p += count;
  }
  if (defs->type == BPP_MSCIDEFS_DEFINE)
  {
    /* DEFINE */

    count = parse_define(p, defs);
    if (!count) goto l_unwind;

    p += count;

  }
  else if (defs->type == BPP_MSCIDEFS_HYBRID)
  {
    /* HYBRIDIZATION */

    count = parse_hybrid(p, defs);
    if (!count) goto l_unwind;

    p += count;

  }
  else if (defs->type == BPP_MSCIDEFS_BIDIR)
  {
    /* BIDIRECTION */

    count = parse_bidir(p, defs);
    if (!count) goto l_unwind;

    p += count;
  }
  else if (defs->type == BPP_MSCIDEFS_SHOWBL)
  {
    /* SHOWBL */

    count = parse_showbl(p);
    if (!count) goto l_unwind;

    p += count;
  }

  if (is_emptyline(p)) ret = 1;

l_unwind:
  if (!ret)
  {
    mscidefs_dealloc(defs);
    defs = NULL;
  }
  free(s);
  if (tag)
    free(tag);

  return defs;
}


list_t * parse_mscifile(const char * mscifile)
{
  long line_count = 0;
  long ret = 1;
  FILE * fp;
  list_t * list;
  mscidefs_t * defs = NULL;
  long treeline = -1;

  fp = xopen(mscifile,"r");

  list = (list_t *)xcalloc(1,sizeof(list_t));

  while (getnextline(fp))
  {
    ++line_count;

    if (is_emptyline(line)) continue;

    if (!(defs = parse_mscidefs(line)))
    {
      ret = 0;
      fprintf(stderr, "Invalid entry in %s (line %ld)\n", mscifile, line_count);
      goto l_unwind;
    }

    defs->lineno = line_count;

    if (defs->type == BPP_MSCIDEFS_TREE)
    {
      if (treeline != -1)
        fatal("Duplicate tree in file %s (lines %ld and %ld)",
              mscifile, treeline, line_count);
      else
        treeline = line_count;
    }

    /* insert item into list */
    list_append(list,(void *)defs);
  }

l_unwind:
  fclose(fp);
  if (ret == 0)
  {
    list_clear(list,mscidefs_dealloc);
    free(list);
    list = NULL;
  }
  return list;
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


static snode_t * find_lca(stree_t * stree, hashtable_t * ht, const char * nodelist)
{
  char * s = xstrdup(nodelist);
  long i,k = 0;

  if (!strip(&s))
    fatal("No taxa found");

  char * p = s;


  while (*s)
  {
    /* get next tip */
    size_t taxon_len = strcspn(s,",");
    char * taxon = xstrndup(s, taxon_len);
    
    /* remove prefix and suffix spaces */
    if (!strip(&taxon))
      fatal("No taxon found");
    
    /* find taxon in the hash list and obtain its index (i) */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    cb_cmp_pairlabel);
    if (!query)
      fatal("Taxon %s does not appear in the tree", taxon);

    i = (unsigned int)(uintptr_t)(query->data);

    assert(i < stree->tip_count);
    assert(!strcmp(stree->nodes[i]->label, taxon));

    stree->nodes[i]->mark[0] = 1;

    s += taxon_len;
    if (*s == ',')
      s += 1;

    ++k;
    free(taxon);
  }
  free(p);

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

  return lca;
}

static void label_inner_nodes(list_t * list, stree_t * stree)
{
  long i;
  snode_t * lca = NULL;
  list_item_t * li;

  /* create hash table */
  hashtable_t * ht = hashtable_create(stree->tip_count);

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

  /* go through list of definitions */
  li = list->head;
  while (li)
  {
    mscidefs_t * def = (mscidefs_t *)(li->data);
    
    if (def->type == BPP_MSCIDEFS_DEFINE)
    {
      /* i. find LCA */
      lca = find_lca(stree, ht, def->label1);

      /* ii. label LCA */
      if (lca)
      {
        if (lca->label)
          fatal("Inner node to be labeled %s (line %ld) already has a label (%s)",
                def->node1_1, def->lineno, lca->label);

        lca->label = xstrdup(def->node1_1);
      }
      else
        fatal("LCA not found for list of taxa (line %ld): %s",
              def->lineno, def->label1);

    }
    li = li->next;
  }

  hashtable_destroy(ht,free);
}

static void destroy_pptable(stree_t * stree)
{
  unsigned int i;
  unsigned int nodes_count;

  if (!stree->pptable)
    return;

  nodes_count = stree->tip_count+stree->inner_count+stree->hybrid_count;

  for (i = 0; i < nodes_count; ++i)
    free(stree->pptable[i]);
  free(stree->pptable);
  stree->pptable = NULL;
}

static snode_t * edge_basenode(stree_t * stree,
                               const char * ep1,
                               const char * ep2,
                               long lineno)
{
  long i,j;

  /* get edge1 first end point (node1_1) */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (!stree->nodes[i]->label) continue;
    if (!strcmp(stree->nodes[i]->label,ep1))
      break;
  }
  if (i == stree->tip_count+stree->inner_count)
    fatal("Node with label %s does not exist (line %ld)",
          ep1, lineno);

  /* get edge1 second end point (node1_2) */
  for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
  {
    if (!stree->nodes[j]->label) continue;
    if (!strcmp(stree->nodes[j]->label,ep2))
      break;
  }
  if (j == stree->tip_count+stree->inner_count)
    fatal("Node with label %s does not exist (line %ld)",
          ep2, lineno);

  /* resolve parent child relationship */
  snode_t * a = stree->nodes[i];
  snode_t * b = stree->nodes[j];

  /* first check for parent-child relation */
  snode_t * p = NULL;  /* parent */
  snode_t * c = NULL;  /* child */
  if (a->hybrid && b->hybrid)
  {
    /* Note: mirror nodes do not have descendants, therefore we do not check
       those cases */
    if (a->parent == b)
    {
      c = a; p = b;
    }
    else if (a->hybrid->parent == b)
    {
      c = a->hybrid; p = b;
    }
    else if (b->parent == a)
    {
      c = b; p = a;
    }
    else if (b->hybrid->parent == a)
    {
      c = b->hybrid; p = a;
    }
    else
      fatal("Nodes %s and %s do not form an edge (line %ld)", ep1, ep2, lineno);
  }
  else if (a->hybrid)
  {
    if (a->parent == b)
    {
      c = a; p = b;
    }
    else if (a->hybrid->parent == b)
    {
      c = a->hybrid; p = b;
    }
    else if (b->parent == a)
    {
      c = b; p = a;
    }
    else
      fatal("Nodes %s and %s do not form an edge (line %ld)", ep1, ep2, lineno);
  }
  else if (b->hybrid)
  {
    if (a->parent == b)
    {
      c = a; p = b;
    }
    else if (b->parent == a)
    {
      c = b; p = a;
    }
    else if (b->hybrid->parent == a)
    {
      c = b->hybrid; p = a;
    }
    else
      fatal("Nodes %s and %s do not form an edge (line %ld)", ep1, ep2, lineno);
  }
  else
  {
    if (a->parent == b)
    {
      c = a; p = b;
    }
    else if (b->parent == a)
    {
      c = b; p = a;
    }
    else
      fatal("Nodes %s and %s do not form an edge (line %ld)", ep1, ep2, lineno);
  }

  assert(p && c);
  assert(!p->hybrid || !node_is_mirror(p));
  
  /* check that the specified edge is not horizontal (introgression) */
  if (c->hybrid && !c->htau)
  {
    fatal("Cannot insert hybridization on a horizontal introgression edge "
          "defined by %s and %s (line %ld)",
          ep1,ep2,lineno); 
  }

  return c;
}

static void process_bidir(stree_t * stree, mscidefs_t * def)
{
  long i;

  destroy_pptable(stree);

  opt_msci = 1;

  snode_t * a = edge_basenode(stree, def->node1_1, def->node1_2, def->lineno);
  snode_t * b = edge_basenode(stree, def->node2_1, def->node2_2, def->lineno);

  if (a == b)
    fatal("Cannot create bidirection on a single edge %s - %s (line %ld)",
          def->node1_1, def->node1_2, def->lineno);

  snode_t * pa = a->parent;
  snode_t * pb = b->parent;

  assert(pa && pb);
  
  /* create four new hybridization nodes (2 inner and 2 mirror) */
  snode_t * sh = (snode_t *)xcalloc(1,sizeof(snode_t));
  snode_t * sm = (snode_t *)xcalloc(1,sizeof(snode_t));
  snode_t * th = (snode_t *)xcalloc(1,sizeof(snode_t));
  snode_t * tm = (snode_t *)xcalloc(1,sizeof(snode_t));

  sh->hybrid = sm; sm->hybrid = sh;
  th->hybrid = tm; tm->hybrid = th;

  sh->mark = (int *)xcalloc(1,sizeof(int));
  sm->mark = (int *)xcalloc(1,sizeof(int));
  th->mark = (int *)xcalloc(1,sizeof(int));
  tm->mark = (int *)xcalloc(1,sizeof(int));

  sh->label = xstrdup(def->label1);
  sm->label = xstrdup(def->label1);
  th->label = xstrdup(def->label2);
  tm->label = xstrdup(def->label2);

  /* relink all edges adjacent to t */
  th->parent = pb;
  if (pb->left == b)
  {
    pb->left = th;
  }
  else
  {
    assert(pb->right == b);
    pb->right = th;
  }
  th->left = b;
  th->right = sm;
  b->parent = th;
  sm->parent = th;

  /* relink all edges adjacent to s */
  sh->parent = pa;
  if (pa->left == a)
  {
    pa->left = sh;
  }
  else
  {
    assert(pa->right == a);
    pa->right = sh;
  }
  sh->left   = a;
  sh->right  = tm;
  a->parent  = sh;
  tm->parent = sh;

  assert(def->phi1 >= 0 && def->phi1 <= 1);
  assert(def->phi2 >= 0 && def->phi2 <= 1);

  /* set phi */
  sm->hphi = def->phi1;
  sh->hphi = 1 - def->phi1;

  tm->hphi = def->phi2;
  th->hphi = 1 - def->phi2;

  /* set parent has taus. We will check correctness later when constructing the
     pptable */
  sh->htau = 1;
  

  /* reconstruct nodes array */
  long tip_count = stree->tip_count;
  long inner_count = stree->inner_count+2;
  long hybrid_count = stree->hybrid_count+2;
  long total_nodes = tip_count+inner_count+hybrid_count;
  snode_t ** nodes = (snode_t **)xmalloc((size_t)total_nodes*sizeof(snode_t *));

  /* first copy tips and inner nodes */
  memcpy(nodes,
         stree->nodes,
         (stree->tip_count+stree->inner_count)*sizeof(snode_t *));
  nodes[stree->tip_count+stree->inner_count] = sh;
  nodes[stree->tip_count+stree->inner_count+1] = th;

  /* now copy hybrid nodes */
  memcpy(nodes+tip_count+inner_count,
         stree->nodes+stree->tip_count+stree->inner_count,
         stree->hybrid_count*sizeof(snode_t *));
  nodes[tip_count+inner_count+stree->hybrid_count] = sm;
  nodes[tip_count+inner_count+stree->hybrid_count+1] = tm;

  free(stree->nodes);
  stree->nodes = nodes;
  stree->tip_count = tip_count;
  stree->inner_count = inner_count;
  stree->hybrid_count = hybrid_count;

  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->node_index = i;

  assert(!stree->pptable);
  stree_init_pptable(stree);

  unsigned int sh_index = sh->node_index;
  unsigned int th_index = th->node_index;

  /* node s cannot be ancestral to node t */
  if (stree->pptable[th_index][sh_index] == 1)
    fatal("Node %s cannot be ancestral to node %s in bidirections (line %ld) ",
          sh->label, th->label, def->lineno);
  /* node t cannot be ancestral to node s */
  if (stree->pptable[sh_index][th_index] == 1)
    fatal("Node %s cannot be ancestral to node %s in bidirections (line %ld)",
          th->label, sh->label, def->lineno);
}

static void process_hybrid(stree_t * stree, mscidefs_t * def)
{
  int parallel = 0;
  int proot = 0;    /* special case where parent is root */
  long i;
  snode_t * a = NULL;
  snode_t * b = NULL;
  snode_t * pa = NULL;
  snode_t * pb = NULL;

  destroy_pptable(stree);

  opt_msci = 1;

  a = edge_basenode(stree, def->node1_1, def->node1_2, def->lineno);
  if (def->node2_1)
  {
    b = edge_basenode(stree, def->node2_1, def->node2_2, def->lineno);
  }

  #if 0
  if (a == b)
    parallel = 1;
  #endif
  if (!b)
    parallel = 1;

  pa = a->parent;
  if (!parallel)
    pb = b->parent;

  if (!pa || (!parallel && !pb))
    fatal("Cannot create hybridization event on root node (line %ld)",
          def->lineno);

  /* create two new hybridization nodes */
  snode_t * hl = (snode_t *)xcalloc(1,sizeof(snode_t));
  snode_t * hr = (snode_t *)xcalloc(1,sizeof(snode_t));
  hl->hybrid = hr;
  hr->hybrid = hl;

  /* create the parent of hr */
  snode_t * t = NULL;
  if (!parallel || pa->parent || (pa->left && pa->right))
  {
    t = (snode_t *)xcalloc(1,sizeof(snode_t));
    proot = 1;
  }

  #if 0
  if (parallel)
  {
    if (pa->left && pa->right)
      fatal("Parent of node %s cannot be bifurcating (line %ld)",
            def->node1, def->lineno);
  }
  #endif

  /* create mark array */
  hl->mark = (int *)xcalloc(1,sizeof(int));
  hr->mark = (int *)xcalloc(1,sizeof(int));
  if (proot)
    t->mark  = (int *)xcalloc(1,sizeof(int));

  hl->label = xstrdup(def->label1);
  hr->label = xstrdup(def->label1);
  if (proot)
    t->label = xstrdup(def->label2);

  if (!parallel)
  {
    /* relink all edges adjacent to t */

    t->parent = pb;
    if (pb->left == b)
    {
      pb->left = t;
    }
    else
    {
      assert(pb->right == b);
      pb->right = t;
    }
    t->left = hr;
    t->right = b;
    hr->parent = b->parent = t;

    /* relink all links adjacent to hl */
    if (pa->left == a)
    {
      pa->left = hl;
    }
    else
    {
      assert(pa->right == a);
      pa->right = hl;
    }
    hl->parent = pa;
    hl->left = a;
    a->parent = hl;
    hl->right = NULL;
  }
  else if (parallel)
  {
    if (!pa->parent && (!pa->left || !pa->right))
    {
      /* relink all edges adjacent to A */

      hl->parent = pa;
      hr->parent = pa;
      pa->left = hl;
      pa->right = hr;
      hl->left = a;
      a->parent = hl;
    }
    else
    {
      /* relink all edges adjacent to T */
      assert(t);

      t->parent = pa;
      if (pa->left == a)
      {
        pa->left = t;
      }
      else
      {
        assert(pa->right == a);
        pa->right = t;
      }
      t->left = hl;
      t->right = hr;
      hl->parent = t;
      hr->parent = t;
      hl->left = a;
      a->parent = hl;
    }
  }

  /* set phi */
  if (def->phi1 != -1)
  {
    hl->hphi = def->phi1;
    hr->hphi = 1 - def->phi1;
  }
  else
  {
    hl->hphi = hr->hphi = -1;
  }

  /* set parent has taus. We will check correctness later when constructing the
     pptable */
  if (parallel)
  {
    hl->htau = 1;
    hr->htau = 1;
  }
  else
  {
    hl->htau = def->has_tau1;
    hr->htau = def->has_tau2;
  }

  long new_inners = 0;
  if (!parallel || t)
    new_inners = 2;
  else
    new_inners = 1;

  /* reconstruct nodes array */
  long tip_count = stree->tip_count;
  long inner_count = stree->inner_count+new_inners;
  long hybrid_count = stree->hybrid_count+1;
  long total_nodes = tip_count+inner_count+hybrid_count;
  snode_t ** nodes = (snode_t **)xmalloc((size_t)total_nodes*sizeof(snode_t *));

  /* first copy tips and inner nodes */
  memcpy(nodes,
         stree->nodes,
         (stree->tip_count+stree->inner_count)*sizeof(snode_t *));
  nodes[stree->tip_count+stree->inner_count] = hl;
  if (new_inners == 2)
    nodes[stree->tip_count+stree->inner_count+1] = t;

  /* now copy hybrid nodes */
  memcpy(nodes+tip_count+inner_count,
         stree->nodes+stree->tip_count+stree->inner_count,
         stree->hybrid_count*sizeof(snode_t *));
  nodes[tip_count+inner_count+stree->hybrid_count] = hr;

  free(stree->nodes);
  stree->nodes = nodes;
  stree->tip_count = tip_count;
  stree->inner_count = inner_count;
  stree->hybrid_count = hybrid_count;

  for (i = 0; i < total_nodes; ++i)
    stree->nodes[i]->node_index = i;

  assert(!stree->pptable);
  stree_init_pptable(stree);

  /* if a has a no tau then it must not be parental to the other */

  if (!parallel)
  {
    unsigned int pa_index = pa->node_index;
    unsigned int t_index = t->node_index;;

    if (def->has_tau1 == 0)
    {
      /* node pa cannot be ancestral to node t */
      if (stree->pptable[t_index][pa_index] == 1)
        fatal("Node %s cannot have tau=no and be ancestral to node %s (line %ld)",
              pa->label, t->label, def->lineno);
    }
    if (def->has_tau2 == 0)
    {
      /* node t cannot be ancestral to node pa */
      if (stree->pptable[pa_index][t_index] == 1)
        fatal("Node %s cannot have tau=no and be ancestral to node %s (line %ld)",
              t->label, pa->label, def->lineno);
    }
  }
}


/*** Ziheng 2020-10-2
*  if (opt_msci_faketree_binarize), each hybrid mirror node is converted to a ghost tip.
***/
extern long opt_msci_faketree_binarize;

static char* msci_export_newick_recursive(const snode_t* root, char* (*cb_serialize)(const snode_t*))
{
  char* newick;
  int size_alloced;
  assert(root != NULL);

  if (!(root->left) && !(root->right))
  {
    /* tip or hybrid mirror node */
    if (root->hybrid)
    {
      /* hybrid mirror node */
      assert(node_is_mirror(root));
      if (cb_serialize)
      {
        if (opt_msci_faketree_binarize)
          size_alloced = xasprintf(&newick, "%s :%f", root->label, root->parent->tau);
        else
        {
          newick = cb_serialize(root);
          size_alloced = strlen(newick);
        }
      }
      else
      {
        if (node_is_hybridization(root))
        {
          /* hybridization event */
          if (root->hphi >= 0)
            size_alloced = xasprintf(&newick, "%s[&phi=%f,tau-parent=%s]", root->label, root->hphi, root->htau ? "yes" : "no");
          else
            size_alloced = xasprintf(&newick, "%s[tau-parent=%s]", root->label, root->htau ? "yes" : "no");
        }
        else
        {
          /* bidirectional introgression */
          if (root->hphi >= 0)
            size_alloced = xasprintf(&newick, "%s[&phi=%f]", root->label, root->hphi);
          else
            size_alloced = xasprintf(&newick, "%s", root->label);
        }
      }
    }
    else
    {
      /* tip node */
      if (cb_serialize)
      {
        newick = cb_serialize(root);
        size_alloced = strlen(newick);
      }
      else
      {
        if (show_branches)
          size_alloced = xasprintf(&newick, "%s:%f", root->label, root->length);
        else
          size_alloced = xasprintf(&newick, "%s", root->label);
      }
    }
  }
  else
  {
    /* inner or hybrid node */
    if (root->hybrid)
    {
      /* hybrid non-mirror node */
      assert(!node_is_mirror(root));
      char* subtree = NULL;
      if (node_is_hybridization(root))
      {
        /* hybridization event */
        assert(root->left && !root->right);
        subtree = msci_export_newick_recursive(root->left, cb_serialize);
        if (!subtree)
          return NULL;

        if (cb_serialize)
        {
          char* temp = cb_serialize(root);
          if (!opt_msci_faketree_binarize)
            size_alloced = xasprintf(&newick, "(%s)%s", subtree, temp);
          else
            size_alloced = xasprintf(&newick, "(%s, ghost_%s :%f)%s", subtree, root->label, root->tau, temp);
          free(temp);
        }
        else
        {
          if (root->hphi >= 0)
            size_alloced = xasprintf(&newick, "(%s)%s[&phi=%f,tau-parent=%s]", subtree, root->label ? root->label : "", root->hphi, root->htau ? "yes" : "no");
          else
            size_alloced = xasprintf(&newick, "(%s)%s[tau-parent=%s]", subtree, root->label ? root->label : "", root->htau ? "yes" : "no");
        }
      }
      else
      {
        /* bidirectional introgression */
        assert(root->left && root->right);
        assert(root->right->hybrid && !root->right->left && !root->right->right);
        subtree = msci_export_newick_recursive(root->left, cb_serialize);
        if (!subtree)
          return NULL;

        if (cb_serialize)
        {
          /* assert(0); */
          char* temp = cb_serialize(root);
          if (!opt_msci_faketree_binarize)
            size_alloced = xasprintf(&newick, "(%s, %s)%s", subtree, root->right->label, temp);
          else
            size_alloced = xasprintf(&newick, "(%s, %s_ghost : %f)%s", subtree, root->label, root->tau, temp);
          free(temp);
        }
        else
        {
          if (root->hphi >= 0)
          {
            if (show_branches)
              size_alloced = xasprintf(&newick, "(%s,%s)%s[&phi=%f]:%f", subtree, root->right->label, root->label, root->hphi, root->length);
            else
              size_alloced = xasprintf(&newick, "(%s,%s)%s[&phi=%f]", subtree, root->right->label, root->label, root->hphi);
          }
          else
          {
            if (show_branches)
              size_alloced = xasprintf(&newick, "(%s,%s)%s:%f", subtree, root->right->label, root->label, root->length);
            else
              size_alloced = xasprintf(&newick, "(%s,%s)%s", subtree, root->right->label, root->label);
          }
        }
      }
      free(subtree);
      if (size_alloced < 0)
        fatal("Memory allocation during newick export failed");
    }
    else
    {
      /* inner node */
      char* subtree1 = msci_export_newick_recursive(root->left, cb_serialize);
      if (!subtree1)
        return NULL;
      char* subtree2 = msci_export_newick_recursive(root->right, cb_serialize);
      if (!subtree2)
        return NULL;

      if (cb_serialize)
      {
        char* temp = cb_serialize(root);
        size_alloced = xasprintf(&newick, "(%s, %s)%s", subtree1, subtree2, temp);
        free(temp);
      }
      else
      {
        if (show_branches)
          size_alloced = xasprintf(&newick, "(%s,%s)%s:%f", subtree1, subtree2, root->label ? root->label : "", root->length);
        else
          size_alloced = xasprintf(&newick, "(%s,%s)%s", subtree1, subtree2, root->label ? root->label : "");
      }
      free(subtree1);
      free(subtree2);
      if (size_alloced < 0)
        fatal("Memory allocation during newick export failed");
    }
  }
  return newick;
}



char * msci_export_newick(const snode_t * root, char * (*cb_serialize)(const snode_t *))
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
      if (show_branches)
        size_alloced = xasprintf(&newick, "%s:%f", root->label ? root->label : "", root->length);
      else
        size_alloced = xasprintf(&newick, "%s", root->label ? root->label : "");
    }
  }
  else
  {
    char* subtree1, * subtree2;
    subtree1 = msci_export_newick_recursive(root->left, cb_serialize);
    if (!subtree1)
      fatal("Unable to allocate enough memory.");
    subtree2 = msci_export_newick_recursive(root->right, cb_serialize);
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
      if (show_branches)
        size_alloced = xasprintf(&newick, "(%s, %s)%s:%f;", subtree1, subtree2,  root->label ? root->label : "", root->length);
      else
        size_alloced = xasprintf(&newick, "(%s, %s)%s;", subtree1, subtree2, root->label ? root->label : "");
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("memory allocation during newick export failed");

  return newick;
}

void cmd_msci_create()
{
  long i;
  list_t * list;
  list_item_t * li;
  stree_t * stree = NULL;

  /* 1. Parse definitions file */
  printf("Parsing definititions file %s ...\n", opt_mscifile);
  list = parse_mscifile(opt_mscifile);
  if (!list)
  {
    fatal("Invalid syntax found on file %s", opt_mscifile);
  }

  /* NOTE: list was already checked in parse_mscifile for duplicate tree
     entries */

  /* 2. Reconstruct the species tree */
  printf("Reconstructing species tree ...\n");
  li = list->head;
  while (li)
  {
    mscidefs_t * def = (mscidefs_t *)(li->data);

    if (def->type == BPP_MSCIDEFS_TREE)
    {
      assert(!stree);

      stree = bpp_parse_newick_string(def->node1_1);
      if (!stree)
        fatal("Could not parse species tree");

      /* create marks */
      for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
        stree->nodes[i]->mark = (int *)xcalloc(1,sizeof(int));
    }
    li = li->next;
  }

  /* 3. Label inner nodes according to DEFINE statements */
  label_inner_nodes(list,stree);

  /* 4. Go through definitions */
  printf("Processing hybridization/introgression events ...\n");
  li = list->head;
  while (li)
  {
    mscidefs_t * def = (mscidefs_t *)(li->data);

    if (def->type == BPP_MSCIDEFS_HYBRID)
    {
      process_hybrid(stree,def);
    }
    else if (def->type == BPP_MSCIDEFS_BIDIR)
    {
      process_bidir(stree,def);
    }
    li = li->next;
  }

  /* X. Deallocate */
  list_clear(list,mscidefs_dealloc);
  free(list);

  char * newick = msci_export_newick(stree->root, NULL);
  printf("Newick tree:\n%s\n", newick);
  free(newick);

  stree_destroy(stree,NULL);
}
