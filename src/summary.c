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

static long ulong_bits;      /* number of bits in unsigned long */
static long bitmask_bits;    /* number of bits in bitmask */
static long bitmask_elms;    /* how many unsigned long a bitmask is made of */

static hashtable_t * ht_trivial = NULL;  /* tip node bitmasks */
static hashtable_t * ht_biparts = NULL;  /* bipartitions */

struct bptrivial_s
{
  unsigned long * bitmask;
  char * label;
};

struct bipartition_s
{
  unsigned long * bitmask;
  long count;
};

/* serialize hashtable of bipartitions to an array of bipartitions */
static struct bipartition_s ** hashtable_serialize1p(hashtable_t * ht)
{
  unsigned long i,k;
  struct bipartition_s ** blist;

  blist = (struct bipartition_s **)xmalloc((ht->entries_count+1) *
                                           sizeof(struct bipartition_s *));

  for (i = 0, k = 0; i < ht->table_size; ++i)
  {
    list_t * list = ht->entries[i];

    list_item_t * head = list->head;
    while (head)
    {
      ht_item_t * hi = (ht_item_t *)(head->data);
      blist[k++] = (struct bipartition_s *)(hi->value);
      head = head->next;
    }
  }

  assert(k == ht->entries_count);

  return blist;
}

static unsigned long hash_fnv_long(unsigned long * l, long count)
{
  unsigned long i;
  char * s = (char *)l;
  unsigned long hash = 14695981039346656037UL;
  unsigned long c;

  for (i = 0; i < (size_t)count * sizeof(unsigned long); ++i)
  {
    c = (unsigned long)*s++;
    hash ^= c;
    hash *= 1099511628211UL;
  }

  return hash;
}

static int cb_cmp_bitmask(void * a, void * b)
{
  long i;
  struct bipartition_s * bp = (struct bipartition_s *)a;
  unsigned long * bitmask = (unsigned long *)b;

  for (i = 0; i < bitmask_elms; ++i)
    if (bp->bitmask[i] != bitmask[i])
      return 0;

  return 1;
}
static int cb_cmp_trivial(void * a, void * b)
{
  struct bptrivial_s * trivial = (struct bptrivial_s *)a;
  char * label = (char *)b;

  return (!strcmp(trivial->label,label));
}

/* returns an unary code with the i-th position set */
static unsigned long * bitencode(long i, long count)
{
  assert(count > 0);

  //long alloc_elms = (count / ulong_bits) + ((count % ulong_bits) ? 1 : 0);

  unsigned long * bitmask = (unsigned long *)xcalloc((size_t)bitmask_elms,
                                                     sizeof(unsigned long));
  //long * bitmask = (long *)xcalloc(alloc_elms,sizeof(long));

  long index = i / ulong_bits;

  bitmask[index] = 1ul << (i - index*ulong_bits);

  return bitmask;
}

static void hash_trivial_bp(char ** species, long count)
{
  long i;
  struct bptrivial_s * trivial;

  assert(ht_trivial);
  assert(count > 0);

  for (i = 0; i < count; ++i)
  {
    trivial = (struct bptrivial_s *)xmalloc(sizeof(struct bptrivial_s));
    trivial->label = xstrdup(species[i]);
    trivial->bitmask = bitencode(i,count);

    if (!hashtable_insert(ht_trivial,
                          (void *)(trivial),
                          hash_fnv(trivial->label),
                          hashtable_strcmp))
    {
      /* this should never happen because duplicate taxa were
         already checked for during tree parsing */
      fatal("Duplicate label (%s)", trivial->label);
    }
  }
}

void bipartitions_init(char ** species, long species_count)
{
  assert(species_count > 0);

  /* compute the size (in bits) of an 'unsigned long' and how many of them we
     need to store a bipartition bitmask (bitmask_elms) */
  ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  bitmask_elms = (species_count / ulong_bits) + !!(species_count % ulong_bits);
  bitmask_bits = species_count;

  assert(bitmask_elms > 0);
  assert(bitmask_bits > 0);
  assert(ulong_bits > 0);

  /* create a hash table of trivial bipartitions; i.e. one bitmask for each
     species label with a single set bit indicating the position of the label
     in the list */
  ht_trivial = hashtable_create((size_t)species_count);
  hash_trivial_bp(species,species_count);

  /* now initialize a hash table where we will be storing unique non-trivial
     bipartitions as bitmasks along with their frequencies */
  ht_biparts = hashtable_create(100*(size_t)species_count);

  /* TODO - think of a better size to initialize hashtable */
}

static void bitmask_update_recursive(snode_t * node)
{
  long i;
  if (!node->left) return;

  bitmask_update_recursive(node->left);
  bitmask_update_recursive(node->right);

  node->bitmask = (unsigned long *)xmalloc((size_t)bitmask_elms *
                                           sizeof(unsigned long));

  for (i = 0; i < bitmask_elms; ++i)
    node->bitmask[i] = node->left->bitmask[i] | node->right->bitmask[i];
}

/* prints a bitmask as a sequence of zeroes and ones */
static void bitmask_print(FILE * fp_out, unsigned long * bitmask)
{
  long i,j;
  long bits_left = bitmask_bits;
  long bits_per_elm = ulong_bits;

  for (i = 0; i < bitmask_elms; ++i)
  {
    unsigned long bits = bitmask[i];
    long bits_avail = MIN(bits_left,bits_per_elm);
    for (j = 0; j < bits_avail; ++j)
    {
      fprintf(stdout, "%c", (char)((bits & 1) + 0x30));
      fprintf(fp_out, "%c", (char)((bits & 1) + 0x30));
      bits >>= 1ul;
    }
    bits_left -= bits_per_elm; 
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");
}

/* assign trivial bipartition bitmasks from hashtable to the tips of stree and
then recursively compute bitmasks for all nodes */
static void assign_bitmasks(stree_t * stree)
{
  long i;
  struct bptrivial_s * trivial;

  /* attach bitmasks at tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    trivial = hashtable_find(ht_trivial,
                             (void *)(stree->nodes[i]->label),
                             hash_fnv(stree->nodes[i]->label),
                             cb_cmp_trivial);
    if (!trivial)
      fatal("Internal error in locating tip label %s (assign_bitmasks())",
             stree->nodes[i]->label);

    stree->nodes[i]->bitmask = trivial->bitmask;
  }

  /* now recursively create bitmasks */

  bitmask_update_recursive(stree->root);
}

/* updates counts in hashtable with bipartitions of current tree */
void bipartitions_update(stree_t * stree)
{
  long i;
  struct bipartition_s * bp;

  assign_bitmasks(stree);

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->parent)
    {
      bp = hashtable_find(ht_biparts,
                          (void *)(stree->nodes[i]->bitmask),
                          hash_fnv_long(stree->nodes[i]->bitmask,bitmask_elms),
                          cb_cmp_bitmask);
      if (bp)
      {
        bp->count++;
      }
      else
      {
        bp = (struct bipartition_s *)xmalloc(sizeof(struct bipartition_s));
        bp->bitmask = (unsigned long *)xmalloc((size_t)bitmask_elms *
                                               sizeof(unsigned long));
        bp->count = 1;
        memcpy(bp->bitmask,
               stree->nodes[i]->bitmask,
               (size_t)bitmask_elms*sizeof(unsigned long));
        hashtable_insert_force(ht_biparts,
                               (void *)bp,
                               hash_fnv_long(stree->nodes[i]->bitmask,
                                             bitmask_elms));
      }
    }
  }

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->bitmask)
      free(stree->nodes[i]->bitmask);
}

static void cb_bptrivial_dealloc(void * data)
{
  struct bptrivial_s * trivial = data;
  free(trivial->bitmask);
  free(trivial->label);
  free(trivial);
}

static void cb_bipartition_dealloc(void * data)
{
  struct bipartition_s * bp = data;
  free(bp->bitmask);
  free(bp);
}


static char * cb_serialize_support(const snode_t * node)
{
  char * s = NULL;

  /* inner node */
  if (node->left)
  {
    if (node->parent)
      xasprintf(&s, " #%f", node->support);
    else
      s = xstrdup("");
  }
  else
    xasprintf(&s, "%s", node->label);

  return s;
}

static void print_stree_with_support(FILE * fp_out,
                                     const char * treestr,
                                     size_t freq,
                                     size_t trees_count)
{
  long i;
  struct bipartition_s * bp;

  stree_t * stree = bpp_parse_newick_string(treestr);
  if (!stree)
    fatal("Internal error while parsing species tree");

  /* assign bitmasks to tree nodes (present bits indicate species in subtree) */
  assign_bitmasks(stree);

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->parent)
    {
      bp = hashtable_find(ht_biparts,
                          (void *)(stree->nodes[i]->bitmask),
                          hash_fnv_long(stree->nodes[i]->bitmask,bitmask_elms),
                          cb_cmp_bitmask);
      if (!bp)
        fatal("Internal error when printing best tree with support values");

      stree->nodes[i]->support = bp->count / (double)trees_count;
    }
  }

  char * newick = stree_export_newick(stree->root, cb_serialize_support);
  fprintf(stdout, "%s   [P = %f]\n", newick, freq / (double)trees_count);
  fprintf(fp_out, "%s   [P = %f]\n", newick, freq / (double)trees_count);
  free(newick);

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->bitmask)
      free(stree->nodes[i]->bitmask);

  stree_destroy(stree,NULL);
}

void summary_dealloc_hashtables()
{
  hashtable_destroy(ht_biparts,cb_bipartition_dealloc);
  hashtable_destroy(ht_trivial,cb_bptrivial_dealloc);
}

static int cb_countcmp(const void * a, const void * b)
{
  const struct bipartition_s * pa = *((const struct bipartition_s **)a);
  const struct bipartition_s * pb = *((const struct bipartition_s **)b);

  if (pa->count < pb->count) return 1;
  else if (pa->count > pb->count) return -1;

  return 0;

}

static int cb_popcntcmp(const void * a, const void * b)
{
  long i;
  long bitsa = 0,bitsb = 0;
  const struct bipartition_s * pa = *((const struct bipartition_s **)a);
  const struct bipartition_s * pb = *((const struct bipartition_s **)b);

  for (i = 0; i < bitmask_elms; ++i)
  {
    bitsa += PLL_POPCOUNTL(pa->bitmask[i]);
    bitsb += PLL_POPCOUNTL(pb->bitmask[i]);
  }

  if (bitsa > bitsb) return 1;
  else if (bitsa < bitsb) return -1;

  return 0;

}

static void bipartitions_finalize(FILE * fp_out, size_t trees_count, char ** species)
{
  unsigned long i,j,k,entry;
  unsigned long majority = 0;

  struct bipartition_s ** x = hashtable_serialize1p(ht_biparts);
  fprintf(stdout, "\n(B) Best splits in the sample of trees (%ld splits in all)\n",
          ht_biparts->entries_count);
  fprintf(fp_out, "\n(B) Best splits in the sample of trees (%ld splits in all)\n",
          ht_biparts->entries_count);

  /* find number of bipartitions appearing in at least 50% of the samples */
  qsort(x,ht_biparts->entries_count,sizeof(struct bipartition_s *), cb_countcmp);
  for (i = 0; i < ht_biparts->entries_count; ++i)
  {
      fprintf(stdout, "%6ld %f  ", x[i]->count, x[i]->count / (double)trees_count);
      fprintf(fp_out, "%6ld %f  ", x[i]->count, x[i]->count / (double)trees_count);
      bitmask_print(fp_out,x[i]->bitmask);
      if (x[i]->count / (double)trees_count >= 0.5)
        majority++;
  }


  /* now sort them according to pop count */

  /* we store an array of pointers to pointers of character strings. If sptr[i] == NULL then the name is taken from
  species[i]. If sptr[i] != NULL then we check if *(sptr[i]) != NULL and if yes, we use that string in output, otherwise
  if *(sptr[i]) == NULL it means that the string was already used in this round from another selected taxon */
  char *** sptr = (char ***)xcalloc((size_t)bitmask_bits,sizeof(char **));
  char ** s = (char **)xmalloc((size_t)bitmask_bits * sizeof(char *));
  unsigned long * index = (unsigned long *)xmalloc((size_t)bitmask_bits *
                                                   sizeof(unsigned long));
  unsigned long * index_all = (unsigned long *)xmalloc((size_t)bitmask_bits *
                                                       sizeof(unsigned long));
  long ptrcount;
  long allcount;


  /* sort bipartitions by number of set bits in bitmask */
  qsort(x, majority, sizeof(struct bipartition_s *), cb_popcntcmp);

  /* add one more bitmask with all bits set */

  x[majority] = (struct bipartition_s *)xmalloc(sizeof(struct bipartition_s));
  struct bipartition_s * debug_free = x[majority];
  x[majority]->count = 0;
  x[majority]->bitmask = (unsigned long *)xmalloc((size_t)bitmask_elms *
                                                  sizeof(unsigned long));
  memset(x[majority]->bitmask,0xff,(size_t)bitmask_elms*sizeof(unsigned long));
  majority++;


  /* now construct the newick string for the majority rule consensus tree */
  char ** newick = (char **)xcalloc(majority,sizeof(char *));
  for (i = 0; i < majority; ++i)
  {
    ptrcount = 0;
    entry = 0;
    allcount = 0;

    for (j = 0; j < (unsigned long)bitmask_elms; ++j)
    {
      for (k = 0; k < (unsigned long)ulong_bits; ++k)
      {
        if (entry == (unsigned long)bitmask_bits) break;

        if ((x[i]->bitmask[j] >> k) & 1)
        {
          if (!sptr[entry])
          {
            s[ptrcount] = xstrdup(species[entry]);
            index[ptrcount++] = entry;
          }
          else if (*(sptr[entry]))
          {
            s[ptrcount] = *(sptr[entry]);
            index[ptrcount++] = entry;
            *(sptr[entry]) = NULL;
          }
          else
          {
            /* already used in this iteration, just skip */

          }
          index_all[allcount++] = entry;
        }
        entry++;
      }
    }
    
    /* assemble */
    char * temp;
    xasprintf(newick+i,"(%s",s[0]);
    free(s[0]);
    for (j = 1; j < (unsigned long)ptrcount; ++j)
    {
      xasprintf(&temp, "%s, %s", newick[i], s[j]);
      free(newick[i]);
      newick[i] = temp;
      free(s[j]);
    }

    if (i == majority - 1)
      xasprintf(&temp, "%s);", newick[i]);
    else
      xasprintf(&temp, "%s) #%f", newick[i], x[i]->count / (double)trees_count);

    free(newick[i]);
    newick[i] = temp;

    for (j = 0; j < (unsigned long)allcount; ++j)
    {
      sptr[index_all[j]] = newick+i;
    }
  }
  fprintf(stdout, "\n(C) Majority-rule consensus tree\n");
  fprintf(fp_out, "\n(C) Majority-rule consensus tree\n");
  fprintf(stdout, "%s\n", newick[majority-1]);
  fprintf(fp_out, "%s\n", newick[majority-1]);
  for (i = 0; i < majority; ++i)
    free(newick[i]);
  free(newick);

  /* cleanup */
  free(debug_free->bitmask);
  free(debug_free);
  free(x);
  free(s);
  free(index);
  free(index_all);

  free(sptr);
}



/* TODO: New summary code */

struct distinct_s 
{
  size_t start;
  size_t count;
};

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static char * cb_serialize_none(const snode_t * snode)
{
  if (!snode->left) return xstrdup(snode->label);

  return xstrdup("");
}

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

static void strip_attributes(char * s)
{
  char * p = s;

  while (*s)
  {
    if (*s == '#' || *s == ':')
    {
      while (*s && *s != ',' && *s != ')' && *s != ';')
        ++s;
    }
    else if (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n')
      ++s;
    else
      *p++ = *s++;
  }
  *p = 0;
}

static int cb_strcmp(const void * a, const void * b)
{
  const char ** pa = (const char **)a;
  const char ** pb = (const char **)b;

  return strcmp(*pa,*pb);
}

static void stree_sort_recursive(snode_t * node)
{
  if (!node->left)
    return;

  stree_sort_recursive(node->left);
  stree_sort_recursive(node->right);

  if (strcmp(node->left->label, node->right->label) > 0)
    SWAP(node->left,node->right);

  node->label = (char *)xmalloc(strlen(node->left->label) +
                                strlen(node->right->label) + 1);

  /* concatenate species labels */
  node->label[0] = 0;
  strcat(node->label,node->left->label);
  strcat(node->label,node->right->label);
}

static void stree_sort(stree_t * stree)
{
  stree_sort_recursive(stree->root); 
}

static int cb_dtree_cmp(const void * a, const void * b)
{
  const struct distinct_s * pa = (const struct distinct_s *)a;
  const struct distinct_s * pb = (const struct distinct_s *)b;

  if (pa->count < pb->count) return 1;
  else if (pa->count > pb->count) return -1;

  return 0;
}

void stree_summary(FILE * fp_out, char ** species_names, long species_count)
{
  size_t i,distinct;
  size_t line_count = 0;
  FILE * fp_mcmc;
  char ** treelist;
  struct distinct_s * dtree;

  /* allocate space for holding all species tree samples */
  treelist = (char **)xmalloc((size_t)(opt_samples+1)*sizeof(char *));

  /* open mcmc file */
  #ifndef DEBUG_MAJORITY
  fp_mcmc = xopen(opt_mcmcfile,"r");
  #else
  fp_mcmc = xopen("test.txt","r");
  #endif

  bipartitions_init(species_names,species_count);

  /* read each line from the file, and strip all thetas and branch lengths
     such that only the tree topology and tip names remain, and store them
     in treelist */
  while (getnextline(fp_mcmc))
  {
    strip_attributes(line);
    stree_t * t = bpp_parse_newick_string(line);
    if (!t)
      fatal("Internal error while parsing species tree");
    stree_sort(t);
    treelist[line_count++] = stree_export_newick(t->root,cb_serialize_none);

    bipartitions_update(t);
    stree_destroy(t,NULL);
  }
  assert(line_count);


  fprintf(stdout, "Species in order:\n");
  fprintf(fp_out, "Species in order:\n");
  for (i = 0; i < (size_t)species_count; ++i)
  {
    fprintf(stdout, " %3ld. %s\n", i+1, species_names[i]);
    fprintf(fp_out, " %3ld. %s\n", i+1, species_names[i]);
  }
  fprintf(stdout, "\n");
  fprintf(fp_out, "\n");

  qsort(treelist,(size_t)line_count,sizeof(char *), cb_strcmp);

  distinct = 1;
  size_t * uniquepos = (size_t *)xmalloc(line_count*sizeof(size_t));
  uniquepos[0] = 0;
  for (i = 1; i < line_count; ++i)
  {
    if (strcmp(treelist[i],treelist[i-1]))
    {
      uniquepos[distinct++] = i;
    }
  }

  assert(distinct > 0);

  dtree = (struct distinct_s *)xmalloc(distinct * sizeof(struct distinct_s));
  for (i = 0; i < distinct; ++i)
  {
    dtree[i].start = uniquepos[i];
    dtree[i].count = (i == distinct - 1) ?
                      line_count - uniquepos[i] : uniquepos[i+1] - uniquepos[i];
  }

  qsort(dtree, distinct, sizeof(struct distinct_s), cb_dtree_cmp);
  fprintf(stdout, "(A) Best trees in the sample (%ld distinct trees in all)\n", distinct);
  fprintf(fp_out, "(A) Best trees in the sample (%ld distinct trees in all)\n", distinct);
  double cdf = 0;
  for (i = 0; i < distinct; ++i)
  {
    double pdf = dtree[i].count / (double)line_count;
    cdf += pdf;
    fprintf(stdout, " %8ld %8.5f %8.5f %s\n",
            dtree[i].count, pdf, cdf, treelist[dtree[i].start]);
    fprintf(fp_out, " %8ld %8.5f %8.5f %s\n",
            dtree[i].count, pdf, cdf, treelist[dtree[i].start]);
  }

  bipartitions_finalize(fp_out,line_count,species_names);

  fprintf(stdout, "\n(D) Best tree (or trees from the mastertree file) "
          "with support values\n");
  fprintf(fp_out, "\n(D) Best tree (or trees from the mastertree file) "
          "with support values\n");
  for (i = 0; i < distinct; ++i)
  {
    if (i && dtree[i].count != dtree[i-1].count)
      break;

    print_stree_with_support(fp_out,
                             treelist[dtree[i].start],
                             dtree[i].count,
                             line_count);
  }

  summary_dealloc_hashtables();
  free(uniquepos);
  free(dtree);
 
  /* deallocate list of trees */
  for (i = 0; i < line_count; ++i)
    free(treelist[i]);
  free(treelist);

  fclose(fp_mcmc);
}

long getlinecount(const char * filename)
{
  long linecount = 0;
  FILE * fp;

  fp = xopen(filename,"r");

  /* read number of lines */
  while (getnextline(fp)) linecount++;

  fclose(fp);

  return linecount;
}
