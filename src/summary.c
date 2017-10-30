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

#include "bpp.h"

static int ** splits;
static long splits_size;
static long splits_maxsize;

static long splits_incrementsize;
static long splits_initialsize;

static long bitmask_elms;
static long bits_count;

static hashtable_t * ht_bm;     /* tip node bitmasks */
static hashtable_t * ht_bp;     /* bipartitions */

struct bitmask_pair_s
{
  long * bitmask;
  char * label;
};

struct bipfreq_s
{
  long * bitmask;
  long freq;
};

static struct bipfreq_s ** hashtable_serialize1p(hashtable_t * ht)
{
  unsigned long i,k;

  struct bipfreq_s ** bf = (struct bipfreq_s **)xmalloc((ht->entries_count+1) *
                                                        sizeof(struct bipfreq_s *));
  for (i = 0, k = 0; i < ht->table_size; ++i)
  {
    list_t * list = ht->entries[i];

    list_item_t * head = list->head;
    while (head)
    {
      ht_item_t * hi = (ht_item_t *)(head->data);
      bf[k++] = (struct bipfreq_s *)(hi->value);
      head = head->next;
    }
  }

  assert(k == ht->entries_count);

  return bf;
}

static unsigned long hash_fnv_long(long * l, long count)
{
  unsigned long i;
  char * s = (char *)l;
  unsigned long hash = 14695981039346656037UL;
  unsigned long c;

  for (i = 0; i < count * sizeof(long); ++i)
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
  struct bipfreq_s * bpf = (struct bipfreq_s *)a;
  long * bitmask = (long *)b;

  for (i = 0; i < bitmask_elms; ++i)
    if (bpf->bitmask[i] != bitmask[i])
      return 0;

  return 1;
}
static int cb_cmp_bmpair(void * a, void * b)
{
  struct bitmask_pair_s * bmpair = (struct bitmask_pair_s *)a;
  char * label = (char *)b;

  return (!strcmp(bmpair->label,label));
}

/* returns an unary code with the i-th position set */
static long * bitencode(long i, long count)
{
  assert(count > 0);

  long lsize_bits = sizeof(long) * CHAR_BIT;
  long alloc_elms = (count / lsize_bits) + ((count % lsize_bits) ? 1 : 0);

  long * bitmask = (long *)xmalloc(alloc_elms * sizeof(long));

  long index = count / lsize_bits;

  bitmask[index] = 1 << i;

  return bitmask;
}

static hashtable_t * bitmask_hash(char ** species, long count)
{
  long i;
  struct bitmask_pair_s * bmpair;

  hashtable_t * ht = hashtable_create(count);

  for (i = 0; i < count; ++i)
  {
    bmpair = (struct bitmask_pair_s *)xmalloc(sizeof(struct bitmask_pair_s));
    bmpair->label = xstrdup(species[i]);
    bmpair->bitmask = bitencode(i,count);

    if (!hashtable_insert(ht,
                          (void *)(bmpair),
                          hash_fnv(bmpair->label),
                          hashtable_strcmp))
    {
      /* this should never happen because duplicate taxa were
         already checked for during tree parsing */
      fatal("Duplicate label (%s)", bmpair->label);
    }
  }

  return ht;
}

void splits_init(long init, long increment, char ** species, long species_count)
{
  splits_initialsize = init;
  splits_incrementsize = increment;
  splits_size = 0;

  splits = (int **)xmalloc(init*sizeof(int *));
  splits_maxsize = init;

  long lsize_bits = sizeof(long) * CHAR_BIT;
  bitmask_elms = (species_count / lsize_bits) +
                 ((species_count % lsize_bits) ? 1 : 0);
  bits_count = species_count;

  ht_bm = bitmask_hash(species,species_count);
  //ht_bp = hashtable_create(4194304);

  /* TODO : CREATE HASHTABLE ACCORDING TO NUMBER OF SPECIES */
  ht_bp = hashtable_create(100);
}

static void bitmask_update_recursive(snode_t * node)
{
  long i;
  if (!node->left) return;

  bitmask_update_recursive(node->left);
  bitmask_update_recursive(node->right);

  node->bitmask = (long *)xmalloc(bitmask_elms * sizeof(long));

  for (i = 0; i < bitmask_elms; ++i)
    node->bitmask[i] = node->left->bitmask[i] | node->right->bitmask[i];
}

static void bitmask_print(long * bitmask)
{
  long i,j;
  long bits_left = bits_count;
  long bits_per_elm = sizeof(long)*CHAR_BIT;

  for (i = 0; i < bitmask_elms; ++i)
  {
    long bits = bitmask[i];
    long bits_avail = MIN(bits_left,bits_per_elm);
    for (j = 0; j < bits_avail; ++j)
    {
      printf("%c", (char)((bits & 1) + 0x30));
      bits >>= 1;
    }
    bits_left -= bits_per_elm; 
  }
  printf("\n");
}

static void bitmasks_init(stree_t * stree)
{
  long i;
  struct bitmask_pair_s * bmpair;

  /* attach bitmasks at tip nodes */
  for (i = 0; i < stree->tip_count; ++i)
  {
    bmpair = hashtable_find(ht_bm,
                            (void *)(stree->nodes[i]->label),
                            hash_fnv(stree->nodes[i]->label),
                            cb_cmp_bmpair);
    if (!bmpair)
      fatal("Internal error in splits_init");

    stree->nodes[i]->bitmask = bmpair->bitmask;
  }

  /* now recursively create bitmasks */

  bitmask_update_recursive(stree->root);
}

void splits_update(stree_t * stree)
{
  long i;
  struct bipfreq_s * bpf;

  bitmasks_init(stree);

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->parent)
    {
      bpf = hashtable_find(ht_bp,
                           (void *)(stree->nodes[i]->bitmask),
                           hash_fnv_long(stree->nodes[i]->bitmask,bitmask_elms),
                           cb_cmp_bitmask);
      if (bpf)
      {
        bpf->freq++;
      }
      else
      {
        bpf = (struct bipfreq_s *)xmalloc(sizeof(struct bipfreq_s));
        bpf->bitmask = (long *)xmalloc(bitmask_elms*sizeof(long));
        bpf->freq = 1;
        memcpy(bpf->bitmask,stree->nodes[i]->bitmask,bitmask_elms*sizeof(long));
        hashtable_insert_force(ht_bp,
                               (void *)bpf,
                               hash_fnv_long(stree->nodes[i]->bitmask,
                                             bitmask_elms));
      }

    }
  }

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->bitmask)
      free(stree->nodes[i]->bitmask);
}
static void cb_bm_dealloc(void * data)
{
  struct bitmask_pair_s * bm = data;
  free(bm->bitmask);
  free(bm->label);
  free(bm);
}

static void cb_bf_dealloc(void * data)
{
  struct bipfreq_s * bf = data;
  free(bf->bitmask);
  free(bf);
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

void print_stree_with_support(const char * treestr, long freq, long trees_count)
{
  long i;
  struct bipfreq_s * bpf;

  stree_t * stree = stree_parse_newick_string(treestr);

  /* assign bitmasks to tree nodes (present bits indicate species in subtree) */
  bitmasks_init(stree);

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    if (stree->nodes[i]->parent)
    {
      bpf = hashtable_find(ht_bp,
                           (void *)(stree->nodes[i]->bitmask),
                           hash_fnv_long(stree->nodes[i]->bitmask,bitmask_elms),
                           cb_cmp_bitmask);
      if (!bpf)
        fatal("Internal error when printing best tree with support values");

      stree->nodes[i]->support = bpf->freq / (double)trees_count;
    }
  }

  char * newick = stree_export_newick(stree->root, cb_serialize_support);
  fprintf(stdout, "%s   [P = %f]\n", newick, freq / (double)trees_count);
  free(newick);

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
    if (stree->nodes[i]->bitmask)
      free(stree->nodes[i]->bitmask);

  stree_destroy(stree,NULL);
}

void summary_dealloc_hashtables()
{
  hashtable_destroy(ht_bp,cb_bf_dealloc);
  hashtable_destroy(ht_bm,cb_bm_dealloc);
}

static int cb_freqcmp(const void * a, const void * b)
{
  const struct bipfreq_s * pa = *((const struct bipfreq_s **)a);
  const struct bipfreq_s * pb = *((const struct bipfreq_s **)b);

  if (pa->freq < pb->freq) return 1;
  else if (pa->freq > pb->freq) return -1;

  return 0;

}

static int cb_popcntcmp(const void * a, const void * b)
{
  long i;
  long bitsa = 0,bitsb = 0;
  const struct bipfreq_s * pa = *((const struct bipfreq_s **)a);
  const struct bipfreq_s * pb = *((const struct bipfreq_s **)b);

  for (i = 0; i < bitmask_elms; ++i)
  {
    bitsa += PLL_POPCOUNT(pa->bitmask[i]);
    //bitsa += __builtin_popcountl(pa->bitmask[i]);
    bitsb += PLL_POPCOUNT(pb->bitmask[i]);
    //bitsb += __builtin_popcountl(pb->bitmask[i]);
  }

  if (bitsa > bitsb) return 1;
  else if (bitsa < bitsb) return -1;

  return 0;

}

void splits_finalize(long trees_count, char ** species)
{
  unsigned long i,j,k,entry;
  unsigned long majority = 0;

  struct bipfreq_s ** x = hashtable_serialize1p(ht_bp);
  printf("\n(B) Best splits in the sample of trees (%ld splits in all)\n",
         ht_bp->entries_count);
  qsort(x,ht_bp->entries_count,sizeof(struct bipfreq_s *), cb_freqcmp);
  for (i = 0; i < ht_bp->entries_count; ++i)
  {
      printf("%6ld %f  ", x[i]->freq, x[i]->freq / (double)trees_count);
      bitmask_print(x[i]->bitmask);
      if (x[i]->freq / (double)trees_count >= 0.5)
        majority++;
  }


  /* now sort them according to pop count */

  /* we store an array of pointers to pointers of character strings. If sptr[i] == NULL then the name is taken from
  species[i]. If sptr[i] != NULL then we check if *(sptr[i]) != NULL and if yes, we use that string in output, otherwise
  if *(sptr[i]) == NULL it means that the string was already used in this round from another selected taxon */
  char *** sptr = (char ***)xcalloc(bits_count,sizeof(char **));
  char ** s = (char **)xmalloc(bits_count * sizeof(char *));
  long * index = (long *)xmalloc(bits_count * sizeof(long));
  long * index_all = (long *)xmalloc(bits_count * sizeof(long));
  long ptrcount;
  long allcount;


  qsort(x, majority, sizeof(struct bipfreq_s *), cb_popcntcmp);

  /* add one more entry of all ones set */

  x[majority] = (struct bipfreq_s *)xmalloc(sizeof(struct bipfreq_s));
  struct bipfreq_s * debug_free = x[majority];
  x[majority]->freq = 0;
  x[majority]->bitmask = (long *)xmalloc(bitmask_elms*sizeof(long));
  memset(x[majority]->bitmask,0xff,bitmask_elms*sizeof(long));
  majority++;


  char ** newick = (char **)xcalloc(majority,sizeof(char *));
  for (i = 0; i < majority; ++i)
  {
    ptrcount = 0;
    entry = 0;
    allcount = 0;

    for (j = 0; j < (unsigned long)bitmask_elms; ++j)
    {
      for (k = 0; k < sizeof(long)*CHAR_BIT && entry < (unsigned long)bits_count; ++k)
      {
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
      xasprintf(&temp, "%s) #%f", newick[i], x[i]->freq / (double)trees_count);

    free(newick[i]);
    newick[i] = temp;

    for (j = 0; j < (unsigned long)allcount; ++j)
    {
      sptr[index_all[j]] = newick+i;
    }
  }
  printf("\n(C) Majority-rule consensus tree\n");
  printf("%s\n", newick[majority-1]);
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
  free(splits);
}
