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

static hashtable_t * sht;
static hashtable_t * mht;

static long * ext_qsort_arg;

static int cb_cmp_indices(const void * a, const void * b)
{
  const long * x = (long *)a;
  const long * y = (long *)b;

  return ((ext_qsort_arg[*x] < ext_qsort_arg[*y]) -
          (ext_qsort_arg[*x] > ext_qsort_arg[*y]));
}

static void cb_dealloc_pairlabel(void * data)
{
  pair_t * pair = data;

  free(pair->label);
  free(pair);
}

/* return an array of ones and zeroes indicating whether sequences in alignment
   msa are diploids (1) or haploids (0) */
static unsigned int * get_diploid_info(msa_t * msa, int msa_id)
{
  long i;
  unsigned int * diploid;

  diploid = (unsigned int *)xmalloc((size_t)(msa->count)*sizeof(unsigned int));

  /* if we have only one species */
  if (sht == NULL && mht == NULL)
  {
    assert(opt_diploid_size == 1);

    for (i = 0; i < msa->count; ++i)
      diploid[i] = opt_diploid[0] ? 1 : 0;

    return diploid;
  }

  for (i = 0; i < msa->count; ++i)
  {
    char * label = msa->label[i];
    label = strchr(label, '^');
    if (!label)
      fatal("Cannot find species tag on sequence %s of locus %d",
            msa->label[i], msa_id);

    /* skip the '^' mark */
    label++;
    if (!(*label))
      fatal("Sequence %s of locus %d contains no label",
            msa->label[i], msa_id);

    pair_t * pair;
    pair = hashtable_find(mht,
                          (void *)label,
                          hash_fnv(label),
                          cb_cmp_pairlabel);

    if (!pair)
      fatal("Cannot find species mapping for sequence %s of locus %d",
            label, msa_id);

    snode_t * node = (snode_t *)(pair->data);

    diploid[i] = node->diploid;
  }
  return diploid;
}

#if 0
static list_t * update_map_list()
{
  list_t * newlist = (list_t *)xmalloc(sizeof(list_t));

  for (i = 0; i < ht->table_size; ++i)
  {
    list_t * list = mht->entries[i];

    list_item_t * head = list->head;
    while(head)
    {
      ht_item_t * hi = (ht_item_t *)(head->data);

      pair_t * pair  = (pair_t *)(hi->value);
      char * label   = pair->label;
      snode_t * node = (snode_t *)(pair->data);



      if (node->diploid)
      {
        ht_item_t 

      }
      else
      {
      }


    }
  }

  return newlist;
}
#endif

static int cb_cmp_nodelabel(void * a, void * b)
{
  snode_t * node = (snode_t *)a;
  char * label = (char * )b;

  return (!strcmp(node->label,label));
}

static void update_map_list(list_t * maplist)
{
  mapping_t * newmap;

  list_t * newlist = (list_t *)xcalloc(1,sizeof(list_t));

  list_item_t * li = maplist->head;
  while (li)
  {
    mapping_t * map = (mapping_t *)(li->data);

    snode_t * node = hashtable_find(sht,
                                    (void *)(map->species),
                                    hash_fnv(map->species),
                                    cb_cmp_nodelabel);
    if (!node)
      fatal("cannot find node with population label %s", map->species);

    printf("%s -> %s (diploid=%u)\n", map->individual, map->species, node->diploid);
    if (node->diploid)
    {
      /* add first label mapping */
      newmap = (mapping_t *)xmalloc(sizeof(mapping_t));
      xasprintf(&(newmap->individual), "%s.1", map->individual);
      newmap->species = xstrdup(map->species);
      newmap->lineno = map->lineno;
      list_append(newlist, (void *)newmap);

      /* add second label mapping */
      newmap = (mapping_t *)xmalloc(sizeof(mapping_t));
      xasprintf(&(newmap->individual), "%s.2", map->individual);
      newmap->species = xstrdup(map->species);
      newmap->lineno = map->lineno;
      list_append(newlist, (void *)newmap);
    }
    else
    {
      newmap = (mapping_t *)xmalloc(sizeof(mapping_t));
      newmap->individual = xstrdup(map->individual);
      newmap->species    = xstrdup(map->species);
      newmap->lineno = map->lineno;
      list_append(newlist,(void *)newmap);
    }
    li = li->next;
  }

  printf("\n\n");

  li = newlist->head;
  while(li)
  {
    mapping_t * map = (mapping_t *)(li->data);

    snode_t * node = hashtable_find(sht,
                                    (void *)(map->species),
                                    hash_fnv(map->species),
                                    cb_cmp_nodelabel);
    if (!node)
      fatal("cannot find node with population label %s", map->species);

    printf("%s -> %s (diploid=%u)\n", map->individual, map->species, node->diploid);

    li = li->next;
  }

  list_clear(maplist,map_dealloc);

  memcpy(maplist,newlist,sizeof(list_t));

  free(newlist);


}

static void diploid_resolution_init(stree_t * stree, list_t * maplist)
{
  sht = species_hash(stree);
  mht = maplist_hash(maplist,sht);
}

static void diploid_resolution_fini()
{
  hashtable_destroy(sht,NULL);
  hashtable_destroy(mht,cb_dealloc_pairlabel);
}

static unsigned int findmax(const unsigned int * map)
{
  int i;
  unsigned int max = 0;

  for (i = 0; i < ASCII_SIZE; ++i)
    if (map[i] > max)
      max = map[i];

  return max;
}

static unsigned long * diploid_resolve_locus(msa_t * msa,
                                             int msa_index,
                                             unsigned int * weight,
                                             int * cleandata,
                                             const unsigned int * map)
{
  int card;
  long i,j,k,n,m;
  long round;
  long unresolved_count = 0;
  unsigned int c;
  unsigned char inv_charmap[ASCII_SIZE];
  unsigned char charmap[ASCII_SIZE];
  int * hmat;
  int * hptr;
  long * resolved;  
  long * seqhets;         /* # heterozygotes at sequence */
  long * sitehets;        /* # heterozygotes at site */
  long * singletons;      /* # singleton sites at sequence */
  long * single_indices;  /* indices of singletons with at least one het */
  unsigned int * diploid;
  char ** newlabel;
  unsigned long * resolution_count;

  *cleandata = 1;

  /* allocate necessary arrays */
  resolved = (long *)xmalloc((size_t)(msa->count) * sizeof(long));
  seqhets = (long *)xmalloc((size_t)(msa->count) * sizeof(long));
  singletons = (long *)xcalloc((size_t)(msa->count),sizeof(long));
  sitehets = (long *)xcalloc((size_t)(msa->length),sizeof(long));
  single_indices = (long *)xmalloc((size_t)(msa->length) * sizeof(long));
  hmat = (int *)xcalloc((size_t)(msa->count*msa->length),sizeof(int));
  hptr = hmat;

  diploid = get_diploid_info(msa,msa_index);


  /* if map states are out of the BYTE range, remap */
  if (findmax(map) >= ASCII_SIZE)
  {
    /* for now only DNA with pll_map_nt */
    assert(0);
  }
  else
  {
    for (i = 0; i < ASCII_SIZE; ++i)
      charmap[i] = (unsigned char)(map[i]);
  }

  /* create inverse charmap to decode states back to characters when finished */
  for (i = 0; i < ASCII_SIZE; ++i)
    if (map[i])
      inv_charmap[charmap[i]] = (unsigned char)i;

  for (i = 0; i < msa->count; ++i)
    resolved[i] = 1;

  /* 1. Fill the h matrix and set cleandata for current locus */
  for (i = 0; i < msa->count; ++i)
  {
    /* Ensures that cleandata is set to 0 in case there exists a haploid with at
       least one ambiguous character */
    if (!diploid[i] && *cleandata)
    {
      for (j = 0; j < msa->length; ++j)
      {
        c = map[(int)(msa->sequence[i][j])];
        if (PLL_POPCOUNT(c) > 1)
        {
          *cleandata = 0;
          break;
        }
      }
    }

    /* if diploid, compute entries of current row in matrix h */
    if (diploid[i])
    {
      for (j = 0; j < msa->length; ++j)
      {
        c = map[(int)(msa->sequence[i][j])];
        
        /* number of states of ambiguity */
        card = PLL_POPCOUNT(c);
        
        if (card == 1)
          continue;
        else if (card == 2)
        {
          hptr[j] = 1;
          sitehets[j]++;
          resolved[i] = 0;
          unresolved_count++;
          if (weight[j] == 1) singletons[i]++;
        }
        else if (card == 3)
        {
          fprintf(stderr, "%c not allowed in sequence %ld at locus %d\n",
                  msa->sequence[i][j], i, msa_index);
        }
        else
          *cleandata = 0;
      }
    }
    
    /* move to the next row of h matrix */
    hptr += msa->length;
  }

  /* 2. Create a list of indices to singleton sites with at least one het */
  for (i=0,k=0; i < msa->length; ++i)
    if (weight[i] == 1 && sitehets[i])
      single_indices[k++] = i;

  /* 3. Resolve one individual at a time */
  for (round = 0; round < msa->count && unresolved_count; ++round)
  {

    /* sort list of indices by number of hets (descending order) */
    ext_qsort_arg = sitehets;
    qsort(single_indices,(size_t)k,sizeof(long),cb_cmp_indices);

    /* find most variable singleton site */
    long chosen = -1;
    for (i = 0; i < k; ++i)
    {
      long site = single_indices[i];

      /* find least variable sequence at site */
      hptr = hmat;
      long y = msa->length+1;
      for (j = 0; j < msa->count; ++j)  
      {
        if (resolved[j] || hptr[site] == 0) 
        {
          hptr += msa->length;
          continue;
        }
        
        if (singletons[j] < y)
        {
          y = singletons[j];
          chosen = j;
        }

        hptr += msa->length;
      }
      
      /* we found a site to resolve */
      if (chosen >= 0)
      {
        hmat[chosen * msa->length + site] = -1;
        sitehets[site]--;
        resolved[chosen] = 1;
        unresolved_count--;

        /* update singleton indices */
        if (sitehets[site] == 0)
        {
          memmove(single_indices+i,single_indices+i+1,(size_t)(k-(i+1)));
          --k;
        }
        break;
      }
    }

    /* if no singletons to resolve, break */
    if (chosen == -1) { break;}
  }

  /* 4. Construct alignment A2 of expanded site patterns. Done in several steps  */

  /* allocate buffer to store number of resolutions for each site in A1 */
  resolution_count = (unsigned long *)xmalloc((size_t)msa->length *
                                              sizeof(unsigned long));

  /* 4a. Determine number of sites for new alignment A2 */
  size_t patterns = (size_t)(msa->length);
  for (i = 0; i < msa->length; ++i)
  {
    assert(sitehets[i] >= 0);
    assert((size_t)sitehets[i] <= sizeof(long) * CHAR_BIT);
    if (sitehets[i])
    {
      size_t temp = patterns;   /* overflow check */
      patterns += (1ul << sitehets[i]) - 1;
      resolution_count[i] = (1ul << sitehets[i]);
      assert(patterns > temp);
    }
    else
      resolution_count[i] = 1;
  }

  /* 4b. Determine number of sequences A2 consists of, allocate space and
     determine mapping of sequence indices from A1 to A2 */
  long newseq_count = 0;
  for (i = 0; i < msa->count; ++i)
    newseq_count += (diploid[i]) ? 2 : 1;
    
  char ** newseq = (char **)xmalloc((size_t)newseq_count*sizeof(char *));
  for (i = 0; i < newseq_count; ++i)
    newseq[i] = (char *)xmalloc((patterns+1)*sizeof(char));

  long * mapping = (long *)xmalloc((size_t)(msa->count)*sizeof(long));
  for (i=0,k=0; i < msa->count; ++i)
  {
    mapping[i] = k++;
    if (diploid[i]) ++k;
  }

  newlabel = (char **)xmalloc((size_t)newseq_count * sizeof(char *));
  for (i=0; i < msa->count; ++i)
  {
    k = mapping[i];
    if (diploid[i])
    {
      char * s = NULL;
      xasprintf(&s, "%s.1", msa->label[i]);
      newlabel[k] = s;
      s = NULL;
      xasprintf(&s, "%s.2", msa->label[i]);
      newlabel[k+1] = s;
    }
    else
    {
      newlabel[k] = xstrdup(msa->label[i]);
    }
  }
  #if 0
  printf("npatt in A2 after expansion = %ld\n", patterns);
  #endif

  /* 4c. loop through sites in heterogenous alignment A1 to generate resolved
     site patterns for alignment A2 */
  long * hets = (long *)xmalloc((size_t)(msa->count) * sizeof(long));
  long q = 0;
  char * newsite = (char *)xmalloc((size_t)newseq_count * sizeof(char));
  for (i = 0; i < msa->length; ++i)
  {
    hptr = hmat;

    for (j=0, n=0; j < msa->count; ++j)
    {
      k = mapping[j];           /* get index of sequence j in A2 */

      if (hptr[i] == 0)         /* haploid or homo */
      {
        newsite[k] = msa->sequence[j][i];
        if (diploid[j])
        {
          newsite[k+1] = msa->sequence[j][i];
        }
      }
      else if (hptr[i] == -1)    /* fixed resolution */
      {
        unsigned int state = map[(int)(msa->sequence[j][i])];

        assert(diploid[j]);
        assert(PLL_POPCOUNT(state) == 2);
        assert(state < 16);  /* currently only DNA */

        unsigned int state1 = 1 << PLL_CTZ(state);
        unsigned int state2 = state & ~state1;

        newsite[k]   = inv_charmap[state1];
        newsite[k+1] = inv_charmap[state2];
      }
      else
        hets[n++] = j;

      /* point to next row */
      hptr += msa->length;
    }

    assert(n == sitehets[i]);

    /* generate resolved site patterns in A2 */
    for (j = 0; j < (1ul << n); ++j)
    {
      for (k=0,m=j; k < n; ++k)
      {
        long index = m & 1;
        m >>= 1;


        long i1 = hets[n-1-k];
        unsigned int state = map[(int)(msa->sequence[i1][i])];

        unsigned int state1 = 1 << PLL_CTZ(state);
        unsigned int state2 = state & ~state1;

        if (index)
          SWAP(state1,state2);

        newsite[mapping[i1]]   = inv_charmap[state1];
        newsite[mapping[i1]+1] = inv_charmap[state2];
      }

      for (k = 0; k < newseq_count; ++k) newseq[k][q] = newsite[k];
      ++q;
    }
  }

  /* add terminating zero to phased sequences */
  for (i = 0; i < newseq_count; ++i)
    newseq[i][patterns] = 0;

  /* 5a. replace alignment a1 with alignment a2 */
  for (i = 0; i < msa->count; ++i)
  {
    free(msa->label[i]);
    free(msa->sequence[i]);
  }
  free(msa->label);
  free(msa->sequence);

  msa->sequence = newseq;
  msa->label = newlabel;
  msa->count = newseq_count;
  msa->length = patterns;

  /* 5b. compress for JC69 */
  /* This is done outside of this function */

  /* free allocated memory */
  free(hets);
  free(mapping);
  free(newsite);
  free(resolved);
  free(diploid);
  free(seqhets);
  free(singletons);
  free(sitehets);
  free(single_indices);
  free(hmat);

  return resolution_count;
}

unsigned long ** diploid_resolve(stree_t * stree,
                                 msa_t ** msa_list,
                                 list_t * maplist,
                                 unsigned int ** weights,
                                 int msa_count,
                                 const unsigned int * map)
{
  long i;
  int * cleandata;
  unsigned long ** resolution_count;

  cleandata = (int *)xcalloc((size_t)msa_count,sizeof(int));

  /* init hash tables */
  if (stree->tip_count == 1)
    sht = mht = NULL;
  else
    diploid_resolution_init(stree,maplist);

  resolution_count = (unsigned long **)xmalloc((size_t)msa_count *
                                               sizeof(unsigned long *));
  for (i = 0; i < msa_count; ++i)
  {
    resolution_count[i] = diploid_resolve_locus(msa_list[i],
                                                (int)i,
                                                weights[i],
                                                cleandata+i,
                                                map);
  }

  /* update map file with new labels */
  if (stree->tip_count > 1)
    update_map_list(maplist);

  /* deallocate hash tables */
  if (stree->tip_count > 1)
    diploid_resolution_fini();

  free(cleandata);

  return resolution_count;
}
