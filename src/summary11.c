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

/* A11 method summary */
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

typedef struct db_bitvector_s
{
  uint64_t * bitvector;
  int64_t bits;
  int64_t count;
} db_bitvector_t;

typedef struct db_stree_s
{
  char * newick;
  int64_t species;
  int64_t count;
} db_stree_t;

typedef struct stringfreq_s
{
  char * label;
  int64_t count;
  int64_t word_count;
} stringfreq_t;

static void cb_stringfreq_dealloc(void * data)
{
  stringfreq_t * sf = data;
  if (sf->label)
    free(sf->label);
  free(sf);
}

/* serialize hashtable of string frequencies into an array */
static stringfreq_t ** hashtable_serialize(hashtable_t * ht)
{
  unsigned long i,k;
  stringfreq_t ** sflist;

  /* allocate the array to hold all items in the hash table */
  sflist = (stringfreq_t **)xmalloc((size_t)(ht->entries_count+1) *
                                   sizeof(stringfreq_t *));

  /* go through all entries of the hashtable and serialize */
  for (i = 0, k = 0; i < ht->table_size; ++i)
  {
    list_t * list = ht->entries[i];

    list_item_t * head = list->head;
    while (head)
    {
      ht_item_t * hi = (ht_item_t *)(head->data);
      sflist[k++] = (stringfreq_t *)(hi->value);
      head = head->next;
    }
  }

  /* serialized items must equal the number of items in hashtable */
  assert(k == ht->entries_count);

  return sflist;
}

/* comparator for strings */
static int cb_delimit_strcmp(const void * s1, const void * s2)
{
  const char * a = *(const char **)s1;
  const char * b = *(const char **)s2;

  return strcmp(a, b);
}

/* recursively fill buf (starting at position index) with all tip nodes of
   subtree rooted at node */
static void snode_getleaves(const snode_t * node, int * index, char ** buf)
{
  if (!(node->left))
  {
    buf[*index] = xstrdup(node->label);
    *index = *index+1;
    return;
  }

  snode_getleaves(node->left,index,buf);
  snode_getleaves(node->right,index,buf);
}

static char * delimit_string(const snode_t * root)
{
  int index = 0;
  long i;
  char ** labels = (char **)xmalloc((size_t)(root->leaves) * sizeof(char *));

  snode_getleaves(root,&index,labels);

  qsort(labels,root->leaves,sizeof(char *),cb_delimit_strcmp);

  long sumsize = 0;
  for (i = 0; i < root->leaves; ++i)
    sumsize += strlen(labels[i]);
  sumsize++;

  char * l = (char *)xmalloc((size_t)sumsize * sizeof(char));
  l[0] = 0;

  char * p = l;
  for (i = 0; i < root->leaves; ++i)
  {
    size_t len = strlen(labels[i]);
    memcpy(p,labels[i],len*sizeof(char));
    p+=len;
  }
  *p = 0;

  for (i = 0; i < root->leaves; ++i)
    free(labels[i]);
  free(labels);
  
  return l;
}

static char * stree_export_delimitation_recursive(const snode_t * root)
{
  char * newick;
  int size_alloced;
  assert(root != NULL);

  if (!(root->left) || !(root->right))
    size_alloced = xasprintf(&newick, "%s", root->label);
  else
  {
    if (root->tau)
    {
      char * subtree1 = stree_export_delimitation_recursive(root->left);
      if (!subtree1)
        fatal("Unable to allocate enough memory.");

      char * subtree2 = stree_export_delimitation_recursive(root->right);
      if (!subtree2)
        fatal("Unable to allocate enough memory.");

      size_alloced = xasprintf(&newick, "(%s, %s)", subtree1, subtree2);
      free(subtree1);
      free(subtree2);
    }
    else
    {
      char * subtree = delimit_string(root);
      size_alloced = xasprintf(&newick, "%s", subtree);
      free(subtree);
    }
  }

  if (size_alloced < 0)
    fatal("memory allocation during newick export failed");

  return newick;
}

static char * stree_export_delimitation(const snode_t * root)
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
    size_alloced = xasprintf(&newick, "%s", root->label);
  else
  {
    if (root->tau)
    {
      char * subtree1 = stree_export_delimitation_recursive(root->left);
      if (!subtree1)
        fatal("Unable to allocate enough memory.");

      char * subtree2 = stree_export_delimitation_recursive(root->right);
      if (!subtree2)
        fatal("Unable to allocate enough memory.");

      size_alloced = xasprintf(&newick, "(%s, %s);", subtree1, subtree2);
      free(subtree1);
      free(subtree2);
    }
    else
    {
      char * subtree = delimit_string(root);
      size_alloced = xasprintf(&newick, "%s;", subtree);
      free(subtree);
    }
  }

  if (size_alloced < 0)
    fatal("memory allocation during newick export failed");

  return newick;
}

static int logint64_len(int64_t x)
{
  return x ? (int)floor(log10(labs(x)))+1 : 1;
}

/* utterly ridiculous way of creating a delimitation string - should redesign
 * at some point */
static char * create_delim_string(stree_t * stree)
{
  long i;

  char ** labels = (char **)xmalloc((size_t)(stree->tip_count)*sizeof(char *));
  for (i = 0; i < stree->tip_count; ++i)
    labels[i] = xstrdup(stree->nodes[i]->label);

  qsort(labels,stree->tip_count,sizeof(char *),cb_delimit_strcmp);

  long sumlen = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    sumlen += strlen(labels[i]);
  }

  /* for spaces */
  sumlen += stree->tip_count - 1;

  char * s = (char *)xmalloc((size_t)(sumlen+1) * sizeof(char));

  sumlen = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    size_t len = strlen(labels[i]);

    memcpy(s+sumlen,labels[i],len*sizeof(char));
    sumlen += len;
    s[sumlen++] = ' ';
  }
  s[--sumlen] = 0;

  for (i = 0; i < stree->tip_count; ++i)
    free(labels[i]);
  free(labels);

  return s;
}

static int64_t get_int64(const char * pline, int64_t * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(pline);
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

  ret = sscanf(start, "%" PRIu64 "%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static void recursive_age(snode_t * node)
{
  if (!node->left)
  {
    node->tau = 0;
    return;
  }

  recursive_age(node->left);
  recursive_age(node->right);

  node->tau = node->left->tau + node->left->length;
}

static int cb_stree_count(const void * pa, const void * pb)
{
  const db_stree_t * a = (const db_stree_t *)pa;
  const db_stree_t * b = (const db_stree_t *)pb;

  if (a->count < b->count) return 1;
  else if (a->count > b->count) return -1;

  return 0;
}
static int cb_stree_species(const void * pa, const void * pb)
{
  db_stree_t * a = (db_stree_t *)pa;
  db_stree_t * b = (db_stree_t *)pb;

  if (a->species > b->species) return 1;
  if (a->species < b->species) return -1;

  return 0;
  //return (int)(a->species - b->species);
}

static int cb_stree_strcmp(const void * pa, const void * pb)
{
  const db_stree_t * sa = (db_stree_t *)pa;
  const db_stree_t * sb = (db_stree_t *)pb;

  const char * a = sa->newick;
  const char * b = sb->newick;

  return strcmp(a, b);
}

static int cb_cmp_label(void * a, void * b)
{
  stringfreq_t * sf = (stringfreq_t *)a;
  char * label =  (char *)b;

  return !strcmp(sf->label,label);
}
static char * get_delimit_string(stree_t * stree)
{
  recursive_age(stree->root);
  return stree_export_delimitation(stree->root);
}

static int cb_countcmp(const void * a, const void * b)
{
  const stringfreq_t * pa = *((const stringfreq_t **)a);
  const stringfreq_t * pb = *((const stringfreq_t **)b);

  if (pa->count < pb->count) return 1;
  else if (pa->count > pb->count) return -1;

  return 0;
}

static void strip_theta_attributes(char * s)
{
  char * p = s;

  while (*s)
  {
    if (*s == '#')
    {
      while (*s && *s != ',' && *s != ')' && *s != ';' && *s != ':')
        ++s;
    }
    else if (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n')
      ++s;
    else
      *p++ = *s++;
  }
  *p = 0;
}

static void replace(db_stree_t * treelist,
                    int64_t index,
                    int64_t start,
                    int64_t freq)
{
  if (index != start)
  {
    if (treelist[index].newick)
      free(treelist[index].newick);
    treelist[index].newick = treelist[start].newick;
    treelist[start].newick = NULL;

    treelist[index].species = treelist[start].species;
  }
  treelist[index].count = freq;
}

void mixed_summary(FILE * fp_out)
{
  int64_t line_count = 0;
  int64_t i,j;
  FILE * fp_mcmc;
  db_stree_t * treelist;
  snode_t ** inner;

  /* open MCMC file for reading */
  fp_mcmc = xopen(opt_mcmcfile,"r");

  /* allocate space for storing inner nodes */
  inner = (snode_t **)xmalloc((size_t)opt_max_species_count*sizeof(snode_t *));

  /* allocate space for reading trees from MCMC file */
  treelist = (db_stree_t *)xmalloc((size_t)(opt_samples+1)*sizeof(db_stree_t));

  /* read trees and species counts from MCMC file */
  while (getnextline(fp_mcmc))
  {
    /* separate line into two zero-terminated strings, the first one (line)
       contains the newick tree string and the second (tmp) holds the species
       count */
    char * tmp =  strchr(line,';');
    tmp++;
    *tmp = 0;
    tmp++;

    /* parse newick string and unambiguously sort tree by its labels */
    if (opt_est_theta)
      strip_theta_attributes(line);
    stree_t * t = stree_parse_newick_string(line);
    stree_sort(t);

    /* convert expanded tree into delimited tree, e.g.:

       ((A:0,B:0):0.02,(C:0.01,D:0.01):0.01); -> (AB,(C,D)); */
    char * dnewick = get_delimit_string(t);

    int64_t species_count;
    if (!get_int64(tmp,&species_count))
      fatal("Cannot read number of species; line %ld of %s",line_count+1,opt_mcmcfile);

    treelist[line_count].newick = dnewick;
    treelist[line_count].species = species_count;
    treelist[line_count].count = 1;
    line_count++;

   // if (species_count > max_species_count)
   //   max_species_count = species_count;

    stree_destroy(t,NULL);
  }

  /* sort by number of species (ascending order - not relevant) */
  qsort(treelist,(size_t)line_count, sizeof(db_stree_t), cb_stree_species);

  /* now sort by newick each range of trees having the same number of species */
  int64_t start = 0;
  int64_t freq = 1; 
  for (i = 1; i < line_count; ++i)
  {
    if (treelist[i].species == treelist[i-1].species)
      ++freq;
    else
    {
      qsort(treelist+start,(size_t)freq,sizeof(db_stree_t),cb_stree_strcmp);
      start = i;
      freq = 1;
    }
  }
  if (freq > 1)
    qsort(treelist+start,(size_t)freq,sizeof(db_stree_t),cb_stree_strcmp);

  /* now keep only the unique trees in treelist */
  int64_t index = 0;
  start = 0; freq = 1;
  for (i = 1; i < line_count; ++i)
  {
    if (!strcmp(treelist[i].newick,treelist[i-1].newick))
    {
      ++freq;
    }
    else
    {
      replace(treelist,index,start,freq);
      ++index;

      start = i;
      freq = 1;

    }
  }
  replace(treelist,index,start,freq);
  ++index;

  /* Print summary statistics (A) with the following columns: 
  
     1. Frequency
     2. Posterior
     3. Cummulative posterior
     4. Number of species
     5. List of delimited species
     6. Newick tree with delimited species
  */
  
  /* create two hashtables for indexing species (ht_species) and delimitations
     (ht_delims) */
  hashtable_t * ht_species = hashtable_create(100*treelist[0].species);
  hashtable_t * ht_delims  = hashtable_create(100*treelist[0].species);

  qsort(treelist,(size_t)index,sizeof(db_stree_t),cb_stree_count);

  int maxlen = logint64_len(treelist[0].count);
  double prob;
  double cum = 0;
  start = 0; freq = 1;
  fprintf(stdout,
          "\n(A) List of best models (count postP #species SpeciesTree)\n");
  fprintf(fp_out,
          "\n(A) List of best models (count postP #species SpeciesTree)\n");
  for (i = 0; i < index; ++i)
  {
    /* TODO: This is an ugly hack to disable checking that the tip labels of
       the loaded tree are equal to the ones given in the control file */
    char * dbgtmp = opt_reorder;
    opt_reorder = NULL;
    stree_t * t = stree_parse_newick_string(treelist[i].newick);
    opt_reorder = dbgtmp;

    /* create delimitation string */
    char * delim = create_delim_string(t);

    /* print summary statistics */
    cum += treelist[i].count / (double)opt_samples;
    prob = treelist[i].count / (double)opt_samples;
    fprintf(stdout,
            "%*ld %f %f %ld ",
            maxlen,treelist[i].count,prob,cum,treelist[i].species);
    fprintf(fp_out,
            "%*ld %f %f %ld ",
            maxlen,treelist[i].count,prob,cum,treelist[i].species);
    fprintf(stdout, " (%s) ", delim);
    fprintf(fp_out, " (%s) ", delim);
    fprintf(stdout, " %s\n", treelist[i].newick);
    fprintf(fp_out, " %s\n", treelist[i].newick);


    /* query delimitation string against hash table */
    stringfreq_t * query = hashtable_find(ht_delims,
                                          (void *)delim,
                                          hash_fnv(delim),
                                          cb_cmp_label);
    if (query)
    {
      /* delimitation found; increase frequency counter and deallocate delim */
      query->count += treelist[i].count;
      free(delim);
    }
    else
    {
      /* delimitation was not found in hashtable; create a new record */
      query = (stringfreq_t *)xmalloc(sizeof(stringfreq_t));
      query->count = treelist[i].count;
      query->label = delim;
      query->word_count = t->tip_count;
      hashtable_insert_force(ht_delims,
                             (void *)query,
                             hash_fnv(query->label));
    }

    /* query each individual species against hash table */
    for (j = 0; j < t->tip_count; ++j)
    {
      stringfreq_t * sf = hashtable_find(ht_species,
                                         (void *)(t->nodes[j]->label),
                                         hash_fnv(t->nodes[j]->label),
                                         cb_cmp_label);
      if (sf)
      {
        sf->count += treelist[i].count;
      }
      else
      {
        sf = (stringfreq_t *)xmalloc(sizeof(stringfreq_t));
        sf->count = treelist[i].count;
        sf->label = xstrdup(t->nodes[j]->label);
        hashtable_insert_force(ht_species,
                               (void *)sf,
                               hash_fnv(sf->label));
      }
    }
    stree_destroy(t,NULL);
  }
  for (i = 0; i < line_count; ++i)
  {
    if (treelist[i].newick)
      free(treelist[i].newick);
  }
  free(treelist);

  /* Print delimitation summary statistics (B) with the following columns:
  
     1. Frequency
     2. Posterior
     3. Number of species
     4. List of delimited species
  */

  /* serialize hashtable of delimitation frequencies into an array */
  stringfreq_t ** dfreqs = hashtable_serialize(ht_delims);

  /* sort delimitations by frequency count (descending order) */
  qsort(dfreqs,ht_delims->entries_count,sizeof(stringfreq_t *),cb_countcmp);

  /* print statistics */
  fprintf(stdout,
          "\n(B) %ld species delimitations & their posterior probabilities\n",
          ht_delims->entries_count);
  fprintf(fp_out,
          "\n(B) %ld species delimitations & their posterior probabilities\n",
          ht_delims->entries_count);
  maxlen = logint64_len(dfreqs[0]->count);
  for (i = 0; i < ht_delims->entries_count; ++i)
  {
    fprintf(stdout, "%*ld %f %3ld (%s)\n", 
            maxlen,
            dfreqs[i]->count,
            dfreqs[i]->count / (double)opt_samples,
            dfreqs[i]->word_count,
            dfreqs[i]->label);
    fprintf(fp_out, "%*ld %f %3ld (%s)\n", 
            maxlen,
            dfreqs[i]->count,
            dfreqs[i]->count / (double)opt_samples,
            dfreqs[i]->word_count,
            dfreqs[i]->label);
  }

  /* compute summary statistics (posterior for # of species) for section (D) */
  double * posterior = (double *)xcalloc((size_t)(opt_max_species_count+1),
                                         sizeof(double));
  for (i = 0; i < ht_delims->entries_count; ++i)
    posterior[dfreqs[i]->word_count] += dfreqs[i]->count;

  for (i = 1; i <= opt_max_species_count; ++i)
    posterior[i] /= opt_samples;
  free(dfreqs);

  /* Print summary statistics (C) for each delimited species:
  
     1. Frequency
     2. Posterior
     3. Putative species 
  */

  /* serialize hashtable of delimited species frequencies into an array */
  dfreqs = hashtable_serialize(ht_species);

  /* sort delimited species labels by frequency count (descending order) */
  qsort(dfreqs,ht_species->entries_count,sizeof(stringfreq_t *),cb_countcmp);

  /* print statistics */
  fprintf(stdout,
          "\n(C) %ld delimited species & their posterior probabilities\n",
          ht_species->entries_count);
  fprintf(fp_out,
          "\n(C) %ld delimited species & their posterior probabilities\n",
          ht_species->entries_count);
  maxlen = logint64_len(dfreqs[0]->count);
  for (i = 0; i < ht_species->entries_count; ++i)
  {
    fprintf(stdout, "%*ld %f %s\n",
            maxlen,
            dfreqs[i]->count,
            dfreqs[i]->count / (double)opt_samples,
            dfreqs[i]->label);
    fprintf(fp_out, "%*ld %f %s\n",
            maxlen,
            dfreqs[i]->count,
            dfreqs[i]->count / (double)opt_samples,
            dfreqs[i]->label);
  }

  free(dfreqs);

  /* Print summary statistics (D) for each delimited species:
  
     1. Frequency
     2. Posterior
     3. Number of species
     4. List of delimited species
  */

  fprintf(stdout, "\n(D) Posterior probability for # of species\n");
  fprintf(fp_out, "\n(D) Posterior probability for # of species\n");
                 
  maxlen = logint64_len(opt_max_species_count);
  double * prior_A11 = getpriorA11();
  if (opt_delimit_prior == BPP_SPECIES_PRIOR_SLH ||
      opt_delimit_prior == BPP_SPECIES_PRIOR_SUNIFORM)
    for (i = 0; i < opt_max_species_count; ++i)
      prior_A11[i] = 1.0 / opt_max_species_count;
  for (i = 1; i <= opt_max_species_count; ++i)
  {
    fprintf(stdout, "P[%*ld] = %f  prior[%*ld] = %f\n",
            maxlen, i, posterior[i], maxlen, i, prior_A11[i-1]);
    fprintf(fp_out, "P[%*ld] = %f  prior[%*ld] = %f\n",
            maxlen, i, posterior[i], maxlen, i, prior_A11[i-1]);
  }
  free(posterior);
                 
                 
  hashtable_destroy(ht_species,cb_stringfreq_dealloc);
  hashtable_destroy(ht_delims,cb_stringfreq_dealloc);
                 
  free(inner);   
  fclose(fp_mcmc);
}                
