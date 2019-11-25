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

#define SWAP_CLV_INDEX(n,i) ((n)+((i)-1)%(2*(n)-2))
#define SWAP_PMAT_INDEX(e,i) (i) = (((e)+(i))%((e)<<1))
#define SWAP_SCALER_INDEX(n,i) (((n)+((i)-1))%(2*(n)-2))

#define GET_HINDEX(t,p) (((node_is_mirror((p)) ? \
                          (p)->node_index : (p)->hybrid->node_index)) - \
                        ((t)->tip_count+(t)->inner_count))

/* TODO: The below two functions are duplicated several times. Move them as
   public functions in gtree.c */
static void all_partials_recursive(gnode_t * node,
                                   unsigned int * trav_size,
                                   gnode_t ** outbuffer)
{
  if (!node->left)
    return;

  all_partials_recursive(node->left,  trav_size, outbuffer);
  all_partials_recursive(node->right, trav_size, outbuffer);

  outbuffer[*trav_size] = node;
  *trav_size = *trav_size + 1;
}

static void gtree_all_partials(gnode_t * root,
                               gnode_t ** travbuffer,
                               unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return;

  all_partials_recursive(root, trav_size, travbuffer);
}

static void dealloc_locus_data(locus_t * locus)
{
  unsigned int i;

  if (!locus) return;

  free(locus->rates);
  free(locus->rate_weights);
  free(locus->eigen_decomp_valid);
  if (!locus->pattern_weights)
    free(locus->pattern_weights);

  if (locus->scale_buffer)
    for (i = 0; i < locus->scale_buffers; ++i)
      free(locus->scale_buffer[i]);
  free(locus->scale_buffer);

  if (locus->tipchars)
    for (i = 0; i < locus->tips; ++i)
      pll_aligned_free(locus->tipchars[i]);
  free(locus->tipchars);

  if (locus->ttlookup)
    pll_aligned_free(locus->ttlookup);

  if (locus->charmap)
    free(locus->charmap);

  if (locus->tipmap)
    free(locus->tipmap);

  if (locus->clv)
  {
    int start = (locus->attributes & PLL_ATTRIB_PATTERN_TIP) ?
                    locus->tips : 0;
    for (i = start; i < locus->clv_buffers + locus->tips; ++i)
      pll_aligned_free(locus->clv[i]);
  }
  free(locus->clv);

  if (locus->pmatrix)
  {
    //for (i = 0; i < partition->prob_matrices; ++i)
      pll_aligned_free(locus->pmatrix[0]);
  }
  free(locus->pmatrix);

  if (locus->subst_params)
    for (i = 0; i < locus->rate_matrices; ++i)
      pll_aligned_free(locus->subst_params[i]);
  free(locus->subst_params);

  if (locus->eigenvecs)
    for (i = 0; i < locus->rate_matrices; ++i)
      pll_aligned_free(locus->eigenvecs[i]);
  free(locus->eigenvecs);

  if (locus->inv_eigenvecs)
    for (i = 0; i < locus->rate_matrices; ++i)
      pll_aligned_free(locus->inv_eigenvecs[i]);
  free(locus->inv_eigenvecs);

  if (locus->eigenvals)
    for (i = 0; i < locus->rate_matrices; ++i)
      pll_aligned_free(locus->eigenvals[i]);
  free(locus->eigenvals);

  if (locus->frequencies)
    for (i = 0; i < locus->rate_matrices; ++i)
      pll_aligned_free(locus->frequencies[i]);
  free(locus->frequencies);

  free(locus->param_indices);
  free(locus->heredity);

  if (locus->pattern_weights)
    free(locus->pattern_weights);

  if (locus->diploid)
  {
    free(locus->diploid_mapping);
    free(locus->diploid_resolution_count);
    free(locus->likelihood_vector);
  }

  free(locus);
}

#if 0
static int mark_ambiguous_sites(msa_t * msa,
                                const unsigned int * map)
{
  int i,j,k;
  unsigned int c;
  int clean_count = 0;

  int len = msa->length;
  int count = msa->count;

  char ** seq = msa->sequence;

  /* iterate sites */
  for (i = 0; i < len; ++i)
  {
    /* iterate sequences */
    for (j = 0; j < count; ++j)
    {
      /* convert character to state code */
      c = (map[(unsigned int)(seq[j][i])]);

      /* get number of bits set */
      k = __builtin_popcount(c);

      /* if no bits are set then invalid state */
      if (!k)
        fatal("Illegal state code \%c\" in tip %s",
              seq[j][i], msa->label[j]);

      /* if ambiguity, set the whole site to an ambiguous state and increase
         the number of sites to be removed */
      if (k > 1)
      {
        c = seq[j][i];
        for (k = 0; k < count; ++k)
          seq[k][i] = c;
        clean_count++;
        break;
      }
    }
  }

  return clean_count;
}

static void set_tipstates(locus_t * locus,
                          msa_t * msa,
                          const unsigned int * map)
{
  int i,j,k;
  unsigned int c;
  int clean_count = 0;
  
  int len = msa->length;
  int count = msa->count;

  /* if --cleandata then get number of sites to be removed */
  if (opt_cleandata)
    clean_count = mark_ambiguous_sites(msa,map);

  /* allocate necessary space for storing tip sequences */
  locus->tipstates = (unsigned char **)xmalloc(count*sizeof(unsigned char *));
  for (i = 0; i < count; ++i)
    locus->tipstates[i] = (unsigned char *)xmalloc((len-clean_count) *
                                                   sizeof(unsigned char));
    
  /* iterate sequences */
  for (i = 0; i < count; ++i)
  {
    char * seq = msa->sequence[i];

    /* go through the current sequence */
    for (j=0, k=0; j < len; ++j)
    {
      /* convert sequence character to state code */
      c = map[(unsigned int)(seq[j])];

      /* seq[j] is an illegal character */
      if (!c)
        fatal("Illegal state code \%c\" in tip %s",
              seq[j], msa->label[i]);

      /* if cleandata enabled and ambiguous state then skip */
      if (opt_cleandata && __builtin_popcount(c) > 1)

      /* store state */
      locus->tipstates[i][k++] = c;
    }

  }
}

void init_locus(locus_t * locus, msa_t * msa)
{
  /* convert alignments from ASCII character representation to state nucleotide
     state representation and store them in the locus structure. If cleandata
     is specified, remove all sites that contain an ambiguity at some sequence
  */

  set_tipstates(locus, msa, pll_map_nt);
}
#endif

void pll_set_pattern_weights(locus_t * locus,
                             const unsigned int * pattern_weights)
{
  unsigned int i;

  memcpy(locus->pattern_weights,
         pattern_weights,
         sizeof(unsigned int)*locus->sites);

  locus->pattern_weights_sum = 0;
  for (i = 0; i < locus->sites; ++i)
    locus->pattern_weights_sum += pattern_weights[i];

}

void pll_set_category_rates(locus_t * locus, const double * rates)
{
  memcpy(locus->rates, rates, locus->rate_cats*sizeof(double));
}

static int update_charmap(locus_t * locus, const unsigned int * map)
{
  unsigned int i,j,k = 0;
  unsigned int new_states_count = 0;
  unsigned int mapcopy[ASCII_SIZE];

  memcpy(mapcopy, map, ASCII_SIZE * sizeof(unsigned int));

  /* find maximum value in charmap table */
  k = 0;
  while (locus->tipmap[k]) ++k;

  /* compute the number of new states in the map */
  for (i = 0; i < ASCII_SIZE; ++i)
  {
    if (mapcopy[i])
    {
      /* check whether state map[i] already exists in the tipmap */
      for (j = 0; j < k; ++j)
        if (mapcopy[i] == locus->tipmap[j])
          break;

      /* if it does not exist */
      if (j == k)
      {
        /* check whether it is the first time we find it */
        for (j = 0; j < i; ++j)
          if (mapcopy[j] == mapcopy[i])
            break;

        /* if first time, increase number of new states */
        if (j == i) new_states_count++;
      }
    }
  }

  /* erase old charmap */
  memset(locus->charmap,0,ASCII_SIZE*sizeof(unsigned char));

  /* using this map we will have more than 256 states, so return an error */
  if (new_states_count + k >= ASCII_SIZE)
    fatal("Cannot specify 256 or more states with PLL_ATTRIB_PATTERN_TIP.");

  /* traverse the new map */
  for (i = 0; i < ASCII_SIZE; ++i)
  {
    if (mapcopy[i])
    {
      unsigned int code = 0;

      /* check whether state map[i] already exists in the tipmap */
      for (j = 0; j < k; ++j)
        if (mapcopy[i] == locus->tipmap[j])
          break;

      if (j == k)
      {
        /* if it does not exist */
        code = k;
        locus->tipmap[code] = mapcopy[i];
        ++k;
      }
      else
      {
        /* if it exists already then j is the index to the tipmap */
        code = j;
      }

      locus->charmap[i] = code;

      /* find all characters with the same state in the map */
      for (j=i+1; j < ASCII_SIZE; ++j)
      {
        if (mapcopy[i] == mapcopy[j])
        {
          locus->charmap[j] = code;
          mapcopy[j] = 0;
        }
      }
    }
  }

  /* set maximum number of states (including ambiguities), its logarithm,
     and the logarithm of states */
  if (new_states_count)
  {
    /* special cases which do not use remapping */
    if (locus->states == 4)
    {
      for (k = 0, i = 0; locus->tipmap[i]; ++i)
        if (locus->tipmap[i] > k)
          k = locus->tipmap[i];

      locus->maxstates = k+1;
    }
    else
      locus->maxstates += new_states_count;

    unsigned int l2_maxstates = (unsigned int)ceil(log2(locus->maxstates));

    /* allocate space for the precomputed tip-tip likelihood vector */
    size_t alloc_size = (1 << (2 * l2_maxstates)) *
                        (locus->states_padded * locus->rate_cats);

    /* for AVX we do not need to reallocate ttlookup as it has fixed size */
    if ((locus->states == 4) && (locus->attributes & PLL_ATTRIB_ARCH_AVX))
      return BPP_SUCCESS;

    free(locus->ttlookup);
    locus->ttlookup = pll_aligned_alloc(alloc_size * sizeof(double),
                                        locus->alignment);
    if (!locus->ttlookup)
      fatal("Cannot allocate space for storing precomputed tip-tip CLVs.");
  }

  return BPP_SUCCESS;
}

/* create a bijective mapping from states to the range <1,maxstates> where
   maxstates is the maximum number of states (including ambiguities). It is
   neeed to index the precomputed conditional likelihoods for each pair of
   states. The sequences are then encoded using this charmap, and we store
   the precomputated CLV for a charmapped pair i and j, at index:

   (i << ceil(log(maxstate)) + j) << log(states) << log(rates) */
static int create_charmap(locus_t * locus, const unsigned int * usermap)
{
  unsigned int i,j,m = 0;
  unsigned char k = 0;
  unsigned int map[ASCII_SIZE];

  /* If ascertainment bias correction attribute is set, CLVs will be allocated
     with additional sites for each state */
  #if 0
  unsigned int sites_alloc = locus->asc_bias_alloc ?
                   locus->sites + locus->states : locus->sites;
  #endif
  unsigned int sites_alloc = locus->sites;

  //memcpy(map, partition->map, PLL_ASCII_SIZE * sizeof(unsigned int));
  memcpy(map, usermap, ASCII_SIZE * sizeof(unsigned int));

  if (!(locus->charmap = (unsigned char *)xcalloc(ASCII_SIZE,
                                                  sizeof(unsigned char))))

  if (!(locus->tipmap = (unsigned int *)xcalloc(ASCII_SIZE,
                                                sizeof(unsigned int))))


  /* create charmap (remapped table of ASCII characters to range 0,|states|)
     and tipmap which is a (1,|states|) -> state */
  for (i = 0; i < ASCII_SIZE; ++i)
  {
    if (map[i])
    {
      if (map[i] > m) m = map[i];

      locus->charmap[i] = k;
      locus->tipmap[(unsigned int)k] = map[i];
      for (j = i+1; j < ASCII_SIZE; ++j)
      {
        if (map[i] == map[j])
        {
          locus->charmap[j] = k;
          map[j] = 0;
        }
      }
      ++k;
    }
  }

  /* For all state settings for which remapping will not be done, we need to
     increment maxstates by one to account for a fictive state 0 which will
     never be used */
  if (locus->states == 4) k = m+1;

  /* set maximum number of states (including ambiguities), its logarithm,
     and the logarithm of states */
  locus->maxstates = (unsigned int)k;

  unsigned int l2_maxstates = (unsigned int)ceil(log2(locus->maxstates));

  /* allocate space for the precomputed tip-tip likelihood vector */
  size_t alloc_size = (1 << (2 * l2_maxstates)) *
                      (locus->states_padded * locus->rate_cats);

  /* dedicated 4x4 function  - if AVX is not used we can allocate less space
     in case not all 16 possible ambiguities are present */
  if ((locus->states == 4) &&
      (locus->attributes & PLL_ATTRIB_ARCH_AVX))
  {
    locus->ttlookup = pll_aligned_alloc(1024 * locus->rate_cats *
                                        sizeof(double),
                                        locus->alignment);
  }
  else
  {
    locus->ttlookup = pll_aligned_alloc(alloc_size * sizeof(double),
                                        locus->alignment);
  }

  /* allocate tip character arrays */
  locus->tipchars = (unsigned char **)xcalloc(locus->tips,
                                              sizeof(unsigned char *));

  for (i = 0; i < locus->tips; ++i)
    locus->tipchars[i] = (unsigned char *)xmalloc(sites_alloc *
                                                  sizeof(unsigned char));

  return BPP_SUCCESS;
}
static int set_tipchars_4x4(locus_t * locus,
                            unsigned int tip_index,
                            const unsigned int * map,
                            const char * sequence)
{
  unsigned int c;
  unsigned int i;

  /* iterate through sites */
  for (i = 0; i < locus->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
      fatal("Illegal state code in tip \"%c\"", sequence[i]);

    /* store states as the remapped characters from charmap */
    locus->tipchars[tip_index][i] = (unsigned char)c;
  }

  /* tipmap is never used in the 4x4 case except create and update_charmap */

  return BPP_SUCCESS;
}

static int set_tipchars(locus_t * locus,
                        unsigned int tip_index,
                        const unsigned int * map,
                        const char * sequence)
{
  unsigned int c;
  unsigned int i;

  /* iterate through sites */
  for (i = 0; i < locus->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
      fatal("Illegal state code in tip \"%c\"", sequence[i]);

    /* store states as the remapped characters from charmap */
    locus->tipchars[tip_index][i] = locus->charmap[(int)(sequence[i])];
  }

  return BPP_SUCCESS;
}

static int set_tipclv(locus_t * locus,
                     unsigned int tip_index,
                     const unsigned int * map,
                     const char * sequence)
{
  unsigned int c;
  unsigned int i,j;
  double * tipclv = locus->clv[tip_index];

  /* iterate through sites */
  for (i = 0; i < locus->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
      fatal("Illegal state code in tip \"%c\"", sequence[i]);

    /* decompose basecall into the encoded residues and set the appropriate
       positions in the tip vector */
    for (j = 0; j < locus->states; ++j)
    {
      tipclv[j] = c & 1;
      c >>= 1;
    }

    /* fill in the entries for the other gamma values */
    tipclv += locus->states_padded;
    for (j = 0; j < locus->rate_cats - 1; ++j)
    {
      memcpy(tipclv, tipclv - locus->states_padded,
             locus->states * sizeof(double));
      tipclv += locus->states_padded;
    }
  }

  return BPP_SUCCESS;
}

int pll_set_tip_states(locus_t * locus,
                       unsigned int tip_index,
                       const unsigned int * map,
                       const char * sequence)
{
  int rc;

  if (locus->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    /* create (or update) character map for tip-tip precomputations */
    if (locus->tipchars)
    {
      update_charmap(locus,map);
    }
    else
    {
      if (!create_charmap(locus,map))
      {
        dealloc_locus_data(locus);
        return BPP_FAILURE;
      }
    }

    if (locus->states == 4)
      rc = set_tipchars_4x4(locus, tip_index, map, sequence);
    else
      rc = set_tipchars(locus, tip_index, map, sequence);
  }
  else
    rc = set_tipclv(locus, tip_index, map, sequence);

  return rc;
}

//TODO: <DOC> We should account for padding before calling this function
int pll_set_tip_clv(locus_t * locus,
                    unsigned int tip_index,
                    const double * clv,
                    int padding)
{
  unsigned int i,j;

  if (locus->attributes & PLL_ATTRIB_PATTERN_TIP)
    fatal("Cannot use pll_set_tip_clv with PLL_ATTRIB_PATTERN_TIP.");

  double * tipclv = locus->clv[tip_index];

  for (i = 0; i < locus->sites; ++i)
  {
    for (j = 0; j < locus->rate_cats; ++j)
    {
      memcpy(tipclv, clv, locus->states*sizeof(double));
      tipclv += locus->states_padded;
    }
    clv += padding ? locus->states_padded : locus->states;
  }

  return BPP_SUCCESS;
}


locus_t * locus_create(unsigned int dtype,
                       unsigned int model,
                       unsigned int tips,
                       unsigned int clv_buffers,
                       unsigned int states,
                       unsigned int sites,
                       unsigned int rate_matrices,
                       unsigned int prob_matrices,
                       unsigned int rate_cats,
                       unsigned int scale_buffers,
                       unsigned int attributes)
{
  unsigned int i;
  unsigned int sites_alloc;

  /* TODO: In the case of ascertainment bias correction, change sites_alloc
     to sites+states */
  sites_alloc = sites;

  /* make sure that multiple ARCH were not specified */
  if (PLL_POPCOUNT(attributes & PLL_ATTRIB_ARCH_MASK) > 1)
    fatal("Internal error in setting locus attributes");

  /* allocate locus partition */
  locus_t * locus = (locus_t *)xcalloc(1,sizeof(locus_t));

  /* extract architecture and set vectorization parameters */
  locus->alignment = PLL_ALIGNMENT_CPU;
  locus->attributes = attributes;
  locus->states_padded = states;

  /* by default we assume the locus does not contain diploid sequences */
  locus->diploid = 0;
  locus->diploid_mapping = NULL;
  locus->diploid_resolution_count = NULL;
  locus->likelihood_vector = NULL;

  if (attributes & PLL_ATTRIB_ARCH_SSE)
  {
    locus->alignment = PLL_ALIGNMENT_SSE;
    locus->states_padded = (states+1) & 0xFFFFFFFE;
  }
  if (attributes & PLL_ATTRIB_ARCH_AVX)
  {
    locus->alignment = PLL_ALIGNMENT_AVX;
    locus->states_padded = (states+3) & 0xFFFFFFFC;
  }
  if (attributes & PLL_ATTRIB_ARCH_AVX2)
  {
    locus->alignment = PLL_ALIGNMENT_AVX;
    locus->states_padded = (states+3) & 0xFFFFFFFC;
  }

  unsigned int states_padded = locus->states_padded;

  /* initialize properties */

  locus->dtype = dtype;
  locus->model = model;

  locus->tips = tips;
  locus->clv_buffers = clv_buffers;
  locus->states = states;
  locus->sites = sites;

  assert(opt_alpha_cats == rate_cats);
  if (rate_cats > 1)
  {
    /* set alpha to mean */
    locus->rates_alpha = opt_alpha_alpha / opt_alpha_beta;
  }
  else
    locus->rates_alpha = 1;

  locus->rate_matrices = rate_matrices;
  locus->prob_matrices = prob_matrices;
  locus->rate_cats = rate_cats;
  locus->scale_buffers = scale_buffers;

  locus->pattern_weights = NULL;

  locus->eigenvecs = NULL;
  locus->inv_eigenvecs = NULL;
  locus->eigenvals = NULL;

  locus->rates = NULL;
  locus->rate_weights = NULL;
  locus->subst_params = NULL;
  locus->scale_buffer = NULL;
  locus->frequencies = NULL;
  locus->eigen_decomp_valid = 0;

  locus->ttlookup = NULL;
  locus->tipchars = NULL;
  locus->charmap = NULL;
  locus->tipmap = NULL;

  /* param indices. By default we use the same frequencies/qmatrix for computing
     the pmatrices for each rate category */
  locus->param_indices = (unsigned int *)xmalloc((size_t)locus->rate_cats *
                                                 sizeof(unsigned int));
  for (i = 0; i < locus->rate_cats; ++i)
    locus->param_indices[i] = 0;

  /* eigen_decomp_valid */
  locus->eigen_decomp_valid = (int *)xcalloc(locus->rate_matrices,
                                             sizeof(int));
  /* clv */
  locus->clv = (double **)xcalloc(locus->tips + locus->clv_buffers,
                                      sizeof(double *));

  /* if tip pattern precomputation is enabled, then do not allocate CLV space
     for the tip nodes */
  int start = (locus->attributes & PLL_ATTRIB_PATTERN_TIP) ? locus->tips : 0;

  for (i = start; i < locus->tips + locus->clv_buffers; ++i)
  {
    locus->clv[i] = pll_aligned_alloc(sites_alloc * states_padded * rate_cats *
                                      sizeof(double),
                                      locus->alignment);
    /* zero-out CLV vectors to avoid valgrind warnings when using odd number of
       states with vectorized code */
    memset(locus->clv[i],
           0,
           (size_t)sites_alloc*states_padded*rate_cats*sizeof(double));
  }

  /* pmatrix */
  locus->pmatrix = (double **)xcalloc(locus->prob_matrices, sizeof(double *));

  /* allocate transition probability matrices in contiguous space, in order
     to save the 'displacement' amount of memory per matrix, which is
     required for updating partials when the number of states is not a multiple
     of states_padded. */
  size_t displacement = (states_padded - states)*(states_padded)*sizeof(double);
  locus->pmatrix[0] = pll_aligned_alloc(locus->prob_matrices * states *
                                        states_padded * rate_cats *
                                        sizeof(double) + displacement,
                                        locus->alignment);

  for (i = 1; i < locus->prob_matrices; ++i)
    locus->pmatrix[i] = locus->pmatrix[i-1] + states*states_padded*rate_cats;

  /* zero-out p-matrices to avoid valgrind warnings when using odd number of
     states with vectorized code */
  memset(locus->pmatrix[0],0,
         locus->prob_matrices * states * states_padded * rate_cats *
         sizeof(double) + displacement);

  /* eigenvecs */
  locus->eigenvecs = (double **)xcalloc(locus->rate_matrices,
                                        sizeof(double *));
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    locus->eigenvecs[i] = pll_aligned_alloc(states*states_padded*sizeof(double),
                                            locus->alignment);
    memset(locus->eigenvecs[i], 0, states * states_padded * sizeof(double));
  }

  /* inv_eigenvecs */
  locus->inv_eigenvecs = (double **)xcalloc(locus->rate_matrices,
                                            sizeof(double *));
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    locus->inv_eigenvecs[i] = pll_aligned_alloc(states*states_padded *
                                                sizeof(double),
                                                locus->alignment);
    memset(locus->inv_eigenvecs[i], 0, states*states_padded*sizeof(double));
  }

  /* eigenvals */
  locus->eigenvals = (double **)xcalloc(locus->rate_matrices,sizeof(double *));
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    locus->eigenvals[i] = pll_aligned_alloc(states_padded*sizeof(double),
                                            locus->alignment);
    memset(locus->eigenvals[i], 0, states_padded * sizeof(double));
  }

  /* subst_params */
  locus->subst_params = (double **)xcalloc(locus->rate_matrices,
                                           sizeof(double *));
  for (i = 0; i < locus->rate_matrices; ++i)
    locus->subst_params[i] = pll_aligned_alloc(((states*states-states)/2) *
                                               sizeof(double),
                                               locus->alignment);

  /* frequencies */
  locus->frequencies = (double **)xcalloc(locus->rate_matrices,
                                          sizeof(double *));
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    locus->frequencies[i] = pll_aligned_alloc(states_padded*sizeof(double),
                                              locus->alignment);
    memset(locus->frequencies[i],0,states_padded*sizeof(double));
  }

  /* heredity scalers */
  locus->heredity = (double *)xcalloc(locus->rate_matrices,sizeof(double));
  for (i = 0; i < locus->rate_matrices; ++i)
    locus->heredity[i] = 1;

  /* rates */
  locus->rates = (double *)xcalloc(locus->rate_cats,sizeof(double));
  pll_compute_gamma_cats(locus->rates_alpha,
                         locus->rates_alpha,
                         locus->rate_cats,
                         locus->rates,
                         PLL_GAMMA_RATES_MEAN);
  #if 0
  for (i = 0; i < locus->rate_cats; ++i)
    locus->rates[i] = 1;
  #endif

  /* rate weights */
  locus->rate_weights = (double *)xcalloc(locus->rate_cats,sizeof(double));
    /* initialize to 1/n_rates */
  for (i = 0; i < locus->rate_cats; ++i)
    locus->rate_weights[i] = 1.0 / locus->rate_cats;

  /* site weights */
  locus->pattern_weights = (unsigned int *)xmalloc(sites_alloc *
                                                   sizeof(unsigned int));
  /* implicitely set all weights to 1 */
  for (i = 0; i < locus->sites; ++i)
    locus->pattern_weights[i] = 1;
  locus->pattern_weights_sum = sites;

  /* scale_buffer */
  locus->scale_buffer = (unsigned int **)xcalloc(locus->scale_buffers,
                                                 sizeof(unsigned int *));
  for (i = 0; i < locus->scale_buffers; ++i)
  {
    size_t scaler_size = (attributes & PLL_ATTRIB_RATE_SCALERS) ?
                             sites_alloc * rate_cats : sites_alloc;
    locus->scale_buffer[i] = (unsigned int *)xcalloc(scaler_size,
                                                     sizeof(unsigned int));
  }

  return locus;
}

void locus_destroy(locus_t * locus)
{
  dealloc_locus_data(locus);
}

void pll_set_subst_params(locus_t * locus,
                          unsigned int param_index,
                          const double * params)
{
  unsigned int count = (locus->states * (locus->states-1)) / 2;

  memcpy(locus->subst_params[param_index], params, count*sizeof(double));
  locus->eigen_decomp_valid[param_index] = 0;

  /* NOTE: For protein models PLL/RAxML do a rate scaling by 10.0/max_rate */
}

void pll_set_frequencies(locus_t * locus,
                         unsigned int freqs_index,
                         const double * frequencies)
{
  memcpy(locus->frequencies[freqs_index],
         frequencies,
         locus->states*sizeof(double));
  locus->eigen_decomp_valid[freqs_index] = 0;
}

void locus_set_frequencies_and_rates(locus_t * locus)
{
  double frequencies[4] = {0.25, 0.25, 0.25, 0.25};
  double jcrates[6] = {1,1,1,1,1,1};
  const double * freqs;
  const double * rates;

  if (locus->dtype == BPP_DATA_DNA)
  {
    assert(locus->states == 4);

    if (locus->model == BPP_DNA_MODEL_GTR)
    {
      int i;
      double sum;

      /****** Ziheng 2019-11-14 ****************/
      for (i = 0, sum = 0; i < 4; i++)
        sum += frequencies[i] = 0.8 + 0.4*legacy_rndu(0);
      for (i = 0; i < 4; i++)
        frequencies[i] /= sum;

      for (i = 0, sum = 0; i < 6; i++)
        sum += jcrates[i] = 0.8 + 0.4*legacy_rndu(0);
      for (i = 0; i < 6; i++)
        jcrates[i] /= sum;
    }



    freqs = frequencies;
    rates = jcrates;

  }
  else
  {
    assert(locus->dtype == BPP_DATA_AA);
    assert(locus->states == 20);

    switch (locus->model)
    {
      case BPP_AA_MODEL_DAYHOFF:
        freqs = pll_aa_freqs_dayhoff;
        rates = pll_aa_rates_dayhoff;
        break;

      case BPP_AA_MODEL_LG:
        freqs = pll_aa_freqs_lg;
        rates = pll_aa_rates_lg;
        break;

      case BPP_AA_MODEL_DCMUT:
        freqs = pll_aa_freqs_dcmut;
        rates = pll_aa_rates_dcmut;
        break;

      case BPP_AA_MODEL_JTT:
        freqs = pll_aa_freqs_jtt;
        rates = pll_aa_rates_jtt;
        break;

      case BPP_AA_MODEL_MTREV:
        freqs = pll_aa_freqs_mtrev;
        rates = pll_aa_rates_mtrev;
        break;

      case BPP_AA_MODEL_WAG:
        freqs = pll_aa_freqs_wag;
        rates = pll_aa_rates_wag;
        break;

      case BPP_AA_MODEL_RTREV:
        freqs = pll_aa_freqs_rtrev;
        rates = pll_aa_rates_rtrev;
        break;

      case BPP_AA_MODEL_CPREV:
        freqs = pll_aa_freqs_cprev;
        rates = pll_aa_rates_cprev;
        break;

      case BPP_AA_MODEL_VT:
        freqs = pll_aa_freqs_vt;
        rates = pll_aa_rates_vt;
        break;

      case BPP_AA_MODEL_BLOSUM62:
        freqs = pll_aa_freqs_blosum62;
        rates = pll_aa_rates_blosum62;
        break;

      case BPP_AA_MODEL_MTMAM:
        freqs = pll_aa_freqs_mtmam;
        rates = pll_aa_rates_mtmam;
        break;

      case BPP_AA_MODEL_MTART:
        freqs = pll_aa_freqs_mtart;
        rates = pll_aa_rates_mtart;
        break;

      case BPP_AA_MODEL_MTZOA:
        freqs = pll_aa_freqs_mtzoa;
        rates = pll_aa_rates_mtzoa;
        break;

      case BPP_AA_MODEL_PMB:
        freqs = pll_aa_freqs_pmb;
        rates = pll_aa_rates_pmb;
        break;

      case BPP_AA_MODEL_HIVB:
        freqs = pll_aa_freqs_hivb;
        rates = pll_aa_rates_hivb;
        break;

      case BPP_AA_MODEL_HIVW:
        freqs = pll_aa_freqs_hivw;
        rates = pll_aa_rates_hivw;
        break;

      case BPP_AA_MODEL_JTTDCMUT:
        freqs = pll_aa_freqs_jttdcmut;
        rates = pll_aa_rates_jttdcmut;
        break;

      case BPP_AA_MODEL_FLU:
        freqs = pll_aa_freqs_flu;
        rates = pll_aa_rates_flu;
        break;

      case BPP_AA_MODEL_STMTREV:
        freqs = pll_aa_freqs_stmtrev;
        rates = pll_aa_rates_stmtrev;
        break;

      default:
        fatal("Internal error while selecting protein substitution model");
    }
  }
  
  assert(locus->rate_matrices == 1);

  pll_set_frequencies(locus,0,freqs);
  pll_set_subst_params(locus,0,rates);
}

void locus_set_heredity_scalers(locus_t * locus, const double * heredity)
{
  /* one mutation rate per substitution matrix available */
  memcpy(locus->heredity, heredity, locus->rate_matrices*sizeof(double));
}

static double update_branchlength_relaxed_clock(stree_t * stree,
                                                gnode_t * node,
                                                long msa_index)
{
  /* relaxed clock */
  double t = node->time;
  double length = 0;
  snode_t * start = node->pop;
  snode_t * end   = node->parent->pop;
  assert(end);

  length = 0;

  while (start != end)
  {
    snode_t * pop = start;
    assert(start && start->parent);
    start = start->parent;

    if (start->hybrid)
    {
      assert(!node_is_mirror(start));

      unsigned int hindex = GET_HINDEX(stree,start);
      assert(hindex >= 0 && hindex < stree->hybrid_count);

      /* find correct parent node according to hpath flag */
      assert(node->hpath[hindex] != BPP_HPATH_NONE);
      assert(start->left);
      if (node->hpath[hindex] == BPP_HPATH_RIGHT)
        start = start->hybrid;
    }

    /* skip using branch rates on horizontal edges in hybridization events */
    if (!(pop->hybrid && start->htau == 0))
      length += (start->tau - t)*pop->brate[msa_index];
    t = start->tau;
  }
  length += (node->parent->time - t) * node->parent->pop->brate[msa_index];

  return length;
}

static void locus_update_all_matrices_generic_recursive(locus_t * locus,
                                                        gtree_t * gtree,
                                                        gnode_t * node,
                                                        stree_t * stree,
                                                        long msa_index,
                                                        double * memexpd,
                                                        double * memtemp)
{
  unsigned int n,j,k,m;
  unsigned int states = locus->states;
  unsigned int states_padded = states;
  unsigned int rate_cats = locus->rate_cats;
  unsigned int attrib = locus->attributes;
  double t;
  const double * rates = locus->rates;
  double * const * eigenvals = locus->eigenvals;
  double * const * eigenvecs = locus->eigenvecs;
  double * const * inv_eigenvecs = locus->inv_eigenvecs;
  double * expd;
  double * temp;

  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;

  expd = memexpd;
  temp = memtemp;

  unsigned int * param_indices = locus->param_indices;

  if (opt_clock == BPP_CLOCK_GLOBAL)
  {
    /* strict clock */
    t = node->length = (node->parent->time - node->time)*gtree->rate_mui;
  }
  else
  {
    /* relaxed clock */
    t = update_branchlength_relaxed_clock(stree,node,msa_index);
  }

  assert(t >= 0);

  /* compute effective pmatrix location */
  for (n = 0; n < rate_cats; ++n)
  {
    pmat = locus->pmatrix[node->pmatrix_index] + n*states*states_padded;
    t *= locus->rates[n];

    evecs = eigenvecs[param_indices[n]];
    inv_evecs = inv_eigenvecs[param_indices[n]];
    evals = eigenvals[param_indices[n]];

    /* if branch length is zero then set the p-matrix to identity matrix */
    if (t < 1e-100)
    {
      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k)
          pmat[j*states_padded + k] = (j == k) ? 1 : 0;
    }
    else
    {
      /* NOTE: in order to deal with numerical issues in cases when Qt -> 0, we
       * use a trick suggested by Ben Redelings and explained here:
       * https://github.com/xflouris/libpll/issues/129#issuecomment-304004005
       * In short, we use expm1() to compute (exp(Qt) - I), and then correct
       * for this by adding an identity matrix I in the very end */

      /* exponentiate eigenvalues */
      for (j = 0; j < states; ++j)
        expd[j] = expm1(evals[j] * rates[n] * t);

      for (j = 0; j < states; ++j)
        for (k = 0; k < states; ++k)
          temp[j*states+k] = inv_evecs[j*states_padded+k] * expd[k];

      for (j = 0; j < states; ++j)
      {
        for (k = 0; k < states; ++k)
        {
          pmat[j*states_padded+k] = (j==k) ? 1.0 : 0;
          for (m = 0; m < states; ++m)
          {
            pmat[j*states_padded+k] +=
                temp[j*states+m] * evecs[m*states_padded+k];
          }
        }
      }
    }
    #ifdef DEBUG
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
        assert(pmat[j*states_padded+k] >= 0);
    #endif
  }

  if (!(node->left)) return;

  locus_update_all_matrices_generic_recursive(locus,
                                              gtree,
                                              node->left,
                                              stree,
                                              msa_index,
                                              memexpd,
                                              memtemp);
  locus_update_all_matrices_generic_recursive(locus,
                                              gtree,
                                              node->right,
                                              stree,
                                              msa_index,
                                              memexpd,
                                              memtemp);

}

static void locus_update_all_matrices_generic(locus_t * locus,
                                              gtree_t * gtree,
                                              stree_t * stree,
                                              long msa_index)
{
  unsigned int n;
  double * expd;
  double * temp;

  unsigned int * param_indices = locus->param_indices;

  for (n = 0; n < locus->rate_cats; ++n)
  {
    unsigned int param_index = param_indices[n];

    if (!locus->eigen_decomp_valid[param_index])
    {
      pll_update_eigen(locus->eigenvecs[param_index],
                       locus->inv_eigenvecs[param_index],
                       locus->eigenvals[param_index],
                       locus->frequencies[param_index],
                       locus->subst_params[param_index],
                       locus->states,
                       locus->states_padded);
      locus->eigen_decomp_valid[param_index] = 1;
    }

  }

  expd = (double *)xmalloc(locus->states * sizeof(double));
  temp = (double *)xmalloc(locus->states*locus->states*sizeof(double));

  locus_update_all_matrices_generic_recursive(locus,
                                              gtree,
                                              gtree->root->left,
                                              stree,
                                              msa_index,
                                              expd,
                                              temp);
  locus_update_all_matrices_generic_recursive(locus,
                                              gtree,
                                              gtree->root->right,
                                              stree,
                                              msa_index,
                                              expd,
                                              temp);

  free(expd);
  free(temp);
}

static void locus_update_all_matrices_jc69_recursive(locus_t * locus,
                                                     gtree_t * gtree,
                                                     gnode_t * root,
                                                     stree_t * stree,
                                                     long msa_index)
{
  long n;
  double t;
  double * pmat;
  unsigned int states = locus->states;
  unsigned int states_padded = locus->states_padded;

  if (opt_clock == BPP_CLOCK_GLOBAL)
  {
    /* strict clock */
    t = root->length = (root->parent->time - root->time)*gtree->rate_mui;
  }
  else
  {
    /* relaxed clock */
    t = root->length = update_branchlength_relaxed_clock(stree,root,msa_index);
  }

  for (n = 0; n < locus->rate_cats; ++n)
  {
    pmat = locus->pmatrix[root->pmatrix_index] + n*states*states_padded;
    t *= locus->rates[n];

    if (t < 1e-100)
    {
      pmat[0]  = 1;
      pmat[1]  = 0;
      pmat[2]  = 0;
      pmat[3]  = 0;

      pmat[4]  = 0;
      pmat[5]  = 1;
      pmat[6]  = 0;
      pmat[7]  = 0;

      pmat[8]  = 0;
      pmat[9]  = 0;
      pmat[10] = 1;
      pmat[11] = 0;

      pmat[12] = 0;
      pmat[13] = 0;
      pmat[14] = 0;
      pmat[15] = 1;
    }
    else
    {
      double a =  (1 + 3*exp(-4*t/3) ) / 4;
      double b = (1 - a) / 3;

      pmat[0]  = a;
      pmat[1]  = b;
      pmat[2]  = b;
      pmat[3]  = b;

      pmat[4]  = b;
      pmat[5]  = a;
      pmat[6]  = b;
      pmat[7]  = b;

      pmat[8]  = b;
      pmat[9]  = b;
      pmat[10] = a;
      pmat[11] = b;

      pmat[12] = b;
      pmat[13] = b;
      pmat[14] = b;
      pmat[15] = a;
    }
  }

  if (!(root->left)) return;

  locus_update_all_matrices_jc69_recursive(locus,
                                           gtree,
                                           root->left,
                                           stree,
                                           msa_index);
  locus_update_all_matrices_jc69_recursive(locus,
                                           gtree,
                                           root->right,
                                           stree,
                                           msa_index);
}

static void locus_update_all_matrices_jc69(locus_t * locus,
                                           gtree_t * gtree,
                                           stree_t * stree,
                                           long msa_index)
{
  locus_update_all_matrices_jc69_recursive(locus,
                                           gtree,
                                           gtree->root->left,
                                           stree,
                                           msa_index);
  locus_update_all_matrices_jc69_recursive(locus,
                                           gtree,
                                           gtree->root->right,
                                           stree,
                                           msa_index);
}

void locus_update_all_matrices(locus_t * locus,
                               gtree_t * gtree,
                               stree_t * stree,
                               long msa_index)
{
  if (locus->dtype == BPP_DATA_DNA && locus->model == BPP_DNA_MODEL_JC69)
  {
    locus_update_all_matrices_jc69(locus,gtree,stree,msa_index);
    return;
  }

  locus_update_all_matrices_generic(locus, gtree, stree, msa_index);

}

static void locus_update_matrices_jc69(locus_t * locus,
                                       gtree_t * gtree,
                                       gnode_t ** traversal,
                                       stree_t * stree,
                                       long msa_index,
                                       unsigned int count)
{
  unsigned int i,n;
  double t;
  double * pmat;
  gnode_t * node;

  unsigned int states = locus->states;
  unsigned int states_padded = locus->states_padded;

  if (!opt_usedata) return;

  for (i = 0; i < count; ++i)
  {
    node = traversal[i];

    
    if (opt_clock == BPP_CLOCK_GLOBAL)
    {
      /* strict clock */
      t = node->length = (node->parent->time - node->time)*gtree->rate_mui;
    }
    else
    {
      /* relaxed clock */
      t = node->length = update_branchlength_relaxed_clock(stree,node,msa_index);
    }

    for (n = 0; n < locus->rate_cats; ++n)
    {
      pmat = locus->pmatrix[node->pmatrix_index] + n*states*states_padded;
      t *= locus->rates[n];

      if (t < 1e-100)
      {
        pmat[0]  = 1;
        pmat[1]  = 0;
        pmat[2]  = 0;
        pmat[3]  = 0;

        pmat[4]  = 0;
        pmat[5]  = 1;
        pmat[6]  = 0;
        pmat[7]  = 0;

        pmat[8]  = 0;
        pmat[9]  = 0;
        pmat[10] = 1;
        pmat[11] = 0;

        pmat[12] = 0;
        pmat[13] = 0;
        pmat[14] = 0;
        pmat[15] = 1;
      }
      else
      {
        double a =  (1 + 3*exp(-4*t/3) ) / 4;
        double b = (1 - a) / 3;

        pmat[0]  = a;
        pmat[1]  = b;
        pmat[2]  = b;
        pmat[3]  = b;

        pmat[4]  = b;
        pmat[5]  = a;
        pmat[6]  = b;
        pmat[7]  = b;

        pmat[8]  = b;
        pmat[9]  = b;
        pmat[10] = a;
        pmat[11] = b;

        pmat[12] = b;
        pmat[13] = b;
        pmat[14] = b;
        pmat[15] = a;
      }
    }
  }
}

void locus_update_matrices(locus_t * locus,
                           gtree_t * gtree,
                           gnode_t ** traversal,
                           stree_t * stree,
                           long msa_index,
                           unsigned int count)
{
  if (!opt_usedata) return;
  if (locus->dtype == BPP_DATA_DNA && locus->model == BPP_DNA_MODEL_JC69)
  {
    locus_update_matrices_jc69(locus,gtree,traversal,stree,msa_index,count);
    return;
  }

  unsigned int n;
  
  unsigned int * param_indices = locus->param_indices;

  for (n = 0; n < locus->rate_cats; ++n)
  {
    unsigned int param_index = param_indices[n];
    if (!locus->eigen_decomp_valid[0])
    {
      pll_update_eigen(locus->eigenvecs[param_index],
                       locus->inv_eigenvecs[param_index],
                       locus->eigenvals[param_index],
                       locus->frequencies[param_index],
                       locus->subst_params[param_index],
                       locus->states,
                       locus->states_padded);
      locus->eigen_decomp_valid[param_index] = 1;
    }
  }

  bpp_core_update_pmatrix(locus,gtree,traversal,stree,msa_index,count);
}


static void locus_update_all_partials_recursive(locus_t * locus, gnode_t * root)
{
  unsigned int * scaler;
  unsigned int * lscaler;
  unsigned int * rscaler;
  gnode_t * lnode;
  gnode_t * rnode;


  if (!(root->left)) return;

  locus_update_all_partials_recursive(locus,root->left);
  locus_update_all_partials_recursive(locus,root->right);

  lnode = root->left;
  rnode = root->right;

  /* check if we use scalers */
  scaler = (root->scaler_index == PLL_SCALE_BUFFER_NONE) ?
             NULL : locus->scale_buffer[root->scaler_index];

  lscaler = (lnode->scaler_index == PLL_SCALE_BUFFER_NONE) ?
              NULL : locus->scale_buffer[lnode->scaler_index];

  rscaler = (rnode->scaler_index == PLL_SCALE_BUFFER_NONE) ?
              NULL : locus->scale_buffer[rnode->scaler_index];

  pll_core_update_partial_ii(locus->states,
                             locus->sites,
                             locus->rate_cats,
                             locus->clv[root->clv_index],
                             scaler,
                             locus->clv[lnode->clv_index],
                             locus->clv[rnode->clv_index],
                             locus->pmatrix[lnode->pmatrix_index],
                             locus->pmatrix[rnode->pmatrix_index],
                             lscaler,
                             rscaler,
                             locus->attributes);
}

void locus_update_all_partials(locus_t * locus, gtree_t * gtree)
{
  if (!opt_usedata) return;

  locus_update_all_partials_recursive(locus,gtree->root);
}

void locus_update_partials(locus_t * locus, gnode_t ** traversal, unsigned int count)
{
  unsigned int i;
  unsigned int * scaler;
  unsigned int * lscaler;
  unsigned int * rscaler;
  gnode_t * node;
  gnode_t * lnode;
  gnode_t * rnode;

  if (!opt_usedata) return;

  for (i = 0; i < count; ++i)
  {
    node  = traversal[i];
    lnode = traversal[i]->left;
    rnode = traversal[i]->right;

    /* check if we use scalers */
    scaler = (node->scaler_index == PLL_SCALE_BUFFER_NONE) ?
               NULL : locus->scale_buffer[node->scaler_index];

    lscaler = (lnode->scaler_index == PLL_SCALE_BUFFER_NONE) ?
                NULL : locus->scale_buffer[lnode->scaler_index];

    rscaler = (rnode->scaler_index == PLL_SCALE_BUFFER_NONE) ?
                NULL : locus->scale_buffer[rnode->scaler_index];

    pll_core_update_partial_ii(locus->states,
                               locus->sites,
                               locus->rate_cats,
                               locus->clv[node->clv_index],
                               scaler,
                               locus->clv[lnode->clv_index],
                               locus->clv[rnode->clv_index],
                               locus->pmatrix[lnode->pmatrix_index],
                               locus->pmatrix[rnode->pmatrix_index],
                               lscaler,
                               rscaler,
                               locus->attributes);
  }
}

double locus_root_loglikelihood(locus_t * locus,
                                gnode_t * root,
                                const unsigned int * freqs_indices,
                                double * persite_lnl)
{
  double logl;
  unsigned int * scaler;

  if (!opt_usedata) return 0;

  scaler = (root->scaler_index == PLL_SCALE_BUFFER_NONE) ?
             NULL : locus->scale_buffer[root->scaler_index];

  if (locus->diploid)
  {
    pll_core_root_likelihood_vector(locus->states,
                                    locus->sites,
                                    locus->rate_cats,
                                    locus->clv[root->clv_index],
                                    scaler,
                                    locus->frequencies,
                                    locus->rate_weights,
                                    locus->pattern_weights,
                                    freqs_indices,
                                    locus->likelihood_vector,
                                    locus->attributes);
    
    long i,j,k=0;
    logl = 0;

    for (i = 0; i < locus->unphased_length; ++i)
    {
      double meanl = 0;

      for (j = 0; j < locus->diploid_resolution_count[i]; ++j)
        meanl += locus->likelihood_vector[locus->diploid_mapping[k++]];

      meanl /= locus->diploid_resolution_count[i];

      logl += log(meanl) * locus->pattern_weights[i];
    }
  }
  else
  {
    logl = pll_core_root_loglikelihood(locus->states,
                                       locus->sites,
                                       locus->rate_cats,
                                       locus->clv[root->clv_index],
                                       scaler,
                                       locus->frequencies,
                                       locus->rate_weights,
                                       locus->pattern_weights,
                                       freqs_indices,
                                       persite_lnl,
                                       locus->attributes);
  }
  return opt_bfbeta * logl;
}

#if 0
static long propose_freqs(stree_t * stree,
                          locus_t * locus,
                          gtree_t * gtree,
                          long msa_index,
                          long thread_index)
{
  unsigned int i,j,k,m,n;
  long accepted = 0;
  double lnacceptance;
  double logl;
  double old_freqj, old_freqk;
  double sum;
  double x,y;
  unsigned int * param_indices = locus->param_indices;
  gnode_t ** gt_nodes;

  /* allocate temporary space for gene tree traversal */
  gt_nodes = (gnode_t **)xmalloc((gtree->tip_count+gtree->inner_count) *
                                 sizeof(gnode_t *));
    
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    double * freqs = locus->frequencies[param_indices[i]];

    /* select two base frequencies j and k at random and save them */
    j = (unsigned int)(legacy_rndu(thread_index) * locus->states);
    k = (unsigned int)(legacy_rndu(thread_index) * (locus->states-1));
    if (k >= j)
      ++k;
    old_freqj = freqs[j];
    old_freqk = freqs[k];

    sum = freqs[j] + freqs[k];

    /* compute ratio */
    x = freqs[j] / sum;

    /* min/max bounds for proposed value */
    double minv = PLL_MISC_EPSILON / sum;
    double maxv = 1 - minv;

    /* propose new ratio */
    y = x + opt_finetune_freqs*legacy_rnd_symmetrical(thread_index);
    y = reflect(y,minv,maxv,thread_index);

    /* set new proposed frequencies */
    freqs[j] = y*sum;
    freqs[k] = sum - freqs[j];

    #if 0
    if (fabs(locus->frequencies[0][0] + locus->frequencies[0][1] + locus->frequencies[0][2] + locus->frequencies[0][3]-1) >= 1e-10)
    {
      printf ("%f %f %f %f = %f\n",
              locus->frequencies[0][0], locus->frequencies[0][1], locus->frequencies[0][2], locus->frequencies[0][3], 
              locus->frequencies[0][0]+ locus->frequencies[0][1]+ locus->frequencies[0][2]+ locus->frequencies[0][3]);
      assert(0);
    }
    #endif
    
    /* swap pmatrix indices to new buffers, and update pmatrices */
    for (m=0,n=0; m < gtree->tip_count+gtree->inner_count; ++m)
    {
      gnode_t * p = gtree->nodes[m];
      if (p->parent)
      {
        SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
        gt_nodes[n++] = p;
      }
    }
    locus->eigen_decomp_valid[i] = 0;
    locus_update_matrices(locus,gtree,gt_nodes,stree,msa_index,n);

    /* get postorder traversal of inner nodes, swap CLV indidces to point to new
       buffer, and update partials */
    gtree_all_partials(gtree->root,gt_nodes,&n);
    for (m = 0; m < n; ++m)
    {
      gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                              gt_nodes[m]->clv_index);
      if (opt_scaling)
        gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                      gt_nodes[m]->scaler_index);
    }
    locus_update_partials(locus,gt_nodes,n);

    /* compute log-likelihood */
    logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

    lnacceptance = logl - gtree->logl;

    #if 0
    /* This is code for dubugging purposes to ensure that we obtain the prior
       when running the program without data */

    double alpha[4] = {10,25,30,35};
    lnacceptance += (alpha[j]-1.0)*log(freqs[j]/old_freqj);
    lnacceptance += (alpha[k]-1.0)*log(freqs[k]/old_freqk);
    #endif

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
    {
      /* accepted */
      accepted = 1;
      gtree->logl = logl;
    }
    else
    {
      /* rejected */

      /* revert pmatrices */
      for (m = 0; m < gtree->tip_count + gtree->inner_count; ++m)
      {
        gnode_t * p = gtree->nodes[m];
        if (p->parent)
          SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
      }

      /* revert CLV */
      for (m = 0; m < n; ++m)
      {
        gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                gt_nodes[m]->clv_index);
        if (opt_scaling)
          gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                        gt_nodes[m]->scaler_index);
      }

      /* revert old frequencies */
      freqs[j] = old_freqj;
      freqs[k] = old_freqk;

      /* update eigen decomposition */
      if (opt_usedata)
        pll_update_eigen(locus->eigenvecs[param_indices[i]],
                         locus->inv_eigenvecs[param_indices[i]],
                         locus->eigenvals[param_indices[i]],
                         locus->frequencies[param_indices[i]],
                         locus->subst_params[param_indices[i]],
                         locus->states,
                         locus->states_padded);
    }
  }
  free(gt_nodes);

  return accepted;
}
#endif

static void propose_freqs(stree_t * stree,
                          locus_t * locus,
                          gtree_t * gtree,
                          long msa_index,
                          long thread_index,
                          long * p_accepted,
                          long * p_candidates)
{
  unsigned int i,j,m,n;
  unsigned int ref = locus->states - 1;
  long accepted = 0;
  long candidates = 0;
  double lnacceptance;
  double logl;
  double old_freq_j, old_freq_ref;
  double old_logfreq_j, new_logfreq_j;
  double sum;
  unsigned int * param_indices = locus->param_indices;
  gnode_t ** gt_nodes;

  /* allocate temporary space for gene tree traversal */
  gt_nodes = (gnode_t **)xmalloc((gtree->tip_count+gtree->inner_count) *
                                 sizeof(gnode_t *));
    
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    double * freqs = locus->frequencies[param_indices[i]];

    for (j = 0; j < locus->states; ++j)
    {
      /* skip reference rate */
      if (j == ref) continue;

      ++candidates;

      /* set bounds for proposing new freq */
      sum = freqs[j] + freqs[ref];
      double minv = log(1e-5);
      double maxv = log(sum);

      old_freq_j   = freqs[j];
      old_freq_ref = freqs[ref];
      old_logfreq_j   = log(old_freq_j);

      /* propose new ratio */
      new_logfreq_j = log(old_freq_j) + opt_finetune_freqs *
                      legacy_rnd_symmetrical(thread_index);
      new_logfreq_j = reflect(new_logfreq_j,minv,maxv,thread_index);

      /* set new proposed rates */
      freqs[j]   = exp(new_logfreq_j);
      freqs[ref] = sum - freqs[j];

      /* swap pmatrix indices to new buffers, and update pmatrices */
      for (m=0,n=0; m < gtree->tip_count+gtree->inner_count; ++m)
      {
        gnode_t * p = gtree->nodes[m];
        if (p->parent)
        {
          SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
          gt_nodes[n++] = p;
        }
      }
      locus->eigen_decomp_valid[i] = 0;
      locus_update_matrices(locus,gtree,gt_nodes,stree,msa_index,n);

      /* get postorder traversal of inner nodes, swap CLV indidces to point to new
         buffer, and update partials */
      gtree_all_partials(gtree->root,gt_nodes,&n);
      for (m = 0; m < n; ++m)
      {
        gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                gt_nodes[m]->clv_index);
        if (opt_scaling)
          gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                        gt_nodes[m]->scaler_index);
      }
      locus_update_partials(locus,gt_nodes,n);

      /* compute log-likelihood */
      logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

      lnacceptance = new_logfreq_j - old_logfreq_j +
                     logl - gtree->logl;

      if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
      {
        /* accepted */
        ++accepted;
        gtree->logl = logl;
      }
      else
      {
        /* rejected */

        /* revert pmatrices */
        for (m = 0; m < gtree->tip_count + gtree->inner_count; ++m)
        {
          gnode_t * p = gtree->nodes[m];
          if (p->parent)
            SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
        }

        /* revert CLV */
        for (m = 0; m < n; ++m)
        {
          gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                  gt_nodes[m]->clv_index);
          if (opt_scaling)
            gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                          gt_nodes[m]->scaler_index);
        }

        /* revert old frequencies */
        freqs[j] = old_freq_j;
        freqs[ref] = old_freq_ref;

        /* update eigen decomposition */
        if (opt_usedata)
          pll_update_eigen(locus->eigenvecs[param_indices[i]],
                           locus->inv_eigenvecs[param_indices[i]],
                           locus->eigenvals[param_indices[i]],
                           locus->frequencies[param_indices[i]],
                           locus->subst_params[param_indices[i]],
                           locus->states,
                           locus->states_padded);
      }
    }
  }
  free(gt_nodes);

  *p_accepted = accepted;
  *p_candidates = candidates;
}

double locus_propose_freqs_serial(stree_t * stree, locus_t ** locus, gtree_t ** gtree)
{
  long i;
  long accepted = 0;
  long candidates = 0;
  long loc_acc = 0;
  long loc_cand = 0;
  long thread_index = 0;

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (locus[i]->dtype == BPP_DATA_DNA && locus[i]->model != BPP_DNA_MODEL_JC69)
    {
      propose_freqs(stree,locus[i],gtree[i],i,thread_index,&loc_acc,&loc_cand);

      accepted += loc_acc;
      candidates += loc_cand;
    }
  }

  if (!candidates)
    return 0;

  return (accepted/(double)candidates);
}

void locus_propose_freqs_parallel(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree,
                                  long locus_start,
                                  long locus_count,
                                  long thread_index,
                                  long * p_proposal_count,
                                  long * p_accepted)
{
  long i;
  long accepted = 0;
  long candidates = 0;
  long loc_acc = 0;
  long loc_cand = 0;

  assert(locus_start >= 0);
  assert(locus_count > 0);

  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    if (locus[i]->dtype == BPP_DATA_DNA && locus[i]->model != BPP_DNA_MODEL_JC69)
    {
      propose_freqs(stree,locus[i],gtree[i],i,thread_index,&loc_acc,&loc_cand);
      candidates += loc_cand;
      accepted += loc_acc;
    }
  }

  *p_proposal_count = candidates;
  *p_accepted = accepted;
}

#if 0
static long propose_qrates(stree_t * stree,
                           locus_t * locus,
                           gtree_t * gtree,
                           long msa_index,
                           long thread_index)
{
  unsigned int i,j,k,m,n;
  long accepted = 0;
  long rates_count = 0;
  double lnacceptance;
  double logl;
  double old_ratej, old_ratek;
  double sum;
  double x,y;
  unsigned int * param_indices = locus->param_indices;
  gnode_t ** gt_nodes;

  /* TODO: Implement amino acids */
  assert(locus->dtype == BPP_DATA_DNA);
  assert(locus->states == 4);

  switch (locus->model)
  {
    case BPP_DNA_MODEL_K80:
      rates_count = 1;
      assert(0);
      break;

    case BPP_DNA_MODEL_F81:
      rates_count = 0;
      assert(0);
      break;

    case BPP_DNA_MODEL_HKY:
    case BPP_DNA_MODEL_T92:
    case BPP_DNA_MODEL_F84:
      rates_count = 1;
      assert(0);
      break;

    case BPP_DNA_MODEL_TN93:
      rates_count = 2;
      break;

    case BPP_DNA_MODEL_GTR:
      rates_count = 6;
      break;

    default:
      assert(0);
  }


  /* allocate temporary space for gene tree traversal */
  gt_nodes = (gnode_t **)xmalloc((gtree->tip_count+gtree->inner_count) *
                                 sizeof(gnode_t *));
    
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    double * rates = locus->subst_params[param_indices[i]];

    /* select two rates j and k at random and save them */
    j = (unsigned int)(legacy_rndu(thread_index) * rates_count);
    k = (unsigned int)(legacy_rndu(thread_index) * (rates_count-1));
    if (k >= j)
      ++k;
    old_ratej = rates[j];
    old_ratek = rates[k];

    sum = rates[j] + rates[k];

    /* compute ratio */
    x = rates[j] / sum;

    /* min/max bounds for proposed value */
    double minv = PLL_MISC_EPSILON / sum;
    double maxv = 1 - minv;

    /* propose new ratio */
    y = x + opt_finetune_qrates*legacy_rnd_symmetrical(thread_index);
    y = reflect(y,minv,maxv,thread_index);

    /* set new proposed rates */
    rates[j] = y*sum;
    rates[k] = sum - rates[j];

    #if 0
    if (fabs(locus->subst_params[0][0] + locus->subst_params[0][1] + locus->subst_params[0][2] +
             locus->subst_params[0][3] + locus->subst_params[0][4] + locus->subst_params[0][5] -6) >= 1e-10)
    {
      printf ("%f %f %f %f %f %f = %f\n",
              locus->subst_params[0][0], locus->subst_params[0][1], locus->subst_params[0][2],
              locus->subst_params[0][3], locus->subst_params[0][4], locus->subst_params[0][5],
              locus->subst_params[0][0]+ locus->subst_params[0][1]+ locus->subst_params[0][2]+
              locus->subst_params[0][3]+ locus->subst_params[0][4]+ locus->subst_params[0][5]);
      assert(0);
    }
    #endif
    
    /* swap pmatrix indices to new buffers, and update pmatrices */
    for (m=0,n=0; m < gtree->tip_count+gtree->inner_count; ++m)
    {
      gnode_t * p = gtree->nodes[m];
      if (p->parent)
      {
        SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
        gt_nodes[n++] = p;
      }
    }
    locus->eigen_decomp_valid[i] = 0;
    locus_update_matrices(locus,gtree,gt_nodes,stree,msa_index,n);

    /* get postorder traversal of inner nodes, swap CLV indidces to point to new
       buffer, and update partials */
    gtree_all_partials(gtree->root,gt_nodes,&n);
    for (m = 0; m < n; ++m)
    {
      gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                              gt_nodes[m]->clv_index);
      if (opt_scaling)
        gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                      gt_nodes[m]->scaler_index);
    }
    locus_update_partials(locus,gt_nodes,n);

    /* compute log-likelihood */
    logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

    lnacceptance = logl - gtree->logl;

    #if 0
    /* This is code for dubugging purposes to ensure that we obtain the prior
       when running the program without data */

    double alpha[6] = {1,2,1,1,2,1};
    lnacceptance += (alpha[j]-1.0)*log(rates[j]/old_ratej);
    lnacceptance += (alpha[k]-1.0)*log(rates[k]/old_ratek);
    #endif

    if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
    {
      /* accepted */
      accepted = 1;
      gtree->logl = logl;
    }
    else
    {
      /* rejected */

      /* revert pmatrices */
      for (m = 0; m < gtree->tip_count + gtree->inner_count; ++m)
      {
        gnode_t * p = gtree->nodes[m];
        if (p->parent)
          SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
      }

      /* revert CLV */
      for (m = 0; m < n; ++m)
      {
        gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                gt_nodes[m]->clv_index);
        if (opt_scaling)
          gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                        gt_nodes[m]->scaler_index);
      }

      /* revert old rates */
      rates[j] = old_ratej;
      rates[k] = old_ratek;

      /* update eigen decomposition */
      if (opt_usedata)
        pll_update_eigen(locus->eigenvecs[param_indices[i]],
                         locus->inv_eigenvecs[param_indices[i]],
                         locus->eigenvals[param_indices[i]],
                         locus->frequencies[param_indices[i]],
                         locus->subst_params[param_indices[i]],
                         locus->states,
                         locus->states_padded);
    }
  }
  free(gt_nodes);
  return accepted;
}
#endif

static void propose_qrates(stree_t * stree,
                             locus_t * locus,
                             gtree_t * gtree,
                             long msa_index,
                             long thread_index,
                             long * p_accepted,
                             long * p_candidates)
{
  unsigned int i,j,m,n;
  unsigned int ref = 0;
  long accepted = 0;
  long candidates = 0;
  long qrates_count = 0;
  double lnacceptance;
  double logl;
  double old_rate_j, old_rate_ref;
  double old_lograte_j, new_lograte_j;
  double sum;
  unsigned int * param_indices = locus->param_indices;
  gnode_t ** gt_nodes;
#if 0
  double gtr_alpha[6] = { 1,1,1,1,1,1 };
#else
  double gtr_alpha[6] = { 1,2,1,1,2,1 };
#endif

  /* TODO: Implement amino acids */
  assert(locus->dtype == BPP_DATA_DNA);
  assert(locus->states == 4);

  switch (locus->model)
  {
    case BPP_DNA_MODEL_K80:
      qrates_count = 1;
      assert(0);
      break;

    case BPP_DNA_MODEL_F81:
      qrates_count = 0;
      assert(0);
      break;

    case BPP_DNA_MODEL_HKY:
    case BPP_DNA_MODEL_T92:
    case BPP_DNA_MODEL_F84:
      qrates_count = 1;
      assert(0);
      break;

    case BPP_DNA_MODEL_TN93:
      qrates_count = 2;
      break;

    case BPP_DNA_MODEL_GTR:
      qrates_count = 6;
      ref = 1;          /* A->G */
      break;

    default:
      assert(0);
  }


  /* allocate temporary space for gene tree traversal */
  /* TODO: Allocate this memory only once inside the gene tree so that allocas/
     deallocs are no longer necessary. This will also simplify gtree.c */
  gt_nodes = (gnode_t **)xmalloc((gtree->tip_count+gtree->inner_count) *
                                 sizeof(gnode_t *));
    
  for (i = 0; i < locus->rate_matrices; ++i)
  {
    double * qrates = locus->subst_params[param_indices[i]];

    for (j = 0; j < qrates_count; ++j)
    {
      /* skip reference rate */
      if (j == ref) continue;

      ++candidates;

      /* set bounds for proposing new rate */
      sum = qrates[j] + qrates[ref];
      double minv = log(1e-5);
      double maxv = log(sum);

      old_rate_j   = qrates[j];
      old_rate_ref = qrates[ref];
      old_lograte_j   = log(old_rate_j);
      
      /* propose new ratio */
      new_lograte_j = old_lograte_j + opt_finetune_qrates *
                      legacy_rnd_symmetrical(thread_index);
      new_lograte_j = reflect(new_lograte_j,minv,maxv,thread_index);

      /* set new proposed rates */
      qrates[j]   = exp(new_lograte_j);
      qrates[ref] = sum - qrates[j];

      /* swap pmatrix indices to new buffers, and update pmatrices */
      for (m=0,n=0; m < gtree->tip_count+gtree->inner_count; ++m)
      {
        gnode_t * p = gtree->nodes[m];
        if (p->parent)
        {
          SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
          gt_nodes[n++] = p;
        }
      }
      locus->eigen_decomp_valid[i] = 0;
      locus_update_matrices(locus,gtree,gt_nodes,stree,msa_index,n);

      /* get postorder traversal of inner nodes, swap CLV indidces to point to new
         buffer, and update partials */
      gtree_all_partials(gtree->root,gt_nodes,&n);
      for (m = 0; m < n; ++m)
      {
        gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                gt_nodes[m]->clv_index);
        if (opt_scaling)
          gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                        gt_nodes[m]->scaler_index);
      }
      locus_update_partials(locus,gt_nodes,n);

      /* compute log-likelihood */
      logl = locus_root_loglikelihood(locus,gtree->root,param_indices,NULL);

      lnacceptance = new_lograte_j - old_lograte_j +
                     logl - gtree->logl;

      /* This is code for dubugging purposes to ensure that we obtain the prior when running the program without data */
      if (gtr_alpha[j] - 1)
        lnacceptance += (gtr_alpha[j] - 1.0)*(new_lograte_j - old_lograte_j);
      if (gtr_alpha[ref] - 1)
        lnacceptance += (gtr_alpha[ref] - 1.0)*log(qrates[ref] / old_rate_ref);

      if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
      {
        /* accepted */
        ++accepted;
        gtree->logl = logl;
      }
      else
      {
        /* rejected */

        /* revert pmatrices */
        for (m = 0; m < gtree->tip_count + gtree->inner_count; ++m)
        {
          gnode_t * p = gtree->nodes[m];
          if (p->parent)
            SWAP_PMAT_INDEX(gtree->edge_count, p->pmatrix_index);
        }

        /* revert CLV */
        for (m = 0; m < n; ++m)
        {
          gt_nodes[m]->clv_index = SWAP_CLV_INDEX(gtree->tip_count,
                                                  gt_nodes[m]->clv_index);
          if (opt_scaling)
            gt_nodes[m]->scaler_index = SWAP_SCALER_INDEX(gtree->tip_count,
                                                          gt_nodes[m]->scaler_index);
        }

        /* revert old qrates */
        qrates[j]   = old_rate_j;
        qrates[ref] = old_rate_ref;

        /* update eigen decomposition */
        if (opt_usedata)
          pll_update_eigen(locus->eigenvecs[param_indices[i]],
                           locus->inv_eigenvecs[param_indices[i]],
                           locus->eigenvals[param_indices[i]],
                           locus->frequencies[param_indices[i]],
                           locus->subst_params[param_indices[i]],
                           locus->states,
                           locus->states_padded);
      }
    }
  }
  free(gt_nodes);

  *p_accepted = accepted;
  *p_candidates = candidates;
}

double locus_propose_qrates_serial(stree_t * stree,
                                 locus_t ** locus,
                                 gtree_t ** gtree)
{
  long i;
  long accepted = 0;
  long candidates = 0;
  long loc_acc = 0;
  long loc_cand = 0;
  long thread_index = 0;

  for (i = 0; i < opt_locus_count; ++i)
  {
    if (locus[i]->dtype == BPP_DATA_DNA && locus[i]->model != BPP_DNA_MODEL_JC69)
    {
      propose_qrates(stree,locus[i],gtree[i],i,thread_index,&loc_acc,&loc_cand);

      accepted += loc_acc;
      candidates += loc_cand;
    }
  }

  if (!candidates)
    return 0;

  return (accepted/(double)candidates);
}

void locus_propose_qrates_parallel(stree_t * stree,
                                  locus_t ** locus,
                                  gtree_t ** gtree,
                                  long locus_start,
                                  long locus_count,
                                  long thread_index,
                                  long * p_proposal_count,
                                  long * p_accepted)
{
  long i;
  long accepted = 0;
  long candidates = 0;
  long loc_acc = 0;
  long loc_cand = 0;

  assert(locus_start >= 0);
  assert(locus_count > 0);

  for (i = locus_start; i < locus_start+locus_count; ++i)
  {
    if (locus[i]->dtype == BPP_DATA_DNA && locus[i]->model != BPP_DNA_MODEL_JC69)
    {
      propose_qrates(stree,locus[i],gtree[i],i,thread_index,&loc_acc,&loc_cand);
      candidates += loc_cand;
      accepted += loc_acc;
    }
  }

  *p_proposal_count = candidates;
  *p_accepted = accepted;
}

