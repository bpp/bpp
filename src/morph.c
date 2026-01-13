//
//  Created by Chi Zhang on 2023/12/5.
//

#include "bpp.h"

#define LABEL_LEN 99

// #define DEBUG_Morph_BM     1

/* get a non-blank character from file */
static int get_nb_char(FILE *fp)
{
  int c;
  while (isspace(c=fgetc(fp)));
  return c;
}

static int parse_header(FILE * fp, int * nrow, int * ncol, int * type)
{
  int t;
  
  if (fscanf(fp, "%d %d", nrow, ncol) != 2) {
    fprintf(stderr, "Expecting two numbers\n");
    return 0;
  }
  else if (*nrow <= 0) {
    fprintf(stderr, "Number of species must be >0(%d)\n", *nrow);
    return 0;
  }
  else if (*ncol <= 0) {
    fprintf(stderr, "Number of traits must be >0 (%d)\n", *ncol);
    return 0;
  }

  t = get_nb_char(fp);
  if (t == 'C' || t == 'c')
  {
    *type = BPP_DATA_CONT;
  }
  else if (t == 'D' || t == 'd')
  {
    *type = BPP_DATA_DISC;
  }
  else
  {
    fprintf(stderr, "Unrecognized data type (%c)\n", t);
    return 0;
  }

  return 1;
}

static int parse_label(FILE * fp, char * name, int len)
{
  int c, j;
  
  c = get_nb_char(fp);
  if (c == EOF)
    return 0;
  
  for (j = 0; j < len; ++j)
  {
    name[j] = c;
    c = fgetc(fp);
    if (c == EOF || isspace(c)) break;
  }
  
  if (j < len)
    name[j+1] = '\0';
  else
    name[len] = '\0';

  return 1;
}

static int parse_value_c(FILE * fp, double * trait, int len)
{
  int c, j;
  
  for (j = 0; j < len; ++j)
  {
    c = get_nb_char(fp);
    
    if (c == EOF)
    {
      return 0;
    }
    else if (c == '-')
    {
      c = fgetc(fp);
      if (isspace(c) || c == EOF)
        trait[j] = NAN;    // inapplicable value
      else
      {
        ungetc(c, fp);
        if(fscanf(fp, "%lf", &trait[j]) != 1)
          return 0;
        trait[j] = -trait[j];  // negative value
      }
    }
    else if (c == '?')
    {
      trait[j] = INFINITY;      // missing value
    }
    else if (c == 'N' || c == 'n')
    {
      c = fgetc(fp);
      if(c != 'A' && c != 'a')
        return 0;
      else
      {
        c = fgetc(fp);
        if (isspace(c) || c == EOF)
          trait[j] = INFINITY;  // missing value
        else if (c == 'N' || c == 'n')
          trait[j] = NAN;  // inapplicable value
        else
          return 0;
      }
    }
    else
    {
      ungetc(c, fp);           // positive value
      if(fscanf(fp, "%lf", &trait[j]) != 1)
        return 0;
    }
  }
  assert(j == len);

  return 1;
}

static int state_bin(int x)
{
  int std_bin[] = {1,2,4,8,16,32,64,128,256,512};
  
  if (x >= 0 && x <= 9)
    return std_bin[x];
  else
  {
    fprintf(stderr, "Unsupported trait value (%d)\n", x);
    return 0;
  }
}

static int parse_value_d(FILE * fp, int * std, int len)
{
  int c, j;
  
  /* store the states in binary to get ready for bit operations */
  for (j = 0; j < len; ++j)
  {
    c = get_nb_char(fp);
    
    if (c == EOF)
      return 0;
    
    if (isdigit(c))  // 0-9
    {
      std[j] = state_bin(c-'0');
    }
    else if (c == '?' || c == '-')  // treat gap as missing
    {
      std[j] = 1023;
    }
    else if (c == '{' || c == '(')  // ambiguity
    {
      std[j] = 0;
      while ((c = get_nb_char(fp)) != '}' && c != ')')
        if (isdigit(c))
          std[j] |= state_bin(c-'0');
        
      if (std[j] == 0) {
        fprintf(stderr, "Misspecified ambiguity state at %d\n", j+1);
        return 0;
      }
    }
    else
    {
      fprintf(stderr, "Unrecognized trait value (%c) at %d\n", c, j+1);
      return 0;
    }
  }
  
  return 1;
}

static void parse_comment(FILE * fp)
{
  int c;
  
  while ((c=get_nb_char(fp)) == '#')
    while ((c=fgetc(fp)) != '\n'); //eat the line

  ungetc(c, fp);
}

static morph_t * parse_trait_part(FILE * fp)
{
  int i, j;
  
  morph_t * morph = (morph_t *)xmalloc(sizeof(morph_t));
  
  /* skip comments, assuming they are at the beginning of the trait block */
  parse_comment(fp);
  
  /* read header */
  if (!parse_header(fp, &(morph->ntaxa), &(morph->length), &(morph->dtype)))
  {
    fprintf(stderr, "Error in header\n");
    return NULL;
  }
  
  /* allocate space */
  morph->label = (char **)xmalloc((morph->ntaxa)*sizeof(char *));
  if (morph->dtype == BPP_DATA_CONT)
  {
    morph->conti = (double **)xmalloc((morph->ntaxa)*sizeof(double *));
    for (i = 0; i < morph->ntaxa; ++i)
    {
      morph->conti[i] = (double *)xmalloc((morph->length)*sizeof(double));
      morph->label[i] = (char *)xmalloc((LABEL_LEN+1)*sizeof(char));
    }
  }
  else // morph->dtype == BPP_DATA_DISC
  {
    morph->discr = (int **)xmalloc((morph->ntaxa)*sizeof(int *));
    for (i = 0; i < morph->ntaxa; ++i)
    {
      morph->discr[i] = (int *)xmalloc((morph->length)*sizeof(int));
      morph->label[i] = (char *)xmalloc((LABEL_LEN+1)*sizeof(char));
    }
  }
  
  /* read morphological traits of each species,
     assuming each line has a label followed by m->length numbers.
     trait matrix has dimension m->ntaxa * m->length */
  for (i = 0; i < morph->ntaxa; ++i)
  {
    /* read the label (species name) */
    if (!parse_label(fp, morph->label[i], LABEL_LEN+1))
    {
      fprintf(stderr, "Failed to read label of species %d\n", i+1);
      return NULL;
    }
    /* read the trait values of this species */
    if (morph->dtype == BPP_DATA_CONT &&
        !parse_value_c(fp, morph->conti[i], morph->length))
    {
      fprintf(stderr, "Failed to read traits of species %s\n", morph->label[i]);
      return NULL;
    }
    else if (morph->dtype == BPP_DATA_DISC &&
             !parse_value_d(fp, morph->discr[i], morph->length))
    {
      fprintf(stderr, "Failed to read traits of species %s\n", morph->label[i]);
      return NULL;
    }
  }
  
  /* TODO: also read in correlation matrix */

#ifdef DEBUG_Morph_Matrix
  printf("\n%d %d  %d\n", morph->ntaxa, morph->length, morph->dtype);
  for (i = 0; i < morph->ntaxa; ++i) {
    printf("%s\t", morph->label[i]);
    if (morph->dtype == BPP_DATA_CONT)
      for (j = 0; j < morph->length; ++j)
        printf("%+.3lf\t", morph->conti[i][j]);
    else
      for (j = 0; j < morph->length; ++j)
        printf("%4d ", morph->discr[i][j]);
    printf("\n");
  }
#endif
  
  return morph;
}

morph_t ** parse_traitfile(const char * traitfile, long * count)
{
  int c, i, m_slotalloc = 10, m_maxcount = 10;
  
  FILE * fp = xopen(traitfile,"r");

  morph_t ** pmorph = (morph_t **)xmalloc(m_maxcount*sizeof(morph_t *));

  *count = 0;

  while ((c=fgetc(fp)) != EOF)
  {
    if (isspace(c)) continue;
    ungetc(c,fp);
    
    if (*count == m_maxcount)
    {
      m_maxcount += m_slotalloc;
      morph_t ** temp = (morph_t **)xmalloc(m_maxcount*sizeof(morph_t *));
      memcpy(temp, pmorph, (*count)*sizeof(morph_t *));
      free(pmorph);
      pmorph = temp;
    }

    pmorph[*count] = parse_trait_part(fp);
    if (pmorph[*count] == NULL)
    {
      for (i = 0; i < *count; ++i)
        morph_destroy(pmorph[i]);
      free(pmorph);

      fatal("Error parsing trait partition %d", *count+1);
    }

    *count += 1;  // count the number of partitions
  }

  fclose(fp);

  return pmorph;
}

void morph_destroy(morph_t * morph)
{
  int i;

  if (morph->label)
  {
    for (i = 0; i < morph->ntaxa; ++i)
      if (morph->label[i])
        free(morph->label[i]);
    free(morph->label);
  }

  if (morph->conti)
  {
    for (i = 0; i < morph->ntaxa; ++i)
      if (morph->conti[i])
        free(morph->conti[i]);
    free(morph->conti);
  }

  if (morph->discr)
  {
    for (i = 0; i < morph->ntaxa; ++i)
      if (morph->discr[i])
        free(morph->discr[i]);
    free(morph->discr);
  }

  free(morph);
}

void trait_destroy(stree_t * stree)
{
  int n, i, j;
  snode_t * snode;
  trait_t * trait;

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    if (snode->trait)
    {
      for (n = 0; n < stree->trait_count; ++n)
      {
        trait = snode->trait[n];
        if (trait)
        {
          if (trait->state_d)
            free(trait->state_d);
          if (trait->condprob)
          {
            for (j = 0; j < stree->trait_dim[n] +54; ++j)
              free(trait->condprob[j]);
            free(trait->condprob);
          }
          if (trait->tranprob)
          {
            for (j = 0; j < 9; ++j)
              free(trait->tranprob[j]);
            free(trait->tranprob);
          }
          if (trait->state_m)
            free(trait->state_m);
          if (trait->contrast)
            free(trait->contrast);
          if (trait->active)
            free(trait->active);
          free(trait);
        }
      }
      free(snode->trait);
    }
  }
  
  if (stree->trait_dim)
    free(stree->trait_dim);
  if (stree->trait_type)
    free(stree->trait_type);
  if (stree->trait_missing)
    free(stree->trait_missing);

  if (stree->trait_nstate)
  {
    for (n = 0; n < stree->trait_count; ++n)
    {
      if (stree->trait_nstate[n])
        free(stree->trait_nstate[n]);
    }
    free(stree->trait_nstate);
  }
  if (stree->trait_v_pop)
    free(stree->trait_v_pop);
  if (stree->trait_ldetRs)
    free(stree->trait_ldetRs);

  if (stree->trait_logl)
    free(stree->trait_logl);
  if (stree->trait_old_logl)
    free(stree->trait_old_logl);
  if (stree->trait_logpr)
    free(stree->trait_logpr);
  if (stree->trait_old_logpr)
    free(stree->trait_old_logpr);
}

static void trait_update_ic(int idx, snode_t * snode, stree_t * stree)
{
  int j;
  double v_pop, v_k, v_k1, v_k2, *m_k1, *m_k2;

  if (snode->left && snode->right)  /* internal node */
  {
    trait_update_ic(idx, snode->left, stree);
    trait_update_ic(idx, snode->right, stree);
    
    if (snode->parent)
      v_k = (snode->parent->tau - snode->tau) * snode->trait[idx]->brate;
    else
      v_k = 0;
    v_k1 = snode->left->trait[idx]->brlen;
    v_k2 = snode->right->trait[idx]->brlen;
    snode->trait[idx]->brlen = v_k + v_k1*v_k2 / (v_k1+v_k2);  // v_k'
    
    m_k1 = snode->left->trait[idx]->state_m;
    m_k2 = snode->right->trait[idx]->state_m;
    for (j = 0; j < stree->trait_dim[idx]; ++j)
    {
      snode->trait[idx]->contrast[j] = m_k1[j] - m_k2[j];  // x_k and m_k'
      snode->trait[idx]->state_m[j] = (v_k2*m_k1[j]+v_k1*m_k2[j]) / (v_k1+v_k2);
    }
  }
  else  /* tip node */
  {
    v_k = (snode->parent->tau - snode->tau) * snode->trait[idx]->brate;
    /* the trait matrix has been standardized so that all characters have the
       same variance and the population noise has variance of v_pop (1.0).
       Alvarez-Carretero et al. 2019. p.970. */
    v_pop = 1.0;  //TODO: deal with stree->trait_v_pop[idx];
    assert(v_k > 0 || v_pop > 0);
    snode->trait[idx]->brlen = v_k + v_pop;
  }

#ifdef DEBUG_Morph_BM
  printf("%s\t v'=%lf\n", snode->label, snode->trait[idx]->brlen);
  for (j = 0; j < stree->trait_dim[idx]; ++j) {
    printf("m'=%lf ", snode->trait[idx]->state_m[j]);
    if (snode->left)
      printf("x=%lf", snode->trait[idx]->contrast[j]);
    printf("\n");
  }
  printf("\n");
#endif
}

static void trait_update_lmr(int idx, snode_t * snode, stree_t * stree)
{
  

}

static void trait_trprob_mk(double ** p, double v, int max_state)
{
  int k;
  
  /* Mk model, Lewis 2001 */
  for (k = 2; k <= max_state; ++k)
  {
    p[k-2][0] = 1.0/k + (k-1.0)/k * exp(-v * k/(k-1.0)); // no change
    p[k-2][1] = 1.0/k - 1.0/k * exp(-v * k/(k-1.0));     // change
  }
}

static void trait_update_cp(int idx, snode_t * snode, stree_t * stree)
{
  int h, j, k, a, x, y, z;
  int * nstate, nchar, max_state;
  double v, prob_l, prob_r, tr_prob;
  
  /* update the branch length */
  if (snode->parent)
    v = (snode->parent->tau - snode->tau) * snode->trait[idx]->brate;
  else
    v = 0.0;
  
  // this is a hack to avoid log likelihood becoming nan
  if (v < 1e-8) v = 1e-8;
  
  snode->trait[idx]->brlen = v;

  nchar = stree->trait_dim[idx];
  nstate = stree->trait_nstate[idx];
  max_state = nstate[nchar];

  /* calculate the transition probabilities */
  trait_trprob_mk(snode->trait[idx]->tranprob, v, max_state);

  /* pruning algorithm */
  /* TODO: this can be optimized by grouping characters by patterns */
  if (snode->left && snode->right)  /* internal node */
  {
    trait_update_cp(idx, snode->left, stree);
    trait_update_cp(idx, snode->right, stree);

    for (h = 0; h < nchar; ++h)
    {
      k = nstate[h];
      for (x = 0; x < k; ++x)
      {
        prob_l = prob_r = 0.0;
        for (y = 0; y < k; ++y)
        {
          if (y == x)
            tr_prob = snode->left->trait[idx]->tranprob[k-2][0];
          else
            tr_prob = snode->left->trait[idx]->tranprob[k-2][1];
          prob_l += tr_prob * snode->left->trait[idx]->condprob[h][y];
        }
        for (z = 0; z < k; ++z)
        {
          if (z == x)
            tr_prob = snode->right->trait[idx]->tranprob[k-2][0];
          else
            tr_prob = snode->right->trait[idx]->tranprob[k-2][1];
          prob_r += tr_prob * snode->right->trait[idx]->condprob[h][z];
        }
        snode->trait[idx]->condprob[h][x] = prob_l * prob_r;
      }
    }
    for (k = 2; k <= max_state; ++k) {
      for (a = 0; a < k; ++a) // dummy constant chars
      {
        j = nchar + k*(k-1)/2 - 1 + a;
        for (x = 0; x < k; ++x)
        {
          prob_l = prob_r = 0.0;
          for (y = 0; y < k; ++y)
          {
            if (y == x)
              tr_prob = snode->left->trait[idx]->tranprob[k-2][0];
            else
              tr_prob = snode->left->trait[idx]->tranprob[k-2][1];
            prob_l += tr_prob * snode->left->trait[idx]->condprob[j][y];
          }
          for (z = 0; z < k; ++z)
          {
            if (z == x)
              tr_prob = snode->right->trait[idx]->tranprob[k-2][0];
            else
              tr_prob = snode->right->trait[idx]->tranprob[k-2][1];
            prob_r += tr_prob * snode->right->trait[idx]->condprob[j][z];
          }
          snode->trait[idx]->condprob[j][x] = prob_l * prob_r;
        }
      }
    }
  }
  else  /* tip node */
  {
    for (h = 0; h < nchar; ++h)
    {
      /* set the conditional probabilities (L) to 1 for any state that is
         compatible with the observed state and to 0 otherwise */
      for (x = 0; x < nstate[h]; ++x)
      {
        if (snode->trait[idx]->state_d[h] & state_bin(x))
          snode->trait[idx]->condprob[h][x] = 1.;
        else
          snode->trait[idx]->condprob[h][x] = 0.;
      }
    }
    for (k = 2; k <= max_state; ++k) {
      for (a = 0; a < k; ++a) // dummy constant chars
      {
        j = nchar + k*(k-1)/2 - 1 + a;
        for (x = 0; x < k; ++x)
        {
          if (x == a)
            snode->trait[idx]->condprob[j][x] = 1.;
          else
            snode->trait[idx]->condprob[j][x] = 0.;
        }
      }
    }
  }

#ifdef DEBUG_Morph_Mkv
  printf("%s\n", snode->label);
  
  double ** condprob = snode->trait[idx]->condprob,
         ** tranprob = snode->trait[idx]->tranprob;
  for (k = 2; k <= max_state; ++k)
    printf("k=%d\tp0=%lf\tp1=%lf\n", k, tranprob[k-2][0], tranprob[k-2][1]);
  for (h = 0; h < nchar; ++h) {
    printf("ch%d\t", h+1);
    for (x = 0; x < nstate[h]; ++x)
      printf("%.2e ", condprob[h][x]);
    printf("\n");
  }
  for (k = 2; k <= max_state; ++k) {
    for (a = 0; a < k; ++a) // dummy constant chars
    {
      printf("dm%d\t", a);
      j = nchar + k*(k-1)/2 - 1 + a;
      for (x = 0; x < k; ++x)
        printf("%.2e ", condprob[j][x]);
      printf("\n");
    }
  }
  printf("\n");
#endif
}

static void trait_update_part(int idx, stree_t * stree)
{
#ifdef DEBUG_Chi
  // stree->nodes[3]->tau = 0.13;
  // stree->nodes[4]->tau = 0.08;
#endif

  if (stree->trait_type[idx] == BPP_DATA_DISC)
    trait_update_cp(idx, stree->root, stree);
  else if (stree->trait_type[idx] == BPP_DATA_CONT &&
           stree->trait_missing[idx])
    trait_update_lmr(idx, stree->root, stree);
  else
    trait_update_ic(idx, stree->root, stree);
}

void trait_update(stree_t * stree)
{
  int n;

  /* TODO: this can be optimized. Instead of updating all nodes,
     update the modified node and its ancestors all the way to the root. */
  for (n = 0; n < stree->trait_count; ++n)
    trait_update_part(n, stree);
}

static int trait_fill_tip(stree_t * stree, morph_t ** morph_list)
{
  int n, i, j, k, l, nchar, state, max_state;
  double value;
  snode_t * snode;
  morph_t * morph;
  
  /* copy the trait values over */
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode = stree->nodes[i];
    for (n = 0; n < stree->trait_count; ++n)
    {
      morph = morph_list[n];
      for (l = 0; l < morph->ntaxa; ++l)
      {
        if (strncmp(snode->label, morph->label[l], LABEL_LEN) == 0)
        {
          for (j = 0; j < morph->length; ++j) {
            if (morph->dtype == BPP_DATA_CONT)
              snode->trait[n]->state_m[j] = morph->conti[l][j];
            else
              snode->trait[n]->state_d[j] = morph->discr[l][j];
          }
          break;
        }
      }
      if (l == morph->ntaxa)
      {
        fprintf(stderr, "Species name %s not found in partition %d\n",
                snode->label, n+1);
        return 0;
      }
    }
  }
  
  for (n = 0; n < stree->trait_count; ++n)
  {
    nchar = stree->trait_dim[n];

    if (morph_list[n]->dtype == BPP_DATA_DISC)
    {
      /* check whether all discrete characters are variable */
      for (j = 0; j < nchar; ++j)
      {
        k = l = max_state = 0;
        for (i = 0; i < stree->tip_count; ++i)
        {
          state = stree->nodes[i]->trait[n]->state_d[j];
          if (state >= 1023) { // ? or -
            stree->trait_missing[n] = 1;
            l++;
            if (k == i) k++;
          }
          else if (state == stree->nodes[k]->trait[n]->state_d[j])
            l++;
          if (state < 1023 && state > max_state)
            max_state = state;
        }
        if (l == stree->tip_count)
        {
          fprintf(stderr, "Constant char at column %d partition %d\n", j, n);
          return 0;
        }
        
        /* record the number of states for each character */
        for (k = 2; state_bin(k) <= max_state; ++k);
        stree->trait_nstate[n][j] = k;
        
        /* record the max number of states of this partition */
        if (stree->trait_nstate[n][nchar] < k)
          stree->trait_nstate[n][nchar] = k;
      }
    }
    else // (morph_list[n]->dtype == BPP_DATA_CONT)
    {
      /* update the vector for the active coordinates */
      for (i = 0; i < stree->tip_count; ++i)
      {
        for (j = 0; j < nchar; ++j)
        {
          value = stree->nodes[i]->trait[n]->state_m[j];
          if (isnan(value))  // inapplicable value
          {
            stree->trait_missing[n] = 1;
            stree->nodes[i]->trait[n]->active[j] = -1;
          }
          else if (isinf(value))  // missing value
          {
            stree->trait_missing[n] = 1;
            stree->nodes[i]->trait[n]->active[j] = 0;
          }
          else {
            stree->nodes[i]->trait[n]->active[j] = 1;
          }
        }
      }
    }
  }
  
#ifdef DEBUG_Morph_Matrix
  for (n = 0; n < stree->trait_count; ++n)
  {
    printf("\n");
    for (i = 0; i < stree->tip_count; ++i)
    {
      snode = stree->nodes[i];
      printf("%s\t", snode->label);
      if (stree->trait_type[n] == BPP_DATA_DISC)
        for (j = 0; j < stree->trait_dim[n]; ++j)
          printf("%4d ", snode->trait[n]->state_d[j]);
      else
        for (j = 0; j < stree->trait_dim[n]; ++j)
          printf("%+.3lf\t", snode->trait[n]->state_m[j]);
      printf("\n");
    }
    if (stree->trait_type[n] == BPP_DATA_DISC)
    {
      printf("#\t");
      for (j = 0; j < stree->trait_dim[n] +1; ++j)
        printf("%4d ", stree->trait_nstate[n][j]);
      printf(" states\n");
    }
    else
    {
      for (i = 0; i < stree->tip_count; ++i)
      {
        snode = stree->nodes[i];
        printf("%s\t", snode->label);
        for (j = 0; j < stree->trait_dim[n]; ++j)
          printf("%4d ", snode->trait[n]->active[j]);
        printf("\n");
      }
      printf("\n");
    }
  }
#endif

  return 1;
}

static void trait_alloc_mem(stree_t * stree, morph_t ** morph_list, int n_part)
{
  int n, i, j, nchar;
  snode_t * snode;
  trait_t * trait;
  
  stree->trait_dim = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_type = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_missing = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_v_pop = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_ldetRs = (double *)xcalloc(n_part, sizeof(double));
  
  stree->trait_nstate = (int **)xcalloc(n_part, sizeof(int *));
  for (n = 0; n < n_part; ++n)
  {
    nchar = morph_list[n]->length;
    if (morph_list[n]->dtype == BPP_DATA_DISC)
    { /* for the number of states of each character
         use the last element to store the max number of states */
      stree->trait_nstate[n] = (int *)xcalloc(nchar +1, sizeof(int));
    }
  }
  
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->trait = (trait_t **)xmalloc(n_part*sizeof(trait_t *));
    for (n = 0; n < n_part; ++n)
    {
      snode->trait[n] = (trait_t *)xmalloc(sizeof(trait_t));
      trait = snode->trait[n];
      
      nchar = morph_list[n]->length;
      if (morph_list[n]->dtype == BPP_DATA_CONT)
      {
        trait->state_m  = (double *)xcalloc(nchar, sizeof(double));
        trait->contrast = (double *)xcalloc(nchar, sizeof(double));
        trait->active   = (int *)xcalloc(nchar, sizeof(int));
      }
      else
      {
        trait->state_d  = (int *)xcalloc(nchar, sizeof(int));
        /* each character has maximally ten states; the last 2+3+...+10=54
           cells are for storing conditional probs of dummy constant chars
           of 2, 3, ..., 10 states */
        trait->condprob = (double **)xcalloc(nchar +54, sizeof(double *));
        for (j = 0; j < nchar +54; ++j)
          trait->condprob[j] = (double *)xcalloc(10, sizeof(double));
        /* store the transition probabilities for characters
           of 2, 3, ..., 10 states */
        trait->tranprob = (double **)xcalloc(9, sizeof(double *));
        for (j = 0; j < 9; ++j)
          trait->tranprob[j] = (double *)xcalloc(j+2, sizeof(double));
      }
    }
  }
  
  stree->trait_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_logpr = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logpr = (double *)xcalloc(n_part, sizeof(double));
}

void trait_init(stree_t * stree, morph_t ** morph_list, int n_part)
{
  int n, i;

  assert(stree != NULL && morph_list != NULL && n_part > 0);
  for (n = 0; n < n_part; ++n) assert(morph_list[n] != NULL);
  
  /* allocate necessary memory */
  trait_alloc_mem(stree, morph_list, n_part);

  /* set up values */
  stree->trait_count = n_part;
  for (n = 0; n < n_part; ++n) {
    stree->trait_dim[n] = morph_list[n]->length;
    stree->trait_type[n] = morph_list[n]->dtype;
    if (morph_list[n]->dtype == BPP_DATA_CONT)
    {
      stree->trait_v_pop[n] = morph_list[n]->v_pop;
      stree->trait_ldetRs[n] = morph_list[n]->ldetRs;
    }
  }
  
  /* fill the trait values for the tip nodes */
  if (!trait_fill_tip(stree, morph_list))
  {
    trait_destroy(stree);
    fatal("Error filling traits");
  }
  
  /* initialize branch rates */
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    for (n = 0; n < n_part; ++n)
      stree->nodes[i]->trait[n]->brate = 1.0;
  
  /* then fill up relevant things */
  trait_update(stree);
}

static void trait_store_part(int idx, stree_t * stree)
{
  int i;
  snode_t * snode;

  stree->trait_old_logl[idx] = stree->trait_logl[idx];
  stree->trait_old_logpr[idx] = stree->trait_logpr[idx];

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->trait[idx]->old_brate = snode->trait[idx]->brate;
  }
}

static void trait_restore_part(int idx, stree_t * stree)
{
  int i;
  snode_t * snode;

  stree->trait_logl[idx] = stree->trait_old_logl[idx];
  stree->trait_logpr[idx] = stree->trait_old_logpr[idx];

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->trait[idx]->brate = snode->trait[idx]->old_brate;
  }
}

void trait_store(stree_t * stree)
{
  int n;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
    trait_store_part(n, stree);
}

void trait_restore(stree_t * stree)
{
  int n;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
    trait_restore_part(n, stree);
}

static double loglikelihood_trait_AC(int idx, stree_t * stree)
{
  int i, j, p;
  double v_k1, v_k2, zz, ldetRs, logl;
  snode_t * snode;

  /* log determinant of shrinkage estimate of the correlation matrix,
      i.e. log(det(R*)) */
  ldetRs = 0.0;  //TODO: stree->trait_ldetRs[idx];
  p = stree->trait_dim[idx];

  logl = 0.0;

  /* loop over the internal species nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    
    v_k1 = snode->left->trait[idx]->brlen;
    v_k2 = snode->right->trait[idx]->brlen;

    for (zz = 0, j = 0; j < p; ++j)
      zz += snode->trait[idx]->contrast[j] * snode->trait[idx]->contrast[j];
    
    /* Alvarez-Carretero et al. 2019. eq.5. */
    logl += -0.5* (p*log(2.0*BPP_PI*(v_k1+v_k2)) + ldetRs + zz/(v_k1+v_k2));
  }
  
  stree->trait_logl[idx] = logl;

#ifdef DEBUG_Morph_BM
  printf("part%d: cur log(like)=%lf, old log(like)=%lf\n\n", idx+1,
         stree->trait_logl[idx], stree->trait_old_logl[idx]);
#endif

  return logl;
}

/* TODO: */
static double loglikelihood_trait_Mitov(int idx, stree_t * stree)
{
  double logl;


  logl = 0.0;

  
  stree->trait_logl[idx] = logl;

  return logl;
}

static double loglikelihood_trait_Mkv(int idx, stree_t * stree)
{
  int h, j, k, a, x;
  int * nstate, nchar, max_state;
  double prob, p_c, p[9]={0}, logl;
  double ** root_prob;
  
  nchar = stree->trait_dim[idx];
  nstate = stree->trait_nstate[idx];
  max_state = nstate[nchar];
  root_prob = stree->root->trait[idx]->condprob;
  
  /* prob of being constant */
  for (k = 2; k <= max_state; ++k) {
    for (a = 0; a < k; ++a)
    {
      j = nchar + k*(k-1)/2 - 1 + a;
      for (x = 0; x < k; ++x)
        p[k-2] += root_prob[j][x] / k;
    }
  }
  
  logl = 0.0;

  /* assuming the characters are independent */
  for (h = 0; h < nchar; ++h)
  {
    /* average the conditional probabilities over the root states */
    prob = 0.0;
    for (x = 0; x < nstate[h]; ++x)
      prob += root_prob[h][x] / nstate[h];
    
    /* correct for the variable coding bias */
    p_c = p[nstate[h]-2];
    /* note: MrBayes uses p[0] for all characters, including those
       with more than two states (is that a good approximation?) */

    logl += log(prob) - log(1 - p_c);
  }

  stree->trait_logl[idx] = logl;

#ifdef DEBUG_Morph_Mkv
  printf("part%d: cur log(like)=%lf, old log(like)=%lf\n\n", idx+1,
         stree->trait_logl[idx], stree->trait_old_logl[idx]);
#endif

  return logl;
}

static double loglikelihood_trait_part(int idx, stree_t * stree)
{
  if (stree->trait_type[idx] == BPP_DATA_DISC)
  {
    /* Mkv model */
    return loglikelihood_trait_Mkv(idx, stree);
  }
  else if (stree->trait_type[idx] == BPP_DATA_CONT &&
           stree->trait_missing[idx])
  {
    /* BM model with missing data; Mitov et al. 2020 */
    return loglikelihood_trait_Mitov(idx, stree);
  }
  else
  {
    /* BM model without missing data; Alvarez-Carretero et al. 2019 */
    return loglikelihood_trait_AC(idx, stree);
  }
}

double loglikelihood_trait(stree_t * stree)
{
  int n;
  double logl_sum = 0.0;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
    logl_sum += loglikelihood_trait_part(n, stree);
  
  return logl_sum;
}

static double logprior_trait_part(int idx, stree_t * stree)
{
  int i;
  double logpr, a, b, x;
  snode_t * snode;

  logpr = 0.0;
  
  /* the branch rates follow i.i.d. gamma distributions
     with parameters opt_brate_m_alpha and opt_brate_m_beta */
  a = opt_brate_m_alpha;
  b = opt_brate_m_beta;

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    
    /* skip the root */
    if (!snode->parent)
      continue;
    
    x = snode->trait[idx]->brate;
 
    logpr += a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;
  }
  
  stree->trait_logpr[idx] = logpr;
  
  return logpr;
}

double logprior_trait(stree_t * stree)
{
  int n;
  double logpr_sum = 0.0;

  for (n = 0; n < stree->trait_count; ++n)
    logpr_sum += logprior_trait_part(n, stree);
  
  return logpr_sum;
}

static double prop_branch_rates_relax(stree_t * stree)
{
  int n, i,  proposed, accepted;
  long thread_index = 0;
  double old_rate, old_lograte, new_rate, new_lograte;
  double lnacceptance, a, b;
  snode_t * snode;

  /* the branch rates follow i.i.d. gamma distributions
     with parameters opt_brate_m_alpha and opt_brate_m_beta */
  a = opt_brate_m_alpha;
  b = opt_brate_m_beta;

  /* morphological rates are independent among partitions,
     and within a partition, each branch has a rate parameter */
  proposed = accepted = 0;
  for (n = 0; n < stree->trait_count; ++n)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode = stree->nodes[i];

      /* skip the root */
      if (!snode->parent)
        continue;
      
      old_rate = snode->trait[n]->brate;
      old_lograte = log(old_rate);
      new_lograte = old_lograte +
           opt_finetune_brate_m * legacy_rnd_symmetrical(thread_index);
      new_lograte = reflect(new_lograte, -99, 99, thread_index);
      snode->trait[n]->brate = new_rate = exp(new_lograte);
      
      lnacceptance = new_lograte - old_lograte;
      
      /* calculate the log prior difference */
      lnacceptance += (a - 1) * log(new_rate / old_rate)
                             - b * (new_rate - old_rate);
      
      /* then calculate the log likelihood difference */
      trait_update_part(n, stree);
      lnacceptance += loglikelihood_trait_part(n, stree)
                       - stree->trait_old_logl[n];
      
      if (lnacceptance >= -1e-10 ||
          legacy_rndu(thread_index) < exp(lnacceptance))
      {
        accepted++;
        trait_store_part(n, stree);
      }
      else  // rejected
      {
        trait_restore_part(n, stree);
      }
      
      proposed++;
    }
  }
  
  return (double)accepted/proposed;
}

static double prop_branch_rates_strict(stree_t * stree)
{
  int n, i,  proposed, accepted;
  long thread_index = 0;
  double old_rate, old_lograte, new_rate, new_lograte;
  double lnacceptance, a, b;
  snode_t * snode;
 
  /* the branch rates follow i.i.d. gamma distributions
     with parameters opt_brate_m_alpha and opt_brate_m_beta */
  a = opt_brate_m_alpha;
  b = opt_brate_m_beta;

  /* morphological rates are independent among partitions (?); within each
     partition, there is a single rate parameter shared across branches */
  proposed = accepted = 0;
  for (n = 0; n < stree->trait_count; ++n)
  {
    old_rate = stree->nodes[0]->trait[n]->brate;
    old_lograte = log(old_rate);
    new_lograte = old_lograte +
         opt_finetune_brate_m * legacy_rnd_symmetrical(thread_index);
    new_lograte = reflect(new_lograte, -99, 99, thread_index);
    new_rate = exp(new_lograte);
    
    lnacceptance = new_lograte - old_lograte;
    
    /* calculate the log prior difference */
    lnacceptance += (a - 1) * log(new_rate / old_rate)
                           - b * (new_rate - old_rate);

    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode = stree->nodes[i];
      
      /* skip the root */
      if (!snode->parent)
        continue;
      
      snode->trait[n]->brate = new_rate;
    }
      
    /* then calculate the log likelihood difference */
    trait_update_part(n, stree);
    lnacceptance += loglikelihood_trait_part(n, stree)
                     - stree->trait_old_logl[n];
    
    if (lnacceptance >= -1e-10 ||
        legacy_rndu(thread_index) < exp(lnacceptance))
    {
      accepted++;
      trait_store_part(n, stree);
    }
    else  // rejected
    {
      trait_restore_part(n, stree);
    }
    
    proposed++;
  }
  
  return (double)accepted/proposed;
}

double prop_branch_rates_trait(stree_t * stree)
{
  if (opt_clock == BPP_CLOCK_GLOBAL)
    return prop_branch_rates_strict(stree);
  else
    return prop_branch_rates_relax(stree);
}
