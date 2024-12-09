//
//  Created by Chi Zhang on 2023/12/5.
//

#include "bpp.h"

#define LABEL_LEN 99

#define DEBUG_Chi 1

/* get a non-blank character from file */
static int get_nb_char(FILE *fp)
{
  int c;
  while (isspace(c=fgetc(fp)));
  return c;
}

static int parse_header(FILE * fp, int * nrow, int * ncol, int * type,
                        double * v_pop, double * ldetRs)
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
    
    /* also read in here v_pop and ldetRs, which have been obtained when
       preprocessing the trait data */
    if (fscanf(fp, "%lf %lf", v_pop, ldetRs) != 2) {
      fprintf(stderr, "Expecting population variance and log(det(R*))\n");
      return 0;
    }
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
    else if (c == '-') // negative value
    {
      if(fscanf(fp, "%lf", &trait[j]) != 1)
        return 0;
      trait[j] = -trait[j];
    }
    else if (c == '?') // missing value
    {
      trait[j] = sqrt(-1); // nan
    }
    else // positive value
    {
      ungetc(c, fp);
      if(fscanf(fp, "%lf", &trait[j]) != 1)
        return 0;
    }
  }
  assert(j == len);

  return 1;
}

static int state_bin(int x)
{
  int std_bin[] = {1,2,4,8,16,32,64,128,256,512,1023,1024};
  
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
    
    if (isdigit(c))
    {
      std[j] = state_bin(c-'0');
    }
    else if (c == '?')
    {
      std[j] = 1023;
    }
    else if (c == '-')
    {
      std[j] = 1024;
    }
    else // TODO: ambiguity
    {
      fprintf(stderr, "Unrecognized trait value (%c)\n", c);
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
  if (!parse_header(fp, &(morph->ntaxa), &(morph->length), &(morph->dtype),
                        &(morph->v_pop), &(morph->ldetRs)))
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
     morph->trait has dimension m->count * m->length */
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
  
#ifdef DEBUG_Chi
  printf("\n%d %d %c\t", morph->ntaxa, morph->length, morph->dtype);
  if (morph->dtype == BPP_DATA_CONT)
    printf("%lf %lf\n", morph->v_pop, morph->ldetRs);
  else
    printf("\n");
  for (i = 0; i < morph->ntaxa; ++i) {
    printf("%s\t", morph->label[i]);
    if (morph->dtype == BPP_DATA_CONT)
      for (j = 0; j < morph->length; ++j)
        printf("%lf\t", morph->conti[i][j]);
    else
      for (j = 0; j < morph->length; ++j)
        printf("%d\t", morph->discr[i][j]);
    printf("\n");
  }
#endif
  
  return morph;
}

morph_t ** parse_traitfile(const char * traitfile, int * count)
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

  if (stree->trait_nstate)
  {
    for (n = 0; n < stree->trait_count; ++i)
      free(stree->trait_nstate[n]);
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

static int trait_fill_tip(stree_t * stree, morph_t ** morph_list)
{
  int n, i, j, k, l, nchar, state, max_state;
  snode_t * snode;
  morph_t * morph;
  
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode = stree->nodes[i];
    
    /* loop over the trait partitions */
    for (n = 0; n < stree->trait_count; ++n)
    {
      morph = morph_list[n];
      for (l = 0; l < morph->ntaxa; ++l)
      {
        if (strncmp(snode->label, morph->label[l], LABEL_LEN) == 0)
        {
          /* copy the trait values over */
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
    if (morph_list[n]->dtype == BPP_DATA_DISC)
    {
      nchar = stree->trait_dim[n];
      for (j = 0; j < nchar; ++j)
      {
        /* check whether all discrete characters are variable */
        l = max_state = 0;
        for (i = 0; i < stree->tip_count; ++i)
        {
          state = stree->nodes[i]->trait[n]->state_d[j];
          if (i == 0 || state >= 1023) // ? or -
            l++;
          else if (state == stree->nodes[0]->trait[n]->state_d[j])
            l++;
          if (state < 1023 && state > max_state)
            max_state = state;
        }
        if (l == stree->tip_count)
        {
          /* this happens when all states are the same (constant) */
          fprintf(stderr, "Constant char at column %d partition %d\n", j, n);
          return 0;
        }
        
        /* record the number of states for each character */
        for (k = 2; state_bin(k) < max_state; ++k);
        stree->trait_nstate[n][j] = k;
        
        /* record the max number of states of this partition */
        if (stree->trait_nstate[n][nchar] < k)
          stree->trait_nstate[n][nchar] = k;
      }
    }
  }
  
  return 1;
}

static void trait_update_pic_part(int idx, snode_t * snode, stree_t * stree)
{
  int j;
  double v_pop, v_k, v_k1, v_k2, *m_k1, *m_k2;

  if (snode->left && snode->right)  /* internal node */
  {
    trait_update_pic_part(idx, snode->left, stree);
    trait_update_pic_part(idx, snode->right, stree);
    
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
       same variance and the population noise has variance of v_pop.
       Alvarez-Carretero et al. 2019. p.970. */
    v_pop = stree->trait_v_pop[idx];
    snode->trait[idx]->brlen = v_k + v_pop;
  }
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

static void trait_update_cpl_part(int idx, snode_t * snode, stree_t * stree)
{
  int h, j, k, a, x, y, z;
  int * nstate, nchar, max_state;
  double v, prob_l, prob_r, tr_prob;
  
  /* update the branch length */
  if (snode->parent)
    v = (snode->parent->tau - snode->tau) * snode->trait[idx]->brate;
  else
    v = 0.0;
  snode->trait[idx]->brlen = v;

  nchar = stree->trait_dim[idx];
  nstate = stree->trait_nstate[idx];
  max_state = nstate[nchar];

  /* calculate the transition probabilities */
  trait_trprob_mk(snode->trait[idx]->tranprob, v, max_state);

  /* pruning algorithm */
  if (snode->left && snode->right)  /* internal node */
  {
    trait_update_cpl_part(idx, snode->left, stree);
    trait_update_cpl_part(idx, snode->right, stree);

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
}

void trait_update(stree_t * stree)
{
  int n;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
  {
    if (stree->trait_type[n] == BPP_DATA_CONT)
      trait_update_pic_part(n, stree->root, stree);
    else
      trait_update_cpl_part(n, stree->root, stree);
  }
}

void trait_init(stree_t * stree, morph_t ** morph_list, int n_part)
{
  int n, i, j;
  snode_t * snode;
  trait_t * trait;

  assert(stree != NULL && morph_list != NULL && n_part > 0);
  for (n = 0; n < n_part; ++n) assert(morph_list[n] != NULL);
  
  stree->trait_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_logpr = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logpr = (double *)xcalloc(n_part, sizeof(double));
  
  stree->trait_count = n_part;
  stree->trait_dim = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_type = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_v_pop = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_ldetRs = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_nstate = (int **)xcalloc(n_part, sizeof(int *));
  for (n = 0; n < n_part; ++n) {
    stree->trait_dim[n] = morph_list[n]->length;
    stree->trait_type[n] = morph_list[n]->dtype;
    if (morph_list[n]->dtype == BPP_DATA_CONT)
    {
      stree->trait_v_pop[n] = morph_list[n]->v_pop;
      stree->trait_ldetRs[n] = morph_list[n]->ldetRs;
    }
    else
    { /* use the last element to store the max number of states */
      stree->trait_nstate[n] =
          (int *)xcalloc(morph_list[n]->length +1, sizeof(int));
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
      if (morph_list[n]->dtype == BPP_DATA_CONT)
      {
        trait->state_m =
          (double *)xcalloc(morph_list[n]->length, sizeof(double));
        trait->contrast =
          (double *)xcalloc(morph_list[n]->length, sizeof(double));
      }
      else
      {
        trait->state_d =
             (int *)xcalloc(morph_list[n]->length, sizeof(int));
        /* each character has maximally ten states; the last 2+3+...+10=54
           cells are for storing conditional probs of dummy constant chars
           of 2, 3, ..., 10 states */
        trait->condprob =
         (double **)xcalloc(morph_list[n]->length +54, sizeof(double *));
        for (j = 0; j < morph_list[n]->length +54; ++j)
          trait->condprob[j] = (double *)xcalloc(10, sizeof(double));
        /* store the transition probabilities for characters
           of 2, 3, ..., 10 states */
        trait->tranprob = (double **)xcalloc(9, sizeof(double *));
        for (j = 0; j < 9; ++j)
          trait->tranprob[j] = (double *)xcalloc(j+2, sizeof(double));
      }
      
      /* initialize branch rates */
      trait->brate = 1.0;
    }
  }
  
  /* fill the trait values for the tip nodes */
  if (!trait_fill_tip(stree, morph_list))
  {
    trait_destroy(stree);
    fatal("Error filling traits");
  }
  
  /* then set up relevant things */
  for (n = 0; n < n_part; ++n)
  {
    if (morph_list[n]->dtype == BPP_DATA_CONT)
    {
      /* update the branch lengths and independent contrasts */
      trait_update_pic_part(n, stree->root, stree);
    }
    else
    {
      /* update the branch lengths and conditional probabilities */
      trait_update_cpl_part(n, stree->root, stree);
    }
  }
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

void debug_print_trait(int idx, stree_t * stree)
{
  int i, j;
  snode_t * snode;
  
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    printf("%s\t", snode->label);
//    printf(" r_k=%lf", snode->pic[idx]->brate);
//    printf(" v_k=%lf", snode->pic[idx]->brlen);
//    printf(" m_k:");
//    for (j = 0; j < stree->trait_dim[idx]; ++j)
//      printf(" %+lf", snode->pic[idx]->trait[j]);
//    printf(" x_k:");
//    for (j = 0; j < stree->trait_dim[idx]; ++j)
//      printf(" %+lf", snode->pic[idx]->contrast[j]);
    printf("\n");
  }
  
  printf("cur log(like)=%lf", stree->trait_logl[idx]);
  printf("\t log(prior)=%lf\n", stree->trait_logpr[idx]);
  printf("old log(like)=%lf", stree->trait_old_logl[idx]);
  printf("\t log(prior)=%lf\n", stree->trait_old_logpr[idx]);
}

static double loglikelihood_trait_c_bm(int idx, stree_t * stree)
{
  int i, j, p;
  double v_k1, v_k2, zz, ldetRs, logl;
  snode_t * snode;

  ldetRs = stree->trait_ldetRs[idx];
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

#ifdef DEBUG_Chi
  debug_print_trait(idx, stree);
#endif

  return logl;
}

static double loglikelihood_trait_d_mkv(int idx, stree_t * stree)
{
  int h, j, a, x;
  int * nstate, nchar;
  double prob, p0, logl;
  
  nchar = stree->trait_dim[idx];
  nstate = stree->trait_nstate[idx];
  
  logl = 0.0;

  /* assuming the characters are independent */
  for (h = 0; h < nchar; ++h)
  {
    /* average the conditional probabilities over the root states */
    prob = 0.0;
    for (x = 0; x < nstate[h]; ++x)
      prob += stree->root->trait[idx]->condprob[h][x] / nstate[h];
    
    /* correct for the variable coding bias */
    p0 = 0.0;  // prob of being constant
    for (a = 0; a < nstate[h]; ++a)
    {
      j = nchar + nstate[h]*(nstate[h]-1)/2 - 1 + a;
      for (x = 0; x < nstate[h]; ++x)
        p0 += stree->root->trait[idx]->condprob[j][x] / nstate[h];
    }

    logl += log(prob) - log(1-p0);
  }

  stree->trait_logl[idx] = logl;

#ifdef DEBUG_Chi
  debug_print_trait(idx, stree);
#endif

  return logl;
}

double loglikelihood_trait(stree_t * stree)
{
  int n;
  double logl_sum = 0.0;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
  {
    if (stree->trait_type[n] == BPP_DATA_CONT)
      logl_sum += loglikelihood_trait_c_bm(n, stree);
    else
      logl_sum += loglikelihood_trait_d_mkv(n, stree);
  }
  
  return logl_sum;
}

double logprior_trait(stree_t * stree)
{
  /* the branch rates follow i.i.d. gamma distributions
     with parameters opt_brate_alpha and opt_brate_beta */

  // TODO: calc prior

  return 0.0;
}

static double prop_branch_rates_relax(stree_t * stree)
{
  int n, i,  proposed, accepted;
  long thread_index = 0;
  double old_rate, old_lograte, new_rate, new_lograte;
  double lnacceptance;
  snode_t * snode;

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
      lnacceptance += (opt_brate_alpha-1) * log(new_rate/old_rate)
                         - opt_brate_beta * (new_rate - old_rate);
      
      if (stree->trait_type[n] == BPP_DATA_CONT)
      {
        /* update the contrasts as branch rate has been changed */
        trait_update_pic_part(n, stree->root, stree);
        
        /* then calculate the log likelihood difference */
        lnacceptance += loglikelihood_trait_c_bm(n, stree)
                      - stree->trait_old_logl[n];
      }
      else
      {
        /* update conditional probs as branch rate has been changed */
        trait_update_cpl_part(n, stree->root, stree);
        
        /* then calculate the log likelihood difference */
        lnacceptance += loglikelihood_trait_d_mkv(n, stree)
                      - stree->trait_old_logl[n];
      }
      
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
  double lnacceptance;
  snode_t * snode;
  
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
    lnacceptance += (opt_brate_alpha-1) * log(new_rate/old_rate)
                       - opt_brate_beta * (new_rate - old_rate);

    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode = stree->nodes[i];
      
      /* skip the root */
      if (!snode->parent)
        continue;
      
      snode->trait[n]->brate = new_rate;
    }
      
    if (stree->trait_type[n] == BPP_DATA_CONT)
    {
      /* update the contrasts as branch rate has been changed */
      trait_update_pic_part(n, stree->root, stree);
      
      /* then calculate the log likelihood difference */
      lnacceptance += loglikelihood_trait_c_bm(n, stree)
                    - stree->trait_old_logl[n];
    }
    else
    {
      /* update conditional probs as branch rate has been changed */
      trait_update_cpl_part(n, stree->root, stree);
      
      /* then calculate the log likelihood difference */
      lnacceptance += loglikelihood_trait_d_mkv(n, stree)
                    - stree->trait_old_logl[n];
    }
    
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
