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

static int parse_header(FILE * fp, int * nrow, int * ncol)
{
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

  /* TODO: read in here v_pop and ldetRs? */

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

static int parse_value(FILE * fp, double * trait, int len)
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

static trait_t * parse_trait_part(FILE * fp)
{
  int i, j;
  
  trait_t * morph = (trait_t *)xmalloc(sizeof(trait_t));
  
  /* read header */
  if (!parse_header(fp, &(morph->count), &(morph->length)))
  {
    fprintf(stderr, "Error in header\n");
    return NULL;
  }
  
  morph->label = (char **)xmalloc((morph->count)*sizeof(char *));
  morph->trait = (double **)xmalloc((morph->count)*sizeof(double *));
  for (i = 0; i < morph->count; ++i)
  {
    morph->trait[i] = (double *)xmalloc((morph->length)*sizeof(double));
    morph->label[i] = (char *)xmalloc((LABEL_LEN+1)*sizeof(char));
  }
  
  /* read morphological traits of each species,
     assuming each line has a label followed by m->length numbers.
     morph->trait has dimension m->count * m->length */
  for (i = 0; i < morph->count; ++i)
  {
    /* read the label (species name) */
    if (!parse_label(fp, morph->label[i], LABEL_LEN+1))
    {
      fprintf(stderr, "Failed to read label of species %d\n", i+1);
      return NULL;
    }
    /* read the trait values of this species */
    if (!parse_value(fp, morph->trait[i], morph->length))
    {
      fprintf(stderr, "Failed to read traits of species %d\n", i+1);
      return NULL;
    }
  }
  
#ifdef DEBUG_Chi
  printf("\n%d %d\n", morph->count, morph->length);
  for (i = 0; i < morph->count; ++i) {
    printf("%s\t", morph->label[i]);
    for (j = 0; j < morph->length; ++j)
      printf("%lf\t", morph->trait[i][j]);
    printf("\n");
  }
#endif
  
  morph->dtype = BPP_DATA_TRAIT;
  
  return morph;
}

trait_t ** parse_traitfile(const char * traitfile, int * count)
{
  int c, i, m_slotalloc = 10, m_maxcount = 0;
  
  FILE * fp = xopen(traitfile,"r");

  trait_t ** pmorph = (trait_t **)xmalloc(m_maxcount*sizeof(trait_t *));

  *count = 0;

  while ((c=fgetc(fp)) != EOF)
  {
    if (isspace(c)) continue;
    ungetc(c,fp);
    
    if (*count == m_maxcount)
    {
      m_maxcount += m_slotalloc;
      trait_t ** temp = (trait_t **)xmalloc(m_maxcount*sizeof(trait_t *));
      memcpy(temp, pmorph, (*count)*sizeof(trait_t *));
      free(pmorph);
      pmorph = temp;
    }

    pmorph[*count] = parse_trait_part(fp);
    if (pmorph[*count] == NULL)
    {
      for (i = 0; i < *count; ++i)
        trait_destroy(pmorph[i]);
      free(pmorph);

      fatal("Error parsing trait partition %d", *count+1);
    }

    *count += 1;  // count the number of partitions
  }

  fclose(fp);

  return pmorph;
}

void trait_destroy(trait_t * morph)
{
  int i;

  if (morph->label)
  {
    for (i = 0; i < morph->count; ++i)
      if (morph->label[i])
        free(morph->label[i]);
    free(morph->label);
  }

  if (morph->trait)
  {
    for (i = 0; i < morph->count; ++i)
      if (morph->trait[i])
        free(morph->trait[i]);
    free(morph->trait);
  }

  free(morph);
}

void pic_destroy(stree_t * stree)
{
  int n, i;
  snode_t * snode;
  contrast_t * ic;

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    if (snode->pic)
    {
      for (n = 0; n < stree->trait_count; ++n)
      {
        ic = snode->pic[n];
        if (ic)
        {
          if(ic->trait)
            free(ic->trait);
          if(ic->contrast)
            free(ic->contrast);
          free(ic);
        }
      }
      free(snode->pic);
    }
  }
  
  if (stree->trait_dim)
    free(stree->trait_dim);
  
  if (stree->trait_logl)
    free(stree->trait_logl);
  if (stree->trait_old_logl)
    free(stree->trait_old_logl);
  if (stree->trait_logpr)
    free(stree->trait_logpr);
  if (stree->trait_old_logpr)
    free(stree->trait_old_logpr);
}

void debug_print_pic(int idx, stree_t * stree)
{
  int i, j;
  snode_t * snode;
  
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    printf("%s\t", snode->label);
    printf(" r_k=%lf", snode->pic[idx]->brate);
    printf(" v_k=%lf", snode->pic[idx]->brlen);
    printf(" m_k:");
    for (j = 0; j < stree->trait_dim[idx]; ++j)
      printf(" %+lf", snode->pic[idx]->trait[j]);
    printf(" x_k:");
    for (j = 0; j < stree->trait_dim[idx]; ++j)
      printf(" %+lf", snode->pic[idx]->contrast[j]);
    printf("\n");
  }
  
  printf("cur log(like)=%lf", stree->trait_logl[idx]);
  printf("\t log(prior)=%lf\n", stree->trait_logpr[idx]);
  printf("old log(like)=%lf", stree->trait_old_logl[idx]);
  printf("\t log(prior)=%lf\n", stree->trait_old_logpr[idx]);
}

static int pic_fill_tip(stree_t * stree, trait_t ** trait_list)
{
  int n, i, j, l;
  snode_t * snode;
  trait_t * morph;
  
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode = stree->nodes[i];
    
    /* loop over the trait partitions */
    for (n = 0; n < stree->trait_count; ++n)
    {
      morph = trait_list[n];
      for (l = 0; l < morph->count; ++l)
      {
        if (strncmp(snode->label, morph->label[l], LABEL_LEN) == 0)
        {
          /* copy the trait values over */
          for (j = 0; j < morph->length; ++j)
            snode->pic[n]->trait[j] = morph->trait[l][j];
          break;
        }
      }
      if (l == morph->count)
      {
        fprintf(stderr, "Species name %s not found in partition %d\n",
                snode->label, n+1);
        return 0;
      }
    }
  }
  
  return 1;
}

static void pic_update_part(int idx, snode_t * snode, int dim)
{
  int j;
  double v_pop, v_k, v_k1, v_k2, *m_k1, *m_k2;

  if (snode->left && snode->right)  /* internal node */
  {
    pic_update_part(idx, snode->left, dim);
    pic_update_part(idx, snode->right, dim);
    
    if (snode->parent)
      v_k = (snode->parent->tau - snode->tau) * snode->pic[idx]->brate;
    else
      v_k = 0;
    v_k1 = snode->left->pic[idx]->brlen;
    v_k2 = snode->right->pic[idx]->brlen;
    snode->pic[idx]->brlen = v_k + v_k1*v_k2/(v_k1+v_k2);  // v_k'
    
    m_k1 = snode->left->pic[idx]->trait;
    m_k2 = snode->right->pic[idx]->trait;
    for (j = 0; j < dim; ++j)
    {
      snode->pic[idx]->contrast[j] = m_k1[j] - m_k2[j];  // x_k and m_k'
      snode->pic[idx]->trait[j] = (v_k2*m_k1[j] + v_k1*m_k2[j]) /(v_k1+v_k2);
    }
  }
  else  /* tip node */
  {
    v_k = (snode->parent->tau - snode->tau) * snode->pic[idx]->brate;
    /* The trait matrix has been standardized so that all characters have the
       same variance and the population noise has variance of v_pop.
       Alvarez-Carretero et al. 2019. p.970. */
    v_pop = 0.01;  /* TODO: //Chi */
    snode->pic[idx]->brlen = v_k + v_pop;
  }
}

void pic_update(stree_t * stree)
{
  int n;

  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
    pic_update_part(n, stree->root, stree->trait_dim[n]);
}

void pic_init(stree_t * stree, trait_t ** trait_list, int n_part)
{
  int n, i;
  snode_t * snode;
  
  assert(stree != NULL && trait_list != NULL && n_part > 0);
  for (n = 0; n < n_part; ++n) assert(trait_list[n] != NULL);
  
  stree->trait_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_logpr = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logpr = (double *)xcalloc(n_part, sizeof(double));

  stree->trait_count = n_part;
  stree->trait_dim = (int *)xcalloc(n_part, sizeof(int));
  for (n = 0; n < n_part; ++n)
    stree->trait_dim[n] = trait_list[n]->length;
  
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->pic = (contrast_t **)xmalloc(n_part*sizeof(contrast_t *));
    for (n = 0; n < n_part; ++n)
    {
      snode->pic[n] = (contrast_t *)xmalloc(sizeof(contrast_t));
      snode->pic[n]->trait = (double *)xcalloc(trait_list[n]->length,
                                               sizeof(double));
      snode->pic[n]->contrast = (double *)xcalloc(trait_list[n]->length,
                                                  sizeof(double));
      /* initialize branch rates */
      snode->pic[n]->brate = 1.0;
    }
  }
  
  /* fill the trait values for the tip nodes */
  if (!pic_fill_tip(stree, trait_list))
  {
    pic_destroy(stree);
    fatal("Error filling traits");
  }
  
  /* then update the independent contrasts */
  pic_update(stree);
}

static void pic_store_part(int idx, stree_t * stree)
{
  int i;
  snode_t * snode;

  stree->trait_old_logl[idx] = stree->trait_logl[idx];
  stree->trait_old_logpr[idx] = stree->trait_logpr[idx];

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->pic[idx]->old_brate = snode->pic[idx]->brate;
  }
}

static void pic_restore_part(int idx, stree_t * stree)
{
  int i;
  snode_t * snode;

  stree->trait_logl[idx] = stree->trait_old_logl[idx];
  stree->trait_logpr[idx] = stree->trait_old_logpr[idx];

  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->pic[idx]->brate = snode->pic[idx]->old_brate;
  }
}

void pic_store(stree_t * stree)
{
  int n;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
    pic_store_part(n, stree);
}

void pic_restore(stree_t * stree)
{
  int n;
  
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
    pic_restore_part(n, stree);
}

static double loglikelihood_trait_part(int idx, stree_t * stree)
{
  int i, j, p;
  double v_k1, v_k2, zz, ldetRs, logl;
  snode_t * snode;

  ldetRs = 0.0; /* TODO: //Chi */
  p = stree->trait_dim[idx];

  logl = 0.0;
  /* loop over the internal species nodes */
  for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    
    v_k1 = snode->left->pic[idx]->brlen;
    v_k2 = snode->right->pic[idx]->brlen;

    for (zz = 0, j = 0; j < p; ++j)
      zz += snode->pic[idx]->contrast[j] * snode->pic[idx]->contrast[j];
    
    /* Alvarez-Carretero et al. 2019. eq.5. */
    logl += -0.5* (p*log(2.0*BPP_PI*(v_k1+v_k2)) + ldetRs + zz/(v_k1+v_k2));
  }
  
  stree->trait_logl[idx] = logl;

#ifdef DEBUG_Chi
  debug_print_pic(idx, stree);
#endif

  return logl;
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


static double prop_branch_rates_relax(stree_t * stree)
{
  int n, i,  proposed, accepted;
  long thread_index = 0;
  double old_rate, old_lograte, new_rate, new_lograte;
  double lnacceptance;
  snode_t * snode;

  proposed = accepted = 0;
  for (n = 0; n < stree->trait_count; ++n)
  {
    /* morphological rates are independent among partitions,
       and within a partition, each branch has a rate parameter */
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode = stree->nodes[i];

      /* skip the root */
      if (!snode->parent)
        continue;
      
      old_rate = snode->pic[n]->brate;
      old_lograte = log(old_rate);
      new_lograte = old_lograte +
           opt_finetune_brate_m * legacy_rnd_symmetrical(thread_index);
      new_lograte = reflect(new_lograte, -99, 99, thread_index);
      snode->pic[n]->brate = new_rate = exp(new_lograte);
      
      lnacceptance = new_lograte - old_lograte;
      
      /* calculate the log prior difference */
      lnacceptance += (opt_brate_alpha-1) * log(new_rate/old_rate)
                         - opt_brate_beta * (new_rate - old_rate);
      
      /* update the contrasts as branch rate has been changed */
      pic_update_part(n, stree->root, stree->trait_dim[n]);
      
      /* then calculate the log likelihood difference */
      lnacceptance += loglikelihood_trait_part(n, stree)
                      - stree->trait_old_logl[n];
      
      if (lnacceptance >= -1e-10 ||
          legacy_rndu(thread_index) < exp(lnacceptance))
      {
        accepted++;
        pic_store_part(n, stree);
      }
      else  // rejected
      {
        pic_restore_part(n, stree);
      }
      
      proposed++;
    }
  }
  
  return (double)accepted/proposed;
}

static double prop_branch_rates_strict(stree_t * stree)
{
  int n, i, accepted;
  long thread_index = 0;
  double old_rate, old_lograte, new_rate, new_lograte;
  double lnacceptance;
  snode_t * snode;
  
  /* there is a single rate parameter shared across branches and among
     partitions, thus we propose one and update them all */
  old_rate = stree->nodes[0]->pic[0]->brate;
  old_lograte = log(old_rate);
  new_lograte = old_lograte +
       opt_finetune_brate_m * legacy_rnd_symmetrical(thread_index);
  new_lograte = reflect(new_lograte, -99, 99, thread_index);
  new_rate = exp(new_lograte);
  
  lnacceptance = new_lograte - old_lograte;

  /* calculate the log prior difference */
  lnacceptance += (opt_brate_alpha-1) * log(new_rate/old_rate)
                     - opt_brate_beta * (new_rate - old_rate);

  for (n = 0; n < stree->trait_count; ++n)
  {
    for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode = stree->nodes[i];
      
      /* skip the root */
      if (!snode->parent)
        continue;
      
      snode->pic[n]->brate = new_rate;
    }
  }
      
  /* update the contrasts as branch rate has been changed */
  pic_update(stree);
  
  /* then calculate the log likelihood difference */
  for (n = 0; n < stree->trait_count; ++n)
    lnacceptance += loglikelihood_trait_part(n, stree)
                  - stree->trait_old_logl[n];
  
  if (lnacceptance >= -1e-10 || legacy_rndu(thread_index) < exp(lnacceptance))
  {
    accepted = 1;
    pic_store(stree);
  }
  else  // rejected
  {
    accepted = 0;
    pic_restore(stree);
  }
  
  return (double)accepted;
}

double prop_branch_rates_trait(stree_t * stree)
{
  if (opt_clock == BPP_CLOCK_GLOBAL)
    return prop_branch_rates_strict(stree);
  else
    return prop_branch_rates_relax(stree);
}
