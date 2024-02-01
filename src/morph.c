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

  /* TODO: read in here ldetRs? */

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
  
#if(DEBUG_Chi)
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
      for (n = 0; n < 2*stree->trait_count; ++n)
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
}

static void debug_print_pic(stree_t * stree)
{
  int n, i, j, dim;
  snode_t * snode;
  
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    printf("%s\t", snode->label);
    for (n = 0; n < 2*stree->trait_count; ++n)
    {
      printf("# %lf |", snode->pic[n]->brlen);
      printf("  %lf |", snode->pic[n]->brate);
      if (n < stree->trait_count)
        dim = stree->trait_dim[n];
      else
        dim = stree->trait_dim[n- stree->trait_count];
      for (j = 0; j < dim; ++j)
        printf(" %+lf", snode->pic[n]->trait[j]);
      printf(" |");
      for (j = 0; j < dim; ++j)
        printf(" %+lf", snode->pic[n]->contrast[j]);
      if (n == stree->trait_count-1)
        printf("\n \t");
    }
    printf("\n");
  }
  printf("\n");
  
  printf("cur log(like):");
  for (n = 0; n < stree->trait_count; ++n)
    printf(" %lf", stree->trait_logl[n]);
  printf("\n");
  printf("old log(like):");
  for (n = 0; n < stree->trait_count; ++n)
    printf(" %lf", stree->trait_old_logl[n]);
  printf("\n");
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

void pic_update(snode_t * snode, int * ncol, int n_part)
{
  int n, j;
  double v_k, v_k1, v_k2, *m_k1, *m_k2;
  contrast_t * ic;

  if (snode->left != NULL && snode->right != NULL)  /* internal node */
  {
    pic_update(snode->left,  ncol, n_part);
    pic_update(snode->right, ncol, n_part);
    
    /* loop over the trait partitions */
    for (n = 0; n < n_part; ++n)
    {
      ic = snode->pic[n];
      ic->brate = 1.0;
      if (snode->parent != NULL)
        v_k = (snode->parent->tau - snode->tau) * ic->brate;
      else
        v_k = 0;
      v_k1 = snode->left->pic[n]->brlen;
      v_k2 = snode->right->pic[n]->brlen;
      ic->brlen = v_k + v_k1*v_k2/(v_k1+v_k2);                     // v_k'
      
      m_k1 = snode->left->pic[n]->trait;
      m_k2 = snode->right->pic[n]->trait;
      for (j = 0; j < ncol[n]; ++j)
      {
        ic->contrast[j] = m_k1[j] - m_k2[j];                       // x_k
        ic->trait[j] = (v_k2*m_k1[j] + v_k1*m_k2[j]) /(v_k1+v_k2); // m_k'
      }
    }
  }
  else  /* tip node */
  {
    for (n = 0; n < n_part; ++n)
    {
      ic = snode->pic[n];
      ic->brate = 1.0;
      v_k = (snode->parent->tau - snode->tau) * ic->brate;
      /* The trait matrix has been standardized so that all characters have the
         same variance and so that the population noise has unit variance. See
         Alvarez-Carretero et al. 2019. p.970. */
      ic->brlen = v_k + 1.0;
    }
  }
}

void pic_init(stree_t * stree, trait_t ** trait_list, int n_part)
{
  int n, i, dim;
  snode_t * snode;
  
  assert(stree != NULL && trait_list != NULL && n_part > 0);
  for (n = 0; n < n_part; ++n) assert(trait_list[n] != NULL);
  
  stree->trait_logl = (double *)xcalloc(n_part, sizeof(double));
  stree->trait_old_logl = (double *)xcalloc(n_part, sizeof(double));
  
  stree->trait_count = n_part;
  stree->trait_dim = (int *)xcalloc(n_part, sizeof(int));
  for (n = 0; n < n_part; ++n)
    stree->trait_dim[n] = trait_list[n]->length;
  
  /* allocate twice the amount of space for store/restore */
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
    snode->pic = (contrast_t **)xmalloc(2*n_part*sizeof(contrast_t *));
    for (n = 0; n < 2*n_part; ++n)
    {
      snode->pic[n] = (contrast_t *)xmalloc(sizeof(contrast_t));
      if (n < n_part)
        dim = stree->trait_dim[n];
      else
        dim = stree->trait_dim[n- n_part];
      snode->pic[n]->trait    = (double *)xcalloc(dim, sizeof(double));
      snode->pic[n]->contrast = (double *)xcalloc(dim, sizeof(double));
    }
  }
  
  /* fill the trait values for the tip nodes */
  if (!pic_fill_tip(stree, trait_list))
  {
    pic_destroy(stree);
    fatal("Error filling traits");
  }
  
  /* update the independent contrasts */
  pic_update(stree->root, stree->trait_dim, n_part);
}

void pic_store(stree_t * stree)
{
  int n, i, j, n_part;
  snode_t * snode;

  n_part = stree->trait_count;
  for (n = 0; n < n_part; ++n)
    stree->trait_old_logl[n] = stree->trait_logl[n];
    
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
      
    for (n = 0; n < n_part; ++n)
    {
      snode->pic[n+ n_part]->brate = snode->pic[n]->brate;
      snode->pic[n+ n_part]->brlen = snode->pic[n]->brlen;
      
      if (i >= stree->tip_count)
      for (j = 0; j < stree->trait_dim[n]; ++j)
      {
        snode->pic[n+ n_part]->trait[j] = snode->pic[n]->trait[j];
        snode->pic[n+ n_part]->contrast[j] = snode->pic[n]->contrast[j];
      }
    }
  }
}

void pic_restore(stree_t * stree)
{
  int n, i, j, n_part;
  snode_t * snode;

  n_part = stree->trait_count;
  for (n = 0; n < n_part; ++n)
    stree->trait_logl[n] = stree->trait_old_logl[n];
    
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    snode = stree->nodes[i];
      
    for (n = 0; n < n_part; ++n)
    {
      snode->pic[n]->brate = snode->pic[n+ n_part]->brate;
      snode->pic[n]->brlen = snode->pic[n+ n_part]->brlen;
      
      if (i >= stree->tip_count)
      for (j = 0; j < stree->trait_dim[n]; ++j)
      {
        snode->pic[n]->trait[j] = snode->pic[n+ n_part]->trait[j];
        snode->pic[n]->contrast[j] = snode->pic[n+ n_part]->contrast[j];
      }
    }
  }
}

double loglikelihood_trait(stree_t * stree)
{
  int n, i, j, p;
  double v_k1, v_k2, zz, ldetRs, logl, logl_sum;
  snode_t * snode;

  logl_sum = 0.0;
  /* loop over the trait partitions */
  for (n = 0; n < stree->trait_count; ++n)
  {
    logl = 0.0;
    
    ldetRs = 0.0; /* TODO: //Chi */
    p = stree->trait_dim[n];

    /* loop over the internal species nodes */
    for (i = stree->tip_count; i < stree->tip_count+stree->inner_count; ++i)
    {
      snode = stree->nodes[i];
      
      v_k1 = snode->left->pic[n]->brlen;
      v_k2 = snode->right->pic[n]->brlen;

      for (zz = 0, j = 0; j < p; ++j)
        zz += snode->pic[n]->contrast[j] * snode->pic[n]->contrast[j];
      
      /* Alvarez-Carretero et al. 2019. eq.5. */
      logl += -0.5* (p*log(2.0*BPP_PI*(v_k1+v_k2)) + ldetRs + zz/(v_k1+v_k2));
    }
    
    stree->trait_logl[n] = logl;
    logl_sum += logl;
  }
  
#if(DEBUG_Chi)
  debug_print_pic(stree);
#endif
  
  return logl_sum;
}
