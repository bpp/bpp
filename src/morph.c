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
  int n = fscanf(fp, "%d %d", nrow, ncol);
  
  if (n != 2)
    fatal("Error in header: expecting two numbers");
  else if (*nrow <= 0)
    fatal("Error in header: number of species must be >0 (%d)", *nrow);
  else if (*ncol <= 0)
    fatal("Error in header: number of traits must be >0 (%d)", *ncol);

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
  int c, j, n;
  
  for (j = 0; j < len; ++j)
  {
    c = get_nb_char(fp);
    
    if (c == EOF)
    {
      return 0;
    }
    else if (c == '-') // negative value
    {
      if((n = fscanf(fp, "%lf", &trait[j])) != 1)
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
      if((n = fscanf(fp, "%lf", &trait[j])) != 1)
        return 0;
    }
  }
  assert(j == len);

  return 1;
}

trait_t * parse_trait_part(FILE * fp)
{
  int i, j;
  
  trait_t * morph = (trait_t *)xmalloc(sizeof(trait_t));
  
  /* read header */
  if (!parse_header(fp, &(morph->count), &(morph->length)))
    return NULL;

  morph->label = (char **)xcalloc((morph->count),sizeof(char *));
  morph->trait = (double **)xcalloc((morph->count),sizeof(double *));
  for (i = 0; i < morph->count; ++i)
  {
    morph->trait[i] = (double *)xmalloc((morph->length)*sizeof(double));
    morph->label[i] = (char *)xmalloc((LABEL_LEN+1)*sizeof(char));
  }
  
  /* read morphological traits of each species,
     assuming each line has a label followed by m->length numbers.
     m->trait has dimension m->count * m->length */
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

trait_t ** parse_traitfile(const char * traitfile, long * count)
{
  int c;
  long m_slotalloc = 10;
  long m_maxcount = 0;
  
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
      fatal("Error parsing trait partition %d", (*count)+1);

    *count += 1;
  }

  fclose(fp);

  return pmorph;
}

static void pic_init_recursive(snode_t * snode, long n_part, int * dim)
{
  int i, j, n;
  double v_k, v_k1, v_k2, *m_k1, *m_k2;
  contrast_t * ic;

  /* assuming the tips already have trait values filled in */
  if (snode->left != NULL && snode->right != NULL)
  {
    pic_init_recursive(snode->left, n_part, dim);
    pic_init_recursive(snode->right, n_part, dim);
    
    /* allocate twice the amount of space for store/restore */
    snode->pic = (contrast_t **)xmalloc(2*n_part*sizeof(contrast_t *));

    /* loop over the trait partitions */
    for (n = 0; n < n_part; ++n)
    {
      ic = snode->pic[n] = (contrast_t *)xmalloc(sizeof(contrast_t));
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
      ic->trait = (double *)xmalloc(dim[n]*sizeof(double));
      ic->contrast = (double *)xmalloc(dim[n]*sizeof(double));
      for (j = 0; j < dim[n]; ++j)
      {
        ic->contrast[j] = m_k1[j] - m_k2[j];                       // x_k
        ic->trait[j] = (v_k2*m_k1[j] + v_k1*m_k2[j]) /(v_k1+v_k2); // m_k'
      }
    }
  }
}

void pic_init(stree_t * stree, trait_t ** trait_list, long n_part)
{
  int i, j, n, l, found;
  snode_t * snode;
  contrast_t * ic;
  trait_t * morph;
  
  assert(stree != NULL && trait_list != NULL && n_part > 0);
  for (n = 0; n < n_part; ++n) assert(trait_list[n] != NULL);
  
  /* tip node(_k) has trait values (m_k) and branch length (v_k) as is */
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode = stree->nodes[i];

    /* allocate twice the amount of space for store/restore */
    snode->pic = (contrast_t **)xmalloc(2*n_part*sizeof(contrast_t *));
    
    /* loop over the trait partitions */
    for (n = 0; n < n_part; ++n)
    {
      ic = snode->pic[n] = (contrast_t *)xmalloc(sizeof(contrast_t));
      ic->brate = 1.0;
      ic->brlen = (snode->parent->tau - snode->tau) * ic->brate;
      /* TODO: how to deal with zero branch length */

      /* find and fill the trait values for this tip node */
      morph = trait_list[n];
      found = 0;
      for (l = 0; l < morph->count; ++l)
      {
        if (strncmp(snode->label, morph->label[l], LABEL_LEN) == 0)
        {
          /* copy the trait values over */
          ic->trait = (double *)xmalloc((morph->length)*sizeof(double));
          for (j = 0; j < morph->length; ++j)
            ic->trait[j] = morph->trait[l][j];
          
          found++;
          break;
        }
      }
      if (found == 0)
      {
        fatal("Species name %s not found in partition %d", snode->label, n+1);
      }
    }
  }
  
  stree->trait_dim = (int *)xmalloc(n_part*sizeof(int));
  for (n = 0; n < n_part; ++n)
    stree->trait_dim[n] = trait_list[n]->length;

  /* internal node: compute contrasts (x_k), ancestral traits (m_k'), and
                    transformed branch length (v_k') */
  pic_init_recursive(stree->root, n_part, stree->trait_dim);
  
#if(DEBUG_Chi)
  for (n = 0; n < n_part; ++n)
  {
    printf("Trait partition %d\n", n+1);
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    {
      snode = stree->nodes[i];
      printf("%s\t", snode->label);
      printf("v'=%lf, ", snode->pic[n]->brlen);
      printf("m':");
      for (j = 0; j < stree->trait_dim[n]; ++j)
        printf(" %lf", snode->pic[n]->trait[j]);
      printf(", ");
      if (i >= stree->tip_count)
      {
        printf("x:");
        for (j = 0; j < stree->trait_dim[n]; ++j)
          printf(" %lf", snode->pic[n]->contrast[j]);
      }
      printf("\n");
    }
    printf("\n");
  }
#endif
}
