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
  
  trait_t * m = (trait_t *)xmalloc(sizeof(trait_t));
  
  /* read header */
  if (!parse_header(fp, &(m->count), &(m->length)))
    return NULL;

  m->label = (char **)xcalloc((m->count),sizeof(char *));
  m->trait = (double **)xcalloc((m->count),sizeof(double *));
  for (i = 0; i < m->count; ++i)
  {
    m->trait[i] = (double *)xmalloc((m->length)*sizeof(double));
    m->label[i] = (char *)xmalloc((LABEL_LEN+1)*sizeof(char));
  }
  
  /* read morphological traits of each species,
     assuming each line has a label followed by m->length numbers.
     m->trait has dimension m->count * m->length */
  for (i = 0; i < m->count; ++i)
  {
    /* read the label (species name) */
    if (!parse_label(fp, m->label[i], LABEL_LEN+1))
    {
      fprintf(stderr, "Failed to read label of species %d\n", i+1);
      return NULL;
    }
    /* read the trait values of this species */
    if (!parse_value(fp, m->trait[i], m->length))
    {
      fprintf(stderr, "Failed to read traits of species %d\n", i+1);
      return NULL;
    }
  }
  
#if(DEBUG_Chi)
  printf("\n%d %d\n", m->count, m->length);
  for (i = 0; i < m->count; ++i) {
    printf("%s\t", m->label[i]);
    for (j = 0; j < m->length; ++j)
      printf("%lf\t", m->trait[i][j]);
    printf("\n");
  }
#endif
  
  m->dtype = BPP_DATA_TRAIT;
  
  return m;
}

trait_t ** parse_traitfile(const char * traitfile, long * count)
{
  int c;
  long m_slotalloc = 10;
  long m_maxcount = 0;
  
  FILE * fp = xopen(traitfile,"r");

  trait_t ** mm = (trait_t **)xmalloc(m_maxcount*sizeof(trait_t *));

  *count = 0;

  while ((c=fgetc(fp)) != EOF)
  {
    if (isspace(c)) continue;
    ungetc(c,fp);
    
    if (*count == m_maxcount)
    {
      m_maxcount += m_slotalloc;
      trait_t ** temp = (trait_t **)xmalloc(m_maxcount*sizeof(trait_t *));
      memcpy(temp, mm, *count * sizeof(trait_t *));
      free(mm);
      mm = temp;
    }

    mm[*count] = parse_trait_part(fp);
    if (mm[*count] == NULL)
      fatal("Error parsing trait partition %d", (*count)+1);

    *count += 1;
  }

  fclose(fp);

  return mm;
}

static void pic_init_recursive(snode_t * snode, long n_part, int * dim)
{
  int i, j, k;
  double v_k, v_k1, v_k2, *m_k1, *m_k2;
  contrast_t * ic;

  /* assuming the tips already have trait values filled in */
  if (snode->left != NULL && snode->right != NULL)
  {
    pic_init_recursive(snode->left, n_part, dim);
    pic_init_recursive(snode->right, n_part, dim);
    
    /* allocate twice the amount of space for store/restore */
    snode->pic = (contrast_t **)xmalloc(2*n_part*sizeof(contrast_t *));
    for (k = 0; k < 2*n_part; ++k)
    {
      ic = snode->pic[k] = (contrast_t *)xmalloc(sizeof(contrast_t));
      ic->trait = (double *)xmalloc(dim[k]*sizeof(double));
      ic->contrast = (double *)xmalloc(dim[k]*sizeof(double));
    }
    
    /* loop over the trait partitions */
    for (k = 0; k < n_part; ++k)
    {
      ic = snode->pic[k];

      if (snode->parent != NULL)
        v_k = snode->parent->tau - snode->tau;
      else
        v_k = 0;
      v_k1 = snode->left->pic[k]->brlen;
      v_k2 = snode->right->pic[k]->brlen;
      ic->brlen = v_k + v_k1*v_k2/(v_k1+v_k2);                     // v_k'
      
      m_k1 = snode->left->pic[k]->trait;
      m_k2 = snode->right->pic[k]->trait;
      for (j = 0; j < dim[k]; ++j)
      {
        ic->contrast[j] = m_k1[j] - m_k2[j];                       // x_k
        ic->trait[j] = (v_k2*m_k1[j] + v_k1*m_k2[j]) /(v_k1+v_k2); // m_k'
      }
    }
  }
}

void pic_init(stree_t * stree, trait_t ** trait_list, long n_part)
{
  assert(stree != NULL && trait_list != NULL);
  assert(n_part > 0);
  
  int i, j, k, l, found;
  snode_t * snode;
  contrast_t * ic;
  trait_t * m;
  
  /* tip node(_k) has trait values (m_k) and branch length (v_k) as is */
  for (i = 0; i < stree->tip_count; ++i)
  {
    snode = stree->nodes[i];
    snode->pic = (contrast_t **)xmalloc(n_part*sizeof(contrast_t *));
    
    /* loop over the trait partitions */
    for (k = 0; k < n_part; ++k)
    {
      m = trait_list[k];
      ic = snode->pic[k] = (contrast_t *)xmalloc(sizeof(contrast_t));

      /* TODO: how to deal with zero branch length */
      ic->brlen = snode->parent->tau - snode->tau;
      
      /* find and fill the trait values for this tip node */
      found = 0;
      for (l = 0; l < m->count; ++l)
      {
        if (strncmp(snode->label, m->label[l], LABEL_LEN) == 0)
        {
          /* copy the trait values over */
          ic->trait = (double *)xmalloc((m->length)*sizeof(double));
          for (j = 0; j < m->length; ++j)
            ic->trait[j] = m->trait[l][j];
          
          found++;
          break;
        }
      }
      if (found == 0)
      {
        fatal("Species name %s not found in partition %d", snode->label, k+1);
      }
    }
  }
  
  stree->trait_dim = (int *)xmalloc(n_part*sizeof(int));
  for (k = 0; k < n_part; ++k)
    stree->trait_dim[k] = trait_list[k]->length;

  /* internal node: compute contrasts (x_k), ancestral traits (m_k'), and
                    transformed branch length (v_k') */
  pic_init_recursive(stree->root, n_part, stree->trait_dim);
  
#if(DEBUG_Chi)
  for (k = 0; k < n_part; ++k)
  {
    printf("Trait partition %d\n", k+1);
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    {
      snode = stree->nodes[i];
      printf("%s\t", snode->label);
      printf("v'=%lf, ", snode->pic[k]->brlen);
      printf("m':");
      for (j = 0; j < stree->trait_dim[k]; ++j)
        printf(" %lf", snode->pic[k]->trait[j]);
      printf(", ");
      if (i >= stree->tip_count)
      {
        printf("x:");
        for (j = 0; j < stree->trait_dim[k]; ++j)
          printf(" %lf", snode->pic[k]->contrast[j]);
      }
      printf("\n");
    }
  }
#endif
}
