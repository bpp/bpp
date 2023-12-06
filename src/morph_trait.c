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

static int parse_trait(FILE * fp, double * trait, int len)
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
  
  trait_t * m = (trait_t *)xcalloc(1,sizeof(trait_t));
  
  /* read header */
  if (!parse_header(fp, &(m->count), &(m->length)))
    return NULL;

  m->label = (char **)xcalloc((m->count),sizeof(char *));
  m->trait = (double **)xcalloc((m->count),sizeof(double *));
  for (i = 0; i < m->count; ++i)
  {
    m->trait[i] = (double *)xmalloc((m->length) * sizeof(double));
    m->label[i] = (char *)xmalloc((LABEL_LEN+1) * sizeof(char));
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
    if (!parse_trait(fp, m->trait[i], m->length))
    {
      fprintf(stderr, "Failed to read traits of species %d\n", i+1);
      return NULL;
    }
  }
  
#if(DEBUG_Chi)
  printf("\n%d %d\n", m->count, m->length);
  for (i = 0; i < m->count; ++i) {
    printf("%s\t", m->label[i]);
    for (j = 0; j < m->length; ++j) {
      printf("%lf\t", m->trait[i][j]);
    }
    printf("\n");
  }
  printf("\n");
#endif
  
  m->dtype = BPP_DATA_TRAIT;
  
  return m;
}

trait_t ** parse_traitfile(const char * traitfile, long * count)
{
  FILE * fp;
  long m_slotalloc = 10;
  long m_maxcount = 0;
  int c;
  
  fp = xopen(traitfile,"r");

  *count = 0;

  trait_t ** mm = (trait_t **)xmalloc(m_maxcount*sizeof(trait_t *));

  while ((c=fgetc(fp)) != EOF) {
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
