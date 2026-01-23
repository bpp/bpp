//
//  Created by Chi Zhang on 2023/12/5.
//

#include "bpp.h"

#define LABEL_LEN 99

#define DEBUG_Chi 1
#define DEBUG_Morph_Matrix 1
#define DEBUG_Morph_BM_M   1
#define BM_AC     1
#define BM_Mitov  0

/* matrix related operations */
static int mat_sub(double *A, double *Asub, int n, int m, int *r, int *c);
static int mat_scale(double *A, double b, double *C, int n, int m);
static int mat_add(double *A, double *B, double *C, int n, int m);
static int mat_multi(double *A, double *B, double *C, int n, int m, int k);
static int mat_trans(double *A, double *At, int n, int m);
static int mat_decom_chol(double *A, double *L, int n);
static int mat_inv(double *L, double *Ainv, int n);
static double mat_logdet(double *L, int n);
static int mat_issym(double *A, int n);

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
    return -1;
  }
  else if (*nrow <= 0) {
    fprintf(stderr, "Number of species must be >0(%d)\n", *nrow);
    return 1;
  }
  else if (*ncol <= 0) {
    fprintf(stderr, "Number of traits must be >0 (%d)\n", *ncol);
    return 1;
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
    return 1;
  }

  return 0;
}

static int parse_label(FILE * fp, char * name, int len)
{
  int c, j;
  
  c = get_nb_char(fp);
  if (c == EOF)
    return -1;
  
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

  return 0;
}

static int parse_value_c(FILE * fp, double * trait, int len)
{
  int c, j;
  
  for (j = 0; j < len; ++j)
  {
    c = get_nb_char(fp);
    
    if (c == EOF)
    {
      return -1;
    }
    else if (c == '-')
    {
      c = fgetc(fp);
      if (isspace(c) || c == EOF)
        trait[j] = NAN;    // inapplicable value
      else
      {
        ungetc(c, fp);
        if (fscanf(fp, "%lf", &trait[j]) != 1)
          return 1;
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
      if (c != 'A' && c != 'a')
        return 1;
      else
      {
        c = fgetc(fp);
        if (isspace(c) || c == EOF)
          trait[j] = INFINITY;  // missing value
        else if (c == 'N' || c == 'n')
          trait[j] = NAN;  // inapplicable value
        else
          return 1;
      }
    }
    else
    {
      ungetc(c, fp);           // positive value
      if (fscanf(fp, "%lf", &trait[j]) != 1)
        return 1;
    }
  }
  assert(j == len);

  return 0;
}

static int state_bin(int x)
{
  int std_bin[] = {1,2,4,8,16,32,64,128,256,512};
  
  if (x >= 0 && x <= 9)
    return std_bin[x];
  else
  {
    fprintf(stderr, "Unsupported trait value (%d)\n", x);
    return -1;
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
      return -1;
    
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
        return 1;
      }
    }
    else
    {
      fprintf(stderr, "Unrecognized trait value (%c) at %d\n", c, j+1);
      return 1;
    }
  }
  
  return 0;
}

static void parse_comment(FILE * fp)
{
  int c;
  
  while ((c=get_nb_char(fp)) == '#')
    while ((c=fgetc(fp)) != '\n'); //eat the line

  ungetc(c, fp);
}

static int parse_matrix_c(FILE * fp, double * mat, int n, int m)
{
  int i, j;
  for (i = 0; i < n; ++i)
    for (j = 0; j < m; ++j)
      if (fscanf(fp, "%lf", &mat[i*m + j]) != 1)
        return 1;
  return 0;
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
  
  if (morph->matRs)
    free(morph->matRs);
  
  free(morph);
}

static morph_t * morph_parse(FILE * fp)
{
  int  i, j, nchar, c;
  char tmpstr[10];
  
  morph_t * morph = (morph_t *)xmalloc(sizeof(morph_t));
  
  /* skip comments, assuming they are at the beginning of the trait block */
  parse_comment(fp);
  
  /* read header */
  if (parse_header(fp, &(morph->ntaxa), &(morph->length), &(morph->dtype)))
  {
    fprintf(stderr, "Error in header\n");
    return NULL;
  }
  nchar = morph->length;
  
  /* allocate space */
  morph->label = (char **)xmalloc((morph->ntaxa)*sizeof(char *));
  for (i = 0; i < morph->ntaxa; ++i)
    morph->label[i] = (char *)xmalloc((LABEL_LEN+1)*sizeof(char));

  if (morph->dtype == BPP_DATA_CONT)
  {
    morph->conti = (double **)xmalloc((morph->ntaxa)*sizeof(double *));
    for (i = 0; i < morph->ntaxa; ++i)
      morph->conti[i] = (double *)xmalloc(nchar*sizeof(double));
  }
  else // morph->dtype == BPP_DATA_DISC
  {
    morph->discr = (int **)xmalloc((morph->ntaxa)*sizeof(int *));
    for (i = 0; i < morph->ntaxa; ++i)
      morph->discr[i] = (int *)xmalloc(nchar*sizeof(int));
  }
  
  /* read morphological traits of each species,
     assuming each line has a label followed by numbers.
     trait matrix has dimension ntaxa * nchar */
  for (i = 0; i < morph->ntaxa; ++i)
  {
    /* read the label (species name) */
    if (parse_label(fp, morph->label[i], LABEL_LEN+1))
    {
      fprintf(stderr, "Failed to read label of species %d\n", i+1);
      morph_destroy(morph);
      return NULL;
    }
    /* read the trait values of this species */
    if (morph->dtype == BPP_DATA_CONT &&
        parse_value_c(fp, morph->conti[i], nchar))
    {
      fprintf(stderr, "Failed to read traits of species %s\n", morph->label[i]);
      morph_destroy(morph);
      return NULL;
    }
    else if (morph->dtype == BPP_DATA_DISC &&
             parse_value_d(fp, morph->discr[i], nchar))
    {
      fprintf(stderr, "Failed to read traits of species %s\n", morph->label[i]);
      morph_destroy(morph);
      return NULL;
    }
  }

  if (morph->dtype == BPP_DATA_CONT)
  {
    c = get_nb_char(fp);
    if (c == 'R')
    {
      /* read correlation matrix (R*) */
      fgetc(fp);  // skip '*'
      morph->matRs = (double *)xmalloc(nchar*nchar*sizeof(double));
      if (parse_matrix_c(fp, morph->matRs, nchar, nchar))
      {
        fprintf(stderr, "Error reading correlation matrix (R*)\n");
        morph_destroy(morph);
        return NULL;
      }
    }
    else if (c == 'l')
    {
      /* read log determinant of R* */
      morph->matRs = (double *)xmalloc(sizeof(double));
      if (fscanf(fp, "%s %lf", tmpstr, morph->matRs) != 2)
      {
        fprintf(stderr, "Error reading log determinant of R*\n");
        morph_destroy(morph);
        return NULL;
      }
    }
    else
    {
      ungetc(c, fp);
      morph->matRs = NULL;  // assume identity matrix
    }
  }

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
  if (c == 'R') {
    printf("\nR* = \n");
    for (i = 0; i < morph->length; ++i) {
      for (j = 0; j < morph->length; ++j)
        printf("%+.3lf\t", morph->matRs[i * morph->length + j]);
      printf("\n");
    }
  }
  else if (c == 'l') {
    printf("logdet(R*) = %.3lf\n", morph->matRs[0]);
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

    pmorph[*count] = morph_parse(fp);
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

static void bm_update_vxm(int idx, snode_t * snode, stree_t * stree)
{
  /* Alvarez-Carretero et al. 2019. p.969. */

  int j;
  double v_k, v_k1, v_k2, *m_k1, *m_k2;

  if (snode->left && snode->right)  /* internal node */
  {
    bm_update_vxm(idx, snode->left, stree);
    bm_update_vxm(idx, snode->right, stree);
    
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
    /* the trait matrix has been standardized so that all characters have
       the same variance and the population noise has unit variance. */
    snode->trait[idx]->brlen = v_k + 1.0;
  }

#ifdef DEBUG_Morph_BM_A
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

static int bm_init_Rs_Phi(stree_t * stree, morph_t ** morph_list)
{
  int n, i, j, nchar;
  morph_t * morph;

  for (n = 0; n < stree->trait_count; ++n)
  { 
    morph = morph_list[n];
    
    if (morph->dtype != BPP_DATA_CONT)
      continue;
    
    nchar = morph->length;

    if (stree->trait_missing[n] || BM_Mitov)
    {
      stree->trait_Rs[n] = (double *)xmalloc(nchar*nchar*sizeof(double));
      if (morph->matRs == NULL)
      {
        /* set up identity matrix for R* */
        for (i = 0; i < nchar; ++i)
          for (j = 0; j < nchar; ++j)
            if (i == j)
              stree->trait_Rs[n][i*nchar + j] = 1.0;
            else
              stree->trait_Rs[n][i*nchar + j] = 0.0;
      }
      else
      {
        /* copy the R* matrix over */
        memcpy(stree->trait_Rs[n], morph->matRs, nchar*nchar*sizeof(double));

        /* check whether R* is symmetric */
        if (mat_issym(stree->trait_Rs[n], nchar))
          fatal("Error: correlation matrix R* is not symmetric");
      }

      /* set up identity matrix for Phi */
      stree->trait_Phi[n] = (double *)xmalloc(nchar*nchar*sizeof(double));
      for (i = 0; i < nchar; ++i)
        for (j = 0; j < nchar; ++j)
          if (i == j)
            stree->trait_Phi[n][i*nchar + j] = 1.0;
          else  
            stree->trait_Phi[n][i*nchar + j] = 0.0;
    }
    else  // no missing data
    {
      if (morph->matRs == NULL)
        stree->trait_ldetRs[n] = 0.0;  // identity matrix
      else
        stree->trait_ldetRs[n] = morph->matRs[0];
    }
  }

  return 0;
}

static void bm_ACEf_Lmr(int idx, snode_t * snode, stree_t * stree,
                        double * L_i, double * m_i, double * r_i)
{
  /* Mitov et al. 2020; BM model (Eq. 2, 10, 11) */

  int    nchar, k_i, k_j, *act, *act_p;
  double t, *R, *I, *A, *C, *E, *L, *m, r,  *V, *T, *x, *y;

  if (snode == stree->root)
    return;
  
  /* find stratch space */
  R = snode->trait[idx]->tranprob[0];
  I = snode->trait[idx]->tranprob[1];
  A = snode->trait[idx]->tranprob[2];
  C = snode->trait[idx]->tranprob[3];
  E = snode->trait[idx]->tranprob[4];
  V = snode->trait[idx]->tranprob[5];
  T = snode->trait[idx]->tranprob[6];
  x = snode->trait[idx]->tranprob[7];
  y = snode->trait[idx]->tranprob[8];

  /* total number of characters */
  nchar = stree->trait_dim[idx];

  /* active coordinates */
  act = snode->trait[idx]->active;
  act_p = snode->parent->trait[idx]->active;
  k_i = act[nchar];   // number of active coordinates at node i
  k_j = act_p[nchar]; // number of active coordinates at parent j

  t = (snode->parent->tau - snode->tau) * snode->trait[idx]->brate;
  /* the trait matrix has been standardized so that all characters have
     the same variance and the population noise has unit variance. */
  if (snode->left == NULL)
    t += 1.0;
  snode->trait[idx]->brlen = t;

  /* Rs is the linear shrinkage estimate of the correlation matrix R,
     which is input along with the morphological data */
  mat_sub(stree->trait_Rs[idx], R, nchar, nchar, act, act);
  mat_scale(R, t, V, k_i, k_i);         // V = t*R
  mat_decom_chol(V, T, k_i);
  mat_inv(T, V, k_i);                   // V: inv(V)

  t = -0.5 * mat_logdet(T, k_i);        // reuse t

  /* E = t(Phi) * inv(V) */
  mat_sub(stree->trait_Phi[idx], I, nchar, nchar, act, act_p);
  mat_trans(I, T, k_i, k_j);            // T: t(Phi)
  mat_multi(T, V, E, k_j, k_i, k_i);

  /* C = -0.5 * E * Phi */
  mat_multi(E, I, C, k_j, k_i, k_j);
  mat_scale(C, -0.5, C, k_j, k_j);

  /* A = -0.5 * inv(V) */
  mat_scale(V, -0.5, A, k_i, k_i);

  if (snode->left == NULL) // tip
  {
    m = snode->trait[idx]->state_m;
    mat_sub(m, x, 1, nchar, NULL, act);

    /* L_i = C */
    mat_scale(C, 1.0, L_i, k_j, k_j);

    /* m_i = E * x */
    mat_multi(E, x, m_i, k_j, k_i, 1);
    
    /* r_i = t(x) * A * x + f */
    mat_multi(x, A, y, 1, k_i, k_i);
    mat_multi(y, x, r_i, 1, k_i, 1);
    *(r_i) += t - 0.5 * k_i * log(2.0 * M_PI);
  }
  else  // internal node
  {
    L = snode->trait[idx]->glinv_L;
    m = snode->trait[idx]->state_m;
    r = snode->trait[idx]->glinv_r;

    mat_add(A, L, V, k_i, k_i);         // V: A+L
    mat_scale(V, -2.0, A, k_i, k_i);    // A: -2*(A+L)
    mat_decom_chol(A, T, k_i);
    t += -0.5 * mat_logdet(T, k_i);

    mat_decom_chol(V, T, k_i);
    mat_inv(T, V, k_i);                 // V: inv(A+L)

    /* L_i = C - 0.25 * E * inv(A+L) * t(E) */
    mat_multi(E, V, A, k_j, k_i, k_i);  // A: E * inv(A+L)
    mat_trans(E, T, k_j, k_i);          // T: t(E)
    mat_multi(A, T, E, k_j, k_i, k_j);  // E: E * inv(A+L) * t(E)
    mat_scale(E, -0.25, L_i, k_j, k_j);
    mat_add(C, L_i, L_i, k_j, k_j);

    /* m_i = -0.5 * E * inv(A+L) * m */
    mat_multi(A, m, m_i, k_j, k_i, 1);
    mat_scale(m_i, -0.5, m_i, k_j, 1);

    /* r_i = r + f + 0.5 * k * log(2pi)
             - 0.5 * logdet(-2 * (A+L))
             - 0.25 * t(m) * inv(A+L) * m */
    mat_multi(m, V, y, 1, k_i, k_i);
    mat_multi(y, m, r_i, 1, k_i, 1);
    *(r_i) = *(r_i) * (-0.25) + r + t;
  }
}

static void bm_update_Lmr(int idx, snode_t * snode, stree_t * stree)
{
  /* Mitov et al. 2020; BM model (Theorem 2) */

  int    nchar, n_act, j;
  double *L1, *L2, *m1, *m2, r1, r2;

  if (snode->left && snode->right)  // internal node
  {
    bm_update_Lmr(idx, snode->left, stree);
    bm_update_Lmr(idx, snode->right, stree);

    /* total number of characters */
    nchar = stree->trait_dim[idx];

    /* active coordinates */
    for (n_act = 0, j = 0; j < nchar; ++j)
    {
      if (snode->left->trait[idx]->active[j] == -1 &&
          snode->right->trait[idx]->active[j] == -1)
      {
        snode->trait[idx]->active[j] = -1;
      }
      else {
        snode->trait[idx]->active[j] = 1;
        n_act += 1;
      }
    }
    snode->trait[idx]->active[nchar] = n_act;
  
    /* find stratch space */
    L1 = snode->trait[idx]->tranprob[0];
    L2 = snode->trait[idx]->tranprob[1];
    m1 = snode->trait[idx]->tranprob[7];
    m2 = snode->trait[idx]->tranprob[8];

    bm_ACEf_Lmr(idx, snode->left,  stree, L1, m1, &r1);
    bm_ACEf_Lmr(idx, snode->right, stree, L2, m2, &r2);

    /* store L, m, r at this node */
    mat_add(L1, L2, snode->trait[idx]->glinv_L, n_act, n_act);
    for (j = 0; j < n_act; ++j)
      snode->trait[idx]->state_m[j] = m1[j] + m2[j];
    snode->trait[idx]->glinv_r = r1 + r2;
  }
}

static void mk_trprob(double ** p, double v, int max_state)
{
  int k;
  
  /* Mk model, Lewis 2001 */
  for (k = 2; k <= max_state; ++k)
  {
    p[k-2][0] = 1.0/k + (k-1.0)/k * exp(-v * k/(k-1.0)); // no change
    p[k-2][1] = 1.0/k - 1.0/k * exp(-v * k/(k-1.0));     // change
  }
}

static void mk_update_cp(int idx, snode_t * snode, stree_t * stree)
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
  mk_trprob(snode->trait[idx]->tranprob, v, max_state);

  /* pruning algorithm */
  /* TODO: this can be optimized by grouping characters by patterns */
  if (snode->left && snode->right)  /* internal node */
  {
    mk_update_cp(idx, snode->left, stree);
    mk_update_cp(idx, snode->right, stree);

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
  stree->nodes[3]->tau = 0.13;
  stree->nodes[4]->tau = 0.08;
#endif

  if (stree->trait_type[idx] == BPP_DATA_DISC)
    mk_update_cp(idx, stree->root, stree);
  else if (stree->trait_missing[idx] || BM_Mitov)
    bm_update_Lmr(idx, stree->root, stree);
  else
    bm_update_vxm(idx, stree->root, stree);
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
  int n, i, j, k, l, nchar, state, max_state, n_act;
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
          if (morph->dtype == BPP_DATA_CONT)
            memcpy(snode->trait[n]->state_m, morph->conti[l], 
                             (morph->length) * sizeof(double));
          else
            memcpy(snode->trait[n]->state_d, morph->discr[l],
                                (morph->length) * sizeof(int));
          break;
        }
      }
      if (l == morph->ntaxa)
      {
        fprintf(stderr, "Species name %s not found in partition %d\n",
                snode->label, n+1);
        return 1;
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
          return 1;
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
      /* update the vector for the active coordinates 
         while checking missing states */
      for (i = 0; i < stree->tip_count; ++i)
      {
        snode = stree->nodes[i];
        for (n_act = 0, j = 0; j < nchar; ++j)
        {
          value = snode->trait[n]->state_m[j];
          if (isnan(value))  // inapplicable value
          {
            stree->trait_missing[n] = 2;
            snode->trait[n]->active[j] = -1;
          }
          else if (isinf(value))  // missing value
          {
            if (stree->trait_missing[n] != 2)
              stree->trait_missing[n] = 1;
            snode->trait[n]->active[j] = 0;
          }
          else {
            snode->trait[n]->active[j] = 1;
            n_act += 1;
          }
        }
        snode->trait[n]->active[nchar] = n_act;
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

  return 0;
}

static void trait_alloc_mem(stree_t * stree, morph_t ** morph_list, int n_part)
{
  int n, i, j, nchar;
  trait_t * trait;
  
  stree->trait_dim = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_type = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_missing = (int *)xcalloc(n_part, sizeof(int));
  stree->trait_ldetRs = (double *)xcalloc(n_part, sizeof(double));
  
  stree->trait_nstate = (int **)xcalloc(n_part, sizeof(int *));
  stree->trait_Rs =  (double **)xcalloc(n_part, sizeof(double *));
  stree->trait_Phi = (double **)xcalloc(n_part, sizeof(double *));
  for (n = 0; n < n_part; ++n)
  {
    nchar = morph_list[n]->length;
    /* the number of states of each character
       use the last element to store the max number of states */
    if (morph_list[n]->dtype == BPP_DATA_DISC)
      stree->trait_nstate[n] = (int *)xcalloc(nchar +1, sizeof(int));
    
    /* allocate trait_Rs[n] and trait_Phi[n] later if needed */
  }
  
  for (i = 0; i < stree->tip_count+stree->inner_count; ++i)
  {
    stree->nodes[i]->trait = (trait_t **)xmalloc(n_part*sizeof(trait_t *));
    for (n = 0; n < n_part; ++n)
    {
      stree->nodes[i]->trait[n] = (trait_t *)xmalloc(sizeof(trait_t));
      trait = stree->nodes[i]->trait[n];
      
      nchar = morph_list[n]->length;
      if (morph_list[n]->dtype == BPP_DATA_CONT)
      {
        trait->state_m = (double *)xcalloc(nchar, sizeof(double));
        trait->active = (int *)xcalloc(nchar +1, sizeof(int));
        if (i >= stree->tip_count)
        {
          trait->contrast = (double *)xcalloc(nchar, sizeof(double));
          trait->glinv_L = (double *)xmalloc(nchar*nchar * sizeof(double));
        }

        /* scratch space to store L_i, m_i, etc */
        trait->tranprob = (double **)xmalloc(9 * sizeof(double *));
        for (j = 0; j < 7; ++j)
          trait->tranprob[j] = (double *)xmalloc(nchar*nchar * sizeof(double));
        for (j = 7; j < 9; ++j)
          trait->tranprob[j] = (double *)xmalloc(nchar * sizeof(double));
      }
      else
      {
        trait->state_d = (int *)xcalloc(nchar, sizeof(int));
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
          if (trait->glinv_L)
            free(trait->glinv_L);  
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
  if (stree->trait_ldetRs)
    free(stree->trait_ldetRs);

  if (stree->trait_nstate)
  {
    for (n = 0; n < stree->trait_count; ++n)
    {
      if (stree->trait_nstate[n])
        free(stree->trait_nstate[n]);
    }
    free(stree->trait_nstate);
  }
  if (stree->trait_Rs)
  {
    for (n = 0; n < stree->trait_count; ++n)
    {
      if (stree->trait_Rs[n])
      {
        free(stree->trait_Rs[n]);
        free(stree->trait_Phi[n]);
      }
    }
    free(stree->trait_Rs);
    free(stree->trait_Phi);
  }

  if (stree->trait_logl)
    free(stree->trait_logl);
  if (stree->trait_old_logl)
    free(stree->trait_old_logl);
  if (stree->trait_logpr)
    free(stree->trait_logpr);
  if (stree->trait_old_logpr)
    free(stree->trait_old_logpr);
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
  }
  
  /* fill the trait values for the tip nodes; */
  if (trait_fill_tip(stree, morph_list))
  {
    trait_destroy(stree);
    fatal("Error filling traits");
  }
  
  /* for continuous traits, set up R and Phi */
  bm_init_Rs_Phi(stree, morph_list);

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

static double loglikelihood_BM_AC(int idx, stree_t * stree)
{
  int i, j, p;
  double v_k1, v_k2, zz, ldetRs, logl;
  snode_t * snode;

  /* log determinant of shrinkage estimate of the correlation matrix,
      i.e. log(det(R*)) */
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

#ifdef DEBUG_Morph_BM_A
  printf("part%d: cur log(like)=%lf, old log(like)=%lf\n\n", idx+1,
         stree->trait_logl[idx], stree->trait_old_logl[idx]);
#endif

  return logl;
}

static double loglikelihood_BM_Mitov(int idx, stree_t * stree)
{
  int nchar, k_0;
  double logl, *L0, *m0, r0, *x0, *T, *L0_1, *y, z;
  snode_t * root = stree->root;

  /* find stratch space */
  T    = root->trait[idx]->tranprob[0];
  L0_1 = root->trait[idx]->tranprob[1];
  x0   = root->trait[idx]->tranprob[7];
  y    = root->trait[idx]->tranprob[8];

  /* number of active coordinates */
  nchar = stree->trait_dim[idx];
  k_0 = root->trait[idx]->active[nchar];

  /* Mitov et al. 2020; Eq. S2
     x0 = -0.5 * inv(L0) * m0 */
  L0 = root->trait[idx]->glinv_L;
  m0 = root->trait[idx]->state_m;
  mat_decom_chol(L0, T, k_0);
  mat_inv(T, L0_1, k_0);  // L0_1: inv(L0)
  mat_multi(L0_1, m0, y, k_0, k_0, 1);
  mat_scale(y, -0.5, x0, k_0, 1);  

  /* logl = t(x0) * L0 * x0 + t(x0) * m0 + r0 */
  mat_multi(x0, L0, y, 1, k_0, k_0);
  mat_multi(y, x0, &z, 1, k_0, 1);
  logl = z;

  mat_multi(x0, m0, &z, 1, k_0, 1);
  r0 = root->trait[idx]->glinv_r;
  logl += z + r0;

  stree->trait_logl[idx] = logl;

#ifdef DEBUG_Morph_BM_M
  printf("part%d: cur log(like)=%lf, old log(like)=%lf\n\n", idx+1,
         stree->trait_logl[idx], stree->trait_old_logl[idx]);
#endif

  return logl;
}

static double loglikelihood_Mkv(int idx, stree_t * stree)
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
    return loglikelihood_Mkv(idx, stree);
  }
  else if (stree->trait_missing[idx] || BM_Mitov)
  {
    /* BM model with missing data; Mitov et al. 2020 */
    return loglikelihood_BM_Mitov(idx, stree);
  }
  else
  {
    /* BM model without missing data; Alvarez-Carretero et al. 2019 */
    return loglikelihood_BM_AC(idx, stree);
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


/* submatrix of A[n*m] by active coordinates */
static int mat_sub(double *A, double *Asub, int n, int m, int *r, int *c)
{
  int i, j, a, b, ncol;
  
  /* get total active columns */
  ncol = c[m];

  /* extract submatrix using active rows and columns */
  for (a = 0, i = 0; i < n; ++i)
  {
    if (n > 1 && r[i] != 1) // n=1 if A is a vector
      continue;
    for (b = 0, j = 0; j < m; ++j)
    {
      if (c[j] != 1) continue;
      Asub[a * ncol + b] = A[i * m + j];
      b++;
    }
    a++;
  }
  
  return 0;
}

/* C = b*A; b is a scalar */
static int mat_scale(double *A, double b, double *C, int n, int m)
{
  int i, j;
  
  for (i = 0; i < n; ++i)
    for (j = 0; j < m; ++j)
      C[i * m + j] = b * A[i * m + j];
  
  return 0;
}

/* C = A+B; A[n*m], B[n*m], C[n*m] */
static int mat_add(double *A, double *B, double *C, int n, int m)
{
  int i, j;
  
  for (i = 0; i < n; ++i)
    for (j = 0; j < m; ++j)
      C[i * m + j] = A[i * m + j] + B[i * m + j];
  
  return 0;
}

/* C = A*B; A[n*m], B[m*k], C[n*k] */
static int mat_multi(double *A, double *B, double *C, int n, int m, int k)
{
  int i, j, r;
  double sum;
  
  for (i = 0; i < n; ++i)
    for (j = 0; j < k; ++j)
    {
      for (sum = 0, r = 0; r < m; ++r)
        sum += A[i * m + r] * B[r * k + j];
      C[i * k + j] = sum;
    }

  return 0;
}

/* A'; A[n*m], At[m*n] */
static int mat_trans(double *A, double *At, int n, int m)
{
  int i, j;
  
  for (i = 0; i < n; ++i)
    for (j = 0; j < m; ++j)
      At[j * n + i] = A[i * m + j];
  
  return 0;
}

/* Cholesky decomposition: A = LL',
   where A is symmetrical and positive definite, and L is lower triangular */
static int mat_decom_chol(double *A, double *L, int n)
{
  int i, j, k;
  double sum;
  
  for (i = 0; i < n; ++i)
    for (j = 0; j <= i; ++j)
    {
      sum = A[i * n + j];
      for (k = 0; k < j; ++k)
        sum -= L[i * n + k] * L[j * n + k];
      if (i == j)
      {
        if (sum <= 0.0) return -1;  // not positive definite
        L[i * n + j] = sqrt(sum);
      }
      else
      {
        L[i * n + j] = sum / L[j * n + j];
      }
    }
  
  /* set upper triangular part to zero */
  for (i = 0; i < n; ++i)
    for (j = i + 1; j < n; ++j)
      L[i * n + j] = 0.0;
  
  return 0;
}

/* inverse A using Cholesky decomposition L */
static int mat_inv(double *L, double *Ainv, int n)
{
  int i, j, k;
  double *Linv, sum, *Linv_t;
  
  /* store inv(L) separately to preserve L */
  Linv   = (double *)xmalloc(n*n*sizeof(double));
  Linv_t = (double *)xmalloc(n*n*sizeof(double));

  for (i = 0; i < n; ++i)
  {
    Linv[i * n + i] = 1.0 / L[i * n + i];
    for (j = 0; j < i; ++j)
    {
      for (sum = 0, k = j; k < i; ++k)
        sum -= L[i * n + k] * Linv[k * n + j];
      Linv[i * n + j] = sum / L[i * n + i];
    }
  }

  /* inv(A) = t(inv(L)) * inv(L) */
  mat_trans(Linv, Linv_t, n, n);
  mat_multi(Linv_t, Linv, Ainv, n, n, n);

  free(Linv);
  free(Linv_t);
  
  return 0;
}

/* determinant A using Cholesky decomposition L */
static double mat_logdet(double *L, int n)
{
  int i;
  double logdet;

  /* logdet(A) = 2 * sum(log(L[i,i])) for A = LL' */
  logdet = 0.0;
  for (i = 0; i < n; ++i)
    logdet += log(L[i * n + i]);
  logdet *= 2.0;

  return logdet;
}

/* check whether A is symmetric */
static int mat_issym(double *A, int n)
{
  int i, j;
  double tol = 1e-5;

  for (i = 0; i < n; ++i)
    for (j = i + 1; j < n; ++j)
      if (fabs(A[i * n + j] - A[j * n + i]) > tol)
        return 1;

  return 0;
}
