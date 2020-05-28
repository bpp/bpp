/*
    Copyright (C) 2016-2019 Tomas Flouri and Ziheng Yang

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

#define GET_HINDEX(t,p) (((node_is_mirror((p)) ? \
                          (p)->node_index : (p)->hybrid->node_index)) - \
                        ((t)->tip_count+(t)->inner_count))

#define DNA_QRATES_COUNT        6
#define DNA_STATES_COUNT        4


/* TODO: Simulation is at the moment single-threaded */
static const long thread_index_zero = 0;

static long * g_order = NULL;

static char charmap_nt_tcag[16] =
{
  '\0', 'T', 'C', 'Y', 'A', 'W', 'M', 'H',
   'G', 'K', 'S', 'B', 'R', 'D', 'V', 'X'
};

static char * cb_serialize_branch(const snode_t * node)
{
  char * s = NULL;

  /* inner node */
  if (node->left)
  {
    if (node->parent)
    {
      if (opt_est_theta && node->theta > 0)
        xasprintf(&s, " #%f: %f", node->theta, node->parent->tau - node->tau);
      else
        xasprintf(&s, ": %f", node->parent->tau - node->tau);
    }
    else
    {
      if (opt_est_theta && node->theta > 0)
        xasprintf(&s, " #%f", node->theta);
      else
      {
        /* replacement to keep GCC warning from showing up when executing
           the following:
           
           xasprintf(&s, "");
        */

        s = (char *)xmalloc(sizeof(char));
        *s = 0;

      }
    }
      
  }
  else
  {
    if (opt_est_theta && node->theta > 0)
      xasprintf(&s, "%s #%f: %f",
               node->label, node->theta, node->parent->tau - node->tau);
    else
      xasprintf(&s, "%s: %f", node->label, node->parent->tau - node->tau);
  }
    

  return s;
}

static void print_settings(stree_t * stree)
{
  long i;

  fprintf(stdout,"%d species:", stree->tip_count);
  for (i = 0; i < stree->tip_count; ++i)
  {
    fprintf(stdout, " %s (%ld)", stree->nodes[i]->label, opt_sp_seqcount[i]);
  }
  fprintf(stdout,"\n");
  if (opt_msci)
  {
    printf("  NETWORK\n");
    print_network_table(stree);
  }
  else
  {
    if (stree->tip_count == 1)
    {
      printf("  %s\n", stree->nodes[0]->label);
    }
    else
    {
      char * newick = stree_export_newick(stree->root, cb_serialize_branch);
      printf("  %s\n", newick);
      free(newick);
    }
  }

  stree_show_pptable(stree, BPP_TRUE);
}

static int MultiNomialAliasSetTable(int ncat, double * prob, double * F, int * L)
{
   /* This sets up the tables F and L for the alias algorithm for generating samples from the
      multinomial distribution MN(ncat, p) (Walker 1974; Kronmal & Peterson 1979).

      F[i] has cutoff probabilities, L[i] has aliases.
      I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.

      Should perhaps check whether prob[] sums to 1.
   */
   signed char *I = (signed char *)xmalloc((size_t)ncat * sizeof(signed char));
   int i, j, k, nsmall;

   for (i = 0; i < ncat; i++)  L[i] = -9;
   for (i = 0; i < ncat; i++)  F[i] = ncat*prob[i];
   for (i = 0, nsmall = 0; i < ncat; i++) {
      if (F[i] >= 1)  I[i] = 1;
      else { I[i] = -1; nsmall++; }
   }
   for (i = 0; nsmall > 0; i++) {
      for (j = 0; j < ncat; j++)  if (I[j] == -1) break;
      for (k = 0; k < ncat; k++)  if (I[k] == 1)  break;
      if (k == ncat)  break;

      L[j] = k;
      F[k] -= 1 - F[j];
      if (F[k] < 1) { I[k] = -1; nsmall++; }
      I[j] = 0;  nsmall--;
   }

   free(I);
   return(0);
}

static int MultiNomialAlias(int n, int ncat, double * F, int * L, int * nobs)
{
   /* This generates multinomial samples using the F and L tables set up before,
      using the alias algorithm (Walker 1974; Kronmal & Peterson 1979).

      F[i] has cutoff probabilities, L[i] has aliases.
      I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.
   */
   int i, k;
   double r;

   for (i = 0; i < ncat; i++)  nobs[i] = 0;
   for (i = 0; i < n; i++) {
      r = legacy_rndu(thread_index_zero)*ncat;
      k = (int)r;
      r -= k;
      if (r <= F[k]) nobs[k]++;
      else        nobs[L[k]]++;
   }
   return (0);
}

static void process_diploid(long species_count)
{
  long i;

  if (!opt_diploid)
    opt_diploid = (long *)xcalloc((size_t)species_count, sizeof(long));

  /* double the number of sequences for species indicated as diploid */
  for (i = 0; i < species_count; ++i)
  {
    if (opt_diploid[i])
      opt_sp_seqcount[i] *= 2;
  }

  opt_cleandata = 1;
  for (i = 0; i < species_count; ++i)
    if (opt_diploid[i])
    {
      opt_cleandata = 0;
      break;
    }
}

static void process_subst_model()
{
  long i;

  if (opt_model == BPP_DNA_MODEL_JC69)
  {
    fprintf(stdout, "Substitution model: JC69\n");
    opt_qrates_fixed = 1;
    opt_basefreqs_fixed = 1;

    if (!opt_qrates_params)
      opt_qrates_params = (double *)xmalloc(DNA_QRATES_COUNT*sizeof(double));
      
    for (i = 0; i < 6; ++i)
      opt_qrates_params[i] = 0;
  }
  else
  {
    /* GTR model */

    if (opt_qrates_fixed)
    {
      fprintf(stdout, "Substitution model: GTR with fixed rates\n");
      opt_qrates_fixed = 1;
        for (i = 0; i < 6; ++i)
          opt_qrates_params[i] /= opt_qrates_params[5];
    }
    else
    {
      fprintf(stdout, "Substitution model: GTR with generated rates from Dirichlet\n");
    }
    fprintf(stdout, "                   ");
    for (i = 0; i < 6; ++i)
      fprintf(stdout, " %f", opt_qrates_params[i]);
    fprintf(stdout, "\n");
  }
}

static void process_basefreqs()
{
  long i;
  char dna[4] = "TCAG";

  if (opt_model == BPP_DNA_MODEL_GTR)
  {
    if (opt_basefreqs_fixed)
    {
      fprintf(stdout, "Base frequencies: fixed\n");
    }
    else
    {
      fprintf(stdout, "Base frequencies: generated from Dirichlet\n");
    }
    fprintf(stdout, "                 ");
    if (opt_basefreqs_fixed)
    {
      for (i = 0; i < 4; ++i)
        fprintf(stdout, " P(%c)=%f", dna[i], opt_basefreqs_params[i]);
    }
    else
    {
      for (i = 0; i < 4; ++i)
        fprintf(stdout, " a_%c=%f", dna[i], opt_basefreqs_params[i]);
    }
    fprintf(stdout, "\n");
  }
}

static double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
   /* returns the incomplete gamma ratio I(x,alpha) where x is the upper
              limit of the integration and alpha is the shape parameter.
      returns (-1) if in error
      ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
      (1) series expansion,     if (alpha>x || x<=1)
      (2) continued fraction,   otherwise
      RATNEST FORTRAN by
      Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
      19: 285-287 (AS32)
   */
   int i;
   double p = alpha, g = ln_gamma_alpha;
   double accurate = 1e-10, overflow = 1e60;
   double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6];

   if (x == 0) return (0);
   if (x < 0 || p <= 0) return (-1);

   factor = exp(p*log(x) - x - g);
   if (x > 1 && x >= p) goto l30;
   /* (1) series expansion */
   gin = 1;  term = 1;  rn = p;
l20:
   rn++;
   term *= x / rn;   gin += term;
   if (term > accurate) goto l20;
   gin *= factor / p;
   goto l50;
l30:
   /* (2) continued fraction */
   a = 1 - p;   b = a + x + 1;  term = 0;
   pn[0] = 1;  pn[1] = x;  pn[2] = x + 1;  pn[3] = x*b;
   gin = pn[2] / pn[3];
l32:
   a++;
   b += 2;
   term++;
   an = a*term;
   for (i = 0; i < 2; i++)
      pn[i + 4] = b*pn[i + 2] - an*pn[i];
   if (pn[5] == 0) goto l35;
   rn = pn[4] / pn[5];
   dif = fabs(gin - rn);
   if (dif > accurate) goto l34;
   if (dif <= accurate*rn) goto l42;
l34:
   gin = rn;
l35:
   for (i = 0; i < 4; i++) pn[i] = pn[i + 2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i = 0; i < 4; i++) pn[i] /= overflow;
   goto l32;
l42:
   gin = 1 - factor*gin;

l50:
   return (gin);
}

static double QuantileNormal(double prob)
{
   /* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
      returns (-9999) if in error
      Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
      Applied Statistics 22: 96-97 (AS70)

      Newer methods:
        Wichura MJ (1988) Algorithm AS 241: the percentage points of the
          normal distribution.  37: 477-484.
        Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
          points of the normal distribution.  26: 118-121.
   */
   double a0 = -.322232431088, a1 = -1, a2 = -.342242088547, a3 = -.0204231210245;
   double a4 = -.453642210148e-4, b0 = .0993484626060, b1 = .588581570495;
   double b2 = .531103462366, b3 = .103537752850, b4 = .0038560700634;
   double y, z = 0, p = prob, p1;

   p1 = (p < 0.5 ? p : 1 - p);
   if (p1 < 1e-20) z = 999;
   else {
      y = sqrt(log(1 / (p1*p1)));
      z = y + ((((y*a4 + a3)*y + a2)*y + a1)*y + a0) / ((((y*b4 + b3)*y + b2)*y + b1)*y + b0);
   }
   return (p < 0.5 ? -z : z);
}

static double QuantileChi2(double prob, double v)
{
   /* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
      returns -1 if in error.   0.000002<prob<0.999998
      RATNEST FORTRAN by
          Best DJ & Roberts DE (1975) The percentage points of the
          Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
      Converted into C by Ziheng Yang, Oct. 1993.
   */
   double e = .5e-6, aa = .6931471805, p = prob, g, smallv = 1e-6;
   double xx, c, ch, a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0, s1, s2, s3, s4, s5, s6;

   if (p < smallv)   return(0);
   if (p > 1 - smallv) return(9999);
   if (v <= 0)      return (-1);

   g = lgamma(v / 2);
   xx = v / 2;   c = xx - 1;
   if (v >= -1.24*log(p)) goto l1;

   ch = pow((p*xx*exp(g + xx*aa)), 1 / xx);
   if (ch - e < 0) return (ch);
   goto l4;
l1:
   if (v > .32) goto l3;
   ch = 0.4;   a = log(1 - p);
l2:
   q = ch;  p1 = 1 + ch*(4.67 + ch);  p2 = ch*(6.73 + ch*(6.66 + ch));
   t = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch*(13.32 + 3 * ch)) / p2;
   ch -= (1 - exp(a + g + .5*ch + c*aa)*p2 / p1) / t;
   if (fabs(q / ch - 1) - .01 <= 0) goto l4;
   else                       goto l2;

l3:
   x = QuantileNormal(p);
   p1 = 0.222222 / v;
   ch = v*pow((x*sqrt(p1) + 1 - p1), 3.0);
   if (ch > 2.2*v + 6)
      ch = -2 * (log(1 - p) - c*log(.5*ch) + g);
l4:
   q = ch;   p1 = .5*ch;
   if ((t = IncompleteGamma(p1, xx, g)) < 0)
      fatal("\nIncompleteGamma");
   p2 = p - t;
   t = p2*exp(xx*aa + g + p1 - c*log(ch));
   b = t / ch;  a = 0.5*t - b*c;

   s1 = (210 + a*(140 + a*(105 + a*(84 + a*(70 + 60 * a))))) / 420;
   s2 = (420 + a*(735 + a*(966 + a*(1141 + 1278 * a)))) / 2520;
   s3 = (210 + a*(462 + a*(707 + 932 * a))) / 2520;
   s4 = (252 + a*(672 + 1182 * a) + c*(294 + a*(889 + 1740 * a))) / 5040;
   s5 = (84 + 264 * a + c*(175 + 606 * a)) / 2520;
   s6 = (120 + c*(346 + 127 * c)) / 5040;
   ch += t*(1 + 0.5*t*s1 - b*c*(s1 - b*(s2 - b*(s3 - b*(s4 - b*(s5 - b*s6))))));
   if (fabs(q / ch - 1) > e) goto l4;

   return (ch);
}

#define QuantileGamma(prob,alpha,beta) QuantileChi2(prob,2.0*(alpha))/(2.0*(beta))

static int DiscreteGamma(double freqK[], double rK[], double alpha, double beta, int K, int UseMedian)
{
   /* discretization of G(alpha, beta) with equal proportions in each category.
   */
   int i;
   double t, mean = alpha / beta, lnga1;

   if (UseMedian) {   /* median */
      for (i = 0; i < K; i++) rK[i] = QuantileGamma((i*2. + 1) / (2.*K), alpha, beta);
      for (i = 0, t = 0; i < K; i++) t += rK[i];
      for (i = 0; i < K; i++) rK[i] *= mean*K / t;   /* rescale so that the mean is alpha/beta. */
   }
   else {            /* mean */
      lnga1 = lgamma(alpha + 1);
      for (i = 0; i < K - 1; i++) /* cutting points, Eq. 9 */
         freqK[i] = QuantileGamma((i + 1.0) / K, alpha, beta);
      for (i = 0; i < K - 1; i++) /* Eq. 10 */
         freqK[i] = IncompleteGamma(freqK[i] * beta, alpha + 1, lnga1);
      rK[0] = freqK[0] * mean*K;
      for (i = 1; i < K - 1; i++)  rK[i] = (freqK[i] - freqK[i - 1])*mean*K;
      rK[K - 1] = (1 - freqK[K - 2])*mean*K;
   }

   for (i = 0; i < K; i++) freqK[i] = 1.0 / K;

   return (0);
}

static double * rates4sites(double locus_siterate_alpha, int cdf)
{
  long i,j,k;
  double * rates = NULL;

  double * rK = (double *)xmalloc((size_t)opt_siterate_cats * sizeof(double));
  double * freqK = (double *)xmalloc((size_t)opt_siterate_cats * sizeof(double));
  double * Falias = (double *)xmalloc((size_t)opt_siterate_cats * sizeof(double));
  int * Lalias = (int *)xmalloc((size_t)opt_siterate_cats * sizeof(int));
  int * counts = (int *)xmalloc((size_t)opt_siterate_cats * sizeof(int));

  assert(locus_siterate_alpha);

  rates = (double *)xmalloc((size_t)opt_locus_simlen * sizeof(double));

  if (opt_siterate_cats > 1)
  {
    double gamma_a = locus_siterate_alpha;
    double gamma_b = gamma_a;

    DiscreteGamma(freqK,rK,gamma_a,gamma_b,opt_siterate_cats,BPP_FALSE);
    MultiNomialAliasSetTable(opt_siterate_cats, freqK, Falias, Lalias);
    MultiNomialAlias(opt_locus_simlen,opt_siterate_cats,Falias,Lalias,counts);

    for (i = 0, k = 0; i < opt_siterate_cats; ++i)
      for (j = 0; j < counts[i]; ++j)
        rates[k++] = rK[i];
  }
  else
  {
    for (i = 0; i < opt_locus_simlen; ++i)
      rates[i] = legacy_rndgamma(thread_index_zero,locus_siterate_alpha) /
                 locus_siterate_alpha;
  }
  if (cdf)
  {
    for (i = 1; i < opt_locus_simlen; ++i)
      rates[i] += rates[i-1];
    for (i = 0; i < opt_locus_simlen; ++i)
      rates[i] /= rates[opt_locus_simlen-1];
  }

  free(rK);
  free(freqK);
  free(Falias);
  free(Lalias);
  free(counts);

  return rates;
}

static void evolve_jc69_recursive(gnode_t * node,
                                  double locus_siterate_alpha,
                                  double * site_rates)
{
  long i,k;
  double r;
  char dna[4] = "TCAG";
  int inverse[9] = { -1, 0, 1, -1, 2, -1, -1, -1, 3};

  /* See notes in make_root_seq(). inverse[9] is used to convert from unary code
     (T=1,C=2,A=4,G=8) back to 0,1,2,3 code for states */

  /* sanity check */
  if (locus_siterate_alpha && fabs(site_rates[opt_locus_simlen-1] -1) > 1e-4)
    fatal("[ERROR] rates CDF: 1 = %.6f?\n", site_rates[opt_locus_simlen-1]);

  assert(node->parent);

  char * xparent = (char *)(node->parent->data);  /* parent sequence */
  char * x = (char *)(node->data);                /* current sequence */

  /* copy parent sequence to current node */
  memcpy(x,xparent,opt_locus_simlen * sizeof(char));
    
  /* generate number of mutations */
  long mut_count = legacy_rndpoisson(thread_index_zero,
                                     node->length * opt_locus_simlen);

  for (i = 0; i < mut_count; ++i)
  {
    /* get a position for the mutation */
    if (locus_siterate_alpha == 0)
      k = (int)(legacy_rndu(thread_index_zero) * opt_locus_simlen);
    else
      for (k = 0, r = legacy_rndu(thread_index_zero); k < opt_locus_simlen; ++k)
        if (r < site_rates[k])
          break;

    /* generate new state */
    int state = (int)(legacy_rndu(thread_index_zero) * 3);
    if (state >= inverse[(int)x[k]])
      state++;

    x[k] = pll_map_nt_tcag[(int)dna[state]];
    assert(x[k] > 0 && x[k] < 16);
  }

  /* recursively process subtree */
  if (node->left)
    evolve_jc69_recursive(node->left,  locus_siterate_alpha, site_rates);
  if (node->right)
    evolve_jc69_recursive(node->right, locus_siterate_alpha, site_rates);
}

static void evolve_gtr_recursive(gnode_t * node,
                                 double locus_siterate_alpha,
                                 double * site_rates,
                                 double * eigenvecs,
                                 double * inv_eigenvecs,
                                 double * eigenvals)
{
  long i,j,k;
  long states = 4;
  double * srates = NULL;
  double constrate[1] = {1};
  char dna[4] = "TCAG";
  int inverse[9] = { -1, 0, 1, -1, 2, -1, -1, -1, 3};

  /* See notes in make_root_seq(). inverse[9] is used to convert from unary code
     (T=1,C=2,A=4,G=8) back to 0,1,2,3 code for states */

  double * pmatrix = (double *)xcalloc(16,sizeof(double));

  if (!site_rates)
    srates = constrate;
  else
    srates = site_rates;
  
  char * xparent = (char *)(node->parent->data);  /* parent sequence */
  char * x = (char *)(node->data);                /* current sequence */

  /* copy parent sequence to current node */
  memcpy(x,xparent,opt_locus_simlen * sizeof(char));
    

  assert(node->parent);
  for (i = 0; i < opt_locus_simlen; ++i)
  {
    unsigned int matrix_indices[1] = {0};
    unsigned int params_indices[1] = {0};
    if (i == 0 || (site_rates && site_rates[i] != site_rates[i-1]))
    {
      /* calculate pmatrix using eigendecomposition */
      pll_core_update_pmatrix(&pmatrix,
                              4,
                              1,
                              srates + ((site_rates) ? i : 0),
                              &(node->length),
                              matrix_indices,
                              params_indices,
                              &eigenvals,
                              &eigenvecs,
                              &inv_eigenvecs,
                              1,
                              opt_arch);
      /* make rows of pmatrix cumulative */
      for (j = 0; j < states; ++j)
        for (k = 1; k < states; ++k)
          pmatrix[j*states+k] += pmatrix[j*states+k-1];
    }
    
    double r = legacy_rndu(thread_index_zero);
    for (j = 0; j < states-1; j++)
      if (r < pmatrix[inverse[(int)x[i]]*states+j])
        break;

    x[i] = pll_map_nt_tcag[(int)dna[j]];
  }

  free(pmatrix);

  /* recursively process subtree */
  if (node->left)
    evolve_gtr_recursive(node->left,  locus_siterate_alpha, site_rates, eigenvecs, inv_eigenvecs, eigenvals);
  if (node->right)
    evolve_gtr_recursive(node->right, locus_siterate_alpha, site_rates, eigenvecs, inv_eigenvecs, eigenvals);
}

static void make_root_seq(gnode_t * root, double * freqs)
{
  long i,j;
  double r;
  double p[4];  /* cumulative prob */

  char * x = (char *)(root->data);
  char dna[4] = "TCAG";

  /* Note: pll_map_nt_tcag uses unary code for states to make consensus
     sequences from diploid data easier to construct (using OR). Therefore
     I use the dna[4] array to convert from 0,1,2,3 -> characters and then
     characters as indices to pll_map_nt_tcag such that t = 1, c=2, g=4, t=8 */

  if (opt_model == BPP_DNA_MODEL_JC69)
  {
    for (i = 0; i < opt_locus_simlen; ++i)
      x[i] = pll_map_nt_tcag[(int)dna[(int)(legacy_rndu(thread_index_zero)*4)]];
  }
  else
  {
    for (i = 0; i < 4; ++i) p[i]  = freqs[i]; 
    for (i = 1; i < 4; ++i) p[i] += p[i-1]; 

    for (i = 0; i < opt_locus_simlen; ++i)
    {
      for (j = 0, r = legacy_rndu(thread_index_zero); j < 4-1; ++j)
        if (r < p[j]) break;
      x[i] = pll_map_nt_tcag[(int)dna[j]];
    }
  }
}

static list_t * create_maplist(stree_t * stree)
{
  long i;

  list_t * list = (list_t *)xcalloc(1,sizeof(list_t));

  for (i = 0; i < stree->tip_count; ++i)
  {
    mapping_t * m = (mapping_t *)xmalloc(sizeof(mapping_t));

    m->individual = xstrdup(stree->nodes[i]->label);
    m->species    = xstrdup(stree->nodes[i]->label);
    m->lineno     = i+1;

    list_append(list,(void *)m);
  }

  return list;
}

/* correlated clock */
static void compute_relaxed_rates_recursive(snode_t * node, gtree_t * gtree)
{
  if (!node) return;

  assert(node->parent);

  if (node->parent->tau == 0)
    node->rate = node->parent->rate;
  else
  {
    if (opt_rate_prior == BPP_BRATE_PRIOR_GAMMA)
    {
      /* gamma prior */

      double a = gtree->rate_mui * gtree->rate_mui / gtree->rate_nui;
      node->rate = legacy_rndgamma(thread_index_zero,a) /
                   a * node->parent->rate;
    }
    else
    {
      /* log-normal prior */

      double nv = node->parent->rate +
                  sqrt(gtree->rate_nui)*rndNormal(thread_index_zero);
      node->rate = exp(nv);
    }
  }

  compute_relaxed_rates_recursive(node->left,gtree);
  compute_relaxed_rates_recursive(node->right,gtree);
    
}


static void relaxed_clock_branch_lengths(stree_t * stree, gtree_t * gtree)
{
  /* TODO: Implement networks */
  long i;
  assert(stree->hybrid_count == 0);
  snode_t * start;
  snode_t * end;

  assert(opt_clock == BPP_CLOCK_IND || opt_clock == BPP_CLOCK_CORR);

  long total_nodes = stree->tip_count + stree->inner_count + stree->hybrid_count;

  if (opt_clock == BPP_CLOCK_IND)
  {
    if (opt_rate_prior == BPP_BRATE_PRIOR_LOGNORMAL)
    {
      /* log-normal distributed branch rates */
      for (i = 0; i < total_nodes; ++i)
      {
        double nv = gtree->rate_mui +
                    sqrt(gtree->rate_nui)*rndNormal(thread_index_zero);
        stree->nodes[i]->rate = exp(nv);
      }
    }
    else
    {
      assert(opt_rate_prior == BPP_BRATE_PRIOR_GAMMA);
      /* gamma distributed branch rates */
      double a = gtree->rate_mui * gtree->rate_mui / gtree->rate_nui;
      double b = gtree->rate_mui / gtree->rate_nui;
      for (i = 0; i < total_nodes; ++i)
        stree->nodes[i]->rate = legacy_rndgamma(thread_index_zero,a) / b;
    }
  }
  else
  {
    assert(opt_clock == BPP_CLOCK_CORR);

    /* correlated clock */
    
    stree->root->rate = gtree->rate_mui;
    compute_relaxed_rates_recursive(stree->root->left,gtree);
    compute_relaxed_rates_recursive(stree->root->right,gtree);
  }

  /* now update gene tree */
  for (i = 0; i < gtree->tip_count + gtree->inner_count; ++i)
  {
    gtree->nodes[i]->length = 0;
    if (!gtree->nodes[i]->parent) continue;

    gnode_t * x = gtree->nodes[i];

    double t = x->time;

    start = x->pop;
    end   = x->parent->pop;
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
        assert(x->hpath[hindex] != BPP_HPATH_NONE);
        assert(start->left);
        if (x->hpath[hindex] == BPP_HPATH_RIGHT)
          start = start->hybrid;
      }

      /* skip using branch rates on horizontal edges in hybridization events */
      if (!(pop->hybrid && pop->htau == 0))
        x->length += (start->tau - t)*pop->rate;
      t = start->tau;
    }
    x->length += (x->parent->time - t) * x->parent->pop->rate;
  }
}

static void randomize_order(long * order, long n)
{
  long i,k;

  if (!g_order)
    g_order = (long *)xmalloc((size_t)n * sizeof(long));

  for (i = 0; i < n; ++i)
    g_order[i] = i;

  for (i = 0; i < n; ++i)
  {
    k = (long)((n-i)*legacy_rndu(thread_index_zero));
    order[i] = g_order[i+k];
    g_order[i+k] = g_order[i];
  }
}

static char * consensus(const char * x, const char * y, long n)
{
  long i;

  char * cons = (char *)xmalloc((size_t)n * sizeof(char));

  for (i = 0; i < n; ++i)
    cons[i] = x[i] | y[i];

  return cons;
}

static void collapse_diploid(stree_t * stree, gtree_t * gtree, msa_t * msa, long hets_count)
{
  long i,j,k,m;
  char ** sequence;
  char ** label;

  long seq_sum = 0;
  for (i = 0; i < stree->tip_count; ++i)
    seq_sum += opt_sp_seqcount[i];
  assert(seq_sum == msa->count);

  /* allocate new storage for sequences and labels */
  sequence = (char **)xmalloc((size_t)hets_count * sizeof(char *));
  label    = (char **)xmalloc((size_t)hets_count * sizeof(char *));

  k = 0; j = 0; m = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    if (opt_diploid[i])
    {
      /* create a lowercase version of species label */
      char * seqname = xstrdup(stree->nodes[i]->label);
      for (j = 0; j < (long)strlen(seqname); ++j)
        seqname[j] = xtolower(seqname[j]);

      #ifndef SIM_OLD_FORMAT
      long index = 0;
      #endif
      for (j = m; j < m+opt_sp_seqcount[i]; j += 2)
      {
        #ifdef SIM_OLD_FORMAT
        xasprintf(label+k, "%s%ld_%ld^%s", seqname, j-m+1, j-m+2, stree->nodes[i]->label);
        #else
        xasprintf(label+k, "%s%ld^%s", seqname, ++index, stree->nodes[i]->label);
        #endif
        sequence[k] = consensus(msa->sequence[j],msa->sequence[j+1],msa->length);

        gtree->nodes[j]->data = sequence[k];
        gtree->nodes[j+1]->data = NULL;

        ++k;

        free(msa->sequence[j]);
        free(msa->sequence[j+1]);
        free(msa->label[j]);
        free(msa->label[j+1]);

      }
      free(seqname);
    }
    else
    {
      for (j = m; j < m+opt_sp_seqcount[i]; ++j)
      {
        label[k]    = msa->label[j];
        sequence[k] = msa->sequence[j];
        ++k;
      }
    }
    m += opt_sp_seqcount[i];
  }

  free(msa->sequence);
  free(msa->label);
  msa->sequence = sequence;
  msa->label = label;

  msa->count = hets_count;
}

static void write_concat_seqs(FILE * fp, msa_t ** msa)
{
  long i,j,k;

  fprintf(fp,"\n\n%d %ld \n\n",msa[0]->count,opt_locus_simlen*opt_locus_count);

  char * x = (char *)xmalloc((size_t)(opt_locus_simlen * opt_locus_count) *
                             sizeof(char));

  for (i = 0; i < msa[0]->count; ++i)
  {
    /* concatenate sequences */
    for (k = 0; k < opt_locus_count; ++k)
      memcpy(x+k*opt_locus_simlen,msa[k]->sequence[i],opt_locus_simlen*sizeof(char));

    fprintf(fp, "%s%-*s ", "", 10, msa[0]->label[i]);
    for (j = 0; j < opt_locus_simlen*opt_locus_count; ++j)
    {
      if (j % 10 == 0) fprintf(fp, " ");
      fprintf(fp,"%c", charmap_nt_tcag[(int)x[j]]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n\n");

  free(x);
}

static void write_seqs(FILE * fp, msa_t * msa, long species_count)
{
  long i,j;

  fprintf(fp, "\n\n%d %ld \n\n", msa->count, opt_locus_simlen);

  long seq_sum = 0;
  for (i = 0; i < species_count; ++i)
    seq_sum += opt_sp_seqcount[i] / (opt_diploid[i] ? 2 : 1);
  //assert(seq_sum == msa->count);

  for (i = 0; i < msa->count; ++i)
  {

    fprintf(fp, "%s%-*s ", "", 10, msa->label[i]);
    for (j = 0; j < opt_locus_simlen; ++j)
    {
      if (j % 10 == 0) fprintf(fp, " ");
      fprintf(fp,"%c", charmap_nt_tcag[(int)msa->sequence[i][j]]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n\n");
}

static void write_diploid_rand_seqs(FILE * fp_seqrand,
                                    stree_t * stree,
                                    msa_t * msa)
{
  long i,j,k,m;
  char ** sequence;

  msa_t * new_msa = (msa_t *)xcalloc(1,sizeof(msa_t));
  new_msa->count = msa->count;
  new_msa->length = msa->length;
  new_msa->label = msa->label;

  long seq_sum = 0;
  for (i = 0; i < stree->tip_count; ++i)
    seq_sum += opt_sp_seqcount[i];
  assert(seq_sum == msa->count);

  /* allocate new storage for sequences and labels */
  sequence = (char **)xmalloc((size_t)(msa->count) * sizeof(char *));
  for (i = 0; i < msa->count; ++i)
    sequence[i] = (char *)xmalloc((size_t)(msa->length) * sizeof(char));

  j = 0; m = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    if (opt_diploid[i])
    {
      for (j = m; j < m+opt_sp_seqcount[i]; j += 2)
      {
        memcpy(sequence[j],   msa->sequence[j],   msa->length * sizeof(char));
        memcpy(sequence[j+1], msa->sequence[j+1], msa->length * sizeof(char));
        for (k = 0; k < msa->length; ++k)

        /* randomly resolve */
        if (sequence[j][k] != sequence[j+1][k] && legacy_rndu(thread_index_zero)<0.5)
          SWAP(sequence[j][k],sequence[j+1][k]);
      }
    }
    else
    {
      for (j = m; j < m+opt_sp_seqcount[i]; ++j)
        memcpy(sequence[j], msa->sequence[j], msa->length * sizeof(char));
    }
    m += opt_sp_seqcount[i];
  }

  new_msa->sequence = sequence;

  write_seqs(fp_seqrand, new_msa, stree->tip_count);
  
  for (i = 0; i < msa->count; ++i)
    free(sequence[i]);
  free(sequence);
  free(new_msa);
}

static void set_migration_rates(stree_t * stree)
{
  long i,j;

  assert(opt_migration_labels);

  assert(!opt_msci);

  assert(opt_migration == stree->tip_count + stree->inner_count);

  long nodes_count = stree->tip_count + stree->inner_count;

  /* check that the order of migration matrix labels corresponds to the order
     of labels in tree->nodes */
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    if (strcmp(stree->nodes[i]->label,opt_migration_labels[i]))
      break;

  if (i != stree->tip_count+stree->inner_count)
  {
    fprintf(stderr, "Please use the following order for migration matrix cells:\n");
    for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
      fprintf(stderr, "  %s\n", stree->nodes[i]->label);
    fatal("Order of migration matrix entries mismatch.");
  }

  /* reset invalid entries in the migration matrix */
  long reset_count = 0;
  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * x = stree->nodes[i];
    for (j = i; j < stree->tip_count + stree->inner_count; ++j)
    {
      snode_t * y = stree->nodes[j];

      if (x == y)
      {
        if (opt_migration_matrix[i*nodes_count + j] > 0) reset_count++;
        opt_migration_matrix[i*nodes_count+j] = 0;
      }
      else if (stree->pptable[x->node_index][y->node_index] ||
               stree->pptable[y->node_index][x->node_index] ||
               x->tau > y->parent->tau ||
               y->tau > x->parent->tau)
      {
        if (opt_migration_matrix[i*nodes_count + j] > 0) reset_count++;
        if (opt_migration_matrix[j*nodes_count + i] > 0) reset_count++;
        if (opt_migration_matrix[i*nodes_count + j] > 0 ||
            opt_migration_matrix[j*nodes_count+i] > 0)
          fprintf(stderr,
                  "\nMigration between %s and %s is impossible...\n",
                  x->label, y->label);
        opt_migration_matrix[i*nodes_count+j] = opt_migration_matrix[j*nodes_count+i] = -1;
      }
    }
  }

  if (reset_count)
    fprintf(stderr,"\nResetting %ld migration rates in the migration matrix\n",
            reset_count);

  /* print migration matrix on screen */
  fprintf(stdout, "\nMigration matrix:\n");
  for (i = 0; i < nodes_count*nodes_count; ++i)
  {
    printf(" %f", opt_migration_matrix[i]);
    if (i && ((i+1) % nodes_count) == 0)
      printf("\n");
  }

  for (i = 0; i < nodes_count; ++i)
    for (j = 0; j < nodes_count; ++j)
      if (opt_migration_matrix[i*nodes_count+j] < 0)
        opt_migration_matrix[i*nodes_count+j] = 0;

  long needtheta;
  long die = 0;
  for (i = 0; i < stree->tip_count; ++i)
  {
    for (j = 0, needtheta=0; j < nodes_count; ++j)
      if (opt_migration_matrix[j*nodes_count+i]) needtheta = 1;
    if (needtheta && stree->nodes[i]->theta <= 0)
    {
      fprintf(stderr, " [Error]: theta for %s should be > 0\n", stree->nodes[i]->label);
      die = 1;
    }
  }
  if (die)
    fatal("Please fix the migration matrix in the control file");
}


static void simulate(stree_t * stree)
{
  long i,j,k,m;
  long hets;
  double tmrca = 0;
  long * siteorder = NULL;
  double * eigenvecs = NULL;
  double * inv_eigenvecs = NULL;
  double * eigenvals = NULL;
  FILE * fp_seq = NULL;
  FILE * fp_concat = NULL;
  FILE * fp_tree = NULL;
  FILE * fp_param = NULL;
  FILE * fp_map = NULL;
  FILE * fp_seqfull = NULL;
  FILE * fp_seqrand = NULL;

  /* open output files */
  if (opt_msafile)
    fp_seq = xopen(opt_msafile, "w");
  if (opt_concatfile)
    fp_concat = xopen(opt_concatfile, "w");
  if (opt_treefile)
    fp_tree = xopen(opt_treefile, "w");
  if (opt_modelparafile)
    fp_param = xopen(opt_modelparafile, "w");
  assert(opt_mapfile);
  if (opt_mapfile)
    fp_map = xopen(opt_mapfile, "w");

  /* print list of output files */
  if (opt_msafile)
    fprintf(stdout, "Sequence data file -> %s\n", opt_msafile);
  if (opt_treefile)
    fprintf(stdout, "Trees -> %s\n", opt_treefile);
  if (opt_concatfile)
    fprintf(stdout, "Concatenated sequence alignment -> %s\n", opt_msafile);
  if (opt_modelparafile)
    fprintf(stdout, "Model parameters for loci -> %s\n", opt_modelparafile);
  if (opt_mapfile)
    fprintf(stdout, "Tags to species mapping (Imap) -> %s\n", opt_mapfile);

  if (opt_migration)
    set_migration_rates(stree);

  hets = 0;
  for (i = 0; i < stree->tip_count; ++i)
    hets += opt_sp_seqcount[i] / (opt_diploid[i] ? 2 : 1);

  /* print model parameter file header */
  if (opt_modelparafile)
  {
    if (opt_model == BPP_DNA_MODEL_JC69)
    {
      if (opt_est_locusrate)
        fprintf(fp_param, "locus\tmu_i");
    }
    else if (opt_model == BPP_DNA_MODEL_GTR)
    {
      if (opt_est_locusrate)
        fprintf(fp_param, "locus\tQrates_abcdef\tpi_TACG\talpha\tmu_i");
      else
        fprintf(fp_param, "locus\tQrates_abcdef\tpi_TACG\talpha");
    }
    else
      assert(0);

    if (opt_clock != BPP_CLOCK_GLOBAL)
      fprintf(fp_param, "\tnu_i");
    
    fprintf(fp_param, "\n");
  }

  /* print imap file */
  for (i = 0; i < stree->tip_count; ++i)
    fprintf(fp_map, "%s\t%s\n", stree->nodes[i]->label, stree->nodes[i]->label);


  /* allocate MSA structures */
  msa_t ** msa = (msa_t **)xmalloc((size_t)opt_locus_count * sizeof(msa_t *));
  for (i = 0; i < opt_locus_count; ++i)
    msa[i] = (msa_t *)xcalloc(1,sizeof(msa_t));

  /* allocate placeholder for gene trees */
  gtree_t ** gtree = (gtree_t **)xmalloc((size_t)opt_locus_count *
                                         sizeof(gtree_t *));

  long locus_seqcount;
  double qrates[6];
  double freqs[4];
  double locus_siterate_alpha = opt_siterate_alpha;
  double * siterates = NULL;
  double * mui_array = NULL;
  double * vi_array = NULL;

  /* allocate array used for shuffling order of sites */
  if (opt_msafile)
    siteorder = (long *)xmalloc((size_t)opt_locus_simlen * sizeof(long));

//  if (opt_siterate_fixed)
//    locus_siterate_alpha = opt_siterate_alpha;

  /* store number of sequences per locus (before collpasing diploid seqs) */
  locus_seqcount = 0;
  for (i = 0; i < stree->tip_count; ++i)
    locus_seqcount += opt_sp_seqcount[i];

  if (opt_msafile)
  {
    if (hets < locus_seqcount)
    {
      char * filename = NULL;
      xasprintf(&filename, "%s.full", opt_msafile);
      fp_seqfull = xopen(filename,"w");
      free(filename);

      filename = NULL;
      xasprintf(&filename, "%s.rand", opt_msafile);
      fp_seqrand = xopen(filename,"w");
      free(filename);
    }
  }


  /* 1. create maplist (Imap) and 2. initialize two hashtables which are used
     when calling gtree_simulate for quick access to a sequence population */
  list_t * maplist = create_maplist(stree);
  gtree_simulate_init(stree,maplist);
  list_clear(maplist,map_dealloc);
  free(maplist);

  /* allocate eigendecomposition structures */
  if (opt_model == BPP_DNA_MODEL_GTR)
  {
    eigenvecs = (double *)xmalloc(16*sizeof(double));
    inv_eigenvecs = (double *)xmalloc(16*sizeof(double));
    eigenvals = (double *)xmalloc(4*sizeof(double));
  }

  /* pre-generate mu_i and v_i */
  if (opt_est_locusrate)
  {
    mui_array = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    {
      /* generate locus rates from Gamma(a_mui,a_mui/mubar) */
      for (i = 0; i < opt_locus_count; ++i)
        mui_array[i] = legacy_rndgamma(thread_index_zero,opt_mui_alpha) /
                       (opt_mui_alpha/opt_locusrate_mubar);
    }
    else
    {
      /* generate locus rate from Dir(a_mui) */
      assert(opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR);
      double tmp_sum = 0;
      for (i = 0; i < opt_locus_count; ++i)
      {
        mui_array[i] = legacy_rndgamma(thread_index_zero,opt_mui_alpha);
        tmp_sum += mui_array[i];
      }
      for (i = 0; i < opt_locus_count; ++i)
      {
        mui_array[i] /= tmp_sum;
        mui_array[i] *= opt_locusrate_mubar*opt_locus_count;
      }
    }
  }

  if (opt_clock != BPP_CLOCK_GLOBAL)
  {
    vi_array = (double *)xmalloc((size_t)opt_locus_count * sizeof(double));
    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    {
      /* generate v_i from Gamma(a_vi,a_vi/vbar) */
      for (i = 0; i < opt_locus_count; ++i)
        vi_array[i] = legacy_rndgamma(thread_index_zero,opt_vi_alpha) /
                      (opt_vi_alpha/opt_clock_vbar);
    }
    else
    {
      /* generate v_i from Dir(a_vi) */
      assert(opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR);
      double tmp_sum = 0;
      for (i = 0; i < opt_locus_count; ++i)
      {
        vi_array[i] = legacy_rndgamma(thread_index_zero,opt_vi_alpha);
        tmp_sum += vi_array[i];
      }
      for (i = 0; i < opt_locus_count; ++i)
      {
        vi_array[i] /= tmp_sum;
        vi_array[i] *= opt_clock_vbar*opt_locus_count;
      }
    }
  }

  /* start generating */
  for (i = 0; i < opt_locus_count; ++i)
  {
    if (opt_modelparafile)
      fprintf(fp_param, "%ld", i+1);


    if (opt_model == BPP_DNA_MODEL_GTR)
    {
      if (!opt_qrates_fixed)
      {
        legacy_rnddirichlet(thread_index_zero,qrates,opt_qrates_params,6);
        for (j = 0; j < 6; ++j)
          qrates[j] /= qrates[5];
      }
      else
        memcpy(qrates,opt_qrates_params,6*sizeof(double));

      if (!opt_basefreqs_fixed)
        legacy_rnddirichlet(thread_index_zero,freqs,opt_basefreqs_params,4);
      else
        memcpy(freqs,opt_basefreqs_params,4*sizeof(double));

      /* print parameters in parameter file */
      assert(opt_modelparafile);
      for (j = 0; j < 6; ++j)
        fprintf(fp_param," %9.6f", qrates[j]);
      for (j = 0; j < 4; ++j)
        fprintf(fp_param," %8.6f", freqs[j]);

      pll_update_eigen(eigenvecs,inv_eigenvecs,eigenvals,freqs,qrates,4,4);
    }

    if (!opt_siterate_fixed)
    {
      locus_siterate_alpha = legacy_rndgamma(thread_index_zero,opt_siterate_alpha) /
                             opt_siterate_beta;
      if (opt_modelparafile)
        fprintf(fp_param, " %9.6f", locus_siterate_alpha);
    }
    if (opt_est_locusrate)
      fprintf(fp_param, " %9.6f", mui_array[i]);
    if (opt_clock != BPP_CLOCK_GLOBAL)
      fprintf(fp_param, " %9.6f", vi_array[i]);

    if (opt_msafile || opt_treefile)
    {
      msa[i]->label    = (char**)xmalloc((size_t)locus_seqcount*sizeof(char *));
      msa[i]->sequence = (char**)xmalloc((size_t)locus_seqcount*sizeof(char *));

      /* create sequence labels and populate msa structure */
      for (j = 0, m = 0; j < stree->tip_count; ++j)
      {
        if (opt_diploid[j])
          for (k = 0; k < opt_sp_seqcount[j]; ++k)
            xasprintf(msa[i]->label+m++,
                      "%s%ld%c^%s",
                      stree->nodes[j]->label, 
                      k / 2 + 1,
                      (char)('a' + k % 2),
                      stree->nodes[j]->label);
        else
          for (k = 0; k < opt_sp_seqcount[j]; ++k)
            xasprintf(msa[i]->label+m++,
                      "%s%ld^%s",
                      stree->nodes[j]->label, 
                      k + 1,
                      stree->nodes[j]->label);

        msa[i]->count += opt_sp_seqcount[j];
      }
      msa[i]->length = opt_locus_simlen;

      /* change all sequence labels to lowercase */
      for (j = 0; j < m; ++j)
        for (k = 0; k < (long)strlen(msa[i]->label[j]) && msa[i]->label[j][k] != '^'; ++k)
          msa[i]->label[j][k] = xtolower(msa[i]->label[j][k]);
    }

    /* simulate gene tree */
    gtree[i] = gtree_simulate(stree,msa[i],i);

    if (opt_est_locusrate)
      gtree[i]->rate_mui = mui_array[i];
    else
      gtree[i]->rate_mui = 1;

    if (opt_clock != BPP_CLOCK_GLOBAL)
      gtree[i]->rate_nui = vi_array[i];

    tmrca += gtree[i]->root->time;

    /* set branch lengths */
    for (j = 0; j < gtree[i]->tip_count+gtree[i]->inner_count; ++j)
    {
      if (!gtree[i]->nodes[j]->parent) continue;

      gtree[i]->nodes[j]->length = gtree[i]->nodes[j]->parent->time -
                                   gtree[i]->nodes[j]->time;
    }
      
    assert(locus_seqcount == gtree[i]->tip_count);

    /* TODO: Count 3S trees */

    /* if clock is assumed, compute species tree branch rates and write them to
       file */
    if (opt_clock == BPP_CLOCK_IND || opt_clock == BPP_CLOCK_CORR)
    {
      relaxed_clock_branch_lengths(stree, gtree[i]);
      assert(opt_modelparafile);
      for (j = 0; j < stree->tip_count + stree->inner_count; ++j)
      {
        fprintf(fp_param, " %.6f", stree->nodes[j]->rate);
      }
    }
    if (opt_modelparafile)
      fprintf(fp_param, "\n");

    /* multiply branches with locus rate */
    if (opt_est_locusrate && opt_clock == BPP_CLOCK_GLOBAL)
    {
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        gtree[i]->nodes[j]->length *= mui_array[i];
    }


    if (opt_treefile)
    {
      char * newick = gtree_export_newick(gtree[i]->root, NULL);
      fprintf(fp_tree, "%s [TH=%.6f]\n", newick, gtree[i]->root->time);
      free(newick);
    }

    if (opt_msafile)
    {
      /* calculate rates for each site */
      if (locus_siterate_alpha)
        siterates = rates4sites(locus_siterate_alpha,
                                (opt_model == BPP_DNA_MODEL_JC69));

      /* allocate space for sequences and map each to a gene tree node */
      char ** x = (char**)xmalloc((size_t)(gtree[i]->tip_count+gtree[i]->inner_count)*
                                  sizeof(char *));
      for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
      {
        x[j] = (char *)xmalloc((size_t)opt_locus_simlen * sizeof(char));
        gtree[i]->nodes[j]->data = (void *)(x[j]);
      }

      /* map also gene tree tip node sequences to the msa alignment structure,
         and free the placeholder */
      for (j = 0; j < gtree[i]->tip_count; ++j)
        msa[i]->sequence[j] = x[j];
      free(x);

      /* generate a sequence at the root */
      make_root_seq(gtree[i]->root, freqs);

      /* recursively generate ancestral sequences and tip sequences */
      if (opt_model == BPP_DNA_MODEL_JC69)
      {
        evolve_jc69_recursive(gtree[i]->root->left, locus_siterate_alpha, siterates);
        evolve_jc69_recursive(gtree[i]->root->right, locus_siterate_alpha, siterates);
      }
      else
      {
        evolve_gtr_recursive(gtree[i]->root->left,locus_siterate_alpha,siterates,
                             eigenvecs,inv_eigenvecs,eigenvals);
        evolve_gtr_recursive(gtree[i]->root->right,locus_siterate_alpha,siterates,
                             eigenvecs,inv_eigenvecs,eigenvals);
      }

      /* shuffle order of sites */
      if (locus_siterate_alpha && opt_siterate_cats > 1)
      {
        randomize_order(siteorder, opt_locus_simlen);
        char * tmpseq = (char *)xmalloc((size_t)opt_locus_simlen * sizeof(char));
        for (j = 0; j < gtree[i]->tip_count + gtree[i]->inner_count; ++j)
        {
          char * seq = (char *)(gtree[i]->nodes[j]->data);   
          memcpy(tmpseq, seq, opt_locus_simlen * sizeof(char));
          for (k = 0; k < opt_locus_simlen; ++k)
            seq[k] = tmpseq[siteorder[k]];
        }
        free(tmpseq);
      }

      /* collapse diploid sequences */
      if (hets < msa[i]->count)
      {
        /* write full data first */
        if (fp_seqfull)
          write_seqs(fp_seqfull, msa[i], stree->tip_count);
        if (fp_seqrand)
          write_diploid_rand_seqs(fp_seqrand,stree,msa[i]);

        /* then collapse sequences */
        collapse_diploid(stree,gtree[i],msa[i],hets);
      }

      /* write sequences */
      write_seqs(fp_seq, msa[i], stree->tip_count);

      /* TODO: Instead of freeing and allocating, create siterates once and
         fill it with ones, and use rates4sites to alter it */
      if (siterates)
        free(siterates);
    }

    if ((i+1) % 1000 == 0 || (opt_locus_count > 1000 && i == opt_locus_count-1))
      printf("%10ld replicates done... mean tMRCA = %9.6f\n", i+1, tmrca/(i+1));

  }  /* end of locus loop */
  
  if (opt_concatfile)
  {
    fprintf(stdout, "Generating concatenated sequence alignment...\n");
    write_concat_seqs(fp_concat, msa);

  }

  if (opt_migration)
  {
    long matrix_size = opt_migration * opt_migration;

    for (i = 0; i < matrix_size; ++i)
      opt_migration_events[i] /= opt_locus_count;

    printf("\nCounts of migration events averaged over replicates: %8.4f\n",
           tmrca / opt_locus_count);
    /* print migration matrix on screen */
    for (i = 0; i < matrix_size; ++i)
    {
      printf(" %f", opt_migration_events[i]);
      if (i && ((i+1) % opt_migration) == 0)
        printf("\n");
    }
  }

  if (mui_array)
    free(mui_array);
  if (vi_array)
    free(vi_array);
  
  /* deallocate hashtables used for mapping sequences to species */
  gtree_simulate_fini();

  /* deallocate alignment structures and gene trees */
  for (i = 0; i < opt_locus_count; ++i)
  {
    gtree_destroy(gtree[i],free);

    for (j = 0; j < msa[i]->count; ++j)
      free(msa[i]->label[j]);
    free(msa[i]->label);

    /* remember sequences were already freed as they are mapped to the gene tree */
    free(msa[i]->sequence);
    free(msa[i]);
  }
  free(msa);
  free(gtree);



  /* deallocate eigendecomposition structures */
  if (opt_model == BPP_DNA_MODEL_GTR)
  {
    free(eigenvecs);
    free(inv_eigenvecs);
    free(eigenvals);
  }

  /* close all open output files */
  if (opt_msafile)
    fclose(fp_seq);
  if (opt_concatfile)
    fclose(fp_concat);
  if (opt_treefile)
    fclose(fp_tree);
  if (opt_modelparafile)
    fclose(fp_param);
  assert(opt_mapfile);
  if (opt_mapfile)
    fclose(fp_map);

  if (fp_seqfull)
    fclose(fp_seqfull);
  if (fp_seqrand)
    fclose(fp_seqrand);

  if (g_order)
    free(g_order);
  if (siteorder)
    free(siteorder);
}

static void assign_thetas(stree_t * stree)
{
  long i;
  assert(opt_msci);
  /* this is copied from stree_init_theta(...) */

  for (i = 0; i < stree->tip_count + stree->inner_count; ++i)
    stree->nodes[i]->has_theta = 1;

  for (i = 0; i < stree->tip_count; ++i)
  {
    if (opt_sp_seqcount[i] < 2)
    {
      stree->nodes[i]->theta = -1;
      stree->nodes[i]->has_theta = 0;
    }
  }

  for (i = stree->tip_count; i < stree->tip_count + stree->inner_count; ++i)
  {
    snode_t * node = stree->nodes[i];

    if (opt_msci && node->hybrid)
    {
      if (!node_is_bidirection(node))
      {
        /* node is hybridization: we assign a theta to the nodes that
           compose it that have a 'tau-parent' (htau) annotation */
        if (node->htau)
        {
          node->has_theta = 1;
          assert(node->theta > 0);
        }
        else
        {
          node->theta = -1;
          node->has_theta = 0;
        }

        if (node->hybrid->htau)
        {
          node->hybrid->has_theta = 1;
          assert(node->hybrid->theta > 0);
        }
        else
        {
          node->hybrid->theta = -1;
          node->hybrid->has_theta = 0;
        }
      }
      else
      {
        /* bidirectional introgression */

        node->has_theta = 1;
        node->hybrid->theta = -1;
        node->hybrid->has_theta = 0;
      }
    }
  }
}

static void validate_and_set_taus(stree_t * stree)
{
  long i;
  long hoffset = stree->tip_count + stree->inner_count;
  double tau = 0;

  for (i = 0; i < stree->hybrid_count; ++i)
  {
    snode_t * mnode = stree->nodes[hoffset+i];
    snode_t * hnode = mnode->hybrid;

    if (!node_is_bidirection(hnode))
    {
      /* hybridization */
      if (mnode->tau > 0 && hnode->tau > 0 && mnode->tau != hnode->tau)
        fatal("Misspecified tau for hybridization node %s", hnode->label);

      if (hnode->tau > 0)
        tau = hnode->tau;
      if (mnode->tau > 0)
        tau = mnode->tau;

      assert(hnode->parent);
      assert(mnode->parent);

      if (tau && !hnode->htau && hnode->parent->tau > 0 &&
          hnode->parent->tau != tau)
        fatal("Conflicting tau for nodes %s and %s",
              hnode->parent->label, hnode->label);

      if (tau && !mnode->htau && mnode->parent->tau > 0 &&
          mnode->parent->tau != tau)
        fatal("Conflicting tau for nodes %s and %s",
              mnode->parent->label, hnode->label);

      if (mnode->tau <= 0 && hnode->tau <= 0)
      {
        if (hnode->htau && mnode->htau)
          fatal("Missing tau for hybridization node %s", hnode->label);
        else if (!hnode->htau && !mnode->htau)
        {
          if (hnode->parent->tau <= 0 && mnode->parent->tau <= 0)
            fatal("Missing tau for hybridization node %s", hnode->label);

          if (hnode->parent->tau > 0 && mnode->parent->tau > 0 &&
              hnode->parent->tau != mnode->parent->tau)
            fatal("Conflicting tau for nodes %s and %s",
                  hnode->parent->label, mnode->parent->label);

          if (hnode->parent->tau <= 0)
            tau = mnode->parent->tau;
          else
            tau = hnode->parent->tau;
        }
        else if (!hnode->htau)
        {
          if (hnode->parent->tau <= 0)
            fatal("Missing tau for hybridization node %s", hnode->label);

          tau = hnode->parent->tau;
        }
        else
        {
          if (mnode->parent->tau <= 0)
            fatal("Missing tau for hybridization node %s", hnode->label);

          tau = mnode->parent->tau;
        }
      }

      assert(tau > 0);

      mnode->tau = hnode->tau = tau;
      if (!mnode->htau)
        mnode->parent->tau = tau;
      if (!hnode->htau)
        hnode->parent->tau = tau;
    }
    else
    {
      /* bidirectional introgression */

      snode_t * ohnode = mnode->parent;
      snode_t * omnode = hnode->right;

      assert(ohnode->hybrid == omnode && omnode != NULL);

      if (mnode->tau > 0 && hnode->tau > 0 && mnode->tau != hnode->tau)
        fatal("Misspecified tau for introgression node %s", hnode->label);

      double atau = 0;
      if (hnode->tau > 0)
        atau = hnode->tau;
      if (mnode->tau > 0)
        atau = mnode->tau;

      if (omnode->tau > 0 && ohnode->tau > 0 && omnode->tau != ohnode->tau)
        fatal("Misspecified tau for introgression node %s", ohnode->label);

      double btau = 0;
      if (ohnode->tau > 0)
        btau = ohnode->tau;
      if (omnode->tau > 0)
        btau = omnode->tau;

      if (atau <= 0 && btau <= 0)
        fatal("Missing tau for introgression node %s", hnode->label);
      else if (atau > 0 && btau > 0 && atau != btau)
        fatal("Conflicting tau for introgression nodes %s and %s",
              hnode->label, ohnode->label);

      if (atau > 0)
        tau = atau;
      else
        tau = btau;

      hnode->tau = mnode->tau = ohnode->tau = omnode->tau = tau;
    }
  }
}

void cmd_simulate()
{
  long i;
  stree_t * stree;

  assert(opt_streenewick);

  stree = bpp_parse_newick_string(opt_streenewick);

  assert(opt_msci == !!stree->hybrid_count);

  for (i = 0; i < stree->tip_count + stree->inner_count + stree->hybrid_count; ++i)
    stree->nodes[i]->tau = stree->nodes[i]->length;

  /* allocate and set pptable */
  stree_init_pptable(stree);
  stree_label(stree);

  if (opt_msci)
  {
    long hoffset = stree->tip_count+stree->inner_count;
    for (i = 0; i < stree->hybrid_count; ++i)
      stree->nodes[hoffset+i]->tau = stree->nodes[hoffset+i]->hybrid->tau;

    assign_thetas(stree);
    validate_and_set_taus(stree);
  }


  /* allocate space for keeping track of coalescent events at each species tree
     node for each locus */
  stree->locus_count = (unsigned int)opt_locus_count;
  for (i = 0; i < stree->tip_count+stree->inner_count+stree->hybrid_count; ++i)
  {
    snode_t * snode = stree->nodes[i];

    snode->event = (dlist_t **)xcalloc(opt_locus_count, sizeof(dlist_t *));
    snode->event_count = (int *)xcalloc(opt_locus_count, sizeof(int));
    snode->seqin_count = (int *)xcalloc(opt_locus_count, sizeof(int));
    snode->gene_leaves = (unsigned int *)xcalloc(opt_locus_count,sizeof(unsigned int));
    /* TODO: The next two allocations might not be necessary when computing
       theta analytically */
    snode->logpr_contrib = (double*)xcalloc(opt_locus_count, sizeof(double));
    snode->old_logpr_contrib = (double *)xcalloc(opt_locus_count, sizeof(double));

    snode->t2h = NULL;
    snode->old_t2h = NULL;
    if (!opt_est_theta)
    {
      snode->t2h = (double*)xcalloc((size_t)opt_locus_count, sizeof(double));
      snode->old_t2h = (double*)xcalloc((size_t)opt_locus_count, sizeof(double));
      snode->t2h_sum = 0;
      snode->event_count_sum = 0;
    }

    long j;
    for (j = 0; j < stree->locus_count; ++j)
      snode->event[j] = dlist_create();
  }

  process_subst_model();
  process_basefreqs();
  print_settings(stree);
  process_diploid(stree->tip_count);

  simulate(stree);

  stree_destroy(stree,free);

  if (opt_diploid)
    free(opt_diploid);
}
