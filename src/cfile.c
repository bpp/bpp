/*
    Copyright (C) 2016-2018 Tomas Flouri and Ziheng Yang

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

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static long is_emptyline(const char * line)
{
  size_t ws = strspn(line, " \t\r\n");
  if (!line[ws] || line[ws] == '*' || line[ws] == '#') return 1;
  return 0;
}


static long get_string(const char * line, char ** value)
{
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star and hash */
  while (*p && *p != '*' && *p != '#') p++;


  /* now go back, skipping all white-space until a character occurs */
  for (--p; *p == ' ' || *p == '\t' || *p == '\r' || *p =='\n'; --p);
  
  char * end = p+1;
  *end = 0;

  *value = xstrdup(start);
  free(s);

  return ws + end - start;
}

static long get_long(const char * line, long * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%ld%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static long get_double(const char * line, double * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%lf%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}

static long get_e(const char * line, long * value)
{
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    *value = 0;
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  if ((*start != 'E' && *start != 'e') || (start+1 != end))
  {
    free(s);
    *value = 2;         /* erroneous value */
    return 0;
  }

  *value = 1;

  free(s);
  return ws + end - start;
}

static int parse_speciestree(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_est_stree);
  if (!count) goto l_unwind;

  p += count;

  if (opt_est_stree == 0 || opt_est_stree == 1) ret = 1;

  /* TODO: At the momemt we ignore the pslider/expandratio/shrinkratio */

l_unwind:
  free(s);
  return ret;
}

static long parse_diploid(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  size_t ws;
  size_t count = 0;

  assert(line);

  while (1)
  {
    /* check if 0 or 1 */
    if (*p == '0' || *p == '1')
    {
      count++;
      p++;
    }

    if ((*p == '*') || (*p == '#') || (*p == '\0')) break;

    /* make sure we have only white-space and skip it all */
    ws = strspn(p, " \t\r\n");
    if (!ws) goto l_unwind;

    p += ws;
  }

  if (!count) goto l_unwind;

  /* allocate memory */
  opt_diploid_size = count;
  opt_diploid = (long *)xmalloc(count*sizeof(long));

  /* go through the string once more */
  p = s;
  count = 0;
  while (1)
  {
    /* check if 0 or 1 */
    if (*p == '0' || *p == '1')
    {
      opt_diploid[count++] = *p - '0';
      p++;
    }

    if ((*p == '*') || (*p == '#') || (*p == '\0')) break;

    /* skip all white-space */
    ws = strspn(p, " \t\r\n");
    if (!ws) goto l_unwind;

    p += ws;
  }
  ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_checkpoint(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_checkpoint_initial);
  if (!count) goto l_unwind;

  p += count;

  ret = 1;

  count = get_long(p, &opt_checkpoint_step);
  if (!count) goto l_unwind;

  p += count;

  ret = 0;

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_speciesdelimitation(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_est_delimit);
  if (!count) goto l_unwind;

  p += count;

  /* if species tree is fixed (opt_est_delimit = 0) finish parsing */
  if (opt_est_delimit == 0) goto l_finish;

  /* if first value other than 0 or 1, error */
  if (opt_est_delimit != 1) goto l_unwind;

  /* now read second token */
  count = get_long(p, &opt_rjmcmc_method);
  if (!count) goto l_unwind;

  p += count;

  /* species delimitation algorithm can be 0 or 1 */
  if (opt_rjmcmc_method != 0 && opt_rjmcmc_method != 1) return 0;

  /* now read third token */
  count = get_double(p, opt_rjmcmc_method ? 
                     &opt_rjmcmc_alpha : &opt_rjmcmc_epsilon);
  if (!count) goto l_unwind;

  p += count;

  /* if algorithm 0 then finish */
  if (opt_rjmcmc_method == 0) goto l_finish;

  /* now read fourth token */
  count = get_double(p, &opt_rjmcmc_mean);
  if (!count) goto l_unwind;

  p += count;

l_finish:
  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_locusrate(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_est_locusrate);
  if (!count) goto l_unwind;

  printf("est_locusrate = %ld\n", opt_est_locusrate);

  p += count;

  if (is_emptyline(p) && !opt_est_locusrate) ret = 1;

  count = 0;
  if (opt_est_locusrate == MUTRATE_ESTIMATE)
  {
    count = get_double(p, &opt_locusrate_alpha);
    if (!count) goto l_unwind;
    printf("alpha = %f\n", opt_locusrate_alpha);
  }
  else if (opt_est_locusrate == MUTRATE_FROMFILE)
  {
    count = get_string(p,&opt_locusrate_filename);
    if (!count) goto l_unwind;
  }
  else
    goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;
  printf("ret= %ld\n", ret);

l_unwind:
  free(s);
  return ret;
}

static long parse_thetaprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_double(p, &opt_theta_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_theta_beta);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;

  count = get_e(p, &opt_est_theta);
  if (opt_est_theta > 1) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;
  
l_unwind:
  free(s);
  return ret;
}

static long parse_tauprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  /* TODO: Add third options for dirichlet */

  long count;

  count = get_double(p, &opt_tau_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_tau_beta);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;
  
l_unwind:
  free(s);
  return ret;
}

static long parse_heredity(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_est_heredity);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p) && !opt_est_heredity) ret = 1;

  count = 0;
  if (opt_est_heredity == HEREDITY_ESTIMATE)
  {
    count = get_double(p, &opt_heredity_alpha);
    if (!count) goto l_unwind;

    p += count;

    count = get_double(p, &opt_heredity_beta);
    if (!count) goto l_unwind;
  }
  else if (opt_est_heredity == HEREDITY_FROMFILE)
  {
    count = get_string(p,&opt_heredity_filename);
    if (!count) goto l_unwind;
  }
  else
    goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_finetune(const char * line)
{
  long ret = 0;
  long count;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (p[ws] != '0' && p[ws] != '1') goto l_unwind;

  if (p[ws] == '1') opt_finetune_reset = 1;

  p += ws+1;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* skip colon */
  if (p[ws] != ':') goto l_unwind;

  p += ws+1;

  /* now read 7 values */

  /* 1. gene tree age finetune param */
  count = get_double(p, &opt_finetune_gtage);
  if (!count) goto l_unwind;

  p += count;
  
  /* 2. gene tree spr finetune param */
  count = get_double(p, &opt_finetune_gtspr);
  if (!count) goto l_unwind;

  p += count;

  /* 3. theta finetune param */
  count = get_double(p, &opt_finetune_theta);
  if (!count) goto l_unwind;

  p += count;

  /* 4. theta finetune param */
  count = get_double(p, &opt_finetune_tau);
  if (!count) goto l_unwind;

  p += count;

  /* 5. theta finetune param */
  count = get_double(p, &opt_finetune_mix);
  if (!count) goto l_unwind;

  p += count;

  /* TODO: The next is not implemented yet */
  double opt_finetune_seqerr;

  /* 6. locusrate finetune param */
  count = get_double(p, &opt_finetune_locusrate);
  if (!count) goto l_unwind;

  p += count;

  /* 7. sequence error finetune param */
  count = get_double(p, &opt_finetune_seqerr);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
  
}

static long parse_print(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  /* samples */
  count = get_long(p, &opt_print_samples);
  if (!count) goto l_unwind;

  p += count;

  /* locusrate */
  count = get_long(p, &opt_print_locusrate);
  if (!count) goto l_unwind;

  p += count;

  /* heredity scalars */
  count = get_long(p, &opt_print_hscalars);
  if (!count) goto l_unwind;

  p += count;

  /* gene trees */
  count = get_long(p, &opt_print_genetrees);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;
  
l_unwind:
  free(s);
  return ret;
}

static long readandvalidatecount(const char * line, long spcount)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  
  long i = 0;
  long count;

  opt_sp_seqcount = (long*)xmalloc((size_t)spcount*sizeof(long));

  while (spcount)
  {
    count = get_long(p, opt_sp_seqcount+i);
    if (!count) break;

    p += count;

    --spcount;
    ++i;
  }

  /* line contains less entries than number of species */
  if (spcount) goto l_unwind;

  /* otherwise check if nothing else is there */
  if (is_emptyline(p)) ret = 1;

  l_unwind:
    free(s);

  return ret;
}

static long parse_speciesandtree(const char * line, long * spcount)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * seqnames; /* white-space separated sequence labels */
  long alloc_size = 1;

  long count;

  long seq_count;
  count = get_long(p, &seq_count);
  *spcount = seq_count;
  if (!count) goto l_unwind;

  p += count;

  /* now read the remaining part of the line and trim comments */
  if (!get_string(p,&seqnames)) goto l_unwind;

  /* now we need to convert it to a comma-separated list for compatibility
     with the opt_reorder functionality */

  /* TODO: Currently spaces are not supported in species names */
  p = seqnames;

  /* assertion that seqnames contains at least one character and that the last
     character is not a whitespace */
  assert(strlen(seqnames) > 0 && 
         strcspn(seqnames+strlen(seqnames)-1," \t\n\r") > 0);

  opt_reorder = (char *)xcalloc(alloc_size,sizeof(char));
  while (*p)
  {
    /* read sequence label */
    count = strcspn(p," \t\n\r");
    if (count == 0 && seq_count != 0) goto l_unwind;

    /* decrease number of remaining sequences */
    --seq_count;

    char * label = xstrndup(p,count);

    long comma = seq_count ? 1 : 0;

    /* append to opt_reorder */
    char * temp = (char *)xmalloc((alloc_size+count+comma)*sizeof(char));
    memcpy(temp,opt_reorder,alloc_size*sizeof(char));
    strcat(temp,label);
    alloc_size += count+comma;
    free(opt_reorder);
    free(label);
    opt_reorder = temp;

    if (comma)
    {
      opt_reorder[alloc_size-2] = ',';
      opt_reorder[alloc_size-1] = '\0';
    }

    p += count;

    /* skip white-space */
    count = strspn(p, " \t\n\r");
    p += count;
  }

  free(seqnames);
  if (!seq_count) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_long(const char * line, long * value)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, value);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long get_token(char * line, char ** token, char ** value)
{
  char * p = line;

  /* here we parse lines which may be in the form:
     
     token = value
  */

  /* skip all white-space */
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

  /* is it a blank line or comment ? */
  if (!*p || *p == '*' || *p == '#') return 0;

  /* store address of token's beginning */
  *token = p;

  /* find occurrence of '=' */
  while (*p && *p != '=') ++p;

  /* if no '=' found return error */
  if (!*p) return -1;

  /* '=' was found, store pointer to value */
  *value = p+1;

  /* '=' was found, move back to the last letter of token (ignore whitespace) */
  for (--p; *p == ' ' || *p == '\t' || *p == '\r' || *p =='\n'; --p);

  /* return length of token */
  return p - *token + 1;
}

static void check_validity()
{
  if (!opt_streenewick)
    fatal("Initial species tree newick format is required in 'species&tree'");

  if (!opt_outfile)
    fatal("Option 'outfile' is required");

  if (!opt_mcmcfile)
    fatal("Option 'mcmcfile' is required");

  if (opt_method < 0 || opt_method > 3)
    fatal("Invalid method");

  if (!opt_usedata && opt_bfbeta != 1)
    fatal("Cannot use option option 'BayesFactorBeta' when usedata=0");

  if (opt_samples < 1)
    fatal("Option 'nsample' must be a positive integer greater than zero");
  
  if (opt_samplefreq < 1)
    fatal("Option 'sampfreq' must be a positive integer greater than zero");

  if (opt_burnin < 0)
    fatal("Option 'burnin' must be a positive integer or zero");

  /* species delimitation specific checks */
  if (opt_method == METHOD_10)          /* species delimitation */
  {
    if (opt_rjmcmc_method == 0)
    {
      if (opt_rjmcmc_epsilon <= 0)
        fatal("rj-MCMC epsilon must be a positive real greater than zero");
    }
    else if (opt_rjmcmc_method == 1)
    {
      if (opt_rjmcmc_alpha <= 0)
        fatal("rj-MCMC alpha must be a positive real greater than zero");

      if (opt_rjmcmc_mean <= 0)
        fatal("rj-MCMC mean must be a positive real greater than zero");
    }
    else
      fatal("Internal error in deciding rjMCMC algorithm");

    if ((opt_delimit_prior == BPP_SPECIES_PRIOR_SLH) ||
        (opt_delimit_prior == BPP_SPECIES_PRIOR_SUNIFORM))
      fatal("Invalid 'speciesmodelprior' value");
  }
}

void load_cfile()
{
  FILE * fp;
  long line_count = 0;

  fp = xopen(opt_cfile,"r");

  while (getnextline(fp))
  {
    int valid = 0;
    char * token;
    char * value;
    long token_len;

    ++line_count;
    token_len = get_token(line,&token,&value);

    if (!token_len) continue;
    if (token_len < 0)
      fatal("Invalid syntax when parsing file %s on line %ld",
            opt_cfile, line_count);
    
    if (token_len == 4)
    {
      if (!strncasecmp(token,"seed",4))
      {
        if (!parse_long(value,&opt_seed))
          fatal("Option 'seed' expects one integer (line %ld)", line_count);

        if (opt_seed == -1)
          opt_seed = (long)time(NULL);

        valid = 1;
      }
      else if (!strncasecmp(token,"arch",4))
      {
        char * temp;
        if (!get_string(value,&temp))
          fatal("Option %s expects a string (line %ld)", token, line_count);

        if (!strcmp(temp,"cpu"))
          opt_arch = PLL_ATTRIB_ARCH_CPU;
        else if (!strcasecmp(temp,"sse"))
          opt_arch = PLL_ATTRIB_ARCH_SSE;
        else if (!strcasecmp(temp,"avx"))
          opt_arch = PLL_ATTRIB_ARCH_AVX;
        else if (!strcasecmp(temp,"avx2"))
          opt_arch = PLL_ATTRIB_ARCH_AVX2;
        else
          fatal("Invalid instruction set (%s) (line %ld)", temp, line_count);

        free(temp);

        valid = 1;
      }
    }
    else if (token_len == 5)
    {
      if (!strncasecmp(token,"nloci",5))
      {
        if (!parse_long(value,&opt_locus_count) || opt_locus_count < 0)
          fatal("Option 'nloci' expects a positive integer or zero (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"print",5))
      {
        if (!parse_print(value))
          fatal("Option 'print' expects four bits (line %ld)", line_count);

        if (opt_print_samples == 0)
          fatal("First bit of 'print' must be set to 1");
        valid = 1;
      }
    }
    else if (token_len == 6)
    {
      if (!strncasecmp(token,"burnin",6))
      {
        if (!parse_long(value,&opt_burnin) || opt_burnin < 0)
          fatal("Option 'burnin' expects positive (or zero) integer (line %ld)",
                 line_count);
        valid = 1;
      }
    }
    else if (token_len == 7)
    {
      if (!strncasecmp(token,"diploid",7))
      {
        if (!parse_diploid(value))
          fatal("Option %s expects values 0 or 1 for each species (line %ld)",
                token,line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"seqfile",7))
      {
        if (!get_string(value, &opt_msafile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"outfile",7))
      {
        if (!get_string(value, &opt_outfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"usedata",7))
      {
        if (!parse_long(value,&opt_usedata) ||
            (opt_usedata != 0 && opt_usedata != 1))
          fatal("Option 'usedata' expects value 0 or 1 (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"nsample",7))
      {
        if (!parse_long(value,&opt_samples) || opt_samples <= 0)
          fatal("Option 'nsample' expects a positive integer (line %ld)",
                 line_count);
        valid = 1;
      }
    }
    else if (token_len == 8)
    {
      if (!strncasecmp(token,"imapfile",8))
      {
        if (!get_string(value, &opt_mapfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"mcmcfile",8))
      {
        if (!get_string(value,&opt_mcmcfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"tauprior",8))
      {
        if (!parse_tauprior(value))
          fatal("Option 'tauprior' expects two doubles (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"heredity",8))
      {
        if (!parse_heredity(value))
          fatal("Invalid format of 'heredity' (line %ld) ", line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"finetune",8))
      {
        if (!parse_finetune(value))
          fatal("Option 'finetune' in wrong format (line %ld)", line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"sampfreq",8))
      {
        if (!parse_long(value,&opt_samplefreq) || opt_samplefreq <= 0)
          fatal("Option 'samplfreq' expects a positive integer (line %ld)",
                line_count);
        valid = 1;
      }
    }
    else if (token_len == 9)
    {
      if (!strncasecmp(token,"cleandata",9))
      {
        if (!parse_long(value,&opt_cleandata) ||
            (opt_cleandata != 0 && opt_cleandata != 1))
          fatal("Option 'cleandata' expects value 0 or 1 (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"locusrate",9))
      {
        if (!parse_locusrate(value))
          fatal("Erroneous format of 'locusrate' (line %ld)", line_count);
        valid = 1;
      }
    }
    else if (token_len == 10)
    {
      if (!strncasecmp(token,"thetaprior",10))
      {
        if (!parse_thetaprior(value))
          fatal("Option 'thetaprior' expects two doubles (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"checkpoint",10))
      {
        if (!parse_checkpoint(value))
          fatal("Erroneous format of 'checkpoint' (line %ld)", line_count);
        opt_checkpoint = 1;
        if (sizeof(BYTE) != 1)
          fatal("Checkpoint does not work on systems with sizeof(char) != 1");
        valid = 1;
      }
    }
    else if (token_len == 11)
    {
      if (!strncasecmp(token,"speciestree",11))
      {
        if (!parse_speciestree(value))
          fatal("Erroneous format of options speciestree (line %ld)",
                line_count);
        valid = 1;
      }
    }
    else if (token_len == 12)
    {
      if (!strncasecmp(token,"species&tree",12))
      {
        /* TODO: Currently only the old BPP format is allowed. Make it also
           accept only the tree in newick format, i.e. one line instead of 3 */

        long spcount = 0;
        if (!parse_speciesandtree(value,&spcount))
          fatal("Erroneous format of 'species&tree' (line %ld)", line_count);

        if (!getnextline(fp))
          fatal("Incomplete 'species&tree' record (line %ld)", line_count);

        ++line_count;
        if (!readandvalidatecount(line,spcount))
          fatal("Erroneous enumeration of species sequences in 'species&tree' "
                "tag (line %ld).\nExpected number of species is %ld.\n",
                line_count, spcount);

        if (!getnextline(fp))
          fatal("Incomplete 'species&tree' record (line %ld)", line_count);

        ++line_count;

        if (!get_string (line,&opt_streenewick))
          fatal("Expected newick tree string in 'species&tree' (line %ld)",
                 line_count);
        valid = 1;
      }
    }
    else if (token_len == 13)
    {
      if (!strncasecmp(token,"sequenceerror",13))
      {
        fatal("Not implemented (%s)", token);
        valid = 1;
      }
    }
    else if (token_len == 15)
    {
      if (!strncasecmp(token,"bayesfactorbeta",15))
      {
        if (!get_double(value,&opt_bfbeta) || opt_bfbeta <= 0)
          fatal("Option 'bayesfactorbeta' expects a positive real (line %ld)",
                line_count);
        valid = 1;
      }
    }
    else if (token_len == 17)
    {
      if (!strncasecmp(token,"speciesmodelprior",17))
      {
        /* TODO: Currently we allow only priors 0 and 1 */
        if (!parse_long(value,&opt_delimit_prior) ||
            (opt_delimit_prior < BPP_SPECIES_PRIOR_MIN ||
             opt_delimit_prior > BPP_SPECIES_PRIOR_MAX))
          fatal("Option 'speciesmodelprior' expects integer between %d and %d "
                "(line %ld)",
                BPP_SPECIES_PRIOR_MIN, BPP_SPECIES_PRIOR_MAX, line_count);

        valid = 1;
      }
    }
    else if (token_len == 19)
    {
      if (!strncasecmp(token,"speciesdelimitation",19))
      {
        if (!parse_speciesdelimitation(value))
          fatal("Erroneous format of option %s (line %ld)", token, line_count);
        valid = 1;
      }
    }

    if (!valid)
      fatal("Invalid syntax when parsing file %s on line %ld",
            opt_cfile, line_count);
  }

  /* set method */
  if (!opt_est_stree && !opt_est_delimit)
    opt_method = METHOD_00;
  else if (!opt_est_stree)
    opt_method = METHOD_10;
  else if (!opt_est_delimit)
    opt_method = METHOD_01;
  else
    opt_method = METHOD_11;

  check_validity();


  fclose(fp);
}

int parsefile_doubles(const char * filename,
                      long n,
                      double * outbuffer,
                      long * errcontext)
{
  long line_count = 0;
  long count;
  FILE * fp = xopen(filename,"r");

  long entry = 0;

  while (getnextline(fp))
  {
    ++line_count;

    char * p = line;

    while (!is_emptyline(p))
    {
      if (entry == n)
      {
        fclose(fp);
        return ERROR_PARSE_MORETHANEXPECTED;
      }
      count = get_double(p, outbuffer+entry++);
      if (!count)
      {
        fclose(fp);
        *errcontext = line_count;
        return ERROR_PARSE_INCORRECTFORMAT;
      }

      p += count;
      
    }
  }
  if (entry != n)
  {
    fclose(fp);
    *errcontext = entry;
    return ERROR_PARSE_LESSTHANEXPECTED;
  }
  
  fclose(fp);
  return 0;
}
