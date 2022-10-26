/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static char errmsg[512];

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static long species_count = 0;

static const long dna_model_count =  8;
static const long aa_model_count  = 19;

/* Important: CUSTOM *MUST* be last in the list */
static const char * dna_model_name[] = 
 {
   "JC69", "K80", "F81", "HKY", "T92", "TN93", "F84", "GTR"
 };

static const char * aa_model_name[] =
 {
   "DAYHOFF", "LG",    "DCMUT", "JTT",      "MTREV",    "WAG",
   "RTREV",   "CPREV", "VT",    "BLOSUM62", "MTMAM",    "MTART",
   "MTZOA",   "PMB",   "HIVB",  "HIVW",     "JTTDCMUT", "FLU", "STMTREV"
 };

static const long dna_model_index[] =
 {
   BPP_DNA_MODEL_JC69, BPP_DNA_MODEL_K80, BPP_DNA_MODEL_F81, BPP_DNA_MODEL_HKY,
   BPP_DNA_MODEL_T92, BPP_DNA_MODEL_TN93, BPP_DNA_MODEL_F84, BPP_DNA_MODEL_GTR
 };

static const long aa_model_index[] =
 {
   BPP_AA_MODEL_DAYHOFF,  BPP_AA_MODEL_LG,       BPP_AA_MODEL_DCMUT, 
   BPP_AA_MODEL_JTT,      BPP_AA_MODEL_MTREV,    BPP_AA_MODEL_WAG,
   BPP_AA_MODEL_RTREV,    BPP_AA_MODEL_CPREV,    BPP_AA_MODEL_VT,
   BPP_AA_MODEL_BLOSUM62, BPP_AA_MODEL_MTMAM,    BPP_AA_MODEL_MTART,
   BPP_AA_MODEL_MTZOA,    BPP_AA_MODEL_PMB,      BPP_AA_MODEL_HIVB,
   BPP_AA_MODEL_HIVW,     BPP_AA_MODEL_JTTDCMUT, BPP_AA_MODEL_FLU,
   BPP_AA_MODEL_STMTREV
 };

static const char * rate_prior_option[] =
 {
   "DIR", "IID"
 };
static const char * rate_prior_name[] =
 {
   "Gamma-Dirichlet", "Conditional iid"
 };

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

static long starts_with_opar(const char * line)
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

  if (p[ws] != '(')
  {
    free(s);
    return 0;
  }

  return 1;
}

/* get string from current position pointed by line until a delimiter. Create
   new string with result, store it in value, and return number of characters
   read */
static long get_delstring(const char * line, const char * del, char ** value)
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

  /* skip all characters until a delimiter is found */
  char * end = start + strcspn(start, del);

  *end = 0;

  if (start==end)
  {
    free(s);
    return 0;
  }

  *value = xstrdup(start);

  free(s);
  return ws + end - start;
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

static long get_doubleordash(const char * line, double * value)
{
  int len=0;
  long ret = 0;
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

  if (!strcmp(start,"-"))
  {
    free(s);
    return ws+end-start;
  }

  ret = sscanf(start, "%lf%n", value, &len);
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

  /*** Ziheng 2020-9-1 ***/
#if (0)  

  count = get_long(p, &opt_est_stree);
  if (opt_est_stree == 0) goto l_unwind;
  if (opt_est_stree != 1) goto l_unwind;
  p += count;

  count = get_double(p, &opt_prob_snl);
  if (!count) goto l_unwind;
  p += count;

  count = get_double(p, &opt_prob_snl_shrink);
  if (!count) goto l_unwind;
  p += count;

  count = get_double(p, &opt_snl_lambda_expand);
  if (!count) goto l_unwind;
  p += count;

  count = get_double(p, &opt_snl_lambda_shrink);
  if (!count) goto l_unwind;
  p += count;

#else
  count = sscanf(p, "%ld%lf%lf%lf%lf", 
    &opt_est_stree, &opt_prob_snl, &opt_prob_snl_shrink, &opt_snl_lambda_expand, &opt_snl_lambda_shrink);
#endif

l_unwind:
  if (opt_est_stree == 0 || opt_est_stree == 1) ret = 1;
  if (opt_snl_lambda_expand < 0 || opt_snl_lambda_expand>1)
    fatal("opt_snl_lambda_expand should be between 0 and 1.");
  if (opt_snl_lambda_shrink < 0 || opt_snl_lambda_shrink>1)
    fatal("opt_snl_lambda_shrink should be between 0 and 1.");
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
  long ones_count = 0;

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
      if (*p-'0' == 1)
        ones_count++;
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
  if (opt_diploid && !ones_count)
  {
    free(opt_diploid);
    opt_diploid = NULL;
  }
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

static long get_dist(const char * line, long * dist)
{
  int ret=0;
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

  /* invalidate return value */
  *dist = BPP_BRATE_PRIOR_MAX+1;

  ret = ws + end - start;
  if (!strcasecmp(start,"G"))
  {
    *dist = BPP_BRATE_PRIOR_GAMMA;
  }
  else if (!strcasecmp(start,"LN"))
  {
    *dist = BPP_BRATE_PRIOR_LOGNORMAL;
  }
  else
  {
    ret = 0;
  }

  free(s);

  return ret;
}

static long get_priordist(const char * line, long * dist)
{
  int ret=0;
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

  /* invalidate return value */
  *dist = BPP_BRATE_PRIOR_MAX+1;

  ret = ws + end - start;
  if (!strcasecmp(start,"DIR"))
  {
    *dist = BPP_LOCRATE_PRIOR_GAMMADIR;
  }
  else if (!strcasecmp(start,"IID"))
  {
    *dist = BPP_LOCRATE_PRIOR_HIERARCHICAL;
  }
  else
  {
    ret = 0;
  }

  free(s);

  return ret;
}

static long parse_clock(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_clock);
  if (!count) goto l_unwind;

  p += count;

  /* if molecular clock then accept no other parameters */
  if (opt_clock == BPP_CLOCK_GLOBAL)
  {
    if (is_emptyline(p)) ret = 1;
    goto l_unwind;
  }

  /* the only other choice is local clock in which mutation rate drifts over
     branches independently among loci */
  if (opt_clock != BPP_CLOCK_IND && opt_clock != BPP_CLOCK_CORR) goto l_unwind;

  /* read vbar alpha and vbar beta */

  /* get a_vbar */
  count = get_double(p, &opt_vbar_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* get b_vbar */
  count = get_double(p, &opt_vbar_beta);
  if (!count) goto l_unwind;

  p += count;

  /* get a_vi */
  count = get_double(p, &opt_vi_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* next two parameters are optional */
  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* get prior */
  #if 0
  long temp = 0;
  count = get_long(p, &temp);
  if (!count) goto l_unwind;

  if (temp < BPP_LOCRATE_PRIOR_MIN || temp > BPP_LOCRATE_PRIOR_MAX)
    fatal("ERROR: Invalid prior value (%ld) in 'locusrate' option", temp);

  /* check whether prior was already specified in the locusrate option */
  if (opt_locusrate_prior != -1)
  {
    if (temp != opt_locusrate_prior)
      fatal("ERROR: prior = %ld (%s) in 'clock' does not match prior = %ld (%s)"
            " in 'locusrate'", temp, rate_prior_name[temp], opt_locusrate_prior,
            rate_prior_name[opt_locusrate_prior]);
  }
  else
  {
    opt_locusrate_prior = temp;
  }
  #else
  long temp = 0;
  count = get_priordist(p,&temp);
  if (!count) goto l_unwind;

  if (temp < BPP_LOCRATE_PRIOR_MIN || temp > BPP_LOCRATE_PRIOR_MAX)
    fatal("ERROR: Invalid prior value (%ld) in 'locusrate' option", temp);

  /* check whether prior was already specified in the locusrate option */
  if (opt_locusrate_prior != -1)
  {
    if (temp != opt_locusrate_prior)
      fatal("ERROR: prior = %s (%s) in 'clock' does not match prior = %s (%s)"
            " in 'locusrate'", rate_prior_option[temp], rate_prior_name[temp],
            rate_prior_option[opt_locusrate_prior],
            rate_prior_name[opt_locusrate_prior]);
  }
  else
  {
    opt_locusrate_prior = temp;
  }
  #endif

  p += count;

  if (opt_locusrate_prior < BPP_LOCRATE_PRIOR_MIN ||
      opt_locusrate_prior > BPP_LOCRATE_PRIOR_MAX)
    goto l_unwind; 

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* get branch rate prior distribution */
  count = get_dist(p,&opt_rate_prior);
  if (!count) goto l_unwind;

  p += count;

  if (opt_rate_prior < BPP_BRATE_PRIOR_MIN ||
      opt_rate_prior > BPP_BRATE_PRIOR_MAX)
    goto l_unwind;

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

  p += count;

  if (is_emptyline(p) && opt_est_locusrate == MUTRATE_CONSTANT) ret = 1;

  count = 0;
  if (opt_est_locusrate == MUTRATE_ESTIMATE)
  {
    /* get a_mubar */
    count = get_double(p, &opt_mubar_alpha);
    if (!count) goto l_unwind;

    p += count;

    if (is_emptyline(p))
      fatal("The syntax for 'locusrate' tag has changed since BPP v4.1.4.\n"
            "Please refer to the BPP manual for the new syntax.\n");

    /* get b_mubar */
    count = get_double(p, &opt_mubar_beta);
    if (!count) goto l_unwind;

    p += count;

    /* get a_mui */
    count = get_double(p, &opt_mui_alpha);
    if (!count) goto l_unwind;

    p += count;

    /* get optional prior */
    if (is_emptyline(p))
    {
      opt_locusrate_prior = BPP_LOCRATE_PRIOR_GAMMADIR;
      ret = 1;
      goto l_unwind;
    }

    /* get prior */
    #if 0
    long temp = 0;
    count = get_long(p, &temp);
    if (!count) goto l_unwind;

    if (temp < BPP_LOCRATE_PRIOR_MIN || temp > BPP_LOCRATE_PRIOR_MAX)
      fatal("ERROR: Invalid prior value (%ld) in 'clock' option", temp);

    /* check whether prior was already specified in the clock option */
    if (opt_locusrate_prior != -1)
    {
      if (temp != opt_locusrate_prior)
        fatal("ERROR: prior = %ld (%s) in 'locusrate' does not match prior = "
              "%ld (%s) in 'clock'", temp, rate_prior_name[temp],
              opt_locusrate_prior, rate_prior_name[opt_locusrate_prior]);
    }
    else
    {
      opt_locusrate_prior = temp;
    }
    #else
    long temp = 0;
    count = get_priordist(p,&temp);
    if (!count) goto l_unwind;

    if (temp < BPP_LOCRATE_PRIOR_MIN || temp > BPP_LOCRATE_PRIOR_MAX)
      fatal("ERROR: Invalid prior value (%ld) in 'clock' option", temp);

    /* check whether prior was already specified in the clock option */
    if (opt_locusrate_prior != -1)
    {
      if (temp != opt_locusrate_prior)
        fatal("ERROR: prior = %s (%s) in 'locusrate' does not match prior = "
              "%s (%s) in 'clock'",
              rate_prior_option[temp],
              rate_prior_name[temp],
              rate_prior_option[opt_locusrate_prior],
              rate_prior_name[opt_locusrate_prior]);
    }
    else
    {
      opt_locusrate_prior = temp;
    }
    #endif

    p += count;

    if (opt_locusrate_prior < BPP_LOCRATE_PRIOR_MIN ||
        opt_locusrate_prior > BPP_LOCRATE_PRIOR_MAX)
      goto l_unwind;

  }
  else if (opt_est_locusrate == MUTRATE_FROMFILE)
  {
    count = get_string(p,&opt_locusrate_filename);
    if (!count) goto l_unwind;
    
    p += count;
  }
  else
    goto l_unwind;


  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_partition_line(const char * line, partition_t * part)
{
  long i;
  long ret = 0;
  long count;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;
  char * tmp = NULL;

  char * datatype = NULL;
  char * model = NULL;

  count = get_delstring(p," \t\r\n*#,-",&tmp);
  if (!count) goto l_unwind;
  p += count;

  count = get_long(tmp,&part->start);
  if (!count) goto l_unwind;
  part->end = part->start;
  free(tmp);
  tmp=NULL;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");
  p += ws;

  if (is_emptyline(p))
    goto l_unwind;

  if (*p == '-')
  {
    /* parse another number (end) */
    p += 1;

    count = get_delstring(p," \t\r\n*#,-",&tmp);
    if (!count) goto l_unwind;
    p += count;

    count = get_long(tmp,&part->end);
    if (!count) goto l_unwind;
    free(tmp);
    tmp=NULL;

    /* skip all white-space */
    ws = strspn(p, " \t\r\n");
    p += ws;
  }

  if (part->end < part->start)
    goto l_unwind;

  /* skip comma */
  if (*p != ',')
    goto l_unwind;
  p += 1;

  /* skip white-space */
  ws = strspn(p, " \t\r\n");
  p += ws;

  if (is_emptyline(p))
    goto l_unwind;

  count = get_delstring(p," \t\r\n*#,",&datatype);
  if (!count) goto l_unwind;

  p += count;

  if (!strcasecmp(datatype,"DNA"))
    part->dtype = BPP_DATA_DNA;
  else if (!strcasecmp(datatype,"AA"))
    part->dtype = BPP_DATA_AA;
  else
    goto l_unwind;

  /* skip white-space */
  ws = strspn(p, " \t\r\n");
  p += ws;

  if (is_emptyline(p))
    goto l_unwind;

  /* skip comma */
  if (*p != ',')
    goto l_unwind;
  p += 1;
  
  /* read model */
  count = get_delstring(p," \t\r\n*#",&model);
  if (!count) goto l_unwind;

  p += count;

  if (part->dtype == BPP_DATA_DNA)
  {
    /* parse dna model (skip last entry which is CUSTOM) */
    for (i = 0; i < dna_model_count; ++i)
      if (!strcasecmp(model,dna_model_name[i]))
        break;
    if (i == dna_model_count)
      goto l_unwind;
    part->model = dna_model_index[i];
  }
  else if (part->dtype == BPP_DATA_AA)
  {
    /* parse aa model */
    for (i = 0; i < aa_model_count; ++i)
      if (!strcasecmp(model,aa_model_name[i]))
        break;
    if (i == aa_model_count)
      goto l_unwind;
    part->model = aa_model_index[i];
  }
  else
    assert(0);

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  if (datatype) free(datatype);
  if (tmp) free(tmp);
  if (model) free(model);
  return ret;
}

static void validate_partitions(list_t * list)
{
  long i,line_count = 0;
  long start = 0;
  long end = 0;
  assert(list);

  list_item_t * item;

  /* find minimum and maximum partition */
  item = list->head;
  while (item)
  {
    partition_t * part = (partition_t *)(item->data);
    assert(part->model >= 0);
    assert((part->dtype == BPP_DATA_DNA && part->model < dna_model_count) ||
           (part->dtype == BPP_DATA_AA && part->model < aa_model_count));

    assert(part->end >= part->start);

    if (start == 0 || part->start < start)
      start = part->start;
    if (part->end > end)
      end = part->end;

    item = item->next;
  }

  if (start != 1)
    fatal("Partitions in partition file %s must start from locus 1",
          opt_partition_file);
  if (end < start)
    fatal("Invalid partition format in file %s", opt_partition_file);

  /* now check that (a) we have a continuous range of loci in the partitions and
     (b) no two partitions overlap */

  long * tmp = (long *)xcalloc((size_t)(end-start+1),sizeof(long));
  item = list->head;
  while (item)
  {
    partition_t * part = (partition_t *)(item->data);
    assert(part);
    line_count++;

    /* check condition (b) */
    for (i = part->start; i <= part->end; ++i)
    {
      if (tmp[i-1])
        fatal("Partition on line %ld contains locus %ld which is already in "
              "partition on line %ld (file %s)",
              line_count, i, tmp[i-1], opt_partition_file);
      tmp[i-1] = line_count;
    }

    item = item->next;
  }

  /* check condition (a) */
  for (i = 0; i < end; ++i)
    if (!tmp[i])
      fatal("Locus %ld not contained in any partition (file %s)",
            i+1, opt_partition_file);

  /* dealloc */
  free(tmp);
}

static list_t * parse_partition_file()
{
  long line_count = 0;
  long ret = 0;
  FILE * fp;
  list_t * list;

  fp = xopen(opt_partition_file, "r");

  list = (list_t *)xcalloc(1,sizeof(list_t));

  while (getnextline(fp))
  {
    ++line_count;

    partition_t * part = (partition_t *)xmalloc(sizeof(partition_t));
    if (!parse_partition_line(line,part))
    {
      sprintf(errmsg,"Invalid entry on line %ld", line_count);
      goto l_unwind;
    }

    list_append(list,(void *)part);
  }
  ret = 1;

l_unwind:
  fclose(fp);
  if (ret == 0)
  {
    list_clear(list,free);
    free(list);
    list = NULL;
  }

  return list;
}

static long parse_model(const char * line)
{
  long i;
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  char * model = NULL;

  long count;

  count = get_delstring(p," \t\r\n*#",&model);
  if (!count) goto l_unwind;

  p += count;

  /* parse model */

  /* check if DNA model */
  for (i = 0; i < dna_model_count; ++i)
  {
    if (!strcasecmp(model,dna_model_name[i]))
      break;
  }
  if (i < dna_model_count)
  {
    opt_model = dna_model_index[i];
    if (is_emptyline(p)) ret = 1;

    goto l_unwind;
  }
  assert(i == dna_model_count);

  /* check if AA model */
  for (i = 0; i < aa_model_count; ++i)
  {
    if (!strcasecmp(model,aa_model_name[i]))
      break;
  }
  if (i < aa_model_count)
  {
    opt_model = aa_model_index[i];
    if (is_emptyline(p)) ret = 1;

    goto l_unwind;
  }
  assert(i == aa_model_count);

  /* check if custom model */
  if (!strcasecmp(model,"custom"))
  {
    opt_model = BPP_DNA_MODEL_CUSTOM;

    count = get_delstring(p," \t\r\n*#",&opt_partition_file);
    if (!count) goto l_unwind;

    p += count;

    if (is_emptyline(p)) ret = 1;
    goto l_unwind;
  }

l_unwind:
  free(s);
  if (model)
    free(model);
  return ret;
}

static long parse_loadbalance(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  char * lb = NULL;

  long count;

  count = get_delstring(p," \t\r\n*#",&lb);
  if (!count) goto l_unwind;

  p += count;

  if (!is_emptyline(p)) goto l_unwind;

  ret = 1;
  if (!strcasecmp(lb, "zigzag"))
    opt_load_balance = BPP_LB_ZIGZAG;
  else if (!strcasecmp(lb, "none"))
    opt_load_balance = BPP_LB_NONE;
  else
    ret = 0;

l_unwind:
  free(s);
  if (lb)
    free(lb);
  return ret;
}

static long parse_alphaprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_double(p, &opt_alpha_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_alpha_beta);
  if (!count) goto l_unwind;

  p += count;

  /* set default categories */
  opt_alpha_cats = 4;

  if (is_emptyline(p)) ret = 1;

  count = get_long(p, &opt_alpha_cats);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;
  
l_unwind:
  free(s);
  return ret;
}

static long parse_thetamodel(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * model = NULL;

  long count;

  count = get_delstring(p," \t\r\n*#,",&model);
  if (!count) goto l_unwind;

  p += count;

  if (!strcasecmp(model, "linked-none")) 
    opt_linkedtheta = BPP_LINKEDTHETA_NONE;
  else if (!strcasecmp(model, "linked-all"))
    opt_linkedtheta = BPP_LINKEDTHETA_ALL;
  else if (!strcasecmp(model, "linked-inner"))
    opt_linkedtheta = BPP_LINKEDTHETA_INNER;
  else if (!strcasecmp(model, "linked-msci"))
    opt_linkedtheta = BPP_LINKEDTHETA_MSCI;
  else
    goto l_unwind;

  ret = 1;

l_unwind:
  if (model)
    free(model);
  free(s);
  return ret;
}


static long parse_thetaprior_args(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;
  double a,b;
  a = b = 0;

  count = get_double(p, &a);
  if (!count) goto l_unwind;

  /* invgamma using default option (i.e. thetaprior = alpha beta e) */
  p += count;

  /* now read second token */
  count = get_double(p, &b);
  if (!count) goto l_unwind;

  p += count;

  if (opt_theta_dist == BPP_THETA_PRIOR_BETA)
  {
    opt_theta_p = a;
    opt_theta_q = b;
  }
  else
  {
    opt_theta_alpha = a;
    opt_theta_beta = b;
  }

  if (is_emptyline(p) && (opt_theta_dist != BPP_THETA_PRIOR_BETA))
    ret = 1;

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
  {

    count = get_e(p, &opt_est_theta);
    if (opt_est_theta > 1) goto l_unwind;

    p += count;

    if (is_emptyline(p)) ret = 1;
  }
  else if (opt_theta_dist == BPP_THETA_PRIOR_BETA)
  {
    count = get_double(p, &opt_theta_min);
    if (!count) goto l_unwind;

    p += count;

    count = get_double(p, &opt_theta_max);
    if (!count) goto l_unwind;

    p += count;

    if (is_emptyline(p))
      ret = 1;
  }

l_unwind:
  free(s);
  return ret;
}

static long parse_thetaprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * dist = NULL;

  long count;

  opt_theta_dist = BPP_THETA_PRIOR_INVGAMMA;

  /* peek at the first argument. If not a double then read distribution type
     otherwise read the arguments of invgamma (default) */
  count = get_double(p, &opt_theta_alpha);
  if (!count)
  {
    count = get_delstring(p," \t\r\n*#,",&dist);
    if (!count) goto l_unwind;

    p += count;

    if (!strcasecmp(dist,"invgamma"))
      opt_theta_dist = BPP_THETA_PRIOR_INVGAMMA;
    else if (!strcasecmp(dist,"gamma"))
      opt_theta_dist = BPP_THETA_PRIOR_GAMMA;
    else if (!strcasecmp(dist,"beta"))
      opt_theta_dist = BPP_THETA_PRIOR_BETA;
    else
      goto l_unwind;
  }

  ret = parse_thetaprior_args(p);
  
l_unwind:
  if (dist)
    free(dist);
  free(s);
  return ret;
}

static long parse_threads(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  /* read number of threads */
  count = get_long(p, &opt_threads);
  if (!count) goto l_unwind;

  p += count;

  if (opt_threads < 1) goto l_unwind;
  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* read starting thread index */
  count = get_long(p, &opt_threads_start);
  if (!count) goto l_unwind;

  p += count;

  if (opt_threads_start < 1) goto l_unwind;
  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* read thread step */
  count = get_long(p, &opt_threads_step);
  if (!count) goto l_unwind;

  p += count;

  if (opt_threads_step < 1) goto l_unwind;
  if (is_emptyline(p))
    ret = 1;


l_unwind:
  free(s);
  return ret;

}

static long parse_tauprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;
  char * dist = NULL;

  /* TODO: Add third options for dirichlet */

  long count;

  opt_tau_dist = BPP_TAU_PRIOR_INVGAMMA;

  /* peak at the first argument. If not a double then read distribution type
     otherwise read the arguments of invgamma (default) */
  count = get_double(p, &opt_tau_alpha);
  if (!count)
  {
    count = get_delstring(p, " \t\r\n*#,", &dist);
    if (!count) goto l_unwind;

    p += count;
    
    if (!strcasecmp(dist,"invgamma"))
      opt_tau_dist = BPP_TAU_PRIOR_INVGAMMA;
    else if (!strcasecmp(dist, "gamma"))
      opt_tau_dist = BPP_TAU_PRIOR_GAMMA;
    else
      goto l_unwind;
  }

  count = get_double(p, &opt_tau_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_tau_beta);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
  }
  else
  {
    /* now read third token if available */
    double ignorevalue;
    count = get_double(p, &ignorevalue);
    if (!count) goto l_unwind;

    p += count;
    if (is_emptyline(p))
      ret = 1;
  }
  
l_unwind:
  free(s);
  if (dist)
    free(dist);
  return ret;
}

static long parse_phiprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  /* TODO: Add third options for dirichlet */

  long count;

  count = get_double(p, &opt_phi_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_phi_beta);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
  }
  else
  {
    /* now read third token if available */
    double ignorevalue;
    count = get_double(p, &ignorevalue);
    if (!count) goto l_unwind;

    p += count;
    if (is_emptyline(p))
      ret = 1;
  }
  
l_unwind:
  free(s);
  return ret;
}

static long parse_migprior(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_double(p, &opt_mig_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_mig_beta);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  count = get_double(p, &opt_pseudo_alpha);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_double(p, &opt_pseudo_beta);
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

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }


  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* skip colon */
  if (p[ws] != ':') goto l_unwind;

  p += ws+1;

  /* now read 7 values */

  /* 1. gene tree age finetune param */
  count = get_doubleordash(p, &opt_finetune_gtage);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  
  /* 2. gene tree spr finetune param */
  count = get_doubleordash(p, &opt_finetune_gtspr);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }


  /* 3. theta finetune param */
  count = get_doubleordash(p, &opt_finetune_theta);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 4. tau finetune param */
  count = get_doubleordash(p, &opt_finetune_tau);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 5. mixing finetune param */
  count = get_doubleordash(p, &opt_finetune_mix);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* TODO: The next is not implemented yet */
  double opt_finetune_seqerr;

  /* 6. locusrate finetune param */
  count = get_doubleordash(p, &opt_finetune_locusrate);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  #if 0
  /* 7. sequence error finetune param */
  count = get_doubleordash(p, &opt_finetune_seqerr);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }
  #endif

  /* 7. phi finetune */
  count = get_doubleordash(p, &opt_finetune_phi);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 8. freqs finetune */
  count = get_doubleordash(p, &opt_finetune_freqs);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 9. qrates finetune */
  count = get_doubleordash(p, &opt_finetune_qrates);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 10. alpha finetune */
  count = get_doubleordash(p, &opt_finetune_alpha);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 11. mubar finetune */
  count = get_doubleordash(p, &opt_finetune_mubar);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 12. nubar finetune */
  count = get_doubleordash(p, &opt_finetune_nubar);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 13. mu_i finetune */
  count = get_doubleordash(p, &opt_finetune_mui);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 14. nu_i finetune */
  count = get_doubleordash(p, &opt_finetune_nui);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 15. branchrates finetune */
  count = get_doubleordash(p, &opt_finetune_branchrate);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

  /* 16. migration rate finetune */
  count = get_doubleordash(p, &opt_finetune_migrates);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p))
  {
    ret = 1;
    goto l_unwind;
  }

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

  opt_print_qmatrix = 0;
  /* samples */
  count = get_long(p, &opt_print_samples);
  if (!count) goto l_unwind;

  p += count;

  /* check whether value is -1 in which case nothing else must follow */
  if (opt_print_samples == -1)
  {
    ret = 1;
    goto l_unwind;
  }

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

  /* optional qmatrix for backwards compatibility */
  count = get_long(p, &opt_print_qmatrix);
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
  species_count = spcount;

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

static long parse_migration(FILE * fp, const char * firstline, long line_count)
{
  long i;
  long ret = 1;
  long count;
  char * s = xstrdup(firstline);
  char * p = s;

  count = get_long(p, &opt_migration);
  if (!count) goto l_unwind;

  p += count;

  if (!opt_migration && is_emptyline(p)) goto l_unwind;

  ret = 0;

  if (!is_emptyline(p)) goto l_unwind;

  opt_mig_source = (char **)xcalloc((size_t)opt_migration, sizeof(char *));
  opt_mig_target = (char **)xcalloc((size_t)opt_migration, sizeof(char *));

  /* start reading potential migration between populations */
  for (i = 0; i < opt_migration; ++i)
  {
    if (!getnextline(fp))
      fatal("Incomplete 'migration' record (line %ld)", line_count+1);

    char * ss = xstrdup(line);
    p = ss;

    count = get_delstring(p, " \t\r\n*#,-", opt_mig_source+i);
    if (!count) goto l_unwind;
    p += count;

    count = get_delstring(p, " \t\r\n*#,-", opt_mig_target+i);
    if (!count) goto l_unwind;
    p += count;
    
    free(ss);
  }
  ret = 1;

  
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

static void update_sp_seqcount()
{
  long i;

  if (opt_diploid_size != species_count)
    fatal("Number of digits in 'diploid' does not match number of species");

  for (i = 0; i < species_count; ++i)
    if (opt_diploid[i])
      opt_sp_seqcount[i] *= 2;
 
}

static int cb_partitioncmp(const void * a, const void * b)
{
  const partition_t * pa = *((const partition_t **)a);
  const partition_t * pb = *((const partition_t **)b);

  if (pa->start < pb->start) return -1;
  else if (pa->start > pb->start) return 1;

  return 0;
}

/* transform list of partitions into an array for easier access */
static partition_t ** linearize_plist(list_t * plist, long * records)
{
  long i = 0;
  list_item_t * item;
  partition_t ** pa;

  *records = plist->count;

  /* allocate array */
  pa = (partition_t **)xmalloc((size_t)plist->count *
                               sizeof(partition_t *));

  item = plist->head;
  while (item)
  {
    partition_t * part = (partition_t *)(item->data);

    pa[i++] = part;

    item = item->next;
  }

  /* Sort partition by 'start' and check if complete range */
  qsort(pa,plist->count,sizeof(partition_t *), cb_partitioncmp);
  if (pa[0]->start != 1)
    fatal("File %s does not specify model for locus 1", opt_partition_file);
  for (i = 1; i < plist->count; ++i)
  {
    if (pa[i]->start != pa[i-1]->end+1)
    {
      fatal("File %s does not specify model(s) for loci %ld-%ld\n",
            opt_partition_file, pa[i-1]->end+1, pa[i]->start-1);
    }
  }


  return pa;
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

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
  {
    if (opt_theta_alpha <= 1)
      fatal("Alpha value of Inv-Gamma(a,b) of thetaprior must be > 1");

    if (opt_theta_beta <= 0)
      fatal("Beta value of Inv-Gamma(a,b) of thetaprior must be > 0");
  }
  else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
  {
    if (opt_theta_alpha <= 0)
      fatal("Alpha value of Gamma(a,b) of thetaprior must be > 0");

    if (opt_theta_beta <= 0)
      fatal("Beta value of Gamma(a,b) of thetaprior must be > 0");
  }
  else
  {
    assert(opt_theta_dist == BPP_THETA_PRIOR_BETA);
    if (opt_theta_p <= 0)
      fatal("Theta prior Beta(p,q,min,max) requires p > 0");
    if (opt_theta_q <= 0)
      fatal("Theta prior Beta(p,q,min,max) requires q > 0");
    if (opt_theta_min < 0)
      fatal("Theta prior Beta(p,q,min,max) requires min >= 0");
    if (opt_theta_max < 0)
      fatal("Theta prior Beta(p,q,min,max) requires max >= 0");
    if (opt_theta_max <= opt_theta_min)
      fatal("Theta prior Beta(p,q,min,max) requires max > min");
  }

  if (opt_migration)
  {
    if (opt_mig_alpha <= 0)
      fatal("Alpha value of Gamma(a,b) of migprior must be > 0");

    if (opt_mig_beta <= 0)
      fatal("Beta value of Gamma(a,b) of migprior must be > 0");
  }

  if (species_count > 1)
  {
    if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
    {
      if (opt_tau_alpha <= 1)
        fatal("Alpha value of Inv-Gamma(a,b) of tauprior must be > 1");

      if (opt_tau_beta <= 0)
        fatal("Beta value of Inv-Gamma(a,b) of tauprior must be > 0");
    }
    else
    {
      assert(opt_tau_dist == BPP_TAU_PRIOR_GAMMA);
      if (opt_tau_alpha <= 0)
        fatal("Alpha value of Gamma(a,b) of tauprior must be > 0");

      if (opt_tau_beta <= 0)
        fatal("Beta value of Gamma(a,b) of tauprior must be > 0");
    }
  }

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

  if (opt_theta_dist < BPP_THETA_PRIOR_MIN || opt_theta_dist > BPP_THETA_PRIOR_MAX)
    fatal("Internal error: invalid theta prior distribution");

  if (opt_theta_dist == BPP_THETA_PRIOR_INVGAMMA)
  {
    double invgammamean = opt_theta_beta / (opt_theta_alpha - 1);
    if (invgammamean > 1)
      fatal("Inverse gamma prior mean for thetas is > 1.\nPlease make sure you "
            "are indeed using Inv-Gamma as prior and not Gamma (bpp versions "
            "<= 3.3)");
  }
  else if (opt_theta_dist == BPP_THETA_PRIOR_GAMMA)
  {
    double gammamean = opt_theta_alpha / opt_theta_beta;
    if (gammamean > 1)
      fatal("Gamma prior mean for thetas is > 1.\nPlease make sure you "
            "are indeed using Gamma as prior and not Inv-Gamma");
  }

  if (opt_tau_dist == BPP_TAU_PRIOR_INVGAMMA)
  {
    double gammamean = opt_tau_beta / (opt_tau_alpha - 1);
    if (gammamean > 1)
      fatal("Inverse gamma prior mean for taus is > 1. Please make sure you "
            "are indeed using Inv-Gamma as prior and not Gamma (bpp versions "
            "<= 3.3)");
  }
  else
  {
    assert(opt_tau_dist == BPP_TAU_PRIOR_GAMMA);
    double gammamean = opt_tau_alpha / opt_tau_beta;
    if (gammamean > 1)
      fatal("Gamma prior mean for taus is > 1. Please make sure you "
            "are indeed using Gamma as prior");
  }

  if (opt_model == BPP_DNA_MODEL_CUSTOM)
  {
    assert(opt_partition_file);

    list_t * partition_list = parse_partition_file();
    if (!partition_list)
      fatal("Error when reading partition file %s:\n  %s",
            opt_partition_file, errmsg);

    
    validate_partitions(partition_list);

    /* transform partition_list into an array for easier access */
    opt_partition_list = linearize_plist(partition_list,&opt_partition_count);

    /* now clear list but do not delete list items as they were placed in
       opt_partition_list */
    list_clear(partition_list,NULL);
    free(partition_list);
  }
  else
  {
    assert(!opt_partition_file);
  }

  /* if cleandata is set to 1 and the sequences for at least one species require
     phasing then print an error */
  if (opt_cleandata && opt_diploid)
    fatal("ERROR: Cannot phase sequences while cleandata flag is set (remove "
          "ambiguities).\nPlease set cleandata=0 when phasing.");

  /* check clock and locusrate/branchrate */
  if (opt_clock < BPP_CLOCK_MIN || opt_clock > BPP_CLOCK_MAX)
    fatal("Invalid 'clock' value");

  if (opt_datefile && (opt_method == METHOD_10 || opt_method == METHOD_11))
    fatal("Cannot use species delimitation models when using tip dating.");


}

static void update_locusrate_information()
{

    if (opt_locusrate_prior == BPP_LOCRATE_PRIOR_GAMMADIR)
    {
      if (opt_mubar_alpha == 0 && opt_mubar_beta == 0)
        opt_locusrate_prior = BPP_LOCRATE_PRIOR_DIR;
    }

    if (opt_est_locusrate == MUTRATE_ESTIMATE &&
        opt_locusrate_prior == BPP_LOCRATE_PRIOR_HIERARCHICAL)
    {
      /* if mubar_alpha = mubar_beta = 0 then fix (do not estimate) mubar */
      if (opt_mubar_alpha == 0 && opt_mubar_beta == 0)
        opt_est_mubar = 0;
      else
        opt_est_mubar = 1;
    }
}


static void set_print_locusfile()
{
  if (opt_print_rates || opt_print_locusrate ||
      opt_print_hscalars || opt_print_qmatrix)
    opt_print_locusfile = 1;
}

static void set_debug_flags()
{
  long * flagptr;
  long i;

  long * flags[] = {
    &opt_debug_sim,   &opt_debug_gage, &opt_debug_gspr, &opt_debug_mui,
    &opt_debug_hs,    &opt_debug_mix,  &opt_debug_rj,   &opt_debug_theta,
    &opt_debug_tau,   &opt_debug_sspr, &opt_debug_br,   &opt_debug_snl,
    &opt_debug_parser, NULL
  };

  if (opt_debug)
  {
    for (i = 0; flags[i]; ++i)
    {
      flagptr = flags[i];

      if (!(*flagptr))
        *flagptr = opt_debug;
    }
  }
  else
  {
    long flag_set = 0;
    for (i = 0; flags[i]; ++i)
    {
      flagptr = flags[i];
      flag_set |= !!(*flagptr);
    }
    opt_debug = flag_set;
  }
}

void load_cfile()
{
  long line_count = 0;
  FILE * fp;

  /* the following variable is used for checking whether we have a newick
     string in the species&tree tag, in the case of 1 species. For species
     trees we do not accept a tree, whereas for network we require a newick
     string. The program always reads a line. If that line is a tree it is
     processed, otherwise this variable is set such that we do not read another
     line */
  long line_not_processed = 0;

  fp = xopen(opt_cfile,"r");

  while (line_not_processed || getnextline(fp))
  {
    int valid = 0;
    char * token;
    char * value;
    long token_len;

    line_not_processed = 0;

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
        valid = 1;
      }
      else if (!strncasecmp(token,"arch",4))
      {
        char * temp;
        if (!get_string(value,&temp))
          fatal("Option %s expects a string (line %ld)", token, line_count);

        if (!strcasecmp(temp,"cpu"))
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
          fatal("Option 'print' expects either four bits or '-1'  and a fifth optional bit(line %ld)",
                line_count);

        if (opt_print_samples == 0)
          fatal("First bit of 'print' must be set to 1");

        if (opt_print_samples == -1)
          opt_onlysummary = 1;

        valid = 1;
      }
      else if (!strncasecmp(token,"model",5))
      {
        if (!parse_model(value))
          fatal("Erroneous format of 'model' option (line %ld)", line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"clock",5))
      {
        if (!parse_clock(value))
          fatal("Erroneous format of 'clock' (line %ld)\n"
                "Syntax:\n"
                "  clock = 1                                # strict clock\n"
                "  clock = 2 a_vbar b_vbar a_vi prior dist  # relaxed clock",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"phase",5))
      {
        if (!parse_diploid(value))
          fatal("Option %s expects values 0 or 1 for each species (line %ld)",
                token,line_count);
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
        fatal("Option 'diploid' was renamed to 'phase' (same syntax).\n"
              "Please update the control file (line %ld).", line_count);
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
      else if (!strncasecmp(token,"scaling",7))
      {
        if (!parse_long(value,&opt_scaling) ||
            (opt_scaling != 0 && opt_scaling != 1))
          fatal("Option 'scaling' expects value 0 or 1 (line %ld)", line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"threads",7))
      {
        #if 0
        if (!parse_long(value,&opt_threads) || (opt_threads <= 0))
          fatal("Option 'threads' requires a positive integer (line %ld)",
                line_count);
        valid = 1;
        #endif
        if (!parse_threads(value))
          fatal("Option 'threads' expects an integer (line %ld)",
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
      else if (!strncasecmp(token,"datefile",8)) {
        if (!get_string(value, &opt_datefile))
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
          fatal("Option 'tauprior' (line %ld) expects the following syntax:\n"
                "  tauprior = invgamma alpha beta     # for inverse gamma prior\n"
                "  tauprior = gamma alpha beta        # for gamma prior\n",
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
      else if (!strncasecmp(token,"phiprior",8))
      {
        if (!parse_phiprior(value))
          fatal("Option 'phiprior' expects two doubles (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"migprior",8))
      {
        if (!parse_migprior(value))
          fatal("Option 'migprior' expects two doubles (line %ld)", line_count);
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
      else if (!strncasecmp(token,"migration",9))
      {
        if (!parse_migration(fp,value,line_count))
          fatal("Erroneous format of 'migration' (line %ld)", line_count);
        valid = 1;
      }
    }
    else if (token_len == 10)
    {
      if (!strncasecmp(token,"thetaprior",10))
      {
        if (!parse_thetaprior(value))
          fatal("Option 'thetaprior' (line %ld) expects the following syntax:\n"
                "  thetaprior = invgamma alpha beta     # for inverse gamma prior\n"
                "  thetaprior = gamma alpha beta        # for gamma prior\n"
                "  thetaprior = beta p q min max        # for beta prior\n",
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
      else if (!strncasecmp(token,"alphaprior",10))
      {
        if (!parse_alphaprior(value))
          fatal("Option 'alphaprior' expects two doubles (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"thetamodel",10))
      {
        if (!parse_thetamodel(value))
          fatal("Erroneous format of 'thetamodel' (line %ld)\n"
                "Possible options:\n"
                " linked-none\n linked-all\n linked-inner\n linked-msci",
                line_count);
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
      else if (!strncasecmp(token,"loadbalance",11))
      {
        if (!parse_loadbalance(value))
          fatal("Invalid load balance option (line %ld)", line_count);
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

        if (spcount > 1)
        {
          if (!getnextline(fp))
            fatal("Incomplete 'species&tree' record (line %ld)", line_count);

          ++line_count;

          if (!get_string(line,&opt_streenewick))
            fatal("Expected newick tree string in 'species&tree' (line %ld)",
                   line_count);
        }
        else if (spcount == 1)
        {
          /* TODO: This is an ugly hack to account for the case where we have 1
             species and a network */
          int reached_eof = 0;
          if (!getnextline(fp))
            reached_eof = 1; 

          ++line_count;
          line_not_processed = 1;

          if (!reached_eof && starts_with_opar(line))
          {
            if (!get_string(line,&opt_streenewick))
              fatal("Expected newick string in 'species&tree' (line %ld)",
                    line_count);
            
            line_not_processed = 0;
          }
          else
          {
            opt_streenewick = (char *)xmalloc((size_t)(strlen(opt_reorder)+2) *
                                              sizeof(char));
            strcpy(opt_streenewick, opt_reorder);
            opt_streenewick[strlen(opt_reorder)] = ';';
            opt_streenewick[strlen(opt_reorder)+1] = '\0';
          }

          if (reached_eof)
            break;
        }
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
    else if (token_len == 14)
    {
      if (!strncasecmp(token,"constraintfile",14))
      {
        if (!get_string(value,&opt_constraintfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
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

  fclose(fp);

  /* set method */
  if (!opt_est_stree && !opt_est_delimit)
    opt_method = METHOD_00;
  else if (!opt_est_stree)
    opt_method = METHOD_10;
  else if (!opt_est_delimit)
    opt_method = METHOD_01;
  else
    opt_method = METHOD_11;

  opt_snl_lambda_expand = log(opt_snl_lambda_expand) / log(1 - opt_snl_lambda_expand);
  opt_snl_lambda_shrink = log(opt_snl_lambda_shrink) / log(1 - opt_snl_lambda_shrink);

  update_locusrate_information();
  check_validity();

  /* decide whether per-locus files for sampling are required */
  set_print_locusfile();
  set_debug_flags();

  if (opt_diploid)
    update_sp_seqcount();

  if (species_count == 1 && opt_method != METHOD_00)
    fatal("You can only use method A00 with one species");
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
