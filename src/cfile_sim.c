/*
    Copyright (C) 2016-2019 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static long species_count = 0;

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

static long get_tree_string_with_thetas(const char * line, char ** value)
{
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* find end of tree string */
  char * tree_end = strchr(p, ';');
  if (!tree_end)
  {
    free(s);
    return 0;
  }

  /* if there are comments after tree string, delete them */
  tree_end = tree_end+1;
  if (*tree_end != 0)
  {
    char * end = strchr(tree_end,'#');
    if (end)
      *end = 0;
    end = strchr(tree_end, '*');
    if (end)
      *end = 0;
  }


  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line */
  if (!p[ws])
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters to find end of string */
  while (*p) p++;


  /* now go back, skipping all white-space until a character occurs */
  for (--p; *p == ' ' || *p == '\t' || *p == '\r' || *p =='\n'; --p);
  
  char * end = p+1;
  *end = 0;

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

#define DNA_QRATES_COUNT        6
#define DNA_STATES_COUNT        4

static long parse_basefreqs(const char * line)
{
  long i;
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_basefreqs_fixed);
  if (!count) goto l_unwind;

  p += count;

  if  (opt_basefreqs_fixed < 0 || opt_basefreqs_fixed > 1)
    fatal("Option 'basefreqs' must start with a '0' or '1'");

  opt_basefreqs_params = (double *)xcalloc(DNA_STATES_COUNT,sizeof(double));
  for (i = 0; i < DNA_STATES_COUNT; ++i)
  {
    count = get_double(p, opt_basefreqs_params+i);
    if (!count) goto l_unwind;

    p += count;
  }
  
  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_qrates(const char * line)
{
  long i;
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_qrates_fixed);
  if (!count) goto l_unwind;

  p += count;

  if  (opt_qrates_fixed < 0 || opt_qrates_fixed > 1)
    fatal("Option 'qrates' must start with a '0' or '1'");

  opt_qrates_params = (double *)xcalloc(DNA_QRATES_COUNT,sizeof(double));
  for (i = 0; i < DNA_QRATES_COUNT; ++i)
  {
    count = get_double(p, opt_qrates_params+i);
    if (!count) goto l_unwind;

    p += count;
  }
  
  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long parse_loci_and_lengths(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_locus_count);
  if (!count) goto l_unwind;

  p += count;

  /* now read second token */
  count = get_long(p, &opt_locus_simlen);
  if (!count) goto l_unwind;

  p += count;

  if (!is_emptyline(p)) goto l_unwind;
  
  if (opt_locus_count > 0 && opt_locus_simlen > 0)
    ret = 1;
  
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

static char ** split_strings(const char * x, long * token_count)
{
  long i;
  size_t ws;
  char * s = xstrdup(x);
  char * p = s;

  /* 1. get the data part of the line, i.e. skip all leading/trailing
     whitespace and potential comments */

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

  char * data = xstrdup(start);
  free(s);

  /* 2. Now split data into strings and use whitespace as delimiter */
  
  
  /* 2a. Go through data and count the number of tokens */
  *token_count = 0;
  p = data;
  while (*p)
  {
    /* skip characters */
    while (*p && *p != '\t' && *p != ' ') p++;

    /* increment token count */
    (*token_count)++;

    /* skip whitespace */
    while (*p && (*p == ' ' || *p == '\t')) p++;
  }

  /* allocate space for storing tokens */
  char ** token = (char **)xmalloc((size_t)(*token_count) * sizeof(char *));

  /* 2b. Now extract the tokens and store them in a list */
  p = data;
  for (i = 0; i < *token_count; ++i)
  {
    start = p;
    while (*p && *p != '\t' && *p != ' ') p++;

    /* copy token into list */
    char temp = *p;
    *p = 0;
    token[i] = xstrdup(start);
    *p = temp;

    /* skip whitespace */
    while (*p && (*p == ' ' || *p == '\t')) p++;
  }

  return token;
}

static void parse_migration_matrix(FILE * fp, long line_count)
{
  long i,j;

  size_t matrix_size = opt_migration * opt_migration;


  opt_migration_matrix = (double *)xmalloc(matrix_size*sizeof(double));
  opt_migration_events = (double *)xcalloc(matrix_size,sizeof(double));

  if (!getnextline(fp))
    fatal("Incomplete 'migration' record (line %ld)", line_count+1);

  long matrix_dim = 0;
  opt_migration_labels = split_strings(line,&matrix_dim);
  if (!opt_migration_labels || matrix_dim == 0)
    fatal("Option 'migration' must be followed by the population labels (line %ld)", line_count+1);

  /* read matrix */
  for (i = 0; i < matrix_dim; ++i)
  {
    long dim;

    if (!getnextline(fp))
      fatal("Incomplete 'migration' record (line %ld)", line_count+2+i);

    char ** data = split_strings(line,&dim);

    if (dim != matrix_dim+1)
      fatal("Wrong number of parameters in migration matrix row (line %ld)"
            " Expected %ld found %ld", line_count+2+i, matrix_dim+1, dim);
    
    if (strcmp(data[0],opt_migration_labels[i]))
      fatal("Migration matrix label of row %ld does not match label of column "
            "%ld (line %ld)", i+1, i+1, line_count+2+i);

    for (j = 0; j < matrix_dim; ++j)
      if (!get_double(data[j+1], opt_migration_matrix+i*matrix_dim+j))
        fatal("Migration matrix cell (%ld,%ld) is not a number (line %ld)",
              i+1, j+1, line_count+2+i);

  }
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

static long parse_siterate(const char * line)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, &opt_siterate_fixed);
  if (!count) goto l_unwind;

  p += count;

  if (opt_siterate_fixed < 0 || opt_siterate_fixed > 1)
    fatal("Option 'alpha_siterate' expects values '0' or '1'");

  if (opt_siterate_fixed == 0)
  {
    /* sampled gamma shape parameter alpha */
    count = get_double(p, &opt_siterate_alpha);
    if (!count) goto l_unwind;

    p += count;

    count = get_double(p, &opt_siterate_beta);
    if (!count) goto l_unwind;

    p += count;

    count = get_long(p, &opt_siterate_cats);
    if (!count) goto l_unwind;

    p += count;
  }
  else
  {
    /* fixed gamma shape parameter alpha */
    count = get_double(p, &opt_siterate_alpha);
    if (!count) goto l_unwind;

    p += count;

    if (is_emptyline(p))
    {
      if (opt_siterate_alpha == 0)
        ret = 1;

      goto l_unwind;
    }

    count = get_long(p, &opt_siterate_cats);
    if (!count) goto l_unwind;

    p += count;
  }

  if (is_emptyline(p)) ret = 1;

l_unwind:
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

  if (opt_clock < 1 || opt_clock > 3)
    fatal("Option 'clock' expects values '1', '2' or '3'");

  if (opt_clock > 1)
  {
    count = get_double(p, &opt_clock_alpha);
    if (!count) goto l_unwind;

    p += count;
  }

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
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

static void check_validity()
{
  long i;

  if (!opt_streenewick)
    fatal("Initial species tree newick format is required in 'species&tree'");

  if (!opt_mapfile)
    fatal("Option 'imapfile' is required");

  if (!opt_locus_count || !opt_locus_simlen)
    fatal("Option 'loci&length' is required");

  if (opt_model == -1)
  {
    fprintf(stdout, "Substitution model not specified, assuming JC69\n");
    opt_model = BPP_DNA_MODEL_JC69;
  }
  assert(opt_model == BPP_DNA_MODEL_JC69 || opt_model == BPP_DNA_MODEL_GTR);

  if (opt_locusrate_alpha < 0)
    fatal("Option 'alpha_locusrate' requires a positive number or 0");

  if (opt_model == BPP_DNA_MODEL_GTR && !opt_modelparafile)
    fatal("Usage of 'GTR' model requires the specification of an output file "
          "with the option 'modelparafile'");

  if (opt_model == BPP_DNA_MODEL_GTR)
  {
    if (opt_basefreqs_fixed == -1)
      fatal("Usage of 'GTR' model requires specifying base frequencies with "
            "the 'basefreqs' option");
    if (opt_basefreqs_fixed == 1)
    {
      double sum = 0;
      for (i = 0; i < 4; ++i)
      {
        if (opt_basefreqs_params[i] < 0)
          fatal("Base frequencies cannot be negative");
        sum += opt_basefreqs_params[i];
      }
      if (sum != 1)
        fatal("Base frequencies must sum to 1");
    }
    else
    {
      for (i = 0; i < 4; ++i)
        if (opt_basefreqs_params[i] <= 0)
          fatal("Alpha parameters for dirichlet distribution must be positive "
                "numbers (option 'basefreqs')");
    }
  }

  if (!opt_siterate_fixed)
  {
    if (opt_siterate_alpha <= 0 || opt_siterate_beta <= 0)
      fatal("Alpha and beta parameters for option 'siterate' must be positive "
            "numbers");
    fprintf(stdout, "Site rates: alpha sampled from Gamma(%f,%f), K = %ld\n",
            opt_siterate_alpha, opt_siterate_beta, opt_siterate_cats);
  }
  else
  {
    if (opt_siterate_cats)
      fprintf(stdout, "Site rates: from discrete Gamma(%f), K = %ld\n",
              opt_siterate_alpha, opt_siterate_cats);
    else
      fprintf(stdout, "Site rates: from discrete Gamma(%f), K = inf\n",
              opt_siterate_alpha);

  }

  if (opt_migration)
  {
    if (opt_migration != species_count*2-1)
      fatal("Option 'migration' must be equal number of nodes in species tree");
  }
}

void load_cfile_sim()
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

  fp = xopen(opt_simulate,"r");

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
            opt_simulate, line_count);
    
    if (token_len == 4)
    {
      if (!strncasecmp(token,"seed",4))
      {
        if (!get_long(value,&opt_seed))
          fatal("Option 'seed' expects one integer (line %ld)", line_count);

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
      if (!strncasecmp(token,"clock",5))
      {
        if (!parse_clock(value))
          fatal("Option 'clock' expects values '1', '2 a' or '3 a' (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"model",5))
      {
        if (!get_long(value,&opt_model) || opt_model < 0)
          fatal("Option 'model' expects value '%d' or '%d' (line %ld)",
                BPP_DNA_MODEL_JC69, BPP_DNA_MODEL_GTR, line_count);
        valid = 1;
      }
    }
    else if (token_len == 6)
    {
      if (!strncasecmp(token,"qrates",6))
      {
        if (!parse_qrates(value))
          fatal("Option 'qrates' expects one switch and 6 values (line %ld)",
                line_count);
        valid = 1;
      }
    }
    else if (token_len == 7)
    {
      if (!strncasecmp(token,"seqfile",7))
      {
        if (!get_string(value, &opt_msafile))
          fatal("Option '%s' expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"diploid",7))
      {
        if (!parse_diploid(value))
          fatal("Option %s expects values 0 or 1 for each species (line %ld)",
                token,line_count);
        valid = 1;
      }
    }
    else if (token_len == 8)
    {
      if (!strncasecmp(token,"treefile",8))
      {
        if (!get_string(value, &opt_treefile))
          fatal("Option 'treefile' expects a string (line %ld)", line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"imapfile",8))
      {
        if (!get_string(value, &opt_mapfile))
          fatal("Option 'imapfile' expects a string (line %ld)", line_count);
        valid = 1;
      }
    }
    else if (token_len == 9)
    {
      if (!strncasecmp(token,"basefreqs",9))
      {
        if (!parse_basefreqs(value))
          fatal("Option 'basefreqs' expects one switch and 4 values (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"migration",9))
      {
        if (!get_long(value,&opt_migration))
          fatal("Option 'migration' expects one integer (line %ld)", line_count);
        
        parse_migration_matrix(fp, line_count);
        valid = 1;
      }
    }
    else if (token_len == 10)
    {
      if (!strncasecmp(token,"concatfile",10))
      {
        if (!get_string(value,&opt_concatfile))
          fatal("Option 'concatfile' expects a string (line %ld)", line_count);
        valid = 1;
      }
    }
    else if (token_len == 11)
    {
      if (!strncasecmp(token,"loci&length",11))
      {
        if (!parse_loci_and_lengths(value))
          fatal("Option 'loci&length' expects two positive integer values (line %ld)",
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

        if (spcount > 1)
        {
          if (!getnextline(fp))
            fatal("Incomplete 'species&tree' record (line %ld)", line_count);

          ++line_count;

          if (!get_tree_string_with_thetas(line,&opt_streenewick))
            fatal("Expected newick tree string in 'species&tree' (line %ld) "
                  "with ending ';' character", line_count);
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
            if (!get_tree_string_with_thetas(line,&opt_streenewick))
              fatal("Expected newick string in 'species&tree' (line %ld) "
                    "with ending ';' character", line_count);
            
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
      if (!strncasecmp(token,"modelparafile",13))
      {
        if (!get_string(value,&opt_modelparafile))
          fatal("Option 'modelparafile' expects a string (line %ld)",
                line_count);
        valid = 1;
      }
    }
    else if (token_len == 14)
    {
      if (!strncasecmp(token,"alpha_siterate",14))
      {
        if (!parse_siterate(value))
          fatal("Erroneous format of 'alpha_siterate' (line %ld)", line_count);
        valid = 1;
      }
    }
    else if (token_len == 15)
    {
      if (!strncasecmp(token,"alpha_locusrate",15))
      {
        if  (!get_double(value, &opt_locusrate_alpha))
          fatal("Option 'alpha_locusrate' expects one value (line %ld)",
                line_count);
        if (opt_locusrate_alpha < 0)
          fatal("Option 'alpha_locusrate' expects a non-negative value (line %ld)",
                line_count);
        if (opt_locusrate_alpha == 0)
          opt_est_locusrate = 0;
        else
          opt_est_locusrate = 1;
        valid = 1;
      }
    }

    if (!valid)
      fatal("Invalid syntax when parsing file %s on line %ld",
            opt_simulate, line_count);
  }

  check_validity();

  #if 0
  if (opt_diploid)
    update_sp_seqcount();
  #endif

  fclose(fp);
}

