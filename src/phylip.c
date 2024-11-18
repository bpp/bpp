/*
    Copyright (C) 2016-2024 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

#define PHYLIP_SEQUENTIAL  1
#define PHYLIP_INTERLEAVED 2

static int dfa_parse(phylip_t * fd,
                     msa_t * msa,
                     char * p,
                     int seqno,
                     int offset)
{
  int j = 0;
  char c,m;

  char * seqdata = msa->sequence[seqno] + offset;

  /* read sequence data */
  while((c = *p++))
  {
    m = (char) fd->chrstatus[(int)c];
    switch(m)
    {
      case 0:
        /* characters to be stripped */
        fd->stripped_count++;
        fd->stripped[(int)c]++;
        break;

      case 1:
        /* legal character */
        if (offset + j >= msa->length)
        {
          bpp_errno = ERROR_PHYLIP_LONGSEQ;
          snprintf(bpp_errmsg, 200, "Sequence %d (%.100s) longer than expected",
                   seqno+1, msa->label[seqno]);
          return -1;
        }
        seqdata[j++] = c;
        break;

      case 2:
        /* fatal character */
        if (c>=32)
        {
          bpp_errno = ERROR_PHYLIP_ILLEGALCHAR;
          snprintf(bpp_errmsg, 200, "illegal character '%c' "
                                    "on line %ld in the fasta file",
                                    c, fd->lineno);
        }
        else
        {
          bpp_errno = ERROR_PHYLIP_UNPRINTABLECHAR;
          snprintf(bpp_errmsg, 200, "illegal unprintable character "
                                    "%#.2x (hexadecimal) on line %ld "
                                    "in the fasta file",
                                    c, fd->lineno);
        }
        return -1;

      case 3:
        /* silently stripped chars */
        break;
    }
  }
  return j;
}

static char * reallocline(phylip_t * fd, size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  memcpy(temp,fd->line,fd->line_size*sizeof(char));
  free(fd->line);
  fd->line = temp;
  fd->line_maxsize = newmaxsize;

  return temp;
}

static char * getnextline(phylip_t * fd)
{
  size_t len = 0;

  fd->line_size = 0;

  /* read from file until newline or eof */
  while (fgets(fd->buffer, LINEALLOC, fd->fp))
  {
    len = strlen(fd->buffer);

    if (fd->line_size + len > fd->line_maxsize)
      if (!reallocline(fd, fd->line_maxsize + LINEALLOC))
        return NULL;

    memcpy(fd->line+fd->line_size,fd->buffer,len*sizeof(char));
    fd->line_size += len;

    if (fd->buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        fd->line[fd->line_size-1] = 0;
      #endif

      return fd->line;
    }
  }

  if (!fd->line_size)
  {
    free(fd->line);
    fd->line = NULL;
    return NULL;
  }

  if (fd->line_size == fd->line_maxsize)
    if (!reallocline(fd,fd->line_maxsize+1))
      return NULL;

  fd->line[fd->line_size] = 0;
  return fd->line;

}

static int args_getint(const char * arg, int * len)
{
  int temp;
  *len = 0;
  
  int ret = sscanf(arg, "%d%n", &temp, len);
  if ((ret == 0) || (!*len)) 
    return 0;

  return temp;
}

static int whitespace(char c)
{
  if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
    return 1;
  return 0;
}


static int parse_header(const char * line,
                        int * seq_count,
                        int * seq_len,
                        int format)
{
  int len;

  /* read number of sequences */
  if (!(*seq_count = args_getint(line,&len)))
  {
    bpp_errno = ERROR_PHYLIP_SYNTAX;
    snprintf(bpp_errmsg, 200, "Invalid number of sequences in header");
    return BPP_FAILURE;
  }

  line += len;

  /* read sequence length */
  if (!(*seq_len = args_getint(line,&len)))
  {
    bpp_errno = ERROR_PHYLIP_SYNTAX;
    snprintf(bpp_errmsg, 200, "Invalid sequence length in header");
    return BPP_FAILURE;
  }

  line += len;

  /* go through all white spaces */
  while (*line && whitespace(*line)) ++line;

  /* if end of line then return successfully */
  if (!*line)
    return 1;

  /* otherwise, continue only if interleaved format specified, otherwise die */
  if (format == PHYLIP_SEQUENTIAL)
    return 0;

  if (*line != 's' && *line != 'S' && *line != 'i' && *line != 'I')
    return 0;

  /* go through all white spaces */
  while (*line && whitespace(*line)) ++line;

  /* if end of line then return successfully */
  if (!*line)
    return 1;

  return 0;
}

static char * parse_oneline_sequence(phylip_t * fd,
                                     msa_t * msa,
                                     char * p,
                                     int seqno,
                                     int offset,
                                     int * aln_len,
                                     int * error)
{
  int j = 0;

  while (p && !j)
  {
    /* read data */
    if ((j = dfa_parse(fd,msa,p,seqno,offset)) == -1)
    {
      *error = 1;
      return NULL;
    }

    if (j)
    {
      if (!(*aln_len))
      {
        *aln_len = j;
      }
      else if (*aln_len != j)
      {
        *error = 1;
        bpp_errno = ERROR_PHYLIP_NONALIGNED;
        snprintf(bpp_errmsg, 200, "Sequence %d (%.100s) data out of alignment",
                 seqno+1, msa->label[seqno]); 
        return NULL;
      }
    }
    else
      p = getnextline(fd);
  }

  return p;
}

phylip_t * phylip_open(const char * filename,
                       const unsigned int * map)
{
  int i;

  phylip_t * fd = (phylip_t *)xmalloc(sizeof(phylip_t));

  /* allocate space */
  fd->line = NULL;
  fd->line_size = 0;
  fd->line_maxsize = 0;

  fd->lineno = 0;

  fd->no = -1;

  fd->chrstatus = map;

  /* open file */
  fd->fp = fopen(filename, "r");
  if (!(fd->fp))
    fatal("Unable to open file (%s)", filename);

  /* get filesize */
  if (fseek(fd->fp, 0, SEEK_END))
    fatal("Unable to seek in file (%s)", filename);

  fd->filesize = ftell(fd->fp);

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for(i=0; i<256; i++)
    fd->stripped[i] = 0;

  /* cache line */
  if (!getnextline(fd))
  {
    if (fd->line)
      free(fd->line);
    fclose(fd->fp);
    free(fd);
    return NULL;
  }

  fd->lineno = 1;

  return fd;
}

int phylip_rewind(phylip_t * fd)
{
  int i;

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for(i=0; i<256; i++)
    fd->stripped[i] = 0;

  if (!getnextline(fd))
    fatal("Unable to rewind and cache data");

  fd->lineno = 1;
  fd->no = -1;

  return BPP_SUCCESS;
}

void phylip_close(phylip_t * fd)
{
  fclose(fd->fp);
  if (fd->line)
    free(fd->line);
  free(fd);
}

static long emptyline(const char * line)
{
  size_t ws = strspn(line, " \t\r\n");
  if (!line[ws]) return 1;
  return 0;
}

msa_t * phylip_parse_interleaved(phylip_t * fd)
{
  int i;
  int aln_len;
  int sumlen;
  int seqno;
  long headerlen;

  msa_t * msa = (msa_t *)xmalloc(sizeof(msa_t));

  while (emptyline(fd->line)) getnextline(fd);

  /* read header */
  if (!parse_header(fd->line,
                    &(msa->count),
                    &(msa->length),
                    PHYLIP_INTERLEAVED))
    return NULL;

  /* allocate msa placeholders */
  msa->sequence = (char **)xcalloc((size_t)(msa->count),sizeof(char *));
  msa->label = (char **)xcalloc((size_t)(msa->count),sizeof(char *));

  /* allocate sequence data placeholders */
  for (i = 0; i < msa->count; ++i)
  {
    msa->sequence[i] = (char *)xmalloc((size_t)(msa->length+1) * sizeof(char));
    msa->sequence[i][msa->length] = 0;
  }

  /* read sequences with headers */
  seqno = 0;
  aln_len = 0;
  int error = 0;
  while (1)
  {
    /* get next line */
    char * p = getnextline(fd);

    /* if no more lines break */
    if (!p) break;

    /* skip whitespace before sequence header */
    while (*p && whitespace(*p)) ++p;

    /* restart loop if blank line */
    if (!*p) continue;

    /* error if there are more sequences than specified */
    if (seqno == msa->count)
    {
      bpp_errno = ERROR_PHYLIP_SYNTAX;
      snprintf(bpp_errmsg, 200, "Found at least %d sequences but expected %d",
               seqno+1, msa->count);
      msa_destroy(msa);
      return NULL;
    }

    /* find first blank after header */
    if (strchr(p,' '))
      headerlen = xstrchrnul(p,' ') - p;
    else if (strchr(p,'\t'))
      headerlen = xstrchrnul(p,'\t') - p;
    else if (strchr(p,'\r'))
      headerlen = xstrchrnul(p,'\r') - p;
    else
      headerlen = xstrchrnul(p,'\n') - p;

    /* headerlen cannot be zero */
    assert(headerlen > 0);

    /* store sequence header */
    msa->label[seqno] = (char *)xmalloc((size_t)(headerlen+1)*sizeof(char));
    memcpy(msa->label[seqno], p, (size_t)headerlen);
    msa->label[seqno][headerlen] = 0;

    p += headerlen;

    /* read (and parse) the first line (starting from p) that contains at
       least one character */
    if (!parse_oneline_sequence(fd,msa,p,seqno,0,&aln_len,&error))
      break;

    ++seqno;

    if (seqno == msa->count)
      break;
  }

  /* was the last block of sequences non-aligned? */
  if (error)
  {
    msa_destroy(msa);
    return NULL;
  }

  if (seqno != msa->count)
  {
    bpp_errno = ERROR_PHYLIP_SYNTAX;
    snprintf(bpp_errmsg, 200, "Found %d sequence(s) but expected %d",
             seqno, msa->count);
    msa_destroy(msa);
    return NULL;
  }

  /* update the length of the alignment read so far, which will be used as the
     offset when appending data to the end of the sequences */
  sumlen = aln_len;

  /* now read the remaining blocks */
  seqno = 0;
  aln_len = 0;
  int block_count = 2;
  while (1)
  {
    char * p = getnextline(fd);

    /* read (and parse) the first line (starting from p) that contains at
       least one character */
    if (!parse_oneline_sequence(fd,msa,p,seqno,sumlen,&aln_len,&error))
      break;
    
    seqno = (seqno+1) % msa->count;

    /* if data for all sequences were read, then append the alignment length
       to the sum, and go for the next block */
    if (!seqno)
    {
      sumlen += aln_len;
      aln_len = 0;
      block_count++;
    }
  }

  /* was the last block of sequences non-aligned? */
  if (error)
  {
    msa_destroy(msa);
    return NULL;
  }

  /* if seqno != 0 then there were more (or less) sequences than expected */
  if (seqno)
  {
    bpp_errno = ERROR_PHYLIP_SYNTAX;
    snprintf(bpp_errmsg, 200, "Found %d sequences in block %d but expected %d",
             seqno, block_count, msa->count);
    msa_destroy(msa);
    return NULL;
  }
  if (sumlen != msa->length)
  {
    snprintf(bpp_errmsg, 200, "Sequence length is %d but expected %d",
             sumlen, msa->length);
    msa_destroy(msa);
    return NULL;
  }

  return msa;
}

msa_t * phylip_parse_sequential(phylip_t * fd)
{
  int i,j;
  long headerlen;

  msa_t * msa = (msa_t *)xcalloc(1,sizeof(msa_t));

  while (emptyline(fd->line)) getnextline(fd);
    
  /* read header */
  if (!parse_header(fd->line,
                    &(msa->count),
                    &(msa->length),
                    PHYLIP_SEQUENTIAL))
    return NULL;

  msa->sequence = (char **)xcalloc((size_t)(msa->count),sizeof(char *));
  msa->label = (char **)xcalloc((size_t)(msa->count),sizeof(char *));

  for (i = 0; i < msa->count; ++i)
  {
    msa->sequence[i] = (char *)xmalloc((size_t)(msa->length+1) * sizeof(char));
    msa->sequence[i][msa->length] = 0;
  }
  
  /* read sequences */
  int seqno = 0;
  while (1)
  {
    /* get next line */
    char * p = getnextline(fd);

    /* if no more lines break */
    if (!p) break;

    /* skip whitespace before sequence header */
    while (*p && whitespace(*p)) ++p;

    /* restart loop if blank line */
    if (!*p) continue;

    /* error if there are more sequences than specified */
    if (seqno == msa->count)
    {
      bpp_errno = ERROR_PHYLIP_SYNTAX;
      snprintf(bpp_errmsg, 200, "Found at least %d sequences but expected %d",
               seqno+1, msa->count);
      msa_destroy(msa);
      return NULL;
    }

    /* find first blank after header */
    if (strchr(p,' '))
      headerlen = xstrchrnul(p,' ') - p;
    else if (strchr(p,'\t'))
      headerlen = xstrchrnul(p,'\t') - p;
    else if (strchr(p,'\r'))
      headerlen = xstrchrnul(p,'\r') - p;
    else
      headerlen = xstrchrnul(p,'\n') - p;

    /* headerlen cannot be zero */
    assert(headerlen > 0);

    /* store sequence header */
    msa->label[seqno] = (char *)xmalloc((size_t)(headerlen+1)*sizeof(char));
    memcpy(msa->label[seqno], p, (size_t)headerlen);
    msa->label[seqno][headerlen] = 0;

    p += headerlen;

    /* go through possibly multiple sequence data lines */
    j=0;
    while (1)
    {
      /* read sequence data */
      int chars_count = dfa_parse(fd,msa,p,seqno,j);
      if (chars_count == -1)
      {
        msa_destroy(msa);
        return NULL;
      }

      j += chars_count;

      /* break if we read all sequence data */
      if (j == msa->length)
        break;

      p = getnextline(fd);

      if (!p)
      {
        bpp_errno = ERROR_PHYLIP_SYNTAX;
        snprintf(bpp_errmsg, 200,
                 "Sequence %d (%.100s) has %d characters but expected %d",
                seqno+1,msa->label[seqno],j,msa->length);
        msa_destroy(msa);
        return NULL;
      }
    }

    ++seqno;
    /* TODO: Updated for BPP */
    if (seqno == msa->count)
      return msa;
  }
  
  if (seqno != msa->count)
  {
    bpp_errno = ERROR_PHYLIP_SYNTAX;
    snprintf(bpp_errmsg, 200, "Found %d sequence(s) but expected %d",
             seqno, msa->count);
    msa_destroy(msa);
    return NULL;
  }

  return msa;
}

msa_t ** phylip_parse_multisequential(phylip_t * fd, long * count)
{
  long msa_slotalloc = 10;
  long msa_maxcount = 0;
  char * p;
  
  *count = 0;

  msa_t ** msa = (msa_t **)xmalloc(msa_maxcount*sizeof(msa_t *));
  
  while (1)
  {
    if (*count == msa_maxcount)
    {
      msa_maxcount += msa_slotalloc;
      msa_t ** temp = (msa_t **)xmalloc(msa_maxcount*sizeof(msa_t *));
      memcpy(temp,msa,*count * sizeof(msa_t *));
      free(msa);
      msa = temp;
    }

    msa[*count] = phylip_parse_sequential(fd);
    if (msa[*count] == NULL)
      fatal("%s",bpp_errmsg);

    *count = *count + 1;

    /* if 'nloci' option was specified, break when the respective number of loci
       was read */
    if (*count == opt_locus_count) break;

    /* skip empty lines */

    while (1)
    {
      /* get next line */
      p = getnextline(fd);

      /* if no more lines break */
      if (!p) break;

      /* skip whitespace before sequence header */
      while (*p && whitespace(*p)) ++p;

      /* restart loop if blank line, otherwise break */
      if (*p) break;
    }

    if (!p) break;
  }

  return msa;
}
