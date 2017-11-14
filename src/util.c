/*
    Copyright (C) 2016 Tomas Flouri and Ziheng Yang

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

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "bpp.h"

static const char * progress_prompt;
static unsigned long progress_next;
static unsigned long progress_size;
static unsigned long progress_chunk;
static const unsigned long progress_granularity = 200;

void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(1);
}

void progress_init(const char * prompt, unsigned long size)
{
  if (!opt_quiet)
  {
    progress_prompt = prompt;
    progress_size = size;
    progress_chunk = size < progress_granularity ?
      1 : size  / progress_granularity;
    progress_next = 0;
    fprintf(stderr, "%s %.0f%%", prompt, 0.0);
  }
}

void progress_update(unsigned int progress)
{
  if (!opt_quiet)
  {
    if (progress >= progress_next)
    {
      fprintf(stderr, "  \r%s %.0f%%", progress_prompt,
              100.0 * progress  / progress_size);
      progress_next = progress + progress_chunk;
    }
  }
}

void progress_done()
{
  if (!opt_quiet)
    fprintf(stderr, "  \r%s %.0f%%\n", progress_prompt, 100.0);
}

void * xmalloc(size_t size)
{
  void * t;
  t = malloc(size);
  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xcalloc(size_t nmemb, size_t size)
{
  void * t;
  t = calloc(nmemb,size);
  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

char * xstrchrnul(char *s, int c)
{
  char * r = strchr(s, c);

  if (r)
    return r;
  else
    return (char *)s + strlen(s);
}

char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)xmalloc(len+1);
  return strcpy(p,s);
}

char * xstrndup(const char * s, size_t len)
{
  char * p = (char *)xmalloc(len+1);
  strncpy(p,s,len);
  p[len] = 0;
  return p;
}

#if 0
long getusec(void)
{
  struct timeval tv;
  if(gettimeofday(&tv,0) != 0) return 0;
  return tv.tv_sec * 1000000 + tv.tv_usec;
}
#endif

FILE * xopen(const char * filename, const char * mode)
{
  FILE * out = fopen(filename, mode);
  if (!out)
    fatal("Cannot open file %s", filename);

  return out;
}

void * pll_aligned_alloc(size_t size, size_t alignment)
{
  void * mem;

#if (defined(_WIN32) || defined(_WIN64))
  mem = _aligned_malloc(size, alignment);
#else
  if (posix_memalign(&mem, alignment, size))
    mem = NULL;
#endif

  return mem;
}

void pll_aligned_free(void * ptr)
{
#if (defined(_WIN32) || defined(_WIN64))
  _aligned_free(ptr);
#else
  free(ptr);
#endif
}

#ifdef _MSC_VER
static int vasprintf(char **strp, const char *fmt, va_list ap)
{
  int len = _vscprintf(fmt, ap);
  if (len == -1) return -1;

  size_t size = (size_t)len+1;

  char * str = (char *)malloc(size);
  if (!str) return -1;

  int r = vsprintf_s(str, len + 1, fmt, ap);
  if (r == -1)
  {
    free(str);
    return -1;
  }

  *strp = str;
  return r;
}

int xasprintf(char ** strp, const char * fmt, ...)
{
  va_list ap;
  va_start(ap,fmt);
  int r = vasprintf(strp,fmt,ap);
  va_end(ap);
  return r;
}
#endif
