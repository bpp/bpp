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
  const size_t alignment = 16;
  void * t = NULL;
  if (posix_memalign(& t, alignment, size) == -1)
    fatal("Unable to allocate enough memory.");

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

long getusec(void)
{
  struct timeval tv;
  if(gettimeofday(&tv,0) != 0) return 0;
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

void show_rusage()
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);

  fprintf(stderr, "Time: %.3fs (user)", r_usage.ru_utime.tv_sec * 1.0 + (double) r_usage.ru_utime.tv_usec * 1.0e-6);
  fprintf(stderr, " %.3fs (sys)", r_usage.ru_stime.tv_sec * 1.0 + r_usage.ru_stime.tv_usec * 1.0e-6);

#if defined __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  fprintf(stderr, " Memory: %.0fMB\n", r_usage.ru_maxrss * 1.0e-6);
#else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  fprintf(stderr, " Memory: %.0fMB\n", r_usage.ru_maxrss * 1.0e-3);
#endif
}

FILE * xopen(const char * filename, const char * mode)
{
  FILE * out = fopen(filename, mode);
  if (!out)
    fatal("Cannot open file %s", opt_outfile);

  return out;
}
