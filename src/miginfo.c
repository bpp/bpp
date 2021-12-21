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

static size_t alloc_initial = 128;
static size_t alloc_extend  = 128;

miginfo_t * miginfo_create(size_t alloc_size)
{
  miginfo_t * mi;
  
  mi = (miginfo_t *)xmalloc(sizeof(miginfo_t));

  mi->alloc_size = alloc_size;
  mi->count = 0;

  mi->time     = (double *)xmalloc((size_t)alloc_size * sizeof(double));
  mi->old_time = (double *)xmalloc((size_t)alloc_size * sizeof(double));
  mi->source   = (snode_t **)xmalloc((size_t)alloc_size * sizeof(snode_t *));
  mi->target   = (snode_t **)xmalloc((size_t)alloc_size * sizeof(snode_t *));

  return mi;
}

miginfo_t * miginfo_create_default()
{
  return miginfo_create(alloc_initial);
}

void miginfo_destroy(miginfo_t * mi)
{
  if (!mi) return;

  if (mi->time)     free(mi->time);
  if (mi->old_time) free(mi->old_time);
  if (mi->source)   free(mi->source);
  if (mi->target)   free(mi->target);
  free(mi);
}

miginfo_t * miginfo_realloc(miginfo_t * mi)
{
  miginfo_t * extmi;

  if (!mi)
    return miginfo_create_default();

  assert(mi->count+alloc_extend > 0);
  assert(alloc_extend > 0);
  extmi = miginfo_create(mi->count+alloc_extend);

  /* copy data */
  extmi->alloc_size = mi->count+alloc_extend;
  extmi->count = mi->count;
  memcpy(extmi->time,     mi->time,     mi->count*sizeof(double));
  memcpy(extmi->old_time, mi->old_time, mi->count*sizeof(double));
  memcpy(extmi->source,   mi->source,   mi->count*sizeof(snode_t *));
  memcpy(extmi->target,   mi->target,   mi->count*sizeof(snode_t *));

  miginfo_destroy(mi);
  return extmi;
}

void miginfo_append(miginfo_t ** miptr, snode_t * s, snode_t * t, double time)
{
  miginfo_t * mi;

  if (!miptr)
    fatal("miginfo_append(): the unthinkable has happened");
  if (!(*miptr))
    *miptr = miginfo_create_default();

  mi = *miptr;


  if (mi->count == mi->alloc_size)
    *miptr = mi = miginfo_realloc(mi);

  mi->time[mi->count]     = time;
  mi->old_time[mi->count] = 0;
  mi->source[mi->count]   = s;
  mi->target[mi->count]   = t;

  mi->count++;
}

void miginfo_clean(miginfo_t * mi)
{
  if (!mi) return;

  mi->count = 0;
}

void miginfo_clone(miginfo_t * mi, miginfo_t ** clonebuf)
{
  miginfo_t * minew = *clonebuf;

  if (minew && minew->alloc_size < mi->count)
  {
    miginfo_destroy(minew);
    minew = NULL;
  }
  if (!minew)
    minew = miginfo_create(mi->alloc_size);

  /* copy data */
  minew->count = mi->count;
  memcpy(minew->time,mi->time,mi->count * sizeof(double));
  memcpy(minew->old_time,mi->old_time,mi->count * sizeof(double));
  memcpy(minew->source,mi->source,mi->count * sizeof(snode_t *));
  memcpy(minew->target,mi->target,mi->count * sizeof(snode_t *));

  *clonebuf= minew;
}
