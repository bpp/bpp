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

static size_t alloc_initial = 128;
static size_t alloc_extend  = 128;

miginfo_t * miginfo_create(size_t alloc_size, int alloc_dlist)
{
  long i;
  miginfo_t * mi;
  assert(alloc_dlist == MI_DLI_ALLOC || alloc_dlist == MI_DLI_NOOP);
  
  mi = (miginfo_t *)xmalloc(sizeof(miginfo_t));

  mi->alloc_size = alloc_size;
  mi->count = 0;

  mi->me = (migevent_t *)xmalloc((size_t)alloc_size * sizeof(migevent_t));

  if (alloc_dlist == MI_DLI_ALLOC)
  {
    for (i = 0; i < alloc_size; ++i)
    {
      mi->me[i].di_src = (dlist_item_t *)xcalloc(1,sizeof(dlist_item_t));
      mi->me[i].di_tgt = (dlist_item_t *)xcalloc(1,sizeof(dlist_item_t));

      mi->me[i].di_src->data = (void *)(mi->me+i);
      mi->me[i].di_tgt->data = (void *)(mi->me+i);
    }
  }

  return mi;
}

miginfo_t * miginfo_create_default()
{
  return miginfo_create(alloc_initial, MI_DLI_ALLOC);
}

void miginfo_destroy(miginfo_t * mi, long msa_index, int dealloc_dlist)
{
  long i;
  if (!mi) return;
  assert(dealloc_dlist == MI_DLI_FREE || dealloc_dlist == MI_DLI_NOOP);

  if (mi->me)
  {
    if (dealloc_dlist == MI_DLI_FREE)
    {
      for (i = 0; i < mi->count; ++i)
        migevent_unlink(mi->me+i,msa_index);

      for (i = 0; i < mi->alloc_size; ++i)
      {
        free(mi->me[i].di_src);
        free(mi->me[i].di_tgt);
      }
    }
    free(mi->me);
  }
  free(mi);
}

/* extends the size of miginfo by "units" chunks of alloc_extend */
miginfo_t * miginfo_extend(miginfo_t * mi, size_t units)
{
  long i;
  miginfo_t * extmi;

  if (!mi)
    return miginfo_create_default();

  assert(mi->count+units*alloc_extend > 0);
  assert(units*alloc_extend > 0);

  /* create a new structure but do not allocate dlist items as we will copy
     them and manually allocate new ones for the new elements */
  extmi = miginfo_create(mi->alloc_size+units*alloc_extend, MI_DLI_NOOP);

  /* copy data */
  extmi->alloc_size = mi->alloc_size+units*alloc_extend;
  extmi->count = mi->count;

  memcpy(extmi->me, mi->me, mi->alloc_size*sizeof(migevent_t));

  for (i = mi->alloc_size; i < mi->alloc_size+units*alloc_extend; ++i)
  {
    extmi->me[i].di_src = (dlist_item_t *)xcalloc(1,sizeof(dlist_item_t));
    extmi->me[i].di_tgt = (dlist_item_t *)xcalloc(1,sizeof(dlist_item_t));
  }

  /* update the pointers within the doubly-linked list items.
     Note: the below loop will also update the contents of the first mi->count 
     dlist_item_t entries which are already linked with populations, thus no
     need for linking/unlinking */
  for (i = 0; i < mi->alloc_size+units*alloc_extend; ++i)
  {
    extmi->me[i].di_src->data = (void *)(extmi->me+i);
    extmi->me[i].di_tgt->data = (void *)(extmi->me+i);
  }

  /* deallocate old miginfo *WITHOUT* deleting *NOR* unlinking dlist items */
  miginfo_destroy(mi,-1,MI_DLI_NOOP);

  return extmi;
}

/* check if *miptr structure can fit 'extra' migrations. If not, extend it */
void miginfo_check_and_extend(miginfo_t ** miptr, size_t extra)
{
  miginfo_t * mi;
  size_t ext_units = 0;

  if (!miptr)
    fatal("%s - the unthinkable has happened", __FUNCTION__);

  mi = *miptr;
  
  if (mi && (mi->alloc_size - mi->count) >= extra) return;

  if (!mi)
  {
    if (extra > alloc_initial)
    {
      size_t diff = extra - alloc_initial;
      ext_units = diff/alloc_extend + !!(diff % alloc_extend);
    }
    mi = miginfo_create(alloc_initial + ext_units*alloc_extend, MI_DLI_ALLOC);
  }
  else
  {
    size_t avail = mi->alloc_size - mi->count;
    assert(extra > avail);
    size_t diff = extra - avail;
    ext_units = diff/alloc_extend + !!(diff % alloc_extend);
    assert(ext_units >= 1);

    mi = miginfo_extend(mi, ext_units);
  }

  *miptr = mi;
}

/* append migration to *miptr, potentially extend it */
void miginfo_append(miginfo_t ** miptr, snode_t * s, snode_t * t, double time, long msa_index)
{
  miginfo_t * mi;

  if (!miptr)
    fatal("%s - the unthinkable has happened", __FUNCTION__);

  if (!(*miptr))
    *miptr = miginfo_create_default();

  mi = *miptr;

  if (mi->count == mi->alloc_size)
    *miptr = mi = miginfo_extend(mi,1); /* extend by one block */

  migevent_t * me = mi->me+mi->count;

  me->time     = time;
  me->old_time = 0;
  me->source   = s;
  me->target   = t;

  dlist_item_append(s->mig_source[msa_index],me->di_src);
  dlist_item_append(t->mig_target[msa_index],me->di_tgt);

  mi->count++;
}

/* Danger: the function corrupts the data structure:

   The migration is copied to a different branch and the
   dlist_item_t pointers are swapped, but at the same time
   it doesn't remove the migration from the original branch
   and assumes the caller will do it by decreasing the count. 

   Another assumption is that once we use miginfo_move to
   move an event, all subsequent events in the list must also
   be moved (i.e. we always move n events from the end).

   Copies migration event me to the *miptr structure */

void miginfo_move(migevent_t * me, miginfo_t ** miptr)
{
  miginfo_t * mi;

  if (!miptr)
    fatal("%s - the unthinkable has happened", __FUNCTION__);
  if (!(*miptr))
    *miptr = miginfo_create_default();

  mi = *miptr;

  if (mi->count == mi->alloc_size)
    *miptr = mi = miginfo_extend(mi,1); /* extend by one block */

  mi->me[mi->count].time = me->time;
  mi->me[mi->count].old_time = 0;
  mi->me[mi->count].source = me->source;
  mi->me[mi->count].target = me->target;

  /* TODO: Be careful about this part */
  SWAP(mi->me[mi->count].di_src, me->di_src);
  SWAP(mi->me[mi->count].di_tgt, me->di_tgt);

  /* update the actual pointers within the dlistitem */
  mi->me[mi->count].di_src->data = mi->me+mi->count;
  mi->me[mi->count].di_tgt->data = mi->me+mi->count;

  me->di_src->data = me;
  me->di_tgt->data = me;

  mi->count++;
}


void migevent_link(migevent_t * me, long msa_index)
{
  dlist_item_append(me->source->mig_source[msa_index],me->di_src);
  dlist_item_append(me->target->mig_target[msa_index],me->di_tgt);
}

void migevent_unlink(migevent_t * me, long msa_index)
{
  /* first re-link the event before the current event with the one after */
  if (me->di_src->prev)
    me->di_src->prev->next = me->di_src->next;
  else
    me->source->mig_source[msa_index]->head = me->di_src->next;

  /*  now re-link the event after the current with the one before */
  if (me->di_src->next)
    me->di_src->next->prev = me->di_src->prev;
  else
    me->source->mig_source[msa_index]->tail = me->di_src->prev;

  /* now do the same for the target */
  if (me->di_tgt->prev)
    me->di_tgt->prev->next = me->di_tgt->next;
  else
    me->target->mig_target[msa_index]->head = me->di_tgt->next;

  if (me->di_tgt->next)
    me->di_tgt->next->prev = me->di_tgt->prev;
  else
    me->target->mig_target[msa_index]->tail = me->di_tgt->prev;
}

/* "clean" mi by resetting count to 0 and unlinking events from populations */
void miginfo_clean(miginfo_t * mi, long msa_index)
{
  long i;

  if (!mi) return;

  /* unlink old migrations from populations */
  for (i = 0; i < mi->count; ++i)
    migevent_unlink(mi->me+i,msa_index); 

  mi->count = 0;
}

void miginfo_clone(miginfo_t * mi,
                   miginfo_t ** clonebuf,
                   stree_t * clone_stree,
                   long msa_index)
{
  long i;

  if (!*clonebuf)
  {
    *clonebuf = miginfo_create(mi->alloc_size, MI_DLI_ALLOC);
  }
  else
  {
    miginfo_clean(*clonebuf, msa_index);

    miginfo_check_and_extend(clonebuf, mi->count);
  }

  for (i = 0; i < mi->count; ++i)
    miginfo_append(clonebuf,
                   clone_stree->nodes[mi->me[i].source->node_index],
                   clone_stree->nodes[mi->me[i].target->node_index],
                   mi->me[i].time,
                   msa_index);
}
