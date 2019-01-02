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

#define DEF_LIST_APPEND   0
#define DEF_LIST_PREPEND  1

static dlist_item_t * dlist_insert(dlist_t * dlist, void * data, int where)
{
  if (!dlist) return 0;

  /* create list item */
  dlist_item_t * item = (dlist_item_t *)xmalloc(sizeof(dlist_item_t));
  item->data = data;

  /* if list is empty */
  if (!dlist->head)
  {
    dlist->head = dlist->tail = item;
    item->prev  = NULL;
    item->next  = NULL;
    return item;
  }

  /* append */
  if (where == DEF_LIST_APPEND)
  {
    dlist->tail->next = item;
    item->prev = dlist->tail;
    item->next = NULL;
    dlist->tail = item;
    return item;
  }

  /* prepend */
  item->next = dlist->head;
  item->prev = NULL;
  dlist->head->prev = item;
  dlist->head = item;

  return item;
}

dlist_item_t * dlist_append(dlist_t * dlist, void * data)
{
  return dlist_insert(dlist, data, DEF_LIST_APPEND);
}

dlist_item_t * dlist_prepend(dlist_t * dlist, void * data)
{
  return dlist_insert(dlist, data, DEF_LIST_PREPEND);
}

void dlist_clear(dlist_t * dlist, void (*cb_dealloc)(void *))
{
  dlist_item_t * head = dlist->head;

  while (head)
  {
    dlist_item_t * temp = head;
    head = head->next;
    if (cb_dealloc)
      cb_dealloc(temp->data);
    free(temp);
  }

  dlist->head  = dlist->tail = NULL;
}

void dlist_item_remove(dlist_item_t * item)
{
  if (item->prev)
    item->prev->next = item->next;

  if (item->next)
    item->next->prev = item->prev;

  item->prev = item->next = NULL;
}

void dlist_item_append(dlist_t * dlist, dlist_item_t * item)
{
  assert(dlist);
  
  /* if list is empty */
  if (!dlist->head)
  {
    dlist->head = dlist->tail = item;
    item->next = item->prev = NULL;
    return;
  }

  dlist->tail->next = item;
  item->prev = dlist->tail;
  item->next = NULL;
  dlist->tail = item;
}

void dlist_item_prepend(dlist_t * dlist, dlist_item_t * item)
{
  assert(dlist);

  /* if list is empty */
  if (!dlist->head)
  {
    dlist->head = dlist->tail = item;
    item->next = item->prev = NULL;
    return;
  }

  item->next = dlist->head;
  item->prev = NULL;
  dlist->head->prev = item;
  dlist->head = item;
}

dlist_t * dlist_create()
{
  return (dlist_t *)xcalloc(1,sizeof(dlist_t));
}

void dlist_destroy(dlist_t * dlist)
{
  free(dlist);
}
