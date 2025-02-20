/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static int list_insert(list_t * list, void * data, int where)
{
  if (!list) return 0;

  /* create list item */
  list_item_t * item = (list_item_t *)xmalloc(sizeof(list_item_t));
  item->data = data;

  /* if list is empty */
  if (!(list->count))
  {
    list->head = list->tail = item;
    list->count = 1;
    item->next = NULL;
    return 1;
  }

  /* append */
  if (where == DEF_LIST_APPEND)
  {
    list->tail->next = item;
    list->tail = item;
    item->next = NULL;
    list->count++;
    return 1;
  }

  /* prepend */
  item->next = list->head;
  list->head = item;
  list->count++;

  return 1;
}

void list_append(list_t * list, void * data)
{
  list_insert(list, data, DEF_LIST_APPEND);
}

void list_prepend(list_t * list, void * data)
{
  list_insert(list, data, DEF_LIST_PREPEND);
}

void list_clear(list_t * list, void (*cb_dealloc)(void *))
{
  list_item_t * head = list->head;

  while (head)
  {
    list_item_t * temp = head;
    head = head->next;
    if (cb_dealloc)
      cb_dealloc(temp->data);
    free(temp);
  }

  list->head = list->tail = NULL;
  list->count = 0;
}

long list_delitem(list_t * list, list_item_t * item, void (*cb_dealloc)(void *))
{
  list_item_t * li = list->head;
  list_item_t * prev = NULL;

  while (li)
  {
    if (li == item)
      break;
    prev = li;
    li = li->next;
  }

  /* if item not found then error */
  if (!li)
    return 0;

  if (prev)
  {
    prev->next = li->next;
  }
  else
  {
    /* if li is the first item set head */
    list->head = li->next;
  }

  if (cb_dealloc)
    cb_dealloc(item->data);
  free(item);

  list->count--;

  return 1;
}

/* move item to the end of the list */
long list_reposition_tail(list_t * list, list_item_t * item)
{
  /* check that item is indeed from list */
  list_item_t * li = list->head;
  list_item_t * prev = NULL;
  while (li)
  {
    if (li == item)
      break;
    prev = li;
    li = li->next;
  }

  /* if item not found then error */
  if (!li)
    return 0;


  /* if item was already at the end return */
  if (list->tail == li)
    return 1;

  /* from here on we know there is at least one more item after li */

  if (prev)
  {
    prev->next = li->next;
  }
  else
  {
    /* if li is the first item set head */
    list->head = li->next;
  }

  list->tail->next = li;
  li->next = NULL;
  list->tail = li;

  return 1;
}
