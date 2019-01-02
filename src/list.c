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
