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


/* Daniel J. Bernstein 2a hash function */
unsigned long hash_djb2a(char * s)
{
  unsigned long hash = 5381;
  unsigned long c;

  while ((c = (unsigned long)*s++))
    hash = ((hash << 5) + hash) ^ c;  /* hash*33 ^ c */

  return hash;
}

/* Fowler–Noll–Vo 1a hash function */
unsigned long hash_fnv(char * s)
{
  unsigned long hash = 14695981039346656037UL;
  unsigned long c;

  while ((c = (unsigned long)*s++))
  {
    hash ^= c;
    hash *= 1099511628211UL;
  }

  return hash;
}

static ht_item_t * hashitem_create(unsigned long key, void * value)
{
  ht_item_t * hi = (ht_item_t *)xmalloc(sizeof(ht_item_t));
  
  hi->key   = key;
  hi->value = value;

  return hi;
}

int hashtable_strcmp(void * x, void * y)
{
  return !strcmp((char *)x, (char *)y);
}

int hashtable_ptrcmp(void * x, void * y)
{
  return (x == y);
}

void * hashtable_find(hashtable_t * ht,
                      void * x,
                      unsigned long hash,
                      int (*cb_cmp)(void *, void *))
{
  unsigned long index = hash & (ht->table_size-1);
  list_item_t * li = (list_item_t *)(ht->entries[index]->head);

  while (li)
  {
    ht_item_t * hi = (ht_item_t *)(li->data);

    if ((hash == hi->key) && cb_cmp(hi->value, x))
      return hi->value;
    
    li = li->next; 
  }

  return NULL;
}
               

hashtable_t * hashtable_create(unsigned long items_count)
{
  unsigned long i;
  unsigned long size = 1;

  if (!items_count) return NULL;

  /* check that items_count has not the most-significant bit set */
  assert(!(items_count & (1ul << (sizeof(unsigned long)*CHAR_BIT - 1))));

  /* compute a size of greater than items count that is a multiple of 2 */
  items_count <<= 1;
  while (size < items_count)
    size <<= 1;

  /* allocate and init hash table */
  hashtable_t * ht = (hashtable_t *)xmalloc(sizeof(hashtable_t));
  ht->table_size = size;
  ht->entries_count = 0;

  /* allocate and init entries array */
  ht->entries = (list_t **)xmalloc(size*sizeof(list_t *));
  for (i = 0; i < size; ++i)
  {
    ht->entries[i] = (list_t *)xmalloc(sizeof(list_t));
    memset(ht->entries[i], 0, sizeof(list_t));
  }

  return ht;
}

int hashtable_insert(hashtable_t * ht,
                     void * x,
                     unsigned long hash,
                     int (*cb_cmp)(void *, void *))
{
  /* size is always a multiple of 2 and greater than 2 */
  unsigned long index = hash & (ht->table_size-1);

  list_t * list = ht->entries[index];


  if (hashtable_find(ht, x, hash, cb_cmp))
    return 0;

  ht_item_t * item = hashitem_create(hash,x);
  list_append(list, item);

  ht->entries_count++;

  return 1;
}

void hashtable_destroy(hashtable_t * ht, void (*cb_dealloc)(void *))
{
  unsigned long i;

  if (cb_dealloc)
  {
    for (i = 0; i < ht->table_size; ++i)
    {
      list_t * list = ht->entries[i];
      
      list_item_t * head = list->head;
      while (head)
      {
        ht_item_t * hi = (ht_item_t *)(head->data);
        cb_dealloc(hi->value);
        head = head->next;
      }
    }
  }

  for (i = 0; i < ht->table_size; ++i)
  {
    list_clear(ht->entries[i], free);
    free(ht->entries[i]);
  }
  free(ht->entries);
  free(ht); 
}

int cb_cmp_pairlabel(void * a, void * b)
{
  pair_t * pair = (pair_t *)a;
  char * label = (char *)b;

  return (!strcmp(pair->label,label));
}

void hashtable_insert_force(hashtable_t * ht,
                            void * x,
                            unsigned long hash)
{
  /* size is always a multiple of 2 and greater than 2 */
  unsigned long index = hash & (ht->table_size-1);

  list_t * list = ht->entries[index];

  ht_item_t * item = hashitem_create(hash,x);
  list_append(list, item);

  ht->entries_count++;
}
