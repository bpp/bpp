/*
    Copyright (C) 2016-2018 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static int cb_cmp_pairpair(void * a, void * b)
{
  pair_t * p1 = (pair_t *)a;
  pair_t * p2 = (pair_t *)b;

  return (!strcmp(p1->label,p2->label));
}

static int cb_cmp_nodelabel(void * a, void * b)
{
  snode_t * node = (snode_t *)a;
  char * label = (char * )b;

  return (!strcmp(node->label,label));
}

hashtable_t * maplist_hash(list_t * maplist, hashtable_t * sht)
{
  if (!maplist) return NULL;

  list_item_t * li = maplist->head;

  hashtable_t * ht = hashtable_create(maplist->count);

  while (li)
  {
    mapping_t * mapping = (mapping_t *)(li->data);

    snode_t * node = hashtable_find(sht,
                                    (void *)(mapping->species),
                                    hash_fnv(mapping->species),
                                    cb_cmp_nodelabel);
    if (!node)
      fatal("Cannot find node with population label %s", mapping->species);

    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));

    pair->label = xstrdup(mapping->individual);
    pair->data = (void *)node;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(pair->label),
                          cb_cmp_pairpair))
      fatal("Duplicate mapping (%s -> %s) found in file %s",
            mapping->individual,
            mapping->species,
            opt_mapfile);

#if 0
    printf("Mapped %s -> %s\n", pair->label, node->label);
#endif
    li = li->next;
  }
  
  return ht;
}

void maplist_print(list_t * maplist)
{
  if (!maplist) return;

  list_item_t * li = maplist->head;
  while (li)
  {
    mapping_t * map = (mapping_t *)(li->data);
    printf("%s -> %s\n", map->individual, map->species);
    li = li->next;
  }
}

/* callback function for list deallocator */
void map_dealloc(void * data)
{
  if (data)
  {
    mapping_t * map = (mapping_t *)data;

    free(map->individual);
    free(map->species);
    free(map);
  }
}
