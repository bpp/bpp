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

void maplist_print(list_t * maplist)
{
  if (!maplist) return;

  list_item_t * li = maplist->head;
  while (li)
  {
    map_t * map = (map_t *)(li->data);
    printf("%s -> %s\n", map->individual, map->species);
    li = li->next;
  }
}

/* callback function for list deallocator */
void map_dealloc(void * data)
{
  if (data)
  {
    map_t * map = (map_t *)data;

    free(map->individual);
    free(map->species);
    free(map);
  }
}
