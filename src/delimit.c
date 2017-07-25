/*
    Copyright (C) 2016-2017 Tomas Flouri and Ziheng Yang

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

static char * dmodel = NULL;

static snode_t ** trav = NULL;
static int trav_size = 0;

static long delimitations_count_recursive(snode_t * node)
{
  long x,y;

  if (!node->left)
    return 1;

  x = delimitations_count_recursive(node->left);
  y = delimitations_count_recursive(node->right);

  return x*y+1;
}

long delimitations_count(stree_t * stree)
{
  return delimitations_count_recursive(stree->root);
}

static void print()
{
  int i;

  printf("  ");
  for (i=0; i<trav_size; ++i)
    printf("%c", trav[i]->mark ? '1' : '0');
  printf("\n");
}

static void explore(snode_t ** start, snode_t ** end)
{
  while (end != start)
  {
    snode_t * curnode = *end;

    if (curnode->parent->mark)
    {
      curnode->mark = 1;
      print();

      explore(end, trav+trav_size-1);

      curnode->mark = 0;
    }

    end--;
  }
}

static void preorder_recursive(snode_t * node,
                               unsigned int * index,
                               snode_t ** buffer)
{
  if (!node->left)
    return;

  buffer[*index] = node;

  *index = *index + 1;

  preorder_recursive(node->left,  index, buffer);
  preorder_recursive(node->right, index, buffer);
}

void delimitations_enumerate(stree_t * stree)
{
  unsigned int index = 0;
  unsigned int i;

  trav = (snode_t **)xmalloc(stree->inner_count * sizeof(snode_t *));

  preorder_recursive(stree->root, &index, trav);
  assert(index == stree->inner_count);
  trav_size = index;

  dmodel = (char *)xmalloc((stree->inner_count+1)*sizeof(char));
  dmodel[stree->inner_count] = 0;
  for (i = 0; i < stree->inner_count; ++i)
  {
    trav[i]->mark = 0;
    dmodel[i] = '0';
  }

  /* print current delimitation model, i.e. all zeros 000..0 */
  print();

  /* print next model with only one population at the root, i.e. 1000...0 */
  stree->root->mark = 1;
  print();

  /* recursively explore all other possibilities */
  explore(trav,trav+trav_size-1);

  stree->root->mark = 0;

  free(trav);
}

void histories(stree_t * stree)
{
  
}

