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
%{
#include "bpp.h"

extern int map_lex();
extern FILE * map_in;
extern void map_lex_destroy();
extern int map_lineno;

void map_error(list_t ** list, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

static list_t * yy_create_list()
{
  list_t * list = xrealloc(0, sizeof(list_t));
  memset(list,0,sizeof(list_t));
  return list;
}

static void yy_dealloc_list(list_t * list)
{
  if (!list) return;

  map_t * map = (map_t *)(list->data);

  free(map->individual);
  free(map->species);
  free(map);
  if (list->next);
    yy_dealloc_list(list->next);
  free(list);
}
%}

%union
{
  char * s;
  char * d;
  struct list_s * list;
  struct map_s * map;
}

%error-verbose
%parse-param {struct list_s ** list}
%destructor { yy_dealloc_list($$); } map_list

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> individual
%type<s> species
%type<list> map_list
%type<map> map
%start input
%%

input: map_list 
{
  *list = $1;
}

map_list: map map_list
{
  $$ = yy_create_list();
  
  $$->data = (void *)$1;
  $$->next = $2;
}
        | map
{
  $$ = yy_create_list();

  $$->data = (void *)$1;
  $$->next = NULL;
};

map: individual species
{
  $$ = (map_t *)calloc(1,sizeof(map_t));
  $$->individual = $1;
  $$->species = $2;
  $$->lineno = map_lineno;
};

individual: STRING { $$ = $1; }
          | NUMBER { $$ = $1; };

species: STRING { $$ = $1; }
       | NUMBER { $$ = $1; };

%%

list_t * yy_parse_map(const char * filename)
{
  struct list_s * list; 

  map_lineno = 1;
  map_in = fopen(filename, "r");
  if (!map_in)
  {
    fatal("Cannot open file %s", filename);
  }
  else if (map_parse(&list))
  {
    fatal("Error while parsing file %s", filename);
  }
  
  if (map_in) fclose(map_in);

  map_lex_destroy();

  return list;
}
