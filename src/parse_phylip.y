/*
    Copyright (C) 2015 Tomas Flouri

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

extern int phylip_lex();
extern FILE * phylip_in;
extern void phylip_lex_destroy();
extern int phylip_lineno;

static int sc = 0;

typedef struct list_s
{
  locus_t * locus;
  struct list_s * next;
} list_t;

struct sequence_s
{
  int len;
  char * label;
  char * sequence;
};

void phylip_error(list_t ** list, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

list_t * yy_create_list()
{
  list_t * list = xrealloc(0, sizeof(list_t));
  memset(list,0,sizeof(list_t));
  return list;
}

void yy_dealloc_list(list_t * list)
{
  if (!list) return;

  int i;

  for (i = 0; i < list->locus->seq_count; ++i)
  {
    free(list->locus->seq[i]);
    free(list->locus->label[i]);
  }
  free(list->locus->seq);
  free(list->locus->label);
  free(list->locus);

  if (list->next);
    yy_dealloc_list(list->next);
  free(list);
}

static locus_t ** linearize(list_t * list, int * loci_count)
{
  int i;
  list_t * temp;
  
  *loci_count = 0;

  /* get number of loci */
  for (temp = list; temp; temp = temp->next)
    (*loci_count)++;

  /* allocate array for loci */
  locus_t ** loci = (locus_t **)xmalloc((*loci_count)*sizeof(locus_t *));

  /* copy loci structures and deallocate list */
  for (i = 0; list; i++)
  {
    loci[i] = list->locus;
    temp = list;
    list = list->next;
    free(temp);
  }
  return loci;
}

%}

%union
{
  char * s;
  char * d;
  struct list_s * list;
  struct locus_s * locus;
  struct sequence_s * sequence;
}

%error-verbose
%parse-param {struct list_s ** list}
%destructor { yy_dealloc_list($$); } loci_list

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<d> number 
%type<list> loci_list
%type<locus> locus
%type<sequence> seq_list
%start input
%%

input: loci_list 
{
  *list = $1;
}

loci_list: locus loci_list
{
  $$ = yy_create_list();
  
  $$->locus = $1;
  $$->next = $2;
}
        | locus
{
  $$ = yy_create_list();

  $$->locus = $1;
  $$->next = NULL;
};

locus: number number seq_list
{
  int i;

  if (atoi($1) != sc)
  {
    fprintf(stderr,
            "Number of sequences read is %d but expected %d\n",
            sc,
            atoi($1));
    exit(1);
  }
  if (atoi($2) != ($3+sc-1)->len)
  {
    fprintf(stderr,
            "Sequence length is %d but expected %d\n",
            ($3+sc-1)->len,
            atoi($2));
    exit(1);
  }

  $$ = (locus_t *)calloc(1,sizeof(locus_t));

  $$->seq_count = atoi($1);
  $$->seq_len = atoi($2);
  $$->seq = (char **)malloc(sc*sizeof(char *));
  $$->label = (char **)malloc(sc*sizeof(char *));

  for (i = 0; i < sc; ++i)
  {
    $$->seq[i] = ($3+sc-i-1)->sequence;
    $$->label[i] = ($3+sc-i-1)->label;
  }

  free($1);
  free($2);
  free($3);
  
  sc = 0;
}

seq_list: STRING STRING seq_list
{
  $$ = (struct sequence_s *)realloc((void *)$3,
                                     (sc+1)*sizeof(struct sequence_s));
  ($$+sc)->label = $1;
  ($$+sc)->sequence = $2;
  ($$+sc)->len = strlen($2);

  if (($$+sc)->len != ($$+sc-1)->len)
  {
    fprintf(stderr, "Sequence length does not match\n");
    exit(0);
  }
  sc++;

}
        | STRING STRING
{
  $$ = (struct sequence_s *)calloc(1,sizeof(struct sequence_s));
  ($$+sc)->label = $1;
  ($$+sc)->sequence = $2;
  ($$+sc)->len = strlen($2);
  sc++;
}

number: NUMBER { $$=$1; };

%%

locus_t ** yy_parse_phylip(const char * filename, int * loci_count)
{
  struct list_s * list; 
  locus_t ** loci;

  *loci_count = 0;

  phylip_lineno = 1;
  phylip_in = fopen(filename, "r");
  if (!phylip_in)
  {
    fatal("Cannot open file %s", filename);
  }
  else if (phylip_parse(&list))
  {
    fatal("Cannot parse tree file %s (maybe non-binary?)", filename);
  }
  
  if (phylip_in) fclose(phylip_in);

  phylip_lex_destroy();

  loci = linearize(list,loci_count);

  return loci;
}
