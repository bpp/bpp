/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;

static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static long is_emptyline(const char * line)
{
  size_t ws = strspn(line, " \t\r\n");
  if (!line[ws] || line[ws] == '*' || line[ws] == '#') return 1;
  return 0;
}

/* get string from current position pointed by line until a delimiter. Create
   new string with result, store it in value, and return number of characters
   read */
static long get_delstring(const char * line, const char * del, char ** value)
{
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters until a delimiter is found */
  char * end = start + strcspn(start, del);

  *end = 0;

  if (start==end)
  {
    free(s);
    return 0;
  }

  *value = xstrdup(start);

  free(s);
  return ws + end - start;
}

static long parse_mapping(const char * line, char ** tag, char ** species)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_delstring(p," \t\r\n", tag);
  if (!count) goto l_unwind;

  p += count;

  count = get_delstring(p," \t\r\n", species);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;
  

l_unwind:
  free(s);
  return ret;
}

list_t * parse_date_mapfile(const char * mapfile)
{
  long line_count = 0;
  long ret = 1;
  FILE * fp;
  list_t * list;
  mappingDate_t * mapping = NULL;
  char * tag = NULL;
  char * date = NULL;

  fp = xopen(mapfile,"r");

  list = (list_t *)xcalloc(1,sizeof(list_t));

  regex_t regDecimal;
  int regCompiled, regMatch;
  char* ptr; /* Needed for the strtod function */
   /* Regex of decimal number */
   regCompiled = regcomp(&regDecimal, "^([0-9]*)((\\.)?([0-9]+))$" , REG_EXTENDED);
   if (regCompiled == 1) {
     fprintf(stderr, "Regular expression did not compile.\n");
     exit(1);
   }

  while (getnextline(fp))
  {
    ++line_count;

    if (is_emptyline(line)) continue;

    if (!parse_mapping(line,&tag,&date))
    {
      ret = 0;
      fprintf(stderr, "Invalid entry in %s (line %ld)\n", mapfile, line_count);
      goto l_unwind;
    }

    /* construct a mapping record */
    mapping = (mappingDate_t *)xmalloc(sizeof(mapping_t));
    mapping->individual = tag;


    /* Need to convert to double*/

     regMatch = regexec(&regDecimal, date, 0, NULL, 0);
      if (regMatch != 0) {
        fprintf(stderr, "Problem reading datefile. Regular expression does not match a decimal number.\n");
        exit(1);
      }

    mapping->date = strtod(date, &ptr);
    free(date);

    /* insert item into list */
    list_append(list,(void *)mapping);

    tag = date = NULL;
  }


l_unwind:
  fclose(fp);
  if (ret == 0)
  {
    list_clear(list,mapDate_dealloc);
    free(list);
    list = NULL;
  }

  regfree(&regDecimal);
  return list;
}

list_t * parse_mapfile(const char * mapfile)
{
  long line_count = 0;
  long ret = 1;
  FILE * fp;
  list_t * list;
  mapping_t * mapping = NULL;
  char * tag = NULL;
  char * species = NULL;

  fp = xopen(mapfile,"r");

  list = (list_t *)xcalloc(1,sizeof(list_t));

  while (getnextline(fp))
  {
    ++line_count;

    if (is_emptyline(line)) continue;

    if (!parse_mapping(line,&tag,&species))
    {
      ret = 0;
      fprintf(stderr, "Invalid entry in %s (line %ld)\n", mapfile, line_count);
      goto l_unwind;
    }

    /* construct a mapping record */
    mapping = (mapping_t *)xmalloc(sizeof(mapping_t));
    mapping->individual = tag;
    mapping->species = species;
    printf("%s %s\n", tag, species);

    /* insert item into list */
    list_append(list,(void *)mapping);

    tag = species = NULL;
  }

l_unwind:
  fclose(fp);
  if (ret == 0)
  {
    list_clear(list,map_dealloc);
    free(list);
    list = NULL;
  }
  return list;
}


static int cmp_tipDates(const void * a, const void * b){
  
  mappingDate_t * const * x = a;
  mappingDate_t * const * y = b;

  if ((*x)->date - (*y)->date > 0) return 1;
  return -1;
} 

mappingDate_t ** prepTipDatesInfer(stree_t * stree, list_t* dateList,
		 double mu_bar) {

	int i = 0;

	/* Converts dates to expected number of substitutions and creates an sorted 
	 * array of the dates */ 
        list_item_t* list = dateList->head;

	int tipDateArrayLen = dateList->count;
        mappingDate_t ** tipDateArray = xmalloc(tipDateArrayLen * sizeof(mappingDate_t *));
        list = dateList->head;

        while (list) {
		((mappingDate_t *) list->data)->date = ((mappingDate_t *) list->data )->date * mu_bar;
                tipDateArray[i] = (mappingDate_t *)list->data;

                i++;
                list = list->next;

        }

        qsort(tipDateArray, tipDateArrayLen, sizeof(mappingDate_t *), cmp_tipDates);


  return tipDateArray;
}

