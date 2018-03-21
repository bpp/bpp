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

%{
#include "parse_stree.h"
#include "bpp.h"

#define YY_USER_ACTION preaction();

static size_t string_length = 0;

int stree_colstart = 0;
int stree_colend = 0;

static void preaction()
{
  stree_colstart = stree_colend;
  stree_colend = stree_colend + stree_leng - 1;
}

static char * append(size_t * dstlen, const char * src, size_t srclen)
{
  char * mem = (char *)xmalloc((*dstlen + srclen + 1)*sizeof(char));
  memcpy(mem,stree_lval.s,*dstlen);
  strncpy(mem+(*dstlen),src,srclen);
  mem[*dstlen+srclen] = 0;
  if (*dstlen)
    free(stree_lval.s);
  stree_lval.s =  mem;
  *dstlen += srclen;
  return stree_lval.s;
}

%}
%option noyywrap
%option prefix="stree_"
%option nounput
%option noinput
%option yylineno
%x apos
%x quot

%%
<quot>{
\\\"              { append(&string_length, "\\\"", 2);                         }
\'                { append(&string_length, "\'", 1);                           }
\"                { BEGIN(INITIAL); return STRING;                             }
}

<apos>{
\\\'              { append(&string_length, "\\\'", 2);                         }
\"                { append(&string_length, "\"", 1);                           }
\'                { BEGIN(INITIAL); return STRING;                             }
}

<apos,quot>{
\\n               { append(&string_length, "\\n", 2);                          }
\\t               { append(&string_length, "\\t", 2);                          }
\\                { append(&string_length, "\\", 1);                           }
\\\\              { append(&string_length, "\\\\", 2);                         }
([^\"\'\\]|\n)+   { append(&string_length, stree_text, stree_leng);            }
}

\:                { return COLON;                                              }
\;                { return SEMICOLON;                                          }
\)                { return CPAR;                                               }
\(                { return OPAR;                                               }
\,                { return COMMA;                                              }
\"                { string_length = 0; BEGIN(quot);                            }
\'                { string_length = 0; BEGIN(apos);                            }
\n                { stree_colstart = stree_colend = 0;                         }
[-+]?[0-9]+       { stree_lval.d = xstrndup(stree_text,stree_leng);
                    return NUMBER;                                             }
[+-]?(([0-9]+[\.]?[0-9]*)|([0-9]*[\.]?[0-9]+))([eE][+-]?[0-9]+)? {
                    stree_lval.d = xstrndup(stree_text,stree_leng); 
                    return NUMBER;                                             }
[^ \'\",\(\):;\[\]\t\n\r][^ \t\n\r\)\(\[\]\,:;]* {
                    stree_lval.s = xstrndup(stree_text,stree_leng);
                    return STRING;                                             }
[ \t\r]           { ;                                                          }
.                 { fatal("Syntax error (%c)", stree_text[0]);                 }
%%
