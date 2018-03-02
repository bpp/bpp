/*
    Copyright (C) 2016-2018 Tomas Flouri and Ziheng Yang

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
#include "parse_phylip.h"
#include "bpp.h"

static size_t string_length = 0;

static char * append(size_t * current_len, const char * s2, size_t append_len)
{
  printf("%ld %s %ld\n", *current_len, s2, append_len);
  phylip_lval.s = (char *)xrealloc(phylip_lval.s, sizeof(char) * (*current_len+append_len+1));
  strncpy(&phylip_lval.s[*current_len], s2, append_len);
  phylip_lval.s[*current_len + append_len]=0;
  *current_len += append_len;
  return phylip_lval.s;
}

%}
%option noyywrap
%option prefix="phylip_"
%option nounput
%option noinput
%option yylineno
%s apos
%s quot

%%
<quot>{
\\\"                                         append(&string_length, "\"", 1);
\'                                           append(&string_length, "\'", 1);
\"                                           {BEGIN(INITIAL);return STRING;}
}

<apos>{
\\\'                                         append(&string_length, "\'", 1);
\"                                           append(&string_length, "\"", 1);
\'                                           {BEGIN(INITIAL);return STRING;}
}

<apos,quot>{
\\n                                          append(&string_length, "\n", 1);
\\t                                          append(&string_length, "\t", 1);
\\\\                                         append(&string_length, "\\", 1);
([^\"\'\\]|\n)+                              append(&string_length, phylip_text, phylip_leng);
}

\:                                           return COLON;
\;                                           return SEMICOLON;
\)                                           return CPAR;
\(                                           return OPAR;
\,                                           return COMMA;
\"                                           {string_length =0; phylip_lval.s = NULL; ;BEGIN(quot);}
\'                                           {string_length =0; phylip_lval.s = NULL; ;BEGIN(apos);}
[\+|\-]?[0-9]+                               {phylip_lval.d = xstrndup(phylip_text, phylip_leng); return NUMBER;}
[\+|\-]?[0-9]+\.[0-9]+([e|E][\+|\-]?[0-9]+)? {phylip_lval.d = xstrndup(phylip_text, phylip_leng); return NUMBER;}
([a-zA-Z_0-9]|(\.|\-|\^|\?))+                {phylip_lval.s = xstrndup(phylip_text, phylip_leng); return STRING;}
[ \t\n\r]*\#[^\n]*                           ;
[ \t\n\r]                                    ;
.                                            {fprintf(stderr, "Syntax error (%c)\n", phylip_text[0]); exit(EXIT_FAILURE);}
%%
