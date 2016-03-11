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
#include "parse_rtree.h"
#include "bpp.h"

static size_t string_length = 0;

static char * append(size_t * dstlen, const char * src, size_t srclen)
{
  char * mem = (char *)realloc(rtree_lval.s, (*dstlen + srclen + 1));
  if (!mem)
  {
    fatal("Unable to allocate enough memory.");
  }
  rtree_lval.s =  mem;
  strncpy(&rtree_lval.s[*dstlen], src, srclen);
  rtree_lval.s[*dstlen + srclen] = 0;
  *dstlen += srclen;
  return rtree_lval.s;
}

%}
%option noyywrap
%option prefix="rtree_"
%option nounput
%option noinput
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
([^\"\'\\]|\n)+                              append(&string_length, rtree_text, rtree_leng);
}

\:                                           return COLON;
\;                                           return SEMICOLON;
\)                                           return CPAR;
\(                                           return OPAR;
\,                                           return COMMA;
\"                                           {string_length = 0; BEGIN(quot);}
\'                                           {string_length = 0; BEGIN(apos);}
[\+|\-]?[0-9]+                               {rtree_lval.d = xstrndup(rtree_text, rtree_leng); return NUMBER;}
[\+|\-]?[0-9]+\.[0-9]+([e|E][\+|\-]?[0-9]+)? {rtree_lval.d = xstrndup(rtree_text, rtree_leng); return NUMBER;}
[A-Za-z0-9\'\.\-\_]*                         {rtree_lval.s = xstrndup(rtree_text, rtree_leng); return STRING;}
[ \t\n\r]                                    ;
.                                            {fatal("Syntax error (%c)\n", rtree_text[0]);}
%%
