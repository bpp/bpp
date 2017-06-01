# Copyright (C) 2016 Tomas Flouri and Ziheng Yang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
# Heidelberg Institute for Theoretical Studies,
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# Profiling
#PROFILING=-g -pg
PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare

AVX2DEF=-DHAVE_AVX2
AVX2OBJ=core_partials_avx2.o
ifdef DISABLE_AVX2
  AVX2DEF=
  AVX2OBJ=
endif

ifndef CC
CC = gcc
endif
CFLAGS = -D_GNU_SOURCE -DHAVE_SSE3 -DHAVE_AVX $(AVX2DEF) -g -msse3 -O3 $(WARN) #-DDEBUG_GTREE_SIMULATE -DDEBUG_STREE_INIT
LINKFLAGS=$(PROFILING)
LIBS=-lm

BISON = bison
FLEX = flex

PROG=bpp

all: $(PROG)

OBJS=bpp.o rtree.o util.o parse_stree.o lex_stree.o arch.o phylip.o msa.o \
     maps.o locus.o parse_map.o lex_map.o mapping.o compress.o hash.o list.o \
     stree.o random.o gtree.o core_partials.o core_pmatrix.o core_likelihood.o \
     output.o method_00.o core_partials_avx.o core_partials_sse.o $(AVX2OBJ)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $+ $(LIBS)

%_avx.o: %_avx.c
	$(CC) $(CFLAGS) -c -mavx -o $@ $<

%_avx2.o: %_avx2.c
	$(CC) $(CFLAGS) -c -mavx2 -mfma -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

parse_%.c: parse_%.y
	$(BISON) -p $*_ -d -o $@ $<

lex_%.c: lex_%.l
	$(FLEX) -P $*_ -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
