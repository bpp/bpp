# Copyright (C) 2015 Tomas Flouri
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

ifndef CC
CC = gcc
CC=i686-w64-mingw32-gcc
endif
#CC=i686-w64-mingw32-gcc
CFLAGS = -D_GNU_SOURCE -g -msse3 -O3 $(WARN)
LINKFLAGS=$(PROFILING)
LIBS=-lm

BISON = bison
FLEX = flex

PROG=bpp

all: $(PROG)

OBJS=bpp.o rtree.o util.o parse_rtree.o lex_rtree.o arch.o parse_phylip.o \
     lex_phylip.o alignment.o maps.o locus.o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $+

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

parse_%.c: parse_%.y
	$(BISON) -p $*_ -d -o $@ $<

lex_%.c: lex_%.l
	$(FLEX) -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
