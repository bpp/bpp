# Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang
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
# Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
# Department of Genetics, Evolution and Environment,
# University College London, Gower Street, London WC1E 6BT, England

# Profiling
#PROFILING=-g -pg -no-pie
#PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare

AVXDEF=-DHAVE_AVX
AVXOBJ=core_partials_avx.o core_likelihood_avx.o

AVX2DEF=-DHAVE_AVX2
AVX2OBJ=core_partials_avx2.o core_likelihood_avx2.o

ifdef DISABLE_AVX2
  AVX2DEF=
  AVX2OBJ=
endif

ifdef DISABLE_AVX
  AVXDEF=
  AVXOBJ=
endif

ifndef CC
CC = gcc-7
endif
CFLAGS = -D_GNU_SOURCE -DHAVE_SSE3 $(AVXDEF) $(AVX2DEF) -g -msse3 -O3 $(WARN) # -pg -no-pie -DDEBUG_GTREE_SIMULATE -DDEBUG_STREE_INIT
LINKFLAGS=$(PROFILING)
LIBS=-lm -lpthread

PROG=bpp

all: $(PROG)

OBJS=bpp.o rtree.o util.o \
     arch.o phylip.o msa.o \
     maps.o locus.o \
     mapping.o compress.o hash.o list.o \
     stree.o random.o gtree.o core_partials.o core_pmatrix.o core_likelihood.o \
     output.o core_partials_sse.o dlist.o allfixed.o core_likelihood_sse.o \
     prop_mixing.o method.o delimit.o prop_rj.o summary.o cfile.o hardware.o \
     revolutionary.o diploid.o dump.o load.o summary11.o simulate.o cfile_sim.o \
     gamma.o prop_gamma.o threads.o treeparse.o parsemap.o msci_gen.o visual.o \
     pdfgen.o constraint.o debug.o lswitch.o miginfo.o ming2.o a1b1.o \
     bfdriver.o $(AVXOBJ) $(AVX2OBJ)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $+ $(LIBS) $(LDFLAGS)

%_avx.o: %_avx.c
	$(CC) $(CFLAGS) -c -mavx -o $@ $<

%_avx2.o: %_avx2.c
	$(CC) $(CFLAGS) -c -mavx2 -mfma -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ $(OBJS) gmon.out $(PROG)
