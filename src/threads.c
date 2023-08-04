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

static thread_info_t * ti = NULL;

void threads_init()
{
  if (opt_threads > opt_locus_count)
    fprintf(stdout, "Warning: The number of threads shouldn't be greater than the number of loci");
}

void threads_pin_master() {}

long * threads_load_balance(msa_t ** msa_list)
{
  ti = (thread_info_t *)xmalloc((size_t)opt_threads * sizeof(thread_info_t));
  return NULL;
}

thread_info_t * threads_ti()
{
  return ti;
}

/* this is only used when resuming from a checkpoint */
void threads_set_ti(thread_info_t * tip)
{
  ti = tip;
}

void threads_exit()
{
  free(ti);
}
