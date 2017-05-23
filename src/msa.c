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

#include "bpp.h"

void msa_print(msa_t * msa)
{
  if (!msa) return;

  int i;

  printf("%d %d\n", msa->count, msa->length);
  for (i = 0; i < msa->count; ++i)
    printf("%s %s\n", msa->label[i], msa->sequence[i]);
}

void msa_destroy(msa_t * msa)
{
  int i;

  if (msa->label)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->label[i])
        free(msa->label[i]);
    free(msa->label);
  }

  if (msa->sequence)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->sequence[i])
        free(msa->sequence[i]);
    free(msa->sequence);
  }

  free(msa);
}
