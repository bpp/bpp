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

/* legacy random number generators */
static unsigned int z_rndu = 666;

double legacy_rndu (void)
{
/* 32-bit integer assumed.
   From Ripley (1987) p. 46 or table 2.4 line 2. 
   This may return 0 or 1, which can be a problem.
*/
   z_rndu = z_rndu*69069 + 1;
   if(z_rndu==0 || z_rndu==4294967295)  z_rndu = 13;
   return z_rndu/4294967295.0;
}
