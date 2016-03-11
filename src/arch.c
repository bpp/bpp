/*
    Copyright (C) 2014-2015 Tomas Flouri, Torbjorn Rognes, Jeff Epler

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

unsigned long arch_get_memused()
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  
#if defined __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return (unsigned long)(r_usage.ru_maxrss);
#else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return (unsigned long)r_usage.ru_maxrss * 1024;
#endif
}

unsigned long arch_get_memtotal()
{
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  long phys_pages = sysconf(_SC_PHYS_PAGES);
  long pagesize = sysconf(_SC_PAGESIZE);

  if ((phys_pages == -1) || (pagesize == -1))
    fatal("Cannot determine amount of RAM");

  // sysconf(3) notes that pagesize * phys_pages can overflow, such as
  // when long is 32-bits and there's more than 4GB RAM.  Since vsearch
  // apparently targets LP64 systems like x86_64 linux, this will not
  // arise in practice on the intended platform.

  if (pagesize > LONG_MAX / phys_pages)
    return LONG_MAX;
  else
    return (unsigned long)pagesize * (unsigned long)phys_pages;

#elif defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  size_t length = sizeof(ram);
  if(-1 == sysctl(mib, 2, &ram, &length, NULL, 0))
    fatal("Cannot determine amount of RAM");
  return ram;

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}
