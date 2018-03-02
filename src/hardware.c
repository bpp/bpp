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

#include "bpp.h"

/*
    Apple machines should always default to assembly code due to
    inconsistent versioning in LLVM/clang, see issue #138

    https://github.com/xflouris/libpll/issues/138

*/
#if (defined(__APPLE__) || defined(_MSC_VER)) || \
    (!defined(__clang__) && defined(__GNUC__) && (__GNUC__ < 4 || \
      (__GNUC__ == 4 && __GNUC_MINOR__ < 8))) || \
    (defined(__clang__) && (__clang_major__ < 3 || \
      (__clang_major__ == 3 && __clang_minor__ < 9)))
  
  #ifndef _MSC_VER
    #if defined(__i386__) && defined(__PIC__)
      #if (defined(__GNUC__) && __GNUC__ < 3)
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("xchgl\t%%ebx, %k1\n\t"                      \
           "cpuid\n\t"                                  \
           "xchgl\t%%ebx, %k1\n\t"                      \
           : "=a" (a), "=&r" (b), "=c" (c), "=d" (d)    \
           : "0" (level), "2" (count))
      #else
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("xchg{l}\t{%%}ebx, %k1\n\t"                  \
           "cpuid\n\t"                                  \
           "xchg{l}\t{%%}ebx, %k1\n\t"                  \
           : "=a" (a), "=&r" (b), "=c" (c), "=d" (d)    \
           : "0" (level), "2" (count))
      #endif
    #elif defined(__x86_64__) && (defined(__code_model_medium__) || \
          defined(__code_model_large__)) && defined(__PIC__)
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("xchg{q}\t{%%}rbx, %q1\n\t"                  \
           "cpuid\n\t"                                  \
           "xchg{q}\t{%%}rbx, %q1\n\t"                  \
           : "=a" (a), "=&r" (b), "=c" (c), "=d" (d)    \
           : "0" (level), "2" (count))
    #else
#define cpuid(level, count, a, b, c, d)                 \
  __asm__ ("cpuid\n\t"                                  \
           : "=a" (a), "=b" (b), "=c" (c), "=d" (d)     \
           : "0" (level), "2" (count))
    #endif
  #else
#define cpuid(level,count,a,b,c,d)                      \
    {                                                   \
      int cpui[4];                                      \
      __cpuidex(cpui,level,count);                      \
      a=cpui[0]; b=cpui[1]; c=cpui[2]; d=cpui[3];       \
    }
  #endif

void cpu_features_detect()
{
  unsigned int a,b,c,d;

  mmx_present = 0;
  sse_present = 0;
  sse2_present = 0;
  sse3_present = 0;
  ssse3_present = 0;
  sse41_present = 0;
  sse42_present = 0;
  popcnt_present = 0;
  avx_present = 0;
  avx2_present = 0;

#if defined(__PPC__)
  altivec_present = 1;
#else

  cpuid(0,0,a,b,c,d);
  unsigned int maxlevel = a & 0xff;

  if (maxlevel >= 1)
  {
    cpuid(1,0,a,b,c,d);
    mmx_present    = (d >> 23) & 1;
    sse_present    = (d >> 25) & 1;
    sse2_present   = (d >> 26) & 1;
    sse3_present   = (c >>  0) & 1;
    ssse3_present  = (c >>  9) & 1;
    sse41_present  = (c >> 19) & 1;
    sse42_present  = (c >> 20) & 1;
    popcnt_present = (c >> 23) & 1;
    avx_present    = (c >> 28) & 1;

    if (maxlevel >= 7)
    {
      cpuid(7,0,a,b,c,d);
      avx2_present = (b >> 5) & 1;
    }
  }
#endif
}

#else


void cpu_features_detect()
{
  mmx_present = 0;
  sse_present = 0;
  sse2_present = 0;
  sse3_present = 0;
  ssse3_present = 0;
  sse41_present = 0;
  sse42_present = 0;
  popcnt_present = 0;
  avx_present = 0;
  avx2_present = 0;

#if defined(__PPC__)
  altivec_present = __builtin_cpu_supports("altivec");
#elif defined(__x86_64__) || defined(__i386__)
  mmx_present     = __builtin_cpu_supports("mmx");
  sse_present     = __builtin_cpu_supports("sse");
  sse2_present    = __builtin_cpu_supports("sse2");
  sse3_present    = __builtin_cpu_supports("sse3");
  ssse3_present   = __builtin_cpu_supports("ssse3");
  sse41_present   = __builtin_cpu_supports("sse4.1");
  sse42_present   = __builtin_cpu_supports("sse4.2");
  popcnt_present  = __builtin_cpu_supports("popcnt");
  avx_present     = __builtin_cpu_supports("avx");
  avx2_present    = __builtin_cpu_supports("avx2");
#endif
}

#endif

void cpu_features_show()
{
  fprintf(stderr, "Detected CPU features:");
  if (altivec_present)
    fprintf(stderr, " altivec");
  if (mmx_present)
    fprintf(stderr, " mmx");
  if (sse_present)
    fprintf(stderr, " sse");
  if (sse2_present)
    fprintf(stderr, " sse2");
  if (sse3_present)
    fprintf(stderr, " sse3");
  if (ssse3_present)
    fprintf(stderr, " ssse3");
  if (sse41_present)
    fprintf(stderr, " sse4.1");
  if (sse42_present)
    fprintf(stderr, " sse4.2");
  if (popcnt_present)
    fprintf(stderr, " popcnt");
  if (avx_present)
    fprintf(stderr, " avx");
  if (avx2_present)
    fprintf(stderr, " avx2");
  fprintf(stderr, "\n");
}

void cpu_setarch()
{
  /* if arch specified by user, leave it be */
  if (opt_arch != -1)
  {
    if (opt_arch == PLL_ATTRIB_ARCH_CPU)
      printf("User specified SIMD ISA: CPU\n\n");
    else if (opt_arch == PLL_ATTRIB_ARCH_SSE)
      printf("User specified SIMD ISA: SSE\n\n");
    else if (opt_arch == PLL_ATTRIB_ARCH_AVX)
      printf("User specified SIMD ISA: AVX\n\n");
    else if (opt_arch == PLL_ATTRIB_ARCH_AVX2)
      printf("User specified SIMD ISA: AVX2\n\n");
    else
      fatal("Internal error when setting arch");

    return;
  }

  /* otherwise set best present SIMD */
  opt_arch = PLL_ATTRIB_ARCH_CPU;

  if (sse2_present)
    opt_arch = PLL_ATTRIB_ARCH_SSE;
#ifdef HAVE_AVX
  if (avx_present)
    opt_arch = PLL_ATTRIB_ARCH_AVX;
#endif
#ifdef HAVE_AVX2
  if (avx2_present)
    opt_arch = PLL_ATTRIB_ARCH_AVX2;
#endif

  if (opt_arch == PLL_ATTRIB_ARCH_CPU)
    printf("Auto-selected SIMD ISA: CPU\n\n");
  else if (opt_arch == PLL_ATTRIB_ARCH_SSE)
    printf("Auto-selected SIMD ISA: SSE\n\n");
  else if (opt_arch == PLL_ATTRIB_ARCH_AVX)
    printf("Auto-selected SIMD ISA: AVX\n\n");
  else if (opt_arch == PLL_ATTRIB_ARCH_AVX2)
    printf("Auto-selected SIMD ISA: AVX2\n\n");
  else
    fatal("Internal error when setting arch");
}

#ifdef _MSC_VER
int pll_ctz(unsigned int x)
{
  unsigned long index;

  int rc = _BitScanForward(&index,x);
  if (rc)
    return (int)index;

  return 32;

}

unsigned int pll_popcount(unsigned int x)
{
  if (popcnt_present)
    return __popcnt(x);

  /* non-vectorized way */
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

unsigned int pll_popcount64(unsigned long x)
{
  if (popcnt_present)
    return __popcnt64(x);

  x = (x & 0x5555555555555555ul) + ((x >> 1) & 0x5555555555555555ul);
  x = (x & 0x3333333333333333ul) + ((x >> 2) & 0x3333333333333333ul);
  x = (x & 0x0F0F0F0F0F0F0F0Ful) + ((x >> 4) & 0x0F0F0F0F0F0F0F0Ful);
  return (x * 0x0101010101010101ul) >> 56;
}

#endif
