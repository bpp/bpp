/*
    Copyright (C) 2016-2025 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

/* Self-contained SHA-1 implementation for checkpoint integrity verification.
   Based on the public domain implementation by Steve Reid <steve@edmweb.com>
   and adapted for BPP. */

#include "bpp.h"

typedef struct
{
  uint32_t state[5];
  uint32_t count[2];
  unsigned char buffer[64];
} sha1_ctx_t;

#define SHA1_ROL(value, bits) (((value) << (bits)) | ((value) >> (32-(bits))))

/* blk0() and blk() perform the initial expand. */
static uint32_t sha1_blk0(const unsigned char * buf, int i)
{
  return (uint32_t)buf[i*4+0] << 24 |
         (uint32_t)buf[i*4+1] << 16 |
         (uint32_t)buf[i*4+2] <<  8 |
         (uint32_t)buf[i*4+3];
}

#define SHA1_BLK(i) (w[i&15] = SHA1_ROL(w[(i+13)&15] ^ w[(i+8)&15] ^ \
                                          w[(i+2)&15] ^ w[i&15], 1))

/* SHA-1 round functions */
#define SHA1_R0(v,w_,x,y,z,i) z += ((w_&(x^y))^y)     + w[i]         + 0x5A827999 + SHA1_ROL(v,5); w_ = SHA1_ROL(w_,30);
#define SHA1_R1(v,w_,x,y,z,i) z += ((w_&(x^y))^y)     + SHA1_BLK(i)  + 0x5A827999 + SHA1_ROL(v,5); w_ = SHA1_ROL(w_,30);
#define SHA1_R2(v,w_,x,y,z,i) z += (w_^x^y)            + SHA1_BLK(i)  + 0x6ED9EBA1 + SHA1_ROL(v,5); w_ = SHA1_ROL(w_,30);
#define SHA1_R3(v,w_,x,y,z,i) z += (((w_|x)&y)|(w_&x)) + SHA1_BLK(i) + 0x8F1BBCDC + SHA1_ROL(v,5); w_ = SHA1_ROL(w_,30);
#define SHA1_R4(v,w_,x,y,z,i) z += (w_^x^y)            + SHA1_BLK(i)  + 0xCA62C1D6 + SHA1_ROL(v,5); w_ = SHA1_ROL(w_,30);

static void sha1_transform(uint32_t state[5], const unsigned char buf[64])
{
  uint32_t a, b, c, d, e;
  uint32_t w[16];
  int i;

  for (i = 0; i < 16; i++)
    w[i] = sha1_blk0(buf, i);

  a = state[0];
  b = state[1];
  c = state[2];
  d = state[3];
  e = state[4];

  /* 4 rounds of 20 operations each */
  SHA1_R0(a,b,c,d,e, 0); SHA1_R0(e,a,b,c,d, 1);
  SHA1_R0(d,e,a,b,c, 2); SHA1_R0(c,d,e,a,b, 3);
  SHA1_R0(b,c,d,e,a, 4); SHA1_R0(a,b,c,d,e, 5);
  SHA1_R0(e,a,b,c,d, 6); SHA1_R0(d,e,a,b,c, 7);
  SHA1_R0(c,d,e,a,b, 8); SHA1_R0(b,c,d,e,a, 9);
  SHA1_R0(a,b,c,d,e,10); SHA1_R0(e,a,b,c,d,11);
  SHA1_R0(d,e,a,b,c,12); SHA1_R0(c,d,e,a,b,13);
  SHA1_R0(b,c,d,e,a,14); SHA1_R0(a,b,c,d,e,15);

  SHA1_R1(e,a,b,c,d,16); SHA1_R1(d,e,a,b,c,17);
  SHA1_R1(c,d,e,a,b,18); SHA1_R1(b,c,d,e,a,19);

  SHA1_R2(a,b,c,d,e,20); SHA1_R2(e,a,b,c,d,21);
  SHA1_R2(d,e,a,b,c,22); SHA1_R2(c,d,e,a,b,23);
  SHA1_R2(b,c,d,e,a,24); SHA1_R2(a,b,c,d,e,25);
  SHA1_R2(e,a,b,c,d,26); SHA1_R2(d,e,a,b,c,27);
  SHA1_R2(c,d,e,a,b,28); SHA1_R2(b,c,d,e,a,29);
  SHA1_R2(a,b,c,d,e,30); SHA1_R2(e,a,b,c,d,31);
  SHA1_R2(d,e,a,b,c,32); SHA1_R2(c,d,e,a,b,33);
  SHA1_R2(b,c,d,e,a,34); SHA1_R2(a,b,c,d,e,35);
  SHA1_R2(e,a,b,c,d,36); SHA1_R2(d,e,a,b,c,37);
  SHA1_R2(c,d,e,a,b,38); SHA1_R2(b,c,d,e,a,39);

  SHA1_R3(a,b,c,d,e,40); SHA1_R3(e,a,b,c,d,41);
  SHA1_R3(d,e,a,b,c,42); SHA1_R3(c,d,e,a,b,43);
  SHA1_R3(b,c,d,e,a,44); SHA1_R3(a,b,c,d,e,45);
  SHA1_R3(e,a,b,c,d,46); SHA1_R3(d,e,a,b,c,47);
  SHA1_R3(c,d,e,a,b,48); SHA1_R3(b,c,d,e,a,49);
  SHA1_R3(a,b,c,d,e,50); SHA1_R3(e,a,b,c,d,51);
  SHA1_R3(d,e,a,b,c,52); SHA1_R3(c,d,e,a,b,53);
  SHA1_R3(b,c,d,e,a,54); SHA1_R3(a,b,c,d,e,55);
  SHA1_R3(e,a,b,c,d,56); SHA1_R3(d,e,a,b,c,57);
  SHA1_R3(c,d,e,a,b,58); SHA1_R3(b,c,d,e,a,59);

  SHA1_R4(a,b,c,d,e,60); SHA1_R4(e,a,b,c,d,61);
  SHA1_R4(d,e,a,b,c,62); SHA1_R4(c,d,e,a,b,63);
  SHA1_R4(b,c,d,e,a,64); SHA1_R4(a,b,c,d,e,65);
  SHA1_R4(e,a,b,c,d,66); SHA1_R4(d,e,a,b,c,67);
  SHA1_R4(c,d,e,a,b,68); SHA1_R4(b,c,d,e,a,69);
  SHA1_R4(a,b,c,d,e,70); SHA1_R4(e,a,b,c,d,71);
  SHA1_R4(d,e,a,b,c,72); SHA1_R4(c,d,e,a,b,73);
  SHA1_R4(b,c,d,e,a,74); SHA1_R4(a,b,c,d,e,75);
  SHA1_R4(e,a,b,c,d,76); SHA1_R4(d,e,a,b,c,77);
  SHA1_R4(c,d,e,a,b,78); SHA1_R4(b,c,d,e,a,79);

  state[0] += a;
  state[1] += b;
  state[2] += c;
  state[3] += d;
  state[4] += e;
}

static void sha1_init(sha1_ctx_t * ctx)
{
  ctx->state[0] = 0x67452301;
  ctx->state[1] = 0xEFCDAB89;
  ctx->state[2] = 0x98BADCFE;
  ctx->state[3] = 0x10325476;
  ctx->state[4] = 0xC3D2E1F0;
  ctx->count[0] = 0;
  ctx->count[1] = 0;
}

static void sha1_update(sha1_ctx_t * ctx,
                         const unsigned char * data,
                         size_t len)
{
  size_t i, j;

  j = ctx->count[0];
  if ((ctx->count[0] += (uint32_t)(len << 3)) < j)
    ctx->count[1]++;
  ctx->count[1] += (uint32_t)(len >> 29);
  j = (j >> 3) & 63;

  if ((j + len) > 63)
  {
    memcpy(&ctx->buffer[j], data, (i = 64 - j));
    sha1_transform(ctx->state, ctx->buffer);
    for (; i + 63 < len; i += 64)
      sha1_transform(ctx->state, &data[i]);
    j = 0;
  }
  else
  {
    i = 0;
  }
  memcpy(&ctx->buffer[j], &data[i], len - i);
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
static void sha1_final(sha1_ctx_t * ctx, unsigned char digest[20])
{
  unsigned char finalcount[8];
  unsigned char c;
  unsigned int i;

  for (i = 0; i < 8; i++)
    finalcount[i] = (unsigned char)((ctx->count[(i >= 4) ? 0 : 1]
                     >> ((3 - (i & 3)) * 8)) & 255);

  c = 0200;
  sha1_update(ctx, &c, 1);
  while ((ctx->count[0] & 504) != 448)
  {
    c = 0000;
    sha1_update(ctx, &c, 1);
  }
  sha1_update(ctx, finalcount, 8);

  for (i = 0; i < 20; i++)
    digest[i] = (unsigned char)((ctx->state[i >> 2]
                 >> ((3 - (i & 3)) * 8)) & 255);
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

void sha1_compute(const unsigned char * data,
                  size_t len,
                  unsigned char digest[20])
{
  sha1_ctx_t ctx;

  sha1_init(&ctx);
  sha1_update(&ctx, data, len);
  sha1_final(&ctx, digest);
}
