/* gcvt with automatic memory allocation.
   Copyright (C) 2002-2004, 2007-2019 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License along
   with this program; if not, see <https://www.gnu.org/licenses/>.  */
  
// fork from glib/gnulib/vasnprintf.c

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gcvt.h"
#include "float+.h"

#if defined GCC_LINT || defined lint
#define IF_LINT(Code) Code
#else
#define IF_LINT(Code) /* empty */
#endif

#ifdef __GNUC__
#define G_GNUC_CHECK_VERSION(major, minor) \
    ((__GNUC__ > (major)) || \
     ((__GNUC__ == (major)) && \
      (__GNUC_MINOR__ >= (minor))))
#else
#define G_GNUC_CHECK_VERSION(major, minor) 0
#endif

#if G_GNUC_CHECK_VERSION(2, 8)
#define G_GNUC_EXTENSION __extension__
#else
#define G_GNUC_EXTENSION
#endif

#define G_GUINT64_CONSTANT(val)	(G_GNUC_EXTENSION (val##ULL))
#define G_MAXUINT64	G_GUINT64_CONSTANT(0xffffffffffffffff)
#define G_MAXSIZE	G_MAXUINT64

#define GMP_LIMB_BITS 32
typedef unsigned int mp_limb_t;

typedef unsigned long long mp_twolimb_t;
#define GMP_TWOLIMB_BITS 64

typedef struct mpn_t
{
  size_t nlimbs;
  mp_limb_t *limbs; /* Bits in little-endian order, allocated with malloc().  */
} mpn_t;

static int floorlog10 (double x)
{
  int exp;
  double y;
  double z;
  double l;

  /* Split into exponential part and mantissa.  */
  y = frexp (x, &exp);
  if (!(y >= 0.0 && y < 1.0))
    abort ();
  if (y == 0.0)
    return INT_MIN;
  if (y < 0.5)
    {
      while (y < (1.0 / (1 << (GMP_LIMB_BITS / 2)) / (1 << (GMP_LIMB_BITS / 2))))
        {
          y *= 1.0 * (1 << (GMP_LIMB_BITS / 2)) * (1 << (GMP_LIMB_BITS / 2));
          exp -= GMP_LIMB_BITS;
        }
      if (y < (1.0 / (1 << 16)))
        {
          y *= 1.0 * (1 << 16);
          exp -= 16;
        }
      if (y < (1.0 / (1 << 8)))
        {
          y *= 1.0 * (1 << 8);
          exp -= 8;
        }
      if (y < (1.0 / (1 << 4)))
        {
          y *= 1.0 * (1 << 4);
          exp -= 4;
        }
      if (y < (1.0 / (1 << 2)))
        {
          y *= 1.0 * (1 << 2);
          exp -= 2;
        }
      if (y < (1.0 / (1 << 1)))
        {
          y *= 1.0 * (1 << 1);
          exp -= 1;
        }
    }
  if (!(y >= 0.5 && y < 1.0))
    abort ();
  /* Compute an approximation for l = log2(x) = exp + log2(y).  */
  l = exp;
  z = y;
  if (z < 0.70710678118654752444)
    {
      z *= 1.4142135623730950488;
      l -= 0.5;
    }
  if (z < 0.8408964152537145431)
    {
      z *= 1.1892071150027210667;
      l -= 0.25;
    }
  if (z < 0.91700404320467123175)
    {
      z *= 1.0905077326652576592;
      l -= 0.125;
    }
  if (z < 0.9576032806985736469)
    {
      z *= 1.0442737824274138403;
      l -= 0.0625;
    }
  /* Now 0.95 <= z <= 1.01.  */
  z = 1 - z;
  /* log2(1-z) = 1/log(2) * (- z - z^2/2 - z^3/3 - z^4/4 - ...)
     Four terms are enough to get an approximation with error < 10^-7.  */
  l -= 1.4426950408889634074 * z * (1.0 + z * (0.5 + z * ((1.0 / 3) + z * 0.25)));
  /* Finally multiply with log(2)/log(10), yields an approximation for
     log10(x).  */
  l *= 0.30102999566398119523;
  /* Round down to the next integer.  */
  return (int) l + (l < 0 ? -1 : 0);
}

static void *multiply(mpn_t src1, mpn_t src2, mpn_t *dest)
{
  const mp_limb_t *p1;
  const mp_limb_t *p2;
  size_t len1;
  size_t len2;

  if (src1.nlimbs <= src2.nlimbs)
  {
    len1 = src1.nlimbs;
    p1 = src1.limbs;
    len2 = src2.nlimbs;
    p2 = src2.limbs;
  }
  else
  {
    len1 = src2.nlimbs;
    p1 = src2.limbs;
    len2 = src1.nlimbs;
    p2 = src1.limbs;
  }
  /* Now 0 <= len1 <= len2.  */
  if (len1 == 0)
  {
    /* src1 or src2 is zero.  */
    dest->nlimbs = 0;
    dest->limbs = (mp_limb_t *)malloc(1);
  }
  else
  {
    /* Here 1 <= len1 <= len2.  */
    size_t dlen;
    mp_limb_t *dp;
    size_t k, i, j;

    dlen = len1 + len2;
    dp = (mp_limb_t *)malloc(dlen * sizeof(mp_limb_t));
    if (dp == NULL)
      return NULL;
    for (k = len2; k > 0;)
      dp[--k] = 0;
    for (i = 0; i < len1; i++)
    {
      mp_limb_t digit1 = p1[i];
      mp_twolimb_t carry = 0;
      for (j = 0; j < len2; j++)
      {
        mp_limb_t digit2 = p2[j];
        carry += (mp_twolimb_t)digit1 * (mp_twolimb_t)digit2;
        carry += dp[i + j];
        dp[i + j] = (mp_limb_t)carry;
        carry = carry >> GMP_LIMB_BITS;
      }
      dp[i + len2] = (mp_limb_t)carry;
    }
    /* Normalise.  */
    while (dlen > 0 && dp[dlen - 1] == 0)
      dlen--;
    dest->nlimbs = dlen;
    dest->limbs = dp;
  }
  return dest->limbs;
}

static void *divide(mpn_t a, mpn_t b, mpn_t *q)
{
  const mp_limb_t *a_ptr = a.limbs;
  size_t a_len = a.nlimbs;
  const mp_limb_t *b_ptr = b.limbs;
  size_t b_len = b.nlimbs;
  mp_limb_t *roomptr;
  mp_limb_t *tmp_roomptr = NULL;
  mp_limb_t *q_ptr;
  size_t q_len;
  mp_limb_t *r_ptr;
  size_t r_len;

  /* Allocate room for a_len+2 digits.
     (Need a_len+1 digits for the real division and 1 more digit for the
     final rounding of q.)  */
  roomptr = (mp_limb_t *)malloc((a_len + 2) * sizeof(mp_limb_t));
  if (roomptr == NULL)
    return NULL;

  /* Normalise a.  */
  while (a_len > 0 && a_ptr[a_len - 1] == 0)
    a_len--;

  /* Normalise b.  */
  for (;;)
  {
    if (b_len == 0)
      /* Division by zero.  */
      abort();
    if (b_ptr[b_len - 1] == 0)
      b_len--;
    else
      break;
  }

  /* Here m = a_len >= 0 and n = b_len > 0.  */

  if (a_len < b_len)
  {
    /* m<n: trivial case.  q=0, r := copy of a.  */
    r_ptr = roomptr;
    r_len = a_len;
    memcpy(r_ptr, a_ptr, a_len * sizeof(mp_limb_t));
    q_ptr = roomptr + a_len;
    q_len = 0;
  }
  else if (b_len == 1)
  {
    /* n=1: single precision division.
       beta^(m-1) <= a < beta^m  ==>  beta^(m-2) <= a/b < beta^m  */
    r_ptr = roomptr;
    q_ptr = roomptr + 1;
    {
      mp_limb_t den = b_ptr[0];
      mp_limb_t remainder = 0;
      const mp_limb_t *sourceptr = a_ptr + a_len;
      mp_limb_t *destptr = q_ptr + a_len;
      size_t count;
      for (count = a_len; count > 0; count--)
      {
        mp_twolimb_t num =
            ((mp_twolimb_t)remainder << GMP_LIMB_BITS) | *--sourceptr;
        *--destptr = num / den;
        remainder = num % den;
      }
      /* Normalise and store r.  */
      if (remainder > 0)
      {
        r_ptr[0] = remainder;
        r_len = 1;
      }
      else
        r_len = 0;
      /* Normalise q.  */
      q_len = a_len;
      if (q_ptr[q_len - 1] == 0)
        q_len--;
    }
  }
  else
  {
    /* n>1: multiple precision division.
       beta^(m-1) <= a < beta^m, beta^(n-1) <= b < beta^n  ==>
       beta^(m-n-1) <= a/b < beta^(m-n+1).  */
    /* Determine s.  */
    size_t s;
    {
      mp_limb_t msd = b_ptr[b_len - 1]; /* = b[n-1], > 0 */
                                        /* Determine s = GMP_LIMB_BITS - integer_length (msd).
                                           Code copied from gnulib's integer_length.c.  */
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
      s = __builtin_clz(msd);
#else
#if defined DBL_EXPBIT0_WORD && defined DBL_EXPBIT0_BIT
      if (GMP_LIMB_BITS <= DBL_MANT_BIT)
      {
        /* Use 'double' operations.
           Assumes an IEEE 754 'double' implementation.  */
#define DBL_EXP_MASK ((DBL_MAX_EXP - DBL_MIN_EXP) | 7)
#define DBL_EXP_BIAS (DBL_EXP_MASK / 2 - 1)
#define NWORDS \
  ((sizeof(double) + sizeof(unsigned int) - 1) / sizeof(unsigned int))
        union
        {
          double value;
          unsigned int word[NWORDS];
        } m;

        /* Use a single integer to floating-point conversion.  */
        m.value = msd;

        s = GMP_LIMB_BITS - (((m.word[DBL_EXPBIT0_WORD] >> DBL_EXPBIT0_BIT) & DBL_EXP_MASK) - DBL_EXP_BIAS);
      }
      else
#undef NWORDS
#endif
      {
        s = 31;
        if (msd >= 0x10000)
        {
          msd = msd >> 16;
          s -= 16;
        }
        if (msd >= 0x100)
        {
          msd = msd >> 8;
          s -= 8;
        }
        if (msd >= 0x10)
        {
          msd = msd >> 4;
          s -= 4;
        }
        if (msd >= 0x4)
        {
          msd = msd >> 2;
          s -= 2;
        }
        if (msd >= 0x2)
        {
          msd = msd >> 1;
          s -= 1;
        }
      }
#endif
    }
    if (s > 0)
    {
      tmp_roomptr = (mp_limb_t *)malloc(b_len * sizeof(mp_limb_t));
      if (tmp_roomptr == NULL)
      {
        free(roomptr);
        return NULL;
      }
      {
        const mp_limb_t *sourceptr = b_ptr;
        mp_limb_t *destptr = tmp_roomptr;
        mp_twolimb_t accu = 0;
        size_t count;
        for (count = b_len; count > 0; count--)
        {
          accu += (mp_twolimb_t)*sourceptr++ << s;
          *destptr++ = (mp_limb_t)accu;
          accu = accu >> GMP_LIMB_BITS;
        }
        /* accu must be zero, since that was how s was determined.  */
        if (accu != 0)
          abort();
      }
      b_ptr = tmp_roomptr;
    }
    r_ptr = roomptr;
    if (s == 0)
    {
      memcpy(r_ptr, a_ptr, a_len * sizeof(mp_limb_t));
      r_ptr[a_len] = 0;
    }
    else
    {
      const mp_limb_t *sourceptr = a_ptr;
      mp_limb_t *destptr = r_ptr;
      mp_twolimb_t accu = 0;
      size_t count;
      for (count = a_len; count > 0; count--)
      {
        accu += (mp_twolimb_t)*sourceptr++ << s;
        *destptr++ = (mp_limb_t)accu;
        accu = accu >> GMP_LIMB_BITS;
      }
      *destptr++ = (mp_limb_t)accu;
    }
    q_ptr = roomptr + b_len;
    q_len = a_len - b_len + 1; /* q will have m-n+1 limbs */
    {
      size_t j = a_len - b_len;            /* m-n */
      mp_limb_t b_msd = b_ptr[b_len - 1];  /* b[n-1] */
      mp_limb_t b_2msd = b_ptr[b_len - 2]; /* b[n-2] */
      mp_twolimb_t b_msdd =                /* b[n-1]*beta+b[n-2] */
          ((mp_twolimb_t)b_msd << GMP_LIMB_BITS) | b_2msd;
      /* Division loop, traversed m-n+1 times.
         j counts down, b is unchanged, beta/2 <= b[n-1] < beta.  */
      for (;;)
      {
        mp_limb_t q_star;
        mp_limb_t c1;
        if (r_ptr[j + b_len] < b_msd) /* r[j+n] < b[n-1] ? */
        {
          /* Divide r[j+n]*beta+r[j+n-1] by b[n-1], no overflow.  */
          mp_twolimb_t num =
              ((mp_twolimb_t)r_ptr[j + b_len] << GMP_LIMB_BITS) | r_ptr[j + b_len - 1];
          q_star = num / b_msd;
          c1 = num % b_msd;
        }
        else
        {
          /* Overflow, hence r[j+n]*beta+r[j+n-1] >= beta*b[n-1].  */
          q_star = (mp_limb_t) ~(mp_limb_t)0; /* q* = beta-1 */
          if (r_ptr[j + b_len] > b_msd || (c1 = r_ptr[j + b_len - 1] + b_msd) < b_msd)
            goto subtract;
        }
        /* q_star = q*,
           c1 = (r[j+n]*beta+r[j+n-1]) - q* * b[n-1] (>=0, <beta).  */
        {
          mp_twolimb_t c2 = /* c1*beta+r[j+n-2] */
              ((mp_twolimb_t)c1 << GMP_LIMB_BITS) | r_ptr[j + b_len - 2];
          mp_twolimb_t c3 = /* b[n-2] * q* */
              (mp_twolimb_t)b_2msd * (mp_twolimb_t)q_star;
          /* While c2 < c3, increase c2 and decrease c3.
             Consider c3-c2.  While it is > 0, decrease it by
             b[n-1]*beta+b[n-2].  Because of b[n-1]*beta+b[n-2] >= beta^2/2
             this can happen only twice.  */
          if (c3 > c2)
          {
            q_star = q_star - 1; /* q* := q* - 1 */
            if (c3 - c2 > b_msdd)
              q_star = q_star - 1; /* q* := q* - 1 */
          }
        }
        if (q_star > 0)
        subtract:
        {
          /* Subtract r := r - b * q* * beta^j.  */
          mp_limb_t cr;
          {
            const mp_limb_t *sourceptr = b_ptr;
            mp_limb_t *destptr = r_ptr + j;
            mp_twolimb_t carry = 0;
            size_t count;
            for (count = b_len; count > 0; count--)
            {
              /* Here 0 <= carry <= q*.  */
              carry =
                  carry + (mp_twolimb_t)q_star * (mp_twolimb_t)*sourceptr++ + (mp_limb_t) ~(*destptr);
              /* Here 0 <= carry <= beta*q* + beta-1.  */
              *destptr++ = ~(mp_limb_t)carry;
              carry = carry >> GMP_LIMB_BITS; /* <= q* */
            }
            cr = (mp_limb_t)carry;
          }
          /* Subtract cr from r_ptr[j + b_len], then forget about
             r_ptr[j + b_len].  */
          if (cr > r_ptr[j + b_len])
          {
            /* Subtraction gave a carry.  */
            q_star = q_star - 1; /* q* := q* - 1 */
            /* Add b back.  */
            {
              const mp_limb_t *sourceptr = b_ptr;
              mp_limb_t *destptr = r_ptr + j;
              mp_limb_t carry = 0;
              size_t count;
              for (count = b_len; count > 0; count--)
              {
                mp_limb_t source1 = *sourceptr++;
                mp_limb_t source2 = *destptr;
                *destptr++ = source1 + source2 + carry;
                carry =
                    (carry
                         ? source1 >= (mp_limb_t)~source2
                         : source1 > (mp_limb_t)~source2);
              }
            }
            /* Forget about the carry and about r[j+n].  */
          }
        }
          /* q* is determined.  Store it as q[j].  */
          q_ptr[j] = q_star;
        if (j == 0)
          break;
        j--;
      }
    }
    r_len = b_len;
    /* Normalise q.  */
    if (q_ptr[q_len - 1] == 0)
      q_len--;
    /* Normalise r.  */
    while (r_len > 0 && r_ptr[r_len - 1] == 0)
      r_len--;
  }
  /* Compare r << 1 with b.  */
  if (r_len > b_len)
    goto increment_q;
  {
    size_t i;
    for (i = b_len;;)
    {
      mp_limb_t r_i =
          (i <= r_len && i > 0 ? r_ptr[i - 1] >> (GMP_LIMB_BITS - 1) : 0) | (i < r_len ? r_ptr[i] << 1 : 0);
      mp_limb_t b_i = (i < b_len ? b_ptr[i] : 0);
      if (r_i > b_i)
        goto increment_q;
      if (r_i < b_i)
        goto keep_q;
      if (i == 0)
        break;
      i--;
    }
  }
  if (q_len > 0 && ((q_ptr[0] & 1) != 0))
  /* q is odd.  */
  increment_q:
  {
    size_t i;
    for (i = 0; i < q_len; i++)
      if (++(q_ptr[i]) != 0)
        goto keep_q;
    q_ptr[q_len++] = 1;
  }
  keep_q:
    if (tmp_roomptr != NULL)
      free(tmp_roomptr);
  q->limbs = q_ptr;
  q->nlimbs = q_len;
  return roomptr;
}

size_t xsum (size_t size1, size_t size2)
{
  size_t sum = size1 + size2;
  return (sum >= size1 ? sum : G_MAXSIZE);
}

static char *convert_to_decimal(mpn_t a, size_t extra_zeroes)
{
  mp_limb_t *a_ptr = a.limbs;
  size_t a_len = a.nlimbs;
  /* 0.03345 is slightly larger than log(2)/(9*log(10)).  */
  size_t c_len = 9 * ((size_t)(a_len * (GMP_LIMB_BITS * 0.03345f)) + 1);
  /* We need extra_zeroes bytes for zeroes, followed by c_len bytes for the
     digits of a, followed by 1 byte for the terminating NUL.  */
  char *c_ptr = (char *)malloc(xsum(xsum(extra_zeroes, c_len), 1));
  if (c_ptr != NULL)
  {
    char *d_ptr = c_ptr;
    for (; extra_zeroes > 0; extra_zeroes--)
      *d_ptr++ = '0';
    while (a_len > 0)
    {
      /* Divide a by 10^9, in-place.  */
      mp_limb_t remainder = 0;
      mp_limb_t *ptr = a_ptr + a_len;
      size_t count;
      for (count = a_len; count > 0; count--)
      {
        mp_twolimb_t num =
            ((mp_twolimb_t)remainder << GMP_LIMB_BITS) | *--ptr;
        *ptr = num / 1000000000;
        remainder = num % 1000000000;
      }
      /* Store the remainder as 9 decimal digits.  */
      for (count = 9; count > 0; count--)
      {
        *d_ptr++ = '0' + (remainder % 10);
        remainder = remainder / 10;
      }
      /* Normalize a.  */
      if (a_ptr[a_len - 1] == 0)
        a_len--;
    }
    /* Remove leading zeroes.  */
    while (d_ptr > c_ptr && d_ptr[-1] == '0')
      d_ptr--;
    /* But keep at least one zero.  */
    if (d_ptr == c_ptr)
      *d_ptr++ = '0';
    /* Terminate the string.  */
    *d_ptr = '\0';
  }
  return c_ptr;
}

static void *decode_double(double x, int *ep, mpn_t *mp)
{
  mpn_t m;
  int exp;
  double y;
  size_t i;

  /* Allocate memory for result.  */
  m.nlimbs = (DBL_MANT_BIT + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS;
  m.limbs = (mp_limb_t *)malloc(m.nlimbs * sizeof(mp_limb_t));
  if (m.limbs == NULL)
    return NULL;
  /* Split into exponential part and mantissa.  */
  y = frexp(x, &exp);
  if (!(y >= 0.0 && y < 1.0))
    abort();
    /* x = 2^exp * y = 2^(exp - DBL_MANT_BIT) * (y * 2^DBL_MANT_BIT), and the
       latter is an integer.  */
    /* Convert the mantissa (y * 2^DBL_MANT_BIT) to a sequence of limbs.
       I'm not sure whether it's safe to cast a 'double' value between
       2^31 and 2^32 to 'unsigned int', therefore play safe and cast only
       'double' values between 0 and 2^16 (to 'unsigned int' or 'int',
       doesn't matter).  */
#if (DBL_MANT_BIT % GMP_LIMB_BITS) != 0
#if (DBL_MANT_BIT % GMP_LIMB_BITS) > GMP_LIMB_BITS / 2
  {
    mp_limb_t hi, lo;
    y *= (mp_limb_t)1 << (DBL_MANT_BIT % (GMP_LIMB_BITS / 2));
    hi = (int)y;
    y -= hi;
    if (!(y >= 0.0 && y < 1.0))
      abort();
    y *= (mp_limb_t)1 << (GMP_LIMB_BITS / 2);
    lo = (int)y;
    y -= lo;
    if (!(y >= 0.0 && y < 1.0))
      abort();
    m.limbs[DBL_MANT_BIT / GMP_LIMB_BITS] = (hi << (GMP_LIMB_BITS / 2)) | lo;
  }
#else
  {
    mp_limb_t d;
    y *= (mp_limb_t)1 << (DBL_MANT_BIT % GMP_LIMB_BITS);
    d = (int)y;
    y -= d;
    if (!(y >= 0.0 && y < 1.0))
      abort();
    m.limbs[DBL_MANT_BIT / GMP_LIMB_BITS] = d;
  }
#endif
#endif
  for (i = DBL_MANT_BIT / GMP_LIMB_BITS; i > 0;)
  {
    mp_limb_t hi, lo;
    y *= (mp_limb_t)1 << (GMP_LIMB_BITS / 2);
    hi = (int)y;
    y -= hi;
    if (!(y >= 0.0 && y < 1.0))
      abort();
    y *= (mp_limb_t)1 << (GMP_LIMB_BITS / 2);
    lo = (int)y;
    y -= lo;
    if (!(y >= 0.0 && y < 1.0))
      abort();
    m.limbs[--i] = (hi << (GMP_LIMB_BITS / 2)) | lo;
  }
  if (!(y == 0.0))
    abort();
  /* Normalise.  */
  while (m.nlimbs > 0 && m.limbs[m.nlimbs - 1] == 0)
    m.nlimbs--;
  *mp = m;
  *ep = exp - DBL_MANT_BIT;
  return m.limbs;
}

static char *scale10_round_decimal_decoded(int e, mpn_t m, void *memory, int n)
{
  int s;
  size_t extra_zeroes;
  unsigned int abs_n;
  unsigned int abs_s;
  mp_limb_t *pow5_ptr;
  size_t pow5_len;
  unsigned int s_limbs;
  unsigned int s_bits;
  mpn_t pow5;
  mpn_t z;
  void *z_memory;
  char *digits;

  if (memory == NULL)
    return NULL;
  /* x = 2^e * m, hence
     y = round (2^e * 10^n * m) = round (2^(e+n) * 5^n * m)
       = round (2^s * 5^n * m).  */
  s = e + n;
  extra_zeroes = 0;
  /* Factor out a common power of 10 if possible.  */
  if (s > 0 && n > 0)
  {
    extra_zeroes = (s < n ? s : n);
    s -= extra_zeroes;
    n -= extra_zeroes;
  }
  /* Here y = round (2^s * 5^n * m) * 10^extra_zeroes.
     Before converting to decimal, we need to compute
     z = round (2^s * 5^n * m).  */
  /* Compute 5^|n|, possibly shifted by |s| bits if n and s have the same
     sign.  2.322 is slightly larger than log(5)/log(2).  */
  abs_n = (n >= 0 ? n : -n);
  abs_s = (s >= 0 ? s : -s);
  pow5_ptr = (mp_limb_t *)malloc(((int)(abs_n * (2.322f / GMP_LIMB_BITS)) + 1 + abs_s / GMP_LIMB_BITS + 1) * sizeof(mp_limb_t));
  if (pow5_ptr == NULL)
  {
    free(memory);
    return NULL;
  }
  /* Initialize with 1.  */
  pow5_ptr[0] = 1;
  pow5_len = 1;
  /* Multiply with 5^|n|.  */
  if (abs_n > 0)
  {
    static mp_limb_t const small_pow5[13 + 1] =
        {
            1, 5, 25, 125, 625, 3125, 15625, 78125, 390625, 1953125, 9765625,
            48828125, 244140625, 1220703125};
    unsigned int n13;
    for (n13 = 0; n13 <= abs_n; n13 += 13)
    {
      mp_limb_t digit1 = small_pow5[n13 + 13 <= abs_n ? 13 : abs_n - n13];
      size_t j;
      mp_twolimb_t carry = 0;
      for (j = 0; j < pow5_len; j++)
      {
        mp_limb_t digit2 = pow5_ptr[j];
        carry += (mp_twolimb_t)digit1 * (mp_twolimb_t)digit2;
        pow5_ptr[j] = (mp_limb_t)carry;
        carry = carry >> GMP_LIMB_BITS;
      }
      if (carry > 0)
        pow5_ptr[pow5_len++] = (mp_limb_t)carry;
    }
  }
  s_limbs = abs_s / GMP_LIMB_BITS;
  s_bits = abs_s % GMP_LIMB_BITS;
  if (n >= 0 ? s >= 0 : s <= 0)
  {
    /* Multiply with 2^|s|.  */
    if (s_bits > 0)
    {
      mp_limb_t *ptr = pow5_ptr;
      mp_twolimb_t accu = 0;
      size_t count;
      for (count = pow5_len; count > 0; count--)
      {
        accu += (mp_twolimb_t)*ptr << s_bits;
        *ptr++ = (mp_limb_t)accu;
        accu = accu >> GMP_LIMB_BITS;
      }
      if (accu > 0)
      {
        *ptr = (mp_limb_t)accu;
        pow5_len++;
      }
    }
    if (s_limbs > 0)
    {
      size_t count;
      for (count = pow5_len; count > 0;)
      {
        count--;
        pow5_ptr[s_limbs + count] = pow5_ptr[count];
      }
      for (count = s_limbs; count > 0;)
      {
        count--;
        pow5_ptr[count] = 0;
      }
      pow5_len += s_limbs;
    }
    pow5.limbs = pow5_ptr;
    pow5.nlimbs = pow5_len;
    if (n >= 0)
    {
      /* Multiply m with pow5.  No division needed.  */
      z_memory = multiply(m, pow5, &z);
    }
    else
    {
      /* Divide m by pow5 and round.  */
      z_memory = divide(m, pow5, &z);
    }
  }
  else
  {
    pow5.limbs = pow5_ptr;
    pow5.nlimbs = pow5_len;
    if (n >= 0)
    {
      /* n >= 0, s < 0.
         Multiply m with pow5, then divide by 2^|s|.  */
      mpn_t numerator;
      mpn_t denominator;
      void *tmp_memory;
      tmp_memory = multiply(m, pow5, &numerator);
      if (tmp_memory == NULL)
      {
        free(pow5_ptr);
        free(memory);
        return NULL;
      }
      /* Construct 2^|s|.  */
      {
        mp_limb_t *ptr = pow5_ptr + pow5_len;
        size_t i;
        for (i = 0; i < s_limbs; i++)
          ptr[i] = 0;
        ptr[s_limbs] = (mp_limb_t)1 << s_bits;
        denominator.limbs = ptr;
        denominator.nlimbs = s_limbs + 1;
      }
      z_memory = divide(numerator, denominator, &z);
      free(tmp_memory);
    }
    else
    {
      /* n < 0, s > 0.
         Multiply m with 2^s, then divide by pow5.  */
      mpn_t numerator;
      mp_limb_t *num_ptr;
      num_ptr = (mp_limb_t *)malloc((m.nlimbs + s_limbs + 1) * sizeof(mp_limb_t));
      if (num_ptr == NULL)
      {
        free(pow5_ptr);
        free(memory);
        return NULL;
      }
      {
        mp_limb_t *destptr = num_ptr;
        {
          size_t i;
          for (i = 0; i < s_limbs; i++)
            *destptr++ = 0;
        }
        if (s_bits > 0)
        {
          const mp_limb_t *sourceptr = m.limbs;
          mp_twolimb_t accu = 0;
          size_t count;
          for (count = m.nlimbs; count > 0; count--)
          {
            accu += (mp_twolimb_t)*sourceptr++ << s_bits;
            *destptr++ = (mp_limb_t)accu;
            accu = accu >> GMP_LIMB_BITS;
          }
          if (accu > 0)
            *destptr++ = (mp_limb_t)accu;
        }
        else
        {
          const mp_limb_t *sourceptr = m.limbs;
          size_t count;
          for (count = m.nlimbs; count > 0; count--)
            *destptr++ = *sourceptr++;
        }
        numerator.limbs = num_ptr;
        numerator.nlimbs = destptr - num_ptr;
      }
      z_memory = divide(numerator, pow5, &z);
      free(num_ptr);
    }
  }
  free(pow5_ptr);
  free(memory);

  /* Here y = round (x * 10^n) = z * 10^extra_zeroes.  */

  if (z_memory == NULL)
    return NULL;
  digits = convert_to_decimal(z, extra_zeroes);
  free(z_memory);
  return digits;
}

static char *scale10_round_decimal_double(double x, int n)
{
  int e IF_LINT(= 0);
  mpn_t m;
  void *memory = decode_double(x, &e, &m);
  return scale10_round_decimal_decoded(e, m, memory, n);
}

int _gfcvt_s(
    char* buffer,
    size_t buffer_size,
    double value,
    int ndigit,
    int* decpt,
    int* sign
) {
  if (signbit(value)) {
    value = -value;
    *sign = 1;
  } else {
    *sign = 0;
  }
  char *buf = scale10_round_decimal_double(value, ndigit);
  int len = strlen(buf);
  if (len > buffer_size) {
    return -1;
  }
  *decpt = len - ndigit;
  int offset = 0;
  for (int i = len - 1;i >= 0; i--) {
    int j = len - i - 1;
    buffer[j] = buf[i];
  }
  free(buf);
  buffer[len] = '\0';
  return 0;
}


int _gecvt_s(
    char* buffer,
    size_t buffer_size,
    double value,
    int p_ndigit,
    int* decpt,
    int* sign
) {
  if (signbit(value)) {
    value = -value;
    *sign = 1;
  } else {
    *sign = 0;
  }
  char *digits = NULL;
  size_t ndigits = 0;
  int exponent = floorlog10(value);
  int precision = p_ndigit;
  int adjusted = 0;
  for (;;) {
    digits = scale10_round_decimal_double(value, precision - exponent);
    if (digits == NULL) {
      return -1;
    }
    ndigits = strlen(digits);
    if (ndigits == precision + 1)
      break;
    free(digits);
    if (adjusted)
      return -2;
    if (ndigits == precision)
      exponent -= 1;
    else
      exponent += 1;
    adjusted = 1;
  }
  *decpt = exponent+1;
  int len = strlen(digits);
  if (len > buffer_size) {
    return -3;
  }
  int offset = 0;
  for (int i = len - 1;i >= 0; i--) {
    int j = len - i - 1;
    buffer[j] = digits[i];
  }
  free(digits);
  buffer[len] = '\0';
  return 0;
}