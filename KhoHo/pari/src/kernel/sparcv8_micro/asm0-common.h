/* Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */

/* This file is common to SuperSparc and MicroSparc */
/*
ASM addll mulll
NOASM bfffo
*/
#ifdef ASMINLINE

#define LOCAL_HIREMAINDER  ulong hiremainder
#define LOCAL_OVERFLOW     ulong overflow

#define addll(a,b) \
__extension__ ({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ( "addcc %2,%3,%0; \
          addx  %%g0,%%g0,%1" \
         : "=r" (__value), "=r" (overflow) \
         : "r" (__arg1), "r" (__arg2) \
         : "cc"); \
__value; })

#define addllx(a,b) \
__extension__ ({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ( "subcc %%g0,%1,%%g0; \
          addxcc %2,%3,%0; \
          addx  %%g0,%%g0,%1" \
         : "=r" (__value), "=r" (overflow) \
         : "r" (__arg1), "r" (__arg2), "1" (overflow) \
         : "cc"); \
__value; })

#define subll(a,b) \
__extension__ ({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ( "subcc %2,%3,%0; \
          addx  %%g0,%%g0,%1" \
         : "=r" (__value), "=r" (overflow) \
         : "r" (__arg1), "r" (__arg2) \
         : "cc"); \
__value; })

#define subllx(a,b) \
__extension__ ({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ( "subcc %%g0,%1,%%g0; \
          subxcc %2,%3,%0; \
          addx  %%g0,%%g0,%1" \
         : "=r" (__value), "=r" (overflow) \
         : "r" (__arg1), "r" (__arg2), "1" (overflow) \
         : "cc"); \
__value; })

#define mulll(a,b) \
__extension__ ({ ulong __value, __arg1 = (a), __arg2 = (b); \
   __asm__ ( "umul %2,%3,%0; \
          rd  %%y,%1" \
         : "=r" (__value), "=r" (hiremainder) \
         : "r" (__arg1), "r" (__arg2));        \
__value;})

#define addmul(a,b) \
__extension__ ({ ulong __value, __arg1 = (a), __arg2 = (b), __tmp; \
   __asm__ ( "umul %3,%4,%0; \
          rd  %%y,%2; \
          addcc %0,%1,%0; \
          addx %%g0,%2,%1" \
         : "=&r" (__value), "=&r" (hiremainder), "=&r" (__tmp) \
         : "r" (__arg1), "r" (__arg2), "1" (hiremainder) \
         : "cc");        \
__value;})
#endif
