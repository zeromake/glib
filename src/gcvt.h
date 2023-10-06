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

#ifndef _GCVT_H
#define _GCVT_H

#include <stdint.h>

int _gfcvt_s(
    char* buffer,
    size_t buffer_size,
    double value,
    int ndigit,
    int* decpt,
    int* sign
);

int _gecvt_s(
    char* buffer,
    size_t buffer_size,
    double value,
    int ndigit,
    int* decpt,
    int* sign
);
#endif
