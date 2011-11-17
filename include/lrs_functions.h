/* Copyright 2011 National and Kapodistrian University of Athens, Greece.
 *
 * This file is part of respol.
 *
 * Respol is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Respol is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * See the file COPYING.LESSER for the text of the GNU Lesser General
 * Public License.  If you did not receive this file along with respol, see
 * <http: *www.gnu.org/licenses/>.
 */

/* This C header file must be included from an extern "C" {} block. The
 * header gmp.h must be included before and outside the extern block to
 * avoid "c linkage" errors. One does not need to include this file, since
 * it is included from lrs_cgal.h, the C++ header file that contains the
 * interface functions. */

/* Moreover, before including lrslib.h, GMP must be defined. */
#include <lrslib.h>

/* These two functions are examples, their implementations are in
 * lrs_functions.c but they are not normally compiled. */
void makecyclic(lrs_dic*,lrs_dat*);
int ch();
