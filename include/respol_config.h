// Copyright 2011-2012 National and Kapodistrian University of Athens,
// Greece.
//
// This file is part of respol.
//
// Respol is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// Respol is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with respol, see
// <http://www.gnu.org/licenses/>.

#ifndef RESPOL_CONFIG_H
#define RESPOL_CONFIG_H

#include <fstream>

namespace ResPol{

struct config{
  int verbose;
  bool read_from_file;
  bool output_f_vector;
  size_t polytope_type; // 0 for resultant or 1 for secondary
  std::ifstream inp;
};

} // namespace ResPol

#endif // RESPOL_CONFIG_H
