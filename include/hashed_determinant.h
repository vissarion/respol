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

#ifndef HASHED_DETERMINANT_H
#define HASHED_DETERMINANT_H

// This file is the dirtiest part of the code. Depending on the compilation
// definitions, a determinant algorithm will be chosen. Then a struct will
// act as a template typedef (which is not allowed in C++).

#ifdef USE_BIRD_DET
  #include "hashed_determinant_bird.h"
  #define _HD HashedDeterminantBird
#elif defined USE_CGAL_DET
  #include "hashed_determinant_cgal.h"
  #define _HD HashedDeterminantCGAL
#elif defined USE_CGAL_DET_2
  #include "hashed_determinant_cgal_2.h"
  #define _HD HashedDeterminantCGAL2
#elif defined USE_EIGEN_DET
  #include "hashed_determinant_eigen.h"
  #define _HD HashedDeterminantEigen
#elif defined USE_LINBOX_DET
  #include "hashed_determinant_linbox.h"
  #define _HD HashedDeterminantLinbox
#elif defined USE_ORIENTATION_DET
  #include "hashed_determinant_orientation.h"
  #define _HD HashedDeterminantOrientation
#else // by default, use Laplace determinants
  #include "hashed_determinant_laplace.h"
  #define _HD HashedDeterminantLaplace
#endif

template <class _NT>
struct HashedDeterminant{
        typedef _HD<_NT>                                Table;
};

#endif // HASHED_DETERMINANT_H
