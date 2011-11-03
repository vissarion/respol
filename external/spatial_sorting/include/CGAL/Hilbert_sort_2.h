// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Triangulation-higher_dimensions-odevil_shornus/Spatial_sorting/include/CGAL/Hilbert_sort_2.h $
// $Id: Hilbert_sort_2.h 63205 2011-04-29 09:23:32Z odevil $
//
// Author(s)     : Olivier Devillers

#ifndef CGAL_HILBERT_SORT_2_H
#define CGAL_HILBERT_SORT_2_H

#include <CGAL/Hilbert_policy_tags.h>
#include <CGAL/Hilbert_sort_median_2.h>
#include <CGAL/Hilbert_sort_middle_2.h>

namespace CGAL {

template <class K,  class Hilbert_policy >
  class Hilbert_sort_2;

template <class K>  
  class Hilbert_sort_2<K, Hilbert_sort_median_policy >
  : public Hilbert_sort_median_2<K>
{
 public:
  Hilbert_sort_2 (const K &k=K() , std::ptrdiff_t limit=1 )
   : Hilbert_sort_median_2<K> (k,limit)
    {}
};

template <class K>
  class Hilbert_sort_2<K, Hilbert_sort_middle_policy >
  : public Hilbert_sort_middle_2<K>
{
 public:
 Hilbert_sort_2 (const K &k=K() , std::ptrdiff_t limit=1 )
   : Hilbert_sort_middle_2<K> (k,limit)
    {}
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_2_H