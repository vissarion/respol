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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Triangulation-higher_dimensions-odevil_shornus/Spatial_sorting/include/CGAL/Hilbert_sort_base.h $
// $Id: Hilbert_sort_base.h 63205 2011-04-29 09:23:32Z odevil $
//
// Author(s)     : Christophe Delage

#ifndef CGAL_HILBERT_SORT_BASE_H
#define CGAL_HILBERT_SORT_BASE_H

#include <CGAL/basic.h>
#include <algorithm>

namespace CGAL {

namespace internal {

    template <class RandomAccessIterator, class Cmp>
    RandomAccessIterator
    hilbert_split (RandomAccessIterator begin, RandomAccessIterator end,
                   Cmp cmp = Cmp ())
    {
        if (begin >= end) return begin;

        RandomAccessIterator middle = begin + (end - begin) / 2;
        std::nth_element (begin, middle, end, cmp);
        return middle;
    }
}

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_BASE_H
