// Copyright (c) 2009 INRIA Sophia-Antipolis (France),
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://vfisikop@scm.gforge.inria.fr/svn/cgal/branches/features/Triangulation-higher_dimensions-odevil_shornus/Triangulation/include/CGAL/TDS_full_cell_default_storage_policy.h $
// $Id: TDS_full_cell_default_storage_policy.h 63441 2011-05-10 14:10:39Z shornus $
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_TDS_FULL_CELL_DEFAULT_STORAGE_POLICY_H
#define CGAL_TDS_FULL_CELL_DEFAULT_STORAGE_POLICY_H

#include <CGAL/Dimension.h>
#include <CGAL/Compact_container.h>
#include <CGAL/internal/Static_or_dynamic_array.h>

namespace CGAL {

// POLICY TAG

struct TDS_full_cell_default_storage_policy {}; // stores no additional data. Uses XOR trick.

template< typename V, typename S, typename D, typename StoragePolicy >
struct TFC_data; // TFC = Triangulation Full Cell

template< typename Vertex_handle, typename Full_cell_handle, typename Dimen >
struct TFC_data< Vertex_handle, Full_cell_handle, Dimen, TDS_full_cell_default_storage_policy >
{
    typedef typename internal::Dimen_plus_one<Dimen>::type Dimen_plus;
    typedef typename internal::S_or_D_array< Vertex_handle, Dimen_plus, true >     Vertex_handle_array;
    typedef typename internal::S_or_D_array< Full_cell_handle, Dimen_plus >    Full_cell_handle_array;

    Vertex_handle_array  vertices_;
    Full_cell_handle_array neighbors_;

    TFC_data(const int dmax)
    : vertices_(dmax+1), neighbors_(dmax+1)
    {}
    void*   for_compact_container() const { return vertices_.for_compact_container(); }
    void* & for_compact_container()       { return vertices_.for_compact_container(); }
    int dimension() const { return ( vertices_.size() - 1 ); }
    void set_mirror_index(const int, const int) {}
    typedef unsigned long Xor_type; // TODO: choose a type that works for sure in 32bits and 64bits machines
    Xor_type xor_of_vertices(const int cur_dim) const
    {
        Xor_type result(0);
        for( int i = 0; i <= cur_dim; ++i )
            result ^= reinterpret_cast<Xor_type>(&(*vertices_[i]));
        return result;
    }
    // ASSUMES |*this| is indeed a neighbor of neighbor(i):
    // NOT correct when the hole (in insert_in_hole) is doubly covered.
    int mirror_index(const int i) const
    {
        int index = 0;
        Full_cell_handle n = neighbors_[i];
        Full_cell_handle o = n->neighbor(index);
        while( &(o->combinatorics_) != this )
            o = n->neighbor(++index);
        return index;
    }
    Vertex_handle mirror_vertex(const int i, const int cur_dim) const
    {
        Xor_type opp_vertex = xor_of_vertices(cur_dim)
            ^ neighbors_[i]->xor_of_vertices(cur_dim)
            ^ reinterpret_cast<Xor_type>(&(*vertices_[i]));
        Vertex_handle mirror;
        typedef typename Vertex_handle::pointer pointer;
        // mirror.set_pointer(reinterpret_cast<pointer>(opp_vertex));
        mirror = Compact_container<typename Vertex_handle::value_type>
            ::s_iterator_to(*(reinterpret_cast<pointer>(opp_vertex)));
        return mirror;
    }
    void swap_vertices(const int d1, const int d2)
    {
        std::swap(vertices_[d1], vertices_[d2]);
        std::swap(neighbors_[d1], neighbors_[d2]);
    }
};

} //namespace CGAL

#endif // CGAL_TDS_FULL_CELL_DEFAULT_STORAGE_POLICY_H
