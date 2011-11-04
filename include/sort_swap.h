// Copyright 2011 National and Kapodistrian University of Athens, Greece.
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

#ifndef SORT_SWAP_H
#define SORT_SWAP_H

#include <cassert>
#include <boost/range/algorithm_ext/is_sorted.hpp>

// This function takes as parameter an index vector v (usually, a
// std::vector<size_t>) and returns a pair, containing the sorted vector s
// and a boolean b. b is true iff s can be obtained from v with an even
// number of swaps.
// TODO: This function is naively implemented, it sorts using a bubble sort
// and computes the number of swaps on the fly. It will be better to
// implement a subquadratic algorithn that counts the swaps.
template <class Index>
std::pair<Index,bool> sort_swap(const Index &v){
        typedef typename Index::value_type                      elt_t;
        Index s(v); // the sorted vector
        size_t swaps=0; // the number of swaps used
        elt_t tmp; // the temporary for swaps
        size_t min; // the temporary to store the minimum element
        // bubble sort, n^2 but easy to implement
        for(size_t i=0;i<s.size()-1;++i){
                // find the minimum in v[i+1...]
                min=i;
                for(size_t j=i+1;j<s.size();++j)
                        if(s[j]<s[min])
                                min=j;
                // swap s[i] and s[min], and update swaps if necessary
                if(min!=i){
                        tmp=s[min];
                        s[min]=s[i];
                        s[i]=tmp;
                        ++swaps;
                }
        }
        assert(boost::is_sorted(s));
        return std::make_pair(s,!(swaps%2));
}

#endif // SORT_SWAP_H
// vim: ts=2
