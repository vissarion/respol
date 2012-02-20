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

#ifndef SORT_SWAP_H
#define SORT_SWAP_H

#include <cassert>
#include <boost/version.hpp>
#if BOOST_VERSION < 104300
#include <boost/detail/algorithm.hpp>
#else
#include <boost/range/algorithm_ext/is_sorted.hpp>
#endif
#include <boost/tuple/tuple.hpp>

// This function takes as parameter an index vector v (usually, a
// std::vector<size_t>) and returns a pair, containing the sorted vector s
// and a boolean b. b is true iff s can be obtained from v with an even
// number of swaps.
// TODO: This function is naively implemented, it sorts using a selection
// sort and computes the number of swaps on the fly. It will be better to
// implement a subquadratic algorithm that counts the swaps.
template <class Index>
boost::tuple<Index,bool> sort_swap(const Index &v){
        typedef typename Index::value_type                      elt_t;
        size_t n=v.size();
        Index s(v); // The sorted vector.
        bool swaps=true; // The number of swaps used is even.
        elt_t tmp; // The temporary for swaps.
        size_t min; // The temporary to store the minimum element.
        // Selection sort, O(n^2) but easy to implement.
        for(size_t i=0;i+1<n;++i){
                min=i;
                // Find the minimum in v[i+1...].
                for(size_t j=i+1;j<n;++j)
                        if(s[j]<s[min])
                                min=j;
                // Swap s[i] and s[min] and update swaps if necessary.
                if(min!=i){
                        tmp=s[min];
                        s[min]=s[i];
                        s[i]=tmp;
                        swaps=!swaps;
                }
        }
        assert(boost::is_sorted(s));
        assert(s.size()==v.size());
        return boost::make_tuple(boost::ref(s),swaps);
}

// This function does the same as the previous one. Additionally, it
// returns a permutation vector showing the permutation of the indices of
// the input vector.
template <class Index>
boost::tuple<Index,bool,Index> sort_swap_permutation(const Index &v){
        typedef typename Index::value_type                      elt_t;
        size_t n=v.size();
        Index s(v),permutation;
        permutation.reserve(n);
        for(elt_t p=0;p<n;++p)
                permutation.push_back(p);
        bool swaps=true; // The number of swaps used is even.
        elt_t tmp;
        size_t min;
        for(size_t i=0;i+1<n;++i){
                min=i;
                for(size_t j=i+1;j<n;++j)
                        if(s[j]<s[min])
                                min=j;
                if(min!=i){
                        tmp=s[min];
                        s[min]=s[i];
                        s[i]=tmp;
                        swaps=!swaps;
                        tmp=permutation[min];
                        permutation[min]=permutation[i];
                        permutation[i]=tmp;
                }
        }
        assert(boost::is_sorted(s));
        assert(s.size()==v.size());
        assert(permutation.size()==v.size());
        return boost::make_tuple(boost::ref(s),swaps,boost::ref(permutation));
}

// This function does the same as sort_swap, but it modifies the input.
template <class Index>
bool sort_swap_inplace(Index &s){
        typedef typename Index::value_type                      elt_t;
        size_t n=s.size();
        bool swaps=true; // The number of swaps used is even.
        elt_t tmp; // The temporary for swaps.
        size_t min; // The temporary to store the minimum element.
        // Selection sort, O(n^2) but easy to implement.
        for(size_t i=0;i+1<n;++i){
                min=i;
                // Find the minimum in v[i+1...].
                for(size_t j=i+1;j<n;++j)
                        if(s[j]<s[min])
                                min=j;
                // Swap s[i] and s[min] and update swaps if necessary.
                if(min!=i){
                        tmp=s[min];
                        s[min]=s[i];
                        s[i]=tmp;
                        swaps=!swaps;
                }
        }
        assert(boost::is_sorted(s));
        return swaps;
}

#endif // SORT_SWAP_H
// vim: ts=2
