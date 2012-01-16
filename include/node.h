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

#include <iostream>
#include <vector>
#include <map>

namespace Trie{

template <class _data_type>
class Node:
public std::map<std::size_t,Node<_data_type> >
{
        public:
        typedef _data_type                                      data_type;
        typedef std::map<std::size_t,Node<data_type> >          base;
        typedef std::vector<std::size_t>                        key_type;

        public:
        Node(); // constructor
        data_type& get();
        unsigned lookups()const;
        void set_lookup_counter(unsigned);
        void decrease_lookup_counter(unsigned);
        unsigned size_subtree()const;
        Node& child(std::size_t);
        std::size_t children()const;
        std::ostream& print_map(std::ostream&)const;

        // the data stored in the node
        private:
        data_type _d; // the data associated with this node
        unsigned _num_lookups; // how many times this node was traversed
};

} // namespace Trie

#include "node_impl.h"
