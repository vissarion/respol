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

#include "node.h"
#include <vector>

namespace Trie{

template <class _data_type>
class Trie{
        private:
        typedef _data_type                                      data_type;
        typedef Node<data_type>                                 node_type;
        typedef typename node_type::key_type                    key_type;

        private:
        node_type root;

        public:
        data_type& operator[](const key_type&)const;
        std::size_t count(const key_type&)const;
        void erase(const key_type&);
        unsigned use_count(const key_type&)const;
        bool empty()const;
        void clear();
};

} // namespace Trie

#include "trie_impl.h"
