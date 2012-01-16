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

namespace Trie{

// This function returns the a reference to the element associated with
// idx. If it does not exist, the path in the tree to access it is created
// and a reference to the place where it will be stored is returned.
template <class _Data>
typename Trie<_Data>::data_type&
Trie<_Data>::operator[](const typename Trie<_Data>::key_type &idx)const{
        // We work here with pointers to avoid copying of the data (doing
        // so would lead to incorrect results).
        node_type *n=const_cast<node_type*>(&root);
        for(std::size_t i=0;i<idx.size();++i)
                n=&(n->child(idx[i]));
        return n->get();
}

// This function returns 1 iff the node identified by idx exists. The trie
// is not modified.
template <class _Data>
std::size_t
Trie<_Data>::count(const typename Trie<_Data>::key_type &idx)const{
        node_type *n=const_cast<node_type*>(&root);
        std::size_t i=0;
        for(;i<idx.size()-1;++i){
                if(n->count(idx[i])==0)
                        return 0;
                n=&(n->child(idx[i]));
        }
        return n->count(idx[i]);
}

// This erases the node pointed to idx and the subtree below it.
template <class _Data>
void
Trie<_Data>::erase(const typename Trie<_Data>::key_type &idx){
        std::size_t size=idx.size()-1;
        node_type *nodes[size];
        nodes[0]=&root;
        for(std::size_t i=1;i<idx.size();++i){
                //if(nodes[i-1]->count(idx[i-1])==0)
                //        // idx is not associated in the trie
                //        return;
                //else
                        nodes[i]=&(nodes[i-1]->child(idx[i-1]));
        }
        // This is the number of times the erased node was accessed.
        unsigned lookups=(nodes[size]->child(idx[size])).lookups();
        nodes[size]->erase(idx[size]);
        // Before returning, we have to update the upstream lookups.
        for(size_t i=0;i<size;++i)
                nodes[i]->decrease_lookup_counter(lookups);
        return;
}

// This function counts the number of times the data pointed to by idx was
// accessed.
template <class _Data>
unsigned
Trie<_Data>::use_count(const typename Trie<_Data>::key_type &idx)const{
        std::size_t size=idx.size()-1;
        node_type *n=const_cast<node_type*>(&root);
        for(std::size_t i=1;i<idx.size();++i){
                if(n->count(idx[i-1])==0)
                        return 0; // idx is not associated in the trie
                else
                        n=&(n->child(idx[i-1]));
        }
        return (n->child(idx[size])).lookups();
}

// This function returns true iff the trie is empty.
template <class _Data>
bool
Trie<_Data>::empty()const{
        return root.empty();
}

// This function removes all elements from the trie.
template <class _Data>
void
Trie<_Data>::clear(){
        root.clear();
}

} // namespace Trie
