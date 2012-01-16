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

// constructor
template <class _Data>
Node<_Data>::Node():_num_lookups(0){}

template <class _Data>
typename Node<_Data>::data_type& Node<_Data>::get(){
        ++_num_lookups;
        return _d;
}

template <class _Data>
unsigned Node<_Data>::lookups()const{
        return _num_lookups;
}

template <class _Data>
void Node<_Data>::set_lookup_counter(unsigned c){
        _num_lookups=c;
        return;
}

template <class _Data>
void Node<_Data>::decrease_lookup_counter(unsigned c){
        _num_lookups-=c;
        return;
}

template <class _Data>
typename Node<_Data>::Node& Node<_Data>::child(std::size_t s){
        ++_num_lookups;
        return (*this)[s];
}

template <class _Data>
std::size_t Node<_Data>::children()const{
        return this->size();
}

template <class _Data>
std::ostream& Node<_Data>::print_map(std::ostream &o)const{
        o<<'{';
        for(typename base::const_iterator i=this->begin();i!=this->end();++i){
                if(i!=this->begin())
                        o<<',';
                o<<'('<<(i->first)<<','<<&(i->second)<<')';
        }
        return o<<'}';
}

} // namespace Trie
