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

#ifndef HASHED_DETERMINANT_LINBOX_H
#define HASHED_DETERMINANT_LINBOX_H

#include "hashed_determinant.h"
#include <CGAL/LinBox/mpq_class_field.h>
#include <CGAL/LinBox/dense_matrix.h>
#include <linbox/solutions/det.h>

template <class _NT>
class HashedDeterminantLinbox:public HashedDeterminant<_NT>{
        private:
        typedef _NT                                     NT;
        typedef HashedDeterminant<NT>                   HD;
        typedef typename HD::Index                      Index;
        public:
        HashedDeterminantLinbox():HD(){};
        HashedDeterminantLinbox(size_t columns):HD(columns){};
        template <class Iterator>
        HashedDeterminantLinbox(Iterator begin,Iterator end):HD(begin,end){}
        const char* algorithm()const{return "Linbox\0";};
        inline NT compute_determinant(const Index&);
};

// This function is extremely inefficient for small matrices, since there
// is a big overhead in constructing the LinBox matrix.
template <class _NT>
inline
typename HashedDeterminantLinbox<_NT>::NT
HashedDeterminantLinbox<_NT>::compute_determinant(
                const typename HashedDeterminantLinbox<_NT>::Index &idx){
        typedef CGAL::Linbox_rational_field<NT>         Field;
        typedef CGAL::Linbox_dense_matrix<Field>        LBMatrix;
        size_t d=this->_points[0].size();
        LBMatrix M((int)d,(int)d);
        for(size_t row=0;row<d;++row)
                for(size_t column=0;column<d;++column)
                        M.setEntry(row,
                                   column,
                                   this->_points[idx[column]][row]);
        NT det(0);
        LinBox::det(det,
                    M,
                    LinBox::RingCategories::RationalTag(),
                    LinBox::Method::Hybrid());
        return det;
}

#endif // HASHED_DETERMINANT_LINBOX_H
