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

#ifndef HASHED_DETERMINANT_EIGEN_H
#define HASHED_DETERMINANT_EIGEN_H

#include "hashed_determinant.h"
#include <Eigen/Eigen>

template <class _NT>
class HashedDeterminantEigen:public HashedDeterminant<_NT>{
        private:
        typedef _NT                                     NT;
        typedef HashedDeterminant<NT>                   HD;
        typedef typename HD::Index                      Index;
        public:
        HashedDeterminantEigen():HD(){};
        HashedDeterminantEigen(size_t columns):HD(columns){};
        template <class Iterator>
        HashedDeterminantEigen(Iterator begin,Iterator end):HD(begin,end){}
        const char* algorithm()const{return "Eigen\0";};
        inline NT compute_determinant(const Index&);
};

template <class _NT>
inline
typename HashedDeterminantEigen<_NT>::NT
HashedDeterminantEigen<_NT>::compute_determinant(
                const typename HashedDeterminantEigen<_NT>::Index &idx){
        size_t n=idx.size();
        Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> mat(n,n);
        for(size_t i=0;i<n;++i)
                mat.col(i)=Eigen::Map<
                                        Eigen::Matrix<NT,
                                                      Eigen::Dynamic,
                                                      1> >
                                (&this->_points[idx[i]][0],n);
        return mat.determinant();
}

#endif // HASHED_DETERMINANT_EIGEN_H
