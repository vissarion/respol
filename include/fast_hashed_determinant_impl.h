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

#ifndef FAST_HASHED_DETERMINANT_IMPL_H
#define FAST_HASHED_DETERMINANT_IMPL_H

#ifdef USE_LINBOX_DET
// this function is extremely inefficient for small matrices, since
// there is a big overhead in constructing the LinBox matrix
template <class _NT>
inline
typename FastHashedDeterminant<_NT>::NT
FastHashedDeterminant<_NT>::compute_determinant(
                const typename FastHashedDeterminant<_NT>::Index &idx)const{
        typedef CGAL::Linbox_rational_field<NT>         Field;
        typedef CGAL::Linbox_dense_matrix<Field>        LBMatrix;
        size_t d=_points[0].size();
        LBMatrix M((int)d,(int)d);
        for(size_t row=0;row<d;++row)
                for(size_t column=0;column<d;++column)
                        M.setEntry(row,
                                   column,
                                   _points[idx[column]][row]);
        NT det(0);
        LinBox::det(det,
                    M,
                    LinBox::RingCategories::RationalTag(),
                    LinBox::Method::Hybrid());
        return det;
}

#elif USE_CGAL_DET

template <class _NT>
inline
typename FastHashedDeterminant<_NT>::NT
FastHashedDeterminant<_NT>::compute_determinant(
                const typename FastHashedDeterminant<_NT>::Index &idx)const{
        int d = idx.size();
        //std::cout << first-last << "|" << d << std::endl;
        typename LA::Matrix M(d);
        //std::vector<CPoint_d>::iterator s = first;
        for( int j = 0; j < d; ++j ){
                //std::cout << *s << std::endl;
                for( int i = 0; i < d; ++i ){
                        //std::cout << i << "," << j;
                        M(i,j) = _points[idx[i]][j];
                        //std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        return LA::determinant(M);
}

#elif USE_CGAL_DET_2

template <class _NT>
inline
typename FastHashedDeterminant<_NT>::NT
FastHashedDeterminant<_NT>::compute_determinant(
                const typename FastHashedDeterminant<_NT>::Index &idx)const{
        int d = idx.size();
        //std::cout << first-last << "|" << d << std::endl;
        typename LA::Matrix A(d/2);
        typename LA::Matrix B(d/2);
        typename LA::Matrix C(d/2);
        typename LA::Matrix D(d/2);
        //std::vector<CPoint_d>::iterator s = first;

        for( int j = d/2; j < d; ++j ){
                //std::cout << *s << std::endl;
                for( int i = 0; i < d/2; ++i ){
                        std::cout << i << "," << j;
                        A(i,j-d/2) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        for( int j = d/2; j < d; ++j ){
                //std::cout << *s << std::endl;
                for( int i = d/2; i < d; ++i ){
                        std::cout << i << "," << j;
                        B(i-d/2,j-d/2) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        for( int j = 0; j < d/2; ++j ){
                //std::cout << *s << std::endl;
                for( int i = 0; i < d/2; ++i ){
                        std::cout << i << "," << j;
                        C(i,j) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        for( int j = 0; j < d/2; ++j ){
                //std::cout << *s << std::endl;
                for( int i = d/2; i < d; ++i ){
                        std::cout << i << "," << j;
                        A(i-d/2,j) = _points[idx[i]][j];
                        std::cout << " -> " << M(i,j) << std::endl;
                }
        }
        exit(0);
        return LA::determinant(M);
}

#elif USE_EIGEN_DET

#include <eigen3/Eigen/Eigen>

template <class _NT>
inline
typename FastHashedDeterminant<_NT>::NT
FastHashedDeterminant<_NT>::eigendet(
                typename FastHashedDeterminant<_NT>::NT **m,
                const typename FastHashedDeterminant<_NT>::Index &idx3)const{
        size_t n=idx3.size();
        Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> mat(n,n);
        for(size_t i=0;i<n;++i)
                mat.col(i)=Eigen::Map<
                                        Eigen::Matrix<NT,
                                        Eigen::Dynamic,
                                        1> >(m[idx3[i]],n);
        return mat.determinant();
}

template <class _NT>
inline
typename FastHashedDeterminant<_NT>::NT
FastHashedDeterminant<_NT>::compute_determinant(
                const typename FastHashedDeterminant<_NT>::Index &idx){
        size_t n=idx.size();
        Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> mat(n,n);
        for(size_t i=0;i<n;++i)
                mat.col(i)=Eigen::Map<
                                        Eigen::Matrix<NT,
                                                      Eigen::Dynamic,
                                                      1> >
                                (&_points[idx[i]][0],n);
        return mat.determinant();
}

#else // USE_LINBOX_DET USE_CGAL_DET USE_CGAL_DET_2 USE_EIGEN_DET not defined

template <class _NT>
inline
typename FastHashedDeterminant<_NT>::NT
FastHashedDeterminant<_NT>::compute_determinant(
                const typename FastHashedDeterminant<_NT>::Index &idx)
#ifndef USE_HASHED_DETERMINANTS
const
#endif
{
        Index idx2;
        size_t n=idx.size();
        for(size_t i=1;i<n;++i)
                idx2.push_back(idx[i]);
        assert(idx2.size()==idx.size()-1);
        NT det(0);
        for(size_t i=0;i<n;++i){
                if(_points[idx[i]][n-1]!=0){
                        if((i+n)%2)
                                det+=(_points[idx[i]][n-1]*
                                      determinant(idx2));
                        else
                                det-=(_points[idx[i]][n-1]*
                                      determinant(idx2));
                }
                // update the index array
                idx2[i]=idx[i];
        }
        return det;
}
#endif // USE_LINBOX_DET USE_CGAL_DET USE_CGAL_DET_2 USE_EIGEN_DET

#endif // FAST_HASHED_DETERMINANT_IMPL_H
// vim: ts=2
