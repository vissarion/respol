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

template <class _NT>
FastHashedDeterminant<_NT>::~FastHashedDeterminant(){
#ifdef HASH_STATISTICS
#ifdef USE_ORIENTATION_DET
        std::cerr<<"hash statistics:\n"<<
                _points.size()<<" hashed points, of max dim "<<dimension<<
                '\n'<<number_of_computed_determinants<<
                " determinants computed, out of "<<
                number_of_determinant_calls<<
                '\n'<<number_of_computed_determinants<<
                " matrices computed in "<<
                (double)matrix_time/CLOCKS_PER_SEC<<" seconds\n"<<
                (double)full_determinant_time/CLOCKS_PER_SEC<<
                " seconds spent in all determinant computations"<<
                std::endl;
#else // USE_ORIENTATION_DET is not defined
        // collect non-homogeneous hash statistics
        size_t nh_collisions=0,nh_bad_buckets=0,nh_biggest_bucket=0;
        for(size_t bucket=0;bucket!=_determinants.bucket_count();++bucket)
                if(_determinants.bucket_size(bucket)>1){
                        nh_collisions+=(_determinants.bucket_size(bucket)-1);
                        if(_determinants.bucket_size(bucket)>nh_biggest_bucket)
                                nh_biggest_bucket=
                                        _determinants.bucket_size(bucket);
                        ++nh_bad_buckets;
                }
        // collect homogeneous hash statistics
        size_t h_collisions=0,h_bad_buckets=0,h_biggest_bucket=0;
        for(size_t bucket=0;
            bucket!=_h_determinants.bucket_count();
            ++bucket)
                if(_h_determinants.bucket_size(bucket)>1){
                        h_collisions+=(_h_determinants.bucket_size(bucket)-1);
                        if(_h_determinants.bucket_size(bucket)>h_biggest_bucket)
                                h_biggest_bucket=
                                        _h_determinants.bucket_size(bucket);
                        ++h_bad_buckets;
                }
        std::cerr<<
                "hash statistics:\n"<<
                _points.size()<<" hashed points, of max dim "<<dimension<<
                "\nnon-hom hash: "<<
                _determinants.bucket_count()<<" buckets, "<<
                nh_collisions<<" collisions (" <<
                ((double)nh_collisions/
                 (double)number_of_computed_determinants)*100<<
                "%) in "<<nh_bad_buckets<<
                " buckets of max size "<<nh_biggest_bucket<<
                "\nnon-hom determinants: computed "<<
                number_of_computed_determinants<<" out of "<<
                number_of_determinant_calls<<
                "\n load_factor()="<<_determinants.load_factor()<<
                ",max_load_factor()="<<_determinants.max_load_factor()<<
                "\ntime in non-hom full-dim determinant computations: "<<
                (double)full_determinant_time/CLOCKS_PER_SEC<<
                " seconds\nhom hash: "<<
                _h_determinants.bucket_count()<<" buckets, "<<
                h_collisions<<" collisions ("<<
                ((double)h_collisions/
                 (double)number_of_computed_hom_determinants)*100<<
                "%) in "<<h_bad_buckets<<
                " buckets of max size "<<h_biggest_bucket<<
                "\nhom determinants: computed "<<
                number_of_computed_hom_determinants<<
                " out of "<<number_of_hom_determinants<<
                "\n load_factor()="<<_h_determinants.load_factor()<<
                ",max_load_factor()="<<_h_determinants.max_load_factor()<<
                std::endl;
#endif // USE_ORIENTATION_DET
#endif // HASH_STATISTICS
        }

template <class _NT>
double FastHashedDeterminant<_NT>::get_determinant_time(){
#if defined(LOG_DET_TIME) && defined(HASH_STATISTICS)
#ifdef USE_ORIENTATION_DET
        return (double)full_determinant_time/CLOCKS_PER_SEC;
#else // USE_ORIENTATION_DET is not defined
        return (double)determinant_time/CLOCKS_PER_SEC;
#endif // USE_ORIENTATION_DET
#else // (defined(LOG_DET_TIME) && defined(HASH_STATISTICS)) is false
        return -1.;
#endif // defined(LOG_DET_TIME) && defined(HASH_STATISTICS)
}

template <class _NT>
void FastHashedDeterminant<_NT>::set_column(
                size_t i,
                const FastHashedDeterminant<_NT>::Column &c){
        assert((i>=0)&&(i<_points.size()));
#ifdef HASH_STATISTICS
        if(c.size()>dimension)
                dimension=c.size();
#endif
        _points[i]=c;
}

template <class _NT>
size_t FastHashedDeterminant<_NT>::add_column(
                const FastHashedDeterminant<_NT>::Column &c){
#ifdef HASH_STATISTICS
        if(c.size()>dimension)
                dimension=c.size();
#endif
        _points.push_back(c);
        return _points.size()-1;
}

template <class _NT>
int FastHashedDeterminant<_NT>::find(
                const FastHashedDeterminant<_NT>::Column &c)const{
        typename Matrix::const_iterator result;
        result=std::find(_points.begin(),_points.end(),c);
        return (result==_points.end()?-1:result-_points.begin());
}

template <class _NT>
std::ostream& FastHashedDeterminant<_NT>::print_matrix(std::ostream &o)const{
        if(!_points.size())
                return o;
        size_t s=_points[0].size();
        for(size_t row=0;row<s;++row){
                o<<"[ ";
                for(size_t i=0;i<_points.size();++i)
                        o<<_points[i][row]<<' ';
                o<<"]\n";
        }
        return o;
}

template <class _NT>
std::ostream& FastHashedDeterminant<_NT>::print_submatrix(
                const FastHashedDeterminant<_NT>::Index &idx,
                std::ostream &o)const{
        for(size_t row=0;row<idx.size();++row){
                o<<"[ ";
                for(Index::const_iterator i=idx.begin();
                    i!=idx.end();
                    ++i)
                        o<<_points[row][*i]<<' ';
                o<<"]\n";
        }
        return o;
}

template <class _NT>
typename FastHashedDeterminant<_NT>::Column&
FastHashedDeterminant<_NT>::operator[](size_t i){
        CGAL_assertion(i>=0&&i<_points.size());
        return _points[i];
}

template <class _NT>
typename FastHashedDeterminant<_NT>::Matrix::iterator
FastHashedDeterminant<_NT>::begin(){
        return _points.begin();
}

template <class _NT>
typename FastHashedDeterminant<_NT>::Matrix::iterator
FastHashedDeterminant<_NT>::end(){
        return  _points.end();
}

template <class _NT>
int FastHashedDeterminant<_NT>::size(){
        return  _points.size();
}

#ifdef USE_LINBOX_DET

#include <CGAL/LinBox/mpq_class_field.h>
#include <CGAL/LinBox/dense_matrix.h>
#include <linbox/solutions/det.h>

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
