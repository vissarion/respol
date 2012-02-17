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

#ifndef HASHED_DETERMINANT_BASE_H
#define HASHED_DETERMINANT_BASE_H

// Check Boost version. We use boost::functional and boost::unordered,
// which appeared in versions 1.16.0 and 1.36.0, respectively. We give an
// error if Boost version is older than 1.36.0.
#include <boost/version.hpp>
#if BOOST_VERSION < 103600
#error This needs Boost 1.36 or newer
#endif

#ifndef USE_HASHED_DETERMINANTS
#undef HASH_STATISTICS
#endif

#include <vector>
#include <set>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <cassert>

#ifdef HASH_STATISTICS
#include <ostream>
#include <ctime>
#endif

#ifdef USE_SORTED_INDICES
#include "sort_swap.h"
#endif

#ifdef USE_ONLY_CAYLEY_DET_HASH
template <class _NT> class HashedDeterminantBase;
#include <CGAL/Cartesian_d.h>
#endif

// HashedDeterminantBase constructs a big matrix of columns and provides
// functions to compute and store determinants of matrices formed by
// columns of this matrix. It takes one template parameter, _NT, which is
// the number type of the matrix elements.Determinants are stored in a hash
// table, and are never recomputed. The member function compute_determinant
// implements the actual determinant algorithm, that is Laplace expansion
// by default. By overriding this function, a new class can implement
// hashed determinants with other determinant algorithm.  In the case of
// Laplace expansion, since it is recursive, the computed minors of the
// submatrices are also stored in the hash table, which speeds up the
// computation of determinants of some matrices.

template <class _NT>
class HashedDeterminantBase{
        public:
        typedef _NT                                     NT;
#ifdef USE_ONLY_CAYLEY_DET_HASH
        typedef CGAL::Cartesian_d<NT>                   CK;
        typedef typename CK::LA                         LA;
#endif

        typedef std::vector<NT>                         Column;
        typedef std::vector<NT>                         Row;
        typedef std::vector<size_t>                     Index;
        private:
        typedef std::vector<Column>                     Matrix;
        public:
        typedef typename Matrix::iterator               iterator;
        private:
        typedef boost::unordered_map<Index,NT>          Determinants;
        typedef boost::unordered_map<Index,NT>          HDeterminants;
#ifdef USE_ORIENTATION_DET
        typedef boost::unordered_map<Index,NT>          ODeterminants;
        typedef std::set<Index>                         CDeterminants;
#endif // USE_ORIENTATION_DET
#ifdef USE_SORTED_INDICES
        // define the pair type to store a sorted index and a swap boolean
        typedef std::pair<Index,bool>                   SS;
#endif // USE_SORTED_INDICES

        public:
        // constructor for incremental det table
        // (for the space of the projection of the Resultant i.e. PD)
        HashedDeterminantBase():
#ifdef LOG_DET_TIME
        determinant_time(0),
#endif
#ifdef HASH_STATISTICS
        dimension(0),
        number_of_max_dimension_calls(0),
        number_of_determinant_calls(0),
        number_of_computed_determinants(0),
        number_of_hom_determinants(0),
        number_of_computed_hom_determinants(0),
        full_determinant_time(0),
#ifdef USE_ORIENTATION_DET
        matrix_time(0),
#endif
#endif
        _points(),
        _determinants(),
        _h_determinants(),
#ifdef USE_ORIENTATION_DET
        _o_determinants(),
        _c_determinants(),
#endif
#ifdef USE_ONLY_CAYLEY_DET_HASH
        _hashed(false),
#endif
        number_of_hashed_determinants(0)
        {};

        //
        HashedDeterminantBase(size_t columns):
#ifdef LOG_DET_TIME
        determinant_time(0),
#endif
#ifdef HASH_STATISTICS
        dimension(0),
        number_of_max_dimension_calls(0),
        number_of_determinant_calls(0),
        number_of_computed_determinants(0),
        number_of_hom_determinants(0),
        number_of_computed_hom_determinants(0),
        full_determinant_time(0),
#ifdef USE_ORIENTATION_DET
        matrix_time(0),
#endif
#endif
        _points(columns),
        _determinants(),
        _h_determinants(),
#ifdef USE_ORIENTATION_DET
        _o_determinants(),
        _c_determinants(),
#endif
        number_of_hashed_determinants(0)
        {};

        // constructor for static det table
        // (for the space of the lifted Cayley pointset i.e. CD)
        template <class Iterator>
        HashedDeterminantBase(Iterator begin,Iterator end):
#ifdef LOG_DET_TIME
        determinant_time(0),
#endif
#ifdef HASH_STATISTICS
        dimension(0),
        number_of_max_dimension_calls(0),
        number_of_determinant_calls(0),
        number_of_computed_determinants(0),
        number_of_hom_determinants(0),
        number_of_computed_hom_determinants(0),
        full_determinant_time(0),
#ifdef USE_ORIENTATION_DET
        matrix_time(0),
#endif
#endif
        _points(),
        _determinants(),
        _h_determinants(),

#ifdef USE_ONLY_CAYLEY_DET_HASH
        _hashed(true),
#endif
#ifdef USE_ORIENTATION_DET
        _o_determinants(),
        _c_determinants(),
#endif
        number_of_hashed_determinants(0){
                for(Iterator i=begin;i!=end;++i){
#ifdef HASH_STATISTICS
                        if(i->size()>dimension)
                                dimension=i->size();
#endif
                        _points.push_back(*i);
                }
        }

        ~HashedDeterminantBase();

        // Returns the total time spent in determinant computations.
        double get_determinant_time();

        // This function sets a column of the matrix. This function must be
        // called before any determinant computation, since it will
        // invalidate the determinants which are hashed. The column i must
        // exist in the matrix.
        void set_column(size_t,const Column&);

        // Push back a column, at the end of the matrix; returns the index
        // of this inserted element. In contrast with set_column, this
        // function will not invalidate hashed values.
        size_t add_column(const Column&);

        // a linear search to Matrix to see if the Column c exists
        int find(const Column&)const;

        // A string describing the determinant algorithm.
        const char* algorithm()const{return "Default\0";};

        // This function prints the full matrix to an output stream.
        std::ostream& print_matrix(std::ostream&)const;

        // This function prints a square submatrix, formed by the first
        // elements of the columns whose indices are in idx, to an output
        // stream.
        std::ostream& print_submatrix(const Index&,std::ostream&)const;

        // The following functions deal with the columns of the points
        // table.
        Column& operator[](size_t);
        typename Matrix::iterator begin();
        typename Matrix::iterator end();

        // This functions returns the number of points stored in the
        // determinant table.
        int size();

        // This function returns the determinant of a submatrix of _points.
        // This submatrix is formed by the columns whose indices are in
        // idx. If this determinant was already computed (i.e., it is in
        // the hash table), it is returned. Otherwise, the private function
        // compute_determinant is called.
        private:
#ifdef USE_HASHED_DETERMINANTS
        NT& determinant(const Index&);
#else
        #ifdef LOG_DET_TIME
        NT determinant(const Index &idx);
        #else
        NT determinant(const Index &idx)const;
        #endif
#endif

        // This function computes the determinant of a submatrix, enlarged
        // with row at the bottom full of ones.
        public:
        NT homogeneous_determinant(const Index&);

        // This function computes the determinant of a submatrix, enlarged
        // with two rows at the bottom: one lifting row and a row full of
        // ones.
        public:
        NT homogeneous_determinant(const Index&,const Row&);

#ifdef USE_ORIENTATION_DET
        private:
        NT det_minor(NT **m,const Index &idx3,const Index &idx4){
                // idx3 contains the n indices of the columns of the matrix
                // m whose determinant is to be computed. idx4 contains the
                // indices of _points that correspond to the columns used
                // for the computation of m. The size of idx4 is n+1. The
                // (n+1)th element is the number of column subtracted from
                // the first n columns.
                size_t n=idx3.size();
                assert(n>0);
                assert(idx4.size()==n+1);
                if(n==1)
                        return m[idx3[0]][0];
                // check whether idx4 is in _o_determinants
                if(_o_determinants.count(idx4)!=0){
                        assert(_o_determinants.count(idx4)==1);
                        return _o_determinants[idx4];
                }
                // idx5 contains the n-1 indices of the columns which
                // determine the minor whose determinant compute. It is
                // updated on each iteration. idx6 contains n indices,
                // which indicate which columns of _points were used to
                // construct the matrix represented by idx5.
                Index idx5,idx6;
                for(size_t i=1;i<n;++i){
                        idx5.push_back(idx3[i]);
                        idx6.push_back(idx4[i]);
                }
                assert(idx5.size()==n-1);
                idx6.push_back(idx4[n]);
                assert(idx6.size()==n);
                NT det(0);
                for(size_t col=0;col<n;++col){
                        if(m[idx3[col]][n-1]!=0){
                                if((col+n)%2)
                                        det+=m[idx3[col]][n-1]*
                                             det_minor(m,idx5,idx6);
                                else
                                        det-=m[idx3[col]][n-1]*
                                             det_minor(m,idx5,idx6);
                        }
                        idx5[col]=idx3[col];
                        idx6[col]=idx4[col];
                }
                _o_determinants[idx4]=det;
                return det;
        }

        private:
        NT orientation(const Index &idx,const Row &r){
#ifdef HASH_STATISTICS
                clock_t start=clock();
                ++number_of_determinant_calls;
#endif // HASH_STATISTICS
                assert(idx.size()==r.size());
                assert(_points.size()>=idx.size());
                // n is the dimension, which is one less than the size of
                // the orientation matrix
                size_t n=idx.size()-1;
                // idx3 contains the indices of the n-1 columns of m that
                // determine the minor whose determinant is to be computed.
                // idx4 contains the indices of _points that correspond to
                // the columns used for the computation of the minor of m.
                // The size of idx4 is n. The n-th element is the number of
                // column subtracted from the first (n-1) columns.
                Index idx3,idx4;
                for(size_t i=1;i<n;++i){
                        idx3.push_back(i);
                        idx4.push_back(idx[i]);
                }
                assert(idx3.size()==n-1);
                idx4.push_back(idx[n]);
                assert(idx4.size()==n);
                // the determinant will be stored in det
                NT det(0);
                // Check first whether this determinant was computed. This
                // means that the determinants of its minors are cached,
                // thus we just need to perform only d subtractions, d
                // multiplications and d additions.
                if(_c_determinants.count(idx)!=0){
                        assert(_c_determinants.count(idx)==1);
                        NT element;
                        for(size_t col=0;col<n;++col){
                                assert(_o_determinants.count(idx4)==1);
                                element=r[col]-r[n];
                                if(element!=0){
                                        if((col+n)%2)
                                                det+=element*
                                                     _o_determinants[idx4];
                                        else
                                                det-=element*
                                                     _o_determinants[idx4];
                                }
                                idx4[col]=idx[col];
                        }
                        return det;
                }
#ifdef HASH_STATISTICS
                ++number_of_computed_determinants;
                clock_t matrix_start=clock();
#endif
                _c_determinants.insert(idx);
                assert(_c_determinants.count(idx)>0);
                // m is a n*n matrix, which is obtained from the original
                // (n+1)*(n+1) matrix by subtracting the (n+1)th column to
                // the first n columns
                NT **m=(NT**)malloc(n*sizeof(NT*));
                for(size_t col=0;col<n;++col){
                        m[col]=new NT[n];
                        for(size_t row=0;row<n-1;++row)
                                m[col][row]=_points[idx[col]][row]-
                                            _points[idx[n]][row];
                        // TODO: it is superfluous to store this row
                        m[col][n-1]=r[col]-r[n];
                }
#ifdef HASH_STATISTICS
                matrix_time+=(clock()-matrix_start);
#endif
                for(size_t col=0;col<n;++col){
                        NT minor=det_minor(m,idx3,idx4);
                        if(m[col][n-1]!=0){
                                if((col+n)%2)
                                        det+=m[col][n-1]*minor;
                                else
                                        det-=m[col][n-1]*minor;
                        }
                        idx3[col]=col;
                        idx4[col]=idx[col];
                }
                for(size_t col=0;col<n;++col)
                        delete[] m[col];
                free(m);
#ifdef HASH_STATISTICS
                full_determinant_time+=(clock()-start);
#endif // HASH_STATISTICS
                return det;
        }
#endif // USE_ORIENTATION_DET

        // This function computes the determinant of a submatrix of
        // _points. The parameter idx is a vector of indices of the indices
        // of the columns which will form the submatrix. Inlining this
        // function is very important for efficiency reasons!
        private:
#ifndef USE_HASHED_DETERMINANTS
        inline NT compute_determinant(const Index&)const;
#else
        inline NT compute_determinant(const Index&);
#endif

        protected:
        Matrix _points;
        private:
        Determinants _determinants;
        HDeterminants _h_determinants;
#ifdef USE_ONLY_CAYLEY_DET_HASH
				bool _hashed;
#endif // USE_ONLY_CAYLEY_DET_HASH
#ifdef USE_ORIENTATION_DET
        ODeterminants _o_determinants;
        CDeterminants _c_determinants;
#endif // USE_ORIENTATION_DET
        unsigned long number_of_hashed_determinants;
#ifdef LOG_DET_TIME
        clock_t determinant_time;
#endif // LOG_DET_TIME
#ifdef HASH_STATISTICS
        public:
        size_t dimension;
        unsigned long number_of_max_dimension_calls;
        unsigned long number_of_determinant_calls;
        unsigned long number_of_computed_determinants;
        unsigned long number_of_hom_determinants;
        unsigned long number_of_computed_hom_determinants;
        clock_t full_determinant_time;
#ifdef USE_ORIENTATION_DET
        clock_t matrix_time;
#endif // USE_ORIENTATION_DET
#endif // HASH_STATISTICS
};

#include "hashed_determinant_base_impl.h"

#endif // HASHED_DETERMINANT_BASE_H
// vim: ts=2
