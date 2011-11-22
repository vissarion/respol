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

#ifndef FAST_HASHED_DETERMINANT_H
#define FAST_HASHED_DETERMINANT_H

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
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <cassert>

#ifdef HASH_STATISTICS
#include <ostream>
#include <ctime>
#endif

#ifdef USE_LINBOX_DET
#include <CGAL/LinBox/mpq_class_field.h>
#include <CGAL/LinBox/dense_matrix.h>
#include <linbox/solutions/det.h>
#endif

#ifdef USE_CGAL_DET
#include <CGAL/Cartesian_d.h>
#endif

#ifdef USE_EIGEN_DET
#include <Eigen/Eigen>
#endif

#ifdef USE_SORTED_INDICES
#include <sort_swap.h>
#endif

// FastHashedDeterminant constructs a big matrix of columns and provides
// methods to compute and store determinants of matrices formed by columns
// of this matrix. It takes one template parameter, _NT, which is the
// number type of the matrix elements.Determinants are stored in a hash
// table, and are never recomputed. Moreover, the computed minors of the
// submatrices are also stored in the hash table, which speeds up the
// computation of determinants of some matrices.

template <class _NT>
class FastHashedDeterminant{
        public:
        typedef _NT                                     NT;
#ifdef  USE_CGAL_DET 
        typedef CGAL::Cartesian_d<NT>                   CK;
        typedef typename CK::LA                         LA;
#endif
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
        typedef typename Matrix::iterator                     iterator;
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
        FastHashedDeterminant():
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
        FastHashedDeterminant(size_t columns):
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
        FastHashedDeterminant(Iterator begin,Iterator end):
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

        ~FastHashedDeterminant(){
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

        // Returns the total time spent in determinant computations.
        double get_determinant_time(){
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

        // This function sets a column of the matrix. This function must be
        // called before any determinant computation, since it will
        // invalidate the determinants which are hashed. The column i must
        // exist in the matrix.
        void set_column(size_t i,const Column &c){
                assert((i>=0)&&(i<_points.size()));
#ifdef HASH_STATISTICS
                if(c.size()>dimension)
                        dimension=c.size();
#endif
                _points[i]=c;
        }

        // Push back a column, at the end of the matrix; returns the index
        // of this inserted element. In contrast with set_column, this
        // function will not invalidate hashed values.
        size_t add_column(const Column &c){
#ifdef HASH_STATISTICS
                if(c.size()>dimension)
                        dimension=c.size();
#endif
                _points.push_back(c);
                return _points.size()-1;
        }

        // a linear search to Matrix to see if the Column c exists
        int find(const Column &c)const{
                typename Matrix::const_iterator result;
                result=std::find(_points.begin(),_points.end(),c);
                return (result==_points.end()?-1:result-_points.begin());
        }

        // This function returns the determinant of a submatrix of _points.
        // This submatrix is formed by the columns whose indices are in
        // idx. If this determinant was already computed (i.e., it is in
        // the hash table), it is returned. Otherwise, the private function
        // compute_determinant is called.
#ifdef USE_HASHED_DETERMINANTS
#ifdef USE_SORTED_INDICES
        NT
#else // USE_SORTED_INDICES is not defined
        NT&
#endif // USE_SORTED_INDICES
        determinant(const Index &idx){
#ifdef USE_CLEAR_DET_HASH
                if(number_of_hashed_determinants==CLEAR_DET_NUMBER){
#ifdef PRINT_INFO
                        std::cout<<"CLEAR HASH!\n\n\n\n"<<std::endl;
#endif // PRINT_INFO
                        number_of_hashed_determinants=0;
                        _determinants.clear();
                }
#endif // USE_CLEAR_DET_HASH
                if(idx.size()==1)
                        return _points[idx[0]][0];
#ifdef LOG_DET_TIME
                clock_t start_all;
                if(idx.size()==_points[0].size())
                        start_all=clock();
#endif // LOG_DET_TIME
#ifdef HASH_STATISTICS
                clock_t start;
                if(idx.size()==dimension)
                        ++number_of_max_dimension_calls;
                ++number_of_determinant_calls;
#endif // HASH_STATISTICS
#ifdef USE_SORTED_INDICES
                SS idxs=sort_swap(idx);
                // compute the determinant if it is not stored in the table
                if(_determinants.count(idxs.first)==0)
#else // USE_SORTED_INDICES is not defined
                if(_determinants.count(idx)==0)
#endif // USE_SORTED_INDICES
                {
                        ++number_of_hashed_determinants;
#ifdef HASH_STATISTICS
                        ++number_of_computed_determinants;
                        if(idx.size()==dimension){
                                start=clock();
                        }
#ifdef LOG_DET_TIME
                        if(idx.size()==_points[0].size())
                                determinant_time+=clock()-start_all;
#endif // LOG_DET_TIME
#else // HASH_STATISTICS is not defined
                        return
#endif // HASH_STATISTICS
#ifdef USE_SORTED_INDICES
                        (_determinants[idxs.first]=
                                compute_determinant(idxs.first));
#else // USE_SORTED_INDICES is not defined
                        (_determinants[idx]=compute_determinant(idx));
#endif // USE_SORTED_INDICES
#ifdef HASH_STATISTICS
                        if(idx.size()==dimension)
                                full_determinant_time+=(clock()-start);
#endif // HASH_STATISTICS
                }
#ifdef LOG_DET_TIME
                if(idx.size()==_points[0].size())
                        determinant_time+=clock()-start_all;
#endif // LOG_DET_TIME
#ifdef USE_SORTED_INDICES
                assert(_determinants.count(idxs.first)>0);
                return (idxs.second?
                        _determinants[idxs.first]:
                        -_determinants[idxs.first]);
#else // USE_SORTED_INDICES is not defined
                return _determinants[idx];
#endif // USE_SORTED_INDICES
        }
#else // USE_HASHED_DETERMINANTS is not defined
        // TODO: we can use here gaussian elimination or any O(n^3)
        // factorization algorithm
        NT determinant(const Index &idx)
        #ifndef LOG_DET_TIME
        const
        #endif
        {
                size_t n=idx.size();
                if(n==1)
                        return _points[idx[0]][0];
#ifdef LOG_DET_TIME
                clock_t start_all;
                if(idx.size()==_points[0].size())
                        start_all=clock();
#endif // LOG_DET_TIME
                Index idx2;
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
#ifdef LOG_DET_TIME
                if(idx.size()==_points[0].size())
                        determinant_time+=clock()-start_all;
#endif // LOG_DET_TIME
                return det;
        }
#endif // USE_HASHED_DETERMINANTS

        // This function computes the determinant of a submatrix, enlarged
        // with row at the bottom full of ones.
        NT homogeneous_determinant(const Index &idx){
#ifdef HASH_STATISTICS
                number_of_hom_determinants+=1;
#endif // HASH_STATISTICS
#ifdef LOG_DET_TIME
                clock_t start_all,det_old;
                start_all=clock();
                det_old=determinant_time;
#endif // LOG_DET_TIME
                if(_h_determinants.count(idx)!=0){
                        assert(_h_determinants.count(idx)==1);
                        return _h_determinants[idx];
                }
#ifdef HASH_STATISTICS
                number_of_computed_hom_determinants+=1;
#endif // HASH_STATISTICS
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i) idx2.push_back(idx[i]);
                assert(idx2.size()==n-1);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det+=determinant(idx2);
                        else
                                det-=determinant(idx2);
                        // update the index array
                        idx2[i]=idx[i];
                }
#ifdef USE_HASHED_DETERMINANTS
                _h_determinants[idx]=det;
#endif // USE_HASHED_DETERMINANTS
#ifdef LOG_DET_TIME
                determinant_time=det_old+(clock()-start_all);
#endif // LOG_DET_TIME
                return det;
        }

        NT homogeneous_determinant(const Index &idx,const Row &r){
#ifdef USE_ONLY_CAYLEY_DET_HASH
          if (_hashed){      
#endif
                assert(idx.size()==r.size());
#ifdef LOG_DET_TIME
                clock_t start_all,det_old;
                start_all=clock();
                det_old=determinant_time;
#endif
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i)
                        idx2.push_back(idx[i]);
                assert(idx2.size()==n-1);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if(r[i]!=0){
                                if((i+n)%2)
                                        det-=r[i]*homogeneous_determinant(idx2);
                                else
                                        det+=r[i]*homogeneous_determinant(idx2);
                        }
                        // update the index array
                        idx2[i]=idx[i];
                }
#ifdef LOG_DET_TIME
                determinant_time=det_old+(clock()-start_all);
#endif
                return det;
#ifdef USE_ONLY_CAYLEY_DET_HASH
          } else {
						int d = idx.size() - 1;
						//std::cout << first-last << "|" << d << std::endl;
						typename LA::Matrix M(d);
						//std::vector<CPoint_d>::iterator s = first;
						for( int j = 0; j < d-1; ++j ){
										//std::cout << *s << std::endl;
										for( int i = 1; i <= d; ++i ){
														//std::cout << i << "," << j;
														M(i-1,j) = _points[idx[i]][j] - _points[idx[0]][j];
														//std::cout << " -> " << M(i,j) << std::endl;
								//						std::cout //<< _points[idx[i]][j] << " " ;
														//std::cout << "(" << _points[idx[i]][j] << "-"
														//          << _points[idx[0]][j] << ") " 
							//							          << M(i-1,j) << " ";
										}
									//	std::cout << "\n" ;
						}
						for( int i = 1; i <= d; ++i ){
										//std::cout << i << "," << j;
										M(i-1,d-1) = r[i] - r[0];
							//			std::cout //<< r[i] << " " ;
										//std::cout << "(" << r[i] << "-" << r[0] << ") =" 
								//		          << r[i] - r[0] << " ";
										//std::cout << " -> " << M(i,j) << std::endl;
						}
						//std::cout << std::endl;
						//exit(0);
						return LA::sign_of_determinant(M);
					}
#endif // USE_ONLY_CAYLEY_DET_HASH
        }

#ifdef USE_EIGEN_DET
        inline NT eigendet(NT **m,const Index &idx3){
                size_t n=idx3.size();
                Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> mat(n,n);
                for(size_t i=0;i<n;++i)
                        mat.col(i)=Eigen::Map<
                                                Eigen::Matrix<NT,
                                                Eigen::Dynamic,
                                                1> >(m[idx3[i]],n);
                return mat.determinant();
        }
#endif // USE_EIGEN_DET

#ifdef USE_ORIENTATION_DET
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
                // show matrix m
                //for(size_t row=0;row<n;++row){
                //        std::cout<<"[ ";
                //        for(size_t col=0;col<n;++col)
                //                std::cout<<m[col][row]<<' ';
                //        std::cout<<"]\n";
                //}
#ifdef USE_EIGEN_DET
                for(size_t col=0;col<n;++col){
                        NT minor;
                        if(_o_determinants.count(idx4)==0){
                                minor=eigendet(m,idx3);
                                _o_determinants[idx4]=minor;
                        }else{
                                minor=_o_determinants[idx4];
                        }
                        if(m[col][n-1]!=0){
                                if((col+n)%2)
                                        det+=m[col][n-1]*minor;
                                else
                                        det-=m[col][n-1]*minor;
                        }
                        idx3[col]=col;
                        idx4[col]=idx[col];
                }
#else // USE_EIGEN_DET is not defined
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
#endif // USE_EIGEN_DET
                for(size_t col=0;col<n;++col)
                        delete[] m[col];
                free(m);
#ifdef HASH_STATISTICS
                full_determinant_time+=(clock()-start);
#endif // HASH_STATISTICS
                return det;
        }
#endif // USE_ORIENTATION_DET

        // This function prints the full matrix to an output stream.
        std::ostream& print_matrix(std::ostream &o)const{
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

        // This function prints a square submatrix, formed by the first
        // elements of the columns whose indices are in idx, to an output
        // stream.
        std::ostream& print_submatrix(const Index &idx,std::ostream &o)const{
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

        private:
        // This function computes the determinant of a submatrix of
        // _points. The parameter idx is a vector of indices of the indices
        // of the columns which will form the submatrix. Inlining this
        // function is very important for efficiency reasons!
#ifdef USE_LINBOX_DET
        // this function is extremely inefficient for small matrices, since
        // there is a big overhead in constructing the LinBox matrix
        inline NT compute_determinant(const Index &idx)const{
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
        inline NT compute_determinant(const Index &idx)const{
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
        inline NT compute_determinant(const Index &idx)const{
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
        inline NT compute_determinant(const Index &idx){
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
        inline NT compute_determinant(const Index &idx)
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
public:
        Column& operator[](size_t i) { 
					CGAL_assertion(i>=0&&i<_points.size());
					return _points[i]; 
			  }
			  
        typename Matrix::iterator begin(){
					return  _points.begin();
			  }
			  
			  typename Matrix::iterator end(){
					return  _points.end();
			  }
			  
			  int size(){
					return  _points.size();
			  }
			  
        private:
        Matrix _points;
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

#endif // FAST_HASHED_DETERMINANT_H
// vim: ts=2
