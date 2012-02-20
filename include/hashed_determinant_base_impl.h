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

#ifndef HASHED_DETERMINANT_BASE_IMPL_H
#define HASHED_DETERMINANT_BASE_IMPL_H

#ifdef USE_SORTED_INDICES
#include "sort_swap.h"
#endif

template <class _NT>
HashedDeterminantBase<_NT>::~HashedDeterminantBase(){
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
                (nh_collisions==0?
                 0.:
                 ((double)nh_collisions/
                  (double)number_of_computed_determinants)*100)<<
                "%) in "<<nh_bad_buckets<<
                " buckets of max size "<<nh_biggest_bucket<<
                "\nnon-hom determinants: computed "<<
                number_of_computed_determinants<<" out of "<<
                number_of_determinant_calls<<
                "\nload_factor()="<<_determinants.load_factor()<<
                ", max_load_factor()="<<_determinants.max_load_factor()<<
                "\ntime in non-hom full-dim determinant computations: "<<
                (double)full_determinant_time/CLOCKS_PER_SEC<<
                " seconds\nhom hash: "<<
                _h_determinants.bucket_count()<<" buckets, "<<
                h_collisions<<" collisions ("<<
                (h_collisions==0?
                 0.:
                 ((double)h_collisions/
                 (double)number_of_computed_hom_determinants)*100)<<
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
double HashedDeterminantBase<_NT>::get_determinant_time(){
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
void HashedDeterminantBase<_NT>::set_column(
                size_t i,
                const HashedDeterminantBase<_NT>::Column &c){
        assert((i>=0)&&(i<_points.size()));
#ifdef HASH_STATISTICS
        if(c.size()>dimension)
                dimension=c.size();
#endif
        _points[i]=c;
}

template <class _NT>
size_t HashedDeterminantBase<_NT>::add_column(
                const HashedDeterminantBase<_NT>::Column &c){
#ifdef HASH_STATISTICS
        if(c.size()>dimension)
                dimension=c.size();
#endif
        _points.push_back(c);
        return _points.size()-1;
}

template <class _NT>
int HashedDeterminantBase<_NT>::find(
                const HashedDeterminantBase<_NT>::Column &c)const{
        typename Matrix::const_iterator result;
        result=std::find(_points.begin(),_points.end(),c);
        return (result==_points.end()?-1:result-_points.begin());
}

template <class _NT>
std::ostream& HashedDeterminantBase<_NT>::print_matrix(std::ostream &o)const{
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
std::ostream& HashedDeterminantBase<_NT>::print_submatrix(
                const HashedDeterminantBase<_NT>::Index &idx,
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
typename HashedDeterminantBase<_NT>::Column&
HashedDeterminantBase<_NT>::operator[](size_t i){
        assert(i>=0&&i<_points.size());
        return _points[i];
}

template <class _NT>
typename HashedDeterminantBase<_NT>::Matrix::iterator
HashedDeterminantBase<_NT>::begin(){
        return _points.begin();
}

template <class _NT>
typename HashedDeterminantBase<_NT>::Matrix::iterator
HashedDeterminantBase<_NT>::end(){
        return  _points.end();
}

template <class _NT>
int HashedDeterminantBase<_NT>::size(){
        return  _points.size();
}

template <class _NT>
typename HashedDeterminantBase<_NT>::NT
#ifdef USE_HASHED_DETERMINANTS
& HashedDeterminantBase<_NT>::determinant(
                const typename HashedDeterminantBase<_NT>::Index &idx){
#ifdef USE_SORTED_INDICES
        assert(boost::is_sorted(idx));
#endif // USE_SORTED_INDICES
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
        if(_determinants.count(idx)==0){
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
                return (_determinants[idx]=compute_determinant(idx));
#endif // HASH_STATISTICS
                _determinants[idx]=compute_determinant(idx);
#ifdef HASH_STATISTICS
                if(idx.size()==dimension)
                        full_determinant_time+=(clock()-start);
#endif // HASH_STATISTICS
        }
#ifdef LOG_DET_TIME
        if(idx.size()==_points[0].size())
                determinant_time+=clock()-start_all;
#endif // LOG_DET_TIME
        return _determinants[idx];
}
#else // USE_HASHED_DETERMINANTS is not defined
HashedDeterminantBase<_NT>::determinant(
                const typename HashedDeterminantBase<_NT>::Index &idx)
#ifndef LOG_DET_TIME
        const
#endif
{
        if(idx.size()==1)
                return _points[idx[0]][0];
#ifdef LOG_DET_TIME
        clock_t start_all;
        if(idx.size()==_points[0].size())
                start_all=clock();
        NT det(compute_determinant(idx));
        if(idx.size()==_points[0].size())
                determinant_time+=clock()-start_all;
        return det;
#else // LOG_DET_TIME
        return compute_determinant(idx);
#endif // LOG_DET_TIME
}
#endif // USE_HASHED_DETERMINANTS

template <class _NT>
typename HashedDeterminantBase<_NT>::NT
HashedDeterminantBase<_NT>::homogeneous_determinant(
                const typename HashedDeterminantBase<_NT>::Index &idx,
                bool idx_is_sorted){
#ifdef HASH_STATISTICS
        number_of_hom_determinants+=1;
#endif // HASH_STATISTICS
#ifdef LOG_DET_TIME
        clock_t start_all,det_old;
        start_all=clock();
        det_old=determinant_time;
#endif // LOG_DET_TIME

#ifdef USE_SORTED_INDICES
        boost::tuple<Index,bool> idxs;
        if(idx_is_sorted){
                assert(boost::is_sorted(idx));
                idxs=boost::make_tuple(boost::ref(idx),true);
        }else{
                idxs=sort_swap(idx);
        }
#ifdef USE_HASHED_DETERMINANTS
        if(_h_determinants.count(boost::get<0>(idxs))!=0){
                assert(_h_determinants.count(boost::get<0>(idxs))==1);
                return (boost::get<1>(idxs)?
                        _h_determinants[boost::get<0>(idxs)]:
                        -_h_determinants[boost::get<0>(idxs)]);
        }
#endif // USE_HASHED_DETERMINANTS
#ifdef HASH_STATISTICS
        number_of_computed_hom_determinants+=1;
#endif // HASH_STATISTICS
        // Compute the determinant of boost::get<0>(idxs).
        Index idx2;
        size_t n=idx.size();
        for(size_t i=1;i<n;++i)
                idx2.push_back(boost::get<0>(idxs)[i]);
        assert(idx2.size()==n-1);
        NT det(0);
        for(size_t i=0;i<n;++i){
                if((i+n)%2)
                        det+=determinant(idx2);
                else
                        det-=determinant(idx2);
                // update the index array
                idx2[i]=boost::get<0>(idxs)[i];
        }
#if defined USE_HASHED_DETERMINANTS
                _h_determinants[boost::get<0>(idxs)]=det;
        #endif // USE_HASHED_DETERMINANTS
        #ifdef LOG_DET_TIME
                determinant_time=det_old+(clock()-start_all);
        #endif // LOG_DET_TIME
                return (boost::get<1>(idxs)?det:-det);

#else // USE_SORTED_INDICES

#ifdef USE_HASHED_DETERMINANTS
        if(_h_determinants.count(idx)!=0){
                assert(_h_determinants.count(idx)==1);
                return _h_determinants[idx];
        }
#endif // defined USE_HASHED_DETERMINANTS
#ifdef HASH_STATISTICS
        number_of_computed_hom_determinants+=1;
#endif // HASH_STATISTICS
        // Compute the determinant of idx.
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
#endif // USE_SORTED_INDICES
}

template <class _NT>
typename HashedDeterminantBase<_NT>::NT
HashedDeterminantBase<_NT>::homogeneous_determinant(
                const typename HashedDeterminantBase<_NT>::Index &idx,
                const typename HashedDeterminantBase<_NT>::Row &r){
#ifdef USE_ONLY_CAYLEY_DET_HASH
  if(_hashed){
#endif
        assert(idx.size()==r.size());
#ifdef LOG_DET_TIME
        clock_t start_all,det_old;
        start_all=clock();
        det_old=determinant_time;
#endif

#ifdef USE_SORTED_INDICES
        boost::tuple<Index,bool,Index> ssp=sort_swap_permutation(idx);
        Index idx2;
        size_t n=idx.size();
        for(size_t i=1;i<n;++i)
                idx2.push_back(boost::get<0>(ssp)[i]);
        assert(idx2.size()==n-1);
        assert(boost::is_sorted(idx2));
        NT det(0);
        for(size_t i=0;i<n;++i){
                if(r[boost::get<2>(ssp)[i]]!=0){
                        if((i+n)%2)
                                det-=r[boost::get<2>(ssp)[i]]*
                                        homogeneous_determinant(idx2,true);
                        else
                                det+=r[boost::get<2>(ssp)[i]]*
                                        homogeneous_determinant(idx2,true);
                }
                // update the index array
                idx2[i]=boost::get<0>(ssp)[i];
                assert(boost::is_sorted(idx2));
        }
#else // USE_SORTED_INDICES
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
#endif // USE_SORTED_INDICES

#ifdef LOG_DET_TIME
        determinant_time=det_old+(clock()-start_all);
#endif
#ifdef USE_SORTED_INDICES
        return (boost::get<1>(ssp)?det:-det);
#else // USE_SORTED_INDICES
        return det;
#endif // USE_SORTED_INDICES
#ifdef USE_ONLY_CAYLEY_DET_HASH
  } else {
        int d=idx.size()-1;
        typename LA::Matrix M(d);
        for(int j=0;j<d-1;++j)
                for(int i=1;i<=d;++i)
                        M(i-1,j)=_points[idx[i]][j]-_points[idx[0]][j];
        for(int i=1;i<=d;++i)
                M(i-1,d-1)=r[i]-r[0];
        return LA::sign_of_determinant(M);
  }
#endif // USE_ONLY_CAYLEY_DET_HASH
}

template <class _NT>
inline
typename HashedDeterminantBase<_NT>::NT
HashedDeterminantBase<_NT>::compute_determinant(
                const typename HashedDeterminantBase<_NT>::Index &idx)
#ifndef USE_HASHED_DETERMINANTS
const
#endif
{
#ifdef USE_SORTED_INDICES
        assert(boost::is_sorted(idx));
#endif
        size_t n=idx.size();
        Index idx2;
        for(size_t i=1;i<n;++i)
                idx2.push_back(idx[i]);
        assert(idx2.size()==idx.size()-1);
        NT det(0);
        for(size_t i=0;i<n;++i){
                if(_points[idx[i]][n-1]!=0){
                        if((i+n)%2)
                                det+=(_points[idx[i]][n-1]*determinant(idx2));
                        else
                                det-=(_points[idx[i]][n-1]*determinant(idx2));
                }
                // update the index array
                idx2[i]=idx[i];
        }
        return det;
}

#endif // HASHED_DETERMINANT_BASE_IMPL_H
// vim: ts=2
