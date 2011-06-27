#ifndef FAST_HASHED_DETERMINANT_H
#define FAST_HASHED_DETERMINANT_H

#include <vector>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <cassert>
#ifdef HASH_STATISTICS
#include <ostream>
#include <ctime>
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
        typedef std::vector<NT>                         Column;
        typedef std::vector<NT>                         Row;
        typedef std::vector<size_t>                     Index;
        private:
        typedef std::vector<Column>                     Matrix;
        typedef boost::unordered_map<Index,NT>          Determinants;
        typedef boost::unordered_map<Index,NT>          HDeterminants;

        public:
        // constructor for incremental det table
        // (for the space of the projection of the Resultant i.e. PD)
        FastHashedDeterminant():
#ifdef HASH_STATISTICS
        number_of_determinant_calls(0),
        number_of_computed_determinants(0),
        number_of_hom_determinants(0),
        number_of_computed_hom_determinants(0),
        full_determinant_time(0),
#endif
        _points(),
        _determinants(),
        _homogeneous_determinants()
        {};
        
        //
        FastHashedDeterminant(size_t columns):
#ifdef HASH_STATISTICS
        number_of_determinant_calls(0),
        number_of_computed_determinants(0),
        number_of_hom_determinants(0),
        number_of_computed_hom_determinants(0),
        full_determinant_time(0),
#endif
        _points(columns),
        _determinants(),
        _homogeneous_determinants()
        {};
				
				// constructor for static det table 
        // (for the space of the lifted Cayley pointset i.e. CD)
        template <class Iterator>
        FastHashedDeterminant(Iterator begin,Iterator end):
#ifdef HASH_STATISTICS
        number_of_determinant_calls(0),
        number_of_computed_determinants(0),
        number_of_hom_determinants(0),
        number_of_computed_hom_determinants(0),
        full_determinant_time(0),
#endif
        _points(),
        _determinants(),
        _homogeneous_determinants(){
                for(Iterator i=begin;i!=end;++i)
                        _points.push_back(*i);
        }

        ~FastHashedDeterminant(){
#ifdef HASH_STATISTICS
        size_t number_of_collisions=0,bad_buckets=0;
        for(size_t bucket=0;bucket!=_determinants.bucket_count();++bucket)
                if(_determinants.bucket_size(bucket)>1){
                        number_of_collisions+=
                                (_determinants.bucket_size(bucket)-1);
                        ++bad_buckets;
                }
        std::cerr<<
                "hash statistics:"<<
                "\nnum of hashed points:"<<
                _points.size()<<
                "\nnon-hom hash: "<<
                _determinants.bucket_count()<<" buckets, "<<
                number_of_collisions<<" collisions in "<<bad_buckets<<
                " buckets\nnon-hom determinants: computed "<<
                number_of_computed_determinants<<" out of "<<
                number_of_determinant_calls<<
                "\ntime in non-hom determinant computations: "<<
                (double)full_determinant_time/CLOCKS_PER_SEC<<
                " seconds\nhom determinants: computed "<<
                number_of_computed_hom_determinants<<
                " out of "<<number_of_hom_determinants<<std::endl;
#endif
        }

        // This function sets a column of the matrix. This function must be
        // called before any determinant computation, since it will
        // invalidate the determinants which are hashed. The column i must
        // exist in the matrix.
        void set_column(size_t i,const Column &c){
                assert((i>=0)&&(i<_points.size()));
                _points[i]=c;
        }

        // Push back a column, at the end of the matrix; returns the index
        // of this inserted element. In contrast with set_column, this
        // function will not invalidate hashed values.
        size_t add_column(const Column &c){
                _points.push_back(c);
                return _points.size()-1;
        }
				
				// a linear search to Matrix to see if the Column c exists
				int find(const Column &c){
								typename Matrix::iterator result;
								result = std::find(_points.begin(),_points.end(),c);
								if (result == _points.end())
									return -1;
								else 
									return result - _points.begin();
				}
        // This function returns the determinant of a submatrix of _points.
        // This submatrix is formed by the columns whose indices are in
        // idx. If this determinant was already computed (i.e., it is in
        // the hash table), it is returned. Otherwise, the private function
        // compute_determinant is called.
        NT& determinant(const Index &idx){
                if(idx.size()==1)
                        return _points[idx[0]][0];
#ifdef HASH_STATISTICS
                ++number_of_determinant_calls;
#endif
                if(_determinants.count(idx)==0){
#ifdef HASH_STATISTICS
                        ++number_of_computed_determinants;
                        clock_t start;
                        start=clock();
#endif
                        _determinants[idx]=compute_determinant(idx);
#ifdef HASH_STATISTICS
                        full_determinant_time+=(clock()-start);
#endif
                }
                return _determinants[idx];
        }

        // This function computes the determinant of a submatrix, enlarged
        // with row at the bottom full of ones.
        NT homogeneous_determinant(const Index &idx){
#ifdef HASH_STATISTICS
                number_of_hom_determinants+=1;
#endif
                if(_homogeneous_determinants.count(idx)!=0){
                        assert(_homogeneous_determinants.count(idx)==1);
                        return _homogeneous_determinants[idx];
                }
#ifdef HASH_STATISTICS
                number_of_computed_hom_determinants+=1;
#endif
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i)
                        idx2.push_back(idx[i]);
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
                _homogeneous_determinants[idx]=det;
                return det;
        }

        // This function is the same homogeneous_determinant(idx), but
        // adding the vector r between the submatrix and the vector of
        // ones.
        NT old_homogeneous_determinant(const Index &idx,const Row &r){
                assert(idx.size()==r.size());
                Index idx2,idxr;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i){
                        idx2.push_back(idx[i]);
                        idxr.push_back(i);
                }
                assert(idx2.size()==n-1);
                assert(idxr.size()==n-1);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det+=enlarged_homogeneous_determinant(idx2,
                                                                      r,
                                                                      idxr);
                        else
                                det-=enlarged_homogeneous_determinant(idx2,
                                                                      r,
                                                                      idxr);
                        // update the index array
                        idx2[i]=idx[i];
                        idxr[i]=i;
                }
                return det;
        }

        NT homogeneous_determinant(const Index &idx,const Row &r){
                assert(idx.size()==r.size());
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i)
                        idx2.push_back(idx[i]);
                assert(idx2.size()==n-1);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det-=r[i]*homogeneous_determinant(idx2);
                        else
                                det+=r[i]*homogeneous_determinant(idx2);
                        // update the index array
                        idx2[i]=idx[i];
                }
                return det;
        }

        // This function prints the full matrix to an output stream.
        std::ostream& print_matrix(std::ostream &o)const{
                for(size_t i=0;i<_points.size();++i){
                        o<<"[ ";
                        for(size_t j=0;j<_points[i].size();++j)
                                o<<_points[i][j]<<' ';
                        o<<"]\n";
                }
                return o;
        }

        // This function prints a submatrix, formed by the columns whose
        // indices are in idx, to an output stream.
        /*std::ostream& print_submatrix(const Index &idx,std::ostream &o)const{
                for(size_t row=0;row<dim;++row){
                        o<<"[ ";
                        for(Index::const_iterator i=idx.begin();
                            i!=idx.end();
                            ++i)
                                o<<_points[row][*i]<<' ';
                        o<<"]\n";
                }
                return o;
        }*/

        private:
        // This function computes the determinant of a submatrix of
        // _points. The parameter idx is a vector of indices of the indices
        // of the columns which will form the submatrix. Inlining this
        // function is very important for efficiency reasons!
        inline NT compute_determinant(const Index &idx){
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i)
                        idx2.push_back(idx[i]);
                assert(idx2.size()==idx.size()-1);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det+=(_points[idx[i]][n-1]*determinant(idx2));
                        else
                                det-=(_points[idx[i]][n-1]*determinant(idx2));
                        // update the index array
                        idx2[i]=idx[i];
                }
                return det;
        }

        // This function computes a determinant of size dim+1, where the
        // matrix is enlarged by adding at the bottom the vector r. There
        // is also a vector of indices of r, idxr, which specifies the
        // elements of r used to compute the determinant. The size of idx
        // and idxr must be dim+1.
        NT enlarged_homogeneous_determinant(const Index &idx,
                                            const Row &r,
                                            const Index &idxr){
                assert(idx.size()==idxr.size());
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i)
                        idx2.push_back(idx[i]);
                assert(idx2.size()==n-1);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det+=(r[idxr[i]]*determinant(idx2));
                        else
                                det-=(r[idxr[i]]*determinant(idx2));
                        // update the index array
                        idx2[i]=idx[i];
                }
                return det;
        }

        private:
        Matrix _points;
        Determinants _determinants;
        HDeterminants _homogeneous_determinants;
#ifdef HASH_STATISTICS
        public:
        unsigned number_of_determinant_calls;
        unsigned number_of_computed_determinants;
        unsigned number_of_hom_determinants;
        unsigned number_of_computed_hom_determinants;
        clock_t full_determinant_time;
#endif
};

#endif // FAST_HASHED_DETERMINANT_H
