#ifndef FAST_HASHED_DETERMINANT_H
#define FAST_HASHED_DETERMINANT_H

#include <vector>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#ifdef HASH_STATISTICS
#include <ostream>
#include <ctime>
#endif

// FastHashedDeterminant constructs a big matrix of columns and provides
// methods to compute and store determinants of matrices formed by columns
// of this matrix. It takes two template parameters: _NT is the number type
// of the matrix elements. dim is the number of rows of the matrix (which
// is equal to the dimension of the points). Determinants are stored in a
// hash table, and are never recomputed. Moreover, the computed minors of
// the submatrices are also stored in the hash table, which speeds up the
// computation of determinants of some matrices.

template <class _NT,size_t dim>
class FastHashedDeterminant{
        private:
        typedef _NT                                     NT;
        typedef std::vector<NT>                         Column;
        typedef std::vector<NT>                         Row;
        typedef std::vector<Column>                     Matrix;
        typedef std::vector<size_t>                     Index;
        typedef boost::unordered_map<Index,NT>          Determinants;

        public:
        FastHashedDeterminant(size_t columns):
#ifdef HASH_STATISTICS
        number_of_full_determinant_calls(0),
        number_of_determinant_calls(0),
        number_of_hashed_determinants(0),
        number_of_computed_determinants(0),
        determinant_time(0),
#endif
        _points(columns,Column(dim)),_determinants()
        {};

        template <class Iterator>
        FastHashedDeterminant(Iterator begin,Iterator end):
#ifdef HASH_STATISTICS
        number_of_full_determinant_calls(0),
        number_of_determinant_calls(0),
        number_of_hashed_determinants(0),
        number_of_computed_determinants(0),
        determinant_time(0),
#endif
        _points(),_determinants(){
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
        std::cerr<<"hash statistics:\nnumber of determinant calls: "<<
                number_of_determinant_calls<<
                " ("<<number_of_full_determinant_calls<<
                " full-dimension)\nnumber of hashed determinants: "<<
                number_of_hashed_determinants<<
                "\nnumber of computed determinants: "<<
                number_of_computed_determinants<<
                "\nnumber of collisions: "<<
                number_of_collisions<<
                " (in "<<bad_buckets<<" buckets)\ndeterminant time: "<<
                (double)determinant_time/CLOCKS_PER_SEC<<
                " seconds"<<std::endl;
#endif
        }

        // This function sets a column of the matrix. This function must be
        // called before any determinant computation, since it will
        // invalidate the determinants which are hashed. The column i must
        // exist in the matrix.
        void set_column(size_t i,const Column &c){
                _points[i]=c;
        }

        // Push back a column, at the end of the matrix; returns the index
        // of this inserted element. In contrast with set_column, this
        // function will not invalidate hashed values.
        size_t add_column(const Column &c){
                _points.push_back(c);
                return _points.size()-1;
        }

        // This function returns the determinant of a submatrix of _points.
        // This submatrix is formed by the columns whose indices are in
        // idx. If this determinant was already computed (i.e., it is in
        // the hash table), it is returned. Otherwise, the private function
        // compute_determinant is called. The size of idx must be dim.
        NT& determinant(const Index &idx){
                if(idx.size()==1)
                        return _points[idx[0]][0];
#ifdef HASH_STATISTICS
                ++number_of_determinant_calls;
                clock_t start;
                if(idx.size()==dim){
                        start=clock();
                        determinant_time+=(clock()-start);
                        ++number_of_full_determinant_calls;
                }
#endif
                if(_determinants.count(idx)==0){
#ifdef HASH_STATISTICS
                        ++number_of_computed_determinants;
#endif
                        _determinants[idx]=compute_determinant(idx);
#ifdef HASH_STATISTICS
                }else{
                        ++number_of_hashed_determinants;
#endif
                }
#ifdef HASH_STATISTICS
                if(idx.size()==dim)
                        determinant_time+=(clock()-start);
#endif
                return _determinants[idx];
        }

        // This computes a determinant of size dim+1, where the matrix
        // _points is enlarged by adding at the bottom the vector r. The
        // size of idx must be dim+1.
        NT determinant(const Index &idx,const Row &r){
                Index idx2;
                size_t n=dim+1;
                for(size_t i=0;i<n;++i)
                        idx2.push_back(idx[i+1]);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det+=(r[idx[i]]*determinant(idx2));
                        else
                                det-=(r[idx[i]]*determinant(idx2));
                        // update the index array
                        idx2[i]=idx[i];
                }
                return det;
        }

        // This function is the same of determinant(idx,r), where r is a
        // vector full of ones. The size of idx must be dim+1.
        NT homogeneous_determinant(const Index &idx){
                Index idx2;
                size_t n=dim+1;
                for(size_t i=0;i<n;++i)
                        idx2.push_back(idx[i+1]);
                NT det(0);
                for(size_t i=0;i<n;++i){
                        if((i+n)%2)
                                det+=determinant(idx2);
                        else
                                det-=determinant(idx2);
                        // update the index array
                        idx2[i]=idx[i];
                }
                return det;
        }

        // This function prints the matrix _points to an output stream.
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
        std::ostream& print_submatrix(const Index &idx,std::ostream &o)const{
                for(size_t row=0;row<dim;++row){
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
        // of the columns which will form the submatrix. This function is
        // private, since it is not supposed to be called from outside the
        // class. Inlining this function is very important for efficiency
        // reasons!
        inline NT compute_determinant(const Index &idx){
                Index idx2;
                size_t n=idx.size();
                for(size_t i=1;i<n;++i)
                        idx2.push_back(idx[i]);
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

        private:
        Matrix _points;
        Determinants _determinants;
#ifdef HASH_STATISTICS
        public:
        unsigned number_of_full_determinant_calls;
        unsigned number_of_determinant_calls;
        unsigned number_of_hashed_determinants;
        unsigned number_of_computed_determinants;
        clock_t determinant_time;
#endif
};

#endif // FAST_HASHED_DETERMINANT_H
