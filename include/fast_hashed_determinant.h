#ifndef FAST_HASHED_DETERMINANT_H
#define FAST_HASHED_DETERMINANT_H

#include <ostream>
#include <vector>
#include <numeric>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#ifdef HASH_STATISTICS
#include <ctime>
#endif

// FastHashedDeterminant constructs a big matrix of columns and provides
// methods to compute and store determinants of matrices formed by columns
// of this matrix. It takes two template parameters: _NT is the number type
// of the matrix elements. _dimensions is the number of rows of the matrix
// (which is equal to dimension of the points). Determinants are stored in
// a hash table, and are never recomputed. Moreover, the computed minors of
// the submatrices are also stored in the hash table, which speeds up the
// computation of determinants of some matrices.

template <class _NT,size_t d>
class FastHashedDeterminant{
        private:
        typedef _NT                                     NT;
        typedef std::vector<NT>                         Column;
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
        _points(columns,Column(d)),_determinants()
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

        void set_column(size_t i,Column c){
                _points[i]=c;
        }

        NT& determinant(const Index &idx){
                if(idx.size()==1)
                        return _points[idx[0]][0];
#ifdef HASH_STATISTICS
                ++number_of_determinant_calls;
                if(idx.size()==d)
                        ++number_of_full_determinant_calls;
#endif
                if(_determinants.count(idx)==0){
#ifdef HASH_STATISTICS
                        ++number_of_computed_determinants;
                        clock_t start=clock();
#endif
                        _determinants[idx]=compute_determinant(idx);
#ifdef HASH_STATISTICS
                        determinant_time+=(clock()-start);
                }else{
                        ++number_of_hashed_determinants;
#endif
                }
                return _determinants[idx];
        }

        // inlining this function is very important for efficiency reasons
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

        std::ostream& print_matrix(std::ostream &o)const{
                for(size_t i=0;i<_points.size();++i){
                        o<<"[ ";
                        for(size_t j=0;j<_points[i].size();++j)
                                o<<_points[i][j]<<' ';
                        o<<"]\n";
                }
                return o;
        }

        std::ostream& print_submatrix(const Index &idx,std::ostream &o)const{
                for(size_t row=0;row<d;++row){
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
