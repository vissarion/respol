#ifndef HASHED_DETERMINANT_H
#define HASHED_DETERMINANT_H

#include <vector>
#include <numeric>
#include <boost/unordered_map.hpp>

// Test equality of two vectors.
template <class _vec>
struct eqvec:std::binary_function<_vec,_vec,bool>{
        bool operator()(const _vec v1,const _vec v2)const{
                return std::equal(v1.begin(),v1.end(),v2.begin());
        }
};

// This gives a perfect hashing function when d is at least the maximum
// element of the vector.
template <class _vec,int d>
struct hashvec:std::unary_function<_vec,size_t>{
        size_t operator()(const _vec v)const{
                size_t hash=0,count=0;
                for(typename _vec::const_iterator i=v.begin();i!=v.end();++i)
                        hash+=(*i)*d^(count++);
                return hash;
        }
};

// HashedDeterminant constructs a big matrix of columns and provides
// methods to compute and store determinants of matrices formed by columns
// of this matrix. It takes two template parameters: _NT is the number type
// of the matrix elements and _dimensions is the number of rows of the
// matrix (which is equal to dimension of the points). The amount of
// columns is determined at construction time.

template <class _NT,int d>
class HashedDeterminant{
        private:
        typedef _NT                                     NT;
        typedef std::vector<NT>                         Column;
        typedef std::vector<Column>                     Matrix;
        typedef std::vector<size_t>                     Index;
        typedef boost::unordered_map<Index,
                                     NT,
                                     hashvec<Index,d>,
                                     eqvec<Index> >     Determinants;

        public:
        HashedDeterminant(int columns):
        _points(columns,Column(d)),_determinants(){};

        template <class Iterator>
        HashedDeterminant(Iterator begin,Iterator end):
        _points(),_determinants(){
                for(Iterator i=begin;i!=end;++i)
                        _points.push_back(*i);
        }

        void set_column(int i,Column c){
                _points[i]=c;
        }

        NT& determinant(Index idx){
                if(_determinants.count(idx)==0){
                        _determinants[idx]=compute_determinant(idx);
                }
                return _determinants[idx];
        }

        // This function is a copy/paste of the TOPCOM determinant.
        NT compute_determinant(Index idx)const{
                Matrix tmp;
                for(Index::iterator ii=idx.begin();ii!=idx.end();++ii)
                        tmp.push_back(_points[*ii]);
                NT scale(1);
                for(size_t i=0;i<d-1;++i){
                        if(tmp[i][i]==NT(0)){
                                for(size_t k=i+1;k<d;++k){
                                        if(tmp[i][k]!=NT(0)){
                                                Column tmpcol=tmp[i];
                                                tmp[i]=tmp[k];
                                                tmp[k]=tmpcol;
                                                scale*=-1;
                                                continue;
                                        }
                                }
                                if(tmp[i][i]==NT(0)){
                                        return NT(0);
                                }
                        }
                        NT eraser(tmp[i][i]);
                        for(size_t j=i+1;j<d;++j){
                                NT delinquent(tmp[i][j]);
                                if(delinquent==NT(0)){
                                        continue;
                                }
                                for(size_t k=i+1;k<d;++k){
                                        tmp[k][j]-=(tmp[k][i]*
                                                    delinquent/eraser);
                                }
                                tmp[i][j]=NT(0);
                        }
                }
                NT result(1);
                for(size_t i=0;i<d;++i){
                        result*=tmp[i][i];
                }
                return result*scale;
        }

        private:
        Matrix _points;
        Determinants _determinants;
};

#endif // HASHED_DETERMINANT_H
