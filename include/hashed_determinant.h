#ifndef HASHED_DETERMINANT_H
#define HASHED_DETERMINANT_H

// Force inlining of template functions is very important to avoid
// recursion in the naive determinant computation.
#ifdef _MSC_VER
#define __FORCE_INLINE __forceinline
#elif defined(__GNUC__)
#define __FORCE_INLINE __inline__ __attribute__((always_inline))
#else
#define __FORCE_INLINE inline
#endif

#include <ostream>
#include <vector>
#include <numeric>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#ifdef USE_LINBOX_DET
#include <CGAL/LinBox/mpq_class_field.h>
#include <CGAL/LinBox/dense_matrix.h>
#include <linbox/solutions/det.h>
#endif
#ifdef HASH_STATISTICS
#include <ctime>
#endif

// This function object computes the determinant of a matrix using the
// naive algorithm. The template parameter d makes the generated code not
// recursive, since the function calls are expanded at compilation time.
template <int d,class _NT>
struct det_naive:
public std::binary_function<size_t*,
                            std::vector<std::vector<_NT> >,
                            _NT>{
        typedef _NT                                             NT;
        typedef std::vector<std::vector<NT> >                   Matrix;
        typedef det_naive<d-1,NT>                               NaiveMinorDet;
        __FORCE_INLINE NT operator()(const size_t *idx,const Matrix &m)const{
                NaiveMinorDet nd;
                NT det(0);
                size_t *idx2=(size_t*)malloc((d-1)*sizeof(size_t));
                // Compute the first index array.
                for(size_t i=1;i<d;++i)
                        idx2[i-1]=idx[i];
                for(size_t i=0;i<d-1;++i){
                        if((i+d)%2)
                                det+=(m[idx[i]][d-1]*nd(idx2,m));
                        else
                                det-=(m[idx[i]][d-1]*nd(idx2,m));
                        // Update the index array.
                        idx2[i]=idx[i];
                }
                // We do the last iteration outside the loop to avoid
                // modifying the array (small speedup).
                det+=m[idx[d-1]][d-1]*nd(idx2,m);
                free(idx2);
                return det;
        }
};

// Base cases for the previous functor, d=2 and d=1.
template <class _NT>
struct det_naive<2,_NT>:
public std::binary_function<size_t*,
                            std::vector<std::vector<_NT> >,
                            _NT>{
        typedef _NT                                             NT;
        typedef std::vector<std::vector<NT> >                   Matrix;
        __FORCE_INLINE NT operator()(const size_t *idx,const Matrix &m)const{
                return m[idx[0]][0]*m[idx[1]][1]-m[idx[0]][1]*m[idx[1]][0];
        }
};

template <class _NT>
struct det_naive<1,_NT>:
public std::binary_function<size_t*,
                            std::vector<std::vector<_NT> >,
                            _NT>{
        typedef _NT                                             NT;
        typedef std::vector<std::vector<NT> >                   Matrix;
        __FORCE_INLINE NT operator()(const size_t *idx,const Matrix &m)const{
                return m[idx[0]][0];
        }
};

// This function object computes the determinant of a matrix using TOPCOM's
// algorithm (copied from TOPCOM sources).
// TODO: this function seems to be buggy!
template <int d,class _NT>
struct det_topcom:
public std::binary_function<size_t*,
                            std::vector<std::vector<_NT> >,
                            _NT>{
        typedef _NT                                             NT;
        typedef std::vector<std::vector<NT> >                   Matrix;
        __FORCE_INLINE NT operator()(const size_t *idx,const Matrix &m)const{
                // Copy the matrix, since it can be modified.
                Matrix tmp;
                for(size_t i=0;i<d;++i)
                        tmp.push_back(m[idx[i]]);
                // In order to swap rows efficiently, we copy the pointers
                // to the columns of tmpmatrix to a pointer array.
                //NT** tmp=(NT**)malloc(d*sizeof(NT*));
                //for(size_t i=0;i<d;++i)
                //        tmp[i]=&(Column(_points[idx[i]])[0]);
                // The remaining code of this function is a copy/paste of
                // the TOPCOM determinant.
                NT scale(1);
                for(size_t i=0;i<d-1;++i){
                        if(tmp[i][i]==NT(0)){
                                for(size_t k=i+1;k<d;++k){
                                        if(tmp[i][k]!=NT(0)){
                                                // option 1: easy swap
                                                const std::vector<NT> tmpcol=
                                                        tmp[i];
                                                tmp[i]=tmp[k];
                                                tmp[k]=tmpcol;
                                                // option 2: STL swap
                                                //tmp[i].swap(tmp[k]);
                                                // option 3: pointer swap
                                                // (if works, it's faster)
                                                //NT *tmpi=tmp[i];
                                                //tmp[i]=tmp[k];
                                                //tmp[i]=tmpi;
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
};

#ifdef USE_LINBOX_DET
// This function object computes the determinant of a matrix using LinBox.
template <class _NT>
struct det_linbox:
public std::binary_function<size_t*,
                            std::vector<std::vector<_NT> >,
                            _NT>{
        typedef _NT                                             NT;
        typedef std::vector<std::vector<NT> >                   Matrix;
        __FORCE_INLINE NT operator()(const size_t *idx,const Matrix &m)const{
                typedef CGAL::Linbox_rational_field<NT>         Field;
                typedef CGAL::Linbox_dense_matrix<Field>        LBMatrix;
                // TODO: check that the constructed matrix is correct!
                size_t d=m[0].size();
                LBMatrix M((int)d,(int)d);
                for(size_t row=0;row<d;++row)
                        for(size_t column=0;column<d;++column)
                                M.setEntry(row,column,m[idx[column]][row]);
                NT det(0);
                LinBox::detin(det,M);
                return det;
        }
};
#endif

// HashedDeterminant constructs a big matrix of columns and provides
// methods to compute and store determinants of matrices formed by columns
// of this matrix. It takes three template parameters: _NT is the number
// type of the matrix elements. _dimensions is the number of rows of the
// matrix (which is equal to dimension of the points; the amount of columns
// is determined at construction time). _DetF is a function object that
// computes a determinant of a submatrix of the given one.
// TODO: describe the _DetF interface

template <class _NT,
          size_t d,
#ifdef USE_LINBOX_DET
          class _DetF=det_linbox<_NT>
#else
          //class _DetF=det_topcom<d,_NT>
          class _DetF=det_naive<d,_NT>
#endif
         >
class HashedDeterminant{
        private:
        typedef _NT                                     NT;
        typedef _DetF                                   DetF;
        typedef std::vector<NT>                         Column;
        typedef std::vector<Column>                     Matrix;
        typedef std::vector<size_t>                     Index;
        typedef boost::unordered_map<Index,NT>          Determinants;

        public:
        HashedDeterminant(size_t columns):
#ifdef HASH_STATISTICS
        number_of_determinant_calls(0),
        number_of_hashed_determinants(0),
        number_of_computed_determinants(0),
        determinant_time(0),
#endif
        _points(columns,Column(d)),_determinants(),_compute_det()
        {};

        template <class Iterator>
        HashedDeterminant(Iterator begin,Iterator end):
#ifdef HASH_STATISTICS
        number_of_determinant_calls(0),
        number_of_hashed_determinants(0),
        number_of_computed_determinants(0),
        determinant_time(0),
#endif
        _points(),_determinants(),_compute_det(){
                for(Iterator i=begin;i!=end;++i)
                        _points.push_back(*i);
        }

        ~HashedDeterminant(){
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
                "\nnumber of hashed determinants: "<<
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
#ifdef HASH_STATISTICS
                ++number_of_determinant_calls;
#endif
                if(_determinants.count(idx)==0){
#ifdef HASH_STATISTICS
                        ++number_of_computed_determinants;
                        clock_t start=clock();
#endif
                        _determinants[idx]=_compute_det(&(idx[0]),_points);
#ifdef HASH_STATISTICS
                        determinant_time+=(clock()-start);
                }else{
                        ++number_of_hashed_determinants;
#endif
                }
                return _determinants[idx];
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
        DetF _compute_det;
#ifdef HASH_STATISTICS
        public:
        unsigned number_of_determinant_calls;
        unsigned number_of_hashed_determinants;
        unsigned number_of_computed_determinants;
        clock_t determinant_time;
#endif
};

template <class Matrix>
std::ostream& print_matrix(std::ostream &o,const Matrix &m){
        for(size_t i=0;i<m.size();++i){
                o<<"[ ";
                for(size_t j=0;j<m[i].size();++j)
                        o<<m[i][j]<<' ';
                o<<"]\n";
        }
        return o;
};

#endif // HASHED_DETERMINANT_H
