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

#include <boost/version.hpp>
#if BOOST_VERSION >= 104800
#include <boost/container/vector.hpp>
#define svector boost::container::vector
#else
#include <vector>
#define svector std::vector
#endif

#ifndef BIRD_H
#define BIRD_H

namespace Bird{

template <class NT>
std::ostream& print_matrix(std::ostream &o,const svector<svector<NT> > &M){
        o<<'[';
        size_t rows=M.size();
        if(rows>0){
                size_t cols=M[0].size();
                for(size_t r=0;r<rows;++r){
                        o<<"[ ";
                        for(size_t c=0;c<cols;++c)
                                o<<M[r][c]<<' ';
                        o<<']';
                }
        }
        return o<<']';
}

template <class NT>
std::ostream& print_vector(std::ostream &o,const svector<NT> &V){
        o<<"[ ";
        size_t d=V.size();
        for(size_t i=0;i<d;++i)
                o<<V[i]<<' ';
        return o<<']';
}

// Given matrices X and A, this function computes F_A(X)=mu(X)*A and stores
// it in the matrix muXA. The input vector muv contains the diagonal of
// mu(X). At exit, the vector muv is updated with the diagonal of muXA.
template <class NT>
svector<svector<NT> > fax(const svector<svector<NT> > &X,
                          const svector<svector<NT> > &A,
                          svector<NT> &muv){
        size_t rows=X.size();
        size_t cols=A.size();
        svector<NT> thisrow;
        thisrow.reserve(cols);
        svector<svector<NT> > muXA;
        muXA.reserve(rows-1);
        for(size_t r=0;r+1<rows;++r){
                thisrow.clear();
                for(size_t c=0;c<cols;++c){
                        if(r+2!=rows||c==r){
                                // Computation of thisrow[c].
                                thisrow.push_back(muv[rows-r-2]*A[r][c]);
                                for(size_t i=r+1;i<cols;++i)
                                        thisrow[c]+=X[r][i]*A[i][c];
                        }else{
                                // This is the last row. The only element
                                // of this row we'll use later is the
                                // diagonal element (which will be used to
                                // compute the mu-vector). Thus, we set all
                                // the other elements to zero.
                                thisrow.push_back(0);
                        }
                }
                muXA.push_back(thisrow);
        }
        // Compute the mu-vector of mu(X).A.
        muv.clear();
        size_t s=muXA.size()-1;
        muv.push_back(-muXA[s][s]);
        for(size_t i=s-1;i>0;--i)
                muv.push_back(muv.back()-muXA[i][i]);
        return muXA;
}

// Given a matrix A, this function computes F_A^t(A) and stores it in FAnA.
template <class NT>
void fata(size_t t,
          const svector<svector<NT> > &A,
          svector<svector<NT> > &FAnA){
        svector<NT> muv;
        size_t d=A.size()-1;
        // Compute the mu-vector of A.
        muv.reserve(d);
        muv.push_back(-A[d][d]);
        for(size_t i=d-1;i>0;--i)
                muv.push_back(muv.back()-A[i][i]);
        // i=1: FAnA=mu(A)*A
        FAnA=fax(A,A,muv);
        // The vector muv contains now the mu-vector of FAnA.
        // i=2..n-1
        for(size_t i=1;i<t;++i){
                // FAnA=mu(FAnA)*A
                FAnA=fax(FAnA,A,muv);
        }
        return;
}

// This function computes the determinant of a given matrix A. For that, it
// applies n-2 times the function F_A to it. The last application is done
// by multiplying the first row of the resulting matrix by the first column
// of A, because only the first element is needed.
template <class NT>
NT determinant(const svector<svector<NT> > &A){
        size_t n=A.size();
        svector<svector<NT> > Fnm2; // This matrix will store F_A^{n-2}(A).
        fata(n-2,A,Fnm2);
        // F_A was applied n-2 times to A. It remains to compute the
        // determinant, which is just a vector-vector product. The first
        // scalar product comes from the mu-vector of Fnm2 times the first
        // element of the first row of A.
        if(n%2){
                // If n is odd, the determinant equals the vector-vector
                // product.
                NT det(-Fnm2[1][1]*A[0][0]);
                for(size_t i=1;i<n;++i)
                        det+=Fnm2[0][i]*A[i][0];
                return det;
        }else{
                // If n is even, the determinant equals the additive
                // inverse of the vector-vector product.
                NT det(Fnm2[1][1]*A[0][0]);
                for(size_t i=1;i<n;++i)
                        det-=Fnm2[0][i]*A[i][0];
                return det;
        }
}

} // namespace Bird

#endif // BIRD_H
