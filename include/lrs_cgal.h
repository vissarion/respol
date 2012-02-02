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

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

#include <gmp.h>
extern "C"{
#include <lrs_functions.h>
}

template <class NT_>
class LRS_CH{
        // A polytope is represented as a set of halfspaces. A halfspace
        // a_0*x_0+a_1*x_1+...+a_(d-1)*x_(d-1)+c>=0 in dimension d is
        // represented as a vector <a_0,a_1,...,a_(d-1),c>. A (outer)
        // normal to this halfspace is the vector <-a_0,-a_1,...,-a_(d-1)>.
        private:
        typedef NT_                                     NT;
        typedef std::vector<NT>                         Point;
        typedef typename std::vector<Point>             Points;
        typedef std::vector<NT>                         Halfspace;
        typedef std::vector<NT>                         Normal;
        typedef std::set<Halfspace>                     Polytope;

        Points v_rep;
        Polytope h_rep;
        NT denominator;
        std::set<Normal> normals;
        bool h_rep_computed,normals_computed;

        public:
        LRS_CH(const Points&);
        ~LRS_CH();
        Points& get_v_rep()const;
        Polytope& get_h_rep();
        // compute_h_rep() is the main function of the class. It calls LRS,
        // and computes a V-representation of the polytope given a set of
        // points (its H-representation). The functions called depend on
        // the number type used to represent the points coordinates.
        int compute_h_rep();
        NT& get_denominator()const;
        void compute_normals();
        std::set<Normal>& get_normals();
        private:
        void construct_matrix(lrs_dic*,lrs_dat*);
        void get_halfspaces(lrs_dic*,lrs_dat*,lrs_mp_vector);
};

template <class NT_>
LRS_CH<NT_>::LRS_CH(const LRS_CH<NT_>::Points &ps){
        v_rep=ps;
        h_rep=Polytope();
        denominator=NT(0);
        normals=std::set<Normal>();
        h_rep_computed=false;
        normals_computed=false;
}

template <class NT_>
LRS_CH<NT_>::~LRS_CH(){}

template <class NT_>
typename LRS_CH<NT_>::Points& LRS_CH<NT_>::get_v_rep()const{
        return v_rep;
}

template <class NT_>
typename LRS_CH<NT_>::Polytope& LRS_CH<NT_>::get_h_rep(){
        if(!h_rep_computed)
                compute_h_rep();
        return h_rep;
}

// This function is based in LRS' chdemo.c, combining its functions main()
// and makecyclic(). Most original LRS comments are left, since they might
// be helpful in debugging.
template <class NT_>
int LRS_CH<NT_>::compute_h_rep(){
        if(h_rep_computed)
                return 0;
        lrs_dic *P; // structure for holding current dictionary and indices
        lrs_dat *Q; // structure for holding static problem data
        lrs_mp_vector output; // one line of output:ray,vertex,facet,linearity
        lrs_mp_matrix Lin; // holds input linearities if any are found

        long col; // output column index for dictionary

        // Global initialization - done once
#ifdef LRS_DEBUG
        if(!lrs_init((char*)"lrs_cgal\n"))
                return -1;
#else
        if(!lrs_init((char*)""))
                return -2;
#endif

        // allocate and init structure for static problem data
        Q=lrs_alloc_dat((char*)"LRS globals");
        if(Q==NULL)
                return -3;
        // now flags in lrs_dat can be set
        Q->m=v_rep.size(); // number of input rows = number of vertices
        Q->n=v_rep[0].size()+1; // number of input columns (dimension+1)
        Q->hull=TRUE; // yes, convex hull problem: facet enumeration
        Q->polytope=TRUE; // yes, input is a polytope -- TODO: check this
        Q->getvolume=TRUE; // yes, compute the volume

        output=lrs_alloc_mp_vector(Q->n);

        P=lrs_alloc_dic(Q); // allocate and initialize lrs_dic
        if(P==NULL)
                return -4;

        // Build polyhedron: constraints and objective

#ifdef LRS_DEBUG
        std::cout<<"\n*input: "<<Q->m<<" vertices in R^"<<Q->n-1;
#endif

        // construct the data matrix
        construct_matrix(P,Q);

        // Pivot to a starting dictionary
        if(!lrs_getfirstbasis(&P,Q,&Lin,FALSE))
                return -5;

        // There may have been column redundancy (although not for this
        // example of cyclic polytopes) If so the linearity space is
        // obtained and redundant columns are removed. User can access
        // linearity space from lrs_mp_matrix Lin dimensions nredundcol x
        // d+1

#ifdef LRS_DEBUG
        for(col=0L;col<Q->nredundcol;col++) // print linearity space
                lrs_printoutput(Q,Lin[col]); // Array Lin[][] holds the coeffs
#endif

        // put the solutions in h_rep
        get_halfspaces(P,Q,output);

#ifdef LRS_DEBUG
        lrs_printtotals(P,Q); // print final totals
#endif

        // free space: do not change order of next 3 lines!
        lrs_clear_mp_vector(output,Q->n);
        lrs_free_dic(P,Q); // deallocate lrs_dic
        lrs_free_dat(Q); // deallocate lrs_dat

#ifdef LRS_DEBUG
        lrs_close((char*)"lrs_cgal\n");
#else
        lrs_close((char*)"");
#endif
        h_rep_computed=true;
        return 0;
}

template <class NT_>
void LRS_CH<NT_>::construct_matrix(lrs_dic *P,lrs_dat *Q){
        std::cout<<
                "construct_matrix not implemented for this number type"<<
                std::endl;
        return;
}

template <>
void LRS_CH<long int>::construct_matrix(lrs_dic *P,lrs_dat *Q){
#ifdef LRS_DEBUG
        std::cout<<"\nconstruct_matrix<long int>";
#endif
        // this function is based on makecyclic(P,Q)
        long num[Q->n]; // TODO: use v_rep instead of copying to num
        long den[Q->n];
        for(size_t j=0;j<Q->n;++j)
                den[j]=1;
        for(size_t row=0;row<Q->m;++row){
                num[0]=1;
                for(size_t j=1;j<Q->n;j++)
                        num[j]=v_rep[row][j-1];
                lrs_set_row(P,Q,row+1,num,den,TRUE);
        }
        return;
}

template <>
void LRS_CH<CGAL::Gmpz>::construct_matrix(lrs_dic *P,lrs_dat *Q){
#ifdef LRS_DEBUG
        std::cout<<"\nconstruct_matrix<Gmpz>";
#endif
        long size=Q->n-1;
        // TODO: avoid copying v_rep to num
        lrs_mp_vector num=lrs_alloc_mp_vector(size);
        lrs_mp_vector den=lrs_alloc_mp_vector(size);
        for(size_t j=0;j<Q->n;++j)
                mpz_set_ui(den[j],1);
        for(size_t row=0;row<Q->m;++row){
                mpz_set_ui(num[0],1);
                for(size_t j=1;j<Q->n;j++)
                        mpz_set(num[j],v_rep[row][j-1].mpz());
                lrs_set_row_mp(P,Q,row+1,num,den,TRUE);
        }
        lrs_clear_mp_vector(num,size);
        lrs_clear_mp_vector(den,size);
        return;
}

template <>
void LRS_CH<CGAL::Gmpq>::construct_matrix(lrs_dic *P,lrs_dat *Q){
#ifdef LRS_DEBUG
        std::cout<<"\nconstruct_matrix<Gmpq>";
#endif
        long size=Q->n-1;
        // TODO: avoid copying v_rep to num and den
        lrs_mp_vector num=lrs_alloc_mp_vector(size);
        lrs_mp_vector den=lrs_alloc_mp_vector(size);
        for(size_t row=0;row<Q->m;++row){
                mpz_set_ui(num[0],1);
                mpz_set_ui(den[0],1);
                for(size_t j=1;j<Q->n;j++){
                        mpz_set(num[j],mpq_numref(v_rep[row][j-1].mpq()));
                        mpz_set(den[j],mpq_denref(v_rep[row][j-1].mpq()));
                }
                lrs_set_row_mp(P,Q,row+1,num,den,TRUE);
        }
        lrs_clear_mp_vector(num,size);
        lrs_clear_mp_vector(den,size);
        return;
}

template <class NT_>
void LRS_CH<NT_>::get_halfspaces(lrs_dic *P,lrs_dat *Q,lrs_mp_vector output){
        std::cout<<
                "get_halfspaces not implemented for this number type"<<
                std::endl;
        return;
}

template <>
void LRS_CH<long int>::get_halfspaces(lrs_dic *P,
                                       lrs_dat *Q,
                                       lrs_mp_vector output){
#ifdef LRS_DEBUG
        std::cout<<"\nget_halfspaces<long int>";
#endif
        // We initiate reverse search from this dictionary getting new
        // dictionaries until the search is complete. User can access each
        // output line from output which is vertex/ray/facet from the
        // lrs_mp_vector output.
        do{
                for(long col=0;col<=P->d;col++){
                        if(lrs_getsolution(P,Q,output,col)){
#ifdef LRS_DEBUG
                                lrs_printoutput(Q,output);
#endif
                                Halfspace h;
                                for(size_t i=1;i<Q->n;++i)
                                        h.push_back(mpz_get_si(output[i]));
                                h.push_back(mpz_get_si(output[0]));
                                h_rep.insert(h);
                        }
                }
                denominator=mpz_get_si(P->det);
        }while(lrs_getnextbasis(&P,Q,FALSE));
#ifdef LRS_DEBUG
        std::cout<<"\ndenominator="<<denominator;
#endif
        return;
}

template <>
void LRS_CH<CGAL::Gmpz>::get_halfspaces(lrs_dic *P,
                                         lrs_dat *Q,
                                         lrs_mp_vector output){
#ifdef LRS_DEBUG
        std::cout<<"\nget_halfspaces<Gmpz>";
#endif
        // We initiate reverse search from this dictionary getting new
        // dictionaries until the search is complete. User can access each
        // output line from output which is vertex/ray/facet from the
        // lrs_mp_vector output.
        do{
                for(long col=0;col<=P->d;col++){
                        if(lrs_getsolution(P,Q,output,col)){
#ifdef LRS_DEBUG
                                lrs_printoutput(Q,output);
#endif
                                Halfspace h;
                                for(size_t i=1;i<Q->n;++i)
                                        h.push_back(CGAL::Gmpz(output[i]));
                                h.push_back(CGAL::Gmpz(output[0]));
                                h_rep.insert(h);
                        }
                }
                denominator=CGAL::Gmpz(P->det);
        }while(lrs_getnextbasis(&P,Q,FALSE));
#ifdef LRS_DEBUG
        std::cout<<"\ndenominator="<<denominator;
#endif
        return;
}

template <>
void LRS_CH<CGAL::Gmpq>::get_halfspaces(lrs_dic *P,
                                         lrs_dat *Q,
                                         lrs_mp_vector output){
#ifdef LRS_DEBUG
        std::cout<<"\nget_halfspaces<Gmpq>";
#endif
        // We initiate reverse search from this dictionary getting new
        // dictionaries until the search is complete. User can access each
        // output line from output which is vertex/ray/facet from the
        // lrs_mp_vector output.
        do{
                for(long col=0;col<=P->d;col++){
                        if(lrs_getsolution(P,Q,output,col)){
#ifdef LRS_DEBUG
                                lrs_printoutput(Q,output);
#endif
                                Halfspace h;
                                for(size_t i=1;i<Q->n;++i)
                                        h.push_back(CGAL::Gmpq(output[i]));
                                h.push_back(CGAL::Gmpq(output[0]));
                                h_rep.insert(h);
                        }
                }
                denominator=CGAL::Gmpq(P->det);
        }while(lrs_getnextbasis(&P,Q,FALSE));
#ifdef LRS_DEBUG
        std::cout<<"\ndenominator="<<denominator;
#endif
        return;
}

template <class NT_>
typename LRS_CH<NT_>::NT& LRS_CH<NT_>::get_denominator()const{
        return const_cast<NT&>(denominator);
}

template <class NT_>
void LRS_CH<NT_>::compute_normals(){
        if(!h_rep_computed){
                compute_h_rep();
                compute_normals();
        }else{
                if(normals_computed)
                        return;
        }
        for(typename Polytope::const_iterator hi=h_rep.begin();
            hi!=h_rep.end();
            ++hi){
                // *hi is a halfspace
                size_t nsize=hi->size()-1;
                Normal n(nsize);
                for(size_t j=0;j<nsize;++j)
                        n[j]=-(*hi)[j];
                normals.insert(n);
        }
        normals_computed=true;
        return;
}

template <class NT_>
std::set<typename LRS_CH<NT_>::Normal>& LRS_CH<NT_>::get_normals(){
        if(!normals_computed)
                compute_normals();
        return normals;
}
