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

#include <vector>
#include <set>
#include <iostream>

#include <gmp.h>
extern "C"{
#include <lrs_functions.h>
}

template <class NT_>
class LRS_CH{
        private:
        typedef NT_                                     Integer;
        typedef std::vector<Integer>                    Point;
        typedef typename std::vector<Point>             Points;
        typedef std::vector<Integer>                    Hyperplane;
        typedef std::set<Hyperplane>                    Polytope;

        Points v_rep;
        Polytope h_rep;

        public:
        LRS_CH(const Points&);
        ~LRS_CH();
        Points& get_v_rep()const;
        Polytope& get_h_rep()const;
        // compute_h_rep() is the main function of the class. It calls LRS,
        // and computes a V-representation of the polytope given a set of
        // points (its H-representation). The functions called depend on
        // the number type used to represent the points coordinates.
        int compute_h_rep();
};

template <class NT_>
LRS_CH<NT_>::LRS_CH(const LRS_CH<NT_>::Points &ps){
        v_rep=ps;
        h_rep=Polytope();
}

template <class NT_>
LRS_CH<NT_>::~LRS_CH(){
}

template <class NT_>
typename LRS_CH<NT_>::Points& LRS_CH<NT_>::get_v_rep()const{
        return v_rep;
}

template <class NT_>
typename LRS_CH<NT_>::Polytope& LRS_CH<NT_>::get_h_rep()const{
        return const_cast<Polytope&>(h_rep);
}

template <class NT_>
int LRS_CH<NT_>::compute_h_rep(){
        std::cerr<<"not implemented for this type"<<std::endl;
        assert(false);
        return;
}

// This function is based in LRS' chdemo.c, combining its functions main()
// and makecyclic(). Most original LRS comments are left, since they might
// be helpful in debugging.
template <>
int LRS_CH<long int>::compute_h_rep(){
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
        Q->polytope=TRUE; // yes, input is a polytope
        Q->getvolume=TRUE; // yes, compute the volume

        output=lrs_alloc_mp_vector(Q->n);

        P=lrs_alloc_dic(Q); // allocate and initialize lrs_dic
        if(P==NULL)
                return -4;

        // Build polyhedron: constraints and objective

        //std::cout<<"\n\n*cyclic polytope: "<<Q->m<<" vertices in R^"<<Q->n-1;
        // construct the data matrix based on makecyclic(P,Q)
        long num[MAXCOL]; // TODO: use v_rep instead of copying to num
        long den[MAXCOL];
        for(size_t j=0;j<Q->n;++j)
                den[j]=1;
        for(size_t row=0;row<Q->m;++row){
                num[0]=1;
                for(size_t j=1;j<Q->n;j++)
                        num[j]=v_rep[row][j-1];
                lrs_set_row(P,Q,row+1,num,den,TRUE);
        }
        // end of matrix construction

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

        // We initiate reverse search from this dictionary getting new
        // dictionaries until the search is complete User can access each
        // output line from output which is vertex/ray/facet from the
        // lrs_mp_vector output
        do{
                for(col=0;col<=P->d;col++){
                        if(lrs_getsolution(P,Q,output,col)){
#ifdef LRS_DEBUG
                                lrs_printoutput(Q,output);
#endif
                                Hyperplane h;
                                for(size_t i=1;i<Q->n;++i)
                                        h.push_back(mpz_get_si(output[i]));
                                h.push_back(mpz_get_si(output[0]));
                                h_rep.insert(h);
                        }
                }
        }while(lrs_getnextbasis(&P,Q,FALSE));

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
        return 0;
}

template <>
int LRS_CH<CGAL::Gmpz>::compute_h_rep(){
        std::cerr<<"compute_h_rep<Gmpz> not implemented yet"<<std::endl;
        return 0;
}

template <>
int LRS_CH<CGAL::Gmpq>::compute_h_rep(){
        std::cerr<<"compute_h_rep<Gmpq> not implemented yet"<<std::endl;
        return 0;
}
