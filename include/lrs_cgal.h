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
std::set<std::vector<NT_> >
lrs_ch(const std::vector<std::vector<NT_> > &pts){
        typedef NT_                                     Integer;
        typedef std::vector<Integer>                    Point;
        typedef typename std::vector<Point>             Points;
        typedef std::set<std::vector<Integer> >         Normals;
        // print input
        /*std::cout<<"input points:\n";
        for(typename Points::const_iterator i=p.begin();i!=p.end();++i){
                std::cout<<"[ ";
                for(typename Point::const_iterator j=i->begin();
                    j!=i->end();
                    ++j){
                        std::cout<<(*j)<<' ';
                }
                std::cout<<"]\n";
        }*/
        lrs_dic *P;
        lrs_dat *Q;
        lrs_mp_vector output;
        lrs_mp_matrix Lin;
        long i;
        long col;

        if(!lrs_init((char*)"cgal_lrs\n"))
                exit(-1);
        Q=lrs_alloc_dat((char*)"LRS globals");
        if(Q==NULL)
                exit(-2);
        Q->m=pts.size(); // number of input rows = number of vertices
        Q->n=pts[0].size()+1; // number of input columns = dimension + 1
        Q->hull = TRUE; // convex hull problem: facet enumeration
        Q->polytope = FALSE;// input is not a polytope
        Q->getvolume= FALSE; // don't compute the volume
        output=lrs_alloc_mp_vector(Q->n);
        P=lrs_alloc_dic(Q); // allocate and initialize lrs_dic
        if(P==NULL)
                exit(-3);
        insert_data(pts,P,Q);
        if(!lrs_getfirstbasis(&P,Q,&Lin,FALSE))
                exit(-4);
        for(col=0L;col<Q->nredundcol;col++) // print linearity space
                lrs_printoutput(Q,Lin[col]); // Array Lin[][] holds the coeffs
        do{
                for(col=0;col<=P->d;col++)
                        if(lrs_getsolution(P,Q,output,col))
                                lrs_printoutput(Q,output);
        }
        while(lrs_getnextbasis(&P,Q,FALSE));
                lrs_printtotals(P,Q); // print final totals
        // free space : do not change order of next 3 lines!
        lrs_clear_mp_vector(output,Q->n);
        lrs_free_dic(P,Q); // deallocate lrs_dic
        lrs_free_dat(Q); // deallocate lrs_dat
        lrs_close((char*)"cgal_lrs\n");
        // TODO: get results back
        Normals nv;
        return nv;
}

template <class NT_>
void insert_data(const std::vector<std::vector<NT_> > &pts,
                 lrs_dic *P,
                 lrs_dat *Q){
        long m=Q->m;
        long n=Q->n;
        for(size_t row=0;row<m;++row){ 
                lrs_mp_vector num=lrs_alloc_mp_vector(pts[row].size());
                lrs_mp_vector den=lrs_alloc_mp_vector(pts[row].size());
                for(size_t j=0;j<pts[row].size();j++){
                        mpz_set(num[j],pts[row][j].mpz());
                        mpz_set_ui(den[j],1);
                }
                lrs_set_row_mp(P,Q,row+1,num,den,GE);
                lrs_clear_mp_vector(num,pts[row].size());
                lrs_clear_mp_vector(den,pts[row].size());
        }
}
