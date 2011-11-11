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

#include <iostream>
#include <print_functions.h>

///////////////////////////////////////////////////////////
// default projections

// project at the first coordinate of each mi
std::vector<int> proj_first_coord(const int &d,
                                   int m,
                                   const std::vector<int>& mi){
        std::vector<int> proj(d);
        proj[0]=0;
        for(size_t i=1;i<d;++i){
                int mm=0;
                for(size_t j=0;j<i;++j)
                        mm+=mi[j];
                proj[i]=mm;
        }
        return proj;
}

// project at all the coordinates of the last support
std::vector<int> proj_all_from_last_support(const int &d,
                                             int &m,
                                             const std::vector<int>& mi){
        std::vector<int> proj(d);
        int mm=0;
        for (typename std::vector<int>::const_iterator mit=mi.begin();
             mit!=mi.end()-1;
             mit++)
                mm+=*mit;
        size_t j=0;
        for (int i=mm; i<mm+*(mi.end()-1); i++){
                proj[j++]=i;
        }
        std::cout << proj << std::endl;
        return proj;
}

// project at the first coordinate of each mi
std::vector<int> proj_more_coord(const int &d,
                                  int m,
                                  const std::vector<int>& mi){
        size_t nplus1 = mi.size();
        std::vector<int> proj_first = proj_first_coord(nplus1,m,mi);
        std::vector<int> proj;
        int a = (d/nplus1)-1;
        int b = d%nplus1;
        for(typename std::vector<int>::iterator pit=proj_first.begin();
            pit!=proj_first.end();
            pit++){
                proj.push_back(*pit);
                int i=0;
                for(; i<a; i++)
                        proj.push_back((*pit)+i+1);
                if (pit-proj_first.begin() < b)
                        proj.push_back((*pit)+i+1);
        }
        return proj;
}

// project at all the coordinates
std::vector<int> proj_all(const int &m){
        std::vector<int> proj;
        for(size_t i=0;i<m;++i)
                proj.push_back(i);
        return proj;
}

///////////////////////////////////////////////////////////
// input functions

template <class NT_>
int read_pointset(std::vector<std::vector<NT_> >& pointset,
                  std::vector<int>& mi,
                  std::vector<int>& proj,
                  int& m){
        typedef NT_                                             Field;
#ifdef RESTRICTED_RES
        std::cin >> restricted_num_Res;
#endif
        int d;
        std::cin >> d;
        //TODO: change them!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        D=d;
        CD = 2*D+1;             // this is the Cayley space + 1 for lifting
        if (d != D){
                std::cout <<
                "Not matching dimensions of input and compiled program!" <<
                std::endl;
                exit(-1);
        }
        m=0;
        for(size_t i=0;i<d+1;++i){
                int mi_temp;
                std::cin >> mi_temp;
                m += mi_temp;
                mi.push_back(mi_temp);
        }
        if (mi.size() != d+1){
                std::cout << "mi.size() != d+1. The number of polynomials must me one more than the dimension!" << std::endl;
                exit(-1);
        }

        char temp;
        std::string num;
        int numInt;
        bool start=false;
        // ignore the first blanks and search for '|'
        // or '\n'
        if (temp != ' '){
                while(temp != '|' && temp != '\n'){
                        temp = std::cin.get();
                }
        }
        // the projection is given by the input file
        if (temp == '|') {
                do{
                        temp = std::cin.get();
                        // start collecting info in the first non-blank
                        // non-'\n' character
                        if (temp!=' ' && temp!='\n')
                                start=true;
                        if ((temp == ' ' || temp == '\n') && start){
                                std::stringstream numStream(num);
                                numStream >> numInt;
                                proj.push_back(numInt);
                                num.clear();
                                while(temp == ' '){
                                        temp = std::cin.get();
                          }
                        }
                        num.push_back(temp);
                }while(temp != '\n');
                // if there is nothing after '|'
                if (!start){
                        proj = proj_all(m);
                }
        }else{ // use a default projection
                proj = proj_first_coord(D+1,m,mi);
        }
        PD = proj.size(); //this is the dimension of the projection
        sort(proj.begin(),proj.end());

        // compute cayley vector to augment pointset
        if (mi.size() != d+1){
                std::cout << "Input error" << std::endl;
                exit(-1);
        }

        // construct pointset
        std::string point;
        while(!std::getline(std::cin, point, ']').eof()) {
                std::vector<Field> ipoint;
                point.erase( std::remove( point.begin(), point.end(), ' ' ),
                             point.end() );
                point.erase( std::remove( point.begin(), point.end(), '[' ),
                             point.end() );
                std::string coord;
                if (point[0]==',')
                        point.erase(0,1);
                std::stringstream stream(point);
                if (!point.empty()){
                        while( getline(stream, coord, ',') ){
                                std::istringstream buffer(coord);
                                Field temp;
                                buffer >> temp;
                                ipoint.push_back(temp);
                        }
                        pointset.push_back(ipoint);
                }
        }
        if (m != pointset.size()){
                std::cout << "Input error" << std::endl;
                exit(-1);
        }
        return 0;
}
