#include "bird.h"

int main(){
        typedef int                                     NT;
        typedef svector<NT>                             row;
        typedef svector<row>                            matrix;

        matrix A,B,C,D;
        row r1,r2,r3,r4,r5,r6;

        r1.push_back(1);r1.push_back(3);r1.push_back(2);
        r2.push_back(1);r2.push_back(2);r2.push_back(4);
        r3.push_back(1);r3.push_back(1);r3.push_back(1);
        A.push_back(r1);A.push_back(r2);A.push_back(r3);

        r1.push_back(5);
        r2.push_back(3);
        r3.push_back(8);
        r4.push_back(3);r4.push_back(4);r4.push_back(5);r4.push_back(6);
        B.push_back(r1);B.push_back(r2);B.push_back(r3);B.push_back(r4);

        r1.push_back(4);
        r2.push_back(2);
        r3.push_back(11);
        r4.push_back(4);
        r5.push_back(7);r5.push_back(4);r5.push_back(3);
        r5.push_back(5);r5.push_back(1);
        C.push_back(r1);C.push_back(r2);C.push_back(r3);
        C.push_back(r4);C.push_back(r5);

        r1.push_back(7);
        r2.push_back(3);
        r3.push_back(2);
        r4.push_back(2);
        r5.push_back(9);
        r6.push_back(5);r6.push_back(4);r6.push_back(3);
        r6.push_back(7);r6.push_back(5);r6.push_back(1);
        D.push_back(r1);D.push_back(r2);D.push_back(r3);
        D.push_back(r4);D.push_back(r5);D.push_back(r6);

        // These determinants should be 5, -72, 75 and -134.
        std::cout<<"det(A) = "<<Bird::determinant(A)<<std::endl;
        std::cout<<"det(B) = "<<Bird::determinant(B)<<std::endl;
        std::cout<<"det(C) = "<<Bird::determinant(C)<<std::endl;
        std::cout<<"det(D) = "<<Bird::determinant(D)<<std::endl;

        return 0;
}
