//TropLi_disc Source Code --- version 0.1


#include <gmpxx.h>
#include <Eigen/Dense>
#include <LEDA/numbers/integer_matrix.h>
#include <LEDA/numbers/integer_vector.h>
#include <LEDA/numbers/integer.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <list>
#include <time.h>
#include <sstream>


using std::cout;
using std::set;
using std::vector;
using std::map;
using std::list;
using std::pair;
using std::ofstream;
using std::string;
using std::stringstream;
using std::ios;
using std::cin;
using leda::integer_vector;
using leda::integer_matrix;
using leda::integer;



void computepreferences (unsigned long k,
                         list<unsigned long> *order,
                         integer_matrix *pDMred,
                         list<unsigned long> pvec [],
                         list<unsigned long> F [],
                         unsigned long sub [],
                         integer_matrix *pT,
                         integer_matrix *pvectors,
                         integer_vector *ptempvec,
                         bool bits [],
                         bool bits2 []);

void computecone (list<unsigned long> *porder,
                  list<unsigned long> pvec [],
                  list<unsigned long> F [],
                  unsigned long sub [],
                  integer_matrix *pT,
                  integer_matrix *pvectors,
                  integer_vector *ptempvec,
                  bool bits [],
                  bool bits2 []);

unsigned long numray ( bool bits [] );


//initial declarations
unsigned long n, m, d, i, j, q, r, s, ray, ncones, nconesbergman, nrays;
bool admissible;

ofstream outfilecones;
ofstream outfilerays;
ofstream outfilevertices;

list<unsigned long> templist;
set<unsigned long> B;
set<unsigned long> picked;
map<string, unsigned long> rays;

map<string, bool> fulldimpartition;

vector< integer_vector > raysvector;
vector< list<unsigned long> >  raysofcone;

string str;
stringstream strpartition;

set<unsigned long>::iterator it2;
list<unsigned long>::iterator lit2;
list<unsigned long>::iterator lit3;
list<unsigned long>::reverse_iterator rlit;
pair< set<unsigned long>::iterator , bool > ret;
pair< map<string, unsigned long>::iterator , bool > ret2;


bool flagrandom = false;





//main
int main(int argc, char* argv[]) {

    //unsigned long d2, n2, m2;
    cin >> d;
    cin >> n;

    typedef Eigen::Matrix<mpq_class, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    MatrixXd M2(d, n);

    for (auto i=0; i<d; i++) {
        for (auto j=0; j<n; j++) {
            cin >> M2(i,j);
        }
    }

    //cout << M2 << std::endl;

  //initial declarations
  unsigned long t, t2;
  time_t tstart, tend;
  unsigned long row;
  string outstring;

  long sizew, sizeu;

  //read matrix M
  integer_matrix M;
  cin >> M;

  //define n and m and d
  //d = M.dim1();
  //n = M.dim2();
  m = n-d;

  //determine flags
  int nws = 1;
  outstring = "disc";
  sizew = 1000; //random w's will have coordinates between -sizew and sizew
  sizeu = 500;  //random perturbations will have coordinates between -sizeu and sizeu

  for (i=1; i<argc; i++) {
    str = argv[i];
    if (str == "-random") {
      flagrandom = true;
      nws = atoi(argv[i+1]);  //nws = number of vectors w
    }
    if (str.substr(0,5) == "-out:") {
      outstring = str.erase(0,5);
    }
    if (str == "-sizew") sizew = atol(argv[i+1]);  //nws = number of vectors w
    if (str == "-sizeperturbation") sizeu = atol(argv[i+1]);  //nws = number of vectors w
  }

  //TESTTESTTESTTEST
  //cout << "sizew = " << sizew << ". sizeu = " << sizeu << "\n";

  //read w
  integer_matrix w(n,nws);
  if (!flagrandom) {
    for (i=0; i<n; i++) {
      cin >> w(i,0);
    }
  }
  MatrixXd w2(n,nws);
  if (!flagrandom) {
    for (i=0; i<n; i++) {
      cin >> w2(i,0);
    }
  }

  //get time
  time (&tstart);

  //verify rank is m
  //if (rank(M)!=d) {
  Eigen::FullPivLU<MatrixXd> lu_decomp(M2);
  if (lu_decomp.rank()!=d) {
    cout << "The rank of the matrix is not equal to the number of rows.\n";
    return 0;
  }

  //define matrices for ray-shooting
  integer_matrix T(n-1,n);
  integer_vector tempvec(n);
  integer_matrix vectors(n,n);

  MatrixXd T2(n-1,n);
  MatrixXd tempvec(n,1);
  MatrixXd vectors(n,n);

  //definitions
  integer_matrix N(d,d);
  integer_matrix INV(d,d);
  integer deter;
  integer_vector vec(d);

  integer_matrix Mred(d,n);
  integer_matrix DMred(m,n);

  MatrixXd N2(d,d);
  MatrixXd INV2(d,d);
  mpq_t deter2;
  MatrixXd vec2(d,1);

  MatrixXd Mred(d,n);
  MatrixXd DMred(m,n);

  unsigned long sub [m];
  list<unsigned long> F [n];
  list<unsigned long> pvec [m];
  list<unsigned long> order;

  bool bits [n];
  bool bits2 [n];

  //open output files
  outfilecones.open ( (outstring + "_maximalcones.out").c_str() );
  outfilerays.open ( (outstring + "_rays.out").c_str() );

  //more definitions
  bool possible, first;

  //fill some part of the matrix T with the matrix M
  for (j=0; j<n; j++) {
    for (i=0; i<d; i++) {
        //T(i+m-1,j) = M(i,j);
        T2(i+m-1,j) = M2(i,j);
    }
  }

  //initializing variables
  t2 = -1;
  t = m;
  ncones=0;
  nconesbergman=0;
  nrays=0;

  //computes the first m-subset of [n]
  for ( j = 0; j < t; j++ ) {
    sub[m-t+j] = t2 + j + 1;
  }
  if ( t2+1 < n-t ) {
    t = 0;
  }
  t++;
  t2 = sub[m-t];

  //translates sub[] into a set B
  B.clear();
  ret = B.insert(sub[0]);
  it2 = ret.first;
  for (i=1; i<m; i++) {
    B.insert(it2, sub[i]);
    it2++;
  }

  for (j=m; j<n; j++) { //fills the initial matrix N with the corresponding columns of M
    for (i=0; i<d; i++) {
      N2(i,j-m) = M2(i,j);
    }
  }


  while (sub[0] == 0) {

    if (inverse(N, INV, deter, vec)) { //if B is a basis

      Mred = INV * M; //Reduces the columns of M corresponding to the complement of basis B

      //will compute F[i]
      possible = true;
      row=0;
      for (i=0; i<n; i++) {

    (F[i]).clear();

    if (B.find(i) == B.end()) { //if i is not in B

          first = true;

          for (j=0; j<m; j++) {
            if (Mred(row,sub[j]) != 0) {
          DMred(j,i) = 1; //Fills the Dual Reduced Matrix with a nonzero entry
          (F[i]).push_back(j);
              if (first && i < sub[j]) {
                possible = false;
                j=m;
                i=n;
              }
              first = false;
            }
            else {
          DMred(j,i)=0; //Fills the Dual Reduced Matrix with a zero entry
        }
      }
      row++;

    }
    else { //if i is in B
      for (j=0; j<m; j++) {
        if (j == i-row) {
          DMred(j,i)=1;
        }
        else {
          DMred(j,i)=0;
        }
      }
      (F[i]).push_back(i-row);
    }

      }

      if (possible) { // will compute the "preferences" (compatible pairs) arising from basis B
        order.clear();
        for (i=0; i<m; i++) {
          (pvec[i]).clear();
        }
        picked.clear();
    computepreferences(0, &order, &DMred, pvec, F, sub, &T, &vectors, &tempvec, bits, bits2);

    fulldimpartition.clear();
      }

    }

    //computes the next m-subset of [n]
    for ( j = 0; j < t; j++ ) {
      sub[m-t+j] = t2 + j + 1;
    }
    if ( t2+1 < n-t ) {
      t = 0;
    }
    t++;
    t2 = sub[m-t];

    //translates sub[] into a set B
    B.clear();
    ret = B.insert(sub[0]);
    it2 = ret.first;
    for (i=1; i<m; i++) {
      B.insert(it2, sub[i]);
      it2++;
    }

    row=0; //fills the matrix N with the corresponding columns of M
    for (j=0; j<n; j++) {
      if (B.find(j) == B.end()) {
    for (i=0; i<d; i++) {
          N(i,row) = M(i,j);
    }
    row++;
      }
    }

  }

  //close files
  outfilecones.close();
  outfilerays.close();

  //print output
  cout << "\nNumber of maximal cones in the dual cyclic Bergman fan = " << nconesbergman << "\n";
  cout << "Number of maximal cones having codim. 1 after adding rowspace(A) = " << ncones << "\n";
  cout << "Number of rays in the subfan they generate = " << nrays << "\n";
  cout << "Coordinates of these rays written in the file: " << outstring + "_rays.out" << "\n";
  cout << "Maximal cones in this subfan written in the file: " << outstring + "_maximalcones.out" << "\n";


  //get time
  time (&tend);
  cout << "Time elapsed so far: " << difftime (tend,tstart) << " seconds\n\n";
  cout << "Performing ray shooting algorithm...\n";

  //open file
  outfilevertices.open (  (outstring + "_vertices.out").c_str() );


  //definitions
  srand( time(NULL) ); //used for random choices

  integer_matrix C(n,n);
  integer_matrix Cinv(n,n);
  integer_matrix expo(n,nws);
  integer_matrix x(n,nws);
  integer_vector vect(n);
  integer_vector u(n);
  integer_vector y(n);
  integer Dint;
  bool intersects [nws];
  bool repeat;

  integer tempinteger;


  //fill last columns of C, which are always constant
  for (i=m; i<n; i++) {
    for (j=0; j<n; j++) {
      C(j,i) = M(i-m,j);
    }
  }

  for (i=0; i<1; i++) { //will try these many blocks of w's //NOT USED SO FAR

    if (flagrandom) { //get the vectors w at random
      for (q=0; q<nws; q++) {
        for (j=0; j<n; j++) {
          w(j,q) = rand() % (2*sizew+1) - sizew;
        }
      }
    }

    //TESTTESTTESTTEST
    //cout << "w : " << transpose(w) << "\n";

    repeat = true;

    while (repeat==true) {

      //get the perturbation u
      for (j=0; j<n; j++) {
    u[j] = rand() % (2*sizeu+1) - sizeu;
      }

      //initialize exponent vectors
      for (q=0; q<nws; q++) {
    for (j=0; j<n; j++) {
      expo(j,q)=0;
    }
      }

      repeat = false;

      for(j=0; j<ncones; j++) { //will test the intersection of    w + epsilon * u + R>0 * e_r   with cone number j

    //fill C with the cone j
    r=1;
    for ( lit2=(raysofcone[j]).begin(); lit2!=(raysofcone[j]).end(); lit2++ ) {
      for (s=0; s<n; s++) {
        C(s,r) = raysvector[*lit2][s];
      }
      r++;
    }

    for (r=0; r<n; r++) { //will update coordinate r of exponent vector determined by w

      //fill first column of C
      C(r,0) = -1;

      if (inverse(C, Cinv, Dint, vect)) { //since w + epsilon * u is generic, they only intersect when C is invertible

        //tests if the ray shot intersects
        for (q=0; q<nws; q++) {
          intersects[q] = true;
        }
        x = Cinv * w;
        y = Cinv * u;

        //change signs if Dint is negative
        if (Dint < 0) {
          for (s=0; s<n; s++) {
        y[s] =  -y[s];
        for (q=0; q<nws; q++) {
          x(s,q) = - x(s,q);
        }
          }
        }

        for (q=0; q<nws; q++) {
          for (s=0; s<m; s++) {
        if ( x(s,q) < 0 || (x(s,q) == 0 && y[s] < 0) ) {
          intersects[q] = false;
          s=m;
        }
        else {
          if (x(s,q) == 0 && y[s] == 0) { //u was not generic!!!
            repeat = true;
            s=m;
            q=nws;
            C(r,0) = 0; //erase first column of C
            r=n;
            j=ncones;
            cout<< "Non-generic perturbation. Restarting...\n";
          }
        }
          }
        }

        if (repeat == false) { //if u was generic
          // if it intersects then add the corresponding coordinate to the exponent vector
          for (q=0; q<nws; q++) {
            if (intersects[q] == true) {
          expo(r,q) = expo(r,q) + abs(determinant(C));
            }
          }
        }

      }

      if (repeat == false) { //if u was generic
        //erase first column of C
        C(r,0) = 0;
      }

    }

      }

    }

    //print output
    cout << "Vertex matrix written in the file: " << outstring + "_vertices.out" << "\n";
    outfilevertices << transpose(expo) << "\n";


    //compute A-degree
    integer_vector Adeg(d);
    Adeg = M * expo.col(0);
    cout << "A-degree of the A-discriminant:\n";
    for (j=0; j<d; j++) {
      cout << Adeg[j] << " ";
    }
    cout << "\n";

  }

  //close file
  outfilevertices.close();


  //print time
  time (&tend);
  cout << "Total time elapsed: " << difftime (tend,tstart) << " seconds\n\n";

  return 0;

}













//computes the "preferences" (compatible pairs) for columns starting at column k
void computepreferences (unsigned long k, list<unsigned long> *porder, integer_matrix *pDMred, list<unsigned long> pvec [], list<unsigned long> F [], unsigned long sub [], integer_matrix *pT, integer_matrix *pvectors, integer_vector *ptempvec, bool bits [], bool bits2 []) {


  list<unsigned long>::iterator lit;
  list<unsigned long>::iterator littemp;
  list<unsigned long>::iterator it;
  set<unsigned long> S;


  if (k < n) {

    if (B.find(k) == B.end()) { //if k is not in the basis B

      //will compute the first element in the list 'order' whose corresponding entry in column k is non-zero
      lit = (*porder).begin();
      while ( lit != (*porder).end() && (*pDMred)(*lit,k) == 0 ) {
        lit++;
      }

      //will compute the set S of columns having a preference with >= weight than *lit
      S.clear();
      for (lit2 = lit; lit2 != (*porder).end(); lit2++) {
        for (lit3 = (pvec[*lit2]).begin(); lit3 != (pvec[*lit2]).end(); lit3++) {
          S.insert(*lit3);
        }
      }

      if (lit != (*porder).end()) { //will compute preferences if column k prefers *lit
        (pvec[*lit]).push_back(k);
    computepreferences (k+1, porder, pDMred, pvec, F, sub, pT, pvectors, ptempvec, bits, bits2); //computepreferences for next column
        (pvec[*lit]).pop_back();
      }

      //will run over all allowed weights already picked
      do {

        for (it = F[k].begin(); it != F[k].end(); it++) { //will compute preferences if column k picks a -new- number *it with weight right below *lit
          if (sub[*it] < k && picked.find(*it)==picked.end()) { //only accepts preferences to basis columns to the left and which are new
            admissible = true;
            for (it2 = S.begin(); it2 != S.end(); it2++) { //tests there is no conflict with other columns
              if ( (*pDMred)(*it, *it2) != 0 ) {
                admissible = false;
                it2 = --S.end();
              }
            }
            if (admissible) {
              (pvec[*it]).push_back(k);
              picked.insert(*it);
              littemp = (*porder).insert (lit, *it);
          computepreferences (k+1, porder, pDMred, pvec, F, sub, pT, pvectors, ptempvec, bits, bits2); //computepreferences for next column
              picked.erase(*it);
              (*porder).erase(littemp);
              (pvec[*it]).pop_back();
            }

          }

        }

        lit--;

        //will update S
        if (lit != --( (*porder).begin() ) ) {
          for (lit3 = (pvec[*lit]).begin(); lit3 != (pvec[*lit]).end(); lit3++) {
            S.insert(*lit3);
          }
        }

      } while ( lit != --( (*porder).begin() ) );

    }
    else { //if k is a basis column
      computepreferences (k+1, porder, pDMred, pvec, F, sub, pT, pvectors, ptempvec, bits, bits2);
    }

  }
  else { //if all preferences are made

    computecone (porder, pvec, F, sub, pT, pvectors, ptempvec, bits, bits2);
    ++nconesbergman;

  }

  return;
}










//computes the cone corresponding to some "preferences" and an order
void computecone (list<unsigned long> *porder, list<unsigned long> pvec [], list<unsigned long> F [], unsigned long sub [], integer_matrix *pT, integer_matrix *pvectors, integer_vector *ptempvec, bool bits [], bool bits2 []) {

  //computes the string encoding the partition
  strpartition.str("");
  for (i=0; i<m; i++) {
    strpartition << sub[i];
    for (lit3 = (pvec[i]).begin(); lit3 != (pvec[i]).end(); lit3++) {
      strpartition << "-" << *lit3;
    }
    strpartition << " ";
  }

  if ( fulldimpartition.find( strpartition.str() ) == fulldimpartition.end() ) { //if the partition has not been tested for codimension 1

    ray=0;
    for(j=0; j<n; j++) {
      (*ptempvec)[j] = 0;
      for (i=0; i<m-1; i++) {
    (*pT)(i,j)=0;
      }
    }

    for (rlit = (*porder).rbegin(); rlit != (*porder).rend(); rlit++) { //will run over all non-singleton sets in our partition

      //will compute the leaves right on top of our set
      for (lit2 = (pvec[*rlit]).begin(); lit2 != (pvec[*rlit]).end(); lit2++) {

        (*ptempvec)[*lit2] = 1; //adds the elements of our set in the partition
        for (lit3 = F[*lit2].begin(); lit3 != F[*lit2].end(); lit3++) {
      if ( (*ptempvec)[sub[*lit3]] == 0 && *lit3 != *rlit ) {
        //will write the corresponding ray
        (*pT)(ray,sub[*lit3])=1;
        ray++;
        (*ptempvec)[sub[*lit3]] = 1;
      }
        }

      }
      (*ptempvec)[sub[*rlit]] = 1; //finishes computing the new big ray

      if (rlit != --(*porder).rend() ) { //checks if we are not the all 1s ray
        (*pT).row(ray) = *ptempvec;
        ray++;
      }
      else { //if we are in the all 1s ray
        if ( homogeneous_linear_solver(*pT, *pvectors) == 1 ) {  //if it has codimension 1
      fulldimpartition.insert ( pair< string, bool >(strpartition.str(), true) );
    }
        else { //if it has greater codimension
      fulldimpartition.insert ( pair< string, bool >(strpartition.str(), false) );
    }
      }

    }
  }

  if ( fulldimpartition[strpartition.str()] == true ) {  //compute rays if it has codimension 1

    templist.clear();

    //will compute the rays of the cone
    for(j=0; j<n; j++) {
      bits[j] = false;
      bits2[j] = false;
    }
    for (rlit = (*porder).rbegin(); rlit != (*porder).rend(); rlit++) { //will run over all non-singleton sets in our partition

      //will compute the leaves right on top of our set
      for (lit2 = (pvec[*rlit]).begin(); lit2 != (pvec[*rlit]).end(); lit2++) {

        bits[*lit2] = true; //adds the elements of our set in the partition
        for (lit3 = F[*lit2].begin(); lit3 != F[*lit2].end(); lit3++) {
      if ( bits[sub[*lit3]] == false && *lit3 != *rlit ) {
        //add the corresponding ray
        bits2[sub[*lit3]] = true;
        templist.push_back( numray(bits2) );
        bits2[sub[*lit3]] = false;
        bits[sub[*lit3]] = true;
      }
        }

      }
      bits[sub[*rlit]] = true; //finishes computing the new big ray

      if (rlit != --(*porder).rend() ) { //if we are not in the the all 1s ray
        templist.push_back( numray(bits) );
      }
      else { //if we are in the all 1s ray
        raysofcone.push_back( templist );

    //outputs rays of this cone
    outfilecones << "{";
    for (lit2 = templist.begin(); lit2 != templist.end(); lit2++ ) {
      outfilecones << *lit2 << " ";
    }
    outfilecones.seekp (-1, ios::cur);
    outfilecones << "}\n";

    ncones++;
      }

    }

  }

}


//returns the number corresponding to the ray given by bits
unsigned long numray ( bool bits [] ) {

  str = "";

  //computes de string corresponding to S
  for (j=0; j<n; j++) {
    if (bits[j]) {
      str += "1";
    }
    else {
      str += "0";
    }
  }

  //returns the number assigned to the ray, or assigns it a new one
  ret2 = rays.insert ( pair<string, unsigned long>(str, nrays) );
  if (ret2.second) { //if it was a new ray
    integer_vector tempint(n);
    for (j=0; j<n; j++) { //computes the ray and its number
      if ( str.at(j) == '1' ) {
    tempint [j] = 1;
    outfilerays << "1 ";
      }
      else {
    tempint [j] = 0;
    outfilerays << "0 ";
      }
    }
    raysvector.push_back (tempint);
    outfilerays << " # " << nrays << "\n"; //outputs the ray
    nrays++;
  }
  return ret2.first->second;

}
