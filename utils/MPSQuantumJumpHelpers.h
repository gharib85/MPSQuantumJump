#ifndef __MPSQUANTUMJUMPHELPERS__
#define __MPSQUANTUMJUMPHELPERS__

# include "itensor/all.h"
using namespace itensor;
using namespace std;

inline int whichsitejump( vector<double> dpj, double& dp)
    {
    double e2 = ((double) rand() / (RAND_MAX));
    double sumpj=0.;
    int jjump=0;
    for (size_t j=0 ; j < dpj.size(); j++){
        sumpj += dpj.at(j)/dp;
        if (sumpj>e2){
            //cout << "operator #" << j << "jump" << endl;
            jjump=j;
            break;
        }
    }
    return jjump;
    }//whichsitejump

template<class T>
MPOt<T> 
HermConj_MPO(MPOt<T> M)
    {
    auto HM = M;
    HM.mapprime(1,2,Site);
    HM.mapprime(0,1,Site);
    HM.mapprime(2,0,Site);
    for (int i=1; i<= HM.N(); i++){
        HM.Aref(i).dag();
    }

    return HM;
    }//HermConj

    
#endif

