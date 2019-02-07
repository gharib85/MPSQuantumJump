#ifndef __MPSQUANTUMJUMPHELPERS__
#define __MPSQUANTUMJUMPHELPERS__

# include "itensor/all.h"
# include "MPSQuantumJump.h"
using namespace itensor;


template<class T> T  
HermConj_Tensor(T M)
    {
    auto HM = dag(M);
    HM.mapprime(1,2,Site);
    HM.mapprime(0,1,Site);
    HM.mapprime(2,0,Site);
    return HM;
    }//HermConj

template<class T> 
void //C = BA
nmultITensor(T A, T B, T& C)
    {
    auto B_ = B;
    B_.prime();
    C = B_*A;
    C.mapprime(2,1);
    }

template<class T>
void 
ApplyTGate(MPSt<T>& psi, TGate<T> gate, Args args)
    {
    if (gate.TwoSite)
        {
        auto b = gate.i1;
        auto& G = gate.G;
        psi.position(b);

        auto AA = psi.A(b)*psi.A(b+1);
        AA = AA*G;
        AA.mapprime(1,0);

        //Normalize AA after applying G
        AA /= norm(AA);
        //SVD AA to restore MPS form
        auto U = psi.A(b);
        T D,V;
        svd(AA,U,D,V,args);
        psi.setA(b,U);
        psi.setA(b+1,D*V);
        }
    else
        {
        auto b = gate.i1;
        auto& G = gate.G;
        psi.position(b);
        auto AA = psi.A(b);
        AA = AA*G;
        AA.mapprime(1,0);
        AA /= norm(AA);
        psi.setA(b,AA);
        }
    }

template <class Iterable, class T>
void 
applyGate(Iterable const& gatelist,
MPSt<T>& psi, 
Args args)
{
    auto g = gatelist.begin();
    while(g != gatelist.end())
    {
        auto i1 = g->i1();
        auto i2 = g->i2();
        auto AA = psi.A(i1)*psi.A(i2)*g->gate();
        AA.mapprime(1,0,Site);

        ++g;
        if(g != gatelist.end())
            {
            //Look ahead to next gate position
            auto ni1 = g->i1();
            auto ni2 = g->i2();
            //SVD AA to restore MPS form
            //before applying current gate
            if(ni1 >= i2)
                {
                psi.svdBond(i1,AA,Fromleft,args);
                psi.position(ni1); //does no work if position already ni1
                }
            else
                {
                psi.svdBond(i1,AA,Fromright,args);
                psi.position(ni2); //does no work if position already ni2
                }
            }
        else
            {
            //No next gate to analyze, just restore MPS form
            psi.svdBond(i1,AA,Fromright,args);
            }
    }
}

template<class T>
double 
overlap_TGate(MPSt<T> psi, TGate<T> gate)
    {
    if (gate.TwoSite)
        {
        auto b = gate.i1;
        psi.position(b);
        auto hh = dag(prime(psi.A(b),Site))*dag(prime(psi.A(b+1),Site))*gate.G*psi.A(b)*psi.A(b+1);
        return (hh.cplx()).real();
        }
    else 
        {
        auto b = gate.i1;
        psi.position(b);
        return (dag(prime(psi.A(b),Site))*gate.G*psi.A(b)).real();
        }
    }




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

