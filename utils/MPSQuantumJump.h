#ifndef MPSQUANTUMJUMP_H
#define MPSQUANTUMJUMP_H

# include "itensor/all.h"
# include "itensor/util/parallel.h"
# include <vector>
# include <mpi.h>

using namespace itensor;
using namespace std;

using UGate = BondGate<ITensor>;
using IQUGate = BondGate<IQTensor>;

template <class T> 
struct TGate;
using OpGate = TGate<ITensor>;
using IQOpGate = TGate<IQTensor>;

template <class T>
class QuantumJump
    {
    public: 
    QuantumJump(double dt, 
                double Tf, 
                AutoMPO h, 
                vector<MPOt<T>> cm, 
                MPSt<T> psi_init, 
                vector<MPOt<T>> OBS, 
                int number_trial, 
                Args arg,
                Environment const& env);
    
    QuantumJump(double dt, 
                double Tf, 
                MPOt<T> U, 
                vector<MPOt<T>> cm, 
                MPSt<T> psi_init, 
                vector<MPOt<T>> OBS, 
                int number_trial, 
                Args arg,
                Environment const& env);
    
    QuantumJump(double dt, 
                double Tf, 
                vector<BondGate<T>> gates, 
                vector<TGate<T>> cm, 
                MPSt<T> psi_init, 
                vector<MPOt<T>> OBS, 
                int number_trial, 
                Args arg,
                Environment const& env);

    vector<vector<double>> getOBSvsT();
        
    void SteadyStateOnly(bool yn){SS_ = yn;};
    
    void Run();

    private:
    //MPS parameters
    Args args_;

    //Quantum Trajectory parameters
    double Tf_; 
    double dt_;
    int number_trial_;
    int Ntstep_;
    bool SS_;
    vector<vector<double>> o_vs_t_;

    // Multicore parameters
    Environment const* env_;
    int MPI_rank_;
    int MPI_nnodes_;

    //system parameters
    MPSt<T> psi_init_;
    int N_site_;
    vector<MPOt<T>> OBS_;
    int num_OBS_;

    // JumpOpType = 1 : MPO
    // JumpOpType = 2 : Bond gate
    int JumpOpType_;

    // If the time evolution is implemented by MPO
    vector<MPOt<T>> cm_;
    vector<MPOt<T>> CC_;
    MPOt<T> U_;
    

    // If the time evolution is implemented by bond gate
    vector<TGate<T>> cm_Gate_;
    vector<TGate<T>> CC_Gate_;
    vector<BondGate<T>> U_Gates_;

    void onetrial();
    };

template <class T> 
struct TGate
    {
    int i1 ;
    int i2 ;
    T G;
    bool TwoSite = false;
    TGate() { }
    TGate(int i1_, int i2_, T G_) 
    : i1(i1_), i2(i2_), G(G_) {TwoSite = true;}
    TGate(int i1_, T G_) 
    : i1(i1_), G(G_) { }
    };


#endif