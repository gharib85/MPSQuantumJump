# include "itensor/all.h"
# include <vector>
# include "MPSQuantumJump.h"
# include <iostream>
#include "MPSQuantumJumpHelpers.h"

using namespace itensor;
using namespace std;

template<typename T>
QuantumJump<T>::QuantumJump( 
    double dt, 
    double Tf, 
    AutoMPO h, 
    vector<MPOt<T>> cm, 
    MPSt<T> psi_init, 
    vector<MPOt<T>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env)
    : env_(&env)
    {
    JumpOpType_ = 1;
    
    //MPS parameters
    args_ = arg;
    args_.add("DoSVDBond", true);

    //Quantum Trajectory parameters
    Tf_ = Tf;
    dt_ = dt;
    number_trial_ = number_trial;
    Ntstep_ = Tf_/dt_;
    SS_ = false;

    // MPI parameters
    MPI_rank_ = env_ -> rank();
    MPI_nnodes_ = env_ -> nnodes();

    //system parameters
    psi_init_ = psi_init; 
    N_site_ = psi_init_.N();
    OBS_ = OBS;
    cm_ = cm;
    num_OBS_ = OBS_.size();
    U_ = toExpH<T>(h,dt_*Cplx_i);

    // C^dagger * C^-
    CC_.resize(cm_.size());
    for (size_t i=0; i<cm_.size(); i++)
        {
        nmultMPO(cm_.at(i), HermConj_MPO<T>(cm_.at(i)), CC_.at(i));
        }
    }

template<typename T>
QuantumJump<T>::QuantumJump(
    double dt, 
    double Tf, 
    MPOt<T> U, 
    vector<MPOt<T>> cm, 
    MPSt<T> psi_init, 
    vector<MPOt<T>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env)
    : env_(&env)
    {
    JumpOpType_ = 1;
    
    //MPS parameters
    args_ = arg;
    args_.add("DoSVDBond", true);

    //Quantum Trajectory parameters
    Tf_ = Tf;
    dt_ = dt;
    number_trial_ = number_trial;
    Ntstep_ = Tf_/dt_;
    SS_ = false;

    // MPI parameters
    MPI_rank_ = env_ -> rank();
    MPI_nnodes_ = env_ -> nnodes();

    //system parameters
    psi_init_ = psi_init; 
    N_site_ = psi_init_.N();
    OBS_ = OBS;
    cm_ = cm;
    num_OBS_ = OBS_.size();
    U_ = U;

    // C^dagger * C^-
    CC_.resize(cm_.size());
    for (size_t i=0; i<cm_.size(); i++)
        {
        nmultMPO(cm_.at(i), HermConj_MPO<T>(cm_.at(i)), CC_.at(i));
        }
    }

template<typename T>
QuantumJump<T>::QuantumJump( 
    double dt, 
    double Tf, 
    vector<TGate<T>> gates, 
    vector<TGate<T>> cm, 
    MPSt<T> psi_init, 
    vector<MPOt<T>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env)
    : env_(&env)
    {
    JumpOpType_ = 2;
    //MPS parameters
    args_ = arg;

    //Quantum Trajectory parameters
    Tf_ = Tf;
    dt_ = dt;
    number_trial_ = number_trial;
    Ntstep_ = Tf_/dt_;
    SS_ = false;

    // MPI parameters
    MPI_rank_ = env_ -> rank();
    MPI_nnodes_ = env_ -> nnodes();

    //system parameters
    psi_init_ = psi_init; 
    N_site_ = psi_init_.N();
    OBS_ = OBS;
    cm_Gate_ = cm;
    num_OBS_ = OBS_.size();

    // Time Evolution Operator
    U_Gates_ = gates;

    // C^dagger * C^-
    CC_Gate_ = cm_Gate_;
    for (size_t i=0; i<cm_Gate_.size(); i++)
        {
        nmultITensor(cm_Gate_.at(i).G, HermConj_Tensor<T>(cm_Gate_.at(i).G), CC_Gate_.at(i).G );
        }
    }


template<typename T>
void 
QuantumJump<T>::Run()
    {
    cout << "Start Quatum Jump"<< endl;
    // initialize o_vs_t
    o_vs_t_.resize(num_OBS_);
    for (int i=0; i< num_OBS_; i++)
        {
        if (!SS_)
            {
            o_vs_t_.at(i).resize(Ntstep_);
            }
        else
            {
            o_vs_t_.at(i).resize(1);
            } 
        }

    // compute quantum trajectory
    srand(MPI_rank_*1000 + time(NULL));
    int c = 0;
    while (c<number_trial_)
        {
        onetrial();
        c++;
        cout <<"Rank = "<< MPI_rank_ << ": Processing... " 
                << (double(c))/number_trial_*100 <<"%\n";
        }

    // collect data
    env_ -> barrier();
    cout <<"Rank = "<< MPI_rank_ << " Start Send Data" << endl;
    if(MPI_rank_ == 0)
        {
        for (int i = 1; i < MPI_nnodes_; ++i)
            {
            MailBox mailbox(*env_,i);
            vector<vector<Cplx>> tmp;
            mailbox.receive(tmp);
            for (int n=0; n<num_OBS_; n++)
                {
                transform (o_vs_t_.at(n).begin(), o_vs_t_.at(n).end(), tmp.at(n).begin(), 
                                                    o_vs_t_.at(n).begin(), plus<Cplx>());}
            }
        }
    else
        {
        MailBox mailbox(*env_,0);
        mailbox.send(o_vs_t_);
        }
    env_ -> barrier();
    }//Run


template<typename T>
void 
QuantumJump<T>::onetrial()
    {
    auto psi = psi_init_;
    for (int step=0; step<Ntstep_; step++)
        {
        //if (step %10 == 0)
        //    {
        //    cout << "Rank = "<< MPI_rank_ << ": step = " << step << "\t" << maxM(psi) << endl;
        //    }
        // Do measurement
        if (!SS_)
            {
            for (int c=0; c< num_OBS_; c++)
                {
                auto O = overlapC(psi, OBS_.at(c), psi);
                o_vs_t_.at(c).at(step)+=O/number_trial_/MPI_nnodes_;
                }
            }
        else
            {
            if (step == Ntstep_-1)
                {
                for (int c=0; c< num_OBS_; c++)
                    {
                    auto O = overlapC(psi, OBS_.at(c), psi);
                    o_vs_t_.at(c).at(0)+=O/number_trial_/MPI_nnodes_;
                    }
                }
            }
        // Time evolution
        if (JumpOpType_==1)
            {
            vector<double> dpj(CC_.size());

            // Calculate dp
            double dp=0.;
            for (size_t j=0; j < CC_.size(); j++)
                {
                dpj.at(j) = dt_ * real(overlapC(psi, CC_.at(j), psi));
                dp += dpj.at(j);
                }
            //if(dp > 0.1)
            //    {
            //      cout << "Warning : dp = " << dp << ". Please decrease the time step dt." << endl;
            //    }

            // Generate randum number
            double e = ((double) rand() / (RAND_MAX));

            // unitary time evolution
            exactApplyMPO(psi,U_,psi,args_);

            if (e<dp) 
                {  
                int siteJ=whichsitejump(dpj,dp);
                exactApplyMPO(psi,cm_[siteJ],psi,args_);
                psi.position(1,args_);
                }
            
            psi.position(1);
            normalize(psi);
            }
        else if(JumpOpType_==2)
            {
            vector<double> dpj(CC_Gate_.size());
            // Calculate dp
            double dp=0.;
            for (size_t j=0; j < CC_Gate_.size(); j++)
                {
                dpj.at(j) = dt_ * overlap_TGate(psi, CC_Gate_.at(j));
                dp += dpj.at(j);
                }
            //if(dp > 0.1)
            //    {
            //    cout << "Warning : dp = " << dp << ". Please decrease the time step dt." << endl;
            //    }
            
            // Generate randum number
            double e = ((double) rand() / (RAND_MAX));
            
            // unitary time evolution
            applyGate(U_Gates_, psi, args_);

            if (e<dp)
                {
                int siteJ=whichsitejump(dpj,dp);
                ApplyTGate<T>(psi, cm_Gate_.at(siteJ),args_);
                }
            }
        psi.position(1);
        normalize(psi);
        }
    }

template<typename T>
vector<vector<Cplx>>
QuantumJump<T>::getOBSvsT()
        {
        if (MPI_rank_==0)
            {
            return o_vs_t_;
            } 
        else
            {
            error("Only root can access OBSvs_t");
            return vector<vector<Cplx>>();
            }
        }


template
QuantumJump<ITensor>::QuantumJump( 
    double dt, 
    double Tf, 
    AutoMPO h, 
    vector<MPOt<ITensor>> cm, 
    MPSt<ITensor> psi_init, 
    vector<MPOt<ITensor>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env);

template
QuantumJump<ITensor>::QuantumJump( 
    double dt, 
    double Tf, 
    MPOt<ITensor> U, 
    vector<MPOt<ITensor>> cm, 
    MPSt<ITensor> psi_init, 
    vector<MPOt<ITensor>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env);

template
QuantumJump<IQTensor>::QuantumJump( 
    double dt, 
    double Tf, 
    AutoMPO h, 
    vector<MPOt<IQTensor>> cm, 
    MPSt<IQTensor> psi_init, 
    vector<MPOt<IQTensor>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env);

template
QuantumJump<IQTensor>::QuantumJump( 
    double dt, 
    double Tf, 
    MPOt<IQTensor> U, 
    vector<MPOt<IQTensor>> cm, 
    MPSt<IQTensor> psi_init, 
    vector<MPOt<IQTensor>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env);

template
QuantumJump<ITensor>::QuantumJump( 
    double dt, 
    double Tf, 
    vector<TGate<ITensor>> gates, 
    vector<TGate<ITensor>> cm, 
    MPSt<ITensor> psi_init, 
    vector<MPOt<ITensor>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env);

template
QuantumJump<IQTensor>::QuantumJump( 
    double dt, 
    double Tf, 
    vector<TGate<IQTensor>> gates, 
    vector<TGate<IQTensor>> cm, 
    MPSt<IQTensor> psi_init, 
    vector<MPOt<IQTensor>> OBS, 
    int number_trial, 
    Args arg,
    Environment const& env);
template void QuantumJump<ITensor>::Run();
template void QuantumJump<IQTensor>::Run();
template void QuantumJump<ITensor>::onetrial();
template void QuantumJump<IQTensor>::onetrial();
template vector<vector<Cplx>> QuantumJump<ITensor>::getOBSvsT();
template vector<vector<Cplx>> QuantumJump<IQTensor>::getOBSvsT();
