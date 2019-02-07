# include "MPSQuantumJump.h"
#include <iostream>
#include <iomanip>
using namespace std;

int main ( int argc, char *argv[] )
{
	
//---------------------------------------------------------
//----------------------initialize MPI---------------------
	Environment env(argc,argv);

//---------------------------------------------------------
//---------parse parameters from input file ---------------
	auto input = InputGroup(argv[1],"input");
	
	//MPS parameters
	auto cutoff = input.getReal("cutoff",1E-8);
	auto Maxm = input.getInt("Maxm",50);
	auto args = Args("Cutoff=",cutoff,"Maxm=",Maxm);
	
	//Quantum Trajectory parameters
	auto Ntrial = input.getInt("Ntrial",100);
	auto tmax = input.getReal("tmax",2);	
	auto dt = input.getReal("dt",1E-3);
	int trial_per_node = Ntrial/env.nnodes();
	
	//outputfile parameters
	auto output = input.getString("output","data.txt");
	
	//system parameters
	auto N_site = input.getInt("N_site");
	auto SS = input.getYesNo("Steady_state_Only", true);

	// transverse Ising model parameters
    auto gamma = input.getReal("gamma",0.0);
    auto Jz = input.getReal("Jz",0.0);
    auto Jx = input.getReal("Jx",-1);

//---------------------------------------------------------
//----setup Hamiltonian, initial psi, jump operator.-------
	
	//Setup Siteset
	auto sites = SpinHalf(N_site);
	
    //-----------------------------------------------------
    //---------Initial Condition---------------------------
    // i.e. state.set(i,"Up") set the state of site i as up
    //      state.set(i,"Dn") set the state of site i as dn
    // The following script set the inital psi = |udud...>
	auto state = InitState(sites);
	for(int i = 1; i <= N_site; ++i)
   	    {
   		if(i%2 == 1) state.set(i,"Up");
        else state.set(i,"Dn");
    	}	
	auto psi_init = MPS(state);
    //-----------------------------------------------------
    //------------------------------------------------------

	
	//-----------------------------------------------------
    //-------- Hamiltonian --------------------------------
    // AutoMPO automatically convert the string and
    //  position label to Matrix Product Operator. 
    //-------------------------------------------
	auto ampoH = AutoMPO(sites);
	// XX interaction
	for(int i = 1; i < N_site; i++)
		{
		ampoH += Jx, "Sx",i, "Sx",i+1;
		}
	// Z interaction
	for(int i = 1; i <= N_site; i++)
        {
        ampoH += Jz, "Sz",i;
        }
	//Imaginary part of Hamiltonian
    for(int i = 1; i <= N_site; i++)
        {
		ampoH += -Cplx_i*gamma,"projUp",i;
		}

    //-----------------------------------------------------
    //-------- Jump operator ------------------------------
	//Jump operators cm : sigma-_i
	vector<MPO> cm;
	for (int i=1; i<=N_site; i++)
        {
        auto ampoSm = AutoMPO(sites);
		ampoSm += sqrt(gamma),"S-",i;
        cm.push_back(MPO(ampoSm)); 
	    }
	
    //-----------------------------------------------------
    //-------- Observable ---------------------------------
	vector<MPO> OBS;
    // Sz
	for (int i=1; i<=N_site; i++)
        {
        auto ampoSz = AutoMPO(sites);
		ampoSz += "Sz",i;
        OBS.push_back(MPO(ampoSz));
	    }

    // XX correlation function
	for (int m=2; m<=N_site; m++)
        {
        auto ampoxx = AutoMPO(sites);
		ampoxx += "Sx",1, "Sx",m;
        OBS.push_back(MPO(ampoxx));
	    }
	
    // ZZ correlation function
	for (int m=2; m<=N_site; m++)
        {
        auto ampozz = AutoMPO(sites);
		ampozz += "Sz",1, "Sz",m;
        OBS.push_back(MPO(ampozz));
	    }

//-----------------------------------------------------------------
//-------------------Run the Quantum trajectory -------------------
//-------------------And output the data to file ------------------

	//implement QuanutmJump
	auto QJ = QuantumJump<ITensor>( dt, tmax, ampoH, cm,  
                                psi_init,OBS, trial_per_node, args, env);
	QJ.SteadyStateOnly(SS);
	QJ.Run();

    // Write result to file
    cout.precision(10);
	if(env.rank()==0)
        {
		auto M = QJ.getOBSvsT();
		//write to file
		ofstream myfile;
  		myfile.open (output);
  		double t=0;
  		myfile << setw(10);

  		if (!SS) myfile << "t" <<"\t"<< setw(10);

		for (size_t k = 0; k < M.size(); k++)
            {
			myfile << "OBS" <<k<<"\t"<< setw(10);
		    }
		myfile << endl;

		for (size_t m = 0; m < M.at(0).size(); m++)
            {
			myfile << setw(10);
			if (!SS){
			myfile << t <<"\t"<< setw(10);
			t+=dt;
			}
			for (size_t k = 0; k < M.size(); k++)
                {
				myfile << M.at(k).at(m) <<"\t"<< setw(10);
			    }
			myfile << endl;
		    }
	    }


	
	return 0;
}
