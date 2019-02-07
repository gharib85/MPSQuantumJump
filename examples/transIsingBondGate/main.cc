# include "QuantumJump.h"
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
    //-------- Jump operator ------------------------------
	//Jump operators cm : sigma-_i
    vector<OpGate> cm;
    for(int j = 1; j <= N_site; ++j)
        {
        auto T = sqrt(gamma)*sites.op("S-",j);
        cm.emplace_back(j, T);
        }

    //-----------------------------------------------------
    //-------- Time evolution operator --------------------
	auto gates = vector<UGate>();
	for(int j = 1; j < N_site; ++j)
        {
        ITensor hh;
		if(j == 1)
            {
            //non hermitian part 
		    hh += -Cplx_i * gamma * sites.op("projUp",j)* sites.op("Id",j+1);
            hh += -Cplx_i/2 * gamma * sites.op("Id",j)* sites.op("projUp",j+1);
            
            //hermitian part
            //Sz
            hh += Jz * sites.op("Sz",j)* sites.op("Id",j+1);
            hh += Jz/2 * sites.op("Id",j)* sites.op("Sz",j+1);
            }
        else if(j == N_site-1)
            {
            //non hermitian part 
		    hh += -Cplx_i/2 * gamma * sites.op("projUp",j)* sites.op("Id",j+1);
            hh += -Cplx_i * gamma * sites.op("Id",j)* sites.op("projUp",j+1);
            
            //hermitian part
            //Sz
            hh += Jz/2 * sites.op("Sz",j)* sites.op("Id",j+1);
            hh += Jz * sites.op("Id",j)* sites.op("Sz",j+1);
            }
        else
            {
            //non hermitian part 
		    hh += -Cplx_i/2 * gamma * sites.op("projUp",j)* sites.op("Id",j+1);
            hh += -Cplx_i/2 * gamma * sites.op("Id",j)* sites.op("projUp",j+1);
            
            //hermitian part
            //Sz
            hh += Jz/2 * sites.op("Sz",j)* sites.op("Id",j+1);
            hh += Jz/2 * sites.op("Id",j)* sites.op("Sz",j+1);
            }
        //Sx*Sx
        hh += Jx * ITensor(sites.op("Sx",j))* ITensor(sites.op("Sx",j+1));
        
        auto g = UGate(sites,j,j+1,UGate::tReal,dt/2.,hh);
    	gates.push_back(g);
        }
        
    for(int j = N_site-1; j >=1; --j)
        {
        ITensor hh;
		if(j == 1)
            {
            //non hermitian part 
		    hh += -Cplx_i * gamma * sites.op("projUp",j)* sites.op("Id",j+1);
            hh += -Cplx_i/2 * gamma * sites.op("Id",j)* sites.op("projUp",j+1);
            
            //hermitian part
            //Sz
            hh += Jz * sites.op("Sz",j)* sites.op("Id",j+1);
            hh += Jz/2 * sites.op("Id",j)* sites.op("Sz",j+1);
            }
        else if(j == N_site-1)
            {
            //non hermitian part 
		    hh += -Cplx_i/2 * gamma * sites.op("projUp",j)* sites.op("Id",j+1);
            hh += -Cplx_i * gamma * sites.op("Id",j)* sites.op("projUp",j+1);
            
            //hermitian part
            //Sz
            hh += Jz/2 * sites.op("Sz",j)* sites.op("Id",j+1);
            hh += Jz * sites.op("Id",j)* sites.op("Sz",j+1);
            }
        else
            {
            //non hermitian part 
		    hh += -Cplx_i/2 * gamma * sites.op("projUp",j)* sites.op("Id",j+1);
            hh += -Cplx_i/2 * gamma * sites.op("Id",j)* sites.op("projUp",j+1);
            
            //hermitian part
            //Sz
            hh += Jz/2 * sites.op("Sz",j)* sites.op("Id",j+1);
            hh += Jz/2 * sites.op("Id",j)* sites.op("Sz",j+1);
            }
        //Sx*Sx
        hh += Jx * ITensor(sites.op("Sx",j))* ITensor(sites.op("Sx",j+1));
        
        auto g = UGate(sites,j,j+1,UGate::tReal,dt/2.,hh);
        
    	gates.push_back(g);
        }

	//-----------------------------------------------------
    //-------- Observable ---------------------------------
	vector<MPO> OBS;
    for(int j = 1; j <= N_site; j++)
        {
        auto ampoSz = AutoMPO(sites);
		ampoSz += "Sz",j;
		OBS.push_back(MPO(ampoSz));
        }

	for (int j = 2; j <= N_site; j++)
        {
		auto ampocorr = AutoMPO(sites);
		ampocorr += "Sx",1,"Sx",j;
		OBS.push_back(MPO(ampocorr));
	    }
	
    for (int j = 2; j <= N_site; j++)
        {
		auto ampocorr = AutoMPO(sites);
		ampocorr += "Sz",1,"Sz",j;
		OBS.push_back(MPO(ampocorr));
	    }
    
//-----------------------------------------------------------------
//-------------------Run the Quantum trajectory -------------------
//-------------------And output the data to file ------------------

	//implement QuanutmJump
	auto QJ = QuantumJump<ITensor>( dt, tmax, gates, cm,  
                                psi_init,OBS, trial_per_node, args, env);
	QJ.SteadyStateOnly(SS);
	clock_t begin = clock();
	QJ.Run();
	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << env.rank() << "\t" << elapsed_secs << endl;
	
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
