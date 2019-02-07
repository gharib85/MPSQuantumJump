# MPSQuantumJump
Quantum Trajectory Simulation with Matrix Product State based on the [iTensor](http://itensor.org/index.html) library. 
Message Passing interface (MPI) is used to run the simulation in parallel. 

# Compile and Link `MPSQuantumJump`
## Step 0: make sure you have already installed iTensor and MPI.
To install iTensor, see http://itensor.org/index.html

To install MPI, see https://www.open-mpi.org/

## Step 1: Set `CCCOM` parameter to compile the code with MPI. 
1. Go to the itensor directory`cd /path/to/itensor/` and copy the `options.mk` as `options_mpi.mk`.
2. Open `options_mpi.mk` and modify the `CCCOM` parameter:
* For GNU GCC compiler
`CCCOM=mpic++ -m64 -std=c++11 -fPIC`
* For Clang compiler 
`CCCOM=mpic++ -std=c++11 -fPIC`
* For GCC On Windows cygwin
`CCCOM=mpic++ -std=c++11 -Wa,-mbig-obj -O2 -fPIC`

## Step 2: Compile the code by using the `Makefile` provided in the example folder.
In the `Makefile`, set the `iTensor` path in the parameter `LIBRARY_DIR` and set the `MPSQuantumJump` path in the parameter `MPSQuantumJump_DIR`. 

## Step 3: Run the code
Run it with
`mpirun -np # APP Input`

replace `#` by the number of core, `APP` by the filename of the compiled code and `Input` by the input parameters file.

# Example
I used the dissipative transverse Ising model as an example. 

Consider a 1D spin-1/2 chain, the interaction between the nearest neighbor spins is Ising interaction in x-direction and the magnetic field is applied in z-direction. Meanwhile, there exists local dissipation described by the jump operator S^-. 

The system consists of three parameters:

`Jx` : strength of the Ising interaction

`Jz` : strength of the magnetic field

`gamma` : decay rate of the local dissipation

One can modify these parameters in the `inputs` file.
