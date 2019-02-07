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
## Step 2: Compile the code by using the `Makefile` 
Use the `Makefile` provided here to compile the code. 

In the `Makefile`, set the `iTensor` path in the parameter `LIBRARY_DIR` and set the `MPSQuantumJump` path in the parameter `MPSQuantumJump_DIR`. 

See MPSQuantumJump/example/ for more details.


