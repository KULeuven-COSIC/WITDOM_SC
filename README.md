# WITDOM_SC
Contains all code that can be made available related to WITDOM_SC and COSIC implementations
## Homomorphic encryption SC
### 1. Getting started
Open the file solution_fraud_detection_witdom.cpp, go to the main part of the file at the bottom and change the necessary parameters (input file, number of covariates, number of training vectors). Take care if you change the number of covariates and the number of training vectors, other parameters might need to be finetuned. 
After saving the file you compile with
g++ -std=c++11 -funroll-loops -Ofast -Wall -g -I nfllib/include test_big_parameters.cpp -o tsolution_fraud_detection_witdom -L nfllib/lib -lnfllib -lmpfr -lgmpxx -lgmp
with the paths adapted to the correct folders.
Then you can run the program with: 
./solution_fraud_detection

### 2. Installation of the prerequisites
For ubuntu use the following commands
apt-get install gcc<br />
apt-get install libgmp-dev<br />
apt-get install libmpfr-dev<br />
apt-get install git<br />
apt-get install cmake<br />
git clone https://github.com/quarkslab/NFLlib.git<br />
cd NFLlib/<br />
mkdir _build<br />
cd build/<br />
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/nfllib<br />
make<br />
make test<br />
make install<br />
export LD_LIBARY_PATH=/../nfllib/lib/(path naar folder met libnfllib.so file)<br />

### 3. Licensing
The license of the FV-NFLlib and the NFLlib is GPLv3. 
## MPC Component
### Prerequisites
Install SPDZ-2. Tested up to version 0.2. <link>https://github.com/bristolcrypto/SPDZ-2<link/>. <br />
Configure SPDZ-2 to run in 256 bits inputs. <br />
Configure for the Number of parties involved in the computation. (As many as needed) <br />
### 1. Getting Started
Open the file gmdh.mpc, go to the main part of the file and change the parameters i.e. number of inputs for training and testing respectively. 
Save the file and compile it following the instructions detailed on the SPDZ-2 repository. <br />
Add secret inputs normalized (values between 0 and 1) with a precision of 6 decimals. Inputs can be introduced by following the instructions on the SPDZ-2 repository and Google group  <link>https://groups.google.com/forum/#!forum/spdz<link/>. <br />
Run the file and compile it following the instructions detailed on the SPDZ-2 repository. <br />

### 3. Licensing
The license of the FV-NFLlib and the NFLlib is GPLv3. 
