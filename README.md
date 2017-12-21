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
