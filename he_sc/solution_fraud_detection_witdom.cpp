#include <fstream>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <math.h>
#include <nfl.hpp>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

int ti=0;
/// include the FV homomorphic encryption library
namespace FV {
namespace params {
  using poly_t = nfl::poly_from_modulus<uint64_t, 4096, 186>;
  template <typename T>
  struct plaintextModulus;
  template <>
  struct plaintextModulus<mpz_class> {
    static mpz_class value() { return mpz_class(std::to_string(ti)); } 
    };
  using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
  using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;

  gauss_t fg_prng_sk(20.0, 80, 1 << 12);
  gauss_t fg_prng_evk(20.0, 80, 1 << 12);
  gauss_t fg_prng_pk(20.0, 80, 1 << 12);
  gauss_t fg_prng_enc(20.0, 80, 1 << 12);
  }
} 

#include <FV.hpp>
#include "CRT.hpp"

typedef vector<double> record_t;
typedef vector<record_t> data_t;
typedef struct {record_t datay; data_t datax;} input_t;
typedef vector<FV::params::poly_p> record_poly_t;
typedef vector<record_poly_t> data_poly_t;
typedef vector<FV::ciphertext_t> record_cipher_t;
typedef vector<record_cipher_t> data_cipher_t;


//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read the input.
// Remember, a record is a list of doubles separated by commas ','.
istream& operator >> ( istream& ins, record_t& record )
  {
  record.clear();

  string line;
  getline( ins, line );

  stringstream ss( line );
  string field;
  while (getline( ss, field, ',' ))
    {
    stringstream fs( field );
    double f = 0; 
    fs >> f;
    record.push_back( f );
    }

  return ins;
  }

istream& operator >> ( istream& ins, input_t& input)
  {
  input.datay.clear();
  input.datax.clear();

  record_t record;
  while (ins >> record)
    {
      input.datay.push_back(record.back());
      record.pop_back();
      input.datax.push_back( record );
    }

  return ins;  
  }

//the sign function
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//convert the input to its w-NIBNAF expansion
FV::params::poly_p ConvertTowNIBNAFPol(long double input, double base, mpz_t mod, double precision)
{
  FV::params::poly_p output;
  const size_t degree=4096;
  std::array<mpz_t,degree> coeff;
  mpz_t rest,sum,sigmampz;
  mpz_inits(rest,sum,sigmampz,nullptr);
  for (size_t i=0; i<degree;i++)
    {
      mpz_inits(coeff[i],nullptr);
    }
  int sigma = sgn(input);
  long double t = sigma*input;
  while(t>precision)
    {
      int r = ceil(log(t)/log(base));
      if ((pow(base,r)-t)>(t-pow(base,r-1)))
	{
	  r=r-1;
	}
      if(r>-1)
      	{ 
          mpz_set_si(sigmampz,sigma);
          mpz_add(sum,sigmampz,mod);
          mpz_mod(rest,sum,mod);
          mpz_set(coeff[r],rest);
	} 
      else 
	{
          mpz_set_si(sigmampz,-sigma);
          mpz_add(sum,sigmampz,mod);
          mpz_mod(rest,sum,mod);
          mpz_set(coeff[degree+r],rest);
	}
      sigma=sigma*sgn(t-pow(base,r));
      t=abs(t-pow(base,r));
    }
  output.mpz2poly(coeff);
  mpz_clears(rest,sum,sigmampz,nullptr);
  for(size_t i=0; i < degree; i++)
    {
       mpz_clears(coeff[i],nullptr);
    }
  return output;
}

//convert from the w-NIBNAF expansion to a double
template<class I>
long double ConvertToNumber(long double base, mpz_t mod, int split, I &input)
{
  long double result=0;
  const size_t dim = sizeof(input)/sizeof(input[0]);
  const size_t neg_dim=split;
  mpz_t test,div,two,min;
  mpz_inits(test,div,two,min,nullptr);
  mpz_set_ui(two,2);
  mpz_cdiv_q(div,mod,two);
  for (unsigned long int i=0; i<dim; i++)
    {
      mpz_mod(test,input[i],mod);
	if (mpz_cmp(test,div)<0)
	{
	  mpz_set(input[i],test);
	}
      else
	{
	  mpz_sub(min,test,mod);
	  mpz_set(input[i],min);
	}
      long double coeffi=mpz_get_d(input[i]);
      if (coeffi !=0)
      {
	if(i<neg_dim)
	  {
	    result=result+coeffi*pow(base,i);
	  }
	else
	  {
	    int neg_power=i-dim;
	    result=result-coeffi*pow(base,neg_power);
	  }
      }
    }
  mpz_clears(test,div,two,min,nullptr);
  return result;
}


namespace FV{
ciphertext_t poly_to_ciphertext(pk_t &pk, params::poly_p const &poly)
  {
    ciphertext_t ct;
    ct.pk = &pk;
    ct.c0 = poly;
    ct.c0.ntt_pow_phi();
    ct.c0 = nfl::shoup(ct.c0 * ct.pk->delta, ct.pk->delta_shoup);
    ct.isnull = false;
    return ct;
  } 
}

//compute the approximation of the hessian
record_cipher_t compute_H(data_cipher_t d1, FV::ciphertext_t lambda, FV::ciphertext_t cte, FV::ciphertext_t denum, FV::ciphertext_t el_H)
{
  record_cipher_t result(d1[0].size(),el_H);
  record_cipher_t col_sum_row(d1.size(),el_H);
  for (size_t i=0;i<d1.size();i++)
    {
      for (size_t k=0;k<d1[0].size();k++)
	{
	  col_sum_row[i]+=d1[i][k];
	}
    }
  for (size_t j=0; j<d1[0].size();j++)
    {
      FV::ciphertext_t temp=el_H;
      for(size_t i=0; i<d1.size();i++)
	{
	  temp+=d1[i][j]*col_sum_row[i];
	}
      result[j]+=(cte*temp-lambda)*denum;
    }
  return result;
}

//compute the approximation of the gradient
record_cipher_t compute_gradient(FV::ciphertext_t zero, FV::ciphertext_t ahalf,FV::ciphertext_t aquarter, record_cipher_t y, data_cipher_t x, record_cipher_t w_new, FV::ciphertext_t lambda)
{ 
  record_cipher_t g(w_new.size(),zero);
  for (size_t i=0; i<x.size(); i++)
    {
      FV::ciphertext_t cte=zero;
      for (size_t j=0;j<w_new.size();j++)
	{
	  cte+=x[i][j]*w_new[j];
	}
      cte=(ahalf-aquarter*y[i]*cte)*y[i];
      for (size_t k=0;k<w_new.size();k++)
	{
	  g[k]+=cte*x[i][k];
	}
    }
   for (size_t k=0;k<w_new.size();k++)
    {
      g[k]-=lambda*w_new[k];
    }
  return g;
}

//functions for the timings
void diff(timespec &start, timespec &stop, timespec &difference)
{
  if ((stop.tv_nsec-start.tv_nsec)<0)
    {
      difference.tv_sec = stop.tv_sec - start.tv_sec -1;
      difference.tv_nsec = 1000000000+stop.tv_nsec-start.tv_nsec;
    }
  else
    {
      difference.tv_sec=stop.tv_sec-start.tv_sec;
      difference.tv_nsec = stop.tv_nsec-start.tv_nsec;
    }
}

void add(timespec &start, timespec &stop)
{
  start.tv_sec=start.tv_sec+stop.tv_sec;
  start.tv_nsec=start.tv_nsec+stop.tv_nsec;
}


//update the vector of parameters
record_cipher_t update_w(record_cipher_t w, record_cipher_t inv_H, record_cipher_t g)
{
  for(size_t i=0; i<w.size();i++)
    {
      w[i]=w[i]-inv_H[i]*g[i];
    }
  return w;
}

//the overall algorithm to perform one iteration of our fixed hessian method
void compute_parameter_vector_w(const char *inputfile, size_t nb_covariates,size_t n_training, long double * w)
{
  size_t d= nb_covariates+1;
  data_t datax;
  record_t datay;
  input_t input={datay,datax};
  ifstream infile( inputfile );
  infile >> input;

  // Complain if something went wrong with reading the input file.
  if (!infile.eof())
    {
    cout << "Problem reading in file\n";
    }

  infile.close();

  size_t Ndata = input.datax.size();

  for (size_t i=0; i<Ndata;i++)
    {
      input.datay[i]=2*input.datay[i]-1;
      input.datax[i].insert(input.datax[i].begin(),1);
    }
  data_t xtrain_set;
  xtrain_set.assign(input.datax.begin(), input.datax.begin()+n_training);
  record_t ytrain(input.datay.begin(), input.datay.begin()+n_training);

  data_t xtrain;
  record_t temp;
  for (size_t i=0;i<n_training;i++)
  { 
    temp.assign(xtrain_set[i].begin(), xtrain_set[i].begin()+d);
    xtrain.push_back(temp);
  }
  input.datax.clear();
  input.datay.clear();
  xtrain_set.clear();
  const size_t degree=4096;
  int split=300;
  long double base=1.02270245002238236963789094952;
  int primes[2]={2237,2239};
  size_t nb_primes=2;
  mpz_t mpz_primes[nb_primes];
  for (size_t j=0;j<nb_primes;j++)
    {
      mpz_init_set_ui(mpz_primes[j],primes[j]);
    } 
  mpz_t *** coeff_message = (mpz_t ***) malloc(sizeof(mpz_t **)*d);
  for (size_t i = 0; i < d; i++)
    {
      coeff_message[i] = (mpz_t **) malloc(sizeof(mpz_t *)*degree);
      for (size_t j = 0; j < degree; j++)
	{
	  coeff_message[i][j] = (mpz_t *) malloc(sizeof(mpz_t)*nb_primes);
	}
    }
  for (size_t k=0;k<d;k++)
    {
      for (size_t i=0; i<degree; i++)
	{
	  for (size_t j=0; j<nb_primes;j++)
	    {
	      mpz_inits(coeff_message[k][i][j],nullptr);
	    }  
	}
    }
  timespec start_encrypt_data,stop_encrypt_data,time_encrypt_data;
  time_encrypt_data.tv_sec=0.0;
  time_encrypt_data.tv_nsec=0.0;
  timespec start_computation,stop_computation,time_computation;
  time_computation.tv_sec=0.0;
  time_computation.tv_nsec=0.0;
  timespec start_decryption,stop_decryption,time_decryption;
  time_decryption.tv_sec=0.0;
  time_decryption.tv_nsec=0.0;
  timespec start_decoding, stop_decoding;
  for (size_t nb_prime=0;nb_prime<nb_primes;nb_prime++)
    {
    ti=primes[nb_prime];

    // KeyGen
    FV::sk_t secret_key;

    FV::evk_t evaluation_key(secret_key, 32);
    
    FV::pk_t public_key(secret_key, evaluation_key);

    mpz_t modi;
    mpz_inits(modi,nullptr);
    mpz_set_ui(modi,ti);

    double lambda=0.01;
    FV::params::poly_p enc_lambda=ConvertTowNIBNAFPol(lambda,base,modi,0.00001);
    FV::ciphertext_t encr_lambda=FV::poly_to_ciphertext(public_key,enc_lambda);
    
    //encrypting the input
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_encrypt_data);
    data_cipher_t encrypted_xtrain;
    record_cipher_t encrypted_ytrain;
    for(size_t i=0;i<n_training;i++)
      {
	record_cipher_t row;
	for (size_t j=0; j<d; j++)
	  { 
	    FV::params::poly_p tempx=ConvertTowNIBNAFPol(xtrain[i][j],base,modi,0.00001);
	    FV::ciphertext_t encrypt_tempx;
	    FV::encrypt_poly(encrypt_tempx,public_key,tempx);
	    row.push_back(encrypt_tempx);
	  }
	encrypted_xtrain.push_back(row);
	FV::params::poly_p tempy=ConvertTowNIBNAFPol(ytrain[i],base,modi,0.00001);
	FV::ciphertext_t encrypt_tempy;
	FV::encrypt_poly(encrypt_tempy,public_key,tempy);
	encrypted_ytrain.push_back(encrypt_tempy);
      }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop_encrypt_data);
    timespec difference;
    diff(start_encrypt_data,stop_encrypt_data,difference);
    add(time_encrypt_data,difference);

    //start the trainig process
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_computation);
    FV::params::poly_p cte=ConvertTowNIBNAFPol(-0.25,base,modi,0.00001);
    FV::ciphertext_t encr_cte=FV::poly_to_ciphertext(public_key,cte);
    FV::params::poly_p denum=ConvertTowNIBNAFPol(1.0/d,base,modi,0.00001);
    FV::ciphertext_t encr_denum=FV::poly_to_ciphertext(public_key,denum);
    FV::params::poly_p zero=ConvertTowNIBNAFPol(0,base,modi,0.00001);
    FV::ciphertext_t encr_zero=FV::poly_to_ciphertext(public_key,zero);
    record_cipher_t sum_H=compute_H(encrypted_xtrain,encr_lambda,encr_cte,encr_denum,encr_zero);

    FV::params::poly_p two=ConvertTowNIBNAFPol(2,base,modi,0.00001);
    FV::ciphertext_t encr_two=FV::poly_to_ciphertext(public_key,two);
    // double stepsized=-1;
    FV::params::poly_p stepsize;
    FV::ciphertext_t encr_stepsize;
    stepsize=ConvertTowNIBNAFPol(-0.1,base,modi,0.00001);
    encr_stepsize=FV::poly_to_ciphertext(public_key,stepsize);
      
    for(size_t i=0;i<sum_H.size();i++)
      {
	sum_H[i]=encr_stepsize*encr_two-(encr_stepsize*encr_stepsize)*sum_H[i];
      }
    FV::params::poly_p ahalf=ConvertTowNIBNAFPol(0.5,base,modi,0.00001);
    FV::ciphertext_t encr_ahalf=FV::poly_to_ciphertext(public_key,ahalf);
    FV::params::poly_p aquarter=ConvertTowNIBNAFPol(0.25,base,modi,0.00001);
    FV::ciphertext_t encr_aquarter=FV::poly_to_ciphertext(public_key,aquarter);

    FV::params::poly_p startvalue=ConvertTowNIBNAFPol(0.001,base,modi,0.00001);
    FV::ciphertext_t encr_startvalue=FV::poly_to_ciphertext(public_key,startvalue);
    record_cipher_t w(d,encr_startvalue);
    record_cipher_t g=compute_gradient(encr_zero, encr_ahalf,encr_aquarter,encrypted_ytrain, encrypted_xtrain, w, encr_lambda);

    w=update_w(w,sum_H,g);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop_computation);
    timespec difference_computation;
    diff(start_computation,stop_computation,difference_computation);
    add(time_computation,difference_computation);

    //start preparing decryption of the result
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_decryption);
    std::array<mpz_t, degree> poly_message;
    for (size_t i = 0; i <degree; i++) 
      {
	mpz_inits(poly_message[i], nullptr);
      }

    for (size_t k=0;k<d;k++)
      {
	FV::decrypt_poly(poly_message,secret_key, public_key,w[k]);
  
	for (size_t i = 0; i <degree; i++) 
	  {
	    mpz_set(coeff_message[k][i][nb_prime],poly_message[i]);
	  }
      }
    
    //actual decryption of the result
    for (size_t i = 0; i < degree; i++) 
      {
	mpz_clears(poly_message[i],nullptr);
      }
    mpz_clears(modi,nullptr);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop_decryption);
    timespec difference_decryption;
    diff(start_decryption,stop_decryption,difference_decryption);
    add(time_decryption,difference_decryption);
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start_decoding);
  mpz_t coeffN[degree];
  for (size_t i=0;i<degree;i++)
  {
    mpz_inits(coeffN[i],nullptr);
  }
  mpz_t prod;
  mpz_inits(prod,nullptr);
  mpz_set_ui(prod,1);
  for (size_t k=0; k<nb_primes; k++)
  { 
    mpz_mul(prod,prod,mpz_primes[k]);
  }
  for (size_t k=0; k<d; k++)
    {
      for(size_t i=0; i<degree; i++)
	{
	  findMinX(coeffN[i],mpz_primes,coeff_message[k][i],nb_primes);
	}
      w[k]=ConvertToNumber(base,prod,split,coeffN);
    }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop_decoding);
  timespec difference_decoding;
  diff(start_decoding, stop_decoding, difference_decoding);
  add(time_decryption,difference_decoding);

  std::cout<<"time_encrypt_data "<<time_encrypt_data.tv_sec<<" sec "<<time_encrypt_data.tv_nsec<<"\n";
  std::cout<<"time_comp_server "<<time_computation.tv_sec<<" sec "<<time_computation.tv_nsec<<"\n";
  std::cout<<"time_decryption "<<time_decryption.tv_sec<<" sec "<<time_decryption.tv_nsec<<"\n";
  for (size_t k=0;k<d;k++)
    {
      for (size_t i=0; i<degree; i++)
	{
	  for (size_t j=0; j<nb_primes;j++)
	    {
	      mpz_clears(coeff_message[k][i][j],nullptr);
	    }  
	}
    }
  for (size_t i=0;i<d;i++)
    {
      if(coeff_message[i]!=NULL)
	{
	for(size_t j=0;j<degree;j++)
	  {
	    free(coeff_message[i][j]);
	  }
	free(coeff_message[i]);
      }
    }
  free(coeff_message);
  mpz_clears(prod,nullptr);  
  for (size_t i=0; i<nb_primes; i++)
  {
    mpz_clears(mpz_primes[i],nullptr);
  }
  for(size_t k=0; k<degree;k++)
  {
     mpz_clears(coeffN[k],nullptr);
  }
}


int main()
{
  const char *input="prepared_financial_data.txt"; //input file
  size_t d=31; //number of covariates
  size_t n_training=20; //number of training records
  long double* w =(long double*) malloc((d+1)*sizeof(long double));
  compute_parameter_vector_w(input,d, n_training,w);
  std::cout<<"w=[";
  for (size_t k=0; k<d+1;k++)
    {
      std::cout.precision(20);
      std::cout<<w[k]<<(k==d ? "":"; ");
    }
  std::cout<<"]; \n";
  free(w);
  return 0;
}
