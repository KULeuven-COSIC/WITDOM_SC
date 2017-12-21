#include<bits/stdc++.h>
using namespace std;
 
// Returns modulo inverse of a with respect to m using extended
// Euclid Algorithm. Refer below post for details:
// http://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
void inv(mpz_t x1, mpz_t a_input, mpz_t m_input)
{
  mpz_t m0, t, q,temp,x0,m,a;
  mpz_inits(m0,t,q,temp,x0,m,a,nullptr);
  mpz_set(m0,m_input);
  mpz_set(m,m_input);
  mpz_set(a,a_input);
  mpz_set_si(x0,0);
  mpz_set_si(x1,1);
 
  if (mpz_cmp_ui(m,1)==0)
    mpz_set_si(x1,0);
 
  // Apply extended Euclid Algorithm
  while (mpz_cmp_ui(a,1)>0)
  {
    // q is quotient
    mpz_fdiv_q(q,a,m);
    mpz_set(t,m);

    // m is remainder now, process same as
    // euclid's algo
    mpz_fdiv_r(m,a,m);
    mpz_set(a,t);
 
    mpz_set(t,x0);
    mpz_mul(temp,q,x0);
    mpz_sub(x0,x1,temp);
    mpz_set(x1,t);
  }
 
  // Make x1 positive
  if (mpz_cmp_ui(x1,0)<0)
    mpz_add(x1,x1,m0);

  mpz_clears(m0,t,q,temp,x0,m,a,nullptr);
}


// k is size of num[] and rem[].  Returns the smallest
// number x such that:
//  x % num[0] = rem[0],
//  x % num[1] = rem[1],
//  ..................
//  x % num[k-2] = rem[k-1]
// Assumption: Numbers in num[] are pairwise coprime
// (gcd for every pair is 1)
void findMinX(mpz_t result, mpz_t num[], mpz_t rem[], int k)
{
  // Compute product of all numbers
  mpz_t prod,pp,inverse,temp;
  mpz_inits(prod,pp,inverse,temp,nullptr);
  mpz_set_ui(prod,1);
  for (int i = 0; i < k; i++)
  {
    mpz_mul(prod,prod,num[i]);
  }

  // Initialize result
  mpz_set_ui(result,0);
 
  // Apply above formula
  for (int i = 0; i < k; i++)
    {// std::cout<<result<<"+";
     // std::cout<<rem[i]<<"\n";
    mpz_divexact(pp,prod,num[i]);
    inv(inverse,pp,num[i]);
    mpz_mul(temp,inverse,pp);
    mpz_mul(temp,temp,rem[i]);
    // std::cout<<temp<<"=";
    mpz_add(result,result,temp);
    //std::cout<<mpz_class(result).get_str()<<"\n";
  }
 
  mpz_mod(result,result,prod);
  mpz_clears(prod,pp,inverse,temp,nullptr);
}
