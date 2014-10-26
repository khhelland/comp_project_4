#include <armadillo>

using namespace arma;

vec trimul(double a, double b, double c, vec v)
{
  //Function for multiplying a tridiagonal matrix with the 
  //diagonals filled with doubles a,b and c.
  //Returns product
  int n = v.n_elem;
  vec u(n);
  
  u(0) = b*v(0) + c*v(1);
  for(int i = 1; i < n-1; i++)
    {
      u(i) = a*v(i-1) + b*v(i) + c*v(i+1);
    }
  u(n-1) = a*v(n-2) + b*v(n-1);
  return u;
}


void trisolve(double a, double b, double c, vec &u, vec v)
{
  //Function for solving equation Au = v for u.
  //A is a tridiagonal matrix with 
  //diagonals filled with doubles a,b and c,
  // v is known.
  // u is modified in place
 
  int n = v.n_elem;
  vec bv(n);
  bv.fill(b);
  double ac = a*c;
  
  //First: make the matrix upper diagonal:
  for(int i=1; i<n; i++)
    {
      bv(i) -= ac/bv(i-1);
      v(i)  -= a*v(i-1)/bv(i-1);
    }
  
  //Backward substitution to obtain u
  u(n-1) = v(n-1)/bv(n-1);
  
  for(int i= (n-2); i>=0 ; i--)
    {
      u(i) = (v(i) - c*u(i+1))/bv(i);
    }
}

