#include <Armadillo>
#include <fstream>

using namespace arma;
using namespace std;


vec trimul(double a, double b, double c, vec v)
{
  //Function for multiplying a tridiagonal matrix with the 
  //diagonals filled with doubles a,b and c.
  //Returns product
  int n = v.n_elem;
  vec u(n);
  
  for(int i = 0; i < n; i++)
    {
      u(i) = a*v(i-1) + b*v(i) + c*v(i-1);
    }
  
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
      bv(i+1) -= ac/bv(i);
      v(i+1)  -= a*v(i)/bv(i);
    }
  
  //Backward substitution to obtain u
  u(n) = v(n)/bv(n);
  
  for(int i= (n-1); i>0 ; i--)
    {
      u(i) = (v(i) - a*u(i+1))/bv(i);
    }
}

//Schemes for solving diffusion equation

void ForwardEuler(double alpha, vec &u)
{
  u = trimul(alpha,1-2*alpha,alpha,u);
}

void BackwardEuler(double alpha, vec &u)
{
  vec v = u;
  trisolve(-alpha, 1+2*alpha, -alpha, u,v);
}

void CrankNicholson(double alpha, vec &u)
{
  ForwardEuler(alpha, u);
  BackwardEuler(alpha,u);
}


void solve(double dt, double dx, double T, vec v,
           void (*method)(double, vec),
           const char* outfile )
{
  ofstream out(outfile);
  out<<v.t();

  double t = 0;
  double alpha = dt/(dx*dx);
  
  
  while (t < T)
    {
      method(alpha,v);
      
      //if ((int)(t/4) == 0){out << v.t();}
      
      t += dt;
    }
}

int main ()
{
  double dt = 1;
  double dx = 1;
  double T = 10;
  vec v();
  //fill v
  
  solve (dt,dx,T,v,*ForwardEuler,"ForwardEuler.dat")
  

  return 0;
}
