#include <armadillo>
#include <fstream>
#include "tridiag.h"

using namespace arma;
using namespace std;


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
  ForwardEuler(alpha/2, u);
  BackwardEuler(alpha/2,u);
}


void solve(double dt, double dx, int N, vec v,
           void (*method)(double, vec&),
           const char* outfile )
{
  ofstream out(outfile);
  out<<v.t();
  double alpha = dt/(dx*dx);
  
  
  for(int n=0; n<N; n++)
    {
      method(alpha,v);
      
      //if (((n+1)%20) == 0)
      out << v.t();
    }
}

int main ()
{
  double dt = 0.01;
  double dx = 0.1;
  int N  = 300;
  int n = 9;
  vec v(n);
  //fill v
  for(int i = 0; i < n ; i++)
    v(i) = -1 + (i+1)*dx;
  
  solve(dt,dx,N,v,*ForwardEuler,"ForwardEuler.dat");
  
  // vec v = randu<vec>(3);
  // double a = v(0);
  // double b = v(1);
  // double c = v(2);
  // mat A = zeros<mat>(3,3);
  // A.diag(-1).fill(a);  
  // A.diag().fill(b);
  // A.diag(1).fill(c);
  // v = randu<vec>(3);
  // vec u = vec(3);
  // trisolve(a,b,c,u,v);
  // vec w = vec(3);
  // solve(w,A,v);

  // cout << u<<endl<<w;
  

  return 0;
}
