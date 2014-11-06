#include <armadillo>
#include <fstream>
#include <time.h>
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
  BackwardEuler(alpha/2, u);
}


void solve(double dt, double dx, double T, vec v,
           void (*method)(double, vec&),
           const char* outfile )
{
  ofstream out(outfile);
  //out<<v.t();
  double alpha = dt/(dx*dx);
  
  bool print1 = false;
  bool print2 = false;
  
  for(double t=0; t<T; t+=dt)
    {
      method(alpha,v);
      
      if (t>0.02 && !print1) 
        {
          out << v.t(); 
          print1 = true;
        }
      else if (t>0.5 && !print2)
        {
          out<<v.t();
          print2 = true;
        }
      //out << v.t();
    }
}

int main ()
{
  
  double dx = 0.1;
  double dt = 0.49*dx*dx;
  double T  = 2;
  int n = 9;
  vec v(n);
  //fill v
  for(int i = 0; i < n ; i++)
    v(i) = -1 + (i+1)*dx;
  //cout<<v;
 
  clock_t start, mid1, mid2, end;
  start = clock();
  solve(dt, dx, T, v, *ForwardEuler, "ForwardEuler.dat");
  mid1 = clock();
  solve(dt, dx, T, v, *BackwardEuler, "BackwardEuler.dat");
  mid2 = clock();
  solve(dt, dx, T, v, *CrankNicholson, "CrankNicholson.dat");
  end = clock();
  
  double cps = CLOCKS_PER_SEC;
  double FEtime = (mid1-start)/cps;
  double BEtime = (mid2-mid1)/cps;
  double CNtime = (end-mid2)/cps;
  
  cout << FEtime <<"\t" << BEtime << "\t" <<CNtime<<endl;

  return 0;
}
