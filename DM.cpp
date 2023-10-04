using namespace std;
#include <cmath>
#include <string>
#include <iostream>
#include<bits/stdc++.h>
#define EPSILON 0.00000000000001
#define  hbar 1.05457182*TMath::Power(10,-27)
#define  G 6.6732*TMath::Power(10,-8)
#define  c 2.9979*TMath::Power(10,10)
#define  Mo 1.987*TMath::Power(10,33)
#define  eo  TMath::Power(10,15)
#define  kappa c*c/(G*eo)
#define  yDM TMath::Power(10,3)
#define  mx 0.1*1.78266192*c*c*TMath::Power(10,-24)
#define  rho0 TMath::Power(10,14.3445)*kappa*G/(c*c)
#define  rho1 TMath::Power(10,14.7)*kappa*G/(c*c)
#define  rho2 TMath::Power(10,15)*kappa*G/(c*c)
#define  p1 TMath::Power(10,34.331)*kappa*G/(c*c*c*c)
#define  Gamma1 3.418
#define  Gamma2 2.835
#define  Gamma3 2.832
#define  K1 p1/TMath::Power(rho1,Gamma1)
#define  K2 p1/TMath::Power(rho1,Gamma2)
#define  K3 K2*TMath::Power(rho2,Gamma2-Gamma3)
#define  p2 K2*TMath::Power(rho2,Gamma2)
#define  guess 0.007

#define Mb  1.66*TMath::Power(10,-24)


double enDM(double x){
   return (mx*mx*mx*mx*((2*x*x*x+x)*sqrt(1+x*x)-TMath::ASinH(x))/(8*TMath::Pi()*TMath::Pi()*TMath::Power(hbar*c,3))+mx*mx*mx*mx*yDM*yDM*x*x*x*x*x*x/(TMath::Power(3*TMath::Pi()*TMath::Pi(),2)*TMath::Power(hbar*c,3)))/(c*c);
}

double func(double p,double x)
{
    return p*(c*c*c*c)/(kappa*G)-(mx*mx*mx*mx*((2*x*x*x-3*x)*sqrt(1+x*x)+3*TMath::ASinH(x))/(24*TMath::Pi()*TMath::Pi()*TMath::Power(hbar*c,3))+mx*mx*mx*mx*yDM*yDM*x*x*x*x*x*x/(TMath::Power(3*TMath::Pi()*TMath::Pi(),2)*TMath::Power(hbar*c,3)));
}


// Derivative of the above function which is 3*x^x - 2*x
double derivFunc(double x)
{
    return mx*mx*mx*mx*((8*x*x*x*x)/sqrt(x*x+1))/(24*TMath::Pi()*TMath::Pi()*TMath::Power(hbar*c,3))+mx*mx*mx*mx*yDM*yDM*x*x*x*x*x*5/(TMath::Power(3*TMath::Pi()*TMath::Pi(),2)*TMath::Power(hbar*c,3));
}


// Prints root of func(x) with error of EPSILON
double newtonRaphson(double p,double a)
{ 
   a=0.;
   double b=100.;
    if (func(p,a) * func(p,b) >= 0)
    {
        cout << "You have not assumed right a and b\n";
        return 0;
    }
 
    double x = a;
    while ((b-a) >= EPSILON)
    {
        // Find middle point
        x = (a+b)/2;
 
        // Check if middle point is root
        if (func(p,x) == 0.0)
            break;
 
        // Decide the side to repeat the steps
        else if (func(p,x)*func(p,a) < 0)
            b = x;
        else
            a = x;
    }
    return x;
}

void DM(){
  int Num=300;
  double p_center=pow(10,38)*kappa*G/(c*c*c*c);
  double p,nx,hx,ex;
  double h=p_center/pow(2,Num);
  string filename="Dummy";
 
  ofstream output;
  filename="eos"+filename;
  output.open(filename,ios::out);
  output<<Num<<endl;
  p=p_center/pow(2,Num);
  cout<<h<<" "<<p_center<<endl;
  for(int i=0;i<Num;i++){
        p=p*2;
        ex=enDM(newtonRaphson(p,guess));
        nx=pow(newtonRaphson(p,guess)*mx/(hbar*c),3)/(3*TMath::Pi()*TMath::Pi());
        hx=-TMath::Log((mx/(c*c))*nx/(ex+p/(kappa*G/(c*c*c*c))/(c*c)))*c*c;
      output<< std::setprecision(8) <<ex<<" "<<p/(kappa*G/(c*c*c*c))<<" "<<hx<<" "<<nx<<" "<<endl;

  }

  
}
