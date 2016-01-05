 
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void func(double* y0, double x);
void RKstep(double *y0, double* yn , double dx ,double x);
int main()
{
	double y0[3], yn[3];
	const int dim=3;
	y0[0]=1.,  y0[1]=1., y0[2]=1.;
	double x_end=100. , dx=0.001, x=0;
	
	ofstream out("runge4.txt");
	while(x<x_end)
	{
		x+=dx;
		RKstep(y0, yn, dx, x);
		for(int i=0; i<dim; i++) y0[i]=yn[i];
		out<<x<<"\t"<<y0[0]<<"\t"<<y0[1]<<"\t"<<y0[2]<<endl;
	}
	out.close();
  return 0;
}
void RKstep(double *y0, double* yn , double dx ,double x)
{
	const int dim=3;
	double k1[dim], k2[dim], k3[dim], k4[dim];

	for(int i=0; i<dim; i++) k1[i]=y0[i];
	func(k1, x);
	for(int i=0; i<dim; i++)
		k2[i]=y0[i]+0.5*dx*k1[i];
	func(k2,x+0.5*dx);
	for(int i=0; i<dim; i++)
		k3[i]=y0[i]+0.5*dx*k2[i];
	func(k3, x+0.5*dx);
	for(int i=0; i<dim; i++)
		k4[i]=y0[i]+dx*k3[i];
	func(k4,x+dx);
	for(int j=0; j<dim; j++)
		yn[j]=y0[j]+dx/6.*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
}
void func(double* y0, double x)
{
	const double a=10.;
	const double b=28.;
	const double c=8./3. ;
	double y[3]={y0[0], y0[1], y0[2]};
	y0[0]=a*(y[1]-y[0]);
	y0[1]=y[0]*(b-y[2])-y[1];
	y0[2]=y[0]*y[1]-c*y[2];
	 

}
