#include "stdafx.h"
#include <iostream>
#include "conio.h"
#include <fstream>

void runge_kutta ( double *x, double *y,double *z, const int n, const double dx, const double a, const double b, const double c);
double funx(double x1,  double y1,const double a);
double funy (double x1, double y1, double z1, const double b);
double funz (double x1, double y1,double z1, const double c);



using namespace std;

int main()
{
	double x1;
	double y1;
	double z1;
	const double dx=0.1;
	const int n=1000;
	
	//double f[2];
	double *x=new double[n];
	double *y=new double[n];
	double *z=new double[n];
	//double *f=new double[2];
	

	const double a=10.0;
	const double b=28.0;
	const double c=2.666;

	runge_kutta ( x,y,z,n,dx, a,b,c);
    //file_make(x, y, z, n, step);
	delete [] x;
	delete [] y;
	delete [] z;

	getche();
    return 0;
}

void runge_kutta ( double *x, double *y,double *z, const int n, const double dx, const double a, const double b, const double c)
{
	x[0]=1;
	y[0]=1;
	z[0]=1;
	double kx[4];
	double ky[4];
	double kz[4];

	ofstream out("runge_kutta.txt");

	for(int i=1; i<n; i++)
	{
		// berechnung von x
		
		kx[0]=funx(x[i-1],y[i-1],a);
		kx[1]=funx(x[i-1]+dx*0.5*kx[0],y[i-1]+dx*0.5*kx[0],a);
		kx[2]=funx(x[i-1]+dx*0.5*kx[1],y[i-1]+dx*0.5*kx[1],a);
		kx[3]=funx(x[i-1]+dx*kx[2],y[i-1]+dx*kx[2],a);
		

		x[i]=x[i-1]+dx*0.1666*(kx[0]+2*kx[1]+2*kx[2]+kx[3]);

		//berechnung von y
		ky[0]=funy(x[i-1],y[i-1],z[i-1],a);
		ky[1]=funy(x[i-1]+dx*0.5*ky[0],y[i-1]+dx*0.5*ky[0],z[i-1]+dx*0.5*ky[0],a);
		ky[2]=funy(x[i-1]+dx*0.5*ky[1],y[i-1]+dx*0.5*ky[1],z[i-1]+dx*0.5*ky[1],a);
		ky[3]=funy(x[i-1]+dx*ky[2],y[i-1]+dx*ky[2],z[i-1]+dx*0.5*ky[2],a);
		
		y[i]=y[i-1]+dx*0.1666*(ky[0]+2*ky[1]+2*ky[2]+ky[3]);

		kz[0]=funz(x[i-1],y[i-1],z[i-1],c);
		kz[1]=funz(x[i-1]+dx*0.5*kz[0],y[i-1]+dx*0.5*kz[0],z[i-1]+dx*0.5*kz[0],c);
		kz[2]=funz(x[i-1]+dx*0.5*kz[1],y[i-1]+dx*0.5*kz[1],z[i-1]+dx*0.5*kz[1],c);
		kz[3]=funz(x[i-1]+dx*kz[2],y[i-1]+dx*kz[2],z[i-1]+dx*kz[2],c);
		
		z[i]=z[i-1]+dx*0.1666*(kz[0]+2*kz[1]+2*kz[2]+kz[3]);

		out<<i*dx<<" \t"<<x[i-1]<<" \t"<<y[i-1]<<" \t"<<z[i-1]<<endl;
		

	}
	out.close();
	
}


double funx(double x1,  double y1,const double a)
{
	return(a*(y1-x1));
}

double funy (double x1, double y1, double z1, const double b)
{
	return(x1*(b-z1)-y1);
}
double funz (double x1, double y1,double z1, const double c)
{
	
	return(x1*y1-c*z1);
}
