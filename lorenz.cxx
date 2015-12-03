#include <iostream>
#include <fstream>


void fun(double* ki, double x,double y,double z,const int a, const int b, const int c );
void delet_temp(double* temp);


using namespace std;

int main()
{
	
	const double dx=0.01;
	const double T=100;
	const int n=10000;
	const double a=10.0;
	const double b=28.0;
	const double c=2.666;
	double *x= new double [n];
	double *y= new double [n];
	double *z= new double [n];
	double k1[3],k2[3],k3[3],k4[3];
	double temp[3];
	x[0]=1;
	y[0]=1;
	z[0]=1;
	ofstream out("runge_kutta3.txt");
	for (int i=1; i<n; i++)
	{
		//k1-bestimmen
		fun(k1, x[i-1],y[i-1],z[i-1], a,  b, c );
		//k2-bestimmen
		for(int j=0; j<3; j++)
		{
			fun(temp,x[i-1]+dx*k1[j]*0.5, y[i-1]+dx*k1[j]*0.5, z[i-1]+dx*k1[j]*0.5,a, b, c);
			k2[j]=temp[j];
			delet_temp(temp);
		}
		//k3-bestimmen
		for(int j=0; j<3; j++)
		{
			fun(temp,x[i-1]+dx*k2[j]*0.5, y[i-1]+dx*k2[j]*0.5, z[i-1]+dx*k2[j]*0.5,a, b, c);
			k3[j]=temp[j];
			delet_temp(temp);
		}
		//k4-bestimmen
		for(int j=0; j<3; j++)
		{
			fun(temp,x[i-1]+dx*k3[j]*0.5, y[i-1]+dx*k3[j]*0.5, z[i-1]+dx*k3[j]*0.5,a, b, c);
			k4[j]=temp[j];
			delet_temp(temp);
		}
		
		x[i]=x[i-1]+dx/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
		y[i]=y[i-1]+dx/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
		z[i]=z[i-1]+dx/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
		out<<(i-1)*dx<<" \t "<< x[i-1] <<" \t "<< y[i-1]<<" \t "<<z[i-1]<<endl;
	}
	delete [] x;
	delete [] y;
	delete [] z;

	
    return 0;
}


void fun(double* ki, double x,double y,double z,const int a, const int b, const int c ) //f[0]->x, f[1]->y, f[2]->z
{
	ki[0]=a*(y-x);
	ki[1]=x*(b-z)-y;
	ki[2]=(x*y-c*z);
}
void delet_temp(double* temp)
{
	for(int i=0; i<3; i++)
		temp[i]=0;
}
