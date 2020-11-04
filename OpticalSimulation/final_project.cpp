// final_project.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <time.h>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include "nr3.h"
#include "rk4.h"
using namespace std;

double cnot=1;
double munot=1;                  //1.2566370614E-6;

//Call functions before main
void Second_order_nonlinear(double w1, double w2, double w3, double a1, double a2, double a3, double d,int time);
void Second_order_nonlinear_alt(double w1, double w2, double w3, double a1, double a2, double a3, double d);

int main()
{
	double w1=0.02;
	double w2=w1;
	double w3=0.04;
	Second_order_nonlinear(w1,w2,w3,1,1,0,0.09,5000);
	
	return 0;
}

void Second_order_nonlinear(double w1, double w2, double w3, double a1, double a2, double a3, double d, int time)
{
	//first determine the constants related to the second harmonic generation. 
	//In this case it would be the constant d and index of refraction. We will use Gaussian units.
	double n=1.5;
	double cn = cnot/n;			//speed of light in the medium
	//simulation speed will be called cprime
	double cprime=5*cn;
	double deltat = 1;

	/*
	Now create variables that will describe the space and time variabls for each of the waves.
	the space variable is the column element, while the time variable is the row element. We will have the
	matrix empty, but there will be a source at the edge of the medium. In other words,
	there is not field until we begin the simulation, than this field moves into the medium as we run 
	the simulation
	*/
	vector< double > a(1000);
	vector< double > b(1000);
	vector< double > c(1000);
	vector< double > apast(1000);
	vector< double > bpast(1000);
	vector< double > cpast(1000);
	vector< double > atemp(1000);
	vector< double > btemp(1000);
	vector< double > ctemp(1000);

	//Determint the initial conditions. The incident wave is a, and the additional two waves are b and c. 
	//b and c are empty, and a is driven by a sine function. Assign zeros to all the other cells.
	for(int i=0;i<1000;i++){
		apast[i]=0;
		a[i]=0;
		bpast[1]=0;
		b[i]=0;
		cpast[i]=0;
		c[i]=0;
		atemp[i]=0;
		btemp[i]=0;
		ctemp[i]=0;
	}

	//Simulate the motion of the wave
	ofstream filea, fileb, filec;
	filea.open("a.dat");
	fileb.open("b.dat");
	filec.open("c.dat");
	for(int j=0;j<time;j++){
		for(int i=0;i<1000;i++){
			if(i==0){
				atemp[i]=a1*sin(w1*j);
				btemp[i]=a2*sin(w2*j);
				ctemp[i]=a3*sin(w3*j);
			}
			else if(!((i+1)<1000)){
				atemp[i]=atemp[i-1];
				btemp[i]=btemp[i-1];
				ctemp[i]=ctemp[i-1];
			}
			else{
				for(int iter=0;iter<5;iter++){
					atemp[i] = (a[i+1]+a[i-1]-2*a[i])*cn*cn/(cprime*cprime) - d*2*munot*cn*cn*(atemp[i]*atemp[i]+btemp[i]*ctemp[i]+apast[i]*apast[i]+bpast[i]*cpast[i]-2*a[i]*a[i]-2*b[i]*c[i]) - apast[i] + 2*a[i];
					btemp[i] = (b[i+1]+b[i-1]-2*b[i])*cn*cn/(cprime*cprime) - d*2*munot*cn*cn*(btemp[i]*btemp[i]+atemp[i]*ctemp[i]+bpast[i]*bpast[i]+apast[i]*cpast[i]-2*b[i]*b[i]-2*a[i]*c[i]) - bpast[i] + 2*b[i];
					ctemp[i] = (c[i+1]+c[i-1]-2*c[i])*cn*cn/(cprime*cprime) - d*2*munot*cn*cn*(ctemp[i]*ctemp[i]+atemp[i]*btemp[i]+cpast[i]*cpast[i]+apast[i]*bpast[i]-2*c[i]*c[i]-2*a[i]*b[i]) - cpast[i] + 2*c[i];
					//btemp[i] = (b[i+1]+b[i-1]-2*b[i])*cn*cn/(cprime*cprime) + d*2*munot*w2*w2*a[i]*c[i]*deltat*deltat*cn*cn - bpast[i] + 2*b[i];
					//ctemp[i] = (c[i+1]+c[i-1]-2*c[i])*cn*cn/(cprime*cprime) + d*2*munot*w3*w3*b[i]*a[i]*deltat*deltat*cn*cn - cpast[i] + 2*c[i];
				}
			}
			
		}
		if(j%10==0){
			for(int i=0;i<1000;i++){
				filea << atemp[i] << " ";
				fileb << btemp[i] << " ";
				filec << ctemp[i] << " ";
			}
			filea << endl;
			fileb << endl;
			filec << endl;
		}

		apast=a;
		a=atemp;
		bpast=b;
		b=btemp;
		cpast=c;
		c=ctemp;
	}
	filea.close();
	fileb.close();
	filec.close();

}


void Second_order_nonlinear_alt(double w1, double w2, double w3, double a1, double a2, double a3, double d)
{
	//first determine the constants related to the second harmonic generation. 
	//In this case it would be the constant d and index of refraction. We will use Gaussian units.
	double n=1.5;
	double cn = cnot/n;			//speed of light in the medium
	//simulation speed will be called cprime
	double cprime=1.2*cn;
	double deltax=0.0000000001;

	cout << "begin" << endl;
	vector<double> a(500);
	vector<double> b(500);
	vector<double> c(500);
	vector<double> k(3);
	k[0] = w1/cn;
	k[1] = w2/cn;
	k[2] = w3/cn;
	cout << "just getting started" << endl;
	/*
	vector<double> atemp(500);
	vector<double> btemp(500);
	vector<double> ctemp(500);
	vector<double> apast(500);
	vector<double> bpast(500);
	vector<double> cpast(500);
	*/

	for(int i=0;i<500;i++){
		/*
		apast[i]=0;
		bpast[i]=0;
		cpast[i]=0;
		*/
		a[i]=0;
		b[i]=0;
		c[i]=0;
	}

	cout << "creation complete" << endl;

	for(int i=0;(i+1)<500;i++){
		cout << i << endl;
		if(i==0){
			a[i]=a1;
			b[i]=a2;
			c[i]=a3;
			cout << "first stage complete" << endl;
		}
		/*
		else if(!((i+1)<500)){
			a[i+1]=a[i];
			a[i+1]=a[i];
			a[i+1]=a[i];
			cout << "last stage complete" << endl;
		}*/
	
		else{
			for(int iterations=0;iterations<3;iterations++){
				a[i+1] = 2*a[i] - a[i-1] - 2*munot*w1*w1*d*b[i]*c[i]*deltax*deltax - k[0]*k[0]*a[i]*deltax*deltax;
				b[i+1] = 2*b[i] - b[i-1] - 2*munot*w2*w2*d*a[i]*c[i]*deltax*deltax - k[1]*k[1]*b[i]*deltax*deltax;
				c[i+1] = 2*c[i] - c[i-1] - 2*munot*w3*w3*d*a[i]*b[i]*deltax*deltax - k[2]*k[2]*c[i]*deltax*deltax;
			}
		}
	}

	ofstream filea, fileb, filec;
	filea.open("a.dat");
	fileb.open("b.dat");
	filec.open("c.dat");
	for(int i=0;i<500;i++){
		filea << a[i] << " ";
		fileb << b[i] << " ";
		filec << c[i] << " ";
	}

	filea.close();
	fileb.close();
	filec.close();
}
