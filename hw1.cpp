// hw1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
using namespace std;
#include <cmath>
#include <fstream>
#include <iomanip>

//forward declarations
void integrate();
void forwrddiff();
void centraldiff();
double epsilon();
void fivepointstencil();
void bisection();

double pi = 3.14159265;



int main()
{
	double y;
	bisection();
	system("PAUSE");
	return 0;
}

void integrate()
{
	//trapezoidal rule for integrating functions
	double y = 0.0;
	double x = 0.0;
	double step = 0.0009765625;
	int counter = 0;
	double upperbound = pi/2;
	while (true)
	{
		if(counter = 0)
		{
			y = step*cos(x)/2;
			counter += 1;
		}
		else if(x<upperbound && (x+step>upperbound) || x==upperbound)
		{
			y = y + cos(x)*step/2;
			break;
		}
		else if(x>upperbound)
			break;
		else
		{
			y = y + cos(x)*step;
			x = x + step;
		}
	}
	cout << y << "\t" << (y-1) << endl;


	//now this is simpson's method of integration
	y = 0.0;
	double b = step;
	double a = 0.0;
	while (a<upperbound)
	{ 
		y = y + (b-a)/6 * (cos(a) + 4*cos((a+b)/2) + cos(b));
		a = b;
		b += step;
	}
	cout << y << "\t" << (y-1) << endl;

	//and finally the Gaussian Quadrature for 3 points
	double c[] = {0.555555556, 0.888888889, 0.555555556};
	double xpar[] = {-0.774596669, 0.000000000, 0.774596669};
	b = pi/2;
	a = 0.0;
	y = 0.0;
	for(int i=0; i<3; i++)
			y = y + c[i]*cos((b-a)/2*xpar[i]+(b+a)/2)*(b-a)/2;
	cout << y << "\t" << (y-1) << endl;

}

double epsilon()
{
	double y, epsilon2;

	y = 1;
	while (y+1>1)
	{
		epsilon2 = y;
		y /= 2;
	}
	return epsilon2;
}

void forwrddiff()
{
	double x = 1.0;
	int counter = 0;
	double pointdiff = 2;
	double h = 1;
	double diff;
	ofstream fwrderror;
	fwrderror.open("forward1.dat");

	while(true)
	{
		if (counter>150)
			break;
		diff = (pow(x+h,2)-pow(x,2))/h;
		if (diff==0)
			break;
		fwrderror << setprecision(20) << diff << "\t" << "\t" << abs(diff-2)/2 << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	fwrderror.close();

	counter = 0;
	h = 1;
	fwrderror.open("forward2.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (pow(x+h,3)-pow(x,3))/h;
		if (diff==0)
			break;
		fwrderror << setprecision(20) << diff << "\t" << "\t" << abs(diff-3)/3 << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	fwrderror.close();

	counter = 0;
	h = 1;
	double actval = exp(1.0);
	fwrderror.open("forward3.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (exp(-x-h)-exp(-x))/h;
		if (diff==0)
			break;
		fwrderror << setprecision(20) << diff << "\t" << "\t" << abs(diff-actval)/actval << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	fwrderror.close();

	counter = 0;
	h = 1;
	x = pi/6.0;
	actval = -0.5;
	fwrderror.open("forward4.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (cos(x+h)-cos(x))/h;
		if (diff==0)
			break;
		fwrderror << setprecision(20) << diff << "\t" << "\t" << abs((diff-actval)/actval) << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	fwrderror.close();


}

void centraldiff()
{
	double x = 1.0;
	int counter = 0;
	double h = 1.0;
	double diff;
	ofstream centralerror;
	centralerror.open("central1.dat");

	while(true)
	{
		if (counter>150)
			break;
		diff = (pow(x+h/2,2)-pow(x-h/2,2))/h;
		if (diff==0)
			break;
		centralerror << setprecision(20) << diff << "\t" << "\t" << abs(diff-2)/2 << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	centralerror.close();

	counter = 0;
	h = 1.0;
	centralerror.open("central2.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (pow(x+h/2,3)-pow(x-h/2,3))/h;
		if (diff==0)
			break;
		centralerror << setprecision(20) << diff << "\t" << "\t" << abs(diff-3)/3 << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	centralerror.close();

	counter = 0;
	h = 1.0;
	double actval = exp(1.0);
	centralerror.open("central3.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (exp(x+h/2)-exp(x-h/2))/h;
		if (diff==0)
			break;
		centralerror << setprecision(20) << diff << "\t" << "\t" << abs(diff-actval)/actval << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	centralerror.close();


	counter = 0;
	h = 1.0;
	x = pi/6.0;
	actval = -0.5;
	centralerror.open("central4.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (cos(x+h/2)-cos(x-h/2))/h;
		if (diff==0)
			break;
		centralerror << setprecision(20) << diff << "\t" << "\t" << abs((diff-actval)/actval) << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	centralerror.close();


}

void fivepointstencil()
{
	double x = 1.0;
	int counter = 0;
	double h = 1.0;
	double diff;
	ofstream stencilerror;
	stencilerror.open("stencil1.dat");

	while(true)
	{
		if (counter>150)
			break;
		diff = (-pow(x+2*h,2) + 8*pow(x+h,2) - 8*pow(x-h,2) + pow(x-2*h,2))/(12*h);
		if (diff==0)
			break;
		stencilerror << setprecision(20) << diff << "\t" << "\t" << abs(diff-2)/2 << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	stencilerror.close();

	counter = 0;
	h = 1.0;
	stencilerror.open("stencil2.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (-pow(x+2*h,3) + 8*pow(x+h,3) - 8*pow(x-h,3) + pow(x-2*h,3))/(12*h);
		if (diff==0)
			break;
		stencilerror << setprecision(20) << diff << "\t" << "\t" << abs(diff-3)/3 << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	stencilerror.close();

	counter = 0;
	h = 1.0;
	double actval = exp(1.0);
	stencilerror.open("stencil3.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (-exp(x+2*h) + 8*exp(x+h) - 8*exp(x-h) + exp(x-2*h))/(12*h);;
		if (diff==0)
			break;
		stencilerror << setprecision(20) << diff << "\t" << "\t" << abs(diff-actval)/actval << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	stencilerror.close();

	counter = 0;
	h = 1.0;
	x = pi/6.0;
	actval = -0.5;
	stencilerror.open("stencil4.dat");
	while(true)
	{
		if (counter>150)
			break;
		diff = (-cos(x+2*h) + 8*cos(x+h) - 8*cos(x-h) + cos(x-2*h))/(12*h);;
		if (diff==0)
			break;
		stencilerror << setprecision(20) << diff << "\t" << "\t" << abs((diff-actval)/actval) << "\t" << h << '\n';
		h = h*0.75;
		counter += 1;
	}
	stencilerror.close();
}

void bisection()
{
	double leftbound = 0.5;
	double rightbound = 2;
	double x;
	double eps = epsilon();
	int count = 0;
	while(true)
	{
		x = (leftbound-rightbound)/2;
		if (abs(cos(x))<eps)
			break;
		else if(cos(x)>0)
			leftbound = x;
		else if(cos(x)<0)
			rightbound = x;
		count += 1;
	}
	cout << x << endl << count <<  endl;
}