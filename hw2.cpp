// hw2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <time.h>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <fstream>
using namespace std;

//some functions that we will use
double randomdouble(double low, double high);
double kmomentrandom(int moment);
double autocorrelation(int separation);
double integrationrand(int p);
double funcsquared(int p);
double error(int p);

int main()
{

	cout << error(6) << endl;
	system("PAUSE");
}

//The program below will return a number between the 
//low and high parameters.
double randomdouble(double low, double high)
{
	//this function will create random numbers of the double type
	double temp;
	//this calculates random number and returns it
	temp = ((double) rand() / (static_cast<double>(RAND_MAX) + 1.0)) * (high-low) + low;
	return temp;
}

//this will return the kth moment for the random list.
double kmomentrandom(int moment)
{
	srand((unsigned)time(0));
	double y = 0;
	double size = 100000.0;
	for (int k = 0; k < size; k++)
	{
		y = y + pow(randomdouble(0,1),moment);
	}
	double value = y/size;
	return value;
}

//This function returns a value representing a correlation between 
//points in a random list. this will depend on the separation of 
//their indicies.
double autocorrelation(int separation)
{
	srand((unsigned)time(0));
	int k = separation;
	double x = 0;
	double y[100000];
	double size = 100000.0;
	for (int k = 0; k < size; k++)
		y[k] = randomdouble(0,1);
	for(int i=0; i+k <size; i++)
		x = x + y[i]*y[i + k];
	double value = x/size;
	return value;
}
//below can be used with the autocorelation function
/*ofstream correlate;
	correlate.open("autocorrelation.dat");
	for(int i = 0; i<6; i++)
		correlate << autocorrelation(i) << endl;
	system("PAUSE");
	correlate.close();
	return 0;*/


double integrationrand(int p)
{
	double size = pow(10.0,p);
	double y=0;
	
	srand((unsigned)time(0));
	double x;
	for(int i = 0; i < size; i++)
	{
		x = randomdouble(0,1)+randomdouble(0,1)+randomdouble(0,1)+randomdouble(0,1)+randomdouble(0,1);
		y = y + pow(x,3);
	}
	double val = y/size;
	return val;
}
/*for(int i=1; i < 9; i++)
		cout << integrationrand(i) << endl;
*/

double funcsquared(int p)
{
	double size = pow(10.0,p);
	double y=0;
	
	srand((unsigned)time(0));
	double x;
	for(int i = 0; i < size; i++)
	{
		//x = pow(randomdouble(0,1),2)+pow(randomdouble(0,1),2)+pow(randomdouble(0,1),2)+pow(randomdouble(0,1),2)+pow(randomdouble(0,1),2);
		x = randomdouble(0,1)+randomdouble(0,1)+randomdouble(0,1)+randomdouble(0,1)+randomdouble(0,1);
		y = y + pow(x,5);
	}
	double val = y/size;
	return val;
}

double error(int p)
{
	double sigmasquared;
	double sigma;
	double size = pow(10.0,p);
	double y = 0;
	double x = 0;
	double funct1,funct2;
	srand((unsigned)time(0));
	double x1,x2,x3,x4,x5;
	for(int i = 0; i<size; i++)
	{
		x1=randomdouble(0,1);
		x2=randomdouble(0,1);
		x3=randomdouble(0,1);
		x4=randomdouble(0,1);
		x5=randomdouble(0,1);
		x = x + pow(x1 + x2 + x3 + x4 + x5,5);
		y = y + pow(x1 + x2 + x3 + x4 + x5,3);
	}
	sigmasquared = x/size - pow(y/size,2);
	sigma = sqrt(abs(sigmasquared));
	return sigma;
}
