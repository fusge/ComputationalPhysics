// hw4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "nr3.h"
#include "rk4.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>

double pi = 3.14159265359;

using namespace std;

//calling functions
void creatematrix17();
void jacobirelaxation17(int iterations);
void creatematrix129();
void jacobirelaxation129(int iterations);
void gausssidel17(int iterations);
void gausssidel129(int iterations);
int gaussianrelaxation(double beta, bool save);
void betacalc();
void surface_conditions();
void surface_oscillations(int iterations);


int main()
{
	surface_conditions();
	surface_oscillations(100);
	return 0;
}

void creatematrix17()
{
	//The following program will create a matrix for the cross inside
	//the zero potential box. the dimensions of the cross are determined
	//when the program is run in relation to the sie of the box and the 
	//number of cells it has.
	ofstream matrix;
	int rows=17;
	int cols=rows;
	int midpoint=rows/2;
	int lengthofcross = rows/8;
	vector<vector<double>> a(17, vector<double>(17));
	matrix.open("matrix17.dat");		
	for(int j=0;j<rows;j++){
		for(int i=0;i<cols;i++){
			if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
				a[j][i]=10;
			else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
				a[j][i]=10;
			else
				a[j][i]=0;
			matrix << a[j][i] << " ";
		}
		matrix << endl;
	}
	matrix.close();
}

void creatematrix129()
{
	//The following program will create a matrix for the cross inside
	//the zero potential box. the dimensions of the cross are determined
	//when the program is run in relation to the sie of the box and the 
	//number of cells it has.
	ofstream matrix;
	int rows=129;
	int cols=rows;
	int midpoint=rows/2;
	int lengthofcross = rows/8;
	vector<vector<double>> a(129, vector<double>(129));
	matrix.open("matrix129.dat");		
	for(int j=0;j<rows;j++){
		for(int i=0;i<cols;i++){
			if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
				a[j][i]=10;
			else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
				a[j][i]=10;
			else
				a[j][i]=0;
			matrix << a[j][i] << " ";
		}
		matrix << endl;
	}
	matrix.close();
}

void jacobirelaxation17(int iterations)
{
	//This program will read the created matrix and then compute
	//the solution using two matricies: the current matrix and
	//a temporary matrix.
	//j are rows and i are columns.
	vector< vector<double> > temp(17, vector<double>(17));
	vector< vector<double> > b(17, vector<double>(17));
	creatematrix17();
	ifstream matrix;
	matrix.open("matrix17.dat");
	while (!matrix.eof()){
        for(int j=0; j<17; j++){
            for(int i=0; i<17; i++){
                matrix >> b[j][i];
            }
        }
    }
	matrix.close();
	
	
	//now that we read the matrix we will now perform 
	//operations on it to find the solution
	int midpoint=17/2;
	int lengthofcross=17/10;
	for(int iter=0;iter<iterations;iter++){
		for(int j=0;j<17;j++){
			for(int i=0;i<17;i++){
				if(i==0 || i==16 || j==0 || j==16)
					temp[j][i]=b[j][i];
				else if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
					temp[j][i]=b[j][i];
				else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
					temp[j][i]=b[j][i];
				else
					temp[j][i]=0.25*(b[j][i+1]+b[j][i-1]+b[j+1][i]+b[j-1][i]);
			}
		}
		b=temp;
	}

	//save the solution onto a file
	ofstream sol;
	sol.open("matrix17sol.dat");
	for(int j=0;j<17;j++){
		for(int i=0;i<17;i++){
			sol << b[j][i] << " ";
		}
		sol << endl;
	}
	sol.close();


}

void jacobirelaxation129(int iterations)
{
	//This program will read the created matrix and then compute
	//the solution using two matricies: the current matrix and
	//a temporary matrix.
	//j are rows and i are columns.
	vector< vector<double> > temp(129, vector<double>(129));
	vector< vector<double> > b(129, vector<double>(129));
	creatematrix129();
	ifstream matrix;
	matrix.open("matrix129.dat");
	while (!matrix.eof()){
        for(int j=0; j<129; j++){
            for(int i=0; i<129; i++){
                matrix >> b[j][i];
            }
        }
    }
	matrix.close();
	
	
	//now that we read the matrix we will now perform 
	//operations on it to find the solution
	int midpoint=129/2;
	int lengthofcross=129/8;
	for(int iter=0;iter<iterations;iter++){
		for(int j=0;j<129;j++){
			for(int i=0;i<129;i++){
				if(i==0 || i==128 || j==0 || j==128)
					temp[j][i]=b[j][i];
				else if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
					temp[j][i]=b[j][i];
				else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
					temp[j][i]=b[j][i];
				else
					temp[j][i]=0.25*(b[j][i+1]+b[j][i-1]+b[j+1][i]+b[j-1][i]);
			}
		}
		b=temp;
	}

	//save the solution onto a file
	ofstream sol;
	sol.open("matrix129sol.dat");
	for(int j=0;j<129;j++){
		for(int i=0;i<129;i++){
			sol << b[j][i] << " ";
		}
		sol << endl;
	}
	sol.close();

}

void gausssidel17(int iterations)
{
	/* This program creates a matrix of length size 17 and 
	calculates the solution using gaussian-sidel method of 
	computation. j is row index and i is column index.
	*/

	//first let us create the matrix
	vector< vector<double> > a(17, vector<double>(17));
	creatematrix17();
	ifstream matrix;
	matrix.open("matrix17.dat");
	while (!matrix.eof()){
        for(int j=0; j<17; j++){
            for(int i=0; i<17; i++){
                matrix >> a[j][i];
            }
        }
    }
	matrix.close();

	//now we compute the solution. we do the calculations as we go on the matrix. 
	//The if statements make sure that we do not comput on the cross or boundaries.
	int midpoint=17/2;
	int lengthofcross=17/8;
	for(int iter=0;iter<iterations;iter++){
		for(int j=0;j<17;j++){
			for(int i=0;i<17;i++){
				if(i==0 || i==16 || j==0 || j==16)
					a[j][i]=0;
				else if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
					a[j][i]=10;
				else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
					a[j][i]=10;
				else
					a[j][i]=0.25*(a[j][i+1]+a[j][i-1]+a[j+1][i]+a[j-1][i]);
			}
		}
	}

	//save the solution onto a file.
	ofstream sol;
	sol.open("matrix17sidelsol.dat");
	for(int j=0;j<17;j++){
		for(int i=0;i<17;i++){
			sol << a[j][i] << " ";
		}
		sol << endl;
	}
	sol.close();

}

void gausssidel129(int iterations)
{
	/* This program creates a matrix of length size 17 and 
	calculates the solution using gaussian-sidel method of 
	computation. j is row index and i is column index.
	*/

	//first let us create the matrix
	vector< vector<double> > a(129, vector<double>(129));
	creatematrix129();
	ifstream matrix;
	matrix.open("matrix129.dat");
	while (!matrix.eof()){
        for(int j=0; j<129; j++){
            for(int i=0; i<129; i++){
                matrix >> a[j][i];
            }
        }
    }
	matrix.close();

	//now we compute the solution. we do the calculations as we go on the matrix. 
	//The if statements make sure that we do not comput on the cross or boundaries.
	int midpoint=129/2;
	int lengthofcross=129/8;
	for(int iter=0;iter<iterations;iter++){
		for(int j=0;j<129;j++){
			for(int i=0;i<129;i++){
				if(i==0 || i==128 || j==0 || j==128)
					a[j][i]=0;
				else if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
					a[j][i]=10;
				else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
					a[j][i]=10;
				else
					a[j][i]=0.25*(a[j][i+1]+a[j][i-1]+a[j+1][i]+a[j-1][i]);
			}
		}
	}

	//save the solution onto a file.
	ofstream sol;
	sol.open("matrix129sidelsol.dat");
	for(int j=0;j<129;j++){
		for(int i=0;i<129;i++){
			sol << a[j][i] << " ";
		}
		sol << endl;
	}
	sol.close();
}

int gaussianrelaxation(double beta, bool save)
{
	/*This program will arrive to the solution using both under and over relaxations.
	We can determine the value of beta from the begining, but we can also make it so that
	it is inputed into the program. The program will write the solution onto a file and 
	it will return the number of iterations neccessary to arrive to a stable solution.
	*/

	//first create a matrix for storing the solution.
	vector< vector<double> > a(129, vector<double>(129));
	vector< vector<double> > b(129, vector<double>(129));
	ifstream matrix;
	matrix.open("matrix129.dat");
	while (!matrix.eof()){
        for(int j=0; j<129; j++){
            for(int i=0; i<129; i++){
                matrix >> a[j][i];
            }
        }
    }
	matrix.close();

	//The following uses a weighted usage of old and new values using beta, which is inputed at the 
	//very beginning of the function. 
	int iterations=0;
	int midpoint=129/2;
	int lengthofcross=129/8;
	b=a;
	while(true){
		for(int j=0;j<129;j++){
			for(int i=0;i<129;i++){
				if(i==0 || i==128 || j==0 || j==128)
					a[j][i]=0;
				else if(j==midpoint && i>=(midpoint)-lengthofcross && i<=(midpoint)+lengthofcross)
					a[j][i]=10;
				else if(i==midpoint && j>=(midpoint)-lengthofcross && j<=(midpoint)+lengthofcross)
					a[j][i]=10;
				else
					a[j][i]=0.25*beta*(a[j][i+1]+a[j][i-1]+a[j+1][i]+a[j-1][i]) + (1-beta)*b[j][i];
			}
		}
		iterations = iterations+1;
		//We want the loop to stop when values inside the matrix do not change much.
		if(a[midpoint+5][midpoint+5]-b[midpoint+5][midpoint+5]<1E-4)
			break;
		b=a;
	}

	//The following piece of code checks to see if it should save it or not. this value is set as 
	//an input to the funciton where this variable is called save. It is either true of false.
	if(save)
	{
		ofstream sol;
		sol.open("matrix129relaxation.dat");
		for(int j=0;j<129;j++){
			for(int i=0;i<129;i++){
				sol << a[j][i] << " ";
			}
			sol << endl;
		}
		sol.close();
	}
	return iterations;
}
 
void betacalc()
{
	int relax;
	int dummy;
	clock_t t1,t2;
	ofstream relaxation;
	relaxation.open("relaxbeta.dat");
	double beta=0.5;
	while(beta<1.95)
	{
		dummy = gaussianrelaxation(beta, false);
		relaxation << dummy << endl;
		beta = beta+0.05;
	}
	relaxation.close();
	system("PAUSE");
}

void surface_conditions()
{
	//This function will create a 1000X1000 matrix.
	//The physical dimentions are in 1m by 1m, so each element is 
	ofstream matrix;
	int rows=100;
	int cols=rows;
	vector<vector<double>> a(100, vector<double>(100));
	matrix.open("matrix100.dat");		
	for(int j=0;j<rows;j++){
		for(int i=0;i<cols;i++){
			a[j][i]=0;
			matrix << a[j][i] << " ";
		}
		matrix << endl;
	}
	matrix.close();
}

void surface_oscillations(int iterations)
{
	/*
	The way this code works is that it first loads the information of the initial conditions on two matrices, then 
	it calculates the future state and records it onto a third matrix called Temp. When all the values of Temp are calculated,
	the temp is loaded onto b and what was on b is loaded into a. information on a is discarded as it is not neccessary for 
	computation. The program will save "snapshots" of the solutions onto files, which will be analyzed and put together into an
	animation.
	*/
	//the following code will load the initial conditions of the matrix onto matrix a.
	ifstream matrix;
	int rows=100;
	int cols=rows;
	vector<vector <double> > a(100,vector<double>(100));
	matrix.open("matrix100.dat");
	while(!matrix.eof())
		for(int j=0;j<rows;j++)
			for(int i=0;i<cols;i++)
				matrix >> a[j][i];
	matrix.close();

	//now apply the initial conditions related to velocity to find a second matrix, b. 
	//we also define x and y to be i*dx and j*dx, respectively.
	vector<vector <double> > b(100, vector<double>(100));
	double length=pi;		//The length of the square membrane
	double dx=length/rows;
	for(int j=0;j<rows;j++){
		for(int i=0;i<cols;i++){
			if(i==0 || i==99 || j==0 || j==99)
				b[j][i]=0;
			else
				b[j][i] = dx*i*(pi-dx*i)*pow((dx*j),2)*(pi-dx*j*2.0/3.0);
		}
	}

	//All that is left is to compute solutions by using the discretized version of the wave equation.
	//This will all be done while taking into account boundary condtions.
	//First we define the constants
	ofstream membrane;
	double csquared=1;			//This is the veolocity of the membrane.
	double cprimexsquared=2;	//this is the velocity of the simulation. it must be faster than the wave itself
	double cprimeysquared=2;
	string filename;
	vector< vector<double> > temp(100, vector<double>(100));
	for (int k=0;k<iterations;k++){
		for(int j=0;j<100;j++){
			for(int i=0;i<100;i++){
				if(i==0 || i==99 || j==0 || j==99)
					temp[j][i] = 0;
				else
					temp[j][i]=(b[j][i+1]+b[j][i-1]-2*b[j][i])/(cprimexsquared/csquared) + (b[j+1][i]+b[j-1][i]-2*b[j][i])/(cprimeysquared/csquared) - a[j][i]+2*b[j][i];
			}
		}
		//The following if statement is there to record snapshots of the membrane vibrating up and down. We wish to only have a few pictures, not all of it.
		//make sure that if you want to rename the file every time, then use itoa() and include the stdio.h and stdlib.h
		string name="membrane";
		a=b;
		b=temp;
		char buffer[30];
		if(k%10==0){
 			itoa (k,buffer,10);
			filename = name.append(buffer);
 			membrane.open(filename.append(".dat"));		
			for(int j=0;j<rows;j++){
				for(int i=0;i<cols;i++){
					membrane << temp[j][i] << " ";
				}
				membrane << endl;
			}
			membrane.close();
		}
	}
}