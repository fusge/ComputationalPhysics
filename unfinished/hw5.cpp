// hw5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "nr3.h"
#include "gaussj.h"
#include <vector>
#include <cmath>

void gaussjordan(int m, int n, double potential);

int main()
{
	gaussjordan(2,2,2);
	system("PAUSE");
	return 0;
}

void gaussjordan(int m, int n, double potential)
{
	int size = m*n + 1;
	MatDoub_IO A(size,size);
	MatDoub_IO b(size,1);

	//first we create a "map" of what the loops look like with a matrix.
	vector< vector<int> > loops(m, vector<int>(n));
	int count=0;
	for(int j=0;j<m;j++){
		for(int i=0;i<n;i++){
			loops[j][i] = count;
			count++;
		}
	}

	//we also make an empty matrix
	for(int j=0;j<size;j++)
		for(int i=0;i<size;i++)
			A[j][i]=0;


	//now the following program, although complicated, should work in theory. this will create a matrix that 
	//describes the relationship between the connections and the voltage by using the map we created earlier.
	for(int il=0;il<n;il++){
		for(int jl=0;jl<m;jl++){
			if(il==0){
				if(jl==0){
					A[loops[jl][il]][loops[jl][il+1]]=-1.0;
					A[loops[jl][il]][loops[jl+1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il]]=4.0;
					A[loops[jl][il]][size-1]=-1.0;
				}
				else if(jl==m-1){
					A[loops[jl][il]][loops[jl-1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il+1]]=-1.0;
					A[loops[jl][il]][loops[jl][il]]=4.0;
				}
				else{
					A[loops[jl][il]][loops[jl-1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il+1]]=-1.0;
					A[loops[jl][il]][loops[jl+1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il]]=4.0;
				}
			}
			else if(il==n-1){
				if(jl==0){
					A[loops[jl][il]][loops[jl][il-1]]=-1.0;
					A[loops[jl][il]][loops[jl+1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il]]=4.0;
					A[loops[jl][il]][size-1]=-1.0;
				}
				else if(jl==m-1){
					A[loops[jl][il]][loops[jl][il-1]]=-1.0;
					A[loops[jl][il]][loops[jl-1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il]]=4.0;
				}
				else{
					A[loops[jl][il]][loops[jl][il-1]]=-1.0;
					A[loops[jl][il]][loops[jl+1][il]]=-1.0;
					A[loops[jl][il]][loops[jl-1][il]]=-1.0;
					A[loops[jl][il]][loops[jl][il]]=4.0;
				}
			}
			else if(jl==0 && !(il==0 || il==n-1)){
				A[loops[jl][il]][loops[jl][il-1]]=-1.0;
				A[loops[jl][il]][loops[jl+1][il]]=-1.0;
				A[loops[jl][il]][loops[jl][il+1]]=-1.0;
				A[loops[jl][il]][size-1]=-1.0;
				A[loops[jl][il]][loops[jl][il]]=4.0;
			}
			else if(jl==m-1 && !(il==0 || il==n-1)){
				A[loops[jl][il]][loops[jl][il-1]]=-1.0;
				A[loops[jl][il]][loops[jl-1][il]]=-1.0;
				A[loops[jl][il]][loops[jl][il+1]]=-1.0;
				A[loops[jl][il]][loops[jl][il]]=4.0;
			}
			else{
				A[loops[jl][il]][loops[jl-1][il]]=-1.0;
				A[loops[jl][il]][loops[jl][il-1]]=-1.0;
				A[loops[jl][il]][loops[jl+1][il]]=-1.0;
				A[loops[jl][il]][loops[jl][il+1]]=-1.0;
				A[loops[jl][il]][loops[jl][il]]=4.0;
			}
		}
	}
	for(int i=0;i<n;i++)
		A[size-1][loops[0][i]]=-1.0;

	A[size-1][size-1]=2.0;

	for(int i=0;i<size-1;i++)
		b[i][0]=0;
	b[size-1][0]=potential;

	for(int j=0;j<size;j++){
		for(int i=0;i<size;i++){
			cout << A[j][i] << " ";
		}
		cout << endl;
	}

	gaussj(A,b);

	cout<< "The solution to this is: \n";
	for(int i=0;i<size;i++)
		cout << i << ":\t" << b[i][0] << endl;

}