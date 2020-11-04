// hw3.cpp : Defines the entry point for the console application.
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

//the declarations of the constants below are for the planetary problems
double mg = 1.3216e8;
double mgj = 317.8*mg;

double earthmass = 5.97E24;
double jupitermass =34343433;
double sunmass = 3.33E5;
double gconstant = 6.67398E-11;     //we should use gmsun instead. remember to replace them in the other functions.
double earthsemimajor = 1.52E11;
double earthminvelocity = 29.29E3;
double jupitersemimajor = 5.204;
double jupiterminvelocity = 0.425;
double mercurymass = 0.0553;
double mercurysemimajor = 0.459*earthsemimajor;
double mercuryminvelocity = 1.327*earthminvelocity;
double au = 1.496E11;
double alpha = 0.0008*pow(au,2);
double gmsun = 1.328E20;
//-pow(0.007297,0.5)
//the following declarations involve the third problem
//the units of mass for these values are in atomic scale("natural" lorentz-heaviside units)
double bigq = 1.0;
double smallq = -1.0;
double electronmass = 100.0;
double c = 1.0;
double edistance = 5.29;		//units for this distance is picometres.
double bx = 0;
double by = 0;
double bz = 0;

//now comes the deriv declarations
double randomdouble(double low, double high)
{
	//this function will create random numbers of the double type
	double temp;
	//this calculates random number and returns it
	temp = ((double) rand() / (static_cast<double>(RAND_MAX) + 1.0)) * (high-low) + low;
	return temp;
}
void derivsearth(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	double val = pow(y[0]*y[0] + y[2]*y[2],-1.5);
	dydx[0] = y[1];
	dydx[1] = -mg*y[0]*val;
	dydx[2] = y[3];
	dydx[3] = -mg*y[2]*val;
}
void derivsjupiter(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	double val = pow(y[0]*y[0] + y[2]*y[2],-1.5);
	dydx[0] = y[1];
	dydx[1] = -mgj*y[0]*val;
	dydx[2] = y[3];
	dydx[3] = -mgj*y[2]*val;
}
void derivsinteract(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	double vale = pow(y[0]*y[0] + y[2]*y[2],-1.5);
	double valearth = pow((y[0]-y[4])*(y[0]-y[4]) + (y[2]-y[6])*(y[2]-y[6]),-1.5);
	double valj = pow(y[4]*y[4] + y[6]*y[6],-1.5);
	double valjupiter = pow((y[4]-y[0])*(y[4]-y[0]) + (y[6]-y[2])*(y[6]-y[2]),-1.5);
	dydx[0] = y[1];
	dydx[1] = -gconstant*sunmass*y[0]*vale - gconstant*jupitermass*(y[0]-y[4])*valearth;
	dydx[2] = y[3];
	dydx[3] = -gconstant*sunmass*y[2]*vale - gconstant*jupitermass*(y[2]-y[6])*valearth;
	dydx[4] = y[5];
	dydx[5] = -gconstant*sunmass*y[4]*valj - gconstant*earthmass*(y[4]-y[0])*valjupiter;
	dydx[6] = y[7];
	dydx[7] = -gconstant*sunmass*y[6]*valj - gconstant*earthmass*(y[6]-y[2])*valjupiter;
}
void derivsmercury(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	double valm = pow(y[0]*y[0] + y[2]*y[2],-1.5);
	double valmalpha = pow(y[0]*y[0] + y[2]*y[2],-1);
	dydx[0] = y[1];
	dydx[1] = -gmsun*y[0]*valm*(1+alpha*valmalpha);
	dydx[2] = y[3];
	dydx[3] = -gmsun*y[2]*valm*(1+alpha*valmalpha);

}
double gamma = 0.00001;
void derivlorentz(const Doub x, VecDoub_I & y, VecDoub_O & dydx)
{
	double xyzr = pow(y[0]*y[0]+y[2]*y[2]+y[4]*y[4],-0.5);
	double xyr = pow(y[0]*y[0]+y[2]*y[2],-0.5);
	double sintheta = y[2]*xyr;
	double costheta = y[0]*xyr;
	double sinphi = xyzr*1/xyr;
	double cosphi = y[4]*xyzr;
	double esubx = bigq*xyzr*sinphi*costheta;
	double esuby = bigq*xyzr*sinphi*sintheta;
	double esubz = bigq*xyzr*cosphi; 
	double qoverm = smallq/electronmass;
	//below is a piece of code that changes the B field to act like B charges exist.
	//I also comented out the magentic field constants. this can be put back in if you 
	//wish to have constant field. this was replaced by the random B field generator introduced
	//below.
	/*
	double bconst = 0.1;
	double bsubx = bconst*xyzr*sinphi*costheta;
	double bsuby = bconst*xyzr*sinphi*sintheta;
	double bsubz = bconst*xyzr*cosphi;
	
	double bsubx = 0;
	double bsuby = -0.1;
	double bsubz = -0.001;
	*/
	
	double bsubx = bx + gamma*randomdouble(-0.5,0.5);
	double bsuby = by + gamma*randomdouble(-0.5,0.5);
	double bsubz = bz + gamma*randomdouble(-0.5,0.5);
	bx = bsubx;
	by = bsuby;
	bz = bsubz;

	dydx[0] = y[1];
	dydx[1] = qoverm*esubx + pow(c,-1)*qoverm*(bsubz*y[3]-bsuby*y[5]);
	dydx[2] = y[3];
	dydx[3] = qoverm*esuby + pow(c,-1)*qoverm*(bsubz*y[1]-bsubx*y[5]);
	dydx[4] = y[5];
	dydx[5] = qoverm*esubz + pow(c,-1)*qoverm*(bsuby*y[1]-bsubx*y[3]);
}


void earthorbit();
void jupiterorbit();
void orbitinteraction();
void mercuryprecession();
void chargedparticlemotion();


//here is main
int main()
{
	chargedparticlemotion();
	return 0;
}

void earthorbit()
{
	VecDoub y(4), dydx(4);
	Doub x, xmin, xmax, kmax = 10000000000, h = 100;
	VecDoub yout(4);
	int k;
	xmin = 1.; xmax = 100000.;
	y[0] = 152;
	y[1] = 0;
	y[2] = 0;
	y[3] = 924;
	derivsearth(xmin, y, dydx);
	ofstream diff;
	diff.open("earth.dat");
	for(k=0;k<kmax;k++)
	{
		x=xmin+k*h;
		rk4(y, dydx, x, h, yout, derivsearth);
		if (k%500 == 0)
			diff << yout[0] << " " << yout[2] << endl;
		y[0]=yout[0];
		y[1]=yout[1];
		y[2]=yout[2];
		y[3]=yout[3];
		derivsearth(x,y,dydx);
	}
	diff.close();

}

void jupiterorbit()
{
	VecDoub y(4), dydx(4);
	Doub x, xmin, xmax, kmax = 100000, h = 0.0001;
	VecDoub yout(4);
	int k;
	xmin = 1.; xmax = 100000.;
	y[0] = -147*4.95;
	y[1] = 0;
	y[2] = 0;
	y[3] = -405.53;
	derivsearth(xmin, y, dydx);
	ofstream diff;
	diff.open("jupiter.dat");
	for(k=0;k<kmax;k++)
	{
		x=xmin+k*h;
		rk4(y, dydx, x, h, yout, derivsearth);
		if (k%100 == 0)
			diff << yout[0] << " " << yout[2] << endl;
		y[0]=yout[0];
		y[1]=yout[1];
		y[2]=yout[2];
		y[3]=yout[3];
		derivsearth(x,y,dydx);
	}
	diff.close();

}

void orbitinteraction()
{
	VecDoub y(8), dydx(8);
	//we need to make sure that h is big enough. This represents the number of seconds in a day.
	Doub x, xmin, xmax, kmax = 1820000, h = 3600;
	VecDoub yout(8);
	int k;
	xmin = 1.; xmax = 100000;
	//This sets the initial conditions for the formulation of the simulation.
	y[0] = earthsemimajor;
	y[1] = 0;
	y[2] = 0;
	y[3] = earthminvelocity;
	y[4] = jupitersemimajor;
	y[5] = 0;
	y[6] = 0;
	y[7] = jupiterminvelocity;
	derivsinteract(xmin, y, dydx);
	ofstream diff;
	diff.open("orbitsinteraction.dat");
	for(k=0;k<kmax;k++)
	{
		x=xmin+k*h;
		rk4(y, dydx, x, h, yout, derivsinteract);
		if (k%10 == 0)
			diff << yout[0] << " " << yout[2] << " " << yout[4] << " " << yout[6] << endl;
		y[0]=yout[0];
		y[1]=yout[1];
		y[2]=yout[2];
		y[3]=yout[3];
		y[4]=yout[4];
		y[5]=yout[5];
		y[6]=yout[6];
		y[7]=yout[7];
		derivsinteract(x,y,dydx);
	}
	diff.close();

}

void mercuryprecession()
{
	VecDoub y(4), dydx(4);
	Doub x, xmin, xmax, kmax = 400000, h = 60;
	VecDoub yout(4);
	int k;
	xmin = 1.; xmax = 100000.;
	y[0] = mercurysemimajor;
	y[1] = 0;
	y[2] = 0;
	y[3] = mercuryminvelocity;
	derivsmercury(xmin, y, dydx);
	ofstream diff, axis;
	diff.open("mercuryprecession.dat");
	for(k=0;k<kmax;k++)
	{
		x=xmin+k*h;
		rk4(y, dydx, x, h, yout, derivsmercury);
		if (k%10 == 0)
			diff << yout[0] << " " << yout[2] << endl;
		y[0]=yout[0];
		y[1]=yout[1];
		y[2]=yout[2];
		y[3]=yout[3];
		derivsmercury(x,y,dydx);
	}
	diff.close();
}

void chargedparticlemotion()
{
	srand((unsigned)time(0) );
	VecDoub y(6), dydx(6);
	Doub x, xmin, xmax, kmax = 100000, h = 0.1;
	VecDoub yout(6);
	int k;
	xmin = 1.; xmax = 100000.;
	y[0] = edistance;
	y[1] = 0;
	y[2] = 0;
	y[3] = 1;
	y[4] = 0;
	y[5] = 0;
	derivlorentz(xmin, y, dydx);
	ofstream diff;
	diff.open("electron.dat");
	for(k=0;k<kmax;k++)
	{
		x=xmin+k*h;
		rk4(y, dydx, x, h, yout, derivlorentz);
		if (k%100 == 0)
			diff << yout[0] << " " << yout[2] << " " << y[4] << endl;
		y[0]=yout[0];
		y[1]=yout[1];
		y[2]=yout[2];
		y[3]=yout[3];
		y[4]=yout[4];
		y[5]=yout[5];
		derivlorentz(x,y,dydx);
	}
	diff.close();

}