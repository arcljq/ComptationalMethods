

#include <iostream>
#include <cmath> 
#include <fstream>
#include <vector>



#include "Methods.h"

using namespace std;

ofstream outfile("result.dat");


//the follow are given data
const double T_sur = 149.0;  
const double T_in = 38.0;
const double L = 31.0;
const double lastt = 0.5;
const double D = 39.0;
const double Pi = 3.14;

// variables
double dx;
double dt;
int orderofmethod;//the number of method (designed for selection function)
double belta1;
double belta2;
double belta3;
double belta4;
double m1;

// This is the part of different methods
void gridInitialization(int rows, int cols, double** A0,double** AMO);
void gridinputDFF(double** A0);
void gridinputRCCE(double** A0);
void gridinputLSIFCI(double** A0);
void gridinputCNT(double** A0);



int main(int argc, char* argv[])
{

	// input dx and dt;
	cout << "please input the dx :" << endl;
	cin >> dx;
	
	cout << "please input the dt :" << endl;
	cin >> dt;		

	//This is Grid Initialization
	double j = L / dx;
	double k = lastt / dt;
	int rows = k, cols = j;
	double** A0 =0 ;//this is for time and space
	double** AMO = 0;//this is for coefficient for Thomas Algorithm

	gridInitialization(rows, cols, A0, AMO);


	

	
	//This is set max m for the point Initialization
	cout << "please input the Upper limit of m :" << endl;
	cin >> m1;

	// chose a kind of different methods
	cout << "please chose a kind of method" <<endl;
	cout << "(pleas input the number of method. "<< endl;
	cout << "0:DuFort_Frankel; 1:Richardson; 2:Laasonen Simple Implict; 3:Crank-Nicholson.)" << endl;
	cin >> orderofmethod;
	
	/*switch (orderofmethod) 
	{
		case 0://chose DFF method

		case 1://chose RCCE method

		case 2://chose LSIFCI method

		case 3://chose CNT method

	}*/

	return 0;
}


// This is Grid Initialization
void gridInitialization(int rows, int cols, double** A0, double** AMO)
{
	A0 = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		A0[i] = new double[cols];                          
	}

	for (int j = 1; j < cols; j++)
	{
		A0[0][j] = T_in;
	}
	
	for (int i = 0; i < rows; i++)
	{
		A0[i][0] = T_sur;
	}

	AMO = new double* [cols];
	for (int i = 0; i < cols; i++)
	{
		AMO[i] = new double[cols];
	}
}

void gridinputDFF(double** A0) 
{
	//Firstly, caculating the t=1*dt by using different methods(for comparing) 
	//after all ,we chose the CNT method to caculate the t=1*dt
}

void gridinputRCCE(double** A0)
{

}

void gridinputLSIFCI(double** A0)
{

}

void gridinputCNT(double** A0)
{

}

