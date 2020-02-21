#pragma once
#include <cmath>

class Methods
{
public:

	void DFF1(double D, double dx, double dt, double** A1, int rows, int cols);
	void DFF2(double D, double dx, double dt, double** A1, int rows, int cols);
	void DFF3(double D, double dx, double dt, double** A1, int rows, int cols);
	void DFF4(double D, double dx, double dt, double** A1, int rows, int cols);
	void DFF5(double D, double dx, double dt, double** A1, int rows, int cols,int m1);

	void RCCE(double D, double** AT, double dx, double dt, double T_sur, double T_in, double** A1, int rows, int cols);

	void LSIFCI(double D, double** AT,  double dx, double dt, double T_sur, double T_in, double** A1, int rows, int cols);

	void CNT(double D, double** AT,   double dx, double dt, double** A0, int rows, int cols);

};

//This Method is using Crank-Nicholson (Trapezoidal) Method to calculate t=1*dt for DuFort-Frankel Method
void Methods::DFF1(double D, double dx, double dt, double** A0, int rows, int cols)
{
	double** A1 = A0;
	double belta1 = D * dt / dx / dx;

	for (int i = 0; i < cols; i++)
	{		
		for (int j = 2; j < rows; j++)
		{
			A1[i][j+1] = (1 - 2 * belta1) / (1 + 2 * belta1) * A1[i][j-1] + 2 * belta1 / (1 + 2 * belta1) * A1[i + 1][j]+2*belta1/(1+2*belta1)*A1[i-1][j];
		}
	}
}

//This Method is using Laasonen Simple Implicit Method to calculate t=1*dt for DuFort-Frankel Method
void Methods::DFF2(double D, double dx, double dt, double** A0, int rows, int cols)
{
	double** A1 = A0;
	double belta1 = D * dt / dx / dx;

	for (int i = 0; i < cols; i++)
	{
		for (int j = 2; j < rows; j++)
		{
			A1[i][j + 1] = (1 - 2 * belta1) / (1 + 2 * belta1) * A1[i][j - 1] + 2 * belta1 / (1 + 2 * belta1) * A1[i + 1][j] + 2 * belta1 / (1 + 2 * belta1) * A1[i - 1][j];
		}
	}
}

//This Method is using Richardson Method to calculate t=1*dt for DuFort-Frankel Method
void Methods::DFF3(double D,  double dx, double dt, double** A0, int rows, int cols)
{
	double** A1 = A0;
	double belta1 = D * dt / dx / dx;

	for (int i = 0; i < cols; i++)
	{
		for (int j = 2; j < rows; j++)
		{
			A1[i][j + 1] = (1 - 2 * belta1) / (1 + 2 * belta1) * A1[i][j - 1] + 2 * belta1 / (1 + 2 * belta1) * A1[i + 1][j] + 2 * belta1 / (1 + 2 * belta1) * A1[i - 1][j];
		}
	}
}

//This Method is using FTCS to calculate t=1*dt for DuFort-Frankel Method
void Methods::DFF4(double D,double dx, double dt, double** A0, int rows, int cols)
{
	double** A1 = A0;
	double belta1 = D * dt / dx / dx;

	for (int i = 0; i < cols; i++)
	{
		for (int j = 2; j < rows; j++)
		{
			A1[i][j + 1] = (1 - 2 * belta1) / (1 + 2 * belta1) * A1[i][j - 1] + 2 * belta1 / (1 + 2 * belta1) * A1[i + 1][j] + 2 * belta1 / (1 + 2 * belta1) * A1[i - 1][j];
		}
	}
}

//This Method is using Analytical solution to calculate t=1*dt for DuFort-Frankel Method
void Methods::DFF5(double D, double dx, double dt, double** A0, int rows, int cols,int m1)
{
	double** A1 = A0;
	double belta1 = D * dt / dx / dx;
	double S1 = 0.0;
	double Tp = 0.0;
	for (int i = 0; i < rows; i++)
	{
		double ic = i;
		for (int m = 1; m < m1; m++)
		{
			double E = exp(-93 * pow(m * 0.10129, 2) * 1 * dt);
			double sincc = sin(m * 0.10129 * ic * dx);
			S1 = S1 + E * (1 - (-1) ^ m) / (m * 3.14159) * sincc;
			//cout << m << endl;
			//cout << E << endl;
			//cout << sincc << endl;
			Tp = 149 - (222 * S1);

			A1[i][1] = Tp;
			S1 = 0.0;
		}
	}
	for (int i = 0; i < cols; i++)
	{
		for (int j = 2; j < rows; j++)
		{
			A1[i][j + 1] = (1 - 2 * belta1) / (1 + 2 * belta1) * A1[i][j - 1] + 2 * belta1 / (1 + 2 * belta1) * A1[i + 1][j] + 2 * belta1 / (1 + 2 * belta1) * A1[i - 1][j];
		}
	}
}

void Methods::RCCE(double D, double** AT, double dx, double dt, double T_sur, double T_in, double** A0, int rows, int cols)
{

	double** A2 = A0;

	double belta1 = D * dt / dx / dx / 2;

	for (int i = 0; i < rows; i++)
	{
		A2[i] = new double[cols];
	};
	AT = A2;

}

void Methods::LSIFCI(double D, double** AMO, double dx, double dt, double T_sur, double T_in, double** A0, int rows, int cols)
{

	double** A3 = A0;
	double belta3 = D * dt / dx / dx;

	//this is the matrix A of Ax=b;
	for (int i = 0; i < rows; i++)
	{
		AMO[i][i - 1] = -belta3;
		AMO[i][i] = (1 + 2 * belta3);
		AMO[i][i + 1] = (-belta3);
	}
	AMO = A3;
}

void Methods::CNT(double D, double** AMO,  double dx, double dt, double** A0, int rows, int cols)
{
	double** A4 = A0;
	double A4C;
	double belta4 = D * dt / 2 / dx / dx;

	AMO[0][0] = (1 + 2 * belta4);
	AMO[cols][cols] = (1 + 2 * belta4);
	
	//This is a part for creating tridiagonal matrix algorithm
	for (int i = 1; i < cols; i++)
	{
		AMO[i - 1][i] = -1 * belta4;
		AMO[i][i] = (1 + 2 * belta4);
		AMO[i + 1][i] = -1 * belta4;
	}

	//This is a part calculating and putting into the grid
	for (int j = 1; j < rows; j++)
	{
		A4[cols][j] = A4[cols][j-1] / AMO[cols][cols];

		for (int i = 1; i < cols; i++)
		{
			A4C=belta4*A4[i][j-1]+(1-2*belta4)*A4[i+1][j-1]+belta4*A4[i-1][j-1];
			A4[i][j]=(A4C-A4[i-1][j])/AMO[i+1][i];
		}
	}
}
