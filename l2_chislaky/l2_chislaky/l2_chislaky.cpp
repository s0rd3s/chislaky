#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

const double eps = 1e-9;
const double M = 0.001;
const int NIT = 110;

double** initial(int n, int m)
{
	double** A = new double* [n];
	for (int i = 0; i < n; i++)
		A[i] = new double[m];
	return A;
}

double function1(double x1, double x2)
{
	return (x1 * x1 * x1 + x2 * x2 * x2 - 6 * x1 + 3);
}

double function2(double x1, double x2)
{
	return (x1 * x1 * x1 - x2 * x2 * x2 - 6 * x2 + 2);
}


double func11(double x1, double x2)
{
	return ((function1(x1 + M * x1, x2) - function1(x1, x2)) / (M * x1));
}


double func12(double x1, double x2)
{
	return ((function1(x1, x2 + M * x2) - function1(x1, x2)) / (M * x2));
}



double func21(double x1, double x2)
{
	return ((function2(x1 + M * x1, x2) - function2(x1, x2)) / (M * x1));
}


double func22(double x1, double x2)
{
	return ((function2(x1, x2 + M * x2) - function2(x1, x2)) / (M * x2));
}


void J(double** a, double x1, double x2)
{
	a[0][0] = func11(x1, x2);
	a[0][1] = func12(x1, x2);
	a[1][0] = func21(x1, x2);
	a[1][1] = func22(x1, x2);
}

void minus_vector_nevjazki(double* F, double x1, double x2)
{
	F[0] = -function1(x1, x2);
	F[1] = -function2(x1, x2);
}

double* gauss(double** matrix, int n, int m)
{
	//prjamoj
	double elem;
	for (int j = 0; j < n; j++)
	{
		double max = 0;
		int coord_str = 0;
		for (int t = j; t < n; t++)
		{
			if (abs(matrix[t][j]) > max)
			{
				max = abs(matrix[t][j]); coord_str = t;
			}
		}
		if (max > abs(matrix[j][j]))
		{
			double* ptr = matrix[j];
			matrix[j] = matrix[coord_str];
			matrix[coord_str] = ptr;
		}
		elem = matrix[j][j];
		for (int c = j; c < m; c++)
		{
			matrix[j][c] /= elem;   //delenije stroki na elem
		}

		for (int i2 = j + 1; i2 < n; i2++)
		{
			elem = matrix[i2][j];
			for (int k = j; k < m; k++)
				matrix[i2][k] -= elem * matrix[j][k];
		}

	}
	//obratnyj
	double* xx = new double[m];
	xx[n - 1] = matrix[n - 1][n];
	for (int i = n - 2; i >= 0; i--)
	{
		xx[i] = matrix[i][n];
		for (int j = i + 1; j < n; j++) xx[i] -= matrix[i][j] * xx[j];
	}

	cout << endl;

	return xx;
}

double* neutone(int n, double x1, double x2)
{
	double** Jako;
	Jako = initial(n, n);

	double** newmatrix;
	newmatrix = initial(n, n + 1);

	double* F = new double[n];
	double* delta = new double[n];
	double* resh = new double[n];
	resh[0] = x1;
	resh[1] = x2;
	double delta1, delta2;
	int k = 1;

	cout << setw(10) << "x1" << setw(10) << "x2" << setw(17) << "delta1" << setw(17) << "delta2" << setw(6) << "k";
	do
	{
		minus_vector_nevjazki(F, x1, x2);
		J(Jako, x1, x2);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				newmatrix[i][j] = Jako[i][j];
			}
			newmatrix[i][n] = F[i];
		}

		delta = gauss(newmatrix, n, n + 1);
		for (int i = 0; i < n; i++)
			resh[i] += delta[i];

		double max1 = 0;
		double max2 = 0;
		for (int i = 0; i < n; i++)
		{
			if (abs(F[i]) > max1)
				max1 = abs(F[i]);

			if (abs(resh[i]) < 1)
			{
				if (abs(delta[i]) > max2)
					max2 = abs(delta[i]);
			}
			if (abs(delta[i] >= 1))
			{
				if (abs(delta[i] / resh[i]) > max2)
					max2 = abs(delta[i]);
			}
		}

		delta1 = max1;
		delta2 = max2;
		cout << endl;


		x1 = resh[0];
		x2 = resh[1];
		for (int i = 0; i < n; i++)
			cout << setw(10) << resh[i] << "   ";
		cout << setw(13) << delta1 << "   " << setw(13) << delta2 << "   " << setw(2) << k;
		cout << endl;

		k++;
		if (k >= NIT)
		{
			cout << "\n   IER = 2 \n";
			return NULL;
			break;
		}

	} while (delta1 > eps || delta2 > eps);

	return resh;
}


int main()
{
	double x1, x2;
	x1 = 0.5; x2 = 0.2;
	int n = 2;
	double* otvet = new double[n];
	otvet = neutone(n, x1, x2);

	if (otvet != NULL)
	{
		cout << "____________________________________________________________" << endl;
		for (int i = 0; i < n; i++)
			cout << otvet[i] << endl;
	}

	return 0;
}