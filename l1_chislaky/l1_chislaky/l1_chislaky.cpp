#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
using namespace std;

double** initial(int n, int m)
{
	double** A = new double* [n];
	for (int i = 0; i < n; i++)
		A[i] = new double[m];
	return A;
}
void out(double** A, int str, int stolb)
{
	cout << endl;
	for (int i = 0; i < str; i++)
	{
		for (int j = 0; j < stolb; j++)
			cout << setw(10) << A[i][j];
		cout << endl;
	}
	cout << endl;
}

double* gauss(double** matrix, int n, int m)
{														//pryamoi metod
	double elem;
	for (int j = 0; j < n; j++)
	{
		int max = 0;
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
		{												//delenie stroki na element
			matrix[j][c] /= elem;  
		}

		for (int i2 = j + 1; i2 < n; i2++)
		{
			elem = matrix[i2][j];
			for (int k = j; k < m; k++)
				matrix[i2][k] -= elem * matrix[j][k];
		}

	}
														//obratniy
	double* xx = new double[m];
	xx[n - 1] = matrix[n - 1][n];
	for (int i = n - 2; i >= 0; i--)
	{
		xx[i] = matrix[i][n];
		for (int j = i + 1; j < n; j++) xx[i] -= matrix[i][j] * xx[j];
	}

	cout << "result:" << endl;
	for (int i = 0; i < n; i++)
		cout << xx[i] << " ";
	cout << endl;
	delete[] matrix;

	return xx;
}

void mult(double** x, int n, int m, double* y, double* mt)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			mt[i] += x[i][j] * y[j];
}

int main()
{
	int n = 4, m = 4;
	m++;
	double** matrix;
	matrix = initial(n, m);

	matrix[0][0] = 7.90;
	matrix[0][1] = 5.60;
	matrix[0][2] = 5.70;
	matrix[0][3] = -7.20;
	matrix[0][4] = 6.68;

	matrix[1][0] = 8.50;
	matrix[1][1] = -4.80;
	matrix[1][2] = 0.80;
	matrix[1][3] = 3.50;
	matrix[1][4] = 9.95;

	matrix[2][0] = 4.30;
	matrix[2][1] = 4.20;
	matrix[2][2] = -3.20;
	matrix[2][3] = 9.30;
	matrix[2][4] = 8.60;

	matrix[3][0] = 3.20;
	matrix[3][1] = -1.40;
	matrix[3][2] = -8.90;
	matrix[3][3] = 3.30;
	matrix[3][4] = 1.00;

	double** copy = initial(n, n);
	double* svoboda = new double[n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			copy[i][j] = matrix[i][j];
			svoboda[i] = matrix[i][n];
		}

	cout << "matrix: " << endl;
	out(matrix, n, m);
	double* resh = new double[n];
	resh = gauss(matrix, n, m);

	double* nevjaz = new double[n];
	for (int i = 0; i < n; i++)
		nevjaz[i] = 0;

	mult(copy, n, n, resh, nevjaz);

	for (int i = 0; i < n; i++)
		nevjaz[i] -= svoboda[i];

	cout << "nevjazka:" << endl;
	for (int i = 0; i < n; i++)
		cout << nevjaz[i] << " ";
	cout << endl;

	return 0;
}