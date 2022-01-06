#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#define PRES double
#define NXB 15
#define NYB 12
#define NX NXB*3+1
#define NY NYB*3+1
#define NYK2 NYB*2
#define REP 3000
#define DEL 100
#define AMAT 1.1f
#define TEM1 5.0f
#define TEM2 15.0f
#define HX 0.2f
#define HY 0.3f
using namespace std;
int main(int argc, char** argv)
{
	int i, j, k;
	int idt = 0;
	int ndt = 0;
	PRES T1 = TEM1, T2 = TEM2, h = HX, r = HY, a = AMAT, t0;
	PRES T[NY][NX], TT[NY][NX];
	PRES rr = min(h, r);
	PRES tau = 0.25f * rr * rr / a;
	PRES alf_1 = -h / r;
	PRES alf_2 = -r / h;
	PRES alf_3 = 0.5f * alf_2;
	PRES alf_4 = 0.5f * alf_1;
	PRES bet_1 = a * tau / (h * r);
	PRES bet_2 = 2.0f * bet_1;
	PRES bet_3 = 4.0f * bet_1;
	PRES gam_1 = -2.f * (alf_1 + alf_2);
	PRES gam_3 = -(alf_1 + alf_2);
	PRES gam_4 = -(alf_3 + alf_4);
	char filename[128];
	for (j = 0; j < NY; j++)
	{
		for (i = 0; i < NX; i++)
		{
			T[j][i] = 0.0f;
		}
	}
	for (j = 0; j <= NYB; j++)
	{
		T[j][0] = T1;
		TT[j][0] = T1;
	}
	for (j = NYK2; j < NY; j++)
	{
		T[j][NX - 1] = T2;
		TT[j][NX - 1] = T2;
	}
	ofstream fout("T1.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	for (j = 0; j < NY; j++)
	{
		for (i = 0; i < NX; i++)
		{
			PRES w = T[j][i];
			fout.write((char*)&w, sizeof w);
		}
	}
	fout.close();
	for (k = 0; k < REP; k++)
	{
		for (j = 0; j < NY; j++)
		{
			for (i = 0; i < NX; i++)
			{
				t0 = T[j][i];
				if (i > 0 && i < NX - 1 && j == 0)
				{
					TT[j][i] = t0 - bet_2 * (alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j + 1][i] + gam_3 * t0);
				}
				else if (i > 0 && i < NX - 1 && j == NY - 1)
				{
					TT[j][i] = t0 - bet_2 * (alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i] + gam_3 * t0);
				}
				else if (i == 0 && j > NYB && j < NY - 1)
				{
					TT[j][i] = t0 - bet_2 * (alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i + 1] + gam_3 * t0);
				}
				else if (i == NX - 1 && j > 0 && j < NYK2)
				{
					TT[j][i] = t0 - bet_2 * (alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i - 1] + gam_3 * t0);
				}
				else if (i == 0 && j == NY - 1)
				{
					TT[j][i] = t0 - bet_3 * (alf_3 * T[j][i + 1] + alf_4 * T[j - 1][i] + gam_4 * t0);
				}
				else if (i == NX - 1 && j == 0)
				{
					TT[j][i] = t0 - bet_3 * (alf_3 * T[j][i - 1] + alf_4 * T[j + 1][i] + gam_4 * t0);
				}
				else if (i > 0 && i < NX - 1 && j>0 && j < NY - 1)
				{
					TT[j][i] = t0 - bet_1 * (alf_1 * (T[j][i - 1] + T[j][i + 1]) + alf_2 * (T[j - 1][i] + T[j + 1][i]) + gam_1 * t0);
				}
			}
		}
		for (j = 0; j < NY; j++)
		{
			for (i = 0; i < NX; i++)
			{
				T[j][i] = TT[j][i];
			}
		}
		idt++;
		if (idt == DEL)
		{
			idt = 0; ndt++;
			sprintf_s(filename, sizeof(filename), "T%d.dat", ndt + 1);
			ofstream fout(filename, ios_base::out | ios_base::trunc | ios_base::binary);
			for (j = 0; j < NY; j++)
			{
				for (i = 0; i < NX; i++)
				{
					PRES w = T[j][i];
					fout.write((char*)&w, sizeof w);
				}
			}
			fout.close();
		}
	}
	int n_x = NX; int n_y = NY; int n_k = ndt;
	ofstream fou("Param.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fou.write((char*)&n_x, sizeof n_x);
	fou.write((char*)&n_y, sizeof n_y);
	fou.write((char*)&n_k, sizeof n_y);
	fou.close();
}