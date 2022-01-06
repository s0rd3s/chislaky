#include "stdafx.h"
#include <stdlib.h> 
#include <math.h> 
#include <fstream>
#include <iostream>
#define PRES  double 
#define NXB 15   //число шагов сетки по оси X в одном сегменте,
#define NX NXB*3+1  //число узлов сетки по оси X во всей структуре
#define NYB 12   //аналогично по оси Y
#define NY NYB*3+1
#define REP 3000  //максимальное число ПВР-итераций
#define EPSL 1.e-5  //точность
#define LL 1.7f    //итерационный параметр
#define TEM1 5.0f   //температура контакта К1
#define TEM2 15.0f  //температура контакта К2
#define HX 0.2f     //шаги сетки
#define HY 0.3f
using namespace std;

void maxpvr(PRES* t1, PRES* del, PRES* maxdel)  //расчёт максимального относительного изменения температуры на текущей итерации ПВР-методa
{
	PRES d = fabs(*del) / fabs(*t1);
	if (d > *maxdel) *maxdel = d;
}

int main(int argc, char** argv)
{
	ofstream foutT("dT1.dat", ios_base::out | ios_base::trunc | ios_base::binary); //текущая максимальная относительная поправка температуры
	int    i1, i2, i3, j1, j2, j3, rp, i, j, k = 0;
	PRES  T1 = TEM1, T2 = TEM2, h = HX, r = HY, tx, t0, t1, del, maxdel = 0.0f;
	PRES T[NY][NX];  //содержит искомое распределение
	PRES lam = LL;  //итерационный параметр ljambda
	PRES eps = EPSL;   //точность решения СЛАУ
	int prz = 1;   //условие окончания ПВР-итераций
	int nT = 0;   //число итераций ПВР
	PRES  alf_1 = -h / r;
	PRES  alf_2 = -r / h;
	PRES  alf_3 = alf_2 * 0.5f;
	PRES  alf_4 = alf_1 * 0.5f;
	PRES  gam_1 = -2.f * (alf_1 + alf_2);
	PRES  gam_2 = -1.5f * (alf_1 + alf_2);
	PRES  gam_3 = -(alf_1 + alf_2);
	PRES  gam_4 = -(alf_3 + alf_4);
	i1 = NXB;
	i2 = i1 + NXB;
	i3 = i2 + NXB;
	j1 = NYB;
	j2 = j1 + NYB;
	j3 = j2 + NYB;
	rp = REP;
	for (j = 0; j <= j3; j++)   //обнуляется массив температур
		for (i = 0; i <= i3; i++)
			T[j][i] = 0.0f;

	for (j = j1; j <= j3; j++) T[j][0] = T1;    //затем в него заносятся температуры контактов
	for (i = 0; i <= i1; i++) T[j3][i] = T1;
	for (j = 0; j <= j2; j++) T[j][i3] = T2;
	for (i = i2; i <= i3; i++) T[0][i] = T2;

	while (k < rp && prz == 1)  //каждый проход цикла выполняется, если число итераций меньше максимально допустимого числа и не достигнута требуемая точность
	{
		k++;
		for (j = 0; j <= j3; j++) //начинается обход узлов двумерной разностной сетки
		{
			for (i = 0; i <= i3; i++)
			{
				t0 = T[j][i];
				if (i == 0 && j > 0 && j < j1)  //AB
				{
					tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i + 1]) / gam_3;
					del = lam * (tx - t0); //поправка температуры
					t1 = t0 + del; //новое значение температуры
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel); //расчёт максимальной относительной поправки температуры
				}
				else if (i == 0 && j == 0)  //A
				{
					tx = -(alf_3 * T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_4;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i > 0 && i < i2 && j == 0)  //AF
				{
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j + 1][i]) / gam_3;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}

				if (i == i1 && j > j2 && j < j3) //CD
				{
					tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i - 1]) / gam_3;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i == i1 && j == j2) //D
				{
					tx = -(alf_2 * T[j][i - 1] + alf_4 * T[j + 1][i] + alf_3 * T[j][i + 1] + alf_1 * T[j - 1][i]) / gam_2;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}

				if (i > i1 && i < i3 && j == j2)  //DE
				{
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i]) / gam_3;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i > 0 && i < i3 && j > 0 && j < j2) //AGEP
				{
					tx = -(alf_1 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * (T[j][i - 1] + T[j][i + 1])) / gam_1;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i > 0 && i < i1 && j>j2 - 1 && j < j3)  //GECD
				{
					tx = -(alf_1 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * (T[j][i - 1] + T[j][i + 1])) / gam_1;
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
			}
		}
		nT++;
		PRES w = maxdel;
		foutT.write((char*)&w, sizeof w);  //дописывается текущая максимальная относительная поправка температуры
		if (maxdel < eps) prz = 0; maxdel = 0.0f; //проверяется условие выхода из цикла
	}
	foutT.close();  //eсли цикл ПВР закончен, закрывается файл “dT.dat”

	ofstream fouT("nT1.dat", ios_base::out | ios_base::trunc | ios_base::binary); //в файл “nT.dat” записывается число итераций ПВР
	fouT.write((char*)&nT, sizeof nT);
	fouT.close();

	ofstream fout("Pole1.dat", ios_base::out | ios_base::trunc | ios_base::binary);  //в файл “Pole.dat” записывается распределение температуры
	for (j = 0; j < NY; j++)
	{
		for (i = 0; i < NX; i++)
		{
			PRES w = T[j][i];
			fout.write((char*)&w, sizeof w);
		}
	}
	fout.close();
	int n_x = NX;
	int n_y = NY;

	ofstream fou("Param1.dat", ios_base::out | ios_base::trunc | ios_base::binary); //в файл “Param.dat” записывается размерность разностной сетки
	fou.write((char*)&n_x, sizeof n_x);
	fou.write((char*)&n_y, sizeof n_y); fou.close();
	cout << k << endl;
}