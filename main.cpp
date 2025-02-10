/*
bnuGRMHD ©️ 2025
Date: 2024/02/02
*/
#include <cmath>
#include <ctime>
#include <fstream>
#include <thread>
#include "Metric.h"
#include "utils.h"
#include "init.h"
#include "metric/mks.h"

using std::thread;
using std::ref;

int main(int argc, char* argv[])
{
	std::ofstream ofs;
	ofs.open("grmhd.log", std::ios::out);
	auto totalTime = 0.;
	auto totalPhysicalTime = 0.;
	// 初始化：1) 度规张量设置；2) 度规场初始化；3) 守恒量赋初值
	{
		init_metric();
		for (int i = 0; i < N1 + 2 * NG; i++)
			for (int j = 0; j < N2 + 2 * NG; j++)
				for (int k = 0; k < N3 + 2 * NG; k++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
							metricFuncField(i, j, k).m(row, col) = metricFunc(row, col)(X1min + i * dx1Ghost, X2min + j * dx2Ghost, X3min + k * dx3Ghost);

		// 利用中心差分计算alpha的导数
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					alphaDiffField[i][j][k][0] = 0;
					alphaDiffField[i][j][k][1] = (metricFuncField(i + NG + 1, j + NG, k + NG).alpha() - metricFuncField(i + NG - 1, j + NG, k + NG).alpha()) / (2 * dx1);
					alphaDiffField[i][j][k][2] = (metricFuncField(i + NG, j + NG + 1, k + NG).alpha() - metricFuncField(i + NG, j + NG - 1, k + NG).alpha()) / (2 * dx2);
					alphaDiffField[i][j][k][3] = (metricFuncField(i + NG, j + NG, k + NG + 1).alpha() - metricFuncField(i + NG, j + NG, k + NG - 1).alpha()) / (2 * dx3);
				}

		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					for (int l = 0; l < 4; l++)
						for (int row = 0; row < 4; row++)
							for (int col = 0; col < 4; col++)
								metricDiffField(i, j, k, l).m(row, col) = metricDiff(row, col, l)(X1min + i * dx1, X2min + j * dx2, X3min + k * dx3);

		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
						{
							metricFuncHalfField1(i, j, k).m(row, col) = metricFunc(row, col)(X1min + (2 * i + 3) * dx1 / 2, X2min + (j + 2) * dx2, X3min + (k + 2) * dx3);
							metricFuncHalfField2(i, j, k).m(row, col) = metricFunc(row, col)(X1min + (i + 2) * dx1, X2min + (2 * i + 3) * dx2 / 2, X3min + (k + 2) * dx3);
							metricFuncHalfField3(i, j, k).m(row, col) = metricFunc(row, col)(X1min + (i + 2) * dx1, X2min + (j + 2) * dx2, X3min + (2 * k + 3) * dx3 / 2);
						}
		init();
		char filename[13];
		sprintf(filename, "./data/data%0.4d.bin", 0);
		write_bin(fopen(filename, "wb"));
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					ksi[i][j][k] = (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]) * (1 + pow(prim[i][j][k][U1], 2) + pow(prim[i][j][k][U2], 2) + pow(prim[i][j][k][U3], 2));

	}
	/*
	2.迭代方程
		1) 鬼化：给主要量P加上鬼格边界条件
		2) 插值：给主要量P插值
		3) 计算流(flux)
		4) 半步长迭代
		5) 受约束运输
		6) 整步长迭代
	*/
	for(int epoch = 1; epoch <= epochNum; epoch++)
	{
		auto start = clock();
		// 时间步长
		double Delta_t = 1;
		ghostify();

		interpolate();
		
		{
			thread th1(basicCalc, primL1, ref(conL1), ref(fluxL1), ref(srcL1), ref(cpL1), ref(cnL1), ref(metricFuncHalfField1), 0);
			thread th2(basicCalc, primL2, ref(conL2), ref(fluxL2), ref(srcL2), ref(cpL2), ref(cnL2), ref(metricFuncHalfField2), 1);
			thread th3(basicCalc, primL3, ref(conL3), ref(fluxL3), ref(srcL3), ref(cpL3), ref(cnL3), ref(metricFuncHalfField3), 2);
			thread th4(basicCalc, primR1, ref(conR1), ref(fluxR1), ref(srcR1), ref(cpR1), ref(cnR1), ref(metricFuncHalfField1), 0);
			thread th5(basicCalc, primR2, ref(conR2), ref(fluxR2), ref(srcR2), ref(cpR2), ref(cnR2), ref(metricFuncHalfField2), 1);
			thread th6(basicCalc, primR3, ref(conR3), ref(fluxR3), ref(srcR3), ref(cpR3), ref(cnR3), ref(metricFuncHalfField3), 2);
			th1.join();
			th2.join();
			th3.join();
			th4.join();
			th5.join();
			th6.join();
		}

		{
			thread th1(fluxCalc, cpL1, cpR1, cnL1, cnR1, conL1, conR1, fluxL1, fluxR1, ref(fluxHLL1), ref(fluxTVDLF1));
			thread th2(fluxCalc, cpL2, cpR2, cnL2, cnR2, conL2, conR2, fluxL2, fluxR2, ref(fluxHLL2), ref(fluxTVDLF2));
			thread th3(fluxCalc, cpL3, cpR3, cnL3, cnR3, conL3, conR3, fluxL3, fluxR3, ref(fluxHLL3), ref(fluxTVDLF3));
			th1.join();
			th2.join();
			th3.join();
		}
		
#pragma omp parallel for
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					for(int l = 0; l < 8; l++)
					{
						fluxLLF1[i][j][k][l] = theta * fluxHLL1[i][j][k][l] + (1 - theta) * fluxTVDLF1[i][j][k][l];
						fluxLLF2[i][j][k][l] = theta * fluxHLL2[i][j][k][l] + (1 - theta) * fluxTVDLF2[i][j][k][l];
						fluxLLF3[i][j][k][l] = theta * fluxHLL3[i][j][k][l] + (1 - theta) * fluxTVDLF3[i][j][k][l];
					}

		// 4.半步长迭代
		prim2con(prim, con);

		prim2src(prim, con, src);

#pragma omp parallel for
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					auto c1max = max(0, cpR1[i][j][k], cpL1[i][j][k]);
					auto c1min = -min(-0, cnR1[i][j][k], cnL1[i][j][k]);
					auto c2max = max(0, cpR2[i][j][k], cpL2[i][j][k]);
					auto c2min = -min(-0, cnR2[i][j][k], cnL2[i][j][k]);
					auto c3max = max(0, cpR3[i][j][k], cpL3[i][j][k]);
					auto c3min = -min(-0, cnR3[i][j][k], cnL3[i][j][k]);
					auto c1 = abs(max(c1max, c1min));
					auto c2 = abs(max(c2max, c2min));
					auto c3 = abs(max(c3max, c3min));
					Delta_t = min(cour * min(dx1 / (2 * c1), dx2 / (2 * c2), dx3 / (2 * c3)), Delta_t);
				}

#pragma omp parallel for
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
					for (int l = 0; l < 8; l++)
						conHalf[i][j][k][l] = con[i][j][k][l] + src[i][j][k][l] * Delta_t / 2
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField1(i + 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1[i + 1][j][k][l] - sqrt(metricFuncHalfField1(i - 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1[i][j][k][l])
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField2(i, j + 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2[i][j + 1][k][l] - sqrt(metricFuncHalfField1(i, j - 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2[i][j][k][l])
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField3(i, j, k + 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3[i][j][k + 1][l] - sqrt(metricFuncHalfField1(i, j, k - 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3[i][j][k][l]);


		con2prim(conHalf);

		ghostify();

		interpolate();

		{
			thread th1(basicCalc, primL1, ref(conL1), ref(fluxL1), ref(srcL1), ref(cpL1), ref(cnL1), ref(metricFuncHalfField1), 0);
			thread th2(basicCalc, primL2, ref(conL2), ref(fluxL2), ref(srcL2), ref(cpL2), ref(cnL2), ref(metricFuncHalfField2), 1);
			thread th3(basicCalc, primL3, ref(conL3), ref(fluxL3), ref(srcL3), ref(cpL3), ref(cnL3), ref(metricFuncHalfField3), 2);
			thread th4(basicCalc, primR1, ref(conR1), ref(fluxR1), ref(srcR1), ref(cpR1), ref(cnR1), ref(metricFuncHalfField1), 0);
			thread th5(basicCalc, primR2, ref(conR2), ref(fluxR2), ref(srcR2), ref(cpR2), ref(cnR2), ref(metricFuncHalfField2), 1);
			thread th6(basicCalc, primR3, ref(conR3), ref(fluxR3), ref(srcR3), ref(cpR3), ref(cnR3), ref(metricFuncHalfField3), 2);
			th1.join();
			th2.join();
			th3.join();
			th4.join();
			th5.join();
			th6.join();
		}

		{
			thread th1(fluxCalc, cpL1, cpR1, cnL1, cnR1, conL1, conR1, fluxL1, fluxR1, ref(fluxHLL1), ref(fluxTVDLF1));
			thread th2(fluxCalc, cpL2, cpR2, cnL2, cnR2, conL2, conR2, fluxL2, fluxR2, ref(fluxHLL2), ref(fluxTVDLF2));
			thread th3(fluxCalc, cpL3, cpR3, cnL3, cnR3, conL3, conR3, fluxL3, fluxR3, ref(fluxHLL3), ref(fluxTVDLF3));
			th1.join();
			th2.join();
			th3.join();
		}
		/*
		5) 平滑化
		*/

#pragma omp parallel for
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
				{
					fluxSmoothLLF1[i][j][k][6] = 0.125 * (2 * fluxLLF1[i][j][k][6] + fluxLLF1[i][j + 1][k][6] + fluxLLF1[i][j - 1][k][6] - fluxLLF2[i][j][k][5] - fluxLLF2[i][j + 1][k][5] - fluxLLF2[i - 1][j][k][5] - fluxLLF2[i - 1][j + 1][k][5]);
					fluxSmoothLLF2[i][j][k][5] = 0.125 * (2 * fluxLLF2[i][j][k][5] + fluxLLF2[i + 1][j][k][5] + fluxLLF2[i - 1][j][k][5] - fluxLLF1[i][j][k][6] - fluxLLF1[i + 1][j][k][6] - fluxLLF1[i][j - 1][k][6] - fluxLLF1[i + 1][j - 1][k][6]);
					fluxSmoothLLF1[i][j][k][7] = 0.125 * (2 * fluxLLF1[i][j][k][7] + fluxLLF1[i][j][k + 1][7] + fluxLLF1[i][j][k - 1][7] - fluxLLF3[i][j][k][5] - fluxLLF3[i][j][k + 1][7] - fluxLLF3[i - 1][j][k][5] - fluxLLF3[i - 1][j][k + 1][7]);
					fluxSmoothLLF3[i][j][k][5] = 0.125 * (2 * fluxLLF3[i][j][k][5] + fluxLLF3[i + 1][j][k][5] + fluxLLF3[i - 1][j][k][5] - fluxLLF1[i][j][k][7] - fluxLLF1[i + 1][j][k][7] - fluxLLF1[i][j][k - 1][7] - fluxLLF1[i + 1][j][k - 1][7]);
					fluxSmoothLLF2[i][j][k][7] = 0.125 * (2 * fluxLLF2[i][j][k][7] + fluxLLF2[i][j][k + 1][7] + fluxLLF2[i][j][k - 1][7] - fluxLLF3[i][j][k][6] - fluxLLF3[i][j][k + 1][7] - fluxLLF3[i][j - 1][k][6] - fluxLLF3[i][j - 1][k + 1][7]);
					fluxSmoothLLF3[i][j][k][6] = 0.125 * (2 * fluxLLF3[i][j][k][6] + fluxLLF3[i][j + 1][k][6] + fluxLLF3[i][j - 1][k][6] - fluxLLF2[i][j][k][7] - fluxLLF2[i][j + 1][k][7] - fluxLLF2[i][j][k - 1][7] - fluxLLF2[i][j + 1][k - 1][7]);
				}
		/*
		6) 整步迭代
		*/

#pragma omp parallel for
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
					for (int l = 0; l < 8; l++)
						con[i][j][k][l] = con[i][j][k][l] + src[i][j][k][l] * Delta_t / 2
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField1(i + 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1[i + 1][j][k][l] - sqrt(metricFuncHalfField1(i - 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1[i][j][k][l])
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField2(i, j + 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2[i][j + 1][k][l] - sqrt(metricFuncHalfField1(i, j - 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2[i][j][k][l])
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField3(i, j, k + 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3[i][j][k + 1][l] - sqrt(metricFuncHalfField1(i, j, k - 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3[i][j][k][l]);

		con2prim(con);

		fix();
		
		totalTime += clock() - start;
		totalPhysicalTime += Delta_t;
		std::cout << "Time(ms): " << clock() - start << "\tPhysical Time: " << Delta_t << "\tTotal Physical Time: " << totalPhysicalTime << std::endl;
		char filename[13];
		sprintf(filename, "./data/data%0.4d.bin", epoch);
		write_bin(fopen(filename, "wb"));
		ofs << "--------Epoch--------" << epoch << std::endl;
		for(int i = 0; i < N1; i++)
			for(int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					for (int l = 0; l < NPRIM; l++)
						ofs << prim[i][j][k][l] << "\t";
					ofs << std::endl;
				}
		ofs << "Time(ms): " << clock() - start << "\tPhysical Time: " << Delta_t << "\tTotal Physical Time: " << totalPhysicalTime << std::endl;
	}
	ofs << "Total times(ms): " << totalTime << std::endl << "Average time(ms): " << totalTime / epochNum << std::endl;
	return 0;
}
