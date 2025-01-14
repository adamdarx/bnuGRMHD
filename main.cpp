﻿/*
bnuGRMHD ©️ 2024
Date: 2024/12/21
本文件是初代GRMHD的示例代码
TODO:
	1) 目前来说10x10x1的网格能达到0.4s迭代一次
	2) 初始化问题没有解决
*/

#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#include <ctime>
#include "Metric.h"
#include "utils.h"
#include "init.h"


using std::thread;
using std::ref;
using std::vector;

int main()
{
	auto totalTime = 0;
	/*
	1.初始化
		1) 度规张量设置
		2) 度规场初始化
		3) 守恒量赋初值
	*/
	
	/*
	1) 度规张量设置
		创建Eigen::Tensor<MetricComponent, 2>对象，每一个分量本质都是一个函数（类型为MetricComponent，数学地写就是R^3->R的映射定义见Metric.h）
		默认全为零分量（注意不是double 0而是预先定义的ZERO_COMPONENT，定义见Metric.h）
	*/
	metricFunc.setConstant(ZERO_COMPONENT);
	metricDiff.setConstant(ZERO_COMPONENT);
	metricFunc(0, 0) = [](double x1, double x2, double x3) {return -1 + 2 * M / x1; };
	metricFunc(1, 1) = [](double x1, double x2, double x3) {return 1 - 2 * M / x1; };
	metricFunc(2, 2) = [](double x1, double x2, double x3) {return x1 * x1; };
	metricFunc(3, 3) = [](double x1, double x2, double x3) {return pow(x1 * sin(x2), 2); };
	metricDiff(0, 0, 0) = [](double x1, double x2, double x3) {return -M / pow(x1, 2); };
	metricDiff(1, 1, 0) = [](double x1, double x2, double x3) {return M / pow(x1, 2); };


	for (int i = 0; i < N1 + 2 * NG; i++)
		for (int j = 0; j < N2 + 2 * NG; j++)
			for (int k = 0; k < N3 + 2 * NG; k++)
				for(int row = 0; row < 4; row++)
					for (int col = 0; col < 4; col++)
						metricFuncField(i, j, k).m(row, col) = metricFunc(row, col)(X1min + i * dx1, X2min + i * dx2, X3min + i * dx3);
	
	for (int i = 0; i < N1 + 2 * NG; i++)
		for (int j = 0; j < N2 + 2 * NG; j++)
			for (int k = 0; k < N3 + 2 * NG; k++)
				for (int l = 0; l < 4; l++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
							metricDiffField(i, j, k, l).m(row, col) = metricDiff(row, col, l)(X1min + i * dx1, X2min + i * dx2, X3min + i * dx3);
	
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
	prim.setZero();
	init();
	for (int i = 0; i < N1; i++)
		for(int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				prim(NG + i, NG + j, NG + k, RHO) = primInit[i][j][k][RHO];
				prim(NG + i, NG + j, NG + k, UU) = primInit[i][j][k][UU];
				prim(NG + i, NG + j, NG + k, U1) = primInit[i][j][k][U1];
				prim(NG + i, NG + j, NG + k, U2) = primInit[i][j][k][U2];
				prim(NG + i, NG + j, NG + k, U3) = primInit[i][j][k][U3];
				prim(NG + i, NG + j, NG + k, B1) = primInit[i][j][k][B1];
				prim(NG + i, NG + j, NG + k, B2) = primInit[i][j][k][B2];
				prim(NG + i, NG + j, NG + k, B3) = primInit[i][j][k][B3];
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

	/*
	1)鬼化
		分成两边分别鬼化
	*/
	for(int epoch = 0; epoch < 100; epoch++)
	{
		auto start = clock();
		for (int i = NG - 1; i >= 0; i--)
			for (int j = NG - 1; j >= 0; j--)
				for (int k = NG - 1; k >= 0; k--)
				{
					prim(i, j, k, 0) = prim(i + 1, j + 1, k + 1, RHO) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i + NG, j + NG, k + NG).m.determinant());
					prim(i, j, k, 1) = prim(i + 1, j + 1, k + 1, UU) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i + NG, j + NG, k + NG).m.determinant());
					prim(i, j, k, 5) = prim(i + 1, j + 1, k + 1, B1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i + NG, j + NG, k + NG).m.determinant());
					prim(i, j, k, 3) = prim(i + 1, j + 1, k + 1, U2) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + (i + 1) * dx1, 2) + pow(X2min + (j + 1) * dx2, 2) + pow(X3min + (k + 1) * dx3, 2)));
					prim(i, j, k, 4) = prim(i + 1, j + 1, k + 1, U3) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + (i + 1) * dx1, 2) + pow(X2min + (j + 1) * dx2, 2) + pow(X3min + (k + 1) * dx3, 2)));
					prim(i, j, k, 6) = prim(i + 1, j + 1, k + 1, B2) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + (i + 1) * dx1, 2) + pow(X2min + (j + 1) * dx2, 2) + pow(X3min + (k + 1) * dx3, 2)));
					prim(i, j, k, 7) = prim(i + 1, j + 1, k + 1, B3) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + (i + 1) * dx1, 2) + pow(X2min + (j + 1) * dx2, 2) + pow(X3min + (k + 1) * dx3, 2)));
					prim(i, j, k, 2) = prim(i + 1, j + 1, k + 1, U1) * (1 + sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + (i + 1) * dx1, 2) + pow(X2min + (j + 1) * dx2, 2) + pow(X3min + (k + 1) * dx3, 2)));
				}

		for (int i = 0; i < NG; i++)
			for (int j = 0; j < NG; j++)
				for (int k = 0; k < NG; k++)
				{
					prim(i + 1, j + 1, k + 1, RHO) = prim(i, j, k, 0) * sqrt(-metricFuncField(i + NG, j + NG, k + NG).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
					prim(i + 1, j + 1, k + 1, UU) = prim(i, j, k, 1) * sqrt(-metricFuncField(i + NG, j + NG, k + NG).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
					prim(i + 1, j + 1, k + 1, B1) = prim(i, j, k, 5) * sqrt(-metricFuncField(i + NG, j + NG, k + NG).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
					prim(i + 1, j + 1, k + 1, U2) = prim(i, j, k, 3) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + i * dx1, 2) + pow(X2min + i * dx2, 2) + pow(X3min + i * dx3, 2)));
					prim(i + 1, j + 1, k + 1, U3) = prim(i, j, k, 4) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + i * dx1, 2) + pow(X2min + i * dx2, 2) + pow(X3min + i * dx3, 2)));
					prim(i + 1, j + 1, k + 1, B2) = prim(i, j, k, 6) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + i * dx1, 2) + pow(X2min + i * dx2, 2) + pow(X3min + i * dx3, 2)));
					prim(i + 1, j + 1, k + 1, B3) = prim(i, j, k, 7) * (1 - sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + i * dx1, 2) + pow(X2min + i * dx2, 2) + pow(X3min + i * dx3, 2)));
					prim(i + 1, j + 1, k + 1, U1) = prim(i, j, k, 2) * (1 + sqrt(pow(dx1, 2) + pow(dx2, 2) + pow(dx3, 2)) / sqrt(pow(X1min + i * dx1, 2) + pow(X2min + i * dx2, 2) + pow(X3min + i * dx3, 2)));
				}
		std::cout << prim << std::endl;
		// 可并发
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
		Eigen::Tensor<double, 4> fluxLLF1 = theta * fluxHLL1 + (1 - theta) * fluxTVDLF1;
		Eigen::Tensor<double, 4> fluxLLF2 = theta * fluxHLL2 + (1 - theta) * fluxTVDLF2;
		Eigen::Tensor<double, 4> fluxLLF3 = theta * fluxHLL3 + (1 - theta) * fluxTVDLF3;
		// 4.半步长迭代
		prim2con(prim, con);
		prim2src(prim, con, src);
#pragma omp parallel num_threads(2)
		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					auto c1max = max(0, cpR1(i, j, k), cpL1(i, j, k));
					auto c1min = -min(0, cnR1(i, j, k), cnL1(i, j, k));
					auto c2max = max(0, cpR2(i, j, k), cpL2(i, j, k));
					auto c2min = -min(0, cnR2(i, j, k), cnL2(i, j, k));
					auto c3max = max(0, cpR3(i, j, k), cpL3(i, j, k));
					auto c3min = -min(0, cnR3(i, j, k), cnL3(i, j, k));
					auto c1 = max(c1max, c1min);
					auto c2 = max(c2max, c2min);
					auto c3 = max(c3max, c3min);
					auto Delta_t = min(dx1 / (2 * c1), dx2 / (2 * c2), dx3 / (2 * c3));
					for (int l = 0; l < NPRIM; l++)
						conHalf(i, j, k, l) = con(i, j, k, l) + src(i, j, k, l)
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField1(i + 1, j, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxLLF1(i + 1, j, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField2(i, j + 1, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxLLF2(i, j + 1, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField3(i, j, k + 1).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxLLF3(i, j, k + 1, l) - fluxLLF1(i, j, k, l));
				}
		}
		con2prim(conHalf, primHalf);
		interpolate(primHalf);
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
		Eigen::Tensor<double, 4> fluxSmoothLLF1 = fluxLLF1;
		Eigen::Tensor<double, 4> fluxSmoothLLF2 = fluxLLF2;
		Eigen::Tensor<double, 4> fluxSmoothLLF3 = fluxLLF3;
		/*
		5) 平滑化
		*/
#pragma omp parallel num_threads(2)
		for (int i = 1; i < N1 - 1; i++)
		{
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
				{
					fluxSmoothLLF1(i, j, k, 6) = 0.125 * (2 * fluxLLF1(i, j, k, 6) + fluxLLF1(i, j + 1, k, 6) + fluxLLF1(i, j - 1, k, 6) - fluxLLF2(i, j, k, 5) - fluxLLF2(i, j + 1, k, 5) - fluxLLF2(i - 1, j, k, 5) - fluxLLF2(i - 1, j + 1, k, 5));
					fluxSmoothLLF2(i, j, k, 5) = 0.125 * (2 * fluxLLF2(i, j, k, 5) + fluxLLF2(i + 1, j, k, 5) + fluxLLF2(i - 1, j, k, 5) - fluxLLF1(i, j, k, 6) - fluxLLF1(i + 1, j, k, 6) - fluxLLF1(i, j - 1, k, 6) - fluxLLF1(i + 1, j - 1, k, 6));
					fluxSmoothLLF1(i, j, k, 7) = 0.125 * (2 * fluxLLF1(i, j, k, 7) + fluxLLF1(i, j, k + 1, B3) + fluxLLF1(i, j, k - 1, 7) - fluxLLF3(i, j, k, 5) - fluxLLF3(i, j, k + 1, B1) - fluxLLF3(i - 1, j, k, 5) - fluxLLF3(i - 1, j, k + 1, B1));
					fluxSmoothLLF3(i, j, k, 5) = 0.125 * (2 * fluxLLF3(i, j, k, 5) + fluxLLF3(i + 1, j, k, 5) + fluxLLF3(i - 1, j, k, 5) - fluxLLF1(i, j, k, 7) - fluxLLF1(i + 1, j, k, 7) - fluxLLF1(i, j, k - 1, 7) - fluxLLF1(i + 1, j, k - 1, 7));
					fluxSmoothLLF2(i, j, k, 7) = 0.125 * (2 * fluxLLF2(i, j, k, 7) + fluxLLF2(i, j, k + 1, B3) + fluxLLF2(i, j, k - 1, 7) - fluxLLF3(i, j, k, 6) - fluxLLF3(i, j, k + 1, B2) - fluxLLF3(i, j - 1, k, 6) - fluxLLF3(i, j - 1, k + 1, B2));
					fluxSmoothLLF3(i, j, k, 6) = 0.125 * (2 * fluxLLF3(i, j, k, 6) + fluxLLF3(i, j + 1, k, 6) + fluxLLF3(i, j - 1, k, 6) - fluxLLF2(i, j, k, 7) - fluxLLF2(i, j + 1, k, 7) - fluxLLF2(i, j, k - 1, 7) - fluxLLF2(i, j + 1, k - 1, 7));
				}
		}
		/*
		6) 整步迭代
		*/
#pragma omp parallel num_threads(2)
		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					auto c1max = max(0, cpR1(i, j, k), cpL1(i, j, k));
					auto c1min = -min(0, cnR1(i, j, k), cnL1(i, j, k));
					auto c2max = max(0, cpR2(i, j, k), cpL2(i, j, k));
					auto c2min = -min(0, cnR2(i, j, k), cnL2(i, j, k));
					auto c3max = max(0, cpR3(i, j, k), cpL3(i, j, k));
					auto c3min = -min(0, cnR3(i, j, k), cnL3(i, j, k));
					auto c1 = max(c1max, c1min);
					auto c2 = max(c2max, c2min);
					auto c3 = max(c3max, c3min);
					auto Delta_t = min(dx1 / (2 * c1), dx2 / (2 * c2), dx3 / (2 * c3));
					for (int l = 0; l < NPRIM; l++)
						con(i, j, k, l) = conHalf(i, j, k, l) + src(i, j, k, l)
						- Delta_t / (2 * dx1) * (sqrt(metricFuncField(i + 1, j, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxSmoothLLF1(i + 1, j, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * dx2) * (sqrt(metricFuncField(i, j + 1, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxSmoothLLF2(i, j + 1, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * dx3) * (sqrt(metricFuncField(i, j, k + 1).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxSmoothLLF3(i, j, k + 1, l) - fluxLLF1(i, j, k, l));
				}
		}
		con2prim(con, prim);
		totalTime += clock() - start;
		std::cout << "Time(ms): " << clock() - start << std::endl;
		if (epoch % 10 == 0) {
			char filename[13];
			sprintf(filename, "data%0.4d.bin", epoch / 10);
			write_bin(fopen(filename, "wb"));
		}
	}
	std::cout << "Total times(ms): " << totalTime << std::endl << "Average time(ms): " << totalTime / 100 << std::endl;
	return 0;
}
