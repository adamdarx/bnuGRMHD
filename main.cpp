/*
bnuGRMHD ©️ 2025
Date: 2024/02/01
*/
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
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
		// 主要量，对应传统GRMHD方程中的P(带鬼格)
		prim.setZero();
		primHalf.setZero();
		primL1.setZero();
		primL2.setZero();
		primL3.setZero();
		primR1.setZero();
		primR2.setZero();
		primR3.setZero();
		// 守恒量，对应传统GRMHD方程中的U(带鬼格)
		con.setZero();
		conHalf.setZero();
		conL1.setZero();
		conL2.setZero();
		conL3.setZero();
		conR1.setZero();
		conR2.setZero();
		conR3.setZero();
		// 流(flux)
		fluxL1.setZero();
		fluxL2.setZero();
		fluxL3.setZero();
		fluxR1.setZero();
		fluxR2.setZero();
		fluxR3.setZero();
		// HHL流
		fluxHLL1.setZero();
		fluxHLL2.setZero();
		fluxHLL3.setZero();
		// TVDLF流
		fluxTVDLF1.setZero();
		fluxTVDLF2.setZero();
		fluxTVDLF3.setZero();
		// 源(source)
		src.setZero();
		srcL1.setZero();
		srcL2.setZero();
		srcL3.setZero();
		srcR1.setZero();
		srcR2.setZero();
		srcR3.setZero();
		// 特征速度(c_+)
		cpL1.setZero();
		cpL2.setZero();
		cpL3.setZero();
		cpR1.setZero();
		cpR2.setZero();
		cpR3.setZero();
		// 特征速度(c_-)
		cnL1.setZero();
		cnL2.setZero();
		cnL3.setZero();
		cnR1.setZero();
		cnR2.setZero();
		cnR3.setZero();
		alphaDiffField.setZero();

		for (int i = 0; i < N1 + 2 * NG; i++)
			for (int j = 0; j < N2 + 2 * NG; j++)
				for (int k = 0; k < N3 + 2 * NG; k++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
							metricFuncField(i, j, k).m(row, col) = metricFunc(row, col)(X1min + i * dx1, X2min + j * dx2, X3min + k * dx3);

		// 利用中心差分计算alpha的导数
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					alphaDiffField(i, j, k, 0) = 0;
					alphaDiffField(i, j, k, 1) = (metricFuncField(i + NG + 1, j + NG, k + NG).alpha() - metricFuncField(i + NG - 1, j + NG, k + NG).alpha()) / (2 * dx1);
					alphaDiffField(i, j, k, 2) = (metricFuncField(i + NG, j + NG + 1, k + NG).alpha() - metricFuncField(i + NG, j + NG - 1, k + NG).alpha()) / (2 * dx2);
					alphaDiffField(i, j, k, 3) = (metricFuncField(i + NG, j + NG, k + NG + 1).alpha() - metricFuncField(i + NG, j + NG, k + NG - 1).alpha()) / (2 * dx3);
				}
		for (int i = 0; i < N1 + 2 * NG; i++)
			for (int j = 0; j < N2 + 2 * NG; j++)
				for (int k = 0; k < N3 + 2 * NG; k++)
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
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
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
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					ksi(i, j, k) = (prim(i + NG, j + NG, k + NG, RHO) + gam / (gam - 1) * prim(i + NG, j + NG, k + NG, UU)) * (1 + pow(prim(i + NG, j + NG, k + NG, U1), 2) + pow(prim(i + NG, j + NG, k + NG, U2), 2) + pow(prim(i + NG, j + NG, k + NG, U3), 2));
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
	for(int epoch = 0; epoch < epochNum; epoch++)
	{
		auto start = clock();
		// 时间步长
		double Delta_t = 1;
		// 鬼化
		{
			for (int i = NG - 1; i >= 0; i--)
			{
				for (int j = NG - 1; j >= 0; j--)
				{
					for (int k = NG - 1; k >= 0; k--)
					{
						prim(i, j, k, RHO) = prim(i + 1, j + 1, k + 1, RHO) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, UU) = prim(i + 1, j + 1, k + 1, UU) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, B1) = prim(i + 1, j + 1, k + 1, B1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

						prim(i, j, k, U2) = prim(i + 1, j + 1, k + 1, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, U3) = prim(i + 1, j + 1, k + 1, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B2) = prim(i + 1, j + 1, k + 1, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B3) = prim(i + 1, j + 1, k + 1, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i, j, k, U1) = prim(i + 1, j + 1, k + 1, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
					for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
					{
						prim(i + 1, j + 1, k + 1, RHO) = prim(i, j, k, RHO) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, UU) = prim(i, j, k, UU) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, B1) = prim(i, j, k, B1) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

						prim(i + 1, j + 1, k + 1, U2) = prim(i, j, k, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, U3) = prim(i, j, k, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B2) = prim(i, j, k, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B3) = prim(i, j, k, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i + 1, j + 1, k + 1, U1) = prim(i, j, k, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
				}
				for (int j = NG + N2; j < 2 * NG + N2 - 1; j++)
				{
					for (int k = NG - 1; k >= 0; k--)
					{
						prim(i, j, k, RHO) = prim(i + 1, j + 1, k + 1, RHO) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, UU) = prim(i + 1, j + 1, k + 1, UU) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, B1) = prim(i + 1, j + 1, k + 1, B1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

						prim(i, j, k, U2) = prim(i + 1, j + 1, k + 1, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, U3) = prim(i + 1, j + 1, k + 1, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B2) = prim(i + 1, j + 1, k + 1, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B3) = prim(i + 1, j + 1, k + 1, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i, j, k, U1) = prim(i + 1, j + 1, k + 1, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
					for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
					{
						prim(i + 1, j + 1, k + 1, RHO) = prim(i, j, k, RHO) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, UU) = prim(i, j, k, UU) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, B1) = prim(i, j, k, B1) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

						prim(i + 1, j + 1, k + 1, U2) = prim(i, j, k, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, U3) = prim(i, j, k, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B2) = prim(i, j, k, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B3) = prim(i, j, k, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i + 1, j + 1, k + 1, U1) = prim(i, j, k, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
				}
			}

			for (int i = NG + N1; i < 2 * NG + N1 - 1; i++)
			{
				for (int j = NG - 1; j >= 0; j--)
				{
					for (int k = NG - 1; k >= 0; k--)
					{
						prim(i, j, k, RHO) = prim(i + 1, j + 1, k + 1, RHO) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, UU) = prim(i + 1, j + 1, k + 1, UU) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, B1) = prim(i + 1, j + 1, k + 1, B1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

						prim(i, j, k, U2) = prim(i + 1, j + 1, k + 1, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, U3) = prim(i + 1, j + 1, k + 1, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B2) = prim(i + 1, j + 1, k + 1, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B3) = prim(i + 1, j + 1, k + 1, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i, j, k, U1) = prim(i + 1, j + 1, k + 1, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
					for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
					{
						prim(i + 1, j + 1, k + 1, RHO) = prim(i, j, k, RHO) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, UU) = prim(i, j, k, UU) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, B1) = prim(i, j, k, B1) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

						prim(i + 1, j + 1, k + 1, U2) = prim(i, j, k, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, U3) = prim(i, j, k, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B2) = prim(i, j, k, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B3) = prim(i, j, k, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i + 1, j + 1, k + 1, U1) = prim(i, j, k, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
				}
				for (int j = NG + N2; j < 2 * NG + N2 - 1; j++)
				{
					for (int k = NG - 1; k >= 0; k--)
					{
						prim(i, j, k, RHO) = prim(i + 1, j + 1, k + 1, RHO) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, UU) = prim(i + 1, j + 1, k + 1, UU) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
						prim(i, j, k, B1) = prim(i + 1, j + 1, k + 1, B1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

						prim(i, j, k, U2) = prim(i + 1, j + 1, k + 1, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, U3) = prim(i + 1, j + 1, k + 1, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B2) = prim(i + 1, j + 1, k + 1, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i, j, k, B3) = prim(i + 1, j + 1, k + 1, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i, j, k, U1) = prim(i + 1, j + 1, k + 1, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
					for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
					{
						prim(i + 1, j + 1, k + 1, RHO) = prim(i, j, k, RHO) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, UU) = prim(i, j, k, UU) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
						prim(i + 1, j + 1, k + 1, B1) = prim(i, j, k, B1) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

						prim(i + 1, j + 1, k + 1, U2) = prim(i, j, k, U2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, U3) = prim(i, j, k, U3) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B2) = prim(i, j, k, B2) * (1 - dx1 / (X1min + (i + 1) * dx1));
						prim(i + 1, j + 1, k + 1, B3) = prim(i, j, k, B3) * (1 - dx1 / (X1min + (i + 1) * dx1));

						prim(i + 1, j + 1, k + 1, U1) = prim(i, j, k, U1) * (1 + dx1 / (X1min + (i + 1) * dx1));
					}
				}
			}
		}
		interpolate(prim);
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
		// 4.半步长迭代(prim中包含鬼格,需要单独使用函数)
		// prim2con
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					Eigen::Vector3d u{ prim(i + NG, j + NG, k + NG, U1) ,prim(i + NG, j + NG, k + NG, U2) ,prim(i + NG, j + NG, k + NG, U3) };
					Eigen::Vector3d B{ prim(i + NG, j + NG, k + NG, B1) ,prim(i + NG, j + NG, k + NG, B2) ,prim(i + NG, j + NG, k + NG, B3) };
					double Gamma = sqrt(1 + square(i, j, k, u));
					con(i, j, k, 0) = Gamma * prim(i + NG, j + NG, k + NG, RHO);
					con(i, j, k, 1) = (prim(i + NG, j + NG, k + NG, RHO) + gam / (gam - 1) * prim(i + NG, j + NG, k + NG, UU)) * pow(Gamma, 2) - prim(i + NG, j + NG, k + NG, UU) + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma) - pow(dot(i, j, k, B, u / Gamma), 2))) - Gamma * prim(i + NG, j + NG, k + NG, RHO);
					con(i, j, k, 2) = (prim(i + NG, j + NG, k + NG, RHO) + gam / (gam - 1) * prim(i + NG, j + NG, k + NG, UU)) * pow(Gamma, 2) * prim(i + NG, j + NG, k + NG, U1) + square(i, j, k, B) * prim(i + NG, j + NG, k + NG, U1) - dot(i, j, k, B, u / Gamma) * prim(i + NG, j + NG, k + NG, B1);
					con(i, j, k, 3) = (prim(i + NG, j + NG, k + NG, RHO) + gam / (gam - 1) * prim(i + NG, j + NG, k + NG, UU)) * pow(Gamma, 2) * prim(i + NG, j + NG, k + NG, U2) + square(i, j, k, B) * prim(i + NG, j + NG, k + NG, U2) - dot(i, j, k, B, u / Gamma) * prim(i + NG, j + NG, k + NG, B2);
					con(i, j, k, 4) = (prim(i + NG, j + NG, k + NG, RHO) + gam / (gam - 1) * prim(i + NG, j + NG, k + NG, UU)) * pow(Gamma, 2) * prim(i + NG, j + NG, k + NG, U3) + square(i, j, k, B) * prim(i + NG, j + NG, k + NG, U3) - dot(i, j, k, B, u / Gamma) * prim(i + NG, j + NG, k + NG, B3);
					con(i, j, k, 5) = prim(i + NG, j + NG, k + NG, B1);
					con(i, j, k, 6) = prim(i + NG, j + NG, k + NG, B2);
					con(i, j, k, 7) = prim(i + NG, j + NG, k + NG, B3);
				}

		// prim2src
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					Eigen::Vector3d u{ prim(i + NG, j + NG, k + NG, U1) ,prim(i + NG, j + NG, k + NG, U2) ,prim(i + NG, j + NG, k + NG, U3) };
					Eigen::Vector3d B{ prim(i + NG, j + NG, k + NG, B1) ,prim(i + NG, j + NG, k + NG, B2) ,prim(i + NG, j + NG, k + NG, B3) };
					Eigen::Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
					Eigen::Matrix3d betaDiff;
					betaDiff << metricDiffField(i + NG, j + NG, k + NG, 1).betaVec()(0), metricDiffField(i + NG, j + NG, k + NG, 2).betaVec()(0), metricDiffField(i + NG, j + NG, k + NG, 3).betaVec()(0),
						metricDiffField(i + NG, j + NG, k + NG, 1).betaVec()(1), metricDiffField(i + NG, j + NG, k + NG, 2).betaVec()(1), metricDiffField(i + NG, j + NG, k + NG, 3).betaVec()(1),
						metricDiffField(i + NG, j + NG, k + NG, 1).betaVec()(2), metricDiffField(i + NG, j + NG, k + NG, 2).betaVec()(2), metricDiffField(i + NG, j + NG, k + NG, 3).betaVec()(2);
					double Gamma = sqrt(1 + square(i, j, k, u));
					// W^{ij}
					Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim(i + NG, j + NG, k + NG, UU) + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma)) - pow(dot(i, j, k, B, u / Gamma), 2))) * metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(i, j, k, B, u / Gamma) * u / Gamma * B.transpose();
					src(i, j, k, 0) = 0;
					src(i, j, k, 1) = 0.5 * contract(W, (metricFuncField(i + NG, j + NG, k + NG).betaVec()(0) * metricDiffField(i + NG, j + NG, k + NG, 1).gamma() + metricFuncField(i + NG, j + NG, k + NG).betaVec()(1) * metricDiffField(i + NG, j + NG, k + NG,2).gamma() + metricFuncField(i + NG, j + NG, k + NG).betaVec()(2) * metricDiffField(i + NG, j + NG, k + NG, 3).gamma()))
						+ contract(W * metricFuncField(i + NG, j + NG, k + NG).gamma(), betaDiff)
						- (metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() * S)(0) * alphaDiffField(i, j, k, 1) - (metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() * S)(1) * alphaDiffField(i, j, k, 2) - (metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() * S)(2) * alphaDiffField(i, j, k, 3);
					src(i, j, k, 2) = 0.5 * metricFuncField(i + NG, j + NG, k + NG).alpha() * contract(W, metricDiffField(i + NG, j + NG, k + NG, 1).gamma()) + dot(i, j, k, S, metricDiffField(i + NG, j + NG, k + NG, 1).betaVec()) - (con(i, j, k, 0) + con(i, j, k, 1)) * alphaDiffField(i, j, k, 1);
					src(i, j, k, 3) = 0.5 * metricFuncField(i + NG, j + NG, k + NG).alpha() * contract(W, metricDiffField(i + NG, j + NG, k + NG, 2).gamma()) + dot(i, j, k, S, metricDiffField(i + NG, j + NG, k + NG, 2).betaVec()) - (con(i, j, k, 0) + con(i, j, k, 1)) * alphaDiffField(i, j, k, 2);
					src(i, j, k, 4) = 0.5 * metricFuncField(i + NG, j + NG, k + NG).alpha() * contract(W, metricDiffField(i + NG, j + NG, k + NG, 3).gamma()) + dot(i, j, k, S, metricDiffField(i + NG, j + NG, k + NG, 3).betaVec()) - (con(i, j, k, 0) + con(i, j, k, 1)) * alphaDiffField(i, j, k, 3);
					src(i, j, k, 5) = 0;
					src(i, j, k, 6) = 0;
					src(i, j, k, 7) = 0;
				}

		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					auto c1max = max(0, cpR1(i, j, k), cpL1(i, j, k));
					auto c1min = -min(-0, cnR1(i, j, k), cnL1(i, j, k));
					auto c2max = max(0, cpR2(i, j, k), cpL2(i, j, k));
					auto c2min = -min(-0, cnR2(i, j, k), cnL2(i, j, k));
					auto c3max = max(0, cpR3(i, j, k), cpL3(i, j, k));
					auto c3min = -min(-0, cnR3(i, j, k), cnL3(i, j, k));
					auto c1 = abs(max(c1max, c1min));
					auto c2 = abs(max(c2max, c2min));
					auto c3 = abs(max(c3max, c3min));
					Delta_t = min(min(dx1 / (2 * c1), dx2 / (2 * c2), dx3 / (2 * c3)), Delta_t);
				}

		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
					for (int l = 0; l < 8; l++)
						conHalf(i, j, k, l) = con(i, j, k, l) + src(i, j, k, l) * Delta_t / 2
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField1(i + 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1(i + 1, j, k, l) - sqrt(metricFuncHalfField1(i - 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1(i, j, k, l))
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField2(i, j + 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2(i, j + 1, k, l) - sqrt(metricFuncHalfField1(i, j - 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2(i, j, k, l))
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField3(i, j, k + 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3(i, j, k + 1, l) - sqrt(metricFuncHalfField1(i, j, k - 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3(i, j, k, l));
		
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
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
				{
					fluxSmoothLLF1(i, j, k, 6) = 0.125 * (2 * fluxLLF1(i, j, k, 6) + fluxLLF1(i, j + 1, k, 6) + fluxLLF1(i, j - 1, k, 6) - fluxLLF2(i, j, k, 5) - fluxLLF2(i, j + 1, k, 5) - fluxLLF2(i - 1, j, k, 5) - fluxLLF2(i - 1, j + 1, k, 5));
					fluxSmoothLLF2(i, j, k, 5) = 0.125 * (2 * fluxLLF2(i, j, k, 5) + fluxLLF2(i + 1, j, k, 5) + fluxLLF2(i - 1, j, k, 5) - fluxLLF1(i, j, k, 6) - fluxLLF1(i + 1, j, k, 6) - fluxLLF1(i, j - 1, k, 6) - fluxLLF1(i + 1, j - 1, k, 6));
					fluxSmoothLLF1(i, j, k, 7) = 0.125 * (2 * fluxLLF1(i, j, k, 7) + fluxLLF1(i, j, k + 1, 7) + fluxLLF1(i, j, k - 1, 7) - fluxLLF3(i, j, k, 5) - fluxLLF3(i, j, k + 1, 7) - fluxLLF3(i - 1, j, k, 5) - fluxLLF3(i - 1, j, k + 1, 7));
					fluxSmoothLLF3(i, j, k, 5) = 0.125 * (2 * fluxLLF3(i, j, k, 5) + fluxLLF3(i + 1, j, k, 5) + fluxLLF3(i - 1, j, k, 5) - fluxLLF1(i, j, k, 7) - fluxLLF1(i + 1, j, k, 7) - fluxLLF1(i, j, k - 1, 7) - fluxLLF1(i + 1, j, k - 1, 7));
					fluxSmoothLLF2(i, j, k, 7) = 0.125 * (2 * fluxLLF2(i, j, k, 7) + fluxLLF2(i, j, k + 1, 7) + fluxLLF2(i, j, k - 1, 7) - fluxLLF3(i, j, k, 6) - fluxLLF3(i, j, k + 1, 7) - fluxLLF3(i, j - 1, k, 6) - fluxLLF3(i, j - 1, k + 1, 7));
					fluxSmoothLLF3(i, j, k, 6) = 0.125 * (2 * fluxLLF3(i, j, k, 6) + fluxLLF3(i, j + 1, k, 6) + fluxLLF3(i, j - 1, k, 6) - fluxLLF2(i, j, k, 7) - fluxLLF2(i, j + 1, k, 7) - fluxLLF2(i, j, k - 1, 7) - fluxLLF2(i, j + 1, k - 1, 7));
				}
		/*
		6) 整步迭代
		*/
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
					for (int l = 0; l < 8; l++)
						conHalf(i, j, k, l) = con(i, j, k, l) + src(i, j, k, l) * Delta_t / 2
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField1(i + 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1(i + 1, j, k, l) - sqrt(metricFuncHalfField1(i - 1, j, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF1(i, j, k, l))
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField2(i, j + 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2(i, j + 1, k, l) - sqrt(metricFuncHalfField1(i, j - 1, k).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF2(i, j, k, l))
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField3(i, j, k + 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3(i, j, k + 1, l) - sqrt(metricFuncHalfField1(i, j, k - 1).gamma().determinant() / metricFuncField(i + NG, j + NG, k + NG).gamma().determinant()) * fluxLLF3(i, j, k, l));

		//con2prim(prim具有鬼格)
		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				for (int k = 0; k < N3; k++)
				{
					Eigen::Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
					Eigen::Vector3d B{ con(i, j, k, 5) ,con(i, j, k, 6) ,con(i, j, k, 7) };
					auto D = con(i, j, k, 0);
					auto tau = con(i, j, k, 1);
					auto x0 = ksi(i, j, k);
					for (int iter = 0; iter < max_iter; iter++)
					{
						auto x1 = x0 - f(i, j, k, D, tau, S, B, x0) / df(i, j, k, D, tau, S, B, x0); // 牛顿迭代公式
						if (abs((x1 - x0) / x0) < tol)
							break;
						x0 = x1;
					}
					if (ksi(i, j, k) <= 0 || isnan(x0))
						continue;
					ksi(i, j, k) = x0;
					auto Gamma = 1 / sqrt(1 - square(i, j, k, S + dot(i, j, k, S, B) * B / ksi(i, j, k)) / pow(ksi(i, j, k) + square(i, j, k, B), 2));
					prim(i + NG, j + NG, k + NG, RHO) = D / Gamma;
					prim(i + NG, j + NG, k + NG, UU) = (gam - 1) / gam * (ksi(i,j,k) - Gamma * D) / pow(Gamma, 2);
					prim(i + NG, j + NG, k + NG, U1) = (S(0) + dot(i, j, k, S, B) * B(0) / ksi(i, j, k)) / (ksi(i, j, k) + square(i, j, k, B)) * Gamma;
					prim(i + NG, j + NG, k + NG, U2) = (S(1) + dot(i, j, k, S, B) * B(1) / ksi(i, j, k)) / (ksi(i, j, k) + square(i, j, k, B)) * Gamma;
					prim(i + NG, j + NG, k + NG, U3) = (S(2) + dot(i, j, k, S, B) * B(2) / ksi(i, j, k)) / (ksi(i, j, k) + square(i, j, k, B)) * Gamma;
					prim(i + NG, j + NG, k + NG, B1) = B(0);
					prim(i + NG, j + NG, k + NG, B2) = B(1);
					prim(i + NG, j + NG, k + NG, B3) = B(2);
				}
			}
		}
		fix(prim);
		if (epoch % 100 == 0)
		{
			ofs << "-----------------------------Epoch: " << epoch << "-----------------------------" << std::endl;
			for (int i = 0; i < N1; i++)
				for (int j = 0; j < N2; j++)
					for (int k = 0; k < N3; k++)
						for (int l = 0; l < NPRIM; l++)
							ofs << "i: " << i << "\tj: " << j << "\tk: " << k << "\tValue: " << prim(i, j, k, l) << std::endl;
			ofs << "Time(ms): " << clock() - start << "\tPhysical Time: " << Delta_t << "\tTotal Physical Time: " << totalPhysicalTime << std::endl;
		}
		totalTime += clock() - start;
		totalPhysicalTime += Delta_t;
		std::cout << "Time(ms): " << clock() - start << "\tPhysical Time: " << Delta_t << "\tTotal Physical Time: " << totalPhysicalTime << std::endl;
		if (epoch % 100 == 0) {
			char filename[13];
			sprintf(filename, "data%0.4d.bin", epoch / 100);
			write_bin(fopen(filename, "wb"));
		}
	}
	std::cout << "Finished! Details can be found in grmhd.log. " << std::endl;
	ofs << "-----------------------------Epoch: " << epochNum << "-----------------------------" << std::endl;
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int l = 0; l < NPRIM; l++)
					ofs << "i: " << i << "\tj: " << j << "\tk: " << k << "\tValue: " << prim(i, j, k, l) << std::endl;
	ofs << "Total times(ms): " << totalTime << std::endl << "Average time(ms): " << totalTime / epochNum << std::endl;
	return 0;
}
