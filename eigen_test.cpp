/*
bnuGRMHD ©️ 2024
Date: 2024/12/21
本文件是初代GRMHD的示例代码
TODO:
	1) 目前来说10x10x1的网格
	2)
	3)
	4)
*/

#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#include <ctime>
#include <unsupported/Eigen/CXX11/Tensor>
#include "Metric.h"
#include "omp.h"

using Eigen::Tensor;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using std::thread;
using std::ref;
using std::vector;

inline double max(double a, double b) { return a > b ? a : b; }
inline double max(double a, double b, double c) { return max(a, max(b, c)); }
inline double min(double a, double b) { return a < b ? a : b; }
inline double min(double a, double b, double c) { return min(a, min(b, c)); }

double MC(double a, double b, double c)
{
	if (abs(a) < abs(b) && abs(a) < abs(c) && b * c > 0)
		return a;
	else if (abs(b) < abs(a) && abs(b) < abs(c) && b * c > 0)
		return b;
	else if (abs(c) < abs(a) && abs(c) < abs(b) && b * c > 0)
		return c;
	else
		return 0;
}

int main()
{
	omp_set_num_threads(20);
	auto totalTime = 0;
	// 物理尺度
	double xStart = 10;		// 起始点x1坐标
	double yStart = 10;		// 起始点x2坐标
	double zStart = 0;		// 起始点x3坐标
	double L1 = 1;			// x1方向物理长度
	double L2 = 1;			// x2方向物理长度
	double L3 = 1;			// x3方向物理长度
	// 分辨率/格子数
	unsigned short n1 = 10;		// x方向格子数
	unsigned short n2 = 10;		// y方向格子数
	unsigned short n3 = 1;		// z方向格子数
	unsigned short nComp = 8;	// 分量个数
	unsigned short nGhost = 2;	// 单边鬼格数量
	// 常用常量
	double adiabaticIndex = 5.0 / 3.0;				// 多方指数/绝热指数
	double theta = 0.5;								// HHL流和TVDLF流混合参数

	Tensor<MetricComponent, 2> metricFunc(4, 4);												// 度规张量(0,2)型
	Tensor<MetricComponent, 3> metricDiff(4, 4, 4);												// 度规张量导数
	Tensor<Metric, 3> metricFuncField(n1 + 2 * nGhost, n2 + 2 * nGhost, n3 + 2 * nGhost);		// 度规场(0,2)型
	Tensor<Metric, 4> metricDiffField(n1 + 2 * nGhost, n2 + 2 * nGhost, n3 + 2 * nGhost, 4);	// 度规导数场
	Tensor<Metric, 3> metricFuncHalfField1(n1, n2, n3);											// 计算流时需要的半步长度规场(0,2)型
	Tensor<Metric, 3> metricFuncHalfField2(n1, n2, n3);											// 计算流时需要的半步长度规场(0,2)型
	Tensor<Metric, 3> metricFuncHalfField3(n1, n2, n3);											// 计算流时需要的半步长度规场(0,2)型
	// 主要量，对应传统GRMHD方程中的P(带鬼格)
	Tensor<double, 4> prim(n1 + 2 * nGhost, n2 + 2 * nGhost, n3 + 2 * nGhost, nComp);
	Tensor<double, 4> primHalf(n1, n2, n3, nComp);
	Tensor<double, 4> primL1(n1, n2, n3, nComp);
	Tensor<double, 4> primL2(n1, n2, n3, nComp);
	Tensor<double, 4> primL3(n1, n2, n3, nComp);
	Tensor<double, 4> primR1(n1, n2, n3, nComp);
	Tensor<double, 4> primR2(n1, n2, n3, nComp);
	Tensor<double, 4> primR3(n1, n2, n3, nComp);
	// 守恒量，对应传统GRMHD方程中的U(带鬼格)
	Tensor<double, 4> con(n1, n2, n3, nComp);
	Tensor<double, 4> conHalf(n1, n2, n3, nComp);
	Tensor<double, 4> conL1(n1, n2, n3, nComp);
	Tensor<double, 4> conL2(n1, n2, n3, nComp);
	Tensor<double, 4> conL3(n1, n2, n3, nComp);
	Tensor<double, 4> conR1(n1, n2, n3, nComp);
	Tensor<double, 4> conR2(n1, n2, n3, nComp);
	Tensor<double, 4> conR3(n1, n2, n3, nComp);
	// 流(flux)
	Tensor<double, 4> fluxL1(n1, n2, n3, nComp);
	Tensor<double, 4> fluxL2(n1, n2, n3, nComp);
	Tensor<double, 4> fluxL3(n1, n2, n3, nComp);
	Tensor<double, 4> fluxR1(n1, n2, n3, nComp);
	Tensor<double, 4> fluxR2(n1, n2, n3, nComp);
	Tensor<double, 4> fluxR3(n1, n2, n3, nComp);
	// HHL流
	Tensor<double, 4> fluxHLL1(n1, n2, n3, nComp);
	Tensor<double, 4> fluxHLL2(n1, n2, n3, nComp);
	Tensor<double, 4> fluxHLL3(n1, n2, n3, nComp);
	// TVDLF流
	Tensor<double, 4> fluxTVDLF1(n1, n2, n3, nComp);
	Tensor<double, 4> fluxTVDLF2(n1, n2, n3, nComp);
	Tensor<double, 4> fluxTVDLF3(n1, n2, n3, nComp);
	// 源(source)
	Tensor<double, 4> src(n1, n2, n3, nComp);
	Tensor<double, 4> srcL1(n1, n2, n3, nComp);
	Tensor<double, 4> srcL2(n1, n2, n3, nComp);
	Tensor<double, 4> srcL3(n1, n2, n3, nComp);
	Tensor<double, 4> srcR1(n1, n2, n3, nComp);
	Tensor<double, 4> srcR2(n1, n2, n3, nComp);
	Tensor<double, 4> srcR3(n1, n2, n3, nComp);
	// 特征速度(c_+)
	Tensor<double, 3> cpL1(n1, n2, n3);
	Tensor<double, 3> cpL2(n1, n2, n3);
	Tensor<double, 3> cpL3(n1, n2, n3);
	Tensor<double, 3> cpR1(n1, n2, n3);
	Tensor<double, 3> cpR2(n1, n2, n3);
	Tensor<double, 3> cpR3(n1, n2, n3);
	// 特征速度(c_-)
	Tensor<double, 3> cnL1(n1, n2, n3);
	Tensor<double, 3> cnL2(n1, n2, n3);
	Tensor<double, 3> cnL3(n1, n2, n3);
	Tensor<double, 3> cnR1(n1, n2, n3);
	Tensor<double, 3> cnR2(n1, n2, n3);
	Tensor<double, 3> cnR3(n1, n2, n3);
	double M = 1;	// 黑洞质量
	/*
	1.初始化
		1) 度规张量设置
		2) 度规场初始化
		3) 守恒量赋初值
	*/
	
	/*
	1) 度规张量设置
		创建Tensor<MetricComponent, 2>对象，每一个分量本质都是一个函数（类型为MetricComponent，数学地写就是R^3->R的映射定义见Metric.h）
		默认全为零分量（注意不是double 0而是预先定义的ZERO_COMPONENT，定义见Metric.h）
	*/
	metricFunc.setConstant(ZERO_COMPONENT);
	metricDiff.setConstant(ZERO_COMPONENT);
	metricFunc(0, 0) = [M](double x1, double x2, double x3) {return -1 + 2 * M / sqrt(x1 * x1 + x2 * x2 + x3 * x3); };
	metricFunc(1, 1) = [M](double x1, double x2, double x3) {return 1 - 2 * M / sqrt(x1 * x1 + x2 * x2 + x3 * x3); };
	metricFunc(2, 2) = [](double x1, double x2, double x3) {return x1 * x1 + x2 * x2 + x3 * x3; };
	metricFunc(3, 3) = [](double x1, double x2, double x3) {return x1 * x1 + x2 * x2; };
	metricDiff(0, 0, 0) = [M](double x1, double x2, double x3) {return -M / x1 * x1 + x2 * x2 + x3 * x3; };
	metricDiff(1, 1, 0) = [M](double x1, double x2, double x3) {return M / x1 * x1 + x2 * x2 + x3 * x3; };


	for (int i = 0; i < n1 + 2 * nGhost; i++)
		for (int j = 0; j < n2 + 2 * nGhost; j++)
			for (int k = 0; k < n3 + 2 * nGhost; k++)
				for(int row = 0; row < 4; row++)
					for (int col = 0; col < 4; col++)
						metricFuncField(i, j, k).m(row, col) = metricFunc(row, col)(xStart + i * L1 / n1, yStart + j * L2 / n2, zStart + k * L3 / n3);
	
	for (int i = 0; i < n1 + 2 * nGhost; i++)
		for (int j = 0; j < n2 + 2 * nGhost; j++)
			for (int k = 0; k < n3 + 2 * nGhost; k++)
				for (int l = 0; l < 4; l++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
							metricDiffField(i, j, k, l).m(row, col) = metricDiff(row, col, l)(xStart + i * L1 / n1, yStart + j * L2 / n2, zStart + k * L3 / n3);
	
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			for (int k = 0; k < n3; k++)
				for (int row = 0; row < 4; row++)
					for (int col = 0; col < 4; col++)
					{
						metricFuncHalfField1(i, j, k).m(row, col) = metricFunc(row, col)(xStart + (2 * i + 3) * L1 / (2 * n1), yStart + (j + 2) * L2 / n2, zStart + (k + 2) * L3 / n3);
						metricFuncHalfField2(i, j, k).m(row, col) = metricFunc(row, col)(xStart + (i + 2) * L1 / n1, yStart + (2 * j + 3) * L2 / (2 * n2), zStart + (k + 2) * L3 / n3);
						metricFuncHalfField3(i, j, k).m(row, col) = metricFunc(row, col)(xStart + (i + 2) * L1 / n1, yStart + (j + 2) * L2 / n2, zStart + (2 * k + 3) * L3 / (2 * n3));
					}

	prim.setZero();

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
		for (int i = nGhost - 1; i >= 0; i--)
			for (int j = nGhost - 1; j >= 0; j--)
				for (int k = nGhost - 1; k >= 0; k--)
				{
					prim(i, j, k, 0) = prim(i + 1, j + 1, k + 1, 0) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i + nGhost, j + nGhost, k + nGhost).m.determinant());
					prim(i, j, k, 1) = prim(i + 1, j + 1, k + 1, 1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i + nGhost, j + nGhost, k + nGhost).m.determinant());
					prim(i, j, k, 5) = prim(i + 1, j + 1, k + 1, 5) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i + nGhost, j + nGhost, k + nGhost).m.determinant());
					prim(i, j, k, 3) = prim(i + 1, j + 1, k + 1, 3) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + (i + 1) * L1 / n1, 2) + pow(yStart + (j + 1) * L2 / n2, 2) + pow(zStart + (k + 1) * L3 / n3, 2)));
					prim(i, j, k, 4) = prim(i + 1, j + 1, k + 1, 4) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + (i + 1) * L1 / n1, 2) + pow(yStart + (j + 1) * L2 / n2, 2) + pow(zStart + (k + 1) * L3 / n3, 2)));
					prim(i, j, k, 6) = prim(i + 1, j + 1, k + 1, 6) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + (i + 1) * L1 / n1, 2) + pow(yStart + (j + 1) * L2 / n2, 2) + pow(zStart + (k + 1) * L3 / n3, 2)));
					prim(i, j, k, 7) = prim(i + 1, j + 1, k + 1, 7) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + (i + 1) * L1 / n1, 2) + pow(yStart + (j + 1) * L2 / n2, 2) + pow(zStart + (k + 1) * L3 / n3, 2)));
					prim(i, j, k, 2) = prim(i + 1, j + 1, k + 1, 2) * (1 + sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + (i + 1) * L1 / n1, 2) + pow(yStart + (j + 1) * L2 / n2, 2) + pow(zStart + (k + 1) * L3 / n3, 2)));
				}

		for (int i = 0; i < nGhost; i++)
			for (int j = 0; j < nGhost; j++)
				for (int k = 0; k < nGhost; k++)
				{
					prim(i + 1, j + 1, k + 1, 0) = prim(i, j, k, 0) * sqrt(-metricFuncField(i + nGhost, j + nGhost, k + nGhost).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
					prim(i + 1, j + 1, k + 1, 1) = prim(i, j, k, 1) * sqrt(-metricFuncField(i + nGhost, j + nGhost, k + nGhost).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
					prim(i + 1, j + 1, k + 1, 5) = prim(i, j, k, 5) * sqrt(-metricFuncField(i + nGhost, j + nGhost, k + nGhost).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
					prim(i + 1, j + 1, k + 1, 3) = prim(i, j, k, 3) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
					prim(i + 1, j + 1, k + 1, 4) = prim(i, j, k, 4) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
					prim(i + 1, j + 1, k + 1, 6) = prim(i, j, k, 6) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
					prim(i + 1, j + 1, k + 1, 7) = prim(i, j, k, 7) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
					prim(i + 1, j + 1, k + 1, 2) = prim(i, j, k, 2) * (1 + sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
				}

		/*
		2.插值
		*/
		auto interpolate = [L1, L2, L3, n1, n2, n3, nComp, &primL1, &primL2, &primL3, &primR1, &primR2, &primR3](Tensor<double, 4> prim) {
			for (int i = 0; i < n1; i++)
				for (int j = 0; j < n2; j++)
					for (int k = 0; k < n3; k++)
						for (int index = 0; index < nComp; index++)
						{
							primL1(i, j, k, index) = prim(i + 2, j + 2, k + 2, index) - MC((prim(i + 3, j + 2, k + 2, index) - prim(i + 1, j + 2, k + 2, index)) / (2 * L1 / n1),
								2 * (prim(i + 3, j + 2, k + 2, index) - prim(i + 2, j + 2, k + 2, index)) / (L1 / n1),
								2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 1, j + 2, k + 2, index)) / (L1 / n1)) * L1 / n1 / 2;
							primR1(i, j, k, index) = prim(i + 1, j + 2, k + 2, index) + MC((prim(i + 2, j + 2, k + 2, index) - prim(i, j + 2, k + 2, index)) / (2 * L1 / n1),
								2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 1, j + 2, k + 2, index)) / (L1 / n1),
								2 * (prim(i + 1, j + 2, k + 2, index) - prim(i, j + 2, k + 2, index)) / (L1 / n1)) * L1 / n1 / 2;

							primL2(i, j, k, index) = prim(i + 2, j + 3, k + 2, index) - MC((prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 1, k + 2, index)) / (2 * L2 / n2),
								2 * (prim(i + 2, j + 3, k + 2, index) - prim(i + 2, j + 2, k + 2, index)) / (L2 / n2),
								2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 1, k + 2, index)) / (L2 / n2)) * L2 / n2 / 2;
							primR2(i, j, k, index) = prim(i + 2, j + 2, k + 2, index) + MC((prim(i + 2, j + 3, k + 2, index) - prim(i + 2, j, k + 2, index)) / (2 * L2 / n2),
								2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 1, k + 2, index)) / (L2 / n2),
								2 * (prim(i + 2, j + 1, k + 2, index) - prim(i + 2, j, k + 2, index)) / (L2 / n2)) * L2 / n2 / 2;

							primL3(i, j, k, index) = prim(i + 2, j + 2, k + 2, index) - MC((prim(i + 2, j + 2, k + 3, index) - prim(i + 2, j + 2, k + 1, index)) / (2 * L3 / n3),
								2 * (prim(i + 2, j + 2, k + 3, index) - prim(i + 2, j + 2, k + 2, index)) / (L3 / n3),
								2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 2, k + 1, index)) / (L3 / n3)) * L3 / n3 / 2;
							primR3(i, j, k, index) = prim(i + 2, j + 2, k + 1, index) + MC((prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 2, k, index)) / (2 * L3 / n3),
								2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 2, k + 1, index)) / (L3 / n3),
								2 * (prim(i + 2, j + 2, k + 1, index) - prim(i + 2, j + 2, k, index)) / (L3 / n3)) * L3 / n3 / 2;
						}
			};
		/*
		3) 计算流(flux)
		*/
		auto prim2con = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField](Tensor<double, 4> prim, Tensor<double, 4>& con) {
			auto threadFunc = [prim, &con, nGhost, adiabaticIndex, &metricFuncField](int i, int j, int k) {
				Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
				Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
				auto dot = [i, j, k, nGhost, &metricFuncField](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
				auto square = [dot](Vector3d vec) { return dot(vec, vec); };
				double Gamma = 1 / sqrt(1 - square(v));
				con(i, j, k, 0) = Gamma * prim(i, j, k, 0);
				con(i, j, k, 1) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) - prim(i, j, k, 1) + 0.5 * (square(B) * (1 + square(v) - pow(dot(B, v), 2))) - Gamma * prim(i, j, k, 0);
				con(i, j, k, 2) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) * prim(i, j, k, 2) + square(B) * prim(i, j, k, 2) - dot(B, v) * prim(i, j, k, 5);
				con(i, j, k, 3) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) * prim(i, j, k, 3) + square(B) * prim(i, j, k, 3) - dot(B, v) * prim(i, j, k, 6);
				con(i, j, k, 4) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) * prim(i, j, k, 4) + square(B) * prim(i, j, k, 4) - dot(B, v) * prim(i, j, k, 7);
				con(i, j, k, 5) = prim(i, j, k, 5);
				con(i, j, k, 6) = prim(i, j, k, 6);
				con(i, j, k, 7) = prim(i, j, k, 7);
				};

			for (int i = 0; i < n1; i++)
				for (int j = 0; j < n2; j++)
					for (int k = 0; k < n3; k++)
						threadFunc(i, j, k);
			};

		auto con2prim = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField](Tensor<double, 4> con, Tensor<double, 4>& prim) {
			auto max_iter = 5;
			auto tol = 0.01;
			for (int i = 0; i < n1; i++)
			{
				for (int j = 0; j < n2; j++)
				{
					for (int k = 0; k < n3; k++)
					{
						Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
						Vector3d B{ con(i, j, k, 5) ,con(i, j, k, 6) ,con(i, j, k, 7) };
						auto D = con(i, j, k, 0);
						auto tau = con(i, j, k, 1);
						auto dot = [i, j, k, nGhost, &metricFuncField](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
						auto square = [dot](Vector3d vec) { return dot(vec, vec); };
						auto f = [D, tau, S, B, dot, square, adiabaticIndex](double x) {
							auto Gamma = 1 / sqrt(1 - square(S + dot(S, B) * B / x) / pow(x + square(B), 2));
							return x - (adiabaticIndex - 1) / adiabaticIndex * (x - Gamma * D) / pow(Gamma, 2) - tau - D + square(B) - 0.5 * (square(B / Gamma) + pow(dot(S, B), 2) / pow(x, 2));
							};
						auto df = [D, tau, S, B, dot, square, adiabaticIndex](double x) {
							auto Gamma = 1 / sqrt(1 - square(S + dot(S, B) * B / x) / pow(x + square(B), 2));
							return 1 + ((2 * pow(dot(S, B), 2) / pow(x, 3) -
								square(B) * (2 * dot(S, B) * dot(B, (S + (B * dot(S, B)) / x)) / x)) /
								(pow(x, 2) * pow(square(B) + x, 2)) +
								(2 * square(S + (B * dot(S, B)) / x)) /
								pow(square(B) + x, 3)) / 2. -
								((-1 + adiabaticIndex) * (1 - square(S + (B * dot(S, B)) / x) /
									pow(square(B) + x, 2))) *
								(1 + (D * ((2 * dot(S, B) * dot(B, (S + (B * dot(S, B)) / x)) /
									(pow(x, 2) * pow(square(B) + x, 2)) +
									(2 * square(S + (B * dot(S, B)) / x)) /
									pow(square(B) + x, 3)) /
									(2. * pow(1 - square(S + (B * dot(S, B)) / x) /
										pow(square(B) + x, 2), 1.5))))) / adiabaticIndex -
								((-1 + adiabaticIndex) * ((2 * dot(S, B) * dot(B, (S + (B * dot(S, B)) / x))) /
									(pow(x, 2) * pow(square(B) + x, 2)) +
									(2 * square(S + (B * dot(S, B)) / x)) /
									pow(square(B) + x, 3)) *
									(x - D /
										sqrt(1 - square(S + (B * dot(S, B)) / x) /
											pow(square(B) + x, 2)))) / adiabaticIndex;
							};
						double x0 = 0.1;
						for (int iter = 0; iter < max_iter; iter++)
						{
							auto x1 = x0 - f(x0) / df(x0); // 牛顿迭代公式
							if (abs(x1 - x0) < tol)
								break;
							x0 = x1;
						}

						auto Gamma = 1 / sqrt(1 - square(S + dot(S, B) * B / x0) / pow(x0 + square(B), 2));
						prim(i, j, k, 0) = D / Gamma;
						prim(i, j, k, 1) = (adiabaticIndex - 1) / adiabaticIndex * (x0 - Gamma * D) / pow(Gamma, 2);
						prim(i, j, k, 2) = (S(0) + dot(S, B) * B(0) / x0) / (x0 + square(B));
						prim(i, j, k, 3) = (S(1) + dot(S, B) * B(1) / x0) / (x0 + square(B));
						prim(i, j, k, 4) = (S(2) + dot(S, B) * B(2) / x0) / (x0 + square(B));
						prim(i, j, k, 5) = B(0);
						prim(i, j, k, 6) = B(1);
						prim(i, j, k, 7) = B(2);
					}
				}
			}
			};

		auto prim2flux = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField](Tensor<double, 4> prim, Tensor<double, 4> con, Tensor<double, 4>& flux, short comp) {

			for (int i = 0; i < n1; i++)
				for (int j = 0; j < n2; j++)
					for (int k = 0; k < n3; k++)
					{
						Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
						Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
						Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
						auto dot = [i, j, k, nGhost, &metricFuncField](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
						auto square = [dot](Vector3d vec) { return dot(vec, vec); };
						double Gamma = 1 / sqrt(1 - square(v));
						auto W = S * (metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() * v).transpose() + (prim(i, j, k, 1) + 0.5 * (square(B) * (1 - square(v)) + pow(dot(B, v), 2))) * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(B, v) * v * B.transpose();
						flux(i, j, k, 0) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, comp + 2) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp)) * con(i, j, k, 0);
						flux(i, j, k, 1) = metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * (con(i, j, k, 2 + comp) - prim(i, j, k, 2 + comp) * con(i, j, k, 0)) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp) * con(i, j, k, 1);
						flux(i, j, k, 2) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * W * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma())(comp, 0) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp) * con(i, j, k, 2);
						flux(i, j, k, 3) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * W * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma())(comp, 1) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp) * con(i, j, k, 3);
						flux(i, j, k, 4) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * W * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma())(comp, 2) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp) * con(i, j, k, 4);
						flux(i, j, k, 5) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, 2 + comp) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp)) * con(i, j, k, 5) - (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, 2) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(0)) * con(i, j, k, 5 + comp);
						flux(i, j, k, 6) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, 2 + comp) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp)) * con(i, j, k, 6) - (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, 3) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(1)) * con(i, j, k, 5 + comp);
						flux(i, j, k, 7) = (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, 2 + comp) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(comp)) * con(i, j, k, 7) - (metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * prim(i, j, k, 4) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(2)) * con(i, j, k, 5 + comp);
					}
			};

		auto prim2src = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField, &metricDiffField](Tensor<double, 4> prim, Tensor<double, 4> con, Tensor<double, 4>& src) {
			for (int i = 0; i < n1; i++)
				for (int j = 0; j < n2; j++)
					for (int k = 0; k < n3; k++)
					{
						Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
						Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
						Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
						auto dot = [i, j, k, nGhost, &metricFuncField](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
						auto square = [dot](Vector3d vec) { return dot(vec, vec); };
						auto contract = [](Matrix3d A, Matrix3d B) {
							double sum = 0;
							for (int i = 0; i < 3; i++)
								for (int j = 0; j < 3; j++)
									sum += A(i, j) * B(i, j);
							return sum;
							};
						Matrix3d betaDiff;
						betaDiff << metricDiffField(i, j, k, 1).beta()(0), metricDiffField(i, j, k, 2).beta()(0), metricDiffField(i, j, k, 3).beta()(0),
							metricDiffField(i, j, k, 1).beta()(1), metricDiffField(i, j, k, 2).beta()(1), metricDiffField(i, j, k, 3).beta()(1),
							metricDiffField(i, j, k, 1).beta()(2), metricDiffField(i, j, k, 2).beta()(2), metricDiffField(i, j, k, 3).beta()(2);
						double Gamma = 1 / sqrt(1 - square(v));
						auto W = S * (metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() * v).transpose() + (prim(i, j, k, 1) + 0.5 * (square(B) * (1 - square(v)) + pow(dot(B, v), 2))) * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(B, v) * v * B.transpose();
						src(i, j, k, 1) = 0.5 * contract(W, (metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(0) * metricDiffField(i, j, k, 1).gamma() + metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(1) * metricDiffField(i, j, k, 2).gamma() + metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(2) * metricDiffField(i, j, k, 3).gamma()))
							+ contract(W * metricFuncField(i, j, k).gamma(), betaDiff)
							- (metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() * S)(0) * metricDiffField(i, j, k, 1).alpha() - (metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() * S)(1) * metricDiffField(i, j, k, 2).alpha() - (metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma().inverse() * S)(2) * metricDiffField(i, j, k, 3).alpha();
						src(i, j, k, 2) = 0.5 * metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * contract(W, metricDiffField(i, j, k, 1).gamma()) + dot(S, metricDiffField(i, j, k, 1).beta()) - (con(i, j, k, 0) + con(i, j, k, 1)) * metricDiffField(i, j, k, 1).alpha();
						src(i, j, k, 3) = 0.5 * metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * contract(W, metricDiffField(i, j, k, 2).gamma()) + dot(S, metricDiffField(i, j, k, 2).beta()) - (con(i, j, k, 0) + con(i, j, k, 1)) * metricDiffField(i, j, k, 2).alpha();
						src(i, j, k, 4) = 0.5 * metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha() * contract(W, metricDiffField(i, j, k, 3).gamma()) + dot(S, metricDiffField(i, j, k, 3).beta()) - (con(i, j, k, 0) + con(i, j, k, 1)) * metricDiffField(i, j, k, 3).alpha();
					}
			};

		auto prim2c = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField](Tensor<double, 4> prim, Tensor<double, 3> c, Tensor<Metric, 3>& metricFuncHalfField, short sign, short comp) {
			for (int i = 0; i < n1; i++)
				for (int j = 0; j < n2; j++)
					for (int k = 0; k < n3; k++)
					{
						Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
						Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
						auto dot = [i, j, k, nGhost, &metricFuncField](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
						auto square = [dot](Vector3d vec) { return dot(vec, vec); };
						double Gamma = 1 / sqrt(1 - square(v));
						auto u0 = Gamma / metricFuncField(i + nGhost, j + nGhost, k + nGhost).alpha();
						Vector3d ui = { Gamma * (prim(i,j,k,2) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(0)),
										Gamma * (prim(i,j,k,3) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(1)),
										Gamma * (prim(i,j,k,4) - metricFuncField(i + nGhost, j + nGhost, k + nGhost).beta()(2)) };
						auto cs_square = adiabaticIndex * prim(i, j, k, 1) / (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1));
						auto cA_square = (square(B) * (1 - square(v)) + pow(dot(B, v), 2)) / (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1) + square(B) * (1 - square(v)) + pow(dot(B, v), 2));
						auto vf_square = cA_square + cs_square - cA_square * cs_square;
						auto sigmaf = (1 - vf_square) / vf_square;
						auto metricInv = metricFuncHalfField(i, j, k).m.inverse();
						c(i, j, k) = (metricInv(0, comp) - pow(sigmaf, 2) * u0 * ui(comp)) / (metricInv(0, 0) - pow(sigmaf, 2) * u0 * u0) + sign * sqrt(
							pow((metricInv(0, comp) - pow(sigmaf, 2) * u0 * ui(comp)) / (metricInv(0, 0) - pow(sigmaf, 2) * u0 * u0), 2)
							- (metricInv(comp, comp) - sigmaf * ui(comp) * ui(comp)) / (metricInv(0, 0) - sigmaf * u0 * u0)
						);
					}
			};

		auto calFluxHHL = [n1, n2, n3, nComp](Tensor<double, 3> cpL, Tensor<double, 3> cpR,
			Tensor<double, 3> cnL, Tensor<double, 3> cnR,
			Tensor<double, 4> conL, Tensor<double, 4> conR,
			Tensor<double, 4> fluxL, Tensor<double, 4> fluxR,
			Tensor<double, 4>& fluxHLL
			) {
				for (int i = 0; i < n1; i++)
					for (int j = 0; j < n2; j++)
						for (int k = 0; k < n3; k++)
							for (int l = 0; l < nComp; l++)
							{
								auto c_max = max(0, cpR(i, j, k), cpL(i, j, k));
								auto c_min = -min(0, cnR(i, j, k), cnL(i, j, k));
								fluxHLL(i, j, k, l) = (c_min * fluxR(i, j, k, l) + c_max * fluxL(i, j, k, l) - c_max * c_min * (conR(i, j, k, l) - conL(i, j, k, l))) / (c_max + c_min);
							}
			};

		auto calFluxTVDLF = [n1, n2, n3, nComp](Tensor<double, 3> cpL, Tensor<double, 3> cpR,
			Tensor<double, 3> cnL, Tensor<double, 3> cnR,
			Tensor<double, 4> conL, Tensor<double, 4> conR,
			Tensor<double, 4> fluxL, Tensor<double, 4> fluxR,
			Tensor<double, 4>& fluxTVDLF
			) {
				for (int i = 0; i < n1; i++)
					for (int j = 0; j < n2; j++)
						for (int k = 0; k < n3; k++)
							for (int l = 0; l < nComp; l++)
							{
								auto c_max = max(0, cpR(i, j, k), cpL(i, j, k));
								auto c_min = -min(0, cnR(i, j, k), cnL(i, j, k));
								auto c = max(c_max, c_min);
								fluxTVDLF(i, j, k, l) = 0.5 * (fluxR(i, j, k, l) + fluxL(i, j, k, l)) - 0.5 * c * (conR(i, j, k, l) - conL(i, j, k, l));
							}
			};

		auto basicCalc = [prim2con, prim2flux, prim2src, prim2c](Tensor<double, 4> prim, Tensor<double, 4>& con, Tensor<double, 4>& flux, Tensor<double, 4>& src, Tensor<double, 3>& cp, Tensor<double, 3>& cn, Tensor<Metric, 3>& metricFuncHalfField, short comp) {
			prim2con(prim, con);
			prim2flux(prim, con, flux, comp);
			prim2src(prim, con, src);
			prim2c(prim, cp, metricFuncHalfField, 1, comp);
			prim2c(prim, cn, metricFuncHalfField, -1, comp);
			};

		auto fluxCalc = [calFluxHHL, calFluxTVDLF](Tensor<double, 3> cpL, Tensor<double, 3> cpR, Tensor<double, 3> cnL, Tensor<double, 3> cnR, Tensor<double, 4> conL, Tensor<double, 4> conR, Tensor<double, 4> fluxL, Tensor<double, 4> fluxR, Tensor<double, 4>& fluxHLL, Tensor<double, 4>& fluxTVDLF) {
			calFluxHHL(cpL, cpR, cnL, cnR, conL, conR, fluxL, fluxR, fluxHLL);
			calFluxTVDLF(cpL, cpR, cnL, cnR, conL, conR, fluxL, fluxR, fluxTVDLF);
			};
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
		Tensor<double, 4> fluxLLF1 = theta * fluxHLL1 + (1 - theta) * fluxTVDLF1;
		Tensor<double, 4> fluxLLF2 = theta * fluxHLL2 + (1 - theta) * fluxTVDLF2;
		Tensor<double, 4> fluxLLF3 = theta * fluxHLL3 + (1 - theta) * fluxTVDLF3;
		// 4.半步长迭代
		prim2con(prim, con);
		prim2src(prim, con, src);
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < n3; k++)
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
					auto Delta_t = min(L1 / (2 * n1 * c1), L2 / (2 * n2 * c2), L3 / (2 * n3 * c3));
					#pragma omp parallel
					for (int l = 0; l < nComp; l++)
						// TODO: gamma存在半步长和整步长问题
						conHalf(i, j, k, l) = con(i, j, k, l) + src(i, j, k, l)
						- Delta_t / (2 * L1 / n1) * (sqrt(metricFuncHalfField1(i + 1, j, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxLLF1(i + 1, j, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * L2 / n2) * (sqrt(metricFuncHalfField2(i, j + 1, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxLLF2(i, j + 1, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * L3 / n3) * (sqrt(metricFuncHalfField3(i, j, k + 1).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxLLF3(i, j, k + 1, l) - fluxLLF1(i, j, k, l));
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
		Tensor<double, 4> fluxSmoothLLF1 = fluxLLF1;
		Tensor<double, 4> fluxSmoothLLF2 = fluxLLF2;
		Tensor<double, 4> fluxSmoothLLF3 = fluxLLF3;
		/*
		5) 平滑化
		*/
		for (int i = 1; i < n1 - 1; i++)
			for (int j = 1; j < n2 - 1; j++)
				for (int k = 1; k < n3 - 1; k++)
				{
					fluxSmoothLLF1(i, j, k, 6) = 0.125 * (2 * fluxLLF1(i, j, k, 6) + fluxLLF1(i, j + 1, k, 6) + fluxLLF1(i, j - 1, k, 6) - fluxLLF2(i, j, k, 5) - fluxLLF2(i, j + 1, k, 5) - fluxLLF2(i - 1, j, k, 5) - fluxLLF2(i - 1, j + 1, k, 5));
					fluxSmoothLLF2(i, j, k, 5) = 0.125 * (2 * fluxLLF2(i, j, k, 5) + fluxLLF2(i + 1, j, k, 5) + fluxLLF2(i - 1, j, k, 5) - fluxLLF1(i, j, k, 6) - fluxLLF1(i + 1, j, k, 6) - fluxLLF1(i, j - 1, k, 6) - fluxLLF1(i + 1, j - 1, k, 6));
					fluxSmoothLLF1(i, j, k, 7) = 0.125 * (2 * fluxLLF1(i, j, k, 7) + fluxLLF1(i, j, k + 1, 7) + fluxLLF1(i, j, k - 1, 7) - fluxLLF3(i, j, k, 5) - fluxLLF3(i, j, k + 1, 5) - fluxLLF3(i - 1, j, k, 5) - fluxLLF3(i - 1, j, k + 1, 5));
					fluxSmoothLLF3(i, j, k, 5) = 0.125 * (2 * fluxLLF3(i, j, k, 5) + fluxLLF3(i + 1, j, k, 5) + fluxLLF3(i - 1, j, k, 5) - fluxLLF1(i, j, k, 7) - fluxLLF1(i + 1, j, k, 7) - fluxLLF1(i, j, k - 1, 7) - fluxLLF1(i + 1, j, k - 1, 7));
					fluxSmoothLLF2(i, j, k, 7) = 0.125 * (2 * fluxLLF2(i, j, k, 7) + fluxLLF2(i, j, k + 1, 7) + fluxLLF2(i, j, k - 1, 7) - fluxLLF3(i, j, k, 6) - fluxLLF3(i, j, k + 1, 6) - fluxLLF3(i, j - 1, k, 6) - fluxLLF3(i, j - 1, k + 1, 6));
					fluxSmoothLLF3(i, j, k, 6) = 0.125 * (2 * fluxLLF3(i, j, k, 6) + fluxLLF3(i, j + 1, k, 6) + fluxLLF3(i, j - 1, k, 6) - fluxLLF2(i, j, k, 7) - fluxLLF2(i, j + 1, k, 7) - fluxLLF2(i, j, k - 1, 7) - fluxLLF2(i, j + 1, k - 1, 7));
				}
		/*
		6) 整步迭代
		*/
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < n3; k++)
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
					auto Delta_t = min(L1 / (2 * n1 * c1), L2 / (2 * n2 * c2), L3 / (2 * n3 * c3));
					#pragma omp parallel
					for (int l = 0; l < nComp; l++)
						con(i, j, k, l) = conHalf(i, j, k, l) + src(i, j, k, l)
						- Delta_t / (2 * L1 / n1) * (sqrt(metricFuncField(i + 1, j, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxSmoothLLF1(i + 1, j, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * L2 / n2) * (sqrt(metricFuncField(i, j + 1, k).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxSmoothLLF2(i, j + 1, k, l) - fluxLLF1(i, j, k, l))
						- Delta_t / (2 * L3 / n3) * (sqrt(metricFuncField(i, j, k + 1).gamma().determinant() / metricFuncField(i, j, k).gamma().determinant()) * fluxSmoothLLF3(i, j, k + 1, l) - fluxLLF1(i, j, k, l));
				}
		con2prim(con, prim);
		totalTime += clock() - start;
		std::cout << "Time(ms): " << clock() - start << std::endl;
	}
	std::cout << "Total times(ms): " << totalTime << std::endl << "Average time(ms): " << totalTime / 100 << std::endl;
	return 0;
}