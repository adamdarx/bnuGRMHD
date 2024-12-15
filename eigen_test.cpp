/*
bnuGRMHD ©️ 2024
Date: 2024/12/13
本文件是不加并行（OpenMP或MPI）GRMHD的示例代码
*/

#include <iostream>
#include <cmath>
#include <unsupported/Eigen/CXX11/Tensor>
#include "Metric.h"

using Eigen::Tensor;
using Eigen::Vector3d;
using Eigen::Matrix3d;

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
	// 物理尺度
	double adiabaticIndex = 5.0 / 3.0;
	double xStart = 10;		// 起始点x1坐标
	double yStart = 10;		// 起始点x2坐标
	double zStart = 0;		// 起始点x3坐标
	double L1 = 1;			// x1方向物理长度
	double L2 = 1;			// x2方向物理长度
	double L3 = 1;			// x3方向物理长度
	// 分辨率/格子数
	unsigned short n1 = 32;		// x方向格子数
	unsigned short n2 = 32;		// y方向格子数
	unsigned short n3 = 1;		// z方向格子数
	unsigned short nComp = 8;	// 分量个数
	unsigned short nGhost = 2;	// 单边鬼格数量
	
	Tensor<MetricComponent, 2> metricFunc(4, 4);												// 度规张量
	Tensor<MetricComponent, 3> metricDiff(4, 4, 4);												// 度规张量导数
	Tensor<Metric, 3> metricFuncField(n1 + 2 * nGhost, n2 + 2 * nGhost, n3 + 2 * nGhost);				// 度规场
	Tensor<Metric, 4> metricDiffField(n1 + 2 * nGhost, n2 + 2 * nGhost, n3 + 2 * nGhost, 4);				// 度规导数场
	// 主要量，对应传统GRMHD方程中的P(带鬼格)
	Tensor<double, 4> prim(n1 + 2 * nGhost, n2 + 2 * nGhost, n3 + 2 * nGhost, nComp);
	Tensor<double, 4> primL1(n1, n2, n3, nComp);
	Tensor<double, 4> primL2(n1, n2, n3, nComp);
	Tensor<double, 4> primL3(n1, n2, n3, nComp);
	Tensor<double, 4> primR1(n1, n2, n3, nComp);
	Tensor<double, 4> primR2(n1, n2, n3, nComp);
	Tensor<double, 4> primR3(n1, n2, n3, nComp);
	// 守恒量，对应传统GRMHD方程中的U(带鬼格)
	Tensor<double, 4> conL1(n1, n2, n3, nComp);
	Tensor<double, 4> conL2(n1, n2, n3, nComp);
	Tensor<double, 4> conL3(n1, n2, n3, nComp);
	Tensor<double, 4> conR1(n1, n2, n3, nComp);
	Tensor<double, 4> conR2(n1, n2, n3, nComp);
	Tensor<double, 4> conR3(n1, n2, n3, nComp);
	// 流(flux)
	Tensor<double, 5> fluxL1(n1, n2, n3, 3, nComp);
	Tensor<double, 5> fluxL2(n1, n2, n3, 3, nComp);
	Tensor<double, 5> fluxL3(n1, n2, n3, 3, nComp);
	Tensor<double, 5> fluxR1(n1, n2, n3, 3, nComp);
	Tensor<double, 5> fluxR2(n1, n2, n3, 3, nComp);
	Tensor<double, 5> fluxR3(n1, n2, n3, 3, nComp);
	// 源(source)
	Tensor<double, 4> srcL1(n1, n2, n3, nComp);
	Tensor<double, 4> srcL2(n1, n2, n3, nComp);
	Tensor<double, 4> srcL3(n1, n2, n3, nComp);
	Tensor<double, 4> srcR1(n1, n2, n3, nComp);
	Tensor<double, 4> srcR2(n1, n2, n3, nComp);
	Tensor<double, 4> srcR3(n1, n2, n3, nComp);
	Tensor<double, 4> xi(n1, n2, n3, 1);														// conv2prim过程需要用到的辅助变量
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


	prim.setZero();
	primL1.setZero();
	primL2.setZero();
	primL3.setZero();
	primR1.setZero();
	primR2.setZero();
	primR3.setZero();
	conL1.setZero();
	conL2.setZero();
	conL3.setZero();
	conR1.setZero();
	conR2.setZero();
	conR3.setZero();
	fluxL1.setZero();
	fluxL2.setZero();
	fluxL3.setZero();
	fluxR1.setZero();
	fluxR2.setZero();
	fluxR3.setZero();
	srcL1.setZero();
	srcL2.setZero();
	srcL3.setZero();
	srcR1.setZero();
	srcR2.setZero();
	srcR3.setZero();
	xi.setZero();
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
	for (int i = nGhost - 1; i >= 0; i--)
		for (int j = nGhost - 1; j >= 0; j--)
			for (int k = nGhost - 1; k >= 0; k--)
			{
				prim(i, j, k, 0) = prim(i + 1, j + 1, k + 1, 0) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				prim(i, j, k, 1) = prim(i + 1, j + 1, k + 1, 1) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				prim(i, j, k, 5) = prim(i + 1, j + 1, k + 1, 5) * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
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
				prim(i + 1, j + 1, k + 1, 0) = prim(i, j, k, 0) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				prim(i + 1, j + 1, k + 1, 1) = prim(i, j, k, 1) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				prim(i + 1, j + 1, k + 1, 5) = prim(i, j, k, 5) * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				prim(i + 1, j + 1, k + 1, 3) = prim(i, j, k, 3) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
				prim(i + 1, j + 1, k + 1, 4) = prim(i, j, k, 4) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
				prim(i + 1, j + 1, k + 1, 6) = prim(i, j, k, 6) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
				prim(i + 1, j + 1, k + 1, 7) = prim(i, j, k, 7) * (1 - sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
				prim(i + 1, j + 1, k + 1, 2) = prim(i, j, k, 2) * (1 + sqrt(pow(L1 / n1, 2) + pow(L2 / n2, 2) + pow(L3 / n3, 2)) / sqrt(pow(xStart + i * L1 / n1, 2) + pow(yStart + j * L2 / n2, 2) + pow(zStart + k * L3 / n3, 2)));
			}

	/*
	2.插值
	*/
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			for (int k = 0; k < n3; k++)
				for(int index = 0; index < nComp; index++)
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

	/*
	3) 计算流(flux)
	*/
	auto prim2con = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField](Tensor<double, 4> prim, Tensor<double, 4> &con) {
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < n3; k++)
				{
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
				}
	};

	auto prim2flux = [n1, n2, n3, nGhost, adiabaticIndex, &metricFuncField](Tensor<double, 4> prim, Tensor<double, 4> con, Tensor<double, 5> &flux) {
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < n3; k++)
					for(int l = 0; l < 3; l++)
					{
						Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
						Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
						Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
						auto dot = [i, j, k, nGhost, &metricFuncField](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metricFuncField(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
						auto square = [dot](Vector3d vec) { return dot(vec, vec); };
						double Gamma = 1 / sqrt(1 - square(v));
						auto W = S * (metricFuncField(i, j, k).gamma().inverse() * v).transpose() + (prim(i, j, k, 1) + 0.5 * (square(B) * (1 + square(v)) - pow(dot(B, v), 2))) * metricFuncField(i, j, k).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(B, v) * v * B.transpose();
						flux(i, j, k, l, 0) = (metricFuncField(i, j, k).alpha() * prim(i, j, k, l + 2) - metricFuncField(i, j, k).beta()(l)) * con(i, j, k, 0);
						flux(i, j, k, l, 1) = metricFuncField(i, j, k).alpha() * (con(i, j, k, 2 + l) - prim(i, j, k, 2 + l) * con(i, j, k, 0)) - metricFuncField(i, j, k).beta()(l) * con(i, j, k, 1);
						flux(i, j, k, l, 2) = (metricFuncField(i, j, k).alpha() * W * metricFuncField(i, j, k).gamma())(l, 0) - metricFuncField(i, j, k).beta()(l) * con(i, j, k, 2);
						flux(i, j, k, l, 3) = (metricFuncField(i, j, k).alpha() * W * metricFuncField(i, j, k).gamma())(l, 1) - metricFuncField(i, j, k).beta()(l) * con(i, j, k, 3);
						flux(i, j, k, l, 4) = (metricFuncField(i, j, k).alpha() * W * metricFuncField(i, j, k).gamma())(l, 2) - metricFuncField(i, j, k).beta()(l) * con(i, j, k, 4);
						flux(i, j, k, l, 5) = (metricFuncField(i, j, k).alpha() * prim(i, j, k, l + 2) - metricFuncField(i, j, k).beta()(l)) * con(i, j, k, 5) - (metricFuncField(i, j, k).alpha() * prim(i, j, k, 2) - metricFuncField(i, j, k).beta()(0)) * con(i, j, k, 5 + l);
						flux(i, j, k, l, 6) = (metricFuncField(i, j, k).alpha() * prim(i, j, k, l + 2) - metricFuncField(i, j, k).beta()(l)) * con(i, j, k, 6) - (metricFuncField(i, j, k).alpha() * prim(i, j, k, 3) - metricFuncField(i, j, k).beta()(1)) * con(i, j, k, 5 + l);
						flux(i, j, k, l, 7) = (metricFuncField(i, j, k).alpha() * prim(i, j, k, l + 2) - metricFuncField(i, j, k).beta()(l)) * con(i, j, k, 7) - (metricFuncField(i, j, k).alpha() * prim(i, j, k, 4) - metricFuncField(i, j, k).beta()(2)) * con(i, j, k, 5 + l);
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
					auto W = S * (metricFuncField(i, j, k).gamma().inverse() * v).transpose() + (prim(i, j, k, 1) + 0.5 * (square(B) * (1 + square(v)) - pow(dot(B, v), 2))) * metricFuncField(i, j, k).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(B, v) * v * B.transpose();
					src(i, j, k, 1) = 0.5 * contract(W, (metricFuncField(i, j, k).beta()(0) * metricDiffField(i, j, k, 1).gamma() + metricFuncField(i, j, k).beta()(1) * metricDiffField(i, j, k, 2).gamma() + metricFuncField(i, j, k).beta()(2) * metricDiffField(i, j, k, 3).gamma()))
						+ contract(W * metricFuncField(i,j,k).gamma(), betaDiff)
						- (metricFuncField(i, j, k).gamma().inverse() * S)(0) * metricDiffField(i, j, k, 1).alpha() - (metricFuncField(i, j, k).gamma().inverse() * S)(1) * metricDiffField(i, j, k, 2).alpha() - (metricFuncField(i, j, k).gamma().inverse() * S)(2) * metricDiffField(i, j, k, 3).alpha();
					src(i, j, k, 2) = 0.5 * metricFuncField(i, j, k).alpha() * contract(W, metricDiffField(i, j, k, 1).gamma()) + dot(S, metricDiffField(i, j, k, 1).beta()) - (con(i, j, k, 0) + con(i, j, k, 1)) * metricDiffField(i, j, k, 1).alpha();
					src(i, j, k, 3) = 0.5 * metricFuncField(i, j, k).alpha() * contract(W, metricDiffField(i, j, k, 2).gamma()) + dot(S, metricDiffField(i, j, k, 2).beta()) - (con(i, j, k, 0) + con(i, j, k, 1)) * metricDiffField(i, j, k, 2).alpha();
					src(i, j, k, 4) = 0.5 * metricFuncField(i, j, k).alpha() * contract(W, metricDiffField(i, j, k, 3).gamma()) + dot(S, metricDiffField(i, j, k, 3).beta()) - (con(i, j, k, 0) + con(i, j, k, 1)) * metricDiffField(i, j, k, 3).alpha();
				}
	};

	prim2con(primL1, conL1);
	prim2con(primL2, conL2);
	prim2con(primL3, conL3);
	prim2con(primR1, conR1);
	prim2con(primR2, conR2);
	prim2con(primR3, conR3);
	prim2flux(primL1, conL1, fluxL1);
	prim2flux(primL2, conL2, fluxL2);
	prim2flux(primL3, conL3, fluxL3);
	prim2flux(primR1, conR1, fluxR1);
	prim2flux(primR2, conR2, fluxR2);
	prim2flux(primR3, conR3, fluxR3);
	prim2src(primL1, conL1, srcL1);
	prim2src(primL2, conL2, srcL2);
	prim2src(primL3, conL3, srcL3);
	prim2src(primR1, conR1, srcR1);
	prim2src(primR2, conR2, srcR2);
	prim2src(primR3, conR3, srcR3);

	std::cout << conL1.slice(Eigen::array<Eigen::DenseIndex, 4>{1,1,0,0}, Eigen::array<Eigen::DenseIndex, 4>{1,1,1,8}) << std::endl;
	
	return 0;
}