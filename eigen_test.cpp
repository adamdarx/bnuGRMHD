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
	// 物理尺度(笛卡尔坐标系)
	double adiabaticIndex = 5 / 3;
	double xStart = 10;		// 起始点x坐标
	double yStart = 10;		// 起始点x坐标
	double zStart = 0;		// 起始点x坐标
	double Lx = 1;			// x方向物理长度
	double Ly = 1;			// y方向物理长度
	double Lz = 1;			// z方向物理长度
	// 分辨率/格子数
	unsigned short nx = 32;		// x方向格子数
	unsigned short ny = 32;		// y方向格子数
	unsigned short nz = 1;		// z方向格子数
	unsigned short nComp = 8;	// 分量个数
	unsigned short nGhost = 2;	// 单边鬼格数量
	
	Tensor<MetricComponent, 2> metricFunc(4, 4);												// 度规张量
	Tensor<Metric, 3> metric(nx + 2 * nGhost, ny + 2 * nGhost, nz + 2 * nGhost);				// 度规场
	// 主要量，对应传统GRMHD方程中的P(带鬼格)
	Tensor<double, 4> prim(nx + 2 * nGhost, ny + 2 * nGhost, nz + 2 * nGhost, nComp);
	Tensor<double, 4> primLx(nx, ny, nz, nComp);
	Tensor<double, 4> primLy(nx, ny, nz, nComp);
	Tensor<double, 4> primLz(nx, ny, nz, nComp);
	Tensor<double, 4> primRx(nx, ny, nz, nComp);
	Tensor<double, 4> primRy(nx, ny, nz, nComp);
	Tensor<double, 4> primRz(nx, ny, nz, nComp);
	// 守恒量，对应传统GRMHD方程中的U(带鬼格)
	Tensor<double, 4> conLx(nx, ny, nz, nComp);
	Tensor<double, 4> conLy(nx, ny, nz, nComp);
	Tensor<double, 4> conLz(nx, ny, nz, nComp);
	Tensor<double, 4> conRx(nx, ny, nz, nComp);
	Tensor<double, 4> conRy(nx, ny, nz, nComp);
	Tensor<double, 4> conRz(nx, ny, nz, nComp);
	// 流(flux)
	Tensor<double, 5> fluxLx(nx, ny, nz, 3, nComp);
	Tensor<double, 5> fluxLy(nx, ny, nz, 3, nComp);
	Tensor<double, 5> fluxLz(nx, ny, nz, 3, nComp);
	Tensor<double, 5> fluxRx(nx, ny, nz, 3, nComp);
	Tensor<double, 5> fluxRy(nx, ny, nz, 3, nComp);
	Tensor<double, 5> fluxRz(nx, ny, nz, 3, nComp);
	// 源(source)
	Tensor<double, 4> sLx(nx, ny, nz, nComp);
	Tensor<double, 4> sLy(nx, ny, nz, nComp);
	Tensor<double, 4> sLz(nx, ny, nz, nComp);
	Tensor<double, 4> sRx(nx, ny, nz, nComp);
	Tensor<double, 4> sRy(nx, ny, nz, nComp);
	Tensor<double, 4> sRz(nx, ny, nz, nComp);
	Tensor<double, 4> xi(nx, ny, nz, 1);														// conv2prim过程需要用到的辅助变量
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
	metricFunc(0, 0) = [M](double x, double y, double z) {return -1 + 2 * M / sqrt(x * x + y * y + z * z); };
	metricFunc(1, 1) = [M](double x, double y, double z) {return 1 - 2 * M / sqrt(x * x + y * y + z * z); };
	metricFunc(2, 2) = [M](double x, double y, double z) {return x * x + y * y + z * z; };
	metricFunc(3, 3) = [M](double x, double y, double z) {return x * x + y * y; };


	for (int i = 0; i < nx + 2 * nGhost; i++)
		for (int j = 0; j < ny + 2 * nGhost; j++)
			for (int k = 0; k < nz + 2 * nGhost; k++)
				for(int row = 0; row < 4; row++)
					for (int col = 0; col < 4; col++)
						metric(i, j, k).m(row, col) = metricFunc(row, col)(xStart + i * Lx / nx, yStart + j * Ly / ny, zStart + k * Lz / nz);
	

	prim.setZero();
	primLx.setZero();
	primLy.setZero();
	primLz.setZero();
	primRx.setZero();
	primRy.setZero();
	primRz.setZero();
	conLx.setZero();
	conLy.setZero();
	conLz.setZero();
	conRx.setZero();
	conRy.setZero();
	conRz.setZero();
	fluxLx.setZero();
	fluxLy.setZero();
	fluxLz.setZero();
	fluxRx.setZero();
	fluxRy.setZero();
	fluxRz.setZero();
	sLx.setZero();
	sLy.setZero();
	sLz.setZero();
	sRx.setZero();
	sRy.setZero();
	sRz.setZero();
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
				prim(i, j, k, 0) = prim(i + 1, j + 1, k + 1, 0) * sqrt(-metric(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metric(i, j, k).m.determinant());
				prim(i, j, k, 1) = prim(i + 1, j + 1, k + 1, 1) * sqrt(-metric(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metric(i, j, k).m.determinant());
				prim(i, j, k, 5) = prim(i + 1, j + 1, k + 1, 5) * sqrt(-metric(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metric(i, j, k).m.determinant());
				prim(i, j, k, 3) = prim(i + 1, j + 1, k + 1, 3) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + (i + 1) * Lx / nx, 2) + pow(yStart + (j + 1) * Ly / ny, 2) + pow(zStart + (k + 1) * Lz / nz, 2)));
				prim(i, j, k, 4) = prim(i + 1, j + 1, k + 1, 4) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + (i + 1) * Lx / nx, 2) + pow(yStart + (j + 1) * Ly / ny, 2) + pow(zStart + (k + 1) * Lz / nz, 2)));
				prim(i, j, k, 6) = prim(i + 1, j + 1, k + 1, 6) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + (i + 1) * Lx / nx, 2) + pow(yStart + (j + 1) * Ly / ny, 2) + pow(zStart + (k + 1) * Lz / nz, 2)));
				prim(i, j, k, 7) = prim(i + 1, j + 1, k + 1, 7) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + (i + 1) * Lx / nx, 2) + pow(yStart + (j + 1) * Ly / ny, 2) + pow(zStart + (k + 1) * Lz / nz, 2)));
				prim(i, j, k, 2) = prim(i + 1, j + 1, k + 1, 2) * (1 + sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + (i + 1) * Lx / nx, 2) + pow(yStart + (j + 1) * Ly / ny, 2) + pow(zStart + (k + 1) * Lz / nz, 2)));
			}

	for (int i = 0; i < nGhost; i++)
		for (int j = 0; j < nGhost; j++)
			for (int k = 0; k < nGhost; k++)
			{
				prim(i + 1, j + 1, k + 1, 0) = prim(i, j, k, 0) * sqrt(-metric(i, j, k).m.determinant()) / sqrt(-metric(i + 1, j + 1, k + 1).m.determinant());
				prim(i + 1, j + 1, k + 1, 1) = prim(i, j, k, 1) * sqrt(-metric(i, j, k).m.determinant()) / sqrt(-metric(i + 1, j + 1, k + 1).m.determinant());
				prim(i + 1, j + 1, k + 1, 5) = prim(i, j, k, 5) * sqrt(-metric(i, j, k).m.determinant()) / sqrt(-metric(i + 1, j + 1, k + 1).m.determinant());
				prim(i + 1, j + 1, k + 1, 3) = prim(i, j, k, 3) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + i * Lx / nx, 2) + pow(yStart + j * Ly / ny, 2) + pow(zStart + k * Lz / nz, 2)));
				prim(i + 1, j + 1, k + 1, 4) = prim(i, j, k, 4) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + i * Lx / nx, 2) + pow(yStart + j * Ly / ny, 2) + pow(zStart + k * Lz / nz, 2)));
				prim(i + 1, j + 1, k + 1, 6) = prim(i, j, k, 6) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + i * Lx / nx, 2) + pow(yStart + j * Ly / ny, 2) + pow(zStart + k * Lz / nz, 2)));
				prim(i + 1, j + 1, k + 1, 7) = prim(i, j, k, 7) * (1 - sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + i * Lx / nx, 2) + pow(yStart + j * Ly / ny, 2) + pow(zStart + k * Lz / nz, 2)));
				prim(i + 1, j + 1, k + 1, 2) = prim(i, j, k, 2) * (1 + sqrt(pow(Lx / nx, 2) + pow(Ly / ny, 2) + pow(Lz / nz, 2)) / sqrt(pow(xStart + i * Lx / nx, 2) + pow(yStart + j * Ly / ny, 2) + pow(zStart + k * Lz / nz, 2)));
			}

	/*
	2.插值
	*/
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int k = 0; k < nz; k++)
				for(int index = 0; index < nComp; index++)
				{
					primLx(i, j, k, index) = prim(i + 2, j + 2, k + 2, index) - MC((prim(i + 3, j + 2, k + 2, index) - prim(i + 1, j + 2, k + 2, index)) / (2 * Lx / nx),
						2 * (prim(i + 3, j + 2, k + 2, index) - prim(i + 2, j + 2, k + 2, index)) / (Lx / nx),
						2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 1, j + 2, k + 2, index)) / (Lx / nx)) * Lx / nx / 2;
					primRx(i, j, k, index) = prim(i + 1, j + 2, k + 2, index) + MC((prim(i + 2, j + 2, k + 2, index) - prim(i, j + 2, k + 2, index)) / (2 * Lx / nx),
						2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 1, j + 2, k + 2, index)) / (Lx / nx),
						2 * (prim(i + 1, j + 2, k + 2, index) - prim(i, j + 2, k + 2, index)) / (Lx / nx)) * Lx / nx / 2;

					primLy(i, j, k, index) = prim(i + 2, j + 3, k + 2, index) - MC((prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 1, k + 2, index)) / (2 * Ly / ny),
						2 * (prim(i + 2, j + 3, k + 2, index) - prim(i + 2, j + 2, k + 2, index)) / (Ly / ny),
						2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 1, k + 2, index)) / (Ly / ny)) * Ly / ny / 2;
					primRy(i, j, k, index) = prim(i + 2, j + 2, k + 2, index) + MC((prim(i + 2, j + 3, k + 2, index) - prim(i + 2, j, k + 2, index)) / (2 * Ly / ny),
						2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 1, k + 2, index)) / (Ly / ny),
						2 * (prim(i + 2, j + 1, k + 2, index) - prim(i + 2, j, k + 2, index)) / (Ly / ny)) * Ly / ny / 2;

					primLz(i, j, k, index) = prim(i + 2, j + 2, k + 2, index) - MC((prim(i + 2, j + 2, k + 3, index) - prim(i + 2, j + 2, k + 1, index)) / (2 * Lz / nz),
						2 * (prim(i + 2, j + 2, k + 3, index) - prim(i + 2, j + 2, k + 2, index)) / (Lz / nz),
						2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 2, k + 1, index)) / (Lz / nz)) * Lz / nz / 2;
					primRz(i, j, k, index) = prim(i + 2, j + 2, k + 1, index) + MC((prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 2, k, index)) / (2 * Lz / nz),
						2 * (prim(i + 2, j + 2, k + 2, index) - prim(i + 2, j + 2, k + 1, index)) / (Lz / nz),
						2 * (prim(i + 2, j + 2, k + 1, index) - prim(i + 2, j + 2, k, index)) / (Lz / nz)) * Lz / nz / 2;
				}

	/*
	3) 计算流(flux)
	*/
	auto prim2con = [nx, ny, nz, nGhost, adiabaticIndex, &metric](Tensor<double, 4> prim, Tensor<double, 4> &con) {
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int k = 0; k < nz; k++)
				{
					Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
					Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
					auto dot = [i, j, k, nGhost, &metric](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metric(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
					auto square = [dot](Vector3d vec) { return dot(vec, vec); };
					double Gamma = 1 / sqrt(1 - square(v));
					con(i, j, k, 0) = Gamma * prim(i, j, k, 0);
					con(i, j, k, 1) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) - prim(i, j, k, 1) + 0.5 * (square(B) * (1 + square(v) - pow(dot(B,v),2))) - Gamma * prim(i, j, k, 0);
					con(i, j, k, 2) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) * prim(i, j, k, 2) + square(B) * prim(i, j, k, 2) - dot(B, v) * prim(i, j, k, 5);
					con(i, j, k, 3) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) * prim(i, j, k, 3) + square(B) * prim(i, j, k, 3) - dot(B, v) * prim(i, j, k, 6);
					con(i, j, k, 4) = (prim(i, j, k, 0) + adiabaticIndex / (adiabaticIndex - 1) * prim(i, j, k, 1)) * pow(Gamma, 2) * prim(i, j, k, 4) + square(B) * prim(i, j, k, 4) - dot(B, v) * prim(i, j, k, 7);
					con(i, j, k, 5) = prim(i, j, k, 5);
					con(i, j, k, 6) = prim(i, j, k, 6);
					con(i, j, k, 7) = prim(i, j, k, 7);
				}
	};
	auto prim2flux = [nx, ny, nz, nGhost, adiabaticIndex, &metric](Tensor<double, 4> prim, Tensor<double, 4> con, Tensor<double, 5> &flux) {
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int k = 0; k < nz; k++)
					for(int l = 0; l < 3; l++)
					{
						Vector3d v{ prim(i, j, k, 2) ,prim(i, j, k, 3) ,prim(i, j, k, 4) };
						Vector3d B{ prim(i, j, k, 5) ,prim(i, j, k, 6) ,prim(i, j, k, 7) };
						auto dot = [i, j, k, nGhost, &metric](Vector3d vecA, Vector3d vecB) {return double(vecA.transpose() * metric(i + nGhost, j + nGhost, k + nGhost).gamma() * vecB); };
						auto square = [dot](Vector3d vec) { return dot(vec, vec); };
						double Gamma = 1 / sqrt(1 - square(v));
						flux(i, j, k, l, 0) = (metric(i, j, k).alpha() * prim(i, j, k, l + 2) - metric(i, j, k).beta()(l)) * con(i, j, k, 0);
						flux(i, j, k, l, 1) = metric(i, j, k).alpha() * (con(i, j, k, 2 + l) - prim(i, j, k, 2 + l) * con(i, j, k, 0)) - metric(i, j, k).beta()(l) * con(i, j, k, 1);
						flux(i, j, k, l, 2) = 

					}
	};
	std::cout << metric(10, 10, 0).m << std::endl;
	std::cout << conLx.slice(Eigen::array<Eigen::DenseIndex, 4>{1,1,0,0}, Eigen::array<Eigen::DenseIndex, 4>{1,1,1,8}) << std::endl;
	std::cout << metric(10, 10, 0).alpha() << std::endl;
	return 0;
}