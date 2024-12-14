/*
bnuGRMHD ©️ 2024
Date: 2024/12/13
本文件是不加并行（OpenMP或MPI）GRMHD的示例代码
*/

#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>

using Eigen::Tensor;

int main()
{
	// 物理尺度
	unsigned short Lx = 1;		// x方向物理长度
	unsigned short Ly = 1;		// y方向物理长度
	unsigned short Lz = 1;		// z方向物理长度
	// 分辨率/格子数
	unsigned short nx = 32;		// x方向格子数
	unsigned short ny = 32;		// y方向格子数
	unsigned short nz = 1;		// z方向格子数
	unsigned short convNum = 8;	// 守恒量个数
	Tensor<double, 4> conv(nx, ny, nz, convNum);	// 守恒量，对应传统GRMHD方程中的U
	//1. 守恒量赋初值（初值条件很大程度上决定了问题形状）
}