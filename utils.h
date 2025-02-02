#pragma once
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "Metric.h"
MetricComponent ZERO_COMPONENT = [](double x, double y, double z) { return 0; };
constexpr auto a = (0.9375);
constexpr auto h = (0.);
constexpr auto theta = 0.;								// HHL����TVDLF����ϲ���
constexpr auto NDIM = (4);
constexpr auto N1 = (32);
constexpr auto N2 = (16);
constexpr auto N3 = (16);
constexpr auto NG = (2);
constexpr auto PI = (3.14159265358979323846);
constexpr auto X1min = (0.19325057145871735);
constexpr auto X1max = (7.824046010856292);
constexpr auto X2min = (1.e-16);
constexpr auto X2max = (1. * PI);
constexpr auto X3min = (1.e-16);
constexpr auto X3max = (2. * PI);
constexpr auto R0 = (0.);
constexpr auto SMALL = (1.e-16);
constexpr auto isX1periodical = false;
constexpr auto isX2periodical = true;
constexpr auto isX3periodical = true;
constexpr auto cour = 0.8;

//FM_torus disk parameter
constexpr auto rin = (6.);
constexpr auto rmax = (12.);
constexpr auto beta = (100.);
constexpr auto gam = (5. / 3.);
constexpr auto kappa = (1.e-3);

unsigned short max_iter = 5;		// maximum of iteration
double tol = 1e-4;					// tolerance of root devation
auto epochNum = 10000;				// number of iteration epoch
//MKS grid
double Xgrid1[N1][N2][N3];
double Xgrid2[N1][N2][N3];
double Xgrid3[N1][N2][N3];
//MKS grid spacing
double dx1 = (X1max - X1min) / (N1 + 2 * NG);
double dx2 = (X2max - X2min) / (N2 + 2 * NG);
double dx3 = (X3max - X3min) / (N3 + 2 * NG);
//KS grid
double KS_coord1[N1][N2][N3];
double KS_coord2[N1][N2][N3];
double KS_coord3[N1][N2][N3];
//BL grid
double BL_coord1[N1][N2][N3];
double BL_coord2[N1][N2][N3];
double BL_coord3[N1][N2][N3];

//metric at grid point
double gdd_bl[N1][N2][N3][NDIM][NDIM];
double guu_bl[N1][N2][N3][NDIM][NDIM];
double gdet_bl[N1][N2][N3];              /*sqrt(-g_bl)*/
double gdd_ks[N1][N2][N3][NDIM][NDIM];
double guu_ks[N1][N2][N3][NDIM][NDIM];
double gdet_ks[N1][N2][N3];              /*sqrt(-g_ks)*/
double gdd_mks[N1][N2][N3][NDIM][NDIM];
double guu_mks[N1][N2][N3][NDIM][NDIM];
double gdet_mks[N1][N2][N3];             /*sqrt(-g_mks)*/

//Jacobian matrix at grid point
double J_bl2ks[N1][N2][N3][NDIM][NDIM];
double J_ks2bl[N1][N2][N3][NDIM][NDIM];
double J_ks2mks[N1][N2][N3][NDIM][NDIM];
double J_mks2ks[N1][N2][N3][NDIM][NDIM];

//primitive variables
constexpr auto NPRIM = (10);
double primInit[N1][N2][N3][NPRIM];
constexpr auto RHO = (0);
constexpr auto UU = (1);
constexpr auto U0 = (2);
constexpr auto U1 = (3);
constexpr auto U2 = (4);
constexpr auto U3 = (5);
constexpr auto B1 = (6);
constexpr auto B2 = (7);
constexpr auto B3 = (8);
constexpr auto BSQ = (9);
double A[N1 + 1][N2 + 1][N3 + 1];

//fix p
constexpr auto RHOMIN = (1.e-6);
constexpr auto UUMIN = (1.e-8);
constexpr auto SIGMAMAX = (50.);

// metric
Eigen::Tensor<MetricComponent, 2> metricFunc(4, 4);												// �ȹ�����(0,2)��
Eigen::Tensor<MetricComponent, 3> metricDiff(4, 4, 4);												// �ȹ���������
Eigen::Tensor<Metric, 3> metricFuncField(N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG);					// �ȹ泡(0,2)��
Eigen::Tensor<Metric, 4> metricDiffField(N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, 4);					// �ȹ浼����
Eigen::Tensor<double, 4> alphaDiffField(N1, N2, N3, 4);												// alpha������
Eigen::Tensor<Metric, 3> metricFuncHalfField1(N1, N2, N3);											// ������ʱ��Ҫ�İ벽���ȹ泡(0,2)��
Eigen::Tensor<Metric, 3> metricFuncHalfField2(N1, N2, N3);											// ������ʱ��Ҫ�İ벽���ȹ泡(0,2)��
Eigen::Tensor<Metric, 3> metricFuncHalfField3(N1, N2, N3);											// ������ʱ��Ҫ�İ벽���ȹ泡(0,2)��

// ��Ҫ������Ӧ��ͳGRMHD�����е�P(������)
Eigen::Tensor<double, 4> prim(N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NPRIM);
Eigen::Tensor<double, 4> primHalf(N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NPRIM);
Eigen::Tensor<double, 4> primL1(N1, N2, N3, NPRIM);
Eigen::Tensor<double, 4> primL2(N1, N2, N3, NPRIM);
Eigen::Tensor<double, 4> primL3(N1, N2, N3, NPRIM);
Eigen::Tensor<double, 4> primR1(N1, N2, N3, NPRIM);
Eigen::Tensor<double, 4> primR2(N1, N2, N3, NPRIM);
Eigen::Tensor<double, 4> primR3(N1, N2, N3, NPRIM);
// �غ�������Ӧ��ͳGRMHD�����е�U(������)
Eigen::Tensor<double, 4> con(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conHalf(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conL1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conL2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conL3(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conR1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conR2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> conR3(N1, N2, N3, 8);
// ��(flux)
Eigen::Tensor<double, 4> fluxL1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxL2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxL3(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxR1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxR2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxR3(N1, N2, N3, 8);
// HHL��
Eigen::Tensor<double, 4> fluxHLL1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxHLL2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxHLL3(N1, N2, N3, 8);
// TVDLF��
Eigen::Tensor<double, 4> fluxTVDLF1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxTVDLF2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxTVDLF3(N1, N2, N3, 8);
// �����
Eigen::Tensor<double, 4> fluxLLF1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxLLF2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxLLF3(N1, N2, N3, 8);
// �⻬�����
Eigen::Tensor<double, 4> fluxSmoothLLF1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxSmoothLLF2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> fluxSmoothLLF3(N1, N2, N3, 8);
// Դ(source)
Eigen::Tensor<double, 4> src(N1, N2, N3, 8);
Eigen::Tensor<double, 4> srcL1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> srcL2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> srcL3(N1, N2, N3, 8);
Eigen::Tensor<double, 4> srcR1(N1, N2, N3, 8);
Eigen::Tensor<double, 4> srcR2(N1, N2, N3, 8);
Eigen::Tensor<double, 4> srcR3(N1, N2, N3, 8);
// �����ٶ�(c_+)
Eigen::Tensor<double, 3> cpL1(N1, N2, N3);
Eigen::Tensor<double, 3> cpL2(N1, N2, N3);
Eigen::Tensor<double, 3> cpL3(N1, N2, N3);
Eigen::Tensor<double, 3> cpR1(N1, N2, N3);
Eigen::Tensor<double, 3> cpR2(N1, N2, N3);
Eigen::Tensor<double, 3> cpR3(N1, N2, N3);
// �����ٶ�(c_-)
Eigen::Tensor<double, 3> cnL1(N1, N2, N3);
Eigen::Tensor<double, 3> cnL2(N1, N2, N3);
Eigen::Tensor<double, 3> cnL3(N1, N2, N3);
Eigen::Tensor<double, 3> cnR1(N1, N2, N3);
Eigen::Tensor<double, 3> cnR2(N1, N2, N3);
Eigen::Tensor<double, 3> cnR3(N1, N2, N3);
// ţ�ٷ����
Eigen::Tensor<double, 3> ksi(N1, N2, N3);

// useful functions
template<typename T> void print(T info) { std::cout << info << std::endl; }
inline double max(double x, double y) { return x > y ? x : y; }
inline double max(double x, double y, double z) { return max(x, max(y, z)); }
inline double min(double x, double y) { return x < y ? x : y; }
inline double min(double x, double y, double z) { return min(x, min(y, z)); }

inline double MC(double x, double y, double z)
{
	if (abs(x) < abs(y) && abs(x) < abs(z) && y * z > 0)
		return x;
	else if (abs(y) < abs(x) && abs(y) < abs(z) && y * z > 0)
		return y;
	else if (abs(z) < abs(x) && abs(z) < abs(y) && y * z > 0)
		return z;
	else
		return 0;
}

double dot(int i, int j, int k, Eigen::Vector3d vecA, Eigen::Vector3d vecB) {
	return double(vecA.transpose() * metricFuncField(i + NG, j + NG, k + NG).gamma() * vecB);
}

double square(int i, int j, int k, Eigen::Vector3d vec) {
	return dot(i, j, k, vec, vec);
}

double contract(Eigen::Matrix3d A, Eigen::Matrix3d B) {
	double sum = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			sum += A(i, j) * B(i, j);
	return sum;
};

void ghostify(Eigen::Tensor<double, 4> prim) {
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

void interpolate(Eigen::Tensor<double, 4> prim) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int index = 0; index < NPRIM; index++)
				{
					primL1(i, j, k, index) = prim(i + NG, j + NG, k + NG, index) - MC((prim(i + NG + 1, j + NG, k + NG, index) - prim(i + NG - 1, j + NG, k + NG, index)) / (2 * dx1),
						2 * (prim(i + NG + 1, j + NG, k + NG, index) - prim(i + NG, j + NG, k + NG, index)) / (dx1),
						2 * (prim(i + NG, j + NG, k + NG, index) - prim(i + NG - 1, j + NG, k + NG, index)) / (dx1)) * dx1 / 2;
					primR1(i, j, k, index) = prim(i + NG - 1, j + NG, k + NG, index) + MC((prim(i + NG, j + NG, k + NG, index) - prim(i, j + NG, k + NG, index)) / (2 * dx1),
						2 * (prim(i + NG, j + NG, k + NG, index) - prim(i + NG - 1, j + NG, k + NG, index)) / (dx1),
						2 * (prim(i + NG - 1, j + NG, k + NG, index) - prim(i, j + NG, k + NG, index)) / (dx1)) * dx1 / 2;

					primL2(i, j, k, index) = prim(i + NG, j + NG + 1, k + NG, index) - MC((prim(i + NG, j + NG, k + NG, index) - prim(i + NG, j + NG - 1, k + NG, index)) / (2 * dx2),
						2 * (prim(i + NG, j + NG + 1, k + NG, index) - prim(i + NG, j + NG, k + NG, index)) / (dx2),
						2 * (prim(i + NG, j + NG, k + NG, index) - prim(i + NG, j + NG - 1, k + NG, index)) / (dx2)) * dx2 / 2;
					primR2(i, j, k, index) = prim(i + NG, j + NG, k + NG, index) + MC((prim(i + NG, j + NG + 1, k + NG, index) - prim(i + NG, j, k + NG, index)) / (2 * dx2),
						2 * (prim(i + NG, j + NG, k + NG, index) - prim(i + NG, j + NG - 1, k + NG, index)) / (dx2),
						2 * (prim(i + NG, j + NG - 1, k + NG, index) - prim(i + NG, j, k + NG, index)) / (dx2)) * dx2 / 2;

					primL3(i, j, k, index) = prim(i + NG, j + NG, k + NG, index) - MC((prim(i + NG, j + NG, k + NG + 1, index) - prim(i + NG, j + NG, k + NG - 1, index)) / (2 * dx3),
						2 * (prim(i + NG, j + NG, k + NG + 1, index) - prim(i + NG, j + NG, k + NG, index)) / (dx3),
						2 * (prim(i + NG, j + NG, k + NG, index) - prim(i + NG, j + NG, k + NG - 1, index)) / (dx3)) * dx3 / 2;
					primR3(i, j, k, index) = prim(i + NG, j + NG, k + NG - 1, index) + MC((prim(i + NG, j + NG, k + NG, index) - prim(i + NG, j + NG, k, index)) / (2 * dx3),
						2 * (prim(i + NG, j + NG, k + NG, index) - prim(i + NG, j + NG, k + NG - 1, index)) / (dx3),
						2 * (prim(i + NG, j + NG, k + NG - 1, index) - prim(i + NG, j + NG, k, index)) / (dx3)) * dx3 / 2;
				}
}

void prim2con(Eigen::Tensor<double, 4> prim, Eigen::Tensor<double, 4>& con) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim(i, j, k, U1) ,prim(i, j, k, U2) ,prim(i, j, k, U3) };
				Eigen::Vector3d B{ prim(i, j, k, B1) ,prim(i, j, k, B2) ,prim(i, j, k, B3) };
				double Gamma = sqrt(1 + square(i, j, k, u));
				con(i, j, k, 0) = Gamma * prim(i, j, k, RHO);
				con(i, j, k, 1) = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) - prim(i, j, k, UU) + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma) - pow(dot(i, j, k, B, u / Gamma), 2))) - Gamma * prim(i, j, k, RHO);
				con(i, j, k, 2) = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) * prim(i, j, k, U1) + square(i, j, k, B) * prim(i, j, k, U1) - dot(i, j, k, B, u / Gamma) * prim(i, j, k, B1);
				con(i, j, k, 3) = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) * prim(i, j, k, U2) + square(i, j, k, B) * prim(i, j, k, U2) - dot(i, j, k, B, u / Gamma) * prim(i, j, k, B2);
				con(i, j, k, 4) = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) * prim(i, j, k, U3) + square(i, j, k, B) * prim(i, j, k, U3) - dot(i, j, k, B, u / Gamma) * prim(i, j, k, B3);
				con(i, j, k, 5) = prim(i, j, k, B1);
				con(i, j, k, 6) = prim(i, j, k, B2);
				con(i, j, k, 7) = prim(i, j, k, B3);
			}
}

double f(int i, int j, int k, double D, double tau, Eigen::Vector3d S, Eigen::Vector3d B, double x) {
	auto Gamma = 1 / sqrt(1 - square(i, j, k, S + dot(i, j, k, S, B) * B / x) / pow(x + square(i, j, k, B), 2));
	return x - (gam - 1) / gam * (x - Gamma * D) / pow(Gamma, 2) - tau - D + square(i, j, k, B) - 0.5 * (square(i, j, k, B / Gamma) + pow(dot(i, j, k, S, B), 2) / pow(x, 2));
}

double df(int i, int j, int k, double D, double tau, Eigen::Vector3d S, Eigen::Vector3d B, double x) {
	return (f(i, j, k, D, tau, S, B, x + SMALL) - f(i, j, k, D, tau, S, B, x - SMALL)) / (2 * SMALL);
}

void con2prim(Eigen::Tensor<double, 4> con, Eigen::Tensor<double, 4>& prim) {
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
					auto x1 = x0 - f(i, j, k, D, tau, S, B, x0) / df(i, j, k, D, tau, S, B, x0); // ţ�ٵ�����ʽ
					if (abs((x1 - x0) / x0) < tol)
						break;
					x0 = x1;
				}
				if (ksi(i, j, k) <= 0 || isnan(x0))
					continue;
				ksi(i, j, k) = x0;
				auto Gamma = 1 / sqrt(1 - square(i, j, k, S + dot(i, j, k, S, B) * B / ksi(i, j, k)) / pow(ksi(i, j, k) + square(i, j, k, B), 2));
				prim(i + NG, j + NG, k + NG, RHO) = D / Gamma;
				prim(i + NG, j + NG, k + NG, UU) = (gam - 1) / gam * (ksi(i, j, k) - Gamma * D) / pow(Gamma, 2);
				prim(i + NG, j + NG, k + NG, U1) = (S(0) + dot(i, j, k, S, B) * B(0) / ksi(i, j, k)) / (ksi(i, j, k) + square(i, j, k, B)) * Gamma;
				prim(i + NG, j + NG, k + NG, U2) = (S(1) + dot(i, j, k, S, B) * B(1) / ksi(i, j, k)) / (ksi(i, j, k) + square(i, j, k, B)) * Gamma;
				prim(i + NG, j + NG, k + NG, U3) = (S(2) + dot(i, j, k, S, B) * B(2) / ksi(i, j, k)) / (ksi(i, j, k) + square(i, j, k, B)) * Gamma;
				prim(i + NG, j + NG, k + NG, B1) = B(0);
				prim(i + NG, j + NG, k + NG, B2) = B(1);
				prim(i + NG, j + NG, k + NG, B3) = B(2);
			}
		}
	}
}


void prim2flux(Eigen::Tensor<double, 4> prim, Eigen::Tensor<double, 4> con, Eigen::Tensor<double, 4>& flux, short comp) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim(i, j, k, U1) ,prim(i, j, k, U2) ,prim(i, j, k, U3) };
				Eigen::Vector3d B{ prim(i, j, k, B1) ,prim(i, j, k, B2) ,prim(i, j, k, B3) };
				Eigen::Vector3d S{ con(i, j, k, 2) ,con(i, j, k, 3) ,con(i, j, k, 4) };
				double Gamma = sqrt(1 + square(i, j, k, u));
				Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim(i + NG, j + NG, k + NG, UU) + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma)) - pow(dot(i, j, k, B, u / Gamma), 2))) * metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(i, j, k, B, u / Gamma) * u / Gamma * B.transpose();
				flux(i, j, k, 0) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(comp) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp)) * con(i, j, k, 0);
				flux(i, j, k, 1) = metricFuncField(i + NG, j + NG, k + NG).alpha() * (con(i, j, k, 2 + comp) - u(comp) / Gamma * con(i, j, k, 0)) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp) * con(i, j, k, 1);
				flux(i, j, k, 2) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * W * metricFuncField(i + NG, j + NG, k + NG).gamma())(comp, 0) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp) * con(i, j, k, 2);
				flux(i, j, k, 3) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * W * metricFuncField(i + NG, j + NG, k + NG).gamma())(comp, 1) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp) * con(i, j, k, 3);
				flux(i, j, k, 4) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * W * metricFuncField(i + NG, j + NG, k + NG).gamma())(comp, 2) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp) * con(i, j, k, 4);
				flux(i, j, k, 5) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(comp) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp)) * con(i, j, k, 5) - (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(0) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(0)) * con(i, j, k, 5 + comp);
				flux(i, j, k, 6) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(comp) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp)) * con(i, j, k, 6) - (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(1) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(1)) * con(i, j, k, 5 + comp);
				flux(i, j, k, 7) = (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(comp) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(comp)) * con(i, j, k, 7) - (metricFuncField(i + NG, j + NG, k + NG).alpha() * u(2) / Gamma - metricFuncField(i + NG, j + NG, k + NG).betaVec()(2)) * con(i, j, k, 5 + comp);
			}
}

void prim2src(Eigen::Tensor<double, 4> prim, Eigen::Tensor<double, 4> con, Eigen::Tensor<double, 4>& src) {
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
				Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim(i + NG, j + NG, k + NG, UU) + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma)) - pow(dot(i, j, k, B, u / Gamma), 2))) * metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(i, j, k, B, u / Gamma) * u / Gamma * B.transpose();
				src(i, j, k, 0) = 0;
				src(i, j, k, 1) = 0.5 * contract(W, (metricFuncField(i + NG, j + NG, k + NG).betaVec()(0) * metricDiffField(i + NG, j + NG, k + NG, 1).gamma() + metricFuncField(i + NG, j + NG, k + NG).betaVec()(1) * metricDiffField(i + NG, j + NG, k + NG, 2).gamma() + metricFuncField(i + NG, j + NG, k + NG).betaVec()(2) * metricDiffField(i + NG, j + NG, k + NG, 3).gamma()))
					+ contract(W * metricFuncField(i + NG, j + NG, k + NG).gamma(), betaDiff)
					- (metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() * S)(0) * alphaDiffField(i, j, k, 1) - (metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() * S)(1) * alphaDiffField(i, j, k, 2) - (metricFuncField(i + NG, j + NG, k + NG).gamma().inverse() * S)(2) * alphaDiffField(i, j, k, 3);
				src(i, j, k, 2) = 0.5 * metricFuncField(i + NG, j + NG, k + NG).alpha() * contract(W, metricDiffField(i + NG, j + NG, k + NG, 1).gamma()) + dot(i, j, k, S, metricDiffField(i + NG, j + NG, k + NG, 1).betaVec()) - (con(i, j, k, 0) + con(i, j, k, 1)) * alphaDiffField(i, j, k, 1);
				src(i, j, k, 3) = 0.5 * metricFuncField(i + NG, j + NG, k + NG).alpha() * contract(W, metricDiffField(i + NG, j + NG, k + NG, 2).gamma()) + dot(i, j, k, S, metricDiffField(i + NG, j + NG, k + NG, 2).betaVec()) - (con(i, j, k, 0) + con(i, j, k, 1)) * alphaDiffField(i, j, k, 2);
				src(i, j, k, 4) = 0.5 * metricFuncField(i + NG, j + NG, k + NG).alpha() * contract(W, metricDiffField(i + NG, j + NG, k + NG, 3).gamma()) + dot(i, j, k, S, metricDiffField(i + NG, j + NG, k + NG, 3).betaVec()) - (con(i, j, k, 0) + con(i, j, k, 1)) * alphaDiffField(i, j, k, 3);
				src(i, j, k, 5) = 0;
				src(i, j, k, 6) = 0;
				src(i, j, k, 7) = 0;
			}
}

void prim2c(Eigen::Tensor<double, 4> prim, Eigen::Tensor<double, 3>& c, Eigen::Tensor<Metric, 3>& metricFuncHalfField, short sign, short comp) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				if (prim(i, j, k, RHO))
				{
					Eigen::Vector3d u{ prim(i, j, k, U1) ,prim(i, j, k, U2) ,prim(i, j, k, U3) };
					Eigen::Vector3d B{ prim(i, j, k, B1) ,prim(i, j, k, B2) ,prim(i, j, k, B3) };
					double Gamma = sqrt(1 + square(i, j, k, u));
					auto u0 = Gamma / metricFuncField(i + NG, j + NG, k + NG).alpha();
					Eigen::Vector3d ui = { Gamma * (prim(i, j, k, U1) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(0)),
									Gamma * (prim(i, j, k, U2) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(1)),
									Gamma * (prim(i, j, k, U3) - metricFuncField(i + NG, j + NG, k + NG).betaVec()(2)) };
					auto cs_square = gam * prim(i, j, k, UU) / (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU));
					auto cA_square = (square(i, j, k, B) * (1 - square(i, j, k, u / Gamma)) + pow(dot(i, j, k, B, u / Gamma), 2)) / (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU) + square(i, j, k, B) * (1 - square(i, j, k, u / Gamma)) + pow(dot(i, j, k, B, u / Gamma), 2));
					auto vf_square = cA_square + cs_square - cA_square * cs_square;
					auto sigmaf = vf_square ? (1 - vf_square) / vf_square : 0;
					auto metricInv = metricFuncHalfField(i, j, k).m.inverse();
					c(i, j, k) = (metricInv(0, comp) - pow(sigmaf, 2) * u0 * ui(comp)) / (metricInv(0, 0) - pow(sigmaf, 2) * u0 * u0) + sign * sqrt(abs(
						pow((metricInv(0, comp) - pow(sigmaf, 2) * u0 * ui(comp)) / (metricInv(0, 0) - pow(sigmaf, 2) * u0 * u0), 2)
						- (metricInv(comp, comp) - sigmaf * ui(comp) * ui(comp)) / (metricInv(0, 0) - sigmaf * u0 * u0)
					));
				}
}

void calFluxHHL(Eigen::Tensor<double, 3> cpL, Eigen::Tensor<double, 3> cpR,
	Eigen::Tensor<double, 3> cnL, Eigen::Tensor<double, 3> cnR,
	Eigen::Tensor<double, 4> conL, Eigen::Tensor<double, 4> conR,
	Eigen::Tensor<double, 4> fluxL, Eigen::Tensor<double, 4> fluxR,
	Eigen::Tensor<double, 4>& fluxHLL
) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int l = 0; l < 8; l++)
				{
					auto c_max = max(0, cpR(i, j, k), cpL(i, j, k));
					auto c_min = -min(0, cnR(i, j, k), cnL(i, j, k));
					fluxHLL(i, j, k, l) = c_max + c_min ? (c_min * fluxR(i, j, k, l) + c_max * fluxL(i, j, k, l) - c_max * c_min * (conR(i, j, k, l) - conL(i, j, k, l))) / (c_max + c_min) : 0;
				}
}

void calFluxTVDLF(Eigen::Tensor<double, 3> cpL, Eigen::Tensor<double, 3> cpR,
	Eigen::Tensor<double, 3> cnL, Eigen::Tensor<double, 3> cnR,
	Eigen::Tensor<double, 4> conL, Eigen::Tensor<double, 4> conR,
	Eigen::Tensor<double, 4> fluxL, Eigen::Tensor<double, 4> fluxR,
	Eigen::Tensor<double, 4>& fluxTVDLF
) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int l = 0; l < 8; l++)
				{
					auto c_max = max(0, cpR(i, j, k), cpL(i, j, k));
					auto c_min = -min(0, cnR(i, j, k), cnL(i, j, k));
					auto c = max(c_max, c_min);
					fluxTVDLF(i, j, k, l) = 0.5 * (fluxR(i, j, k, l) + fluxL(i, j, k, l)) - 0.5 * c * (conR(i, j, k, l) - conL(i, j, k, l));
				}
}

// ��������
void basicCalc(Eigen::Tensor<double, 4> prim, Eigen::Tensor<double, 4>& con, Eigen::Tensor<double, 4>& flux, Eigen::Tensor<double, 4>& src, Eigen::Tensor<double, 3>& cp, Eigen::Tensor<double, 3>& cn, Eigen::Tensor<Metric, 3>& metricFuncHalfField, short comp) {
	prim2con(prim, con);
	prim2flux(prim, con, flux, comp);
	prim2src(prim, con, src);
	prim2c(prim, cp, metricFuncHalfField, 1, comp);
	prim2c(prim, cn, metricFuncHalfField, -1, comp);
}

void fluxCalc(Eigen::Tensor<double, 3> cpL, Eigen::Tensor<double, 3> cpR, Eigen::Tensor<double, 3> cnL, Eigen::Tensor<double, 3> cnR, Eigen::Tensor<double, 4> conL, Eigen::Tensor<double, 4> conR, Eigen::Tensor<double, 4> fluxL, Eigen::Tensor<double, 4> fluxR, Eigen::Tensor<double, 4>& fluxHLL, Eigen::Tensor<double, 4>& fluxTVDLF) {
	calFluxHHL(cpL, cpR, cnL, cnR, conL, conR, fluxL, fluxR, fluxHLL);
	calFluxTVDLF(cpL, cpR, cnL, cnR, conL, conR, fluxL, fluxR, fluxTVDLF);
}

void check(Eigen::Tensor<double, 4> arr) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int l = 0; l < 8; l++)
					std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tValue: " << arr(i, j, k, l) << std::endl;
	return;
}

void check(Eigen::Tensor<double, 4> arr, int index) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tValue: " << arr(i, j, k, index) << std::endl;
	return;
}

void check(Eigen::Tensor<double, 3> arr) {
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tValue: " << arr(i, j, k) << std::endl;
	return;
}