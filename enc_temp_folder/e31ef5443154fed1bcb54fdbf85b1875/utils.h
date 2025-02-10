#pragma once
#include <omp.h>
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "Metric.h"
MetricComponent ZERO_COMPONENT = [](double x, double y, double z) { return 0; };
constexpr auto a = (0.9375);
constexpr auto h = (0.);
constexpr auto SMALL = (1.e-16);
constexpr auto theta = 0.;								// HHL����TVDLF����ϲ���
constexpr auto NDIM = (4);
constexpr auto N1 = (128);
constexpr auto N2 = (128);
constexpr auto N3 = (4);
constexpr auto NG = (2);
constexpr auto PI = (3.14159265358979323846);
constexpr auto X1min = (0.19325057145871735);
constexpr auto X1max = (6.824046010856292);
constexpr auto X2min = (SMALL);
constexpr auto X2max = (2 * PI - SMALL);
constexpr auto X3min = (SMALL);
constexpr auto X3max = (2. * PI - SMALL);
constexpr auto R0 = (0.);
constexpr auto isX1periodical = false;
constexpr auto isX2periodical = false;
constexpr auto isX3periodical = false;
constexpr auto cour = 0.7;

//FM_torus disk parameter
constexpr auto rin = (6.);
constexpr auto rmax = (10.);
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
double dx1 = (X1max - X1min) / (N1);
double dx2 = (X2max - X2min) / (N2);
double dx3 = (X3max - X3min) / (N3);
double dx1Ghost = (X1max - X1min) / (N1 + 2 * NG);
double dx2Ghost = (X2max - X2min) / (N2 + 2 * NG);
double dx3Ghost = (X3max - X3min) / (N3 + 2 * NG);
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
Eigen::Tensor<Metric, 4> metricDiffField(N1, N2, N3, 4);											// �ȹ浼����
double alphaDiffField[N1][N2][N3][4];																// alpha������
Eigen::Tensor<Metric, 3> metricFuncHalfField1(N1, N2, N3);											// ������ʱ��Ҫ�İ벽���ȹ泡(0,2)��
Eigen::Tensor<Metric, 3> metricFuncHalfField2(N1, N2, N3);											// ������ʱ��Ҫ�İ벽���ȹ泡(0,2)��
Eigen::Tensor<Metric, 3> metricFuncHalfField3(N1, N2, N3);											// ������ʱ��Ҫ�İ벽���ȹ泡(0,2)��

// ��Ҫ������Ӧ��ͳGRMHD�����е�P(������)
double prim[N1][N2][N3][NPRIM];
double primGhost[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NPRIM];
double primL1[N1][N2][N3][NPRIM];
double primL2[N1][N2][N3][NPRIM];
double primL3[N1][N2][N3][NPRIM];
double primR1[N1][N2][N3][NPRIM];
double primR2[N1][N2][N3][NPRIM];
double primR3[N1][N2][N3][NPRIM];
// �غ�������Ӧ��ͳGRMHD�����е�U(������)
double con[N1][N2][N3][8];
double conHalf[N1][N2][N3][8];
double conL1[N1][N2][N3][8];
double conL2[N1][N2][N3][8];
double conL3[N1][N2][N3][8];
double conR1[N1][N2][N3][8];
double conR2[N1][N2][N3][8];
double conR3[N1][N2][N3][8];
// ��(flux)
double fluxL1[N1][N2][N3][8];
double fluxL2[N1][N2][N3][8];
double fluxL3[N1][N2][N3][8];
double fluxR1[N1][N2][N3][8];
double fluxR2[N1][N2][N3][8];
double fluxR3[N1][N2][N3][8];
// HHL��
double fluxHLL1[N1][N2][N3][8];
double fluxHLL2[N1][N2][N3][8];
double fluxHLL3[N1][N2][N3][8];
// TVDLF��
double fluxTVDLF1[N1][N2][N3][8];
double fluxTVDLF2[N1][N2][N3][8];
double fluxTVDLF3[N1][N2][N3][8];
// �����
double fluxLLF1[N1][N2][N3][8];
double fluxLLF2[N1][N2][N3][8];
double fluxLLF3[N1][N2][N3][8];
// �⻬�����
double fluxSmoothLLF1[N1][N2][N3][8];
double fluxSmoothLLF2[N1][N2][N3][8];
double fluxSmoothLLF3[N1][N2][N3][8];
// Դ(source)
double src[N1][N2][N3][8];
double srcL1[N1][N2][N3][8];
double srcL2[N1][N2][N3][8];
double srcL3[N1][N2][N3][8];
double srcR1[N1][N2][N3][8];
double srcR2[N1][N2][N3][8];
double srcR3[N1][N2][N3][8];
// �����ٶ�(c_+)
double cpL1[N1][N2][N3];
double cpL2[N1][N2][N3];
double cpL3[N1][N2][N3];
double cpR1[N1][N2][N3];
double cpR2[N1][N2][N3];
double cpR3[N1][N2][N3];
// �����ٶ�(c_-)
double cnL1[N1][N2][N3];
double cnL2[N1][N2][N3];
double cnL3[N1][N2][N3];
double cnR1[N1][N2][N3];
double cnR2[N1][N2][N3];
double cnR3[N1][N2][N3];
// ţ�ٷ����
double ksi[N1][N2][N3];

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

void ghostify() {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				primGhost[NG + i][NG + j][NG + k][RHO] = prim[i][j][k][RHO];
				primGhost[NG + i][NG + j][NG + k][UU] = prim[i][j][k][UU];
				primGhost[NG + i][NG + j][NG + k][U1] = prim[i][j][k][U1];
				primGhost[NG + i][NG + j][NG + k][U2] = prim[i][j][k][U2];
				primGhost[NG + i][NG + j][NG + k][U3] = prim[i][j][k][U3];
				primGhost[NG + i][NG + j][NG + k][B1] = prim[i][j][k][B1];
				primGhost[NG + i][NG + j][NG + k][B2] = prim[i][j][k][B2];
				primGhost[NG + i][NG + j][NG + k][B3] = prim[i][j][k][B3];
			}

#pragma omp parallel for
	for (int i = NG - 1; i >= 0; i--)
	{
		for (int j = NG - 1; j >= 0; j--)
		{
			for (int k = NG - 1; k >= 0; k--)
			{
				primGhost[i][j][k][RHO] = primGhost[i+1][j+1][k+1][RHO] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][UU] = primGhost[i+1][j+1][k+1][UU] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][B1] = primGhost[i+1][j+1][k+1][B1] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

				primGhost[i][j][k][U2] = primGhost[i+1][j+1][k+1][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][U3] = primGhost[i+1][j+1][k+1][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B2] = primGhost[i+1][j+1][k+1][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B3] = primGhost[i+1][j+1][k+1][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i][j][k][U1] = primGhost[i+1][j+1][k+1][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}

			for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
			{
				primGhost[i+1][j+1][k+1][RHO] = primGhost[i][j][k][RHO] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][UU] = primGhost[i][j][k][UU] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][B1] = primGhost[i][j][k][B1] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

				primGhost[i+1][j+1][k+1][U2] = primGhost[i][j][k][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][U3] = primGhost[i][j][k][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B2] = primGhost[i][j][k][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B3] = primGhost[i][j][k][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i+1][j+1][k+1][U1] = primGhost[i][j][k][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}
		}

		for (int j = NG + N2; j < 2 * NG + N2 - 1; j++)
		{
			for (int k = NG - 1; k >= 0; k--)
			{
				primGhost[i][j][k][RHO] = primGhost[i+1][j+1][k+1][RHO] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][UU] = primGhost[i+1][j+1][k+1][UU] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][B1] = primGhost[i+1][j+1][k+1][B1] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

				primGhost[i][j][k][U2] = primGhost[i+1][j+1][k+1][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][U3] = primGhost[i+1][j+1][k+1][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B2] = primGhost[i+1][j+1][k+1][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B3] = primGhost[i+1][j+1][k+1][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i][j][k][U1] = primGhost[i+1][j+1][k+1][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}

			for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
			{
				primGhost[i+1][j+1][k+1][RHO] = primGhost[i][j][k][RHO] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][UU] = primGhost[i][j][k][UU] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][B1] = primGhost[i][j][k][B1] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

				primGhost[i+1][j+1][k+1][U2] = primGhost[i][j][k][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][U3] = primGhost[i][j][k][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B2] = primGhost[i][j][k][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B3] = primGhost[i][j][k][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i+1][j+1][k+1][U1] = primGhost[i][j][k][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}
		}
	}

#pragma omp parallel for
	for (int i = NG + N1; i < 2 * NG + N1 - 1; i++)
	{
		for (int j = NG - 1; j >= 0; j--)
		{
			for (int k = NG - 1; k >= 0; k--)
			{
				primGhost[i][j][k][RHO] = primGhost[i+1][j+1][k+1][RHO] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][UU] = primGhost[i+1][j+1][k+1][UU] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][B1] = primGhost[i+1][j+1][k+1][B1] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

				primGhost[i][j][k][U2] = primGhost[i+1][j+1][k+1][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][U3] = primGhost[i+1][j+1][k+1][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B2] = primGhost[i+1][j+1][k+1][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B3] = primGhost[i+1][j+1][k+1][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i][j][k][U1] = primGhost[i+1][j+1][k+1][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}

			for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
			{
				primGhost[i+1][j+1][k+1][RHO] = primGhost[i][j][k][RHO] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][UU] = primGhost[i][j][k][UU] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][B1] = primGhost[i][j][k][B1] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

				primGhost[i+1][j+1][k+1][U2] = primGhost[i][j][k][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][U3] = primGhost[i][j][k][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B2] = primGhost[i][j][k][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B3] = primGhost[i][j][k][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i+1][j+1][k+1][U1] = primGhost[i][j][k][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}
		}
		for (int j = NG + N2; j < 2 * NG + N2 - 1; j++)
		{
			for (int k = NG - 1; k >= 0; k--)
			{
				primGhost[i][j][k][RHO] = primGhost[i+1][j+1][k+1][RHO] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][UU] = primGhost[i+1][j+1][k+1][UU] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][B1] = primGhost[i+1][j+1][k+1][B1] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

				primGhost[i][j][k][U2] = primGhost[i+1][j+1][k+1][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][U3] = primGhost[i+1][j+1][k+1][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B2] = primGhost[i+1][j+1][k+1][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i][j][k][B3] = primGhost[i+1][j+1][k+1][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i][j][k][U1] = primGhost[i+1][j+1][k+1][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}

			for (int k = NG + N3; k < 2 * NG + N3 - 1; k++)
			{
				primGhost[i+1][j+1][k+1][RHO] = primGhost[i][j][k][RHO] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][UU] = primGhost[i][j][k][UU] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());
				primGhost[i+1][j+1][k+1][B1] = primGhost[i][j][k][B1] * sqrt(-metricFuncField(i, j, k).m.determinant()) / sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant());

				primGhost[i+1][j+1][k+1][U2] = primGhost[i][j][k][U2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][U3] = primGhost[i][j][k][U3] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B2] = primGhost[i][j][k][B2] * (1 - dx1 / (X1min + (i + 1) * dx1));
				primGhost[i+1][j+1][k+1][B3] = primGhost[i][j][k][B3] * (1 - dx1 / (X1min + (i + 1) * dx1));

				primGhost[i+1][j+1][k+1][U1] = primGhost[i][j][k][U1] * (1 + dx1 / (X1min + (i + 1) * dx1));
			}
		}
	}
}

void interpolate() {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int index = 0; index < NPRIM; index++)
				{
					primL1[i][j][k][index] = primGhost[i+NG][j+NG][k+NG][index] - MC((primGhost[i + NG + 1][j + NG][k + NG][index] - primGhost[i + NG - 1][j + NG][k + NG][index]) / (2 * dx1),
						2 * (primGhost[i + NG + 1][j + NG][k + NG][index] - primGhost[i+NG][j+NG][k+NG][index]) / (dx1),
						2 * (primGhost[i+NG][j+NG][k+NG][index] - primGhost[i + NG - 1][j + NG][k + NG][index]) / (dx1)) * dx1 / 2;
					primR1[i][j][k][index] = primGhost[i + NG - 1][j + NG][k + NG][index] + MC((primGhost[i+NG][j+NG][k+NG][index] - primGhost[i][j + NG][k + NG][index]) / (2 * dx1),
						2 * (primGhost[i+NG][j+NG][k+NG][index] - primGhost[i + NG - 1][j + NG][k + NG][index]) / (dx1),
						2 * (primGhost[i + NG - 1][j + NG][k + NG][index] - primGhost[i][j + NG][k + NG][index]) / (dx1)) * dx1 / 2;

					primL2[i][j][k][index] = primGhost[i + NG][j + NG + 1][k + NG][index] - MC((primGhost[i+NG][j+NG][k+NG][index] - primGhost[i + NG][j + NG - 1][k + NG][index]) / (2 * dx2),
						2 * (primGhost[i + NG][j + NG + 1][k + NG][index] - primGhost[i+NG][j+NG][k+NG][index]) / (dx2),
						2 * (primGhost[i+NG][j+NG][k+NG][index] - primGhost[i + NG][j + NG - 1][k + NG][index]) / (dx2)) * dx2 / 2;
					primR2[i][j][k][index] = primGhost[i+NG][j+NG][k+NG][index] + MC((primGhost[i + NG][j + NG + 1][k + NG][index] - primGhost[i + NG][j][k + NG][index]) / (2 * dx2),
						2 * (primGhost[i+NG][j+NG][k+NG][index] - primGhost[i + NG][j + NG - 1][k + NG][index]) / (dx2),
						2 * (primGhost[i + NG][j + NG - 1][k + NG][index] - primGhost[i + NG][j][k + NG][index]) / (dx2)) * dx2 / 2;

					primL3[i][j][k][index] = primGhost[i + NG][j + NG][k + NG][index] - MC((primGhost[i + NG][j + NG][k + NG + 1][index] - primGhost[i + NG][j + NG][k + NG - 1][index]) / (2 * dx3),
						2 * (primGhost[i + NG][j + NG][k + NG + 1][index] - primGhost[i + NG][j + NG][k + NG][index]) / (dx3),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG][k + NG - 1][index]) / (dx3)) * dx3 / 2;
					primR3[i][j][k][index] = primGhost[i+NG][j+NG][k+NG-1][index] + MC((primGhost[i+NG][j+NG][k+NG][index] - primGhost[i + NG][j + NG][k][index]) / (2 * dx3),
						2 * (primGhost[i+NG][j+NG][k+NG][index] - primGhost[i+NG][j+NG][k+NG-1][index]) / (dx3),
						2 * (primGhost[i+NG][j+NG][k+NG-1][index] - primGhost[i + NG][j + NG][k][index]) / (dx3)) * dx3 / 2;
				}
}

void prim2con(double prim[N1][N2][N3][NPRIM], double con[N1][N2][N3][8]) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim[i][j][k][U1] ,prim[i][j][k][U2] ,prim[i][j][k][U3] };
				Eigen::Vector3d B{ prim[i][j][k][B1] ,prim[i][j][k][B2] ,prim[i][j][k][B3] };
				double Gamma = sqrt(1 + square(i, j, k, u));
				con[i][j][k][0] = Gamma * prim[i][j][k][RHO];
				con[i][j][k][1] = (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]) * pow(Gamma, 2) - prim[i][j][k][UU] + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma) - pow(dot(i, j, k, B, u / Gamma), 2))) - Gamma * prim[i][j][k][RHO];
				con[i][j][k][2] = (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]) * pow(Gamma, 2) * prim[i][j][k][U1] + square(i, j, k, B) * prim[i][j][k][U1] - dot(i, j, k, B, u / Gamma) * prim[i][j][k][B1];
				con[i][j][k][3] = (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]) * pow(Gamma, 2) * prim[i][j][k][U2] + square(i, j, k, B) * prim[i][j][k][U2] - dot(i, j, k, B, u / Gamma) * prim[i][j][k][B2];
				con[i][j][k][4] = (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]) * pow(Gamma, 2) * prim[i][j][k][U3] + square(i, j, k, B) * prim[i][j][k][U3] - dot(i, j, k, B, u / Gamma) * prim[i][j][k][B3];
				con[i][j][k][5] = prim[i][j][k][B1];
				con[i][j][k][6] = prim[i][j][k][B2];
				con[i][j][k][7] = prim[i][j][k][B3];
			}
}

double f(int i, int j, int k, double D, double tau, Eigen::Vector3d S, Eigen::Vector3d B, double x) {
	auto Gamma = 1 / sqrt(1 - square(i, j, k, S + dot(i, j, k, S, B) * B / x) / pow(x + square(i, j, k, B), 2));
	return x - (gam - 1) / gam * (x - Gamma * D) / pow(Gamma, 2) - tau - D + square(i, j, k, B) - 0.5 * (square(i, j, k, B / Gamma) + pow(dot(i, j, k, S, B), 2) / pow(x, 2));
}

double df(int i, int j, int k, double D, double tau, Eigen::Vector3d S, Eigen::Vector3d B, double x) {
	return (f(i, j, k, D, tau, S, B, x + SMALL) - f(i, j, k, D, tau, S, B, x - SMALL)) / (2 * SMALL);
}

void con2prim(double con[N1][N2][N3][8]) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d S{ con[i][j][k][2] ,con[i][j][k][3] ,con[i][j][k][4] };
				Eigen::Vector3d B{ con[i][j][k][5] ,con[i][j][k][6] ,con[i][j][k][7] };
				auto D = con[i][j][k][0];
				auto tau = con[i][j][k][1];
				auto x0 = ksi[i][j][k];
				for (int iter = 0; iter < max_iter; iter++)
				{
					auto x1 = x0 - f(i, j, k, D, tau, S, B, x0) / df(i, j, k, D, tau, S, B, x0);
					if (abs((x1 - x0) / x0) < tol)
						break;
					x0 = x1;
				}
				if (x0 <= SMALL || isnan(x0) || isinf(x0))
					continue;
				ksi[i][j][k] = x0;
				auto Gamma = 1 / sqrt(1 - square(i, j, k, S + dot(i, j, k, S, B) * B / ksi[i][j][k]) / pow(ksi[i][j][k] + square(i, j, k, B), 2));
				prim[i][j][k][RHO] = D / Gamma;
				prim[i][j][k][UU] = (gam - 1) / gam * (ksi[i][j][k] - Gamma * D) / pow(Gamma, 2);
				prim[i][j][k][U1] = (S(0) + dot(i, j, k, S, B) * B(0) / ksi[i][j][k]) / (ksi[i][j][k] + square(i, j, k, B)) * Gamma;
				prim[i][j][k][U2] = (S(1) + dot(i, j, k, S, B) * B(1) / ksi[i][j][k]) / (ksi[i][j][k] + square(i, j, k, B)) * Gamma;
				prim[i][j][k][U3] = (S(2) + dot(i, j, k, S, B) * B(2) / ksi[i][j][k]) / (ksi[i][j][k] + square(i, j, k, B)) * Gamma;
				prim[i][j][k][B1] = B(0);
				prim[i][j][k][B2] = B(1);
				prim[i][j][k][B3] = B(2);
			}
}


void prim2flux(double prim[N1][N2][N3][NPRIM], double con[N1][N2][N3][8], double flux[N1][N2][N3][8], short comp) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim[i][j][k][U1] ,prim[i][j][k][U2] ,prim[i][j][k][U3] };
				Eigen::Vector3d B{ prim[i][j][k][B1] ,prim[i][j][k][B2] ,prim[i][j][k][B3] };
				Eigen::Vector3d S{ con[i][j][k][2] ,con[i][j][k][3] ,con[i][j][k][4] };
				double Gamma = sqrt(1 + square(i, j, k, u));
				Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim[i][j][k][UU] + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma)) - pow(dot(i, j, k, B, u / Gamma), 2))) * metricFuncField(i, j, k).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(i, j, k, B, u / Gamma) * u / Gamma * B.transpose();
				flux[i][j][k][0] = (metricFuncField(i, j, k).alpha() * u(comp) / Gamma - metricFuncField(i, j, k).betaVec()(comp)) * con[i][j][k][0];
				flux[i][j][k][1] = metricFuncField(i, j, k).alpha() * (con[i][j][k][2 + comp] - u(comp) / Gamma * con[i][j][k][0]) - metricFuncField(i, j, k).betaVec()(comp) * con[i][j][k][1];
				flux[i][j][k][2] = (metricFuncField(i, j, k).alpha() * W * metricFuncField(i, j, k).gamma())(comp, 0) - metricFuncField(i, j, k).betaVec()(comp) * con[i][j][k][2];
				flux[i][j][k][3] = (metricFuncField(i, j, k).alpha() * W * metricFuncField(i, j, k).gamma())(comp, 1) - metricFuncField(i, j, k).betaVec()(comp) * con[i][j][k][3];
				flux[i][j][k][4] = (metricFuncField(i, j, k).alpha() * W * metricFuncField(i, j, k).gamma())(comp, 2) - metricFuncField(i, j, k).betaVec()(comp) * con[i][j][k][4];
				flux[i][j][k][5] = (metricFuncField(i, j, k).alpha() * u(comp) / Gamma - metricFuncField(i, j, k).betaVec()(comp)) * con[i][j][k][5] - (metricFuncField(i, j, k).alpha() * u(0) / Gamma - metricFuncField(i, j, k).betaVec()(0)) * con[i][j][k][5 + comp];
				flux[i][j][k][6] = (metricFuncField(i, j, k).alpha() * u(comp) / Gamma - metricFuncField(i, j, k).betaVec()(comp)) * con[i][j][k][6] - (metricFuncField(i, j, k).alpha() * u(1) / Gamma - metricFuncField(i, j, k).betaVec()(1)) * con[i][j][k][5 + comp];
				flux[i][j][k][7] = (metricFuncField(i, j, k).alpha() * u(comp) / Gamma - metricFuncField(i, j, k).betaVec()(comp)) * con[i][j][k][7] - (metricFuncField(i, j, k).alpha() * u(2) / Gamma - metricFuncField(i, j, k).betaVec()(2)) * con[i][j][k][5 + comp];
			}
}

void prim2src(double prim[N1][N2][N3][NPRIM], double con[N1][N2][N3][8], double src[N1][N2][N3][8]) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim[i][j][k][U1] ,prim[i][j][k][U2] ,prim[i][j][k][U3] };
				Eigen::Vector3d B{ prim[i][j][k][B1] ,prim[i][j][k][B2] ,prim[i][j][k][B3] };
				Eigen::Vector3d S{ con[i][j][k][2] ,con[i][j][k][3] ,con[i][j][k][4] };
				Eigen::Matrix3d betaDiff;
				betaDiff << metricDiffField(i, j, k, 1).betaVec()(0), metricDiffField(i, j, k, 2).betaVec()(0), metricDiffField(i, j, k, 3).betaVec()(0),
					metricDiffField(i, j, k, 1).betaVec()(1), metricDiffField(i, j, k, 2).betaVec()(1), metricDiffField(i, j, k, 3).betaVec()(1),
					metricDiffField(i, j, k, 1).betaVec()(2), metricDiffField(i, j, k, 2).betaVec()(2), metricDiffField(i, j, k, 3).betaVec()(2);
				double Gamma = sqrt(1 + square(i, j, k, u));
				Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim[i][j][k][UU] + 0.5 * (square(i, j, k, B) * (1 + square(i, j, k, u / Gamma)) - pow(dot(i, j, k, B, u / Gamma), 2))) * metricFuncField(i, j, k).gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - dot(i, j, k, B, u / Gamma) * u / Gamma * B.transpose();
				src[i][j][k][0] = 0;
				src[i][j][k][1] = 0.5 * contract(W, (metricFuncField(i, j, k).betaVec()(0) * metricDiffField(i, j, k, 1).gamma() + metricFuncField(i, j, k).betaVec()(1) * metricDiffField(i, j, k, 2).gamma() + metricFuncField(i, j, k).betaVec()(2) * metricDiffField(i, j, k, 3).gamma()))
					+ contract(W * metricFuncField(i, j, k).gamma(), betaDiff)
					- (metricFuncField(i, j, k).gamma().inverse() * S)(0) * alphaDiffField[i][j][k][1] - (metricFuncField(i, j, k).gamma().inverse() * S)(1) * alphaDiffField[i][j][k][2] - (metricFuncField(i, j, k).gamma().inverse() * S)(2) * alphaDiffField[i][j][k][3];
				src[i][j][k][2] = 0.5 * metricFuncField(i, j, k).alpha() * contract(W, metricDiffField(i, j, k, 1).gamma()) + dot(i, j, k, S, metricDiffField(i, j, k, 1).betaVec()) - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffField[i][j][k][1];
				src[i][j][k][3] = 0.5 * metricFuncField(i, j, k).alpha() * contract(W, metricDiffField(i, j, k, 2).gamma()) + dot(i, j, k, S, metricDiffField(i, j, k, 2).betaVec()) - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffField[i][j][k][2];
				src[i][j][k][4] = 0.5 * metricFuncField(i, j, k).alpha() * contract(W, metricDiffField(i, j, k, 3).gamma()) + dot(i, j, k, S, metricDiffField(i, j, k, 3).betaVec()) - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffField[i][j][k][3];
				src[i][j][k][5] = 0;
				src[i][j][k][6] = 0;
				src[i][j][k][7] = 0;
			}
}

void prim2c(double prim[N1][N2][N3][NPRIM], double c[N1][N2][N3], Eigen::Tensor<Metric, 3>& metricFuncHalfField, short sign, short comp) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				if (prim[i][j][k][RHO])
				{
					Eigen::Vector3d u{ prim[i][j][k][U1] ,prim[i][j][k][U2] ,prim[i][j][k][U3] };
					Eigen::Vector3d B{ prim[i][j][k][B1] ,prim[i][j][k][B2] ,prim[i][j][k][B3] };
					double Gamma = sqrt(1 + square(i, j, k, u));
					auto u0 = Gamma / metricFuncField(i, j, k).alpha();
					Eigen::Vector3d ui = { Gamma * (prim[i][j][k][U1] - metricFuncField(i, j, k).betaVec()(0)),
									Gamma * (prim[i][j][k][U2] - metricFuncField(i, j, k).betaVec()(1)),
									Gamma * (prim[i][j][k][U3] - metricFuncField(i, j, k).betaVec()(2)) };
					auto cs_square = gam * prim[i][j][k][UU] / (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]);
					auto cA_square = (square(i, j, k, B) * (1 - square(i, j, k, u / Gamma)) + pow(dot(i, j, k, B, u / Gamma), 2)) / (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU] + square(i, j, k, B) * (1 - square(i, j, k, u / Gamma)) + pow(dot(i, j, k, B, u / Gamma), 2));
					auto vf_square = cA_square + cs_square - cA_square * cs_square;
					auto sigmaf = vf_square ? (1 - vf_square) / vf_square : 0;
					auto metricInv = metricFuncHalfField(i, j, k).m.inverse();
					c[i][j][k] = (metricInv(0, comp) - pow(sigmaf, 2) * u0 * ui(comp)) / (metricInv(0, 0) - pow(sigmaf, 2) * u0 * u0) + sign * sqrt(abs(
						pow((metricInv(0, comp) - pow(sigmaf, 2) * u0 * ui(comp)) / (metricInv(0, 0) - pow(sigmaf, 2) * u0 * u0), 2)
						- (metricInv(comp, comp) - sigmaf * ui(comp) * ui(comp)) / (metricInv(0, 0) - sigmaf * u0 * u0)
					));
				}
}

void calFluxHHL(double cpL[N1][N2][N3], double cpR[N1][N2][N3],
	double cnL[N1][N2][N3], double cnR[N1][N2][N3],
	double conL[N1][N2][N3][8], double conR[N1][N2][N3][8],
	double fluxL[N1][N2][N3][8], double fluxR[N1][N2][N3][8],
	double fluxHLL[N1][N2][N3][8]
) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int l = 0; l < 8; l++)
				{
					auto c_max = max(0, cpR[i][j][k], cpL[i][j][k]);
					auto c_min = -min(0, cnR[i][j][k], cnL[i][j][k]);
					fluxHLL[i][j][k][l] = c_max + c_min ? (c_min * fluxR[i][j][k][l] + c_max * fluxL[i][j][k][l] - c_max * c_min * (conR[i][j][k][l] - conL[i][j][k][l])) / (c_max + c_min) : 0;
				}
}

void calFluxTVDLF(double cpL[N1][N2][N3], double cpR[N1][N2][N3],
	double cnL[N1][N2][N3], double cnR[N1][N2][N3],
	double conL[N1][N2][N3][8], double conR[N1][N2][N3][8],
	double fluxL[N1][N2][N3][8], double fluxR[N1][N2][N3][8],
	double fluxTVDLF[N1][N2][N3][8]
) {
#pragma omp parallel for
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int l = 0; l < 8; l++)
				{
					auto c_max = max(0, cpR[i][j][k], cpL[i][j][k]);
					auto c_min = -min(0, cnR[i][j][k], cnL[i][j][k]);
					auto c = max(c_max, c_min);
					fluxTVDLF[i][j][k][l] = 0.5 * (fluxR[i][j][k][l] + fluxL[i][j][k][l]) - 0.5 * c * (conR[i][j][k][l] - conL[i][j][k][l]);
				}
}

void basicCalc(double prim[N1][N2][N3][NPRIM], double con[N1][N2][N3][8], double flux[N1][N2][N3][8], double src[N1][N2][N3][8], double cp[N1][N2][N3], double cn[N1][N2][N3], Eigen::Tensor<Metric, 3>& metricFuncHalfField, short comp) {
	prim2con(prim, con);
	prim2flux(prim, con, flux, comp);
	prim2src(prim, con, src);
	prim2c(prim, cp, metricFuncHalfField, 1, comp);
	prim2c(prim, cn, metricFuncHalfField, -1, comp);
}

void fluxCalc(double cpL[N1][N2][N3], double cpR[N1][N2][N3], double cnL[N1][N2][N3], double cnR[N1][N2][N3], double conL[N1][N2][N3][8], double conR[N1][N2][N3][8], double fluxL[N1][N2][N3][8], double fluxR[N1][N2][N3][8], double fluxHLL[N1][N2][N3][8], double fluxTVDLF[N1][N2][N3][8]) {
	calFluxHHL(cpL, cpR, cnL, cnR, conL, conR, fluxL, fluxR, fluxHLL);
	calFluxTVDLF(cpL, cpR, cnL, cnR, conL, conR, fluxL, fluxR, fluxTVDLF);
}


