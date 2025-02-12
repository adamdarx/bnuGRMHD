#pragma once
#include <omp.h>
#include <iostream>
#include <AMReX_MultiFab.H> //For the method most common at time of writing
#include <AMReX_MFParallelFor.H> //For the second newer method
#include <unsupported/Eigen/CXX11/Tensor>
#include "Metric.h"

const MetricComponent ZERO_COMPONENT = [](double x, double y, double z) { return 0; };
constexpr auto a = (0.9375);
constexpr auto h = (0.);
constexpr auto SMALL = (1.e-16);
constexpr auto theta = 0.;
constexpr auto NDIM = (4);
constexpr auto N1 = (16);
constexpr auto N2 = (16);
constexpr auto N3 = (4);
constexpr auto NG = (2);
constexpr auto max_grid_size = (16);
constexpr auto PI = (3.14159265358979323846);
constexpr auto X1min = (0.19325057145871735);
constexpr auto X1max = (7.824046010856292);
constexpr auto X2min = (SMALL);
constexpr auto X2max = (PI - SMALL);
constexpr auto X3min = (SMALL);
constexpr auto X3max = (2. * PI - SMALL);
constexpr auto R0 = (0.);
constexpr auto isX1periodical = false;
constexpr auto isX2periodical = false;
constexpr auto isX3periodical = false;
constexpr auto cour = 1e-2;
constexpr auto LEFT = 0;
constexpr auto RIGHT = 1;
constexpr auto NEG = 0;
constexpr auto POS = 1;

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
double dx1 = (X1max - X1min) / (N1);
double dx2 = (X2max - X2min) / (N2);
double dx3 = (X3max - X3min) / (N3);
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
Eigen::Tensor<MetricComponent, 2> metricFunc(4, 4);												
Eigen::Tensor<MetricComponent, 3> metricDiff(4, 4, 4);			

Eigen::Tensor<Metric, 3> metricFuncField(N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG);	
Eigen::Tensor<Metric, 3> metricFuncHalfField1(N1 + 1, N2 + 1, N3 + 1);
Eigen::Tensor<Metric, 3> metricFuncHalfField2(N1 + 1, N2 + 1, N3 + 1);
Eigen::Tensor<Metric, 3> metricFuncHalfField3(N1 + 1, N2 + 1, N3 + 1);
Eigen::Tensor<Metric, 3> metricFuncHalfField[3] = { metricFuncHalfField1, metricFuncHalfField2, metricFuncHalfField3 };

Eigen::Tensor<Metric, 4> metricDiffField(N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, 4);
Eigen::Tensor<Metric, 4> metricDiffHalfField1(N1 + 1, N2 + 1, N3 + 1, 4);
Eigen::Tensor<Metric, 4> metricDiffHalfField2(N1 + 1, N2 + 1, N3 + 1, 4);
Eigen::Tensor<Metric, 4> metricDiffHalfField3(N1 + 1, N2 + 1, N3 + 1, 4);
Eigen::Tensor<Metric, 4> metricDiffHalfField[3] = { metricDiffHalfField1, metricDiffHalfField2, metricDiffHalfField3 };

amrex::MultiFab alphaDiffField;
amrex::MultiFab alphaDiffHalfField[3];
amrex::MultiFab prim;
amrex::MultiFab con;
amrex::MultiFab src;
amrex::MultiFab ksi;
amrex::MultiFab primGhost;
amrex::MultiFab primLR[2][3];				// 2是左右(0为左，1为右)，3是方向
amrex::MultiFab conLR[2][3];
amrex::MultiFab srcLR[2][3];
amrex::MultiFab fluxLR[2][3];
amrex::MultiFab fluxHLL[3];					// 3是方向
amrex::MultiFab fluxTVDLF[3];
amrex::MultiFab fluxLLF[3];
amrex::MultiFab fluxSmoothLLF[3];						// 分别是正负，左右，方向
amrex::MultiFab c[2][2][3];
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

double dot(int i, int j, int k, Eigen::Vector3d vecA, Eigen::Vector3d vecB, Eigen::Tensor<Metric, 3> metric) {
	return double(vecA.transpose() * metric(i, j, k).gamma() * vecB);
}

double square(int i, int j, int k, Eigen::Vector3d vec, Eigen::Tensor<Metric, 3> metric) {
	return dot(i, j, k, vec, vec, metric);
}

double contract(Eigen::Matrix3d A, Eigen::Matrix3d B) {
	double sum = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			sum += A(i, j) * B(i, j);
	return sum;
};

void ghostify() {
	for (MFIter mfi(prim, TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
		Array4<Real> const& a = primGhost[mfi].array();
		Array4<Real const> const& b = prim[mfi].array()
		ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, int l)
		{
			a(i, j, k, l) = b(i, j, k, l);
		});
	}

	// 1. 鬼化X1方向

	for (int i = NG - 1; i >= 0; i--)
	{

		for (int j = NG; j < N2 + NG; j++)
		{

			for (int k = NG; k < N3 + NG; k++)
			{
				primGhost(i, j, k, RHO) = primGhost[i + 1][j + 1][k + 1][RHO] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost(i, j, k, UU) = primGhost[i + 1][j + 1][k + 1][UU] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][B1] = primGhost[i + 1][j + 1][k + 1][B1] * sqrt(-metricFuncField(i + 1, j + 1, k + 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

				primGhost(i, j, k, U2) = primGhost[i + 1][j + 1][k + 1][U2] * (1 - dx1);
				primGhost(i, j, k, U3) = primGhost[i + 1][j + 1][k + 1][U3] * (1 - dx1);
				primGhost[i][j][k][B2] = primGhost[i + 1][j + 1][k + 1][B2] * (1 - dx1);
				primGhost[i][j][k][B3] = primGhost[i + 1][j + 1][k + 1][B3] * (1 - dx1);

				primGhost(i, j, k, U1) = primGhost[i + 1][j + 1][k + 1][U1] * (1 + dx1);
			}
		}
	}


	for (int i = NG + N1; i < 2 * NG + N1; i++)
	{

		for (int j = NG; j < N2 + NG; j++)
		{

			for (int k = NG; k < N3 + NG; k++)
			{
				primGhost(i, j, k, RHO) = primGhost[i - 1][j - 1][k - 1][RHO] * sqrt(-metricFuncField(i - 1, j - 1, k - 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost(i, j, k, UU) = primGhost[i - 1][j - 1][k - 1][UU] * sqrt(-metricFuncField(i - 1, j - 1, k - 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());
				primGhost[i][j][k][B1] = primGhost[i - 1][j - 1][k - 1][B1] * sqrt(-metricFuncField(i - 1, j - 1, k - 1).m.determinant()) / sqrt(-metricFuncField(i, j, k).m.determinant());

				primGhost(i, j, k, U2) = primGhost[i - 1][j - 1][k - 1][U2] * (1 - dx1);
				primGhost(i, j, k, U3) = primGhost[i - 1][j - 1][k - 1][U3] * (1 - dx1);
				primGhost[i][j][k][B2] = primGhost[i - 1][j - 1][k - 1][B2] * (1 - dx1);
				primGhost[i][j][k][B3] = primGhost[i - 1][j - 1][k - 1][B3] * (1 - dx1);

				primGhost(i, j, k, U1) = primGhost[i - 1][j - 1][k - 1][U1] * (1 + dx1);
			}
		}
	}

	// 2. 鬼化X2方向

	for(int i = 0; i < N1 + 2 * NG; i++)
	{

		for (int j = 0; j < NG; j++)
		{

			for (int k = NG; k < N3 + NG; k++)
			{
				// 1) 把上面的格子移动到下面
				primGhost(i, j, k, RHO) = primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][RHO];
				primGhost(i, j, k, UU) = primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][UU];
				primGhost(i, j, k, U1) = primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][U1];
				primGhost(i, j, k, U2) = -primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][U2];
				primGhost(i, j, k, U3) = primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][U3];
				primGhost[i][j][k][B1] = primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][B1];
				primGhost[i][j][k][B2] = -primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][B2];
				primGhost[i][j][k][B3] = primGhost[i][NG + 1 - j][(k + N3 / 2) % N3][B3];
				// 2) 把下面的格子移动到上面
				primGhost[i][j + N2 + NG][k][RHO] = primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][RHO];
				primGhost[i][j + N2 + NG][k][UU] = primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][UU];
				primGhost[i][j + N2 + NG][k][U1] = primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][U1];
				primGhost[i][j + N2 + NG][k][U2] = -primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][U2];
				primGhost[i][j + N2 + NG][k][U3] = primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][U3];
				primGhost[i][j + N2 + NG][k][B1] = primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][B1];
				primGhost[i][j + N2 + NG][k][B2] = -primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][B2];
				primGhost[i][j + N2 + NG][k][B3] = primGhost[i][N2 + NG - 1 - j][(k + N3 / 2) % N3][B3];
			}
		}
	}

	// 3. 鬼化X3方向

	for (int i = 0; i < N1 + 2 * NG; i++)
	{

		for (int j = 0; j < N2 + 2 * NG; j++)
		{

			for (int k = 0; k < NG; k++)
			{
				// 1) 把上面的格子移动到下面
				primGhost(i, j, k, RHO) = primGhost[i][j][N3 + NG - 2 + k][RHO];
				primGhost(i, j, k, UU) = primGhost[i][j][N3 + NG - 2 + k][UU];
				primGhost(i, j, k, U1) = primGhost[i][j][N3 + NG - 2 + k][U1];
				primGhost(i, j, k, U2) = -primGhost[i][j][N3 + NG - 2 + k][U2];
				primGhost(i, j, k, U3) = primGhost[i][j][N3 + NG - 2 + k][U3];
				primGhost[i][j][k][B1] = primGhost[i][j][N3 + NG - 2 + k][B1];
				primGhost[i][j][k][B2] = -primGhost[i][j][N3 + NG - 2 + k][B2];
				primGhost[i][j][k][B3] = primGhost[i][j][N3 + NG - 2 + k][B3];
				// 2) 把下面的格子移动到上面
				primGhost(i, j, k, RHO) = primGhost[i][j][k + N3 + NG][RHO];
				primGhost(i, j, k, UU) = primGhost[i][j][k + N3 + NG][UU];
				primGhost(i, j, k, U1) = primGhost[i][j][k + N3 + NG][U1];
				primGhost(i, j, k, U2) = -primGhost[i][j][k + N3 + NG][U2];
				primGhost(i, j, k, U3) = primGhost[i][j][k + N3 + NG][U3];
				primGhost[i][j][k][B1] = primGhost[i][j][k + N3 + NG][B1];
				primGhost[i][j][k][B2] = -primGhost[i][j][k + N3 + NG][B2];
				primGhost[i][j][k][B3] = primGhost[i][j][k + N3 + NG][B3];
			}
		}
	}

}

void interpolate() {

	for (int i = 0; i < N1 + 1; i++)
	{

		for (int j = 0; j < N2 + 1; j++)
		{

			for (int k = 0; k < N3 + 1; k++)
				for (int index = 0; index < NPRIM; index++)
				{
					primLR[LEFT][0][i][j][k][index] = primGhost[i + NG][j + NG][k + NG][index] - MC((primGhost[i + NG + 1][j + NG][k + NG][index] - primGhost[i + NG - 1][j + NG][k + NG][index]) / (2 * dx1),
						2 * (primGhost[i + NG + 1][j + NG][k + NG][index] - primGhost[i + NG][j + NG][k + NG][index]) / (dx1),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG - 1][j + NG][k + NG][index]) / (dx1)) * dx1 / 2;
					primLR[RIGHT][1][i][j][k][index] = primGhost[i + NG - 1][j + NG][k + NG][index] + MC((primGhost[i + NG][j + NG][k + NG][index] - primGhost[i][j + NG][k + NG][index]) / (2 * dx1),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG - 1][j + NG][k + NG][index]) / (dx1),
						2 * (primGhost[i + NG - 1][j + NG][k + NG][index] - primGhost[i][j + NG][k + NG][index]) / (dx1)) * dx1 / 2;

					primLR[LEFT][1][i][j][k][index] = primGhost[i + NG][j + NG + 1][k + NG][index] - MC((primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG - 1][k + NG][index]) / (2 * dx2),
						2 * (primGhost[i + NG][j + NG + 1][k + NG][index] - primGhost[i + NG][j + NG][k + NG][index]) / (dx2),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG - 1][k + NG][index]) / (dx2)) * dx2 / 2;
					primLR[RIGHT][1][i][j][k][index] = primGhost[i + NG][j + NG][k + NG][index] + MC((primGhost[i + NG][j + NG + 1][k + NG][index] - primGhost[i + NG][j][k + NG][index]) / (2 * dx2),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG - 1][k + NG][index]) / (dx2),
						2 * (primGhost[i + NG][j + NG - 1][k + NG][index] - primGhost[i + NG][j][k + NG][index]) / (dx2)) * dx2 / 2;

					primLR[LEFT][2][i][j][k][index] = primGhost[i + NG][j + NG][k + NG][index] - MC((primGhost[i + NG][j + NG][k + NG + 1][index] - primGhost[i + NG][j + NG][k + NG - 1][index]) / (2 * dx3),
						2 * (primGhost[i + NG][j + NG][k + NG + 1][index] - primGhost[i + NG][j + NG][k + NG][index]) / (dx3),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG][k + NG - 1][index]) / (dx3)) * dx3 / 2;
					primLR[RIGHT][2][i][j][k][index] = primGhost[i + NG][j + NG][k + NG - 1][index] + MC((primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG][k][index]) / (2 * dx3),
						2 * (primGhost[i + NG][j + NG][k + NG][index] - primGhost[i + NG][j + NG][k + NG - 1][index]) / (dx3),
						2 * (primGhost[i + NG][j + NG][k + NG - 1][index] - primGhost[i + NG][j + NG][k][index]) / (dx3)) * dx3 / 2;
				}
		}
	}
}

void primLR2conLR() {

	for (int i = 0; i < N1; i++)
	{

		for (int j = 0; j < N2; j++)
		{

			for (int k = 0; k < N3; k++)
			{

				for (int LR = 0; LR < 2; LR++)
				{

					for (int comp = 0; comp < 3; comp++)
					{
						Eigen::Vector3d u{ primLR[LR][comp](i, j, k, U1) ,primLR[LR][comp](i, j, k, U2) ,primLR[LR][comp](i, j, k, U3) };
						Eigen::Vector3d B{ primLR[LR][comp][i][j][k][B1] ,primLR[LR][comp][i][j][k][B2] ,primLR[LR][comp][i][j][k][B3] };
						auto usq = square(i, j, k, u, metricFuncHalfField[comp]);
						auto Bsq = square(i, j, k, B, metricFuncHalfField[comp]);
						double Gamma = sqrt(1 + usq);
						auto vsq = square(i, j, k, u / Gamma, metricFuncHalfField[comp]);
						auto Bv = dot(i, j, k, u / Gamma, B, metricFuncHalfField[comp]);
						conLR[LR][comp][i][j][k][0] = Gamma * primLR[LR][comp](i, j, k, RHO);
						conLR[LR][comp][i][j][k][1] = (primLR[LR][comp](i, j, k, RHO) + gam / (gam - 1) * primLR[LR][comp](i, j, k, UU)) * pow(Gamma, 2) - primLR[LR][comp](i, j, k, UU) + 0.5 * (Bsq * (1 + vsq) - pow(Bv, 2)) - Gamma * primLR[LR][comp](i, j, k, RHO);
						conLR[LR][comp][i][j][k][2] = (primLR[LR][comp](i, j, k, RHO) + gam / (gam - 1) * primLR[LR][comp](i, j, k, UU)) * pow(Gamma, 2) * primLR[LR][comp](i, j, k, U1) + Bsq * primLR[LR][comp](i, j, k, U1) - Bv * primLR[LR][comp][i][j][k][B1];
						conLR[LR][comp][i][j][k][3] = (primLR[LR][comp](i, j, k, RHO) + gam / (gam - 1) * primLR[LR][comp](i, j, k, UU)) * pow(Gamma, 2) * primLR[LR][comp](i, j, k, U2) + Bsq * primLR[LR][comp](i, j, k, U2) - Bv * primLR[LR][comp][i][j][k][B2];
						conLR[LR][comp][i][j][k][4] = (primLR[LR][comp](i, j, k, RHO) + gam / (gam - 1) * primLR[LR][comp](i, j, k, UU)) * pow(Gamma, 2) * primLR[LR][comp](i, j, k, U3) + Bsq * primLR[LR][comp](i, j, k, U3) - Bv * primLR[LR][comp][i][j][k][B3];
						conLR[LR][comp][i][j][k][5] = primLR[LR][comp][i][j][k][B1];
						conLR[LR][comp][i][j][k][6] = primLR[LR][comp][i][j][k][B2];
						conLR[LR][comp][i][j][k][7] = primLR[LR][comp][i][j][k][B3];
					}
				}
			}
		}
	}
}

void primLR2srcLR() {

	for (int i = 0; i < N1; i++)
	{

		for (int j = 0; j < N2; j++)
		{

			for (int k = 0; k < N3; k++)
			{

				for (int LR = 0; LR < 2; LR++)
				{

					for (int comp = 0; comp < 3; comp++)
					{
						Eigen::Vector3d u{ primLR[LR][comp](i, j, k, U1) ,primLR[LR][comp](i, j, k, U2) ,primLR[LR][comp](i, j, k, U3) };
						Eigen::Vector3d B{ primLR[LR][comp][i][j][k][B1] ,primLR[LR][comp][i][j][k][B2] ,primLR[LR][comp][i][j][k][B3] };
						Eigen::Vector3d S{ conLR[LR][comp][i][j][k][2] ,conLR[LR][comp][i][j][k][3] ,conLR[LR][comp][i][j][k][4] };
						auto usq = square(i, j, k, u, metricFuncHalfField[comp]);
						auto Bsq = square(i, j, k, B, metricFuncHalfField[comp]);
						double Gamma = sqrt(1 + usq);
						auto vsq = square(i, j, k, u / Gamma, metricFuncHalfField[comp]);
						auto Bv = dot(i, j, k, u / Gamma, B, metricFuncHalfField[comp]);
						auto metric = metricFuncHalfField[comp](i, j, k);
						auto metricDiff = metricDiffHalfField[comp];
						auto Sb = dot(i, j, k, S, metricDiff(i, j, k, 1).betaVec(), metricFuncHalfField[comp]);
						// W^{ij}
						Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim(i, j, k, UU) + 0.5 * (Bsq * (1 + vsq) - pow(Bv, 2))) * metric.gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - Bv * u / Gamma * B.transpose();
						Eigen::Matrix3d betaDiff;
						betaDiff << metricDiff(i, j, k, 1).betaVec()(0), metricDiff(i, j, k, 2).betaVec()(0), metricDiff(i, j, k, 3).betaVec()(0),
							metricDiff(i, j, k, 1).betaVec()(1), metricDiff(i, j, k, 2).betaVec()(1), metricDiff(i, j, k, 3).betaVec()(1),
							metricDiff(i, j, k, 1).betaVec()(2), metricDiff(i, j, k, 2).betaVec()(2), metricDiff(i, j, k, 3).betaVec()(2);
						srcLR[LR][comp][i][j][k][0] = 0;
						srcLR[LR][comp][i][j][k][1] = 0.5 * contract(W, (metric.betaVec()(0) * metricDiff(i, j, k, 1).gamma() + metric.betaVec()(1) * metricDiff(i, j, k, 2).gamma() + metric.betaVec()(2) * metricDiff(i, j, k, 3).gamma()))
							+ contract(W * metric.gamma(), betaDiff)
							- (metric.gamma().inverse() * S)(0) * alphaDiffHalfField[comp][i][j][k][1] - (metric.gamma().inverse() * S)(1) * alphaDiffHalfField[comp][i][j][k][2] - (metric.gamma().inverse() * S)(2) * alphaDiffHalfField[comp][i][j][k][3];
						srcLR[LR][comp][i][j][k][2] = 0.5 * metric.alpha() * contract(W, metricDiff(i, j, k, 1).gamma()) + Sb - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffHalfField[comp][i][j][k][1];
						srcLR[LR][comp][i][j][k][3] = 0.5 * metric.alpha() * contract(W, metricDiff(i, j, k, 2).gamma()) + Sb - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffHalfField[comp][i][j][k][2];
						srcLR[LR][comp][i][j][k][4] = 0.5 * metric.alpha() * contract(W, metricDiff(i, j, k, 3).gamma()) + Sb - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffHalfField[comp][i][j][k][3];
						srcLR[LR][comp][i][j][k][5] = 0;
						srcLR[LR][comp][i][j][k][6] = 0;
						srcLR[LR][comp][i][j][k][7] = 0;
					}
				}
			}
		}
	}
}

void primLR2fluxLR() {

	for (int i = 0; i < N1 + 1; i++)
		for (int j = 0; j < N2 + 1; j++)
			for (int k = 0; k < N3 + 1; k++)
				for (int LR = 0; LR < 2; LR++)
					for (int comp = 0; comp < 3; comp++)
					{
						Eigen::Vector3d u{ primLR[LR][comp](i, j, k, U1) ,primLR[LR][comp](i, j, k, U2) ,primLR[LR][comp](i, j, k, U3) };
						Eigen::Vector3d B{ primLR[LR][comp][i][j][k][B1] ,primLR[LR][comp][i][j][k][B2] ,primLR[LR][comp][i][j][k][B3] };
						Eigen::Vector3d S{ conLR[LR][comp][i][j][k][2] ,conLR[LR][comp][i][j][k][3] ,conLR[LR][comp][i][j][k][4] };
						auto usq = square(i, j, k, u, metricFuncHalfField[comp]);
						auto Bsq = square(i, j, k, B, metricFuncHalfField[comp]);
						double Gamma = sqrt(1 + usq);
						auto vsq = square(i, j, k, u / Gamma, metricFuncHalfField[comp]);
						auto Bv = dot(i, j, k, u / Gamma, B, metricFuncHalfField[comp]);
						auto metric = metricFuncHalfField[comp](i, j, k);
						auto V = metric.alpha() * u / Gamma - metric.betaVec();
						// W^{ij}
						Eigen::Matrix3d W = S * (u / Gamma).transpose() + (primLR[LR][comp](i, j, k, UU) + 0.5 * (Bsq * (1 + vsq) - pow(Bv, 2))) * metric.gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - Bv * u / Gamma * B.transpose();
						fluxLR[LR][comp][i][j][k][0] = V(comp) * conLR[LR][comp][i][j][k][0];
						fluxLR[LR][comp][i][j][k][1] = metric.alpha() * (conLR[LR][comp][i][j][k][2 + comp] - u(comp) / Gamma * conLR[LR][comp][i][j][k][0]) - metric.betaVec()(comp) * conLR[LR][comp][i][j][k][1];
						fluxLR[LR][comp][i][j][k][2] = (metric.alpha() * W * metric.gamma())(comp, 0) - metric.betaVec()(comp) * conLR[LR][comp][i][j][k][2];
						fluxLR[LR][comp][i][j][k][3] = (metric.alpha() * W * metric.gamma())(comp, 1) - metric.betaVec()(comp) * conLR[LR][comp][i][j][k][3];
						fluxLR[LR][comp][i][j][k][4] = (metric.alpha() * W * metric.gamma())(comp, 2) - metric.betaVec()(comp) * conLR[LR][comp][i][j][k][4];
						fluxLR[LR][comp][i][j][k][5] = V(comp) * conLR[LR][comp][i][j][k][5] - V(0) * conLR[LR][comp][i][j][k][5 + comp];
						fluxLR[LR][comp][i][j][k][6] = V(comp) * conLR[LR][comp][i][j][k][6] - V(1) * conLR[LR][comp][i][j][k][5 + comp];
						fluxLR[LR][comp][i][j][k][7] = V(comp) * conLR[LR][comp][i][j][k][7] - V(2) * conLR[LR][comp][i][j][k][5 + comp];
					}
}

void primLR2cLR() {

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int PN = 0; PN < 2; PN++)
					for (int LR = 0; LR < 2; LR++)
						for (int comp = 0; comp < 3; comp++)
							if (prim(i, j, k, RHO))
							{
								Eigen::Vector3d u{ primLR[LR][comp](i, j, k, U1) ,primLR[LR][comp](i, j, k, U2) ,primLR[LR][comp](i, j, k, U3) };
								Eigen::Vector3d B{ primLR[LR][comp][i][j][k][B1] ,primLR[LR][comp][i][j][k][B2] ,primLR[LR][comp][i][j][k][B3] };
								auto usq = square(i, j, k, u, metricFuncHalfField[comp]);
								auto Bsq = square(i, j, k, B, metricFuncHalfField[comp]);
								double Gamma = sqrt(1 + usq);
								auto vsq = square(i, j, k, u / Gamma, metricFuncHalfField[comp]);
								auto Bv = dot(i, j, k, u / Gamma, B, metricFuncHalfField[comp]);
								auto metric = metricFuncHalfField[comp](i, j, k);
								auto u0 = Gamma / metric.alpha();
								Eigen::Vector3d ui = { Gamma * (prim(i, j, k, U1) - metric.betaVec()(0)),
												Gamma * (prim(i, j, k, U2) - metric.betaVec()(1)),
												Gamma * (prim(i, j, k, U3) - metric.betaVec()(2)) };
								auto cs_square = gam * prim(i, j, k, UU) / (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU));
								auto cA_square = (Bsq * (1 - vsq) + pow(Bv, 2)) / (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU) + Bsq * (1 - vsq) + pow(Bv, 2));
								auto vf_square = cA_square + cs_square - cA_square * cs_square;
								auto metricInv = metric.m.inverse();
								c[PN][LR][comp][i][j][k] = (pow(vf_square, 2) * metricInv(0, comp) - pow(1 - vf_square, 2) * u0 * ui(comp)) / (pow(vf_square, 2) * metricInv(0, 0) - pow(1 - vf_square, 2) * u0 * u0) + (2 * PN - 1) * sqrt(abs(
									pow((pow(vf_square, 2) * metricInv(0, comp) - pow(1 - vf_square, 2) * u0 * ui(comp)) / (pow(vf_square, 2) * metricInv(0, 0) - pow(1 - vf_square, 2) * u0 * u0), 2)
									- (vf_square * metricInv(comp, comp) - (1 - vf_square) * ui(comp) * ui(comp)) / (vf_square * metricInv(0, 0) - (1 - vf_square) * u0 * u0)
								));
							}
}

void calFluxHHL() {

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int LR = 0; LR < 2; LR++)
					for (int comp = 0; comp < 3; comp++)
						for (int l = 0; l < 8; l++)
						{
							auto c_max = max(0, c[POS][RIGHT][comp][i][j][k], c[POS][LEFT][comp][i][j][k]);
							auto c_min = -min(0, c[NEG][LEFT][comp][i][j][k], c[NEG][LEFT][comp][i][j][k]);
							fluxHLL[comp][i][j][k][l] = c_max + c_min ? (c_min * fluxLR[RIGHT][comp][i][j][k][l] + c_max * fluxLR[LEFT][comp][i][j][k][l] - c_max * c_min * (conLR[RIGHT][comp][i][j][k][l] - conLR[LEFT][comp][i][j][k][l])) / (c_max + c_min) : 0;
						}
}

void calFluxTVDLF() {

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int LR = 0; LR < 2; LR++)
					for (int comp = 0; comp < 3; comp++)
						for (int l = 0; l < 8; l++)
						{
							auto c_max = max(0, c[POS][RIGHT][comp][i][j][k], c[POS][LEFT][comp][i][j][k]);
							auto c_min = -min(0, c[NEG][LEFT][comp][i][j][k], c[NEG][LEFT][comp][i][j][k]);
							auto c = max(c_max, c_min);
							fluxTVDLF[comp][i][j][k][l] = 0.5 * (fluxLR[RIGHT][comp][i][j][k][l] + fluxLR[LEFT][comp][i][j][k][l]) - 0.5 * c * (conLR[RIGHT][comp][i][j][k][l] - conLR[LEFT][comp][i][j][k][l]);
						}
}

void prim2con()
{

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim(i, j, k, U1) ,prim(i, j, k, U2) ,prim(i, j, k, U3) };
				Eigen::Vector3d B{ prim[i][j][k][B1] ,prim[i][j][k][B2] ,prim[i][j][k][B3] };
				auto usq = square(i + NG, j + NG, k + NG, u, metricFuncField);
				auto Bsq = square(i + NG, j + NG, k + NG, B, metricFuncField);
				double Gamma = sqrt(1 + usq);
				auto vsq = square(i + NG, j + NG, k + NG, u / Gamma, metricFuncField);
				auto Bv = dot(i + NG, j + NG, k + NG, u / Gamma, B, metricFuncField);
				con[i][j][k][0] = Gamma * prim(i, j, k, RHO);
				con[i][j][k][1] = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) - prim(i, j, k, UU) + 0.5 * (Bsq * (1 + vsq) - pow(Bv, 2)) - Gamma * prim(i, j, k, RHO);
				con[i][j][k][2] = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) * prim(i, j, k, U1) + Bsq * prim(i, j, k, U1) - Bv * prim[i][j][k][B1];
				con[i][j][k][3] = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) * prim(i, j, k, U2) + Bsq * prim(i, j, k, U2) - Bv * prim[i][j][k][B2];
				con[i][j][k][4] = (prim(i, j, k, RHO) + gam / (gam - 1) * prim(i, j, k, UU)) * pow(Gamma, 2) * prim(i, j, k, U3) + Bsq * prim(i, j, k, U3) - Bv * prim[i][j][k][B3];
				con[i][j][k][5] = prim[i][j][k][B1];
				con[i][j][k][6] = prim[i][j][k][B2];
				con[i][j][k][7] = prim[i][j][k][B3];
			}
}

void prim2src() {

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d u{ prim(i, j, k, U1) ,prim(i, j, k, U2) ,prim(i, j, k, U3) };
				Eigen::Vector3d B{ prim[i][j][k][B1] ,prim[i][j][k][B2] ,prim[i][j][k][B3] };
				Eigen::Vector3d S{ con[i][j][k][2] ,con[i][j][k][3] ,con[i][j][k][4] };
				auto usq = square(i + NG, j + NG, k + NG, u, metricFuncField);
				auto Bsq = square(i + NG, j + NG, k + NG, B, metricFuncField);
				double Gamma = sqrt(1 + usq);
				auto vsq = square(i + NG, j + NG, k + NG, u / Gamma, metricFuncField);
				auto Bv = dot(i + NG, j + NG, k + NG, u / Gamma, B, metricFuncField);
				auto metric = metricFuncField(i + NG, j + NG, k + NG);
				auto metricDiff = metricDiffField;
				auto Sb = dot(i + NG, j + NG, k + NG, S, metricDiff(i, j, k, 1).betaVec(), metricFuncField);
				// W^{ij}
				Eigen::Matrix3d W = S * (u / Gamma).transpose() + (prim(i, j, k, UU) + 0.5 * (Bsq * (1 + vsq) - pow(Bv, 2))) * metric.gamma().inverse() - B * B.transpose() / pow(Gamma, 2) - Bv * u / Gamma * B.transpose();
				Eigen::Matrix3d betaDiff;
				betaDiff << metricDiff(i, j, k, 1).betaVec()(0), metricDiff(i, j, k, 2).betaVec()(0), metricDiff(i, j, k, 3).betaVec()(0),
					metricDiff(i, j, k, 1).betaVec()(1), metricDiff(i, j, k, 2).betaVec()(1), metricDiff(i, j, k, 3).betaVec()(1),
					metricDiff(i, j, k, 1).betaVec()(2), metricDiff(i, j, k, 2).betaVec()(2), metricDiff(i, j, k, 3).betaVec()(2);
				src[i][j][k][0] = 0;
				src[i][j][k][1] = 0.5 * contract(W, (metric.betaVec()(0) * metricDiff(i, j, k, 1).gamma() + metric.betaVec()(1) * metricDiff(i, j, k, 2).gamma() + metric.betaVec()(2) * metricDiff(i, j, k, 3).gamma()))
					+ contract(W * metric.gamma(), betaDiff)
					- (metric.gamma().inverse() * S)(0) * alphaDiffField[i][j][k][1] - (metric.gamma().inverse() * S)(1) * alphaDiffField[i][j][k][2] - (metric.gamma().inverse() * S)(2) * alphaDiffField[i][j][k][3];
				src[i][j][k][2] = 0.5 * metric.alpha() * contract(W, metricDiff(i, j, k, 1).gamma()) + Sb - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffField[i][j][k][1];
				src[i][j][k][3] = 0.5 * metric.alpha() * contract(W, metricDiff(i, j, k, 2).gamma()) + Sb - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffField[i][j][k][2];
				src[i][j][k][4] = 0.5 * metric.alpha() * contract(W, metricDiff(i, j, k, 3).gamma()) + Sb - (con[i][j][k][0] + con[i][j][k][1]) * alphaDiffField[i][j][k][3];
				src[i][j][k][5] = 0;
				src[i][j][k][6] = 0;
				src[i][j][k][7] = 0;
			}
}

double f(int i, int j, int k, double D, double tau, Eigen::Vector3d S, Eigen::Vector3d B, double x) {
	auto Bsq = square(i + NG, j + NG, k + NG, B, metricFuncField);
	auto SB = dot(i + NG, j + NG, k + NG, S, B, metricFuncField);
	auto Gamma = 1 / sqrt(1 - square(i, j, k, S + SB * B / ksi[i][j][k], metricFuncField) / pow(ksi[i][j][k] + Bsq, 2));
	return x - (gam - 1) / gam * (x - Gamma * D) / pow(Gamma, 2) - tau - D + Bsq - 0.5 * (Bsq / pow(Gamma, 2) + pow(SB, 2) / pow(x, 2));
}

double df(int i, int j, int k, double D, double tau, Eigen::Vector3d S, Eigen::Vector3d B, double x) {
	return (f(i, j, k, D, tau, S, B, x + SMALL) - f(i, j, k, D, tau, S, B, x)) / SMALL;
}

void con2prim() {

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
			{
				Eigen::Vector3d S{ con[i][j][k][2] ,con[i][j][k][3] ,con[i][j][k][4] };
				Eigen::Vector3d B{ con[i][j][k][5] ,con[i][j][k][6] ,con[i][j][k][7] };
				auto Bsq = square(i + NG, j + NG, k + NG, B, metricFuncField);
				auto SB = dot(i + NG, j + NG, k + NG, S, B, metricFuncField);
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
				auto Gamma = 1 / sqrt(1 - square(i, j, k, S + SB * B / ksi[i][j][k], metricFuncField) / pow(ksi[i][j][k] + Bsq, 2));
				prim(i, j, k, RHO) = D / Gamma;
				prim(i, j, k, UU) = (gam - 1) / gam * (ksi[i][j][k] - Gamma * D) / pow(Gamma, 2);
				prim[i][j][k][U0] = Gamma / metricFuncField(i + NG, j + NG, k + NG).alpha();
				prim(i, j, k, U1) = (S(0) + SB * B(0) / ksi[i][j][k]) / (ksi[i][j][k] + Bsq) * Gamma;
				prim(i, j, k, U2) = (S(1) + SB * B(1) / ksi[i][j][k]) / (ksi[i][j][k] + Bsq) * Gamma;
				prim(i, j, k, U3) = (S(2) + SB * B(2) / ksi[i][j][k]) / (ksi[i][j][k] + Bsq) * Gamma;
				prim[i][j][k][B1] = B(0);
				prim[i][j][k][B2] = B(1);
				prim[i][j][k][B3] = B(2);
				prim[i][j][k][BSQ] = Bsq;
			}
}
