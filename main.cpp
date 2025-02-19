/*
bnuGRMHD ©️ 2025
Date: 2024/02/02
*/
#include <cmath>
#include <ctime>
#include "Metric.h"
#include "utils.h"
#include "init.h"
#include "metric/mks.h"

int main(int argc, char* argv[])
{
	auto totalTime = 0.;
	auto totalPhysicalTime = 0.;

	{
		init_metric();
		for (int i = 0; i < N1 + 2 * NG; i++)
			for (int j = 0; j < N2 + 2 * NG; j++)
				for (int k = 0; k < N3 + 2 * NG; k++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
							metricFuncField[i][j][k].m(row, col) = metricFunc[row][col](X1min + (i - NG) * dx1, X2min + (j - NG) * dx2, X3min + (k - NG) * dx3);

		for (int i = 0; i < N1 + 1; i++)
			for (int j = 0; j < N2 + 1; j++)
				for (int k = 0; k < N3 + 1; k++)
					for (int row = 0; row < 4; row++)
						for (int col = 0; col < 4; col++)
						{
							metricFuncHalfField[0][i][j][k].m(row, col) = metricFunc[row][col](X1min + (2 * i - 1) * dx1 / 2, X2min + j * dx2, X3min + k * dx3);
							metricFuncHalfField[1][i][j][k].m(row, col) = metricFunc[row][col](X1min + i * dx1, X2min + (2 * j - 1) * dx2 / 2, X3min + k * dx3);
							metricFuncHalfField[2][i][j][k].m(row, col) = metricFunc[row][col](X1min + i * dx1, X2min + j * dx2, X3min + (2 * k - 1) * dx3 / 2);
						}

		// 利用中心差分计算alpha的导数
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					{
						alphaDiffField[i][j][k][0] = 0;
						alphaDiffField[i][j][k][1] = (metricFuncField[i + NG + 1][j + NG][k + NG].alpha() - metricFuncField[i + NG - 1][j + NG][k + NG].alpha()) / (2 * dx1);
						alphaDiffField[i][j][k][2] = (metricFuncField[i + NG][j + NG + 1][k + NG].alpha() - metricFuncField[i + NG][j + NG - 1][k + NG].alpha()) / (2 * dx2);
						alphaDiffField[i][j][k][3] = (metricFuncField[i + NG][j + NG][k + NG + 1].alpha() - metricFuncField[i + NG][j + NG][k + NG - 1].alpha()) / (2 * dx3);
					}

		for (int i = 0; i < N1 + 1; i++)
			for (int j = 0; j < N2 + 1; j++)
				for (int k = 0; k < N3 + 1; k++)
					for (int comp = 0; comp < 3; comp++)
					{
						alphaDiffHalfField[comp][i][j][k][0] = 0;
						alphaDiffHalfField[comp][i][j][k][1] = (metricFuncField[i + 1][j][k].alpha() - metricFuncField[i][j][k].alpha()) / (dx1);
						alphaDiffHalfField[comp][i][j][k][2] = (metricFuncField[i][j + 1][k].alpha() - metricFuncField[i][j][k].alpha()) / (dx2);
						alphaDiffHalfField[comp][i][j][k][3] = (metricFuncField[i][j][k + 1].alpha() - metricFuncField[i][j][k].alpha()) / (dx3);
					}

		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					for(int l = 0; l < 4; l++)
						for (int row = 0; row < 4; row++)
							for (int col = 0; col < 4; col++)
								metricDiffField[i][j][k][l].m(row, col) = metricDiff[row][col][l](X1min + (i - NG) * dx1, X2min + (j - NG) * dx2, X3min + (k - NG) * dx3);

		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					for (int l = 0; l < 4; l++)
						for (int row = 0; row < 4; row++)
							for (int col = 0; col < 4; col++)
							{
								metricDiffHalfField[0][i][j][k][l].m(row, col) = metricDiff[row][col][l](X1min + (2 * i - 1) * dx1 / 2, X2min + j * dx2, X3min + k * dx3);
								metricDiffHalfField[1][i][j][k][l].m(row, col) = metricDiff[row][col][l](X1min + i * dx1, X2min + (2 * j - 1) * dx2 / 2, X3min + k * dx3);
								metricDiffHalfField[2][i][j][k][l].m(row, col) = metricDiff[row][col][l](X1min + i * dx1, X2min + j * dx2, X3min + (2 * k - 1) * dx3 / 2);
							}
		init();
		char filename[32];
		sprintf(filename, "./data/data%0.4d.bin", 0);
		write_bin(fopen(filename, "wb"));
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					ksi[i][j][k] = (prim[i][j][k][RHO] + gam / (gam - 1) * prim[i][j][k][UU]) * (1 + pow(prim[i][j][k][U1], 2) + pow(prim[i][j][k][U2], 2) + pow(prim[i][j][k][U3], 2));

	}

	for(int epoch = 1; epoch <= epochNum; epoch++)
	{
		auto start = clock();

		double Delta_t = 1;

		ghostify();

		interpolate();
		
		primLR2conLR();
		
		primLR2srcLR();

		primLR2fluxLR();

		primLR2cLR();

		calFluxHHL();

		calFluxTVDLF();
		
#pragma omp parallel for
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
					for(int l = 0; l < 8; l++)
						for(int comp = 0; comp < 3; comp++)
							fluxLLF[comp][i][j][k][l] = theta * fluxHLL[comp][i][j][k][l] + (1 - theta) * fluxTVDLF[comp][i][j][k][l];

		prim2con();

		prim2src();

#pragma omp parallel for
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
				for (int k = 0; k < N3; k++)
				{
					auto c1max = max(0, c[POS][RIGHT][0][i][j][k], c[POS][LEFT][0][i][j][k]);
					auto c1min = abs(min(-0, c[NEG][RIGHT][0][i][j][k], c[NEG][LEFT][0][i][j][k]));
					auto c2max = max(0, c[POS][RIGHT][1][i][j][k], c[POS][LEFT][1][i][j][k]);
					auto c2min = abs(min(-0, c[NEG][RIGHT][1][i][j][k], c[NEG][LEFT][1][i][j][k]));
					auto c3max = max(0, c[POS][RIGHT][2][i][j][k], c[POS][LEFT][2][i][j][k]);
					auto c3min = abs(min(-0, c[NEG][RIGHT][2][i][j][k], c[NEG][LEFT][2][i][j][k]));
					auto c1 = max(c1max, c1min);
					auto c2 = max(c2max, c2min);
					auto c3 = max(c3max, c3min);
					Delta_t = min(cour * min(dx1 / (2 * c1), dx2 / (2 * c2), dx3 / (2 * c3)), Delta_t);
				}

#pragma omp parallel for
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
					for (int l = 0; l < 8; l++)
						con[i][j][k][l] += src[i][j][k][l] * Delta_t / 2
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField[0][i + 1][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxLLF[0][i + 1][j][k][l] - sqrt(metricFuncHalfField[0][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxLLF[0][i][j][k][l])
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField[1][i][j + 1][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxLLF[1][i][j + 1][k][l] - sqrt(metricFuncHalfField[1][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxLLF[1][i][j][k][l])
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField[2][i][j][k + 1].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxLLF[2][i][j][k + 1][l] - sqrt(metricFuncHalfField[2][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxLLF[2][i][j][k][l]);


		con2prim();

		ghostify();

		interpolate();

		primLR2conLR();

		primLR2srcLR();

		primLR2fluxLR();

		primLR2cLR();

		calFluxHHL();

		calFluxTVDLF();

#pragma omp parallel for
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
				{
					fluxSmoothLLF[0][i][j][k][6] = 0.125 * (2 * fluxLLF[0][i][j][k][6] + fluxLLF[0][i][j + 1][k][6] + fluxLLF[0][i][j - 1][k][6] - fluxLLF[1][i][j][k][5] - fluxLLF[1][i][j + 1][k][5] - fluxLLF[1][i - 1][j][k][5] - fluxLLF[1][i - 1][j + 1][k][5]);
					fluxSmoothLLF[1][i][j][k][5] = 0.125 * (2 * fluxLLF[1][i][j][k][5] + fluxLLF[1][i + 1][j][k][5] + fluxLLF[1][i - 1][j][k][5] - fluxLLF[0][i][j][k][6] - fluxLLF[0][i + 1][j][k][6] - fluxLLF[0][i][j - 1][k][6] - fluxLLF[0][i + 1][j - 1][k][6]);
					fluxSmoothLLF[0][i][j][k][7] = 0.125 * (2 * fluxLLF[0][i][j][k][7] + fluxLLF[0][i][j][k + 1][7] + fluxLLF[0][i][j][k - 1][7] - fluxLLF[2][i][j][k][5] - fluxLLF[2][i][j][k + 1][7] - fluxLLF[2][i - 1][j][k][5] - fluxLLF[2][i - 1][j][k + 1][7]);
					fluxSmoothLLF[2][i][j][k][5] = 0.125 * (2 * fluxLLF[2][i][j][k][5] + fluxLLF[2][i + 1][j][k][5] + fluxLLF[2][i - 1][j][k][5] - fluxLLF[0][i][j][k][7] - fluxLLF[0][i + 1][j][k][7] - fluxLLF[0][i][j][k - 1][7] - fluxLLF[0][i + 1][j][k - 1][7]);
					fluxSmoothLLF[1][i][j][k][7] = 0.125 * (2 * fluxLLF[1][i][j][k][7] + fluxLLF[1][i][j][k + 1][7] + fluxLLF[1][i][j][k - 1][7] - fluxLLF[2][i][j][k][6] - fluxLLF[2][i][j][k + 1][7] - fluxLLF[2][i][j - 1][k][6] - fluxLLF[2][i][j - 1][k + 1][7]);
					fluxSmoothLLF[2][i][j][k][6] = 0.125 * (2 * fluxLLF[2][i][j][k][6] + fluxLLF[2][i][j + 1][k][6] + fluxLLF[2][i][j - 1][k][6] - fluxLLF[1][i][j][k][7] - fluxLLF[1][i][j + 1][k][7] - fluxLLF[1][i][j][k - 1][7] - fluxLLF[1][i][j + 1][k - 1][7]);
				}

#pragma omp parallel for
		for (int i = 1; i < N1 - 1; i++)
			for (int j = 1; j < N2 - 1; j++)
				for (int k = 1; k < N3 - 1; k++)
					for (int l = 0; l < 8; l++)
						con[i][j][k][l] += src[i][j][k][l] * Delta_t / 2
						- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField[0][i + 1][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxSmoothLLF[0][i + 1][j][k][l] - sqrt(metricFuncHalfField[0][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxSmoothLLF[0][i][j][k][l])
						- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField[1][i][j + 1][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxSmoothLLF[1][i][j + 1][k][l] - sqrt(metricFuncHalfField[1][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxSmoothLLF[1][i][j][k][l])
						- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField[2][i][j][k + 1].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxSmoothLLF[2][i][j][k + 1][l] - sqrt(metricFuncHalfField[2][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * fluxSmoothLLF[2][i][j][k][l]);

		con2prim();

		fix();
		
		totalTime += clock() - start;
		totalPhysicalTime += Delta_t;
		std::cout << "Epoch: " << epoch << "\tTime(ms): " << clock() - start << "\tPhysical Time: " << Delta_t << "\tTotal Physical Time: " << totalPhysicalTime << std::endl;
		if (int(totalPhysicalTime) % 10 == 0 && int(totalPhysicalTime) / 10 > 0)
		{
			char filename[32];
			sprintf(filename, "./data/data%0.4d.bin", int(totalPhysicalTime) / 10);
			write_bin(fopen(filename, "wb"));
		}
	}
	return 0;
}
