/*
bnuGRMHD ©️ 2025
Date: 2024/02/02
*/
#include <ctime>
#include <fstream>
#include "Metric.h"
#include "utils.h"
#include "init.h"
#include "metric/mks.h"

int main(int argc, char* argv[])
{
	amrex::Initialize(argc,argv);
	// amrex relevant settings
	// settings for normal meshgrid
	amrex::IntVect dom_lo(0, 0, 0);
	amrex::IntVect dom_hi(N1 - 1, N2 - 1, N3 - 1);
	amrex::Box domain(dom_lo, dom_hi);
	amrex::BoxArray ba(domain);
	ba.maxSize(max_grid_size);
	amrex::DistributionMapping dm(ba);
	amrex::RealBox real_box ({X1min, X2min, X3min}, {X1max, X2max, X3max});
    amrex::Geometry geom(domain, &real_box);
	amrex::GpuArray<amrex::Real,3> dx = geom.CellSizeArray();
	// settings for ghost meshgrid
	amrex::IntVect dom_loGhost(0, 0, 0);
	amrex::IntVect dom_hiGhost(N1 + 2 * NG - 1, N2 + 2 * NG - 1, N3 + 2 * NG - 1);
	amrex::Box domainGhost(dom_loGhost, dom_hiGhost);
	amrex::BoxArray baGhost(domainGhost);
	baGhost.maxSize(max_grid_size + NG);
	amrex::DistributionMapping dmGhost(baGhost);
	amrex::RealBox real_boxGhost ({X1min - NG * dx[0], X2min - NG * dx[1], X3min - NG * dx[2]}, {X1max + NG * dx[0], X2max + NG * dx[1], X3max + NG * dx[2]});
    amrex::Geometry geomGhost(domainGhost, &real_boxGhost);
	// settings for half meshgrid
	amrex::IntVect dom_loHalf(0, 0, 0);
	amrex::IntVect dom_hiHalf(N1 + 1 - 1, N2 + 1 - 1, N3 + 1 - 1);
	amrex::Box domainHalf(dom_loHalf, dom_hiHalf);
	amrex::BoxArray baHalf(domainHalf);
	baHalf.maxSize(max_grid_size + 1);
	amrex::DistributionMapping dmHalf(baHalf);
	amrex::RealBox real_boxHalf ({X1min - 0.5 * dx[0], X2min - 0.5 * dx[1], X3min - 0.5 * dx[2]}, {X1max + 0.5 * dx[0], X2max + 0.5 * dx[1], X3max + 0.5 * dx[2]});
    amrex::Geometry geomHalf(domainHalf, &real_boxHalf);

	alphaDiffField.define(ba, dm, 4, 0);
	for(int dim = 0; dim < 3; dim++)
		alphaDiffHalfField[dim].define(baHalf, dmHalf, 4, 0);
	prim.define(ba, dm, NPRIM, 0);
	con.define(ba, dm, 8, 0);
	src.define(ba, dm, 8, 0);
	ksi.define(ba, dm, 1, 0);
	primGhost.define(baGhost, dmGhost, NPRIM, 2);
	for(int LR = 0; LR < 2; LR++)
		for(int comp = 0; comp < 3; comp++)
		{
			primLR[LR][comp].define(baHalf, dmHalf, NPRIM, 1);
			conLR[LR][comp].define(baHalf, dmHalf, 8, 1);
			srcLR[LR][comp].define(baHalf, dmHalf, 8, 1);
			fluxLR[LR][comp].define(baHalf, dmHalf, 8, 1);
		}
	for(int dim = 0; dim < 3; dim++)
	{
		fluxHLL[dim].define(baHalf, dmHalf, 8, 0);
		fluxTVDLF[dim].define(baHalf, dmHalf, 8, 0);
		fluxLLF[dim].define(baHalf, dmHalf, 8, 1);
		fluxSmoothLLF[dim].define(baHalf, dmHalf, 8, 1);
	}
	for(int PN = 0; PN < 2; PN++)
		for(int LR = 0; LR < 2; LR++)
			for(int comp = 0; comp < 3; comp++)
				c[PN][LR][comp].define(baHalf, dmHalf, 1, 0);

	auto totalTime = 0.;
	auto totalPhysicalTime = 0.;
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
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(alphaDiffField); mfi.isValid(); ++mfi)
	{
		const amrex::Box& bx = mfi.validbox();
		amrex::Array4<amrex::Real> const& a = alphaDiffField[mfi].array();
		const auto lo = lbound(bx);
		const auto hi = ubound(bx);
		for (int i = lo.x; i <= hi.x; i++) {
			for (int j = lo.y; j <= hi.y; j++) {
				for (int k = lo.z; k <= hi.z; k++) {
					a(i, j, k, 0) = 0;
					a(i, j, k, 1) = (metricFuncField[i + NG + 1][j + NG][k + NG].alpha() - metricFuncField[i + NG - 1][j + NG][k + NG].alpha()) / (2 * dx1);
					a(i, j, k, 2) = (metricFuncField[i + NG][j + NG + 1][k + NG].alpha() - metricFuncField[i + NG][j + NG - 1][k + NG].alpha()) / (2 * dx2);
					a(i, j, k, 3) = (metricFuncField[i + NG][j + NG][k + NG + 1].alpha() - metricFuncField[i + NG][j + NG][k + NG - 1].alpha()) / (2 * dx3);
				}
			}
		}
	}

	for (int comp = 0; comp < 3; comp++)
	{
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
		for (amrex::MFIter mfi(alphaDiffHalfField[comp]); mfi.isValid(); ++mfi)
		{
			const amrex::Box& bx = mfi.validbox();
			amrex::Array4<amrex::Real> const& a = alphaDiffHalfField[comp][mfi].array();
			const auto lo = lbound(bx);
			const auto hi = ubound(bx);
			for (int i = lo.x; i <= hi.x; i++) {
				for (int j = lo.y; j <= hi.y; j++) {
					for (int k = lo.z; k <= hi.z; k++) {
						a(i, j, k, 0) = 0;
						a(i, j, k, 1) = (metricFuncField[i + 1][j][k].alpha() - metricFuncField[i][j][k].alpha()) / (dx1);
						a(i, j, k, 2) = (metricFuncField[i][j + 1][k].alpha() - metricFuncField[i][j][k].alpha()) / (dx2);
						a(i, j, k, 3) = (metricFuncField[i][j][k + 1].alpha() - metricFuncField[i][j][k].alpha()) / (dx3);
					}
				}
			}
		}
	}

	for (int i = 0; i < N1 + 2 * NG; i++)
		for (int j = 0; j < N2 + 2 * NG; j++)
			for (int k = 0; k < N3 + 2 * NG; k++)
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
	WriteSingleLevelPlotfile("plt000", prim, {"RHO", "UU", "U0", "U1", "U2", "U3", "B1", "B2", "B3", "BSQ"}, geom, 0., 0);
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
	for (amrex::MFIter mfi(ksi); mfi.isValid(); ++mfi)
	{
		const amrex::Box& bx = mfi.validbox();
		amrex::Array4<amrex::Real> const& a = ksi[mfi].array();
		amrex::Array4<amrex::Real const> const& b = prim[mfi].array();
		const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int i = lo.x; i <= hi.x; i++) {
            for (int j = lo.y; j <= hi.y; j++) {
                for (int k = lo.z; k <= hi.z; k++) {
					a(i, j, k) = (b(i, j, k, RHO) + gam / (gam - 1) * b(i, j, k, UU)) * (1 + pow(b(i, j, k, U1), 2) + pow(b(i, j, k, U2), 2) + pow(b(i, j, k, U3), 2));
				}
			}
		}
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
		
		for(int comp = 0; comp < 3; comp++)
		{
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
			for (amrex::MFIter mfi(fluxLLF[comp]); mfi.isValid(); ++mfi)
			{
				const amrex::Box& bx = mfi.validbox();
				amrex::Array4<amrex::Real> const& a = fluxLLF[comp][mfi].array();
				amrex::Array4<amrex::Real const> const& b = fluxHLL[comp][mfi].const_array();
				amrex::Array4<amrex::Real const> const& c = fluxTVDLF[comp][mfi].const_array();
				const auto lo = lbound(bx);
				const auto hi = ubound(bx);
				for (int i = lo.x; i <= hi.x; i++) {
					for (int j = lo.y; j <= hi.y; j++) {
						for (int k = lo.z; k <= hi.z; k++) {
							for(int l = 0; l < 8; l++)
								a(i,j,k,l) = theta * b(i,j,k,l) + (1 - theta) * c(i,j,k,l);
						}
					}
				}
			}
		}

		prim2con();

		prim2src();
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
		for (amrex::MFIter mfi(c[POS][RIGHT][0]); mfi.isValid(); ++mfi)
		{
			const amrex::Box& bx = mfi.validbox();

			amrex::Array4<amrex::Real const> const& cPR0 = c[POS][RIGHT][0][mfi].array();
			amrex::Array4<amrex::Real const> const& cPL0 = c[POS][LEFT][0][mfi].array();
			amrex::Array4<amrex::Real const> const& cNR0 = c[NEG][RIGHT][0][mfi].array();
			amrex::Array4<amrex::Real const> const& cNL0 = c[NEG][LEFT][0][mfi].array();

			amrex::Array4<amrex::Real const> const& cPR1 = c[POS][RIGHT][1][mfi].array();
			amrex::Array4<amrex::Real const> const& cPL1 = c[POS][LEFT][1][mfi].array();
			amrex::Array4<amrex::Real const> const& cNR1 = c[NEG][RIGHT][1][mfi].array();
			amrex::Array4<amrex::Real const> const& cNL1 = c[NEG][LEFT][1][mfi].array();

			amrex::Array4<amrex::Real const> const& cPR2 = c[POS][RIGHT][2][mfi].array();
			amrex::Array4<amrex::Real const> const& cPL2 = c[POS][LEFT][2][mfi].array();
			amrex::Array4<amrex::Real const> const& cNR2 = c[NEG][RIGHT][2][mfi].array();
			amrex::Array4<amrex::Real const> const& cNL2 = c[NEG][LEFT][2][mfi].array();

			const auto lo = lbound(bx);
			const auto hi = ubound(bx);
			for (int i = lo.x; i <= hi.x; i++) {
				for (int j = lo.y; j <= hi.y; j++) {
					for (int k = lo.z; k <= hi.z; k++) {
						auto c1max = max(0, cPR0(i, j, k), cPL0(i, j, k));
						auto c1min = abs(min(-0, cNR0(i, j, k), cNL0(i, j, k)));
						auto c2max = max(0, cPR1(i, j, k), cPL1(i, j, k));
						auto c2min = abs(min(-0, cNR1(i, j, k), cNL1(i, j, k)));
						auto c3max = max(0, cPR2(i, j, k), cPL2(i, j, k));
						auto c3min = abs(min(-0, cNR2(i, j, k), cNL2(i, j, k)));
						auto c1 = max(c1max, c1min);
						auto c2 = max(c2max, c2min);
						auto c3 = max(c3max, c3min);
						Delta_t = min(cour * min(dx1 / (2 * c1), dx2 / (2 * c2), dx3 / (2 * c3)), Delta_t);
					}
				}
			}
		}

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
		for (amrex::MFIter mfi(con); mfi.isValid(); ++mfi)
		{
			const amrex::Box& bx = mfi.validbox();
			
			amrex::Array4<amrex::Real> const& a = con[mfi].array();
			amrex::Array4<amrex::Real const> const& b = src[mfi].array();
			amrex::Array4<amrex::Real const> const& c1 = fluxLLF[0][mfi].array();
			amrex::Array4<amrex::Real const> const& c2 = fluxLLF[1][mfi].array();
			amrex::Array4<amrex::Real const> const& c3 = fluxLLF[2][mfi].array();
			const auto lo = lbound(bx);
			const auto hi = ubound(bx);
			for (int i = lo.x; i <= hi.x; i++) {
				for (int j = lo.y; j <= hi.y; j++) {
					for (int k = lo.z; k <= hi.z; k++) {
						for(int l = 0; l < 8; l++)
							a(i, j, k, l) += b(i, j, k, l) * Delta_t / 2
							- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField[0][i + 1][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c1(i + 1, j, k, l) - sqrt(metricFuncHalfField[0][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c1(i, j, k, l))
							- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField[1][i][j + 1][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c2(i, j + 1, k, l) - sqrt(metricFuncHalfField[1][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c2(i, j, k, l))
							- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField[2][i][j][k + 1].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c3(i, j, k + 1, l) - sqrt(metricFuncHalfField[2][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c3(i, j, k, l));
					}
				}
			}
		}

		con2prim();

		ghostify();

		interpolate();

		primLR2conLR();

		primLR2srcLR();

		primLR2fluxLR();

		primLR2cLR();

		calFluxHHL();

		calFluxTVDLF();

		amrex::MultiFab::Copy(fluxLLF[0], fluxSmoothLLF[0], 0, 0, 8, 0);
		amrex::MultiFab::Copy(fluxLLF[1], fluxSmoothLLF[1], 0, 0, 8, 0);
		amrex::MultiFab::Copy(fluxLLF[2], fluxSmoothLLF[2], 0, 0, 8, 0);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
		for (amrex::MFIter mfi(fluxSmoothLLF[0]); mfi.isValid(); ++mfi)
		{
			const amrex::Box& bx = mfi.validbox();
			amrex::Array4<amrex::Real> const& a1 = fluxSmoothLLF[0][mfi].array();
			amrex::Array4<amrex::Real> const& a2 = fluxSmoothLLF[1][mfi].array();
			amrex::Array4<amrex::Real> const& a3 = fluxSmoothLLF[2][mfi].array();
			amrex::Array4<amrex::Real const> const& c1 = fluxLLF[0][mfi].array();
			amrex::Array4<amrex::Real const> const& c2 = fluxLLF[1][mfi].array();
			amrex::Array4<amrex::Real const> const& c3 = fluxLLF[2][mfi].array();
			const auto lo = lbound(bx);
			const auto hi = ubound(bx);
			for (int i = lo.x; i <= hi.x; i++) {
				for (int j = lo.y; j <= hi.y; j++) {
					for (int k = lo.z; k <= hi.z; k++) {
						a1(i, j, k, 6) = 0.125 * (2 * c1(i, j, k, 6) + c1(i, j + 1, k, 6) + c1(i, j - 1, k, 6) - c2(i, j, k, 5) - c2(i, j + 1, k, 5) - c2(i - 1, j, k, 5) - c2(i - 1, j + 1, k, 5));
						a2(i, j, k, 5) = 0.125 * (2 * c2(i, j, k, 5) + c2(i + 1, j, k, 5) + c2(i - 1, j, k, 5) - c1(i, j, k, 6) - c1(i + 1, j, k, 6) - c1(i, j - 1, k, 6) - c1(i + 1, j - 1, k, 6));
						a1(i, j, k, 7) = 0.125 * (2 * c1(i, j, k, 7) + c1(i, j, k + 1, 7) + c1(i, j, k - 1, 7) - c3(i, j, k, 5) - c3(i, j, k + 1, 7) - c3(i - 1, j, k, 5) - c3(i - 1, j, k + 1, 7));
						a3(i, j, k, 5) = 0.125 * (2 * c3(i, j, k, 5) + c3(i + 1, j, k, 5) + c3(i - 1, j, k, 5) - c1(i, j, k, 7) - c1(i + 1, j, k, 7) - c1(i, j, k - 1, 7) - c1(i + 1, j, k - 1, 7));
						a2(i, j, k, 7) = 0.125 * (2 * c2(i, j, k, 7) + c2(i, j, k + 1, 7) + c2(i, j, k - 1, 7) - c3(i, j, k, 6) - c3(i, j, k + 1, 7) - c3(i, j - 1, k, 6) - c3(i, j - 1, k + 1, 7));
						a3(i, j, k, 6) = 0.125 * (2 * c3(i, j, k, 6) + c3(i, j + 1, k, 6) + c3(i, j - 1, k, 6) - c2(i, j, k, 7) - c2(i, j + 1, k, 7) - c2(i, j, k - 1, 7) - c2(i, j + 1, k - 1, 7));
					}
				}
			}
		}

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
		for (amrex::MFIter mfi(con); mfi.isValid(); ++mfi)
		{
			const amrex::Box& bx = mfi.validbox();
			
			amrex::Array4<amrex::Real> const& a = con[mfi].array();
			amrex::Array4<amrex::Real const> const& b = src[mfi].array();
			amrex::Array4<amrex::Real const> const& c1 = fluxSmoothLLF[0][mfi].array();
			amrex::Array4<amrex::Real const> const& c2 = fluxSmoothLLF[1][mfi].array();
			amrex::Array4<amrex::Real const> const& c3 = fluxSmoothLLF[2][mfi].array();
			const auto lo = lbound(bx);
			const auto hi = ubound(bx);
			for (int i = lo.x; i <= hi.x; i++) {
				for (int j = lo.y; j <= hi.y; j++) {
					for (int k = lo.z; k <= hi.z; k++) {
						for(int l = 0; l < 8; l++)
							a(i, j, k, l) += b(i, j, k, l) * Delta_t / 2
							- Delta_t / (2 * dx1) * (sqrt(metricFuncHalfField[0][i + 1][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c1(i + 1, j, k, l) - sqrt(metricFuncHalfField[0][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c1(i, j, k, l))
							- Delta_t / (2 * dx2) * (sqrt(metricFuncHalfField[1][i][j + 1][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c2(i, j + 1, k, l) - sqrt(metricFuncHalfField[1][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c2(i, j, k, l))
							- Delta_t / (2 * dx3) * (sqrt(metricFuncHalfField[2][i][j][k + 1].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c3(i, j, k + 1, l) - sqrt(metricFuncHalfField[2][i][j][k].gamma().determinant() / metricFuncField[i + NG][j + NG][k + NG].gamma().determinant()) * c3(i, j, k, l));
					}
				}
			}
		}

		con2prim();

		fix();
		
		totalTime += clock() - start;
		totalPhysicalTime += Delta_t;
		if(int(totalPhysicalTime - Delta_t) / 10 != int(totalPhysicalTime) / 10)
		{
			char filename[16];
			sprintf(filename, "plt%0.3d", int(totalPhysicalTime) / 10);
			WriteSingleLevelPlotfile(filename, prim, {"RHO", "UU", "U0", "U1", "U2", "U3", "B1", "B2", "B3", "BSQ"}, geom, totalPhysicalTime, int(totalPhysicalTime) / 10);
		}
		amrex::Print() << "Epoch: " << epoch << "\tTime(ms): " << clock() - start << "\tPhysical Time: " << Delta_t << "\tTotal Physical Time: " << totalPhysicalTime << std::endl;
	}
	amrex::Finalize();
	return 0;
}
