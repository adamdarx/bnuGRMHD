#pragma once
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <complex.h>
#include <time.h>
#include "utils.h"

#define CLOCKS_PER_SEC ((clock_t)1000)

//given a Jacobian matrix, transform a vector from coord1 to coord2
void vect_trans(double vect1[NDIM], double J[NDIM][NDIM]) {
    double vect2[NDIM];
    for (int i = 0; i < NDIM; i++) {
        vect2[i] = 0.0;
        for (int j = 0; j < NDIM; j++) {
            vect2[i] += J[i][j] * vect1[j];
        }
    }
    for (int i = 0; i < NDIM; i++) {
        vect1[i] = vect2[i];
    }
}

//caculate u^t from u^i
double ut_cal(double g[NDIM][NDIM], double u1, double u2, double u3) {
    double AA = g[0][0];
    double BB = 2. * (g[0][1] * u1 + g[0][2] * u2 + g[0][3] * u3);
    double CC = 1. + g[1][1] * u1 * u1 + g[2][2] * u2 * u2 + g[3][3] * u3 * u3 + 2. * (g[1][2] * u1 * u2 + g[1][3] * u1 * u3 + g[2][3] * u2 * u3);
    double DD = BB * BB - 4. * AA * CC;
    return -(BB + sqrt(DD)) / (2. * AA);
}

//convert 4-velocity u^i from bl coord to mks coord
void convert_ui(int i, int j, int k) {
    double uu[NDIM];
    uu[1] = primInit[i][j][k][U1];
    uu[2] = primInit[i][j][k][U2];
    uu[3] = primInit[i][j][k][U3];
    uu[0] = ut_cal(gdd_bl[i][j][k], uu[1], uu[2], uu[3]);
    vect_trans(uu, J_bl2ks[i][j][k]);
    vect_trans(uu, J_ks2mks[i][j][k]);
    primInit[i][j][k][U0] = uu[0];
    primInit[i][j][k][U1] = uu[1];
    primInit[i][j][k][U2] = uu[2];
    primInit[i][j][k][U3] = uu[3];
}

//caculate the constant l
double lfish_calc(double r)
{
    return(
        ((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
            ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
                sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
                ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) * (2. + r))) /
                sqrt(1 + (2. * a) / pow(r, 1.5) - 3. / r))) /
        (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) * (pow(a, 2) + (-2. + r) * r))
        );
}

//compute bsq
double bsq_cal(int i, int j, int k) {
    double b[NDIM], bsq;  // b^mu and bsq
    b[0] = 0.;
    for (int m = 1; m < 4; m++) {
        for (int n = 0; n < 4; n++) {
            b[0] += gdd_mks[i][j][k][m][n] * primInit[i][j][k][B1 + m - 1] * primInit[i][j][k][U0 + n - 1];
        }
    }
    for (int m = 1; m < 4; m++) {
        b[m] = (primInit[i][j][k][B1 + m - 1] + b[0] * primInit[i][j][k][U1 + m - 1]) / (SMALL + primInit[i][j][k][U0]);
    }
    bsq = 0.;
    for (int m = 0; m < 4; m++) {
        for (int n = 0; n < 4; n++) {
            bsq += gdd_mks[i][j][k][m][n] * b[m] * b[n];
        }
    }
    bsq = fabs(bsq);
    return bsq;
}

//compute B from A_phi in MKS coord, then return bsq_max
double compute_B_from_A(double A[N1 + 1][N2 + 1][N3 + 1]) {
    double bsq_max = 0.;
    double bsq, r;
    for (int i = 1; i < N1; i++) {
        for (int j = 1; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                primInit[i][j][k][B1] = (A[i - 1][j][k] - A[i - 1][j - 1][k] +
                    A[i][j][k] - A[i][j - 1][k]) / (2. * dx2 * gdet_mks[i][j][k]);
                primInit[i][j][k][B2] = -(A[i][j - 1][k] - A[i - 1][j - 1][k] +
                    A[i][j][k] - A[i - 1][j][k]) / (2. * dx1 * gdet_mks[i][j][k]);
                primInit[i][j][k][B3] = 0.;
                r = BL_coord1[i][j][k];
                if (r >= rin) {
                    bsq = bsq_cal(i, j, k);
                    if (bsq > bsq_max) bsq_max = bsq;
                }
            }
        }
    }
    return bsq_max;
}

//fix primitive variable by adding density rho
void fix(double primInit[N1][N2][N3][NPRIM]) {
    double r, rho_floor, ug_floor, bsq, sigma;
    for (int i = 1; i < N1; i++) {
        for (int j = 1; j < N2; j++) {
            for (int k = 1; k < N3; k++) {
                r = BL_coord1[i][j][k];
                rho_floor = RHOMIN * pow(r, -3. / 2.);
                ug_floor = UUMIN * pow(r, -3. / 2. * gam);
                if (primInit[i][j][k][RHO] < rho_floor) primInit[i][j][k][RHO] = rho_floor;
                if (primInit[i][j][k][UU] < ug_floor) primInit[i][j][k][UU] = ug_floor;
                bsq = bsq_cal(i, j, k);
                sigma = bsq / primInit[i][j][k][RHO];
                if (sigma > SIGMAMAX) primInit[i][j][k][RHO] = bsq / SIGMAMAX;
            }
        }
    }
}

void fix(Eigen::Tensor<double, 4> prim) {
    double r, rho_floor, ug_floor, bsq, sigma;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                r = BL_coord1[i][j][k];
                rho_floor = RHOMIN * pow(r, -3. / 2.);
                ug_floor = UUMIN * pow(r, -3. / 2. * gam);
                if (prim(i + NG, j + NG, k + NG, RHO) < rho_floor) prim(i + NG, j + NG, k + NG, RHO) = rho_floor;
                if (prim(i + NG, j + NG, k + NG, UU) < ug_floor) prim(i + NG, j + NG, k + NG, UU) = ug_floor;
                bsq = bsq_cal(i, j, k);
                sigma = bsq / prim(i, j, k, RHO);
                if (sigma > SIGMAMAX) prim(i, j, k, RHO) = bsq / SIGMAMAX;
            }
        }
    }
}

//write bin
int write_bin(FILE* fp) {
    size_t total_elements = N1 * N2 * N3 * NPRIM;
    if (fwrite(&prim, sizeof(double), total_elements, fp) != total_elements) {
        return -1;
    }
    return 0;
}

//read bin
int read_bin(Eigen::Tensor<double, 4> prim, FILE* fp) {
    size_t total_elements = N1 * N2 * N3 * NPRIM;
    if (fread(&prim, sizeof(double), total_elements, fp) != total_elements) {
        return -1;
    }
    return 0;
}

void init() {
    //intermediate quantity
    double r, theta, phi, r2, sinth, sinth2, costh, costh2, sigma, delta, AA;
    double a2 = a * a;
    double rhor = 1 + sqrt(1 - a2);
    //disk parameter and primary variable
    double rho, u, ur, uh, up, lnh, hm1, expm2chi, up1, rancval;
    double rhomax = 0.;
    double umax = 0.;
    double l = lfish_calc(rmax);
    double thin = PI / 2.;          //r = rin, theta = thin
    double sthin = sin(thin);
    double cthin = cos(thin);
    double DDin = rin * rin - 2. * rin + a * a;
    double AAin = (rin * rin + a2) * (rin * rin + a2) - DDin * a2 * sthin * sthin;
    double SSin = rin * rin + a2 * cthin * cthin;
    //the factor from KS metric to MKS metric
    double tfac = 1;
    double rfac, hfac;
    double pfac = 1;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                /**************************************************************************************************
                (1) get the BL, KS and MKS coords at the grid points
                ***************************************************************************************************/

                //MKS grid
                Xgrid1[i][j][k] = X1min + i * dx1;
                Xgrid2[i][j][k] = X2min + j * dx2;
                Xgrid3[i][j][k] = X3min + k * dx3;
                //KS grid
                KS_coord1[i][j][k] = exp(Xgrid1[i][j][k]) + R0;
                KS_coord2[i][j][k] = Xgrid2[i][j][k] + h / 2. * sin(2. * Xgrid2[i][j][k]);
                KS_coord3[i][j][k] = Xgrid3[i][j][k];
                //BL grid
                BL_coord1[i][j][k] = KS_coord1[i][j][k];
                BL_coord2[i][j][k] = KS_coord2[i][j][k];
                //to save computing power, we caculate these quantity at once
                r = BL_coord1[i][j][k];
                theta = BL_coord2[i][j][k];
                r2 = r * r;
                sinth = sin(theta);
                costh = cos(theta);
                sinth2 = sinth * sinth;
                costh2 = costh * costh;
                sigma = r2 + a2 * costh2;
                delta = r2 - 2. * r + a2;
                AA = (r2 + a2) * (r2 + a2) - delta * a2 * sinth2;
                rfac = r;
                hfac = 1. + h * cos(2. * Xgrid2[i][j][k]);

                BL_coord3[i][j][k] = KS_coord3[i][j][k] - a * KS_coord1[i][j][k] / delta;
                BL_coord3[i][j][k] = fmod(BL_coord3[i][j][k], 2. * PI);
                if (BL_coord3[i][j][k] < 0) {
                    BL_coord3[i][j][k] += 2. * PI;
                }
                phi = BL_coord3[i][j][k];
                /**************************************************************************************************
                (2) get BL, KS and MKS metric at the grid points
                ***************************************************************************************************/

                //bl metric
                gdd_bl[i][j][k][0][0] = -1. + 2. * r / sigma;
                gdd_bl[i][j][k][0][1] = 0.;
                gdd_bl[i][j][k][0][2] = 0.;
                gdd_bl[i][j][k][0][3] = -2. * r * a * sinth2 / sigma;

                gdd_bl[i][j][k][1][0] = gdd_bl[i][j][k][0][1];
                gdd_bl[i][j][k][1][1] = sigma / delta;
                gdd_bl[i][j][k][1][2] = 0.;
                gdd_bl[i][j][k][1][3] = 0.;

                gdd_bl[i][j][k][2][0] = gdd_bl[i][j][k][0][2];
                gdd_bl[i][j][k][2][1] = gdd_bl[i][j][k][1][2];
                gdd_bl[i][j][k][2][2] = sigma;
                gdd_bl[i][j][k][2][3] = 0.;

                gdd_bl[i][j][k][3][0] = gdd_bl[i][j][k][0][3];
                gdd_bl[i][j][k][3][1] = gdd_bl[i][j][k][1][3];
                gdd_bl[i][j][k][3][2] = gdd_bl[i][j][k][2][3];
                gdd_bl[i][j][k][3][3] = (r2 + a2 + 2. * r * a2 * sinth2 / sigma) * sinth2;

                //ks metric
                gdd_ks[i][j][k][0][0] = gdd_bl[i][j][k][0][0];//tt component: the same as bl metric
                gdd_ks[i][j][k][0][1] = 2. * r / sigma;
                gdd_ks[i][j][k][0][2] = 0.;
                gdd_ks[i][j][k][0][3] = -2. * a * r * sinth2 / sigma;

                gdd_ks[i][j][k][1][0] = gdd_ks[i][j][k][0][1];
                gdd_ks[i][j][k][1][1] = 1. + 2. * r / sigma;
                gdd_ks[i][j][k][1][2] = 0.;
                gdd_ks[i][j][k][1][3] = -a * sinth2 * (1. + 2. * r / sigma);

                gdd_ks[i][j][k][2][0] = gdd_ks[i][j][k][0][2];
                gdd_ks[i][j][k][2][1] = gdd_ks[i][j][k][1][2];
                gdd_ks[i][j][k][2][2] = sigma;
                gdd_ks[i][j][k][2][3] = 0.;

                gdd_ks[i][j][k][3][0] = gdd_ks[i][j][k][0][3];
                gdd_ks[i][j][k][3][1] = gdd_ks[i][j][k][1][3];
                gdd_ks[i][j][k][3][2] = gdd_ks[i][j][k][2][3];
                gdd_ks[i][j][k][3][3] = sinth2 * (sigma + a2 * sinth2 * (1. + 2. * r / sigma));

                //mks metric
                gdd_mks[i][j][k][0][0] = gdd_ks[i][j][k][0][0] * tfac * tfac;
                gdd_mks[i][j][k][0][1] = gdd_ks[i][j][k][0][1] * tfac * rfac;
                gdd_mks[i][j][k][0][2] = gdd_ks[i][j][k][0][2] * tfac * hfac;
                gdd_mks[i][j][k][0][3] = gdd_ks[i][j][k][0][3] * tfac * pfac;

                gdd_mks[i][j][k][1][0] = gdd_mks[i][j][k][0][1];
                gdd_mks[i][j][k][1][1] = gdd_ks[i][j][k][1][1] * rfac * rfac;
                gdd_mks[i][j][k][1][2] = gdd_ks[i][j][k][1][2] * rfac * hfac;
                gdd_mks[i][j][k][1][3] = gdd_ks[i][j][k][1][3] * rfac * pfac;

                gdd_mks[i][j][k][2][0] = gdd_mks[i][j][k][0][2];
                gdd_mks[i][j][k][2][1] = gdd_mks[i][j][k][1][2];
                gdd_mks[i][j][k][2][2] = gdd_ks[i][j][k][2][2] * hfac * hfac;
                gdd_mks[i][j][k][2][3] = gdd_ks[i][j][k][2][3] * hfac * pfac;

                gdd_mks[i][j][k][3][0] = gdd_mks[i][j][k][0][3];
                gdd_mks[i][j][k][3][1] = gdd_mks[i][j][k][1][3];
                gdd_mks[i][j][k][3][2] = gdd_mks[i][j][k][2][3];
                gdd_mks[i][j][k][3][3] = gdd_ks[i][j][k][3][3] * pfac * pfac;

                //bl guu
                guu_bl[i][j][k][0][0] = (-(r2 + a2) * (r2 + a2) / (SMALL + delta) + a2 * sinth2) / sigma;
                guu_bl[i][j][k][0][1] = 0.;
                guu_bl[i][j][k][0][2] = 0.;
                guu_bl[i][j][k][0][3] = -2. * a * r / ((SMALL + delta) * sigma);

                guu_bl[i][j][k][1][0] = guu_bl[i][j][k][0][1];
                guu_bl[i][j][k][1][1] = delta / sigma;
                guu_bl[i][j][k][1][2] = 0.;
                guu_bl[i][j][k][1][3] = 0.;

                guu_bl[i][j][k][2][0] = guu_bl[i][j][k][0][2];
                guu_bl[i][j][k][2][1] = guu_bl[i][j][k][1][2];
                guu_bl[i][j][k][2][2] = 1. / sigma;
                guu_bl[i][j][k][2][3] = 0.;

                guu_bl[i][j][k][3][0] = guu_bl[i][j][k][0][3];
                guu_bl[i][j][k][3][1] = guu_bl[i][j][k][1][3];
                guu_bl[i][j][k][3][2] = guu_bl[i][j][k][2][3];
                guu_bl[i][j][k][3][3] = (1. / sinth2 - a2 / delta) / sigma;

                //ks guu
                guu_ks[i][j][k][0][0] = -1. - 2. * r / sigma;
                guu_ks[i][j][k][0][1] = 2. * r / sigma;
                guu_ks[i][j][k][0][2] = 0.;
                guu_ks[i][j][k][0][3] = 0.;

                guu_ks[i][j][k][1][0] = guu_ks[i][j][k][0][1];
                guu_ks[i][j][k][1][1] = (r2 - 2. * r + a2) / sigma;
                guu_ks[i][j][k][1][2] = 0.;
                guu_ks[i][j][k][1][3] = a / sigma;

                guu_ks[i][j][k][2][0] = guu_ks[i][j][k][0][2];
                guu_ks[i][j][k][2][1] = guu_ks[i][j][k][1][2];
                guu_ks[i][j][k][2][2] = 1. / sigma;
                guu_ks[i][j][k][2][3] = 0.;

                guu_ks[i][j][k][3][0] = guu_ks[i][j][k][0][3];
                guu_ks[i][j][k][3][1] = guu_ks[i][j][k][1][3];
                guu_ks[i][j][k][3][2] = guu_ks[i][j][k][2][3];
                guu_ks[i][j][k][3][3] = 1. / (sigma * sinth2);

                //mks guu
                guu_mks[i][j][k][0][0] = guu_ks[i][j][k][0][0] / (tfac * tfac);
                guu_mks[i][j][k][0][1] = guu_ks[i][j][k][0][1] / (tfac * rfac);
                guu_mks[i][j][k][0][2] = guu_ks[i][j][k][0][2] / (tfac * hfac);
                guu_mks[i][j][k][0][3] = guu_ks[i][j][k][0][3] / (tfac * pfac);

                guu_mks[i][j][k][1][0] = guu_mks[i][j][k][0][1];
                guu_mks[i][j][k][1][1] = guu_ks[i][j][k][1][1] / (rfac * rfac);
                guu_mks[i][j][k][1][2] = guu_ks[i][j][k][1][2] / (rfac * hfac);
                guu_mks[i][j][k][1][3] = guu_ks[i][j][k][1][3] / (rfac * pfac);

                guu_mks[i][j][k][2][0] = guu_mks[i][j][k][0][2];
                guu_mks[i][j][k][2][1] = guu_mks[i][j][k][1][2];
                guu_mks[i][j][k][2][2] = guu_ks[i][j][k][2][2] / (hfac * hfac);
                guu_mks[i][j][k][2][3] = guu_ks[i][j][k][2][3] / (hfac * pfac);

                guu_mks[i][j][k][3][0] = guu_mks[i][j][k][0][3];
                guu_mks[i][j][k][3][1] = guu_mks[i][j][k][1][3];
                guu_mks[i][j][k][3][2] = guu_mks[i][j][k][2][3];
                guu_mks[i][j][k][3][3] = guu_ks[i][j][k][3][3] / (pfac * pfac);

                //sqrt(-g)
                gdet_bl[i][j][k] = sinth * sigma;
                gdet_ks[i][j][k] = gdet_bl[i][j][k];
                gdet_mks[i][j][k] = gdet_ks[i][j][k] * (tfac * rfac * hfac * pfac);

                /**************************************************************************************************
                (3) get Jacobian matrix at the grid points
                ***************************************************************************************************/

                //BL coord to KS coord
                J_bl2ks[i][j][k][0][0] = 1.;
                J_bl2ks[i][j][k][0][1] = 2. * r / (SMALL + delta);
                J_bl2ks[i][j][k][0][2] = 0.;
                J_bl2ks[i][j][k][0][3] = 0.;

                J_bl2ks[i][j][k][1][0] = 0.;
                J_bl2ks[i][j][k][1][1] = 1.;
                J_bl2ks[i][j][k][1][2] = 0.;
                J_bl2ks[i][j][k][1][3] = 0.;

                J_bl2ks[i][j][k][2][0] = 0.;
                J_bl2ks[i][j][k][2][1] = 0.;
                J_bl2ks[i][j][k][2][2] = 1.;
                J_bl2ks[i][j][k][2][3] = 0.;

                J_bl2ks[i][j][k][3][0] = 0.;
                J_bl2ks[i][j][k][3][1] = a / (SMALL + delta);
                J_bl2ks[i][j][k][3][2] = 0.;
                J_bl2ks[i][j][k][3][3] = 1.;

                //KS coord to BL coord
                J_ks2bl[i][j][k][0][0] = 1.;
                J_ks2bl[i][j][k][0][1] = -2. * r / (SMALL + delta);
                J_ks2bl[i][j][k][0][2] = 0.;
                J_ks2bl[i][j][k][0][3] = 0.;

                J_ks2bl[i][j][k][1][0] = 0.;
                J_ks2bl[i][j][k][1][1] = 1.;
                J_ks2bl[i][j][k][1][2] = 0.;
                J_ks2bl[i][j][k][1][3] = 0.;

                J_ks2bl[i][j][k][2][0] = 0.;
                J_ks2bl[i][j][k][2][1] = 0.;
                J_ks2bl[i][j][k][2][2] = 1.;
                J_ks2bl[i][j][k][2][3] = 0.;

                J_ks2bl[i][j][k][3][0] = 0.;
                J_ks2bl[i][j][k][3][1] = -a / (SMALL + delta);
                J_ks2bl[i][j][k][3][2] = 0.;
                J_ks2bl[i][j][k][3][3] = 1.;

                //KS coord to MKS coord
                J_ks2mks[i][j][k][0][0] = 1.;
                J_ks2mks[i][j][k][0][1] = 0.;
                J_ks2mks[i][j][k][0][2] = 0.;
                J_ks2mks[i][j][k][0][3] = 0.;

                J_ks2mks[i][j][k][1][0] = 0.;
                J_ks2mks[i][j][k][1][1] = rfac;
                J_ks2mks[i][j][k][1][2] = 0.;
                J_ks2mks[i][j][k][1][3] = 0.;

                J_ks2mks[i][j][k][2][0] = 0.;
                J_ks2mks[i][j][k][2][1] = 0.;
                J_ks2mks[i][j][k][2][2] = hfac;
                J_ks2mks[i][j][k][2][3] = 0.;

                J_ks2mks[i][j][k][3][0] = 0.;
                J_ks2mks[i][j][k][3][1] = 0.;
                J_ks2mks[i][j][k][3][2] = 0.;
                J_ks2mks[i][j][k][3][3] = 1.;

                //MKS coord to KS coord
                J_mks2ks[i][j][k][0][0] = 1.;
                J_mks2ks[i][j][k][0][1] = 0.;
                J_mks2ks[i][j][k][0][2] = 0.;
                J_mks2ks[i][j][k][0][3] = 0.;

                J_mks2ks[i][j][k][1][0] = 0.;
                J_mks2ks[i][j][k][1][1] = 1. / rfac;
                J_mks2ks[i][j][k][1][2] = 0.;
                J_mks2ks[i][j][k][1][3] = 0.;

                J_mks2ks[i][j][k][2][0] = 0.;
                J_mks2ks[i][j][k][2][1] = 0.;
                J_mks2ks[i][j][k][2][2] = 1. / hfac;
                J_mks2ks[i][j][k][2][3] = 0.;

                J_mks2ks[i][j][k][3][0] = 0.;
                J_mks2ks[i][j][k][3][1] = 0.;
                J_mks2ks[i][j][k][3][2] = 0.;
                J_mks2ks[i][j][k][3][3] = 1.;

                /**************************************************************************************************
                (4) get primitive variables at the grid points
                ***************************************************************************************************/
                primInit[i][j][k][RHO] = 0.;
                primInit[i][j][k][UU] = 0.;
                primInit[i][j][k][U1] = 0.;
                primInit[i][j][k][U2] = 0.;
                primInit[i][j][k][U3] = 0.;
                primInit[i][j][k][B1] = 0.;
                primInit[i][j][k][B2] = 0.;
                primInit[i][j][k][B3] = 0.;

                if (r >= rin) {
                    lnh = 0.5 * log((1. + sqrt(1. + 4. * (l * l * sigma * sigma) * delta / (AA * AA * sinth2))) / (sigma * delta / AA))
                        - 0.5 * sqrt(1. + 4. * (l * l * sigma * sigma) * delta / (AA * AA * sinth2))
                        - 2. * a * r * l / AA
                        - 0.5 * log((1. + sqrt(1. + 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin * sthin * sthin))) / (SSin * DDin / AAin))
                        + 0.5 * sqrt(1. + 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin * sthin * sthin))
                        + 2. * a * rin * l / AAin;
                }
                else
                    lnh = 1.;

                /*inside torus*/
                if (r >= rin && lnh >= 0.) {
                    hm1 = exp(lnh) - 1.;
                    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));
                    rho = rho * (1. + 0.1 * sin(phi));
                    u = kappa * pow(rho, gam) / (gam - 1.);
                    ur = 0.;
                    uh = 0.;
                    expm2chi = sigma * sigma * delta / (AA * AA * sinth2);
                    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
                    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * sigma * delta) + sqrt(sigma / AA) * up1 / sinth;
                    primInit[i][j][k][RHO] = rho;
                    if (rho > rhomax) rhomax = rho;
                    //rancval = ranc(0);
                    primInit[i][j][k][UU] = u;// * (1. + 4.e-2 * (rancval - 0.5));
                    if (u > umax && r > rin) umax = u;
                    primInit[i][j][k][U1] = ur;
                    primInit[i][j][k][U2] = uh;
                    primInit[i][j][k][U3] = up;
                }
                if (r >= rin) {
                    convert_ui(i, j, k);      //transfrom U1, U2, U3 from bl coord to mks coord
                }
            }
        }
    }
    /*Normalize the density rho*/
    printf("rhomax before normalization: %lf \n", rhomax);
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                primInit[i][j][k][RHO] /= rhomax;
                primInit[i][j][k][UU] /= rhomax;
            }
        }
    }
    /*Aphi at the corner*/
    double rho_ave, q, bsq_max;
    A[N1][N2][N3] = 0.;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                A[i][j][k] = 0.;
                r = BL_coord1[i][j][k];
                if (r >= rin) {
                    rho_ave = 0.25 * (primInit[i][j][k][RHO] + primInit[i - 1][j][k][RHO] + primInit[i][j - 1][k][RHO] + primInit[i - 1][j - 1][k][RHO]);
                    q = rho_ave - 0.2;
                    if (q > 0.) {
                        A[i][j][k] = q;
                        //printf("A: %lf \n", A[i][j][k]);
                    }
                }
            }
        }
    }
    /*caculate B, then get bsq_max and pg_max*/
    bsq_max = compute_B_from_A(A);
    printf("bsq_max before normalization: %lf \n", bsq_max);
    double pg, pg_max = 0.;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                r = BL_coord1[i][j][k];
                if (r >= rin) {
                    pg = (gam - 1.) * primInit[i][j][k][UU];
                    if (pg > pg_max) pg_max = pg;
                }
            }
        }
    }
    printf("pg_max: %lf \n", pg_max);
    /*normalize with beta*/
    double beta_min, norm;
    beta_min = pg_max / (0.5 * bsq_max);
    printf("beta_min before normalization: %.10lf \n", beta_min);
    printf("target beta: %lf \n", beta);
    norm = sqrt(beta_min / beta);
    printf("normalization factor: %lf \n", norm);
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                primInit[i][j][k][B1] *= norm;
                primInit[i][j][k][B2] *= norm;
                primInit[i][j][k][B3] *= norm;
                r = BL_coord1[i][j][k];
                if (r >= rin) {
                    primInit[i][j][k][BSQ] = bsq_cal(i, j, k);
                }
                else {
                    primInit[i][j][k][BSQ] = 0.;
                }
            }
        }
    }
    /*fix p by adding rho*/
    fix(primInit);
}
