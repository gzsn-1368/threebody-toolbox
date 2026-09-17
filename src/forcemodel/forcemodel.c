/**
 * FORCEMODEL USAGE
 *
 * SPICE utilities:
 *    ET = FORCEMODEL('str2et', DATE)
 *    [PTARG, LT] = FORCEMODEL('spkpos', TARG, ET, REF, ABCORR, OBS)
 *    [STATE, LT] = FORCEMODEL('spkezr', TARGET, ET, REF, ABCORR, OBS)
 *    ROTATE = FORCEMODEL('pxform', FROM, TO, ET)       % 3x3 or 3x3xN
 *    XFORM = FORCEMODEL('sxform', FROM, TO, ET)        % 6x6 or 6x6xN
 *    [ROTATE, AV] = FORCEMODEL('xf2rav', XFORM)
 *
 * Rotation utilities:
 *    ROTATE = FORCEMODEL('q2R', Q)
 *    DROTATEDQ = FORCEMODEL('q2Rjac', Q)
 *
 * Accelerations, forces, and magnetic field:
 *    [A, JAC, HESS] = FORCEMODEL('gravityfield', P, DEG, ORD, ROTATE, GM, RE, C, S, DOJAC, DOHESS)
 *    [A, JAC] = FORCEMODEL('fourthbody', ET, P, DOBODIES, DOJAC)
 *    [A, JAC] = FORCEMODEL('relativity', P, V, DOJAC)
 *    [A, JAC] = FORCEMODEL('srp', P, PSUN, RPCM, DOJAC)
 *    [FORCE, JAC, JERK, JERKJAC] = FORCEMODEL('srpnplate_q', P, PSUN, V, VSUN, OMEGA, Q, PLATES, DOJAC, DOJERK)
 *    [FORCE, JAC] = FORCEMODEL('srpnplate_R', P, PSUN, V, VSUN, ROTATE, PLATES, DOJAC)
 *    [A, JAC] = FORCEMODEL('earthalbedo', P, PEARTH, PSUN, RPCM, DOJAC)
 *    [B, JAC] = FORCEMODEL('magneticfield', P, DEG, ORD, ROTATE, RE, G, H, DOJAC)
 *
 * Jerks:
 *    JERK = FORCEMODEL('jerkgravityfield', P, V, A, OMEGA, JAC)
 *    [JERK, JAC] = FORCEMODEL('jerkfourthbody', ET, P, V, DOBODIES, DOJAC)
 *    JERK = FORCEMODEL('jerkrelativity', P, V, A, JAC)
 *    [JERK, JAC] = FORCEMODEL('jerksrp', P, V, PSUN, VSUN, RPCM, DOJAC)
 *    JERK = FORCEMODEL('jerkearthalbedo', P, V, PEARTH, VEARTH, PSUN, VSUN, RPCM)
 *
 * Torques:
 *    [TORQUE, JAC_P, JAC_Q] = FORCEMODEL('gravitytorque_q', JAC_A, Q, INERTIA, DOJAC, HESS_A)
 *    [TORQUE, JAC_P, JAC_R] = FORCEMODEL('gravitytorque_R', JAC_A, ROTATE, INERTIA, DOJAC, HESS_A)
 *    [TORQUE, JAC_P, JAC_Q] = FORCEMODEL('srpnplatetorque_q', FORCE, Q, DISPLACEMENTS, DOJAC, JAC_FORCE)
 *    [TORQUE, JAC_P, JAC_R] = FORCEMODEL('srpnplatetorque_R', FORCE, ROTATE, DISPLACEMENTS, DOJAC, JAC_FORCE)
 *
 * Eclipse model:
 *    [NU, GRAD] = FORCEMODEL('dualcone', P, PSUN, PEARTH)
 *    [NU, GRAD, NUDOT, NUDOTJAC] = FORCEMODEL('dualcone', P, PSUN, PEARTH, V, VSUN, VEARTH, DOJAC)
 * 
 *    Copyright 2025-2026 Guo Zisen.
 */

#include "mex.h"
#include "SpiceUsr.h"
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define IDX(row, col, nRows) ((row) + (col) * (nRows))
#define IDX3(row, col, page, nRows, nCols) ((row) + (col) * (nRows) + (page) * (nRows) * (nCols))
#define PUT3(i, j, k, val, U)         \
    do {                              \
        U[IDX3(i, j, k, 3, 3)] = val; \
        U[IDX3(i, k, j, 3, 3)] = val; \
        U[IDX3(j, i, k, 3, 3)] = val; \
        U[IDX3(j, k, i, 3, 3)] = val; \
        U[IDX3(k, i, j, 3, 3)] = val; \
        U[IDX3(k, j, i, 3, 3)] = val; \
    } while (0)

#define EPS                 1e-15                    // Numerical tolerance
#define PI                  3.141592653589793238     // Pi
#define MOONGM              4.902800066163796e+03    // Lunar gravitational parameter [km^2/s^2]
#define MOONRADIUS          1738                     // Lunar radius [km]
#define LIGHTSPEED          299792.458               // Speed of light [km/s]
#define SOLARLUMINOSITY     3.839e26                 // Solar luminosity [W]
#define EARTHRADIUS         6371                     // Earth radius [km]
#define EARTHREFLECTION     0.3                      // Earth reflectivity [dimensionless]
#define SUNRADIUS           696000                   // Solar radius [km]

static bool isKernelLoaded = false;

static const char* KERNELS[9] = {
    "de440.bsp",
    "earth_1962_240827_2124_combined.bpc",
    "EM_EarthCenteredRotation.fk",
    "EM_MoonCenteredRotation.fk",
    "MoonCenteredInertial.fk",
    "moon_pa_de440_200625.bpc",
    "moon_de440_220930.tf",
    "naif0012.tls",
    "pck00010.tpc"
};

static const char* BODYNAMES[9] = {
    "SUN",                  // Sun
    "MERCURY",              // Mercury
    "VENUS",                // Venus
    "MARS BARYCENTER",      // Mars barycenter
    "JUPITER BARYCENTER",   // Jupiter barycenter
    "SATURN BARYCENTER",    // Saturn barycenter
    "URANUS BARYCENTER",    // Uranus barycenter
    "NEPTUNE BARYCENTER",   // Neptune barycenter
    "PLUTO BARYCENTER"      // Pluto barycenter
};

static const double GM[9] = {
    1.327124400419393e11,   // Sun
    2.203178000000002e4,    // Mercury
    3.248585920000000e5,    // Venus
    4.282837521400002e4,    // Mars barycenter
    1.267127648000002e8,    // Jupiter barycenter
    3.794058520000000e7,    // Saturn barycenter
    5.794548600000008e6,    // Uranus barycenter
    6.836527100580023e6,    // Neptune barycenter
    9.770000000000007e2     // Pluto barycenter
};

typedef enum {
    STR2ET,
    SPKPOS,
    SPKEZR,
    PXFORM,
    SXFORM,
    XF2RAV,
    Q2R,
    Q2RJAC,
    GRAVITYFIELD,
    FOURTHBODY,
    RELATIVITY,
    SRP,
    SRPNPLATE_Q,
    SRPNPLATE_R,
    EARTHALBEDO,
    JERKGRAVITYFIELD,
    JERKFOURTHBODY,
    JERKRELATIVITY,
    JERKSRP,
    JERKEARTHALBEDO,
    MAGNETICFIELD,
    GRAVITYTORQUE_Q,
    GRAVITYTORQUE_R,
    SRPNPLATETORQUE_Q,
    SRPNPLATETORQUE_R,
    DUALCONE
} CMDID;

static const struct {
    const char* fname;
    CMDID id;
} _cmds[] = {
    { "str2et", STR2ET },
    { "spkpos", SPKPOS },
    { "spkezr", SPKEZR },
    { "pxform", PXFORM },
    { "sxform", SXFORM },
    { "xf2rav", XF2RAV },
    { "q2R", Q2R },
    { "q2Rjac", Q2RJAC },
    { "gravityfield", GRAVITYFIELD },
    { "fourthbody", FOURTHBODY },
    { "relativity", RELATIVITY },
    { "srp", SRP },
    { "srpnplate_q", SRPNPLATE_Q },
    { "srpnplate_R", SRPNPLATE_R },
    { "earthalbedo", EARTHALBEDO },
    { "jerkgravityfield", JERKGRAVITYFIELD },
    { "jerkfourthbody", JERKFOURTHBODY },
    { "jerkrelativity", JERKRELATIVITY },
    { "jerksrp", JERKSRP },
    { "jerkearthalbedo", JERKEARTHALBEDO },
    { "magneticfield", MAGNETICFIELD },
    { "gravitytorque_q", GRAVITYTORQUE_Q },
    { "gravitytorque_R", GRAVITYTORQUE_R },
    { "srpnplatetorque_q", SRPNPLATETORQUE_Q },
    { "srpnplatetorque_R", SRPNPLATETORQUE_R },
    { "dualcone", DUALCONE }
};

static const size_t _nCmds = 26;

// Compute the cross product of two three-vectors.
static void cross(const double* a, const double* b, double* out) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

// Convert a quaternion to a rotation matrix.
static void q2R(const double* q, double* R) {
    double q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];

    double q11, q22, q33, q44, q12, q13, q14, q23, q24, q34;
    q11 = q1 * q1;
    q22 = q2 * q2;
    q33 = q3 * q3;
    q44 = q4 * q4;
    q12 = 2 * q1 * q2;
    q13 = 2 * q1 * q3;
    q14 = 2 * q1 * q4;
    q23 = 2 * q2 * q3;
    q24 = 2 * q2 * q4;
    q34 = 2 * q3 * q4;

    R[0] = q11 - q22 - q33 + q44;
    R[1] = q12 - q34;
    R[2] = q13 + q24;
    R[3] = q12 + q34;
    R[4] = -q11 + q22 - q33 + q44;
    R[5] = q23 - q14;
    R[6] = q13 - q24;
    R[7] = q23 + q14;
    R[8] = -q11 - q22 + q33 + q44;
}

// Compute the Jacobian tensor of the quaternion-to-matrix mapping.
static void q2Rjac(const double*q, double* U) {
    double twoq1 = 2.0 * q[0], twoq2 = 2.0 * q[1], twoq3 = 2.0 * q[2], twoq4 = 2.0 * q[3];

    // Derivative with respect to q1.
    U[0] = twoq1; U[3] =  twoq2; U[6] =  twoq3;
    U[1] = twoq2; U[4] = -twoq1; U[7] =  twoq4;
    U[2] = twoq3; U[5] = -twoq4; U[8] = -twoq1;

    // Derivative with respect to q2.
    U[9] = -twoq2; U[12] = twoq1; U[15] = -twoq4;
    U[10] = twoq1; U[13] = twoq2; U[16] =  twoq3;
    U[11] = twoq4; U[14] = twoq3; U[17] = -twoq2;

    // Derivative with respect to q3.
    U[18] = -twoq3; U[21] =  twoq4; U[24] = twoq1;
    U[19] = -twoq4; U[22] = -twoq3; U[25] = twoq2;
    U[20] =  twoq1; U[23] =  twoq2; U[26] = twoq3;

    // Derivative with respect to q4.
    U[27] =  twoq4; U[30] =  twoq3; U[33] = -twoq2;
    U[28] = -twoq3; U[31] =  twoq4; U[34] =  twoq1;
    U[29] =  twoq2; U[32] = -twoq1; U[35] =  twoq4;
}

// Compute associated Legendre polynomials and their derivatives.
static void legendre(mwSize N, mwSize M, double phi, double* P, double* dP, double* ddP, double* dddP) {
    // Precompute trigonometric terms.
    mwSize rows = N + 1;

    double sinp = sin(phi);
    double cosp = cos(phi);

    // Base terms.
    P[0] = 1.0;
    dP[0] = 0.0;
    ddP[0] = 0.0;
    if (dddP) dddP[0] = 0.0;

    if (N) {
        P[IDX(1, 1, rows)] = sqrt(3) * cosp;
        dP[IDX(1, 1, rows)] = -sqrt(3) * sinp;
        ddP[IDX(1, 1, rows)] = -sqrt(3) * cosp;
        if (dddP) dddP[IDX(1, 1, rows)] = sqrt(3) * sinp;
    }

    // Diagonal terms.
    for (size_t i = 2; i <= N; i++) {
        double factor = sqrt((2.0 * i + 1.0) / (2.0 * i));
        size_t idPrev = IDX(i - 1, i - 1, rows);
        size_t idCur = IDX(i, i, rows);

        P[idCur] = factor * cosp * P[idPrev];
        dP[idCur] = factor * (cosp * dP[idPrev] - sinp * P[idPrev]);
        ddP[idCur] = factor * (cosp * (ddP[idPrev] - P[idPrev]) - 2.0 * sinp * dP[idPrev]);
        if (dddP) dddP[idCur] = factor * (-sinp * (ddP[idPrev] - P[idPrev]) + cosp * (dddP[idPrev] - dP[idPrev]) - 2 * cosp * dP[idPrev] - 2 * sinp * ddP[idPrev]);
    }

    // First off-diagonal terms.
    for (size_t i = 1; i <= N; i++) {
        double factor = sqrt(2.0 * i + 1.0);
        size_t idPrev = IDX(i - 1, i - 1, rows);
        size_t idCur = IDX(i, i - 1, rows);

        P[idCur] = factor * sinp * P[idPrev];
        dP[idCur] = factor * (cosp * P[idPrev] + sinp * dP[idPrev]);
        ddP[idCur] = factor * (sinp * (ddP[idPrev] - P[idPrev]) + 2.0 * cosp * dP[idPrev]);
        if (dddP) dddP[idCur] = factor * (cosp * (ddP[idPrev] - P[idPrev]) + sinp * (dddP[idPrev] - dP[idPrev]) - 2 * sinp * dP[idPrev] + 2 * cosp * ddP[idPrev]);
    }

    // Remaining off-diagonal recurrence terms.
    size_t j = 0, k = 2;
    while (1) {
        for (size_t i = k; i <= N; i++) {
            double f1 = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            double f2 = sqrt(2 * i - 1);
            double f3 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));
            P[IDX(i, j, rows)] = f1 * (f2 * sinp * P[IDX(i - 1, j, rows)] - 
                f3 * P[IDX(i - 2, j, rows)]);
        }
        j++;
        k++;
        if (j > M) break;
    }
    j = 0, k = 2;
    while (1) {
        for (size_t i = k; i <= N; i++) {
            double f1 = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            double f2 = sqrt(2 * i - 1);
            double f3 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));
            dP[IDX(i, j, rows)] = f1 * (f2 * (sinp * dP[IDX(i - 1, j, rows)] + 
                cosp * P[IDX(i - 1, j, rows)]) - f3 * dP[IDX(i - 2, j, rows)]);
        }
        j++;
        k++;
        if (j > M) break;
    }
    j = 0, k = 2;
    while (1) {
        for (size_t i = k; i <= N; ++i) {
            double f1 = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            double f2 = sqrt(2 * i - 1);
            double f3 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));
            ddP[IDX(i, j, rows)] = f1 * (f2 * (sinp * (ddP[IDX(i - 1, j, rows)] - 
                P[IDX(i - 1, j, rows)]) + 2 * cosp * dP[IDX(i - 1, j, rows)]) - 
                f3 * ddP[IDX(i - 2, j, rows)]);
        }
        j++;
        k++;
        if (j > M) break;
    }
    if (dddP) {
        j = 0, k = 2;
        while (1) {
            for (size_t i = k; i <= N; ++i) {
                double f1 = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
                double f2 = sqrt(2 * i - 1);
                double f3 = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0));
                dddP[IDX(i, j, rows)] = f1 * (f2 * (sinp * (dddP[IDX(i - 1, j, rows)] - 3 * dP[IDX(i - 1, j, rows)]) +
                                                    cosp * (3 * ddP[IDX(i - 1, j, rows)] - P[IDX(i - 1, j, rows)])) -
                                              f3 * dddP[IDX(i - 2, j, rows)]);
            }
            j++;
            k++;
            if (j > M) break;
        }
    }
}

// Compute derivatives of the Cartesian-to-spherical coordinate mapping.
static void coordinateDerivatives(const double x, const double y, const double z,
                                  const double r, const double rho, mxLogical doJac, mxLogical doHess,
                                  double* rGrad, double* rHess, double* rT3,
                                  double* phiGrad, double* phiHess, double* phiT3,
                                  double* lamGrad, double* lamHess, double* lamT3) {
    // Precompute repeated terms.
    double r2 = r * r, rho2 = rho * rho;
    double r3 = r2 * r, rho3 = rho2 * rho;
    double r4 = r3 * r, rho4 = rho3 * rho;
    double r5 = r4 * r, rho5 = rho4 * rho;
    double r6 = r5 * r, rho6 = rho5 * rho;
    double x2 = x * x, y2 = y * y, z2 = z * z;
    double x3 = x2 * x, y3 = y2 * y, z3 = z2 * z;
    double x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;

    double rinv = 1.0 / r, rhoinv = 1.0 / rho;
    double r2inv = 1.0 / r2, rho2inv = 1.0 / rho2;
    double r3inv = 1.0 / r3, rho3inv = 1.0 / rho3;
    double r4inv = 1.0 / r4, rho4inv = 1.0 / rho4;
    double r5inv = 1.0 / r5, rho5inv = 1.0 / rho5;
    double r6inv = 1.0 / r6, rho6inv = 1.0 / rho6;

    // Derivatives of r.
    rGrad[0] = x * rinv;
    rGrad[1] = y * rinv;
    rGrad[2] = z * rinv;

    if (doJac || doHess) {
        rHess[0] = (y2 + z2) * r3inv;
        rHess[1] = -x * y * r3inv;
        rHess[2] = -x * z * r3inv;
        rHess[3] = rHess[1];
        rHess[4] = (x2 + z2) * r3inv;
        rHess[5] = -y * z * r3inv;
        rHess[6] = rHess[2];
        rHess[7] = rHess[5];
        rHess[8] = (x2 + y2) * r3inv;
    }

    if (doHess) {
        double X[3] = {x, y, z};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    rT3[IDX3(i, j, k, 3, 3)] = -((double)(i == j) * X[k] + (double)(i == k) * X[j] + (double)(j == k) * X[i]) * r3inv +
                        3 * X[i] * X[j] * X[k] * r5inv;
                }
            }
        }
    }

    // Derivatives of rho.
    double rho_i[3] = {x * rhoinv, y * rhoinv, 0};
    double rho_ij[9] = { 0 };
    rho_ij[0] = y2 * rho3inv;
    rho_ij[3] = -x * y * rho3inv;
    rho_ij[1] = rho_ij[3];
    rho_ij[4] = x2 * rho3inv;

    double rho_ijk[27] = { 0 };
    double a = -3 * x * y2 * rho5inv;
    double b = 2 * y * rho3inv - 3 * y3 * rho5inv;
    double c = -x * rho3inv + 3 * x * y2 * rho5inv;
    double d = -3 * x2 * y * rho5inv;
    rho_ijk[0] = a;
    rho_ijk[9] = b;
    rho_ijk[3] = b;
    rho_ijk[1] = b;
    rho_ijk[12] = c;
    rho_ijk[10] = c;
    rho_ijk[4] = c;
    rho_ijk[13] = d;

    // Derivatives of longitude lambda.
    lamGrad[0] = -y * rho2inv;
    lamGrad[1] = x * rho2inv;
    lamGrad[2] = 0.0;

    if (doJac || doHess) {
        lamHess[0] = 2 * x * y * rho4inv;
        lamHess[1] = (y2 - x2) * rho4inv;
        lamHess[3] = lamHess[1];
        lamHess[4] = -2 * x * y * rho4inv;
    }

    if (doHess) {
        lamT3[0] = 2 * y * rho4inv - 8 * x2 * y * rho6inv;
        lamT3[9] = 2 * x * rho4inv - 8 * x * y2 * rho6inv;
        lamT3[3] = lamT3[9];
        lamT3[1] = lamT3[9];
        lamT3[12] = 2 * y * rho4inv - 4 * y * (y2 - x2) * rho6inv;
        lamT3[10] = lamT3[12];
        lamT3[4] = lamT3[12];
        lamT3[13] = -lamT3[9];
    }

    // Derivatives of latitude phi.
    double Fu, Fv, Fuu, Fuv, Fvv, Fuuu, Fuuv, Fuvv, Fvvv;
    Fu = rho * r2inv;
    Fv = -z * r2inv;
    Fuu = -2 * z * rho * r4inv;
    Fuv = (z2 - rho2) * r4inv;
    Fvv = -Fuu;

    Fuuu = -2 * rho * r4inv + 8 * z2 * rho * r6inv;
    Fuuv = -2 * z * r4inv + 8 * z * rho2 * r6inv;
    Fuvv = -2 * rho * r4inv - 4 * rho * (z2 - rho2) * r6inv;
    Fvvv = -Fuuv;

    phiGrad[0] = Fv * rho_i[0];
    phiGrad[1] = Fv * rho_i[1];
    phiGrad[2] = Fu + Fv * rho_i[2];

    if (doJac || doHess) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double a = (double)(i == 2), b = (double)(j == 2);
                phiHess[IDX(i, j, 3)] = Fv * rho_ij[IDX(i, j, 3)] +
                    Fvv * rho_i[i] * rho_i[j] + Fuv * (a * rho_i[j] + b * rho_i[i]) + Fuu * a * b;
            }
        }
    }

    if (doHess) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    int ijk = IDX3(i, j, k, 3, 3), ij = IDX(i, j, 3), ik = IDX(i, k, 3), jk = IDX(j, k, 3);
                    double i3 = (double)(i == 2), j3 = (double)(j == 2), k3 = (double)(k == 2);
                    phiT3[ijk] = Fv * rho_ijk[ijk] +
                        Fvv * (rho_ij[ij] * rho_i[k] + rho_ij[ik] * rho_i[j] + rho_ij[jk] * rho_i[i]) +
                        Fuv * (i3 * rho_ij[jk] + j3 * rho_ij[ik] + k3 * rho_ij[ij]) +
                        Fvvv * rho_i[i] * rho_i[j] * rho_i[k] +
                        Fuvv * (i3 * rho_i[j] * rho_i[k] + j3 * rho_i[i] * rho_i[k] + k3 * rho_i[i] * rho_i[j]) +
                        Fuuv * (i3 * j3 * rho_i[k] + i3 * k3 * rho_i[j] + j3 * k3 * rho_i[i]) +
                        Fuuu * i3 * j3 * k3;
                }
            }
        }
    }
}

// Compute spherical-harmonic gravity acceleration and derivatives.
static void gravityfield(const double* posInertial, mwSize maxDegree, mwSize maxOrder,
                         const double* tm, double gm, double req, const double* C, 
                         const double* S, mwSize coefRows, mxLogical doJac, mxLogical doHess,
                         double* accelOut, double* jacOut, double* hessOut) {
    // Transform to the body-fixed frame.
    double posBody[3] = { 0, 0, 0 };
    for (mwSize i = 0; i < 3; ++i) {
        for (mwSize j = 0; j < 3; ++j) {
            posBody[i] += tm[IDX(i, j, 3)] * posInertial[j];
        }
    }

    double x = posBody[0], y = posBody[1], z = posBody[2];

    // Spherical-coordinate terms.
    double r = sqrt(x * x + y * y + z * z);
    double phi = asin(z / r);
    double lam = atan2(y, x);

    double rho = sqrt(x * x + y * y);

    // Compute the Legendre arrays.
    mwSize dimP = (maxDegree + 1) * (maxDegree + 1);
    double* P = (double*)mxCalloc(dimP, sizeof(double));
    double* dP = (double*)mxCalloc(dimP, sizeof(double));
    double* ddP = (double*)mxCalloc(dimP, sizeof(double));
    double* dddP = (double*)mxCalloc(dimP, sizeof(double));

    legendre(maxDegree, maxOrder, phi, P, dP, ddP, dddP);

    // Accumulate the potential derivatives.
    double dUdr = 0, dUdp = 0, dUdl = 0;
    double d2U[6] = { 0 };
    double d3U[10] = { 0 };

    for (mwSize n = 0; n <= maxDegree;  ++n) {
        // Radial factor and its derivatives.
        double f0 = (gm / r) * pow(req / r, (double)n);
        double f1 = -((double)n + 1) / r * f0;
        double f2 = ((double)n + 1) * ((double)n + 2) / (r * r) * f0;
        double f3 = -((double)n + 1) * ((double)n + 2) * ((double)n + 3) / (r * r * r) * f0;

        // Angular sums.
        double q = 0, q_phi = 0, q_lam = 0;
        double q_phi2 = 0, q_philam = 0, q_lam2 = 0;
        double q_phi3 = 0, q_phi2lam = 0, q_philam2 = 0, q_lam3 = 0;

        for (mwSize m = 0; m <= maxOrder; ++m) {
            double cosML = cos(m * lam);
            double sinML = sin(m * lam);
            size_t idx = IDX(n, m, maxDegree + 1);
            size_t idxCS = IDX(n, m, coefRows);

            double pCpS = C[idxCS] * cosML + S[idxCS] * sinML;
            double pSnC = S[idxCS] * cosML - C[idxCS] * sinML;

            q += P[idx] * pCpS;
            q_phi += dP[idx] * pCpS;
            q_lam += (double)m * P[idx] * pSnC;

            if (doJac || doHess) {
                q_phi2 += ddP[idx] * pCpS;
                q_philam += (double)m * dP[idx] * pSnC;
                q_lam2 -= (double)(m * m) * P[idx] * pCpS;
            }

            if (doHess) {
                q_phi3 += dddP[idx] * pCpS;
                q_phi2lam += (double)m * ddP[idx] * pSnC;
                q_philam2 -= (double)(m * m) * dP[idx] * pCpS;
                q_lam3 -= (double)(m * m * m) * P[idx] * pSnC;
            }
        }

        // First derivatives.
        dUdr += f1 * q;
        dUdp += f0 * q_phi;
        dUdl += f0 * q_lam;

        if (doJac || doHess) {
            d2U[0] += f2 * q;           // rr 
            d2U[1] += f1 * q_phi;       // rφ
            d2U[2] += f1 * q_lam;       // rλ
            d2U[3] += f0 * q_phi2;      // φφ
            d2U[4] += f0 * q_philam;    // φλ
            d2U[5] += f0 * q_lam2;      // λλ
        }

        // Third derivatives.
        if (doHess) {
            d3U[0] += f3 * q;           // rrr
            d3U[1] += f2 * q_phi;       // rrφ
            d3U[2] += f2 * q_lam;       // rrλ
            d3U[3] += f1 * q_phi2;      // rφφ
            d3U[4] += f1 * q_philam;    // rφλ
            d3U[5] += f1 * q_lam2;      // rλλ
            d3U[6] += f0 * q_phi3;      // φφφ
            d3U[7] += f0 * q_phi2lam;   // φφλ
            d3U[8] += f0 * q_philam2;   // φλλ
            d3U[9] += f0 * q_lam3;      // λλλ
        }
    }

    // Derivatives of the coordinate mapping.
    double rGrad[3] = { 0 }, rHess[9] = { 0 }, rT3[27] = { 0 };
    double phiGrad[3] = { 0 }, phiHess[9] = { 0 }, phiT3[27] = { 0 };
    double lamGrad[3] = { 0 }, lamHess[9] = { 0 }, lamT3[27] = { 0 };
    coordinateDerivatives(x, y, z, r, rho, doJac, doHess, rGrad, rHess, rT3, phiGrad, phiHess, phiT3, lamGrad, lamHess, lamT3);

    // Compute acceleration.
    double accelBody[3];
    for (size_t i = 0; i < 3; ++i) {
        accelBody[i] = dUdr * rGrad[i] + dUdp * phiGrad[i] + dUdl * lamGrad[i];
    }
    for (size_t i = 0; i < 3; ++i) {
        accelOut[i] = 0.0;
        for (size_t j = 0; j < 3; ++j) {
            accelOut[i] += tm[IDX(j, i, 3)] * accelBody[j];
        }
    }
    if (!doJac && !doHess) {
        mxFree(P);
        mxFree(dP);
        mxFree(ddP);
        mxFree(dddP);
        return;
    }

    // Compute the Jacobian.
    if (doJac || doHess) {
        double jacBody[9] = { 0 };
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                size_t idx = IDX(i, j, 3);
                jacBody[idx] =
                    d2U[0] * rGrad[i] * rGrad[j] +
                    d2U[3] * phiGrad[i] * phiGrad[j] +
                    d2U[5] * lamGrad[i] * lamGrad[j] +
                    d2U[1] * (rGrad[i] * phiGrad[j] + rGrad[j] * phiGrad[i]) +
                    d2U[2] * (rGrad[i] * lamGrad[j] + rGrad[j] * lamGrad[i]) +
                    d2U[4] * (phiGrad[i] * lamGrad[j] + phiGrad[j] * lamGrad[i]) +
                    dUdr * rHess[idx] +
                    dUdp * phiHess[idx] +
                    dUdl * lamHess[idx];
            }
        }

        double tmp[9] = { 0 };
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    tmp[IDX(i, j, 3)] += tm[IDX(k, i, 3)] * jacBody[IDX(k, j, 3)];
                }
            }
        }
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                double val = 0;
                for (size_t k = 0; k < 3; ++k) {
                    val += tmp[IDX(i, k, 3)] * tm[IDX(k, j, 3)];
                }
                jacOut[IDX(i, j, 3)] = val;
                jacOut[IDX(i, j + 3, 3)] = 0.0;
            }
        }

        if (doHess) {
            // Derivatives of the potential.
            double U1[3] = {dUdr, dUdp, dUdl};
            double U2[9] = {d2U[0], d2U[1], d2U[2], d2U[1], d2U[3], d2U[4], d2U[2], d2U[4], d2U[5]};
            double U3[27] = { 0 };
            PUT3(0, 0, 0, d3U[0], U3);
            PUT3(0, 0, 1, d3U[1], U3);
            PUT3(0, 0, 2, d3U[2], U3);
            PUT3(0, 1, 1, d3U[3], U3);
            PUT3(0, 1, 2, d3U[4], U3);
            PUT3(0, 2, 2, d3U[5], U3);
            PUT3(1, 1, 1, d3U[6], U3);
            PUT3(1, 1, 2, d3U[7], U3);
            PUT3(1, 2, 2, d3U[8], U3);
            PUT3(2, 2, 2, d3U[9], U3);

            // Coordinate-derivative tensors.
            double* q1[3] = {rGrad, phiGrad, lamGrad};
            double* q2[3] = {rHess, phiHess, lamHess};
            double* q3[3] = {rT3, phiT3, lamT3};
            double hessBody[27];

            // Third-order tensor in the body-fixed frame.
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        double s = 0.0;
                        for (int a = 0; a < 3; ++a) {
                            for (int b = 0; b < 3; ++b) {
                                for (int c = 0; c < 3; ++c) {
                                    s += U3[IDX3(a, b, c, 3, 3)] * q1[a][i] * q1[b][j] * q1[c][k];
                                }
                            }
                        }
                        for (int a = 0; a < 3; ++a) {
                            for (int b = 0; b < 3; ++b) {
                                s += U2[IDX(a, b, 3)] * (q2[a][IDX(i, j, 3)] * q1[b][k] + q2[a][IDX(i, k, 3)] * q1[b][j] + q2[a][IDX(j, k, 3)] * q1[b][i]);
                            }
                        }
                        for (int a = 0; a < 3; ++a) {
                            s += U1[a] * q3[a][IDX3(i, j, k, 3, 3)];
                        }
                        hessBody[IDX3(i, j, k, 3, 3)] = s;
                    }
                }
            }

            // Transform from the body-fixed frame to the inertial frame.
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        double s = 0.0;
                        for (int p1 = 0; p1 < 3; ++p1) {
                            for (int p2 = 0; p2 < 3; ++p2) {
                                for (int p3 = 0; p3 < 3; ++p3) {
                                    s += tm[IDX(p1, i, 3)] * tm[IDX(p2, j, 3)] * tm[IDX(p3, k, 3)] * hessBody[IDX3(p1, p2, p3, 3, 3)];
                                }
                            }
                        }
                        hessOut[IDX3(i, j, k, 3, 3)] = s;
                    }
                }
            }
        }
        mxFree(P);
        mxFree(dP);
        mxFree(ddP);
        mxFree(dddP);
    }
}

// Compute gravity-field jerk.
static void jerkGravityfield(const double* posInertial, const double* velInertial,
                             const double* accelInertial, const double* omegaInertial,
                             const double* jacGravityfield, double* jerkOut) {
    // Auxiliary terms.
    double wxr[3], vplus[3], wxa[3];
    cross(omegaInertial, posInertial, wxr);
    cross(omegaInertial, accelInertial, wxa);

    // Jerk.
    for (int i = 0; i < 3; i++) {
        jerkOut[i] = -wxa[i];
        for (int j = 0; j < 3; j++) {
            jerkOut[i] += jacGravityfield[IDX(i, j, 3)] * (velInertial[j] - wxr[j]);
        }
    }
}

// Compute the Jacobian of one third-body gravitational acceleration term.
static void jacobiForthbody(const double* posInertial, double gm, double* jacOut) {
    double x = posInertial[0], y = posInertial[1], z = posInertial[2];
    double r = sqrt(x * x + y * y + z * z);
    double r3 = r * r * r, r5 = r3 * r * r;
    double threeGM = 3.0 * gm;

    jacOut[0] = -gm / r3 + threeGM * x * x / r5;
    jacOut[1] = threeGM * x * y / r5;
    jacOut[2] = threeGM * x * z / r5;

    jacOut[3] = jacOut[1];
    jacOut[4] = -gm / r3 + threeGM * y * y / r5;
    jacOut[5] = threeGM * y * z / r5;

    jacOut[6] = jacOut[2];
    jacOut[7] = jacOut[5];
    jacOut[8] = -gm / r3 + threeGM * z * z / r5;
}

// Compute third-body gravitational acceleration.
static void fourthbody(const double epoch, const double* posInertial, const mxLogical* doBodies, mxLogical doJac, 
                      double* accelOut, double* jacOut) {
    for (int i = 0; i < 9; ++i) {
        if (doBodies[i]) {
            double gm = GM[i];

            // Obtain the body position from SPICE.
            double posPert[3], lt;
            spkpos_c(BODYNAMES[i], epoch, "J2000", "NONE", "MOON", posPert, &lt);
            
            // Compute acceleration.
            double rVec[3] = { posInertial[0] - posPert[0], posInertial[1] - posPert[1], posInertial[2] - posPert[2] };
            double rNorm = sqrt(rVec[0] * rVec[0] + rVec[1] * rVec[1] + rVec[2] * rVec[2]);
            double rNorm3 = rNorm * rNorm * rNorm;

            double pertNorm = sqrt(posPert[0] * posPert[0] + posPert[1] * posPert[1] + posPert[2] * posPert[2]);
            double pertNorm3 = pertNorm * pertNorm * pertNorm;

            for (int j = 0; j < 3; ++j) {
                accelOut[j] -= gm * (rVec[j] / rNorm3 + posPert[j] / pertNorm3);
            }
            
            // Compute the Jacobian.
            if (doJac) {
                double jac[9];
                jacobiForthbody(rVec, gm, jac);
                for (int j = 0; j < 9; ++j) {
                    jacOut[j] += jac[j];
                }
            }
        }
    }
}

// Compute third-body gravitational jerk.
static void jerkFourthbody(const double epoch, const double* posInertial, const double* velInertial,
                          const mxLogical* doBodies, const mxLogical doJac, double* jerkOut, double* jacOut) {
    // Initialize outputs to zero.
    jerkOut[0] = 0.0; jerkOut[1] = 0.0; jerkOut[2] = 0.0;

    for (int i = 0; i < 9; i++) {
        if (doBodies[i]) {
            double gm = GM[i];

            // Obtain the body position and velocity from SPICE.
            double stPert[6], lt;
            spkezr_c(BODYNAMES[i], epoch, "J2000", "NONE", "MOON", stPert, &lt);
            const double* rBody = &stPert[0];
            const double* vBody = &stPert[3];

            // Compute jerk.
            double rVec[3] = { posInertial[0] - rBody[0], posInertial[1] - rBody[1], posInertial[2] - rBody[2] };
            double vVec[3] = { velInertial[0] - vBody[0], velInertial[1] - vBody[1], velInertial[2] - vBody[2] };
            double rNorm = sqrt(rVec[0] * rVec[0] + rVec[1] * rVec[1] + rVec[2] * rVec[2]);
            double bodyNorm = sqrt(rBody[0] * rBody[0] + rBody[1] * rBody[1] + rBody[2] * rBody[2]);

            double rNorm3inv = 1.0 / (rNorm * rNorm * rNorm), rNorm5inv = rNorm3inv / (rNorm * rNorm);
            double bodyNorm3inv = 1.0 / (bodyNorm * bodyNorm * bodyNorm), bodyNorm5inv = bodyNorm3inv / (bodyNorm * bodyNorm);

            double term1 = rVec[0] * vVec[0] + rVec[1] * vVec[1] + rVec[2] * vVec[2];
            double term2 = rBody[0] * vBody[0] + rBody[1] * vBody[1] + rBody[2] * vBody[2];

            for (int j = 0; j < 3; ++j) {
                double j1 = vVec[j] * rNorm3inv - 3.0 * term1 * rVec[j] * rNorm5inv;
                double j2 = vBody[j] * bodyNorm3inv - 3.0 * term2 * rBody[j] * bodyNorm5inv;
                jerkOut[j] -= gm * (j1 + j2);
            }

            // Jacobian.
            if (doJac) {
                double rNorm7inv = rNorm5inv / (rNorm * rNorm);
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        double Jr = gm * (3 * rNorm5inv * (rVec[j] * vVec[k] + rVec[k] * vVec[j] + ((j == k) ? term1 : 0.0)) -
                            15 * term1 * rNorm7inv * (rVec[j] * rVec[k]));
                        double Jv = -gm * (((j == k) ? rNorm3inv : 0.0) - 3 * rNorm5inv * rVec[j] * rVec[k]);
                        jacOut[IDX(j, k, 3)] += Jr;
                        jacOut[IDX(j, k + 3, 3)] += Jv;
                    }
                }
            }
        }
    }
}

// Compute the relativistic perturbation acceleration.
static void relativity(const double* posInertial, const double* velInertial, mxLogical doJac, double* accelOut, double* jacOut) {
    double r = sqrt(posInertial[0] * posInertial[0] + posInertial[1] * posInertial[1] + posInertial[2] * posInertial[2]);
    double s = sqrt(velInertial[0] * velInertial[0] + velInertial[1] * velInertial[1] + velInertial[2] * velInertial[2]);

    double r3 = r * r * r;
    double ls2 = LIGHTSPEED * LIGHTSPEED;
    double coef = MOONGM / (ls2 * r3);
    double term1 = 4.0 * MOONGM / r - s * s;
    double term2 = 4.0 * (posInertial[0] * velInertial[0] + posInertial[1] * velInertial[1] + posInertial[2] * velInertial[2]);
    
    // Compute acceleration.
    for (int i = 0; i < 3; ++i) {
        accelOut[i] = coef * (term1 * posInertial[i] + term2 * velInertial[i]);
    }

    // Compute the Jacobian.
    if (doJac) {
        double x = posInertial[0], y = posInertial[1], z = posInertial[2];
        double vx = velInertial[0], vy = velInertial[1], vz = velInertial[2];

        double c1 = coef;
        double c2 = 3 * MOONGM / (ls2 * r3 * r * r);
        double c3 = 4 * MOONGM / r3;
        double c4 = s * s - (4 * MOONGM) / r;
        double c5 = term2;

        double d1 = -c2 * (vx * c5 - x * c4);
        double d2 = -c2 * (vy * c5 - y * c4);
        double d3 = -c2 * (vz * c5 - z * c4);

        jacOut[0] = -c1 * (c4 - 4 * vx * vx + c3 * x * x) + x * d1;
        jacOut[1] = c1 * (4 * vx * vy - c3 * x * y) + x * d2;
        jacOut[2] = c1 * (4 * vx * vz - c3 * x * z) + x * d3;

        jacOut[3] = c1 * (4 * vx * vy - c3 * x * y) + y * d1;
        jacOut[4] = -c1 * (c4 - 4 * vy * vy + c3 * y * y) + y * d2;
        jacOut[5] = c1 * (4 * vy * vz - c3 * y * z) + y * d3;

        jacOut[6] = c1 * (4 * vx * vz - c3 * x * z) + z * d1;
        jacOut[7] = c1 * (4 * vy * vz - c3 * y * z) + z * d2;
        jacOut[8] = -c1 * (c4 - 4 * vz * vz + c3 * z * z) + z * d3;

        jacOut[9] = c1 * (2 * vx * x + c5);
        jacOut[10] = c1 * (4 * vy * x - 2 * vx * y);
        jacOut[11] = c1 * (4 * vz * x - 2 * vx * z);

        jacOut[12] = -c1 * (2 * vy * x - 4 * vx * y);
        jacOut[13] = c1 * (2 * vy * y + c5);
        jacOut[14] = c1 * (4 * vz * y - 2 * vy * z);

        jacOut[15] = -c1 * (2 * vz * x - 4 * vx * z);
        jacOut[16] = -c1 * (2 * vz * y - 4 * vy * z);
        jacOut[17] = c1 * (2 * vz * z + c5);
    }
}

// Compute the relativistic perturbation jerk.
static void jerkRelativity(const double* posInertial, const double* velInertial, const double* accelInertial,
                           const double* jacRelativity, double* jerkOut) {
    for (int i = 0; i < 3; i++) {
        double tmp = 0.0;
        tmp += jacRelativity[IDX(i, 0, 3)] * velInertial[0];
        tmp += jacRelativity[IDX(i, 1, 3)] * velInertial[1];
        tmp += jacRelativity[IDX(i, 2, 3)] * velInertial[2];
        tmp += jacRelativity[IDX(i, 3, 3)] * accelInertial[0];
        tmp += jacRelativity[IDX(i, 4, 3)] * accelInertial[1];
        tmp += jacRelativity[IDX(i, 5, 3)] * accelInertial[2];
        jerkOut[i] = tmp;
    }
}

// Compute solar-radiation-pressure acceleration using the cannonball model.
static void srp(const double* posInertial, const double* posSun, const double rpcm, mxLogical doJac,
                double* accelOut, double* jacOut) {
    double u[3], r2 = 0, rL2 = 0, dL2 = 0;

    for (int i = 0; i < 3; ++i) {
        u[i] = posSun[i] - posInertial[i];
        dL2 += u[i] * u[i];
    }
    double dL = sqrt(dL2);
    
    // Compute acceleration.
    double P = SOLARLUMINOSITY / (4 * PI * LIGHTSPEED * dL2);
    double coef1 = (P * rpcm * 1e-12) / dL;
    for (int i = 0; i < 3; ++i) {
        accelOut[i] = -coef1 * u[i];
    }

    // Compute the Jacobian.
    if (doJac) {
        double coef2 = SOLARLUMINOSITY * rpcm * 1e-12 / (4 * PI * LIGHTSPEED);
        double x = posInertial[0] - posSun[0], y = posInertial[1] - posSun[1], z = posInertial[2] - posSun[2];
        double dL3 = dL2 * dL, dL5 = dL3 * dL2;

        jacOut[0] = coef2 * (1.0 / dL3 - 3.0 * x * x / dL5);
        jacOut[1] = -3.0 * coef2 * x * y / dL5;
        jacOut[2] = -3.0 * coef2 * x * z / dL5;
        jacOut[3] = jacOut[1];
        jacOut[4] = coef2 * (1.0 / dL3 - 3.0 * y * y / dL5);
        jacOut[5] = -3.0 * coef2 * y * z / dL5;
        jacOut[6] = jacOut[2];
        jacOut[7] = jacOut[5];
        jacOut[8] = coef2 * (1.0 / dL3 - 3.0 * z * z / dL5);
    }
}

// Compute solar-radiation-pressure jerk using the cannonball model.
static void jerkSrp(const double* posInertial, const double* velInertial, const double* posSun,
                    const double* velSun, const double rpcm, const mxLogical doJac, double* jerkOut, double* jacOut) {
    double u[3], r2 = 0, rL2 = 0, dL2 = 0;

    for (int i = 0; i < 3; ++i) {
        u[i] = posSun[i] - posInertial[i];
        r2 += posSun[i] * posSun[i];
        rL2 += posInertial[i] * posInertial[i];
        dL2 += u[i] * u[i];
    }

    double dL = sqrt(dL2);

    double r[3] = { posInertial[0] - posSun[0], posInertial[1] - posSun[1], posInertial[2] - posSun[2] };
    double v[3] = { velInertial[0] - velSun[0], velInertial[1] - velSun[1], velInertial[2] - velSun[2] };
    double rNorm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    double rNorm3inv = 1.0 / (rNorm * rNorm * rNorm), rNorm5inv = rNorm3inv / (rNorm * rNorm);
    double dotTerm = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
    double K = SOLARLUMINOSITY * rpcm * 1e-12 / (4.0 * PI * LIGHTSPEED);
    double coef = -3.0 * rNorm5inv * dotTerm;

    // Compute jerk.
    jerkOut[0] = K * (v[0] * rNorm3inv + coef * r[0]);
    jerkOut[1] = K * (v[1] * rNorm3inv + coef * r[1]);
    jerkOut[2] = K * (v[2] * rNorm3inv + coef * r[2]);

    if (doJac) {
        double rNorm7inv = rNorm5inv / (rNorm * rNorm);
        double coef2 = -3.0 * dotTerm * rNorm5inv;

        // Compute the jerk Jacobian.
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double g = -3.0 * r[j] * rNorm5inv;
                double w = -3.0 * v[j] * rNorm5inv + 15.0 * dotTerm * r[j] * rNorm7inv;
                double Jr = g * v[i] + w * r[i] + ((i == j) ? coef2 : 0.0);
                double Jv = rNorm3inv * ((i == j) ? 1.0 : 0.0) - 3.0 * rNorm5inv * r[i] * r[j];
                jacOut[IDX(i, j, 3)] = K * Jr;
                jacOut[IDX(i, j + 3, 3)] = K * Jv;
            }
        }
    }
}

// Compute N-plate solar-radiation-pressure force using a quaternion attitude state.
static void srpnplate_q(const double* posInertial, const double* posSun,
                        const double* velInertial, const double* velSun,
                        const double* omegaBody, const double* q, const mxArray* plates, 
                        const mxLogical doJac, const mxLogical doJerk, const mwSize numPlates,
                        double* forceOut, double* jacOut, double* jerkOut, double* jerkJacOut) {
    double tm[9], Uq[36], u[3], s[3], du[3], d2 = 0, d, dinv, sdu, Q;
    q2R(q, tm); // Inertial to body-fixed.

    if (doJac || doJerk) q2Rjac(q, Uq);

    for (int i = 0; i < 3; ++i) {
        u[i] = posSun[i] - posInertial[i];
        d2 += u[i] * u[i];
    }
    d = sqrt(d2);
    dinv = 1.0 / d;
    
    for (int i = 0; i < 3; ++i) {
        s[i] = u[i] * dinv;
    }

    double P = 1e-9 * SOLARLUMINOSITY / (4 * PI * LIGHTSPEED * d2);

    double dsdp[9], dPdp[3], dsdt[3], dPdt = 0;
    if (doJac) {
        dsdp[0] = (s[0] * s[0] - 1) * dinv;
        dsdp[1] = s[0] * s[1] * dinv;
        dsdp[2] = s[0] * s[2] * dinv;
        dsdp[3] = dsdp[1];
        dsdp[4] = (s[1] * s[1] - 1) * dinv;
        dsdp[5] = s[1] * s[2] * dinv;
        dsdp[6] = dsdp[2];
        dsdp[7] = dsdp[5];
        dsdp[8] = (s[2] * s[2] - 1) * dinv;

        dPdp[0] = 2 * P * s[0] * dinv;
        dPdp[1] = 2 * P * s[1] * dinv;
        dPdp[2] = 2 * P * s[2] * dinv;
    }
    if (doJerk) {
        for (int i = 0; i < 3; ++i) {
            du[i] = velSun[i] - velInertial[i];
        }
        dsdt[0] = -((s[0] * s[0] - 1) * du[0] + s[0] * s[1] * du[1] + s[0] * s[2] * du[2]) * dinv;
        dsdt[1] = -(s[0] * s[1] * du[0] + (s[1] * s[1] - 1) * du[1] + s[1] * s[2] * du[2]) * dinv;
        dsdt[2] = -(s[0] * s[2] * du[0] + s[1] * s[2] * du[1] + (s[2] * s[2] - 1) * du[2]) * dinv;

        sdu = s[0] * du[0] + s[1] * du[1] + s[2] * du[2];
        Q = -2.0 * P * dinv;

        dPdt = Q * sdu;
    }

    // Accumulate force plate by plate.
    for (int k = 0; k < numPlates; ++k) {
        // Plate normal in the inertial frame.
        double* nBody = mxGetDoubles(mxGetField(plates, k, "UnitNormal"));
        double n[3];
        n[0] = -(tm[0] * nBody[0] + tm[1] * nBody[1] + tm[2] * nBody[2]);
        n[1] = -(tm[3] * nBody[0] + tm[4] * nBody[1] + tm[5] * nBody[2]);
        n[2] = -(tm[6] * nBody[0] + tm[7] * nBody[1] + tm[8] * nBody[2]);

        // Optical properties.
        double A = mxGetScalar(mxGetField(plates, k, "ReferenceArea"));
        double sr = mxGetScalar(mxGetField(plates, k, "SpecularReflectivity"));
        double dr = mxGetScalar(mxGetField(plates, k, "DiffuseReflectivity"));
        double alpha = fmax(0.0, 1 - sr);

        // Compute force.
        double mu = -(n[0] * s[0] + n[1] * s[1] + n[2] * s[2]);
        mu = fmax(0.0, mu);

        double coef1 = P * A * mu * 2 * (dr / 3.0 + sr * mu);
        double coef2 = -P * A * mu * alpha;

        forceOut[IDX(k, 0, numPlates)] = coef1 * n[0] + coef2 * s[0];
        forceOut[IDX(k, 1, numPlates)] = coef1 * n[1] + coef2 * s[1];
        forceOut[IDX(k, 2, numPlates)] = coef1 * n[2] + coef2 * s[2];

        if (mu > 0.0) {
            double dmudp[3], dc1dp[3], dc2dp[3];
            // Compute the Jacobian.
            if (doJac) {
                // Jacobian with respect to position.
                for (int j = 0; j < 3; ++j) {
                    double sum = 0.0;
                    for (int i = 0; i < 3; ++i) {
                        sum += n[i] * dsdp[IDX(i, j, 3)];
                    }
                    dmudp[j] = -sum;
                }
                double g1 = 2 * mu * (dr / 3 + sr * mu), g2 = 2 * P * (dr / 3 + 2 * sr * mu);
                for (int i = 0; i < 3; ++i) {
                    dc1dp[i] = A * (g1 * dPdp[i] + g2 * dmudp[i]);
                    dc2dp[i] = -alpha * A * (mu * dPdp[i] + P * dmudp[i]);
                }
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        jacOut[IDX3(i, j, k, 3, 13)] = n[i] * dc1dp[j] + s[i] * dc2dp[j] + coef2 * dsdp[IDX(i, j, 3)];
                    }
                }

                // Jacobian with respect to the quaternion.
                for (int ell = 0; ell < 4; ++ell) {
                    double* Uell = Uq + 9 * ell;
                    double dndq[3];
                    dndq[0] = -(Uell[0] * nBody[0] + Uell[1] * nBody[1] + Uell[2] * nBody[2]);
                    dndq[1] = -(Uell[3] * nBody[0] + Uell[4] * nBody[1] + Uell[5] * nBody[2]);
                    dndq[2] = -(Uell[6] * nBody[0] + Uell[7] * nBody[1] + Uell[8] * nBody[2]);

                    double dmudq = -(dndq[0] * s[0] + dndq[1] * s[1] + dndq[2] * s[2]);

                    double dc1dq = 2.0 * P * A * (dr / 3.0 + 2.0 * sr * mu) * dmudq;
                    double dc2dq = -P * A * alpha * dmudq;

                    jacOut[IDX3(0, 9 + ell, k, 3, 13)] = dc1dq * n[0] + coef1 * dndq[0] + dc2dq * s[0];
                    jacOut[IDX3(1, 9 + ell, k, 3, 13)] = dc1dq * n[1] + coef1 * dndq[1] + dc2dq * s[1];
                    jacOut[IDX3(2, 9 + ell, k, 3, 13)] = dc1dq * n[2] + coef1 * dndq[2] + dc2dq * s[2];
                }
            }
            
            // Compute force rate.
            if (doJerk) {
                double dmudt, dc1dt, dc2dt, omegaInertial[3], ndot[3] = {0, 0, 0};
                omegaInertial[0] = tm[0] * omegaBody[0] + tm[1] * omegaBody[1] + tm[2] * omegaBody[2];
                omegaInertial[1] = tm[3] * omegaBody[0] + tm[4] * omegaBody[1] + tm[5] * omegaBody[2];
                omegaInertial[2] = tm[6] * omegaBody[0] + tm[7] * omegaBody[1] + tm[8] * omegaBody[2];
                cross(omegaInertial, n, ndot);
                dmudt = -(n[0] * dsdt[0] + n[1] * dsdt[1] + n[2] * dsdt[2] + ndot[0] * s[0] + ndot[1] * s[1] + ndot[2] * s[2]);
                dc1dt = 2 * A * ((mu * dr / 3 + sr * mu * mu) * dPdt + (P * dr / 3 + 2 * sr * P * mu) * dmudt);
                dc2dt = -alpha * A * (mu * dPdt + P * dmudt);

                jerkOut[IDX(k, 0, numPlates)] = n[0] * dc1dt + s[0] * dc2dt + coef1 * ndot[0] + coef2 * dsdt[0];
                jerkOut[IDX(k, 1, numPlates)] = n[1] * dc1dt + s[1] * dc2dt + coef1 * ndot[1] + coef2 * dsdt[1];
                jerkOut[IDX(k, 2, numPlates)] = n[2] * dc1dt + s[2] * dc2dt + coef1 * ndot[2] + coef2 * dsdt[2];

                // Compute the force-rate Jacobian.
                if (doJac) {
                    double d2inv = dinv * dinv;
                    double M[9];
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            M[IDX(i, j, 3)] = (((i == j) ? 1.0 : 0.0) - s[i] * s[j]) * dinv;
                        }
                    }

                    double B1 = 2.0 * A * mu * (dr / 3.0 + sr * mu);
                    double C1 = 2.0 * A * P * (dr / 3.0 + 2.0 * sr * mu);
                    double B2 = -alpha * A * mu;
                    double C2 = -alpha * A * P;

                    // Jacobian with respect to position.
                    for (int j = 0; j < 3; ++j) {
                        double dsdu = dsdp[IDX(0, j, 3)] * du[0] + dsdp[IDX(1, j, 3)] * du[1] + dsdp[IDX(2, j, 3)] * du[2];
                        double dPdtdp = -(2.0 * dinv) * dPdp[j] * sdu - (2.0 * P * dinv) * dsdu - (2.0 * P * d2inv) * sdu * s[j];

                        double dsdtdp[3];
                        for (int i = 0; i < 3; ++i) {
                            dsdtdp[i] = -(dsdp[IDX(i, j, 3)] * sdu + s[i] * dsdu) * dinv + (s[j] * dinv) * dsdt[i];
                        }

                        double ndsdtdp = n[0] * dsdtdp[0] + n[1] * dsdtdp[1] + n[2] * dsdtdp[2];
                        double ndotdsdp = ndot[0] * dsdp[IDX(0, j, 3)] + ndot[1] * dsdp[IDX(1, j, 3)] + ndot[2] * dsdp[IDX(2, j, 3)];
                        double dmudtdp = -(ndsdtdp + ndotdsdp);

                        double dB1dp = 2.0 * A * ((dr / 3.0) + 2.0 * sr * mu) * dmudp[j];
                        double dC1dp = 2.0 * A * (((dr / 3.0) + 2.0 * sr * mu) * dPdp[j] + 2.0 * sr * P * dmudp[j]);
                        double dB2dp = -alpha * A * dmudp[j];
                        double dC2dp = -alpha * A * dPdp[j];

                        double ddc1dtdp = dB1dp * dPdt + B1 * dPdtdp + dC1dp * dmudt + C1 * dmudtdp;
                        double ddc2dtdp = dB2dp * dPdt + B2 * dPdtdp + dC2dp * dmudt + C2 * dmudtdp;

                        for (int i = 0; i < 3; ++i) {
                            jerkJacOut[IDX3(i, j, k, 3, 13)] = n[i] * ddc1dtdp + s[i] * ddc2dtdp + dsdp[IDX(i, j, 3)] * dc2dt + 
                                dc1dp[j] * ndot[i] + coef2 * dsdtdp[i] + dc2dp[j] * dsdt[i];
                        }
                    }

                    // Jacobian with respect to velocity.
                    double tmp = 2.0 * P * dinv;
                    double nTM[3] = {
                        n[0] * M[0] + n[1] * M[1] + n[2] * M[2],
                        n[0] * M[3] + n[1] * M[4] + n[2] * M[5],
                        n[0] * M[6] + n[1] * M[7] + n[2] * M[8]
                    };
                    for (int j = 0; j < 3; ++j) {
                        double ddc1dtdv = B1 * tmp * s[j] + C1 * nTM[j];
                        double ddc2dtdv = B2 * tmp * s[j] + C2 * nTM[j];
                        for (int i = 0; i < 3; ++i) {
                            jerkJacOut[IDX3(i, j + 3, k, 3, 13)] = n[i] * ddc1dtdv + s[i] * ddc2dtdv - coef2 * M[IDX(i, j, 3)];
                        }
                    }

                    // Jacobian with respect to angular velocity.
                    double e1 = P * (dr / 3.0 + 2.0 * sr * mu);
                    for (int j = 0; j < 3; ++j) {
                        double dndotdwj[3] = {
                            tm[j + 3] * n[2] - tm[j + 6] * n[1],
                            tm[j + 6] * n[0] - tm[j] * n[2],
                            tm[j] * n[1] - tm[j + 3] * n[0]  
                        };
                        double dmudtdwj = -(dndotdwj[0] * s[0] + dndotdwj[1] * s[1] + dndotdwj[2] * s[2]);
                        double ddc1dtdwj = 2.0 * A * e1 * dmudtdwj;
                        double ddc2dtdwj = -alpha * A * P * dmudtdwj;
                        jerkJacOut[IDX3(0, j + 6, k, 3, 13)] = n[0] * ddc1dtdwj + s[0] * ddc2dtdwj + coef1 * dndotdwj[0];
                        jerkJacOut[IDX3(1, j + 6, k, 3, 13)] = n[1] * ddc1dtdwj + s[1] * ddc2dtdwj + coef1 * dndotdwj[1];
                        jerkJacOut[IDX3(2, j + 6, k, 3, 13)] = n[2] * ddc1dtdwj + s[2] * ddc2dtdwj + coef1 * dndotdwj[2];
                    }

                    // Jacobian with respect to the quaternion.
                    double c1mu = dr / 3.0 + 2.0 * sr * mu;
                    for (int ell = 0; ell < 4; ++ell) {
                        double* Uell = Uq + 9 * ell;

                        double dndq[3];
                        dndq[0] = -(Uell[0] * nBody[0] + Uell[1] * nBody[1] + Uell[2] * nBody[2]);
                        dndq[1] = -(Uell[3] * nBody[0] + Uell[4] * nBody[1] + Uell[5] * nBody[2]);
                        dndq[2] = -(Uell[6] * nBody[0] + Uell[7] * nBody[1] + Uell[8] * nBody[2]);

                        double dwdq[3];
                        dwdq[0] = Uell[0] * omegaBody[0] + Uell[1] * omegaBody[1] + Uell[2] * omegaBody[2];
                        dwdq[1] = Uell[3] * omegaBody[0] + Uell[4] * omegaBody[1] + Uell[5] * omegaBody[2];
                        dwdq[2] = Uell[6] * omegaBody[0] + Uell[7] * omegaBody[1] + Uell[8] * omegaBody[2];

                        double dndotdq[3];
                        dndotdq[0] = dwdq[1] * n[2] - dwdq[2] * n[1] + omegaInertial[1] * dndq[2] - omegaInertial[2] * dndq[1];
                        dndotdq[1] = dwdq[2] * n[0] - dwdq[0] * n[2] + omegaInertial[2] * dndq[0] - omegaInertial[0] * dndq[2];
                        dndotdq[2] = dwdq[0] * n[1] - dwdq[1] * n[0] + omegaInertial[0] * dndq[1] - omegaInertial[1] * dndq[0];

                        double dmudq = -(dndq[0] * s[0] + dndq[1] * s[1] + dndq[2] * s[2]);
                        double dmudtdq = -(dndq[0] * dsdt[0] + dndq[1] * dsdt[1] + dndq[2] * dsdt[2] + dndotdq[0] * s[0] + dndotdq[1] * s[1] + dndotdq[2] * s[2]);

                        double e1q = P * c1mu;
                        double ddc1dtdq = 2.0 * A * (c1mu * dmudq * dPdt + e1q * dmudtdq + 2.0 * sr * P * dmudq * dmudt);
                        double ddc2dtdq = -alpha * A * (dPdt * dmudq + P * dmudtdq);

                        double dc1dq = 2.0 * P * A * c1mu * dmudq;
                        double dc2dq = -P * A * alpha * dmudq;

                        jerkJacOut[IDX3(0, ell + 9, k, 3, 13)] = dndq[0] * dc1dt + n[0] * ddc1dtdq + s[0] * ddc2dtdq + dc1dq * ndot[0] + coef1 * dndotdq[0] + dc2dq * dsdt[0];
                        jerkJacOut[IDX3(1, ell + 9, k, 3, 13)] = dndq[1] * dc1dt + n[1] * ddc1dtdq + s[1] * ddc2dtdq + dc1dq * ndot[1] + coef1 * dndotdq[1] + dc2dq * dsdt[1];
                        jerkJacOut[IDX3(2, ell + 9, k, 3, 13)] = dndq[2] * dc1dt + n[2] * ddc1dtdq + s[2] * ddc2dtdq + dc1dq * ndot[2] + coef1 * dndotdq[2] + dc2dq * dsdt[2];
                    }
                }
            }
        }
    }
}

// Compute N-plate solar-radiation-pressure force using a rotation-matrix attitude state.
static void srpnplate_R(const double* posInertial, const double* posSun,
                        const double* velInertial, const double* velSun,
                        const double* tm, const mxArray* plates, 
                        const mxLogical doJac, const mwSize numPlates,
                        double* forceOut, double* jacOut) {
    double u[3], s[3], d2 = 0, d, dinv;

    for (int i = 0; i < 3; ++i) {
        u[i] = posSun[i] - posInertial[i];
        d2 += u[i] * u[i];
    }
    d = sqrt(d2);
    dinv = 1.0 / d;
    
    for (int i = 0; i < 3; ++i) {
        s[i] = u[i] * dinv;
    }

    double P = 1e-9 * SOLARLUMINOSITY / (4 * PI * LIGHTSPEED * d2);

    double dsdp[9], dPdp[3];
    if (doJac) {
        dsdp[0] = (s[0] * s[0] - 1) * dinv;
        dsdp[1] = s[0] * s[1] * dinv;
        dsdp[2] = s[0] * s[2] * dinv;
        dsdp[3] = dsdp[1];
        dsdp[4] = (s[1] * s[1] - 1) * dinv;
        dsdp[5] = s[1] * s[2] * dinv;
        dsdp[6] = dsdp[2];
        dsdp[7] = dsdp[5];
        dsdp[8] = (s[2] * s[2] - 1) * dinv;

        dPdp[0] = 2 * P * s[0] * dinv;
        dPdp[1] = 2 * P * s[1] * dinv;
        dPdp[2] = 2 * P * s[2] * dinv;
    }

    // Accumulate force plate by plate.
    for (int k = 0; k < numPlates; ++k) {
        // Plate normal in the inertial frame.
        double* nBody = mxGetDoubles(mxGetField(plates, k, "UnitNormal"));
        double n[3];
        n[0] = -(tm[0] * nBody[0] + tm[1] * nBody[1] + tm[2] * nBody[2]);
        n[1] = -(tm[3] * nBody[0] + tm[4] * nBody[1] + tm[5] * nBody[2]);
        n[2] = -(tm[6] * nBody[0] + tm[7] * nBody[1] + tm[8] * nBody[2]);

        // Optical properties.
        double A = mxGetScalar(mxGetField(plates, k, "ReferenceArea"));
        double sr = mxGetScalar(mxGetField(plates, k, "SpecularReflectivity"));
        double dr = mxGetScalar(mxGetField(plates, k, "DiffuseReflectivity"));
        double alpha = fmax(0.0, 1 - sr);

        // Compute force.
        double mu = -(n[0] * s[0] + n[1] * s[1] + n[2] * s[2]);
        mu = fmax(0.0, mu);

        double coef1 = P * A * mu * 2 * (dr / 3.0 + sr * mu);
        double coef2 = -P * A * mu * alpha;

        forceOut[IDX(k, 0, numPlates)] = coef1 * n[0] + coef2 * s[0];
        forceOut[IDX(k, 1, numPlates)] = coef1 * n[1] + coef2 * s[1];
        forceOut[IDX(k, 2, numPlates)] = coef1 * n[2] + coef2 * s[2];

        if (mu > 0.0) {
            double dmudp[3], dc1dp[3], dc2dp[3];
            // Compute the Jacobian.
            if (doJac) {
                // Jacobian with respect to position.
                for (int j = 0; j < 3; ++j) {
                    double sum = 0.0;
                    for (int i = 0; i < 3; ++i) {
                        sum += n[i] * dsdp[IDX(i, j, 3)];
                    }
                    dmudp[j] = -sum;
                }
                double g1 = 2 * mu * (dr / 3 + sr * mu), g2 = 2 * P * (dr / 3 + 2 * sr * mu);
                for (int i = 0; i < 3; ++i) {
                    dc1dp[i] = A * (g1 * dPdp[i] + g2 * dmudp[i]);
                    dc2dp[i] = -alpha * A * (mu * dPdp[i] + P * dmudp[i]);
                }
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        jacOut[IDX3(i, j, k, 3, 18)] = n[i] * dc1dp[j] + s[i] * dc2dp[j] + coef2 * dsdp[IDX(i, j, 3)];
                    }
                }

                // Jacobian with respect to the rotation matrix.
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        int idxR = 3 * b + a;

                        double dndR[3] = {0.0, 0.0, 0.0};
                        dndR[a] = -nBody[b];

                        double dmudR = -(dndR[0] * s[0] + dndR[1] * s[1] + dndR[2] * s[2]);
                        double dcoef1dR = 2.0 * P * A * (dr / 3.0 + 2.0 * sr * mu) * dmudR;
                        double dcoef2dR = -P * A * alpha * dmudR;
                        jacOut[IDX3(0, 9 + idxR, k, 3, 18)] = dcoef1dR * n[0] + coef1 * dndR[0] + dcoef2dR * s[0];
                        jacOut[IDX3(1, 9 + idxR, k, 3, 18)] = dcoef1dR * n[1] + coef1 * dndR[1] + dcoef2dR * s[1];
                        jacOut[IDX3(2, 9 + idxR, k, 3, 18)] = dcoef1dR * n[2] + coef1 * dndR[2] + dcoef2dR * s[2];
                    }
                }
            }
        }
    }
}

// Compute Earth-albedo perturbation acceleration.
static void earthalbedo(const double* posInertial, const double* posEarth, const double* posSun, const double rpcm, mxLogical doJac,
                        double* accelOut, double* jacOut) {
    double u[3], r2 = 0, rL2 = 0, dL2 = 0, rS2 = 0;

    for (int i = 0; i < 3; ++i) {
        u[i] = posEarth[i] - posInertial[i];
        rS2 += (posSun[i] - posEarth[i]) * (posSun[i] - posEarth[i]);
        r2 += posEarth[i] * posEarth[i];
        rL2 += posInertial[i] * posInertial[i];
        dL2 += u[i] * u[i];
    }

    // Evaluate the illumination condition.
    double rS = sqrt(rS2), r = sqrt(r2), rL = sqrt(rL2), d = sqrt(r2 - MOONRADIUS * MOONRADIUS), dL = sqrt(dL2);
    double cosEta = d / r, cosEtaL = (r2 + dL2 - rL2) / (2 * r * dL);
    mxLogical doAlb = (cosEtaL <= cosEta) || (dL <= d);

    if (doAlb) {
        // Compute acceleration.
        double PE = SOLARLUMINOSITY / (4 * PI * LIGHTSPEED * rS2);
        double P = 0.25 * PE * EARTHREFLECTION * (EARTHRADIUS * EARTHRADIUS) / dL2;
        double coef1 = (P * rpcm * 1e-12) / dL;
        for (int i = 0; i < 3; ++i) {
            accelOut[i] = -coef1 * u[i];
        }

        // Compute the Jacobian.
        if (doJac) {
            double coef2 = 0.25 * PE * EARTHREFLECTION * EARTHRADIUS * EARTHRADIUS * rpcm * 1e-12;
            double x = posInertial[0] - posEarth[0], y = posInertial[1] - posEarth[1], z = posInertial[2] - posEarth[2];
            double dL3 = dL2 * dL, dL5 = dL3 * dL2;

            jacOut[0] = coef2 * (1.0 / dL3 - 3.0 * x * x / dL5);
            jacOut[1] = -3.0 * coef2 * x * y / dL5;
            jacOut[2] = -3.0 * coef2 * x * z / dL5;
            jacOut[3] = jacOut[1];
            jacOut[4] = coef2 * (1.0 / dL3 - 3.0 * y * y / dL5);
            jacOut[5] = -3.0 * coef2 * y * z / dL5;
            jacOut[6] = jacOut[2];
            jacOut[7] = jacOut[5];
            jacOut[8] = coef2 * (1.0 / dL3 - 3.0 * z * z / dL5);
        }
    }
}

// Compute Earth-albedo perturbation jerk.
static void jerkEarthalbedo(const double* posInertial, const double* velInertial, const double* posEarth, const double* velEarth,
                            const double* posSun, const double* velSun, const double rpcm, double* jerkOut) {
    double u[3], r2 = 0, rL2 = 0, dL2 = 0, rS2 = 0;

    for (int i = 0; i < 3; ++i) {
        u[i] = posEarth[i] - posInertial[i];
        rS2 += (posSun[i] - posEarth[i]) * (posSun[i] - posEarth[i]);
        r2 += posEarth[i] * posEarth[i];
        rL2 += posInertial[i] * posInertial[i];
        dL2 += u[i] * u[i];
    }

    // Evaluate the illumination condition.
    double rS = sqrt(rS2), r = sqrt(r2), rL = sqrt(rL2), d = sqrt(r2 - MOONRADIUS * MOONRADIUS), dL = sqrt(dL2);
    double cosEta = d / r, cosEtaL = (r2 + dL2 - rL2) / (2 * r * dL);
    mxLogical doAlb = (cosEtaL <= cosEta) || (dL <= d);

    if (doAlb) {
        const double PE = SOLARLUMINOSITY / (4.0 * PI * LIGHTSPEED * rS2);
        const double coef2 = 0.25 * PE * EARTHREFLECTION * EARTHRADIUS * EARTHRADIUS * rpcm * 1e-12;

        double rSvec[3] = {posSun[0] - posEarth[0], posSun[1] - posEarth[1], posSun[2] - posEarth[2]};
        double vSvec[3] = {velSun[0] - velEarth[0], velSun[1] - velEarth[1], velSun[2] - velEarth[2]};
        double rSdotvS = rSvec[0] * vSvec[0] + rSvec[1] * vSvec[1] + rSvec[2] * vSvec[2];

        double coef2dot = 0.0;
        if (rS2 > 0.0) {
            coef2dot = -2.0 * coef2 * (rSdotvS / rS2);
        }

        double r[3] = { posInertial[0] - posEarth[0], posInertial[1] - posEarth[1], posInertial[2] - posEarth[2] };
        double v[3] = { velInertial[0] - velEarth[0], velInertial[1] - velEarth[1], velInertial[2] - velEarth[2] };

        const double eps = 1e-12;
        if (dL < eps)
            return;

        double dL3 = dL * dL2;
        double dL5 = dL3 * dL2;
        double xvx = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

        jerkOut[0] = coef2 * (v[0] / dL3 - 3.0 * xvx * r[0] / dL5) + coef2dot * (r[0] / dL3);
        jerkOut[1] = coef2 * (v[1] / dL3 - 3.0 * xvx * r[1] / dL5) + coef2dot * (r[1] / dL3);
        jerkOut[2] = coef2 * (v[2] / dL3 - 3.0 * xvx * r[2] / dL5) + coef2dot * (r[2] / dL3);
    }
}

// Compute magnetic flux density from a spherical-harmonic field model.
static void magneticfield(const double* posInertial, mwSize maxDegree, mwSize maxOrder,
                          const double* tm, double req, const double* G, 
                          const double* H, mwSize coefRows, mxLogical doJac,
                          double* fluxOut, double* jacOut) {
    // Transform to the body-fixed frame.
    double posBody[3] = { 0, 0, 0 };
    for (mwSize i = 0; i < 3; ++i) {
        for (mwSize j = 0; j < 3; ++j) {
            posBody[i] += tm[IDX(i, j, 3)] * posInertial[j];
        }
    }

    double x = posBody[0], y = posBody[1], z = posBody[2];

    // Spherical-coordinate terms.
    double r = sqrt(x * x + y * y + z * z);
    double phi = asin(z / r);
    double lam = atan2(y, x);

    double rho = sqrt(x * x + y * y);

    // Compute the Legendre arrays.
    mwSize dimP = (maxDegree + 1) * (maxDegree + 1);
    double* P = (double*)mxCalloc(dimP, sizeof(double));
    double* dP = (double*)mxCalloc(dimP, sizeof(double));
    double* ddP = (double*)mxCalloc(dimP, sizeof(double));

    legendre(maxDegree, maxOrder, phi, P, dP, ddP, NULL);

    // Accumulate the magnetic potential derivatives.
    double dUdr = 0, dUdp = 0, dUdl = 0;
    double d2U[6] = { 0, 0, 0, 0, 0, 0 };

    for (mwSize n = 0; n <= maxDegree;  ++n) {
        double Fn = req * pow(req / r, (double)n + 1);        
        double b1 = -((double)n + 1) * Fn / r;
        double b2 = Fn;
        double b3 = b2;

        double bJ[6];
        if (doJac) {
            bJ[0] = (n + 1) * (n + 2) * Fn / (r * r);
            bJ[1] = b1;
            bJ[2] = b1;
            bJ[3] = b2;
            bJ[4] = b2;
            bJ[5] = -b2;
        }

        double q1 = 0, q2 = 0, q3 = 0, qJ[6] = { 0, 0, 0, 0, 0, 0 };

        for (mwSize m = 0; m <= maxOrder; ++m) {
            double cosML = cos(m * lam);
            double sinML = sin(m * lam);
            size_t idx = IDX(n, m, maxDegree + 1);
            size_t idxCS = IDX(n, m, coefRows);

            double Smn = G[idxCS] * cosML + H[idxCS] * sinML;
            double Tmn = H[idxCS] * cosML - G[idxCS] * sinML;

            q1 += P[idx] * Smn;
            q2 += dP[idx] * Smn;
            q3 += (double)m * P[idx] * Tmn;

            if (doJac) {
                qJ[0] += P[idx] * Smn;
                qJ[1] += dP[idx] * Smn;
                qJ[2] += (double)m * P[idx] * Tmn;
                qJ[3] += ddP[idx] * Smn;
                qJ[4] += (double)m * dP[idx] * Tmn;
                qJ[5] += (double)(m * m) * P[idx] * Smn;
            }
        }

        dUdr += q1 * b1;
        dUdp += q2 * b2;
        dUdl += q3 * b3;

        if (doJac) {
            for (int k = 0; k < 6; ++k) {
                d2U[k] += qJ[k] * bJ[k];
            }
        }
    }

    // Convert to the body-fixed Cartesian frame.
    double rGrad[3] = { x / r, y / r, z / r };
    double phiGrad[3] = { -x * z / (rho * r * r), -y * z / (rho * r * r), rho / (r * r) };
    double lamGrad[3] = { -y / (rho * rho), x / (rho * rho), 0.0 };

    double rHess[9], phiHess[9], lamHess[9];
    if (doJac) {
        double r3 = r * r * r, rho3 = rho * rho * rho;
        double r4 = r3 * r, rho4 = rho3 * rho;
        double x2 = x * x, y2 = y * y, z2 = z * z;
        double x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;

        rHess[0] = (y2 + z2) / r3;
        rHess[1] = -x * y / r3;
        rHess[2] = -x * z / r3;
        rHess[3] = rHess[1];
        rHess[4] = (x2 + z2) / r3;
        rHess[5] = -y * z / r3;
        rHess[6] = rHess[2];
        rHess[7] = rHess[5];
        rHess[8] = (x2 + y2) / r3;

        double rho3r4 = rho3 * r3 * r;
        double rho1r4 = rho * r3 * r;

        phiHess[0] = z * (2 * x4 + x2 * y2 - y4 - y2 * z2) / rho3r4;
        phiHess[1] = x * y * z * (3 * x2 + 3 * y2 + z2) / rho3r4;
        phiHess[2] = -x * (x2 + y2 - z2) / rho1r4;
        phiHess[3] = phiHess[1];
        phiHess[4] = z * (-x4 + x2 * y2 - x2 * z2 + 2 * y4) / rho3r4;
        phiHess[5] = -y * (x2 + y2 - z2) / rho1r4;
        phiHess[6] = phiHess[2];
        phiHess[7] = phiHess[5];
        phiHess[8] = -2 * z * rho / r4;

        lamHess[0] = 2 * x * y / rho4;
        lamHess[1] = (y2 - x2) / rho4;
        lamHess[2] = 0.0;
        lamHess[3] = lamHess[1];
        lamHess[4] = -lamHess[0];
        lamHess[5] = 0.0;
        lamHess[6] = 0.0;
        lamHess[7] = 0.0;
        lamHess[8] = 0.0;
    }

    // Compute the magnetic flux density.
    double fluxBody[3];
    for (size_t i = 0; i < 3; ++i) {
        fluxBody[i] = dUdr * rGrad[i] + dUdp * phiGrad[i] + dUdl * lamGrad[i];
    }
    for (size_t i = 0; i < 3; ++i) {
        fluxOut[i] = 0.0;
        for (size_t j = 0; j < 3; ++j) {
            fluxOut[i] -= tm[IDX(j, i, 3)] * 1e-9 * fluxBody[j];
        }
    }
    
    if (!doJac) {
        mxFree(P);
        mxFree(dP);
        mxFree(ddP);
        return;
    }

    // Compute the Jacobian.
    if (doJac) {
        double jacBody[9] = { 0 };
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                size_t idx = IDX(i, j, 3);
                jacBody[idx] =
                    d2U[0] * rGrad[i] * rGrad[j] +
                    d2U[3] * phiGrad[i] * phiGrad[j] +
                    d2U[5] * lamGrad[i] * lamGrad[j] +
                    d2U[1] * (rGrad[i] * phiGrad[j] + rGrad[j] * phiGrad[i]) +
                    d2U[2] * (rGrad[i] * lamGrad[j] + rGrad[j] * lamGrad[i]) +
                    d2U[4] * (phiGrad[i] * lamGrad[j] + phiGrad[j] * lamGrad[i]) +
                    dUdr * rHess[idx] +
                    dUdp * phiHess[idx] +
                    dUdl * lamHess[idx];
            }
        }

        double tmp[9] = { 0 };
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    tmp[IDX(i, j, 3)] -= tm[IDX(k, i, 3)] * 1e-9 * jacBody[IDX(k, j, 3)];
                }
            }
        }
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                double val = 0;
                for (size_t k = 0; k < 3; ++k) {
                    val += tmp[IDX(i, k, 3)] * tm[IDX(k, j, 3)];
                }
                jacOut[IDX(i, j, 3)] = val;
                jacOut[IDX(i, j + 3, 3)] = 0.0;
            }
        }
        mxFree(P);
        mxFree(dP);
        mxFree(ddP);
    }
}

// Compute gravity-gradient torque using a quaternion attitude state.
static void gravitytorque_q(const double* jacInertial, const double* q, const double* inertiaMatrix, const mxLogical doJac, 
                            const double* hessInertial, double* torqueOut, double* jacOut_r, double* jacOut_q) {
    double tm[9], tmp[9], jacBody[9];
    q2R(q, tm); // Inertial to body-fixed.

    // tmp = R * JI
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            tmp[IDX(i, j, 3)] = jacInertial[3 * j] * tm[i] + jacInertial[1 + 3 * j] * tm[i + 3] + jacInertial[2 + 3 * j] * tm[i + 6];
        }
    }
    // jacBody = tmp * R.'
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            jacBody[IDX(i, j, 3)] = tm[j] * tmp[i] + tm[j + 3] * tmp[i + 3] + tm[j + 6] * tmp[i + 6];
        }
    }
    torqueOut[0] = -jacBody[6] * inertiaMatrix[3] + jacBody[3] * inertiaMatrix[6] +
                    jacBody[4] * inertiaMatrix[7] - jacBody[7] * inertiaMatrix[4] +
                    jacBody[7] * inertiaMatrix[8] - jacBody[8] * inertiaMatrix[7];
    torqueOut[1] = -jacBody[0] * inertiaMatrix[6] + jacBody[6] * inertiaMatrix[0] -
                    jacBody[3] * inertiaMatrix[7] + jacBody[7] * inertiaMatrix[3] -
                    jacBody[6] * inertiaMatrix[8] + jacBody[8] * inertiaMatrix[6];
    torqueOut[2] = -jacBody[3] * inertiaMatrix[0] + jacBody[0] * inertiaMatrix[3] +
                    jacBody[3] * inertiaMatrix[4] - jacBody[4] * inertiaMatrix[3] +
                    jacBody[6] * inertiaMatrix[7] - jacBody[7] * inertiaMatrix[6];
    
    if (!doJac) return;

    // Jacobian with respect to position.
    if (jacOut_r && hessInertial) {
        for (int k = 0; k < 3; ++k) {
            const double* Hk = hessInertial + 9 * k;
            double tmpH[9], dJB[9];

            // tmpH = R * Hk
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    tmpH[IDX(i, j, 3)] = Hk[3 * j] * tm[i] + Hk[1 + 3 * j] * tm[i + 3] + Hk[2 + 3 * j] * tm[i + 6];
                }
            }

            // dJB = tmpH * R.'
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    dJB[IDX(i, j, 3)] = tm[j] * tmpH[i] + tm[j + 3] * tmpH[i + 3] + tm[j + 6] * tmpH[i + 6];
                }
            }

            jacOut_r[IDX(0, k, 3)] = -dJB[6] * inertiaMatrix[3] + dJB[3] * inertiaMatrix[6] +
                                      dJB[4] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[4] + 
                                      dJB[7] * inertiaMatrix[8] - dJB[8] * inertiaMatrix[7];
            jacOut_r[IDX(1, k, 3)] = -dJB[0] * inertiaMatrix[6] + dJB[6] * inertiaMatrix[0] -
                                      dJB[3] * inertiaMatrix[7] + dJB[7] * inertiaMatrix[3] -
                                      dJB[6] * inertiaMatrix[8] + dJB[8] * inertiaMatrix[6];
            jacOut_r[IDX(2, k, 3)] = -dJB[3] * inertiaMatrix[0] + dJB[0] * inertiaMatrix[3] +
                                      dJB[3] * inertiaMatrix[4] - dJB[4] * inertiaMatrix[3] +
                                      dJB[6] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[6];
        }
    }

    // Jacobian with respect to the quaternion.
    if (jacOut_q && q) {
        double TMjac[36];
        q2Rjac(q, TMjac);

        for (int k = 0; k < 4; ++k) {
            double* U = TMjac + 9 * k;

            // A = U * R.'
            double A[9];
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    A[IDX(i, j, 3)] = tm[j] * U[i] + tm[j + 3] * U[i + 3] + tm[j + 6] * U[i + 6];
                }
            }

            // dJB = A * Jb + Jb * A.'
            double AtJ[9], JbA[9], dJB[9];
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    int idcur = IDX(i, j, 3);
                    // AtJ = A * Jb
                    AtJ[idcur] = A[i] * jacBody[3 * j] + A[i + 3] * jacBody[1 + 3 * j] + A[i + 6] * jacBody[2 + 3 * j];
                    // JbA = Jb * A.'
                    JbA[idcur] = jacBody[i] * A[j] + jacBody[i + 3] * A[j + 3] + jacBody[i + 6] * A[j + 6];
                    dJB[idcur] = AtJ[idcur] + JbA[idcur];
                }
            }

            jacOut_q[IDX(0, k, 3)] = -dJB[6] * inertiaMatrix[3] + dJB[3] * inertiaMatrix[6] +
                                      dJB[4] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[4] + 
                                      dJB[7] * inertiaMatrix[8] - dJB[8] * inertiaMatrix[7];
            jacOut_q[IDX(1, k, 3)] = -dJB[0] * inertiaMatrix[6] + dJB[6] * inertiaMatrix[0] -
                                      dJB[3] * inertiaMatrix[7] + dJB[7] * inertiaMatrix[3] -
                                      dJB[6] * inertiaMatrix[8] + dJB[8] * inertiaMatrix[6];
            jacOut_q[IDX(2, k, 3)] = -dJB[3] * inertiaMatrix[0] + dJB[0] * inertiaMatrix[3] +
                                      dJB[3] * inertiaMatrix[4] - dJB[4] * inertiaMatrix[3] +
                                      dJB[6] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[6];
        }
    }
}

// Compute gravity-gradient torque using a rotation-matrix attitude state.
static void gravitytorque_R(const double* jacInertial, const double* tm, const double* inertiaMatrix, const mxLogical doJac, 
                            const double* hessInertial, double* torqueOut, double* jacOut_r, double* jacOut_R) {
    double tmp[9], jacBody[9];

    // tmp = R * JI
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            tmp[IDX(i, j, 3)] = jacInertial[3 * j] * tm[i] + jacInertial[1 + 3 * j] * tm[i + 3] + jacInertial[2 + 3 * j] * tm[i + 6];
        }
    }
    // jacBody = tmp * R.'
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            jacBody[IDX(i, j, 3)] = tm[j] * tmp[i] + tm[j + 3] * tmp[i + 3] + tm[j + 6] * tmp[i + 6];
        }
    }
    torqueOut[0] = -jacBody[6] * inertiaMatrix[3] + jacBody[3] * inertiaMatrix[6] +
                    jacBody[4] * inertiaMatrix[7] - jacBody[7] * inertiaMatrix[4] +
                    jacBody[7] * inertiaMatrix[8] - jacBody[8] * inertiaMatrix[7];
    torqueOut[1] = -jacBody[0] * inertiaMatrix[6] + jacBody[6] * inertiaMatrix[0] -
                    jacBody[3] * inertiaMatrix[7] + jacBody[7] * inertiaMatrix[3] -
                    jacBody[6] * inertiaMatrix[8] + jacBody[8] * inertiaMatrix[6];
    torqueOut[2] = -jacBody[3] * inertiaMatrix[0] + jacBody[0] * inertiaMatrix[3] +
                    jacBody[3] * inertiaMatrix[4] - jacBody[4] * inertiaMatrix[3] +
                    jacBody[6] * inertiaMatrix[7] - jacBody[7] * inertiaMatrix[6];
    
    if (!doJac) return;

    // Jacobian with respect to position.
    if (jacOut_r && hessInertial) {
        for (int k = 0; k < 3; ++k) {
            const double* Hk = hessInertial + 9 * k;
            double tmpH[9], dJB[9];

            // tmpH = R * Hk
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    tmpH[IDX(i, j, 3)] = Hk[3 * j] * tm[i] + Hk[1 + 3 * j] * tm[i + 3] + Hk[2 + 3 * j] * tm[i + 6];
                }
            }

            // dJB = tmpH * R.'
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    dJB[IDX(i, j, 3)] = tm[j] * tmpH[i] + tm[j + 3] * tmpH[i + 3] + tm[j + 6] * tmpH[i + 6];
                }
            }

            jacOut_r[IDX(0, k, 3)] = -dJB[6] * inertiaMatrix[3] + dJB[3] * inertiaMatrix[6] +
                                      dJB[4] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[4] + 
                                      dJB[7] * inertiaMatrix[8] - dJB[8] * inertiaMatrix[7];
            jacOut_r[IDX(1, k, 3)] = -dJB[0] * inertiaMatrix[6] + dJB[6] * inertiaMatrix[0] -
                                      dJB[3] * inertiaMatrix[7] + dJB[7] * inertiaMatrix[3] -
                                      dJB[6] * inertiaMatrix[8] + dJB[8] * inertiaMatrix[6];
            jacOut_r[IDX(2, k, 3)] = -dJB[3] * inertiaMatrix[0] + dJB[0] * inertiaMatrix[3] +
                                      dJB[3] * inertiaMatrix[4] - dJB[4] * inertiaMatrix[3] +
                                      dJB[6] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[6];
        }
    }

    // Jacobian with respect to the rotation matrix.
    if (jacOut_R && tm) {
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                int idxR = 3 * b + a;

                // Construct E_ab.
                double E[9] = {0.0};
                E[IDX(a, b, 3)] = 1.0;

                // term1 = E * jacInertial * R.'
                double EJ[9], term1[9];
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        EJ[IDX(i, j, 3)] =
                            E[IDX(i, 0, 3)] * jacInertial[IDX(0, j, 3)] +
                            E[IDX(i, 1, 3)] * jacInertial[IDX(1, j, 3)] +
                            E[IDX(i, 2, 3)] * jacInertial[IDX(2, j, 3)];
                    }
                }
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        term1[IDX(i, j, 3)] =
                            EJ[IDX(i, 0, 3)] * tm[IDX(j, 0, 3)] +
                            EJ[IDX(i, 1, 3)] * tm[IDX(j, 1, 3)] +
                            EJ[IDX(i, 2, 3)] * tm[IDX(j, 2, 3)];
                    }
                }

                // term2 = R * jacInertial * E.'
                double RE[9], term2[9];
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        RE[IDX(i, j, 3)] =
                            tm[IDX(i, 0, 3)] * jacInertial[IDX(0, j, 3)] +
                            tm[IDX(i, 1, 3)] * jacInertial[IDX(1, j, 3)] +
                            tm[IDX(i, 2, 3)] * jacInertial[IDX(2, j, 3)];
                    }
                }
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        term2[IDX(i, j, 3)] =
                            RE[IDX(i, 0, 3)] * E[IDX(j, 0, 3)] +
                            RE[IDX(i, 1, 3)] * E[IDX(j, 1, 3)] +
                            RE[IDX(i, 2, 3)] * E[IDX(j, 2, 3)];
                    }
                }

                // dJB = term1 + term2
                double dJB[9];
                for (int i = 0; i < 9; ++i) {
                    dJB[i] = term1[i] + term2[i];
                }

                jacOut_R[IDX(0, idxR, 3)] = -dJB[6] * inertiaMatrix[3] + dJB[3] * inertiaMatrix[6] +
                                             dJB[4] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[4] +
                                             dJB[7] * inertiaMatrix[8] - dJB[8] * inertiaMatrix[7];
                jacOut_R[IDX(1, idxR, 3)] = -dJB[0] * inertiaMatrix[6] + dJB[6] * inertiaMatrix[0] -
                                             dJB[3] * inertiaMatrix[7] + dJB[7] * inertiaMatrix[3] -
                                             dJB[6] * inertiaMatrix[8] + dJB[8] * inertiaMatrix[6];
                jacOut_R[IDX(2, idxR, 3)] = -dJB[3] * inertiaMatrix[0] + dJB[0] * inertiaMatrix[3] +
                                             dJB[3] * inertiaMatrix[4] - dJB[4] * inertiaMatrix[3] +
                                             dJB[6] * inertiaMatrix[7] - dJB[7] * inertiaMatrix[6];
            }
        }
    }
}

// Compute N-plate solar-radiation-pressure torque using a quaternion attitude state.
static void srpnplatetorque_q(const double* forceInertial, const double* q, const double* displacements, const mxLogical doJac,
                              const double* jacInertial, int numPlates, double* torqueOut, double* jacOut_r, double* jacOut_q) {
    double tm[9]; // Inertial to body-fixed.
    q2R(q, tm);

    for (int i = 0; i < numPlates; ++i) {
        double forceBody[3];
        // Transform to the body-fixed frame.
        forceBody[0] = tm[0] * forceInertial[IDX(i, 0, numPlates)] + tm[3] * forceInertial[IDX(i, 1, numPlates)] + tm[6] * forceInertial[IDX(i, 2, numPlates)];
        forceBody[1] = tm[1] * forceInertial[IDX(i, 0, numPlates)] + tm[4] * forceInertial[IDX(i, 1, numPlates)] + tm[7] * forceInertial[IDX(i, 2, numPlates)];
        forceBody[2] = tm[2] * forceInertial[IDX(i, 0, numPlates)] + tm[5] * forceInertial[IDX(i, 1, numPlates)] + tm[8] * forceInertial[IDX(i, 2, numPlates)];

        // Accumulate torque.
        torqueOut[0] += displacements[IDX(i, 1, numPlates)] * forceBody[2] - displacements[IDX(i, 2, numPlates)] * forceBody[1];
        torqueOut[1] += displacements[IDX(i, 2, numPlates)] * forceBody[0] - displacements[IDX(i, 0, numPlates)] * forceBody[2];
        torqueOut[2] += displacements[IDX(i, 0, numPlates)] * forceBody[1] - displacements[IDX(i, 1, numPlates)] * forceBody[0];
    }

    if (!doJac) return;

    double q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];

    for (int i = 0; i < numPlates; ++i) {
        // Moment arm.
        double dx = displacements[IDX(i, 0, numPlates)];
        double dy = displacements[IDX(i, 1, numPlates)];
        double dz = displacements[IDX(i, 2, numPlates)];

        // Jacobian with respect to position.
        if (jacOut_r && jacInertial) {
            // Accumulate the three position components.
            for (int k = 0; k < 3; ++k) {
                double gI0 = jacInertial[IDX3(0, k, i, 3, 13)];
                double gI1 = jacInertial[IDX3(1, k, i, 3, 13)];
                double gI2 = jacInertial[IDX3(2, k, i, 3, 13)];

                double gB0 = tm[0] * gI0 + tm[3] * gI1 + tm[6] * gI2;
                double gB1 = tm[1] * gI0 + tm[4] * gI1 + tm[7] * gI2;
                double gB2 = tm[2] * gI0 + tm[5] * gI1 + tm[8] * gI2;

                jacOut_r[IDX(0, k, 3)] += dy * gB2 - dz * gB1;
                jacOut_r[IDX(1, k, 3)] += dz * gB0 - dx * gB2;
                jacOut_r[IDX(2, k, 3)] += dx * gB1 - dy * gB0;
            }
        }

        // Jacobian with respect to attitude.
        if (jacOut_q) {
            double TMjac[36];
            q2Rjac(q, TMjac);

            for (int k = 0; k < 4; ++k) {
                double* U = TMjac + 9 * k;

                double fI0 = forceInertial[IDX(i, 0, numPlates)];
                double fI1 = forceInertial[IDX(i, 1, numPlates)];
                double fI2 = forceInertial[IDX(i, 2, numPlates)];

                // Contribution from the rotation-matrix variation.
                double dFBT0 = U[0] * fI0 + U[3] * fI1 + U[6] * fI2;
                double dFBT1 = U[1] * fI0 + U[4] * fI1 + U[7] * fI2;
                double dFBT2 = U[2] * fI0 + U[5] * fI1 + U[8] * fI2;

                // Contribution from the inertial-force variation.
                double gIq0 = 0.0, gIq1 = 0.0, gIq2 = 0.0;
                if (jacInertial) {
                    int col = 9 + k;
                    gIq0 = jacInertial[IDX3(0, col, i, 3, 13)];
                    gIq1 = jacInertial[IDX3(1, col, i, 3, 13)];
                    gIq2 = jacInertial[IDX3(2, col, i, 3, 13)];
                }
                double dFBF0 = tm[0] * gIq0 + tm[3] * gIq1 + tm[6] * gIq2;
                double dFBF1 = tm[1] * gIq0 + tm[4] * gIq1 + tm[7] * gIq2;
                double dFBF2 = tm[2] * gIq0 + tm[5] * gIq1 + tm[8] * gIq2;

                double dFB0 = dFBT0 + dFBF0;
                double dFB1 = dFBT1 + dFBF1;
                double dFB2 = dFBT2 + dFBF2;

                jacOut_q[IDX(0, k, 3)] += dy * dFB2 - dz * dFB1;
                jacOut_q[IDX(1, k, 3)] += dz * dFB0 - dx * dFB2;
                jacOut_q[IDX(2, k, 3)] += dx * dFB1 - dy * dFB0;
            }
        }
    }
}

// Compute N-plate solar-radiation-pressure torque using a rotation-matrix attitude state.
static void srpnplatetorque_R(const double* forceInertial, const double* tm, const double* displacements, const mxLogical doJac,
                              const double* jacInertial, int numPlates, double* torqueOut, double* jacOut_r, double* jacOut_R) {
    for (int i = 0; i < numPlates; ++i) {
        double forceBody[3];
        // Transform to the body-fixed frame.
        forceBody[0] = tm[0] * forceInertial[IDX(i, 0, numPlates)] + tm[3] * forceInertial[IDX(i, 1, numPlates)] + tm[6] * forceInertial[IDX(i, 2, numPlates)];
        forceBody[1] = tm[1] * forceInertial[IDX(i, 0, numPlates)] + tm[4] * forceInertial[IDX(i, 1, numPlates)] + tm[7] * forceInertial[IDX(i, 2, numPlates)];
        forceBody[2] = tm[2] * forceInertial[IDX(i, 0, numPlates)] + tm[5] * forceInertial[IDX(i, 1, numPlates)] + tm[8] * forceInertial[IDX(i, 2, numPlates)];

        // Accumulate torque.
        torqueOut[0] += displacements[IDX(i, 1, numPlates)] * forceBody[2] - displacements[IDX(i, 2, numPlates)] * forceBody[1];
        torqueOut[1] += displacements[IDX(i, 2, numPlates)] * forceBody[0] - displacements[IDX(i, 0, numPlates)] * forceBody[2];
        torqueOut[2] += displacements[IDX(i, 0, numPlates)] * forceBody[1] - displacements[IDX(i, 1, numPlates)] * forceBody[0];
    }

    if (!doJac) return;

    for (int i = 0; i < numPlates; ++i) {
        // Moment arm.
        double dx = displacements[IDX(i, 0, numPlates)];
        double dy = displacements[IDX(i, 1, numPlates)];
        double dz = displacements[IDX(i, 2, numPlates)];

        // Jacobian with respect to position.
        if (jacOut_r && jacInertial) {
            // Accumulate the three position components.
            for (int k = 0; k < 3; ++k) {
                double gI0 = jacInertial[IDX3(0, k, i, 3, 18)];
                double gI1 = jacInertial[IDX3(1, k, i, 3, 18)];
                double gI2 = jacInertial[IDX3(2, k, i, 3, 18)];

                double gB0 = tm[0] * gI0 + tm[3] * gI1 + tm[6] * gI2;
                double gB1 = tm[1] * gI0 + tm[4] * gI1 + tm[7] * gI2;
                double gB2 = tm[2] * gI0 + tm[5] * gI1 + tm[8] * gI2;

                jacOut_r[IDX(0, k, 3)] += dy * gB2 - dz * gB1;
                jacOut_r[IDX(1, k, 3)] += dz * gB0 - dx * gB2;
                jacOut_r[IDX(2, k, 3)] += dx * gB1 - dy * gB0;
            }
        }

        // Jacobian with respect to rotation matrix R.
        if (jacOut_R) {
            double fI0 = forceInertial[IDX(i, 0, numPlates)];
            double fI1 = forceInertial[IDX(i, 1, numPlates)];
            double fI2 = forceInertial[IDX(i, 2, numPlates)];

            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    int idxR = 3 * b + a;

                    // d(R*fI)/dR_ab = E_ab * fI
                    double dFBT0 = 0.0, dFBT1 = 0.0, dFBT2 = 0.0;
                    if (a == 0) dFBT0 = (b == 0 ? fI0 : (b == 1 ? fI1 : fI2));
                    if (a == 1) dFBT1 = (b == 0 ? fI0 : (b == 1 ? fI1 : fI2));
                    if (a == 2) dFBT2 = (b == 0 ? fI0 : (b == 1 ? fI1 : fI2));

                    // R * d fI / dR_ab
                    double gIR0 = 0.0, gIR1 = 0.0, gIR2 = 0.0;
                    if (jacInertial) {
                        int col = 9 + idxR;
                        gIR0 = jacInertial[IDX3(0, col, i, 3, 18)];
                        gIR1 = jacInertial[IDX3(1, col, i, 3, 18)];
                        gIR2 = jacInertial[IDX3(2, col, i, 3, 18)];
                    }

                    double dFBF0 = tm[0] * gIR0 + tm[3] * gIR1 + tm[6] * gIR2;
                    double dFBF1 = tm[1] * gIR0 + tm[4] * gIR1 + tm[7] * gIR2;
                    double dFBF2 = tm[2] * gIR0 + tm[5] * gIR1 + tm[8] * gIR2;

                    // Total derivative.
                    double dFB0 = dFBT0 + dFBF0;
                    double dFB1 = dFBT1 + dFBF1;
                    double dFB2 = dFBT2 + dFBF2;

                    // Torque derivative: d tau = r_arm x d f_B.
                    jacOut_R[IDX(0, idxR, 3)] += dy * dFB2 - dz * dFB1;
                    jacOut_R[IDX(1, idxR, 3)] += dz * dFB0 - dx * dFB2;
                    jacOut_R[IDX(2, idxR, 3)] += dx * dFB1 - dy * dFB0;
                }
            }
        }
    }
}

// Compute the overlap area of two circles.
static double area2circles(double r1, double r2, double d, mxLogical doHess, double* gradOut, double* hessOut) {
    double eps = 1e-15;
    r1 = fmax(r1, 0.0);
    r2 = fmax(r2, 0.0);
    d = fmax(d, 0.0);

    double r1s = r1 * r1, r2s = r2 * r2, ds = d * d;

    // Disjoint or externally tangent circles.
    if (d >= r1 + r2 - eps) {
        return 0.0;
    }

    // Contained or internally tangent circles.
    if (d <= fabs(r1 - r2) + eps) {
        if (r1 <= r2) {
            gradOut[0] = 2 * PI * r1;
            if (hessOut) hessOut[0] = 2 * PI;
            return PI * r1s;
        }
        else {
            gradOut[1] = 2 * PI * r2;
            if (hessOut) hessOut[4] = 2 * PI;
            return PI * r2s;
        }
    }

    // Coincident circles.
    if (d <= eps && fabs(r1 - r2) < eps) {
        gradOut[0] = 2 * PI * r1;
        gradOut[1] = 2 * PI * r2;
        if (hessOut) {
            hessOut[0] = 2 * PI;
            hessOut[4] = 2 * PI;
        }
        return PI * r1s;
    }

    double x1 = (ds + r1s - r2s) / (2 * d * r1), x2 = (ds + r2s - r1s) / (2 * d * r2);
    x1 = fmax(-1.0, fmin(1.0, x1));
    x2 = fmax(-1.0, fmin(1.0, x2));
    double alpha = acos(x1), beta = acos(x2);

    double f1, f2, f3, f4, delta, invf1, invf2, invf3, invf4;
    f1 = -d + r1 + r2;
    f2 = d + r1 - r2;
    f3 = d - r1 + r2;
    f4 = d + r1 + r2;
    delta = sqrt(fmax(0.0, f1 * f2 * f3 * f4));

    invf1 = 1.0 / f1;
    invf2 = 1.0 / f2;
    invf3 = 1.0 / f3;
    invf4 = 1.0 / f4;

    double gradX1[3], gradX2[3], gradDelta[3];
    gradX1[0] = (r1s + r2s - ds) / (2 * d * r1s);
    gradX1[1] = -r2 / (d * r1);
    gradX1[2] = (ds - r1s + r2s) / (2 * ds * r1);
    gradX2[0] = -r1 / (d * r2);
    gradX2[1] = (r2s + r1s - ds) / (2 * d * r2s);
    gradX2[2] = (ds - r2s + r1s) / (2 * ds * r2);
    gradDelta[0] = 0.5 * delta * (invf1 + invf2 - invf3 + invf4);
    gradDelta[1] = 0.5 * delta * (invf1 - invf2 + invf3 + invf4);
    gradDelta[2] = 0.5 * delta * (-invf1 + invf2 + invf3 + invf4);

    double s1inv = 1.0 / sqrt(1 - x1 * x1), s2inv = 1.0 / sqrt(1 - x2 * x2);
    double coef1 = r1s * s1inv, coef2 = r2s * s2inv;
    gradOut[0] = 2 * r1 * alpha - coef1 * gradX1[0] - coef2 * gradX2[0] - 0.5 * gradDelta[0];
    gradOut[1] = 2 * r2 * beta - coef1 * gradX1[1] - coef2 * gradX2[1] - 0.5 * gradDelta[1];
    gradOut[2] = -coef1 * gradX1[2] - coef2 * gradX2[2] - 0.5 * gradDelta[2];

    if (doHess) {
        double dinv = 1.0 / d, dinv2 = dinv * dinv, dinv3 = dinv2 * dinv;
        double r1inv = 1.0 / r1, r1inv2 = r1inv * r1inv, r1inv3 = r1inv2 * r1inv;
        double r2inv = 1.0 / r2, r2inv2 = r2inv * r2inv, r2inv3 = r2inv2 * r2inv;
        double Hx1[9], Hx2[9];

        // Hessian with respect to x1.
        Hx1[0] = (ds - r2s) * dinv * r1inv3;
        Hx1[1] = r2 * dinv * r1inv2;
        Hx1[2] = -0.5 * (ds + r1s + r2s) * dinv2 * r1inv2;
        Hx1[3] = Hx1[1];
        Hx1[4] = -dinv * r1inv;
        Hx1[5] = r2 * dinv2 * r1inv;
        Hx1[6] = Hx1[2];
        Hx1[7] = Hx1[5];
        Hx1[8] = (r1s - r2s) * dinv3 * r1inv;

        // Hessian with respect to x2.
        Hx2[0] = -dinv * r2inv;
        Hx2[1] = r1 * dinv * r2inv2;
        Hx2[2] = r1 * dinv2 * r2inv;
        Hx2[3] = Hx2[1];
        Hx2[4] = (ds - r1s) * dinv * r2inv3;
        Hx2[5] = -0.5 * (ds + r1s + r2s) * dinv2 * r2inv2;
        Hx2[6] = Hx2[2];
        Hx2[7] = Hx2[5];
        Hx2[8] = (-r1s + r2s) * dinv3 * r2inv;

        // Gradients of alpha and beta.
        double gradAlpha[3] = {-gradX1[0] * s1inv, -gradX1[1] * s1inv, -gradX1[2] * s1inv};
        double gradBeta[3] = {-gradX2[0] * s2inv, -gradX2[1] * s2inv, -gradX2[2] * s2inv};

        // Hessians of alpha and beta.
        double s1inv3 = s1inv * s1inv * s1inv, s2inv3 = s2inv * s2inv * s2inv;
        double k1 = x1 * s1inv3, k2 = x2 * s2inv3;
        double Halpha[9], Hbeta[9];
        for (int i = 0; i < 9; ++i) {
            Halpha[i] = -Hx1[i] * s1inv;
            Hbeta[i] = -Hx2[i] * s2inv;
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Halpha[IDX(i, j, 3)] += -k1 * gradX1[i] * gradX1[j];
                Hbeta[IDX(i, j, 3)] += -k2 * gradX2[i] * gradX2[j];
            }
        }

        // Hessian of delta.
        double Svec[3] = {
            invf1 + invf2 - invf3 + invf4,
            invf1 - invf2 + invf3 + invf4,
            -invf1 + invf2 + invf3 + invf4
        };
        double Hdelta[9] = { 0 };
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Hdelta[IDX(i, j, 3)] += 0.25 * delta * Svec[i] * Svec[j];
            }
        }
        double GG[4][3] = {
            {1.0, 1.0, -1.0},
            {1.0, -1.0, 1.0},
            {-1.0, 1.0, 1.0},
            {1.0, 1.0, 1.0}
        };
        double WI[4] = {invf1 * invf1, invf2 * invf2, invf3 * invf3, invf4 * invf4};
        for (int t = 0; t < 4; ++t) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    Hdelta[IDX(i, j, 3)] += -0.5 * delta * WI[t] * GG[t][i] * GG[t][j];
                }
            }
        }

        // Final Hessian.
        for (int i = 0; i < 9; ++i) {
            hessOut[i] += r1s * Halpha[i] + r2s * Hbeta[i];
        }
        for (int j = 0; j < 3; ++j) {
            hessOut[IDX(0, j, 3)] += 2.0 * r1 * gradAlpha[j];
            hessOut[IDX(j, 0, 3)] += 2.0 * r1 * gradAlpha[j];
        }
        hessOut[IDX(0, 0, 3)] += 2.0 * alpha;
        for (int j = 0; j < 3; ++j) {
            hessOut[IDX(1, j, 3)] += 2.0 * r2 * gradBeta[j];
            hessOut[IDX(j, 1, 3)] += 2.0 * r2 * gradBeta[j];
        }
        hessOut[IDX(1, 1, 3)] += 2.0 * beta;
        for (int i = 0; i < 9; ++i) {
            hessOut[i] += -0.5 * Hdelta[i];
        }
    }

    return r1s * alpha + r2s * beta - 0.5 * delta;
}

// Compute the fraction of sunlight blocked by Earth and the Moon using a dual-cone model.
static double dualcone(const double* posInertial, const double* posSun, const double* posEarth,
                       const double* velInertial, const double* velSun, const double* velEarth,
                       const mxLogical doJac, double* gradOut, double* diffOut, double* diffGradOut) {
    const double Rs = 696000, Re = 6371, Rm = 1738;

    double s[3], e[3], m[3], snorm = 0.0, enorm = 0.0, mnorm = 0.0;
    for (int i = 0; i < 3; ++i) {
        s[i] = posSun[i] - posInertial[i];
        e[i] = posEarth[i] - posInertial[i];
        m[i] = -posInertial[i];
        snorm += s[i] * s[i];
        enorm += e[i] * e[i];
        mnorm += m[i] * m[i];
    }
    snorm = sqrt(snorm);
    enorm = sqrt(enorm);
    mnorm = sqrt(mnorm);
    for (int i = 0; i < 3; ++i) {
        s[i] /= snorm;
        e[i] /= enorm;
        m[i] /= mnorm;
    }

    double a = asin(Rs / snorm);  // Apparent solar radius.
    double bE = asin(Re / enorm); // Apparent Earth radius.
    double bM = asin(Rm / mnorm); // Apparent lunar radius.

    double uSE = 0.0, uSM = 0.0, uEM = 0.0;
    for (int i = 0; i < 3; ++i) {
        uSE += s[i] * e[i];
        uSM += s[i] * m[i];
        uEM += e[i] * m[i];
    }
    double cE = acos(uSE), cM = acos(uSM), dEM = acos(uEM);
    double cosphi = (cE * cE + cM * cM - dEM * dEM) / (2 * cE * cM), phi = acos(cosphi);

    double gradA_SE[3] = {0, 0, 0}, gradA_SM[3] = {0, 0, 0};
    double hessA_SE[9] = { 0 }, hessA_SM[9] = { 0 };

    double A_SE = area2circles(a, bE, cE, doJac, gradA_SE, hessA_SE);
    double A_SM = area2circles(a, bM, cM, doJac, gradA_SM, hessA_SM);

    double A = A_SE + A_SM;
    double occ = A / (PI * a * a);

    // Compute the gradient.
    double denSE = sqrt(1 - uSE * uSE), denSM = sqrt(1 - uSM * uSM);

    // Apparent-radius gradients.
    double snorm2 = snorm * snorm, enorm2 = enorm * enorm, mnorm2 = mnorm * mnorm;
    double Rs2 = Rs * Rs, Re2 = Re * Re, Rm2 = Rm * Rm;
    double gradR0[3], gradR1[3], gradR2[3];
    double rootS = sqrt(snorm2 - Rs2), rootE = sqrt(enorm2 - Re2), rootM = sqrt(mnorm2 - Rm2);
    double ts = Rs / (snorm * rootS);
    double te = Re / (enorm * rootE);
    double tm = Rm / (mnorm * rootM);
    for (int i = 0; i < 3; ++i) {
        gradR0[i] = ts * s[i];
        gradR1[i] = te * e[i];
        gradR2[i] = tm * m[i];
    }

    // Angular-separation gradients.
    double gradD01[3], gradD02[3];
    for (int i = 0; i < 3; ++i) {
        gradD01[i] = ((e[i] - uSE * s[i]) / snorm + (s[i] - uSE * e[i]) / enorm) / denSE;
        gradD02[i] = ((m[i] - uSM * s[i]) / snorm + (s[i] - uSM * m[i]) / mnorm) / denSM;
    }

    // Assemble the final gradient.
    double cc1 = -1.0 / (PI * a * a), cc2 = 2 * A / (PI * a * a * a);
    for (int i = 0; i < 3; ++i) {
        double gradA = (gradA_SE[0] + gradA_SM[0]) * gradR0[i] + gradA_SE[1] * gradR1[i] + gradA_SM[1] * gradR2[i] +
                        gradA_SE[2] * gradD01[i] + gradA_SM[2] * gradD02[i];
        gradOut[i] = cc1 * gradA + cc2 * gradR0[i];
    }

    // Compute the time derivative.
    if (velInertial && velSun && velEarth) {
        double vS[3], vE[3], vM[3];
        for (int i = 0; i < 3; ++i) {
            vS[i] = velSun[i] - velInertial[i];
            vE[i] = velEarth[i] - velInertial[i];
            vM[i] = -velInertial[i];
        }

        // Time derivatives of the apparent radii.
        double svS = 0.0, evE = 0.0, mvM = 0.0;
        for (int i = 0; i < 3; ++i) {
            svS += s[i] * vS[i];
            evE += e[i] * vE[i];
            mvM += m[i] * vM[i];
        }
        double adot = -ts * svS;
        double bEdot = -te * evE;
        double bMdot = -tm * mvM;

        // Time derivatives of the angular separations.
        double dotSE1 = 0.0, dotSE2 = 0.0, dotSM1 = 0.0, dotSM2 = 0.0;
        for (int i = 0; i < 3; ++i) {
            dotSE1 += (e[i] - uSE * s[i]) * vS[i];
            dotSE2 += (s[i] - uSE * e[i]) * vE[i];
            dotSM1 += (m[i] - uSM * s[i]) * vS[i];
            dotSM2 += (s[i] - uSM * m[i]) * vM[i];
        }
        double cEdot = -(dotSE1 / snorm + dotSE2 / enorm) / denSE;
        double cMdot = -(dotSM1 / snorm + dotSM2 / mnorm) / denSM;

        // Assemble the time derivative.
        double Adot = (gradA_SE[0] + gradA_SM[0]) * adot + gradA_SE[1] * bEdot + gradA_SM[1] * bMdot + 
                       gradA_SE[2] * cEdot + gradA_SM[2] * cMdot;
        *diffOut = cc1 * Adot + cc2 * adot;

        // Compute the Jacobian of the time derivative.
        if (doJac) {
            // Gradient with respect to velocity.
            diffGradOut[3] = gradOut[0];
            diffGradOut[4] = gradOut[1];
            diffGradOut[5] = gradOut[2];

            // Prepare geometric chain-rule derivatives.
            double Ds[9], De[9], Dm[9];
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    double Iij = (i == j) ? 1.0 : 0.0;
                    Ds[IDX(i, j, 3)] = -(Iij - s[i] * s[j]) / snorm;
                    De[IDX(i, j, 3)] = -(Iij - e[i] * e[j]) / enorm;
                    Dm[IDX(i, j, 3)] = -(Iij - m[i] * m[j]) / mnorm;
                }
            }

            // Position gradients of adot, bEdot, and bMdot.
            double dadotdp[3] = { 0 }, dbEdotdp[3] = { 0 }, dbMdotdp[3] = { 0 };
            double tsdp = Rs * (2.0 * snorm2 - Rs2) / (snorm2 * rootS * rootS * rootS);
            double tedp = Re * (2.0 * enorm2 - Re2) / (enorm2 * rootE * rootE * rootE);
            double tmdp = Rm * (2.0 * mnorm2 - Rm2) / (mnorm2 * rootM * rootM * rootM);

            for (int j = 0; j < 3; ++j) {
                double dsvSdp = 0.0, devEdp = 0.0, dmvMdp = 0.0;
                for (int i = 0; i < 3; ++i) {
                    dsvSdp += vS[i] * Ds[IDX(i, j, 3)];
                    devEdp += vE[i] * De[IDX(i, j, 3)];
                    dmvMdp += vM[i] * Dm[IDX(i, j, 3)];
                }

                dadotdp[j] = -(tsdp * s[j] * svS + ts * dsvSdp);
                dbEdotdp[j] = -(tedp * e[j] * evE + te * devEdp);
                dbMdotdp[j] = -(tmdp * m[j] * mvM + tm * dmvMdp);
            }

            // Position gradients of cEdot and cMdot.
            double duSEdp[3] = { 0 }, duSMdp[3] = { 0 };
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    duSEdp[j] += Ds[IDX(i, j, 3)] * e[i] + De[IDX(i, j, 3)] * s[i];
                    duSMdp[j] += Ds[IDX(i, j, 3)] * m[i] + Dm[IDX(i, j, 3)] * s[i];
                }
            }

            double ddenSEdp[3], ddenSMdp[3];
            for (int i = 0; i < 3; ++i) {
                ddenSEdp[i] = -(uSE / denSE) * duSEdp[i];
                ddenSMdp[i] = -(uSM / denSM) * duSMdp[i];
            }

            double dq1dp[3], dq1Mdp[3];
            for (int j = 0; j < 3; ++j) {
                double dDotSE1dp = 0.0, dDotSE2dp = 0.0, dDotSM1dp = 0.0, dDotSM2dp = 0.0;;

                for (int i = 0; i < 3; ++i) {
                    double dsidp = Ds[IDX(i, j, 3)];
                    double deidp = De[IDX(i, j, 3)];
                    double dmidp = Dm[IDX(i, j, 3)];

                    dDotSE1dp += (deidp - duSEdp[j] * s[i] - uSE * dsidp) * vS[i];
                    dDotSE2dp += (dsidp - duSEdp[j] * e[i] - uSE * deidp) * vE[i];
                    dDotSM1dp += (dmidp - duSMdp[j] * s[i] - uSM * dsidp) * vS[i];
                    dDotSM2dp += (dsidp - duSMdp[j] * m[i] - uSM * dmidp) * vM[i];
                }
                
                double dsnormdp = -s[j];
                double denormdp = -e[j];
                double dmnormdp = -m[j];

                double dA1dp = dDotSE1dp / snorm - dotSE1 / (snorm * snorm) * dsnormdp;
                double dA2dp = dDotSE2dp / enorm - dotSE2 / (enorm * enorm) * denormdp;
                double dA1Mdp = dDotSM1dp / snorm - dotSM1 / (snorm * snorm) * dsnormdp;
                double dA2Mdp = dDotSM2dp / mnorm - dotSM2 / (mnorm * mnorm) * dmnormdp;

                dq1dp[j] = dA1dp + dA2dp;
                dq1Mdp[j] = dA1Mdp + dA2Mdp;
            }

            double dcEdotdp[3], dcMdotdp[3];;
            for (int i = 0; i < 3; ++i) {
                double num = dq1dp[i] * denSE - (dotSE1 / snorm + dotSE2 / enorm) * ddenSEdp[i];
                double numM = dq1Mdp[i] * denSM - (dotSM1 / snorm + dotSM2 / mnorm) * ddenSMdp[i];
                dcEdotdp[i] = -num / (denSE * denSE);
                dcMdotdp[i] = -numM / (denSM * denSM);
            }

            // Position gradients of the area and the c1/c2 coefficients.
            double Ap[3];
            for (int i = 0; i < 3; ++i) {
                Ap[i] = (gradA_SE[0] + gradA_SM[0]) * gradR0[i] + gradA_SE[1] * gradR1[i] + gradA_SM[1] * gradR2[i] +
                    gradA_SE[2] * gradD01[i] + gradA_SM[2] * gradD02[i];
            }

            double invpi = 1.0 / PI, a3inv = 1.0 / (a * a * a), a4inv = a3inv / a;
            double dcc1dp[3], dcc2dp[3];
            for (int i = 0; i < 3; ++i) {
                dcc1dp[i] = (2.0 * invpi * a3inv) * gradR0[i];
                dcc2dp[i] = (2.0 * invpi * a3inv) * Ap[i] - (6.0 * A * invpi * a4inv) * gradR0[i];
            }

            // Assemble the final Jacobian.
            double* TSE[3] = {gradR0, gradR1, gradD01};
            double* TSM[3] = {gradR0, gradR2, gradD02};

            double USE[9] = { 0 }, USM[9] = { 0 };
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        USE[IDX(i, j, 3)] += hessA_SE[IDX(i, k, 3)] * TSE[k][j];
                        USM[IDX(i, j, 3)] += hessA_SM[IDX(i, k, 3)] * TSM[k][j];
                    }
                }
            }

            // dAdotdp
            double dAdotdp[3] = { 0 };
            double thetaDotSE[3] = {adot, bEdot, cEdot};
            double thetaDotSM[3] = {adot, bMdot, cMdot};

            for (int j = 0; j < 3; ++j) {
                double term1 = 0.0;
                for (int k = 0; k < 3; ++k) {
                    term1 += USE[IDX(k, j, 3)] * thetaDotSE[k] + USM[IDX(k, j, 3)] * thetaDotSM[k];
                }

                double term2 = 0.0;
                term2 += gradA_SE[0] * dadotdp[j] + gradA_SE[1] * dbEdotdp[j] + gradA_SE[2] * dcEdotdp[j];
                term2 += gradA_SM[0] * dadotdp[j] + gradA_SM[1] * dbMdotdp[j] + gradA_SM[2] * dcMdotdp[j];

                dAdotdp[j] = term1 + term2;
            }

            for (int i = 0; i < 3; ++i) {
                diffGradOut[i] = dcc1dp[i] * Adot + cc1 * dAdotdp[i] + dcc2dp[i] * adot + cc2 * dadotdp[i];
            }
        }
    }

    return 1 - occ;
}

// Parse a command name.
static CMDID parseCmd(const char* s) {
    for (size_t i = 0; i < _nCmds; ++i) {
        if (strcmp(s, _cmds[i].fname) == 0) {
            return _cmds[i].id;
        }
    }
    return -1;
}

// MEX entry point.
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (!isKernelLoaded) {
        const char* path = getenv("SPICE_KERNEL_PATH");
        char fullpath[512];
        SpiceChar errAction[] = "RETURN";
        SpiceChar errReport[] = "NONE";

        erract_c("SET", 0, errAction);
        errprt_c("SET", 0, errReport);
        reset_c();

        if (path == NULL || path[0] == '\0') {
            mexErrMsgIdAndTxt(
                "forcemodel:MissingKernelPath",
                "SPICE_KERNEL_PATH is not set. Call dynamics.ephem or set "
                "the variable to the repository kernels directory."
            );
        }
        for (int i = 0; i < 9; ++i) {
            FILE* kernelFile;
            snprintf(fullpath, sizeof(fullpath), "%s/%s", path, KERNELS[i]);
            kernelFile = fopen(fullpath, "rb");
            if (kernelFile == NULL) {
                mexErrMsgIdAndTxt(
                    "forcemodel:MissingKernel",
                    "Required SPICE kernel not found: %s",
                    fullpath
                );
            }
            fclose(kernelFile);
            furnsh_c(fullpath);
            if (failed_c()) {
                SpiceChar message[1841];
                getmsg_c("LONG", sizeof(message), message);
                reset_c();
                mexErrMsgIdAndTxt(
                    "forcemodel:KernelLoadFailed",
                    "Unable to load SPICE kernel %s: %s",
                    fullpath,
                    message
                );
            }
        }
        isKernelLoaded = true;
    }
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("forcemodel:NotEnoughInputs", "输入参数的数目不足。");
    }
    // Read the command name first.
    char* fname = mxArrayToString(prhs[0]);
    switch (parseCmd(fname)) {
        case STR2ET: { // STR2ET
            // Read inputs.
            const char* date = mxArrayToString(prhs[1]);

            if (date == NULL) {
                mexErrMsgIdAndTxt(
                    "forcemodel:STR2ET:InvalidInput",
                    "转换输入日期时出现错误。"
                );
            }

            {
                SpiceChar errAction[] = "RETURN";
                SpiceChar errReport[] = "NONE";

                erract_c("SET", 0, errAction);
                errprt_c("SET", 0, errReport);
            }

            // Clear any residual CSPICE error state.
            reset_c();

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            double* et = mxGetDoubles(plhs[0]);

            // Call CSPICE.
            str2et_c(date, et);

            // Check for a CSPICE error.
            if (failed_c()) {
                SpiceChar shortMsg[1024];
                SpiceChar longMsg[4096];

                shortMsg[0] = '\0';
                longMsg[0] = '\0';

                getmsg_c("SHORT", sizeof(shortMsg), shortMsg);
                getmsg_c("LONG",  sizeof(longMsg),  longMsg);

                reset_c();

                // Release memory.
                mxFree((void*)date);

                mexErrMsgIdAndTxt(
                    "forcemodel:STR2ET:CSPICEError",
                    "CSPICE str2et_c 失败。\n\nSHORT:\n%s\n\nLONG:\n%s",
                    shortMsg,
                    longMsg
                );
            }

            // Release memory.
            mxFree((void*)date);

            break;
        }
        case SPKPOS: { // SPKPOS
            // Read inputs.
            const char* targ = mxArrayToString(prhs[1]);
            const double et = (double)mxGetScalar(prhs[2]);
            const char* ref = mxArrayToString(prhs[3]);
            const char* abcorr = mxArrayToString(prhs[4]);
            const char* obs = mxArrayToString(prhs[5]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* ptarg = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            double* lt = mxGetDoubles(plhs[1]);

            spkpos_c(targ, et, ref, abcorr, obs, ptarg, lt);

            // Release memory.
            mxFree((void*)targ);
            mxFree((void*)ref);
            mxFree((void*)abcorr);
            mxFree((void*)obs);
            break;
        }
        case SPKEZR: { // SPKEZR
            // Read inputs.
            const char* target = mxArrayToString(prhs[1]);
            const double epoch = (double)mxGetScalar(prhs[2]);
            const char* frame = mxArrayToString(prhs[3]);
            const char* abcorr = mxArrayToString(prhs[4]);
            const char* observer = mxArrayToString(prhs[5]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(6, 1, mxREAL);
            double* state = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            double* lt = mxGetDoubles(plhs[1]);

            spkezr_c(target, epoch, frame, abcorr, observer, state, lt);

            // Release memory.
            mxFree((void*)target);
            mxFree((void*)frame);
            mxFree((void*)abcorr);
            mxFree((void*)observer);
            break;
        }
        case PXFORM: { // PXFORM
            // Read inputs.
            if (nrhs != 4) {
                mexErrMsgIdAndTxt(
                    "forcemodel:PXFORM:InvalidInputCount",
                    "pxform 需要 FROM、TO 和 ET 三个输入参数。"
                );
            }
            if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
                    mxIsSparse(prhs[3])) {
                mexErrMsgIdAndTxt(
                    "forcemodel:PXFORM:InvalidEpoch",
                    "ET 必须为实数 double 标量或向量。"
                );
            }

            char* from = mxArrayToString(prhs[1]);
            char* to = mxArrayToString(prhs[2]);

            const mwSize nET = mxGetNumberOfElements(prhs[3]);
            const double* etv = mxGetDoubles(prhs[3]);

            if (!from || !to) {
                if (from) mxFree(from);
                if (to) mxFree(to);
                mexErrMsgIdAndTxt(
                    "forcemodel:PXFORM:InvalidFrame",
                    "坐标系名称必须为字符向量。"
                );
            }
            if (nET > 1 && (mxGetNumberOfDimensions(prhs[3]) > 2 ||
                    (mxGetM(prhs[3]) != 1 && mxGetN(prhs[3]) != 1))) {
                mxFree(from);
                mxFree(to);
                mexErrMsgIdAndTxt(
                    "forcemodel:PXFORM:InvalidEpoch",
                    "ET 必须为标量或向量。"
                );
            }

            if (nET == 0) {
                // Empty input.
                mwSize dims[3] = {3, 3, 0};
                plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
                mxFree(from);
                mxFree(to);
                break;
            }

            if (nET == 1) {
                // Scalar epoch.
                plhs[0] = mxCreateUninitNumericMatrix(
                    3, 3, mxDOUBLE_CLASS, mxREAL
                );
                double* rotate = mxGetDoubles(plhs[0]);

                SpiceDouble xr[3][3];
                pxform_c(from, to, (SpiceDouble)etv[0], xr);

                xpose_c(xr, (SpiceDouble(*)[3])rotate);
            }
            else {
                // Vector epochs.
                mwSize dims[3] = {3, 3, nET};
                plhs[0] = mxCreateUninitNumericArray(
                    3, dims, mxDOUBLE_CLASS, mxREAL
                );
                double* rotate = mxGetDoubles(plhs[0]);

                SpiceDouble xr[3][3];
                for (mwSize k = 0; k < nET; ++k) {
                    pxform_c(from, to, (SpiceDouble)etv[k], xr);

                    xpose_c(xr, (SpiceDouble(*)[3])(rotate + 9 * k));
                }
            }

            mxFree(from);
            mxFree(to);
            break;
        }
        case SXFORM: { // SXFORM
            // Read inputs.
            if (nrhs != 4) {
                mexErrMsgIdAndTxt(
                    "forcemodel:SXFORM:InvalidInputCount",
                    "sxform 需要 FROM、TO 和 ET 三个输入参数。"
                );
            }
            if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
                    mxIsSparse(prhs[3])) {
                mexErrMsgIdAndTxt(
                    "forcemodel:SXFORM:InvalidEpoch",
                    "ET 必须为实数 double 标量或向量。"
                );
            }

            char* from = mxArrayToString(prhs[1]);
            char* to = mxArrayToString(prhs[2]);
            const mwSize nET = mxGetNumberOfElements(prhs[3]);
            const double* etv = mxGetDoubles(prhs[3]);

            if (!from || !to) {
                if (from) mxFree(from);
                if (to) mxFree(to);
                mexErrMsgIdAndTxt(
                    "forcemodel:SXFORM:InvalidFrame",
                    "坐标系名称必须为字符向量。"
                );
            }
            if (nET > 1 && (mxGetNumberOfDimensions(prhs[3]) > 2 ||
                    (mxGetM(prhs[3]) != 1 && mxGetN(prhs[3]) != 1))) {
                mxFree(from);
                mxFree(to);
                mexErrMsgIdAndTxt(
                    "forcemodel:SXFORM:InvalidEpoch",
                    "ET 必须为标量或向量。"
                );
            }

            if (nET == 0) {
                // Empty input.
                mwSize dims[3] = {6, 6, 0};
                plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
                mxFree(from);
                mxFree(to);
                break;
            }

            if (nET == 1) {
                // Preserve the existing 6-by-6 output for a scalar epoch.
                plhs[0] = mxCreateUninitNumericMatrix(
                    6, 6, mxDOUBLE_CLASS, mxREAL
                );
                double* xform = mxGetDoubles(plhs[0]);

                SpiceDouble xf[6][6];
                sxform_c(from, to, (SpiceDouble)etv[0], xf);
                xpose6_c(xf, (SpiceDouble(*)[6])xform);
            }
            else {
                // Allocate one contiguous 6-by-6-by-N output for vector epochs.
                mwSize dims[3] = {6, 6, nET};
                plhs[0] = mxCreateUninitNumericArray(
                    3, dims, mxDOUBLE_CLASS, mxREAL
                );
                double* xform = mxGetDoubles(plhs[0]);

                SpiceDouble xf[6][6];
                for (mwSize k = 0; k < nET; ++k) {
                    sxform_c(from, to, (SpiceDouble)etv[k], xf);
                    xpose6_c(xf, (SpiceDouble(*)[6])(xform + 36 * k));
                }
            }

            // Release memory.
            mxFree(from);
            mxFree(to);
            break;
        }
        case XF2RAV: { // XF2RAV
            // Read inputs.
            const double* xform = mxGetDoubles(prhs[1]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* rot = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* av = mxGetDoubles(plhs[1]);

            double xf[6][6];
            double xr[3][3];

            xpose6_c(xform, (SpiceDouble(*)[6])xf);
            xf2rav_c(xf, (SpiceDouble(*)[3])xr, av);
            xpose_c(xr, (SpiceDouble(*)[3])rot);
            break;
        }
        case Q2R: { // Q2R
            // Read inputs.
            const double* q = mxGetDoubles(prhs[1]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* R = mxGetDoubles(plhs[0]);

            q2R(q, R);
            break;
        }
        case Q2RJAC: { // Q2RJAC
            // Read inputs.
            const double* q = mxGetDoubles(prhs[1]);

            // Create outputs.
            mwSize dims[3] = {3, 3, 4};
            plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            double* U = mxGetDoubles(plhs[0]);

            q2Rjac(q, U);
            break;
        }
        case GRAVITYFIELD: { // GRAVITYFIELD
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            mwSize maxDegree = (mwSize)mxGetScalar(prhs[2]);
            mwSize maxOrder = (mwSize)mxGetScalar(prhs[3]);
            const double *tm = mxGetDoubles(prhs[4]);
            double gm = mxGetScalar(prhs[5]);
            double req = mxGetScalar(prhs[6]);
            const double *C = mxGetDoubles(prhs[7]);
            const double *S = mxGetDoubles(prhs[8]);
            mwSize coefRows = mxGetM(prhs[7]);
            mxLogical doJac = *mxGetLogicals(prhs[9]);
            mxLogical doHess = *mxGetLogicals(prhs[10]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* accelOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            mwSize dims[3] = {3, 3, 3};
            plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            double* hessOut = mxGetDoubles(plhs[2]);

            gravityfield(posInertial, maxDegree, maxOrder, tm, gm, req, C, S, coefRows, doJac, doHess, accelOut, jacOut, hessOut);
            break;
        }
        case FOURTHBODY: { // FOURTHBODY
            // Read inputs.
            const double epoch = mxGetScalar(prhs[1]);
            const double* posInertial = mxGetDoubles(prhs[2]);
            mxLogical* doBodies = mxGetLogicals(prhs[3]);
            mxLogical doJac = *mxGetLogicals(prhs[4]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* accelOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            fourthbody(epoch, posInertial, doBodies, doJac, accelOut, jacOut);
            break;
        }
        case RELATIVITY: { // RELATIVITY
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* velInertial = mxGetDoubles(prhs[2]);
            mxLogical doJac = *mxGetLogicals(prhs[3]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* accelOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            relativity(posInertial, velInertial, doJac, accelOut, jacOut);
            break;
        }
        case SRP: { // SRP
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* posSun = mxGetDoubles(prhs[2]);
            const double rpcm = mxGetScalar(prhs[3]);
            mxLogical doJac = *mxGetLogicals(prhs[4]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* accelOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            srp(posInertial, posSun, rpcm, doJac, accelOut, jacOut);
            break;
        }
        case SRPNPLATE_Q: { // SRPNPLATE_Q
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* posSun = mxGetDoubles(prhs[2]);
            const double* velInertial = mxGetDoubles(prhs[3]);
            const double* velSun = mxGetDoubles(prhs[4]);
            const double* omegaBody = mxGetDoubles(prhs[5]);
            const double* q = mxGetDoubles(prhs[6]);
            const mxArray* plates = prhs[7];
            const mxLogical doJac = *mxGetLogicals(prhs[8]);
            const mxLogical doJerk = *mxGetLogicals(prhs[9]);
            const mwSize numPlates = mxGetNumberOfElements(plates);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(numPlates, 3, mxREAL);
            double* forceOut = mxGetDoubles(plhs[0]);

            mwSize dims[3] = {3, 13, numPlates};
            plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            plhs[2] = mxCreateDoubleMatrix(numPlates, 3, mxREAL);
            double* jerkOut = mxGetDoubles(plhs[2]);

            plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            double* jerkJacOut = mxGetDoubles(plhs[3]);

            srpnplate_q(posInertial, posSun, velInertial, velSun, omegaBody, q, plates, 
                        doJac, doJerk, numPlates, forceOut, jacOut, jerkOut, jerkJacOut);
            break;
        }
        case SRPNPLATE_R: { // SRPNPLATE_R
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* posSun = mxGetDoubles(prhs[2]);
            const double* velInertial = mxGetDoubles(prhs[3]);
            const double* velSun = mxGetDoubles(prhs[4]);
            const double* tm = mxGetDoubles(prhs[5]);
            const mxArray* plates = prhs[6];
            const mxLogical doJac = *mxGetLogicals(prhs[7]);
            const mwSize numPlates = mxGetNumberOfElements(plates);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(numPlates, 3, mxREAL);
            double* forceOut = mxGetDoubles(plhs[0]);

            mwSize dims[3] = {3, 18, numPlates};
            plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            srpnplate_R(posInertial, posSun, velInertial, velSun, tm, plates, doJac, numPlates, forceOut, jacOut);
            break;
        }
        case EARTHALBEDO: { // EARTHALBEDO
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* posEarth = mxGetDoubles(prhs[2]);
            const double* posSun = mxGetDoubles(prhs[3]);
            const double rpcm = mxGetScalar(prhs[4]);
            mxLogical doJac = *mxGetLogicals(prhs[5]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* accelOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            earthalbedo(posInertial, posEarth, posSun, rpcm, doJac, accelOut, jacOut);
            break;
        }
        case JERKGRAVITYFIELD: { // JERKGRAVITYFIELD
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* velInertial = mxGetDoubles(prhs[2]);
            const double* accelInertial = mxGetDoubles(prhs[3]);
            const double* omegaInertial = mxGetDoubles(prhs[4]);
            const double* jacGravityfield = mxGetDoubles(prhs[5]);
            
            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* jerkOut = mxGetDoubles(plhs[0]);

            jerkGravityfield(posInertial, velInertial, accelInertial, omegaInertial, jacGravityfield, jerkOut);
            break;
        }
        case JERKFOURTHBODY: { // JERKFOURTHBODY
            // Read inputs.
            const double epoch = mxGetScalar(prhs[1]);
            const double* posInertial = mxGetDoubles(prhs[2]);
            const double* velInertial = mxGetDoubles(prhs[3]);
            mxLogical* doBodies = mxGetLogicals(prhs[4]);
            mxLogical* doJac = mxGetLogicals(prhs[5]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* jerkOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            jerkFourthbody(epoch, posInertial, velInertial, doBodies, doJac, jerkOut, jacOut);
            break;
        }
        case JERKRELATIVITY: { // JERKRELATIVITY
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* velInertial = mxGetDoubles(prhs[2]);
            const double* accelInertial = mxGetDoubles(prhs[3]);
            const double* jacRelativity = mxGetDoubles(prhs[4]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* jerkOut = mxGetDoubles(plhs[0]);

            jerkRelativity(posInertial, velInertial, accelInertial, jacRelativity, jerkOut);
            break;
        }
        case JERKSRP: { // JERKSRP
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* velInertial = mxGetDoubles(prhs[2]);
            const double* posSun = mxGetDoubles(prhs[3]);
            const double* velSun = mxGetDoubles(prhs[4]);
            const double rpcm = mxGetScalar(prhs[5]);
            const mxLogical doJac = *mxGetLogicals(prhs[6]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* jerkOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            jerkSrp(posInertial, velInertial, posSun, velSun, rpcm, doJac, jerkOut, jacOut);
            break;
        }
        case JERKEARTHALBEDO: { // JERKEARTHALBEDO
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* velInertial = mxGetDoubles(prhs[2]);
            const double* posEarth = mxGetDoubles(prhs[3]);
            const double* velEarth = mxGetDoubles(prhs[4]);
            const double* posSun = mxGetDoubles(prhs[5]);
            const double* velSun = mxGetDoubles(prhs[6]);
            const double rpcm = mxGetScalar(prhs[7]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* jerkOut = mxGetDoubles(plhs[0]);

            jerkEarthalbedo(posInertial, velInertial, posEarth, velEarth, posSun, velSun, rpcm, jerkOut);
            break;
        }
        case MAGNETICFIELD: { // MAGNETICFIELD
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            mwSize maxDegree = (mwSize)mxGetScalar(prhs[2]);
            mwSize maxOrder = (mwSize)mxGetScalar(prhs[3]);
            const double *tm = mxGetDoubles(prhs[4]);
            double req = mxGetScalar(prhs[5]);
            const double *G = mxGetDoubles(prhs[6]);
            const double *H = mxGetDoubles(prhs[7]);
            mwSize coefRows = mxGetM(prhs[6]);
            mxLogical doJac = *mxGetLogicals(prhs[8]);

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* fluxOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 6, mxREAL);
            double* jacOut = mxGetDoubles(plhs[1]);

            magneticfield(posInertial, maxDegree, maxOrder, tm, req, G, H, coefRows, doJac, fluxOut, jacOut);
            break;
        }
        case GRAVITYTORQUE_Q: { // GRAVITYTORQUE_Q
            // Read inputs.
            const double* jacInertial = mxGetDoubles(prhs[1]);
            const double* q = mxGetDoubles(prhs[2]);
            const double* inertiaMatrix = mxGetDoubles(prhs[3]);
            const mxLogical doJac = *mxGetLogicals(prhs[4]);
            double* hessInertial = NULL;
            if (doJac && nrhs > 5) {
                hessInertial = mxGetDoubles(prhs[5]);
            }

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* torqueOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* jacOut_r = mxGetDoubles(plhs[1]);

            plhs[2] = mxCreateDoubleMatrix(3, 4, mxREAL);
            double* jacOut_q = mxGetDoubles(plhs[2]);

            gravitytorque_q(jacInertial, q, inertiaMatrix, doJac, hessInertial, torqueOut, jacOut_r, jacOut_q);
            break;
        }
        case GRAVITYTORQUE_R: { // GRAVITYTORQUE_R
            // Read inputs.
            const double* jacInertial = mxGetDoubles(prhs[1]);
            const double* tm = mxGetDoubles(prhs[2]);
            const double* inertiaMatrix = mxGetDoubles(prhs[3]);
            const mxLogical doJac = *mxGetLogicals(prhs[4]);
            double* hessInertial = NULL;
            if (doJac && nrhs > 5) {
                hessInertial = mxGetDoubles(prhs[5]);
            }

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* torqueOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* jacOut_r = mxGetDoubles(plhs[1]);

            plhs[2] = mxCreateDoubleMatrix(3, 9, mxREAL);
            double* jacOut_R = mxGetDoubles(plhs[2]);

            gravitytorque_R(jacInertial, tm, inertiaMatrix, doJac, hessInertial, torqueOut, jacOut_r, jacOut_R);
            break;
        }
        case SRPNPLATETORQUE_Q: { // SRPNPLATETORQUE_Q
            // Read inputs.
            const double* forceInertial = mxGetDoubles(prhs[1]);
            const double* q = mxGetDoubles(prhs[2]);
            const double* displacements = mxGetDoubles(prhs[3]);
            const mwSize numPlates = mxGetM(prhs[3]);
            const mxLogical doJac = *mxGetLogicals(prhs[4]);
            double* jacInertial = NULL;
            if (doJac && nrhs > 5) {
                jacInertial = mxGetDoubles(prhs[5]);
            }

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* torqueOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* jacOut_r = mxGetDoubles(plhs[1]);

            plhs[2] = mxCreateDoubleMatrix(3, 4, mxREAL);
            double* jacOut_q = mxGetDoubles(plhs[2]);

            srpnplatetorque_q(forceInertial, q, displacements, doJac, jacInertial, numPlates, torqueOut, jacOut_r, jacOut_q);
            break;
        }
        case SRPNPLATETORQUE_R: { // SRPNPLATETORQUE_R
            // Read inputs.
            const double* forceInertial = mxGetDoubles(prhs[1]);
            const double* tm = mxGetDoubles(prhs[2]);
            const double* displacements = mxGetDoubles(prhs[3]);
            const mwSize numPlates = mxGetM(prhs[3]);
            const mxLogical doJac = *mxGetLogicals(prhs[4]);
            double* jacInertial = NULL;
            if (doJac && nrhs > 5) {
                jacInertial = mxGetDoubles(prhs[5]);
            }

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* torqueOut = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* jacOut_r = mxGetDoubles(plhs[1]);

            plhs[2] = mxCreateDoubleMatrix(3, 9, mxREAL);
            double* jacOut_R = mxGetDoubles(plhs[2]);

            srpnplatetorque_R(forceInertial, tm, displacements, doJac, jacInertial, numPlates, torqueOut, jacOut_r, jacOut_R);
            break;
        }
        case DUALCONE: { // DUALCONE
            // Read inputs.
            const double* posInertial = mxGetDoubles(prhs[1]);
            const double* posSun = mxGetDoubles(prhs[2]);
            const double* posEarth = mxGetDoubles(prhs[3]);
            double* velInertial = NULL;
            double* velSun = NULL;
            double* velEarth = NULL;
            mxLogical doJac = false;
            if (nlhs > 2 && nrhs > 4) {
                velInertial = mxGetDoubles(prhs[4]);
                velSun = mxGetDoubles(prhs[5]);
                velEarth = mxGetDoubles(prhs[6]);
                doJac = *mxGetLogicals(prhs[7]);
            }

            // Create outputs.
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            double* nu = mxGetDoubles(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL);
            double* gradOut = mxGetDoubles(plhs[1]);

            plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
            double* diffOut = mxGetDoubles(plhs[2]);
            
            plhs[3] = mxCreateDoubleMatrix(1, 6, mxREAL);
            double* diffGradOut = mxGetDoubles(plhs[3]);

            *nu = dualcone(posInertial, posSun, posEarth, velInertial, velSun, velEarth, doJac, gradOut, diffOut, diffGradOut);
            break;
        }
        default: {
            mxFree((void*)fname);
            mexErrMsgIdAndTxt("forcemodel:UnknownFuncName", "未知的函数名称。");
        }
    }
    mxFree((void*)fname);
}
