package org.oracul;

import static java.lang.Math.*;
import static org.oracul.util.OraculConstants.*;
import static org.oracul.util.Utils.*;

import org.oracul.data.OraculData;

/**
 * User: Pavlo_Ivanenko Date: 9/23/13 Time: 1:52 PM
 */
public abstract class AbstractOracul implements Oracul {

    private final OraculData data;

    public AbstractOracul(OraculData data) {
        this.data = data;
    }

    public OraculData getData() {
        return data;
    }

    protected void setMu(double[] UU, double[] VV, double[] M1, double[] M2) {
        int i, j;
        double dd, dc, tg, cs, ss = S / fRo(SRF_ALT), k0 = 50000.;

        for (j = 1; j <= (SZ_FI - 2); j++) {
            tg = tan(j * dh);
            cs = dh * cos(j * dh);
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                dc = (UU[j * SZ_LAM + i + 1] - UU[j * SZ_LAM + i - 1]) / (2. * cs)
                        - (VV[(j + 1) * SZ_LAM + i] - VV[(j - 1) * SZ_LAM + i]) / (2. * dh) - VV[j * SZ_LAM + i] * tg;
                dd = (VV[j * SZ_LAM + i + 1] - VV[j * SZ_LAM + i - 1]) / (2. * cs)
                        + (UU[(j + 1) * SZ_LAM + i] - UU[(j - 1) * SZ_LAM + i]) / (2. * dh) + UU[j * SZ_LAM + i] * tg;
                M2[j * SZ_LAM + i] = ss * (k0 + 0.08 * ERH_RDS * (cs * cs + dh * dh) * sqrt(dd * dd + dc * dc));
                M1[j * SZ_LAM + i] = M2[j * SZ_LAM + i];
            }
        }
    }

    protected void setBoundaryCnd(double r, double[] UU1, double[] UU2, double[] VV1, double[] VV2) {
        int i, j;
        // lam
        for (i = 0; i <= (SZ_FI - 1); i++)
            for (j = 0; j <= (SZ_LAM - 1); j += (SZ_LAM - 1))
                UU1[i * SZ_LAM + j] = quadrInterpol(r, getData().lamBndCnd[0][0][j / (SZ_LAM - 1)][i],
                        getData().lamBndCnd[0][1][j / (SZ_LAM - 1)][i], getData().lamBndCnd[0][2][j / (SZ_LAM - 1)][i]);
        for (i = 0; i <= (SZ_FI - 1); i++)
            for (j = 0; j <= (SZ_LAM - 1); j += (SZ_LAM - 1))
                VV1[i * SZ_LAM + j] = quadrInterpol(r, getData().lamBndCnd[1][0][j / (SZ_LAM - 1)][i],
                        getData().lamBndCnd[1][1][j / (SZ_LAM - 1)][i], getData().lamBndCnd[1][2][j / (SZ_LAM - 1)][i]);
        // fi
        for (j = 0; j <= (SZ_LAM - 1); j++)
            for (i = 0; i <= (SZ_FI - 1); i += (SZ_FI - 1))
                UU2[i * SZ_LAM + j] = quadrInterpol(r, getData().fiBndCnd[0][0][i / (SZ_FI - 1)][j],
                        getData().fiBndCnd[0][1][i / (SZ_FI - 1)][j], getData().fiBndCnd[0][2][i / (SZ_FI - 1)][j]);
        for (j = 0; j <= (SZ_LAM - 1); j++)
            for (i = 0; i <= (SZ_FI - 1); i += (SZ_FI - 1))
                VV2[i * SZ_LAM + j] = quadrInterpol(r, getData().fiBndCnd[1][0][i / (SZ_FI - 1)][j],
                        getData().fiBndCnd[1][1][i / (SZ_FI - 1)][j], getData().fiBndCnd[1][2][i / (SZ_FI - 1)][j]);
    }

    protected void setActualH(double r, double[] HH) {
        for (int i = 0; i <= (SZ_FI * SZ_LAM - 1); i++)
            HH[i] = quadrInterpol(r, getData().H1[i], getData().H2[i], getData().H3[i]);

    }

    protected void splitLam(double[] R, double[] CC, double[] MMu, double[] FF1, double[] R1) {
        boolean flg1, flg2;
        int i, j, beg = 0;
        double p, q, pm, qm, g, z, h;

        for (j = 1; j <= (SZ_FI - 2); j++) {
            flg1 = true;
            z = dt / (2 * dh * dh * ERH_RDS_2 * cos(j * dh) * cos(j * dh));
            h = dh * ERH_RDS * cos(j * dh);
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                if (CC[j * SZ_LAM + i] > 0.) {
                    flg2 = true;
                } else if (CC[j * SZ_LAM + i] == 0.) {
                    flg2 = flg1;
                } else {
                    flg2 = false;
                }

                if (flg2) {
                    p = z * (h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[j * SZ_LAM + i - 1]);
                    q = z * (h * CC[j * SZ_LAM + i] - MMu[j * SZ_LAM + i] - MMu[j * SZ_LAM + i + 1]);

                    if (!flg1) {// from "-" (i-1) to "+" (i)
                        pm = z * (-h * CC[j * SZ_LAM + i - 1] + MMu[j * SZ_LAM + i - 1] + MMu[j * SZ_LAM + i]);
                        qm = z * (h * CC[j * SZ_LAM + i - 1] + MMu[j * SZ_LAM + i - 1] + MMu[j * SZ_LAM + i - 2]);

                        g = 1. + p + pm;
                        q = R[j * SZ_LAM + i] + FF1[j * SZ_LAM + i] + q * (R[j * SZ_LAM + i] - R[j * SZ_LAM + i + 1]);
                        qm = R[j * SZ_LAM + i - 1] + FF1[j * SZ_LAM + i - 1] + qm
                                * (R[j * SZ_LAM + i - 2] - R[j * SZ_LAM + i - 1]);
                        p = pm * q + p * qm;

                        R1[j * SZ_LAM + i - 1] = (qm + p) / g;
                        R1[j * SZ_LAM + i] = (q + p) / g;

                        minusX1(beg, i - 2, j, R, CC, MMu, FF1, R1);
                    } else {// "+", "+" (pCC[i]>0)
                        R1[j * SZ_LAM + i] = R1[j * SZ_LAM + i - 1]
                                + (R[j * SZ_LAM + i] - R1[j * SZ_LAM + i - 1] + q
                                        * (R[j * SZ_LAM + i] - R[j * SZ_LAM + i + 1]) + FF1[j * SZ_LAM + i]) / (1. + p);
                    }
                } else if (flg1) {
                    beg = i;// from "+" (i-1) to "-" (i)
                }

                flg1 = flg2;
            }

            if (!flg1) {
                minusX1(beg, SZ_LAM - 2, j, R, CC, MMu, FF1, R1);
            }

        }
    }

    protected void minusX1(int i1, int i2, int j, double[] R, double[] CC, double[] MMu, double[] FF1, double[] R1) {
        double p, q, z = dt / (2. * dh * dh * ERH_RDS_2 * cos(j * dh) * cos(j * dh)), h = dh * ERH_RDS * cos(j * dh);

        for (int i = i2; i >= i1; i--) {
            p = z * (-h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[j * SZ_LAM + i + 1]);
            q = z * (h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[j * SZ_LAM + i - 1]);

            R1[j * SZ_LAM + i] = R1[j * SZ_LAM + i + 1]
                    + (R[j * SZ_LAM + i] - R1[j * SZ_LAM + i + 1] + q * (R[j * SZ_LAM + i - 1] - R[j * SZ_LAM + i]) + FF1[j
                            * SZ_LAM + i]) / (1. + p);
        }
    }

    protected void splitFi(double[] R, double[] CC, double[] MMu, double[] FF2, double[] R2) {
        boolean flg1, flg2;
        int i, j, beg = 0;
        double p, q, pm, qm, g, z = dt / (2 * dh * dh * ERH_RDS_2), h = dh * ERH_RDS;

        for (i = 1; i <= (SZ_LAM - 2); i++) {
            flg1 = true;
            for (j = 1; j <= (SZ_FI - 2); j++) {
                if (CC[j * SZ_LAM + i] > 0.) {
                    flg2 = true;
                } else if (CC[j * SZ_LAM + i] == 0.)
                    flg2 = flg1;
                else {
                    flg2 = false;
                }

                if (flg2) {
                    p = z * (h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[(j - 1) * SZ_LAM + i]);
                    q = z * (h * CC[j * SZ_LAM + i] - MMu[j * SZ_LAM + i] - MMu[(j + 1) * SZ_LAM + i]);

                    if (!flg1) {// from "-" (j-1) to "+" (j)
                        pm = z * (-h * CC[(j - 1) * SZ_LAM + i] + MMu[(j - 1) * SZ_LAM + i] + MMu[j * SZ_LAM + i]);
                        qm = z * (h * CC[(j - 1) * SZ_LAM + i] + MMu[(j - 1) * SZ_LAM + i] + MMu[(j - 2) * SZ_LAM + i]);

                        g = 1. + p + pm;
                        q = R[j * SZ_LAM + i] + FF2[j * SZ_LAM + i] + q * (R[j * SZ_LAM + i] - R[(j + 1) * SZ_LAM + i]);
                        qm = R[(j - 1) * SZ_LAM + i] + FF2[(j - 1) * SZ_LAM + i] + qm
                                * (R[(j - 2) * SZ_LAM + i] - R[(j - 1) * SZ_LAM + i]);
                        p = pm * q + p * qm;

                        R2[(j - 1) * SZ_LAM + i] = (qm + p) / g;
                        R2[j * SZ_LAM + i] = (q + p) / g;

                        minusX2(beg, j - 2, i, R, CC, MMu, FF2, R2);
                    } else // "+", "+" (pCC[j]>0)
                    {
                        R2[j * SZ_LAM + i] = R2[(j - 1) * SZ_LAM + i]
                                + (R[j * SZ_LAM + i] - R2[(j - 1) * SZ_LAM + i] + q
                                        * (R[j * SZ_LAM + i] - R[(j + 1) * SZ_LAM + i]) + FF2[j * SZ_LAM + i])
                                / (1. + p);
                    }
                } else if (flg1)
                    beg = j; // from "+" (j-1) to "-" (j)

                flg1 = flg2;
            }

            if (!flg1) {
                minusX2(beg, SZ_FI - 2, i, R, CC, MMu, FF2, R2);
            }
        }
    }

    protected void minusX2(int j1, int j2, int i, double[] R, double[] CC, double[] MMu, double[] FF2, double[] R2) {
        int j;
        double p, q, z = dt / (2. * dh * dh * ERH_RDS_2), h = dh * ERH_RDS;

        for (j = j2; j >= j1; j--) {
            p = z * (-h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[(j + 1) * SZ_LAM + i]);
            q = z * (h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[(j - 1) * SZ_LAM + i]);

            R2[j * SZ_LAM + i] = R2[(j + 1) * SZ_LAM + i]
                    + (R[j * SZ_LAM + i] - R2[(j + 1) * SZ_LAM + i] + q * (R[(j - 1) * SZ_LAM + i] - R[j * SZ_LAM + i]) + FF2[j
                            * SZ_LAM + i]) / (1. + p);
        }
    }

    protected void setFV(double[] UU, double[] VV, double[] HH, double[] FF1, double[] FF2) {
        int i, j;
        double cs, sn, g;

        for (j = 1; j <= (SZ_FI - 2); j++) {
            cs = ERH_RDS * cos(j * dh);
            sn = dt * sin(j * dh);
            g = dt * fG(SRF_ALT, j * dh) / (2. * dh * ERH_RDS);
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                // TODO maybe here is mistake and VV data should be used
                FF1[j * SZ_LAM + i] = -sn * UU[j * SZ_LAM + i] * (UU[j * SZ_LAM + i] / cs + 2. * ERH_OMEGA);
                FF2[j * SZ_LAM + i] = -g * (HH[(j + 1) * SZ_LAM + i] - HH[(j - 1) * SZ_LAM + i]);
            }
        }
    }

    protected void setFU(double[] UU, double[] VV, double[] HH, double[] FF1, double[] FF2) {
        int i, j;
        double cs, sn, g;

        for (j = 1; j <= (SZ_FI - 2); j++) {
            cs = ERH_RDS * cos(j * dh);
            sn = dt * sin(j * dh);
            g = dt * fG(SRF_ALT, j * dh) / (2. * cs * dh);
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                FF1[j * SZ_LAM + i] = -g * (HH[j * SZ_LAM + i + 1] - HH[j * SZ_LAM + i - 1]);
                FF2[j * SZ_LAM + i] = sn * VV[j * SZ_LAM + i] * (UU[j * SZ_LAM + i] / cs + 2. * ERH_OMEGA);
            }
        }

    }

    protected void splitSum(double[] UU, double[] UU1, double[] UU2) {
        int i, j;

        for (i = 1; i <= (SZ_FI - 2); i++)
            for (j = 1; j <= (SZ_LAM - 2); j++)
                UU[i * SZ_LAM + j] = UU1[i * SZ_LAM + j] + UU2[i * SZ_LAM + j] - UU[i * SZ_LAM + j];

        // Boundary condition
        for (j = 0; j <= (SZ_LAM - 1); j += (SZ_LAM - 1))
            for (i = 0; i <= (SZ_FI - 1); i++)
                UU[i * SZ_LAM + j] = UU1[i * SZ_LAM + j];

        for (i = 0; i <= (SZ_FI - 1); i += (SZ_FI - 1))
            for (j = 0; j <= (SZ_LAM - 1); j++)
                UU[i * SZ_LAM + j] = UU2[i * SZ_LAM + j];

    }

    protected void сalcNorma(double r, double[] pU2, double[] pV2, double[] pU, double[] pV) {
        int i, j;
        double maxU = 0., L2U = 0., maxV = 0., L2V = 0., z;

        loadOriginalValues(r, pU2, pV2);
        // For u
        for (j = 1; j <= (SZ_FI - 2); j++)
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                z = pU2[j * SZ_LAM + i] - pU[j * SZ_LAM + i];
                z *= z;
                L2U += z;
                if (z > maxU)
                    maxU = z;
            }
        maxU = sqrt(maxU);
        L2U = dh * sqrt(L2U);

        // For v
        for (j = 1; j <= (SZ_FI - 2); j++)
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                z = pV2[j * SZ_LAM + i] - pV[j * SZ_LAM + i];
                z *= z;
                L2V += z;
                if (z > maxV)
                    maxV = z;
            }
        maxV = sqrt(maxV);
        L2V = dh * sqrt(L2V);

        String statistics = String.format("maxU=%9.2e L2U=%9.2e\nmaxV=%9.2e L2V=%9.2e\n", maxU, L2U, maxV, L2V);
        System.out.println(statistics);

    }

    private void loadOriginalValues(double r, double[] pUU1, double[] pVV1) {
        // char buf[];
        // int i, j, k;
        // double[][] pArOut = new double[2][];
        // pArOut[0] = pUU1;
        // pArOut[1] = pVV1;
        // double[] z = new double[3];
        // File * pFile[3];
        //
        // for (i = 0; i <= 1; i++) {
        // for (k = 0; k <= 2; k++) {
        // sprintf(buf, "IN\\%s", flNames[3 * k + 1 + i]);
        // pFile[k] = fopen(buf, "r");
        // }
        //
        // for (j = 0; j <= (SZ_FI * SZ_LAM - 1); j++) {
        // for (k = 0; k <= 2; k++) {
        // fgets(buf, 64, pFile[k]);
        // z[k] = atof(buf);
        // }
        // pArOut[i][j] = quadrInterpol(r, z[0], z[1], z[2]);
        // }
        //
        // for (k = 0; k <= 2; k++) {
        // fclose(pFile[k]);
        // }
        // }
        //
    }

}
