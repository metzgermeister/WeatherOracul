package org.oracul;

import static java.lang.Math.*;
import static org.oracul.OraculConstants.*;
import static org.oracul.Utils.*;


import org.oracul.data.ForecastInfo;
import org.oracul.data.OraculData;

import java.io.File;

/**
 * User: Pavlo_Ivanenko Date: 9/20/13 Time: 3:01 PM
 */
public class Oracul {
    private final OraculData data;

    public Oracul(OraculData data) {
        this.data = data;
    }

    public void predict(ForecastInfo forecastInfo) {
        int i;

        for (i = 0; (i * forecastInfo.getTimeStep()) <= (forecastInfo.getTimeLimit() - forecastInfo.getTimeStep()); i++) {
            setMu(data.U, data.V, data.Mu1, data.Mu2);
            setBoundaryCnd(2. * (i + 1) * dt / T, data.U1, data.U2, data.V1, data.V2);
            setActualH(2. * i * dt / T, data.H);

            // Calc u
            setFU(data.U, data.V, data.H, data.F1, data.F2);
            System.arraycopy(data.U, 0, data.C, 0, SZ_LAM * SZ_FI);
            splitLam(data.U, data.C, data.Mu1, data.F1, data.U1);

            System.arraycopy(data.V, 0, data.C, 0, SZ_LAM * SZ_FI);
            splitFi(data.U, data.C, data.Mu2, data.F2, data.U2);

            // Calc v
            // SetFV(data.U, data.V, data.H, data.F1, data.F2);

            System.arraycopy(data.U, 0, data.C, 0, SZ_LAM * SZ_FI);
            splitLam(data.V, data.C, data.Mu1, data.F1, data.V1);

            System.arraycopy(data.V, 0, data.C, 0, SZ_LAM * SZ_FI);
            splitFi(data.V, data.C, data.Mu2, data.F2, data.V2);

            // Final
            // SplitSum(data.U, data.U1, data.U2);
            // SplitSum(data.V, data.V1, data.V2);
        }

        // CalcNorma(2. * (i - 1) * dt / T);
    }

    private void setFU(double[] UU, double[] VV, double[] HH, double[] FF1, double[] FF2) {
        int i, j;
        double cs, sn, g;

        for (j = 1; j <= (SZ_FI - 2); j++) {
            cs = ERH_RDS * cos(j * dh);
            sn = dt * sin(j * dh);
            g = dt * fG(SRF_ALT, j * dh) / (2. * cs * dh);
            for (i = 1; i <= (SZ_LAM - 2); i++) {
                FF1[j * SZ_LAM + i] = -g * (HH[j * SZ_LAM + i + 1] - HH[j * SZ_LAM + i - 1]);
                FF2[j * SZ_LAM + i] = sn * VV[j * SZ_LAM + i] * (UU[j * SZ_LAM + i] / cs +
                        2. * ERH_OMEGA);
            }
        }

    }


    private void setMu(double[] UU, double[] VV, double[] M1, double[] M2) {
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

    private void setBoundaryCnd(double r, double[] UU1, double[] UU2, double[] VV1, double[] VV2) {
        int i, j;
        // lam
        for (i = 0; i <= (SZ_FI - 1); i++)
            for (j = 0; j <= (SZ_LAM - 1); j += (SZ_LAM - 1))
                UU1[i * SZ_LAM + j] = quadrInterpol(r, data.lamBndCnd[0][0][j / (SZ_LAM - 1)][i],
                        data.lamBndCnd[0][1][j / (SZ_LAM - 1)][i], data.lamBndCnd[0][2][j / (SZ_LAM - 1)][i]);
        for (i = 0; i <= (SZ_FI - 1); i++)
            for (j = 0; j <= (SZ_LAM - 1); j += (SZ_LAM - 1))
                VV1[i * SZ_LAM + j] = quadrInterpol(r, data.lamBndCnd[1][0][j / (SZ_LAM - 1)][i],
                        data.lamBndCnd[1][1][j / (SZ_LAM - 1)][i], data.lamBndCnd[1][2][j / (SZ_LAM - 1)][i]);
        // fi
        for (j = 0; j <= (SZ_LAM - 1); j++)
            for (i = 0; i <= (SZ_FI - 1); i += (SZ_FI - 1))
                UU2[i * SZ_LAM + j] = quadrInterpol(r, data.fiBndCnd[0][0][i / (SZ_FI - 1)][j], data.fiBndCnd[0][1][i
                        / (SZ_FI - 1)][j], data.fiBndCnd[0][2][i / (SZ_FI - 1)][j]);
        for (j = 0; j <= (SZ_LAM - 1); j++)
            for (i = 0; i <= (SZ_FI - 1); i += (SZ_FI - 1))
                VV2[i * SZ_LAM + j] = quadrInterpol(r, data.fiBndCnd[1][0][i / (SZ_FI - 1)][j], data.fiBndCnd[1][1][i
                        / (SZ_FI - 1)][j], data.fiBndCnd[1][2][i / (SZ_FI - 1)][j]);
    }

    private void setActualH(double r, double[] HH) {
        for (int i = 0; i <= (SZ_FI * SZ_LAM - 1); i++)
            HH[i] = quadrInterpol(r, data.H1[i], data.H2[i], data.H3[i]);

    }


    private void splitLam(double[] R, double[] CC, double[] MMu, double[] FF1, double[] R1) {
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
                        qm = R[j * SZ_LAM + i - 1] + FF1[j * SZ_LAM + i - 1] + qm * (R[j * SZ_LAM + i - 2] - R[j * SZ_LAM + i - 1]);
                        p = pm * q + p * qm;

                        R1[j * SZ_LAM + i - 1] = (qm + p) / g;
                        R1[j * SZ_LAM + i] = (q + p) / g;

                        minusX1(beg, i - 2, j, R, CC, MMu, FF1, R1);
                    } else {// "+", "+" (pCC[i]>0)
                        R1[j * SZ_LAM + i] = R1[j * SZ_LAM + i - 1] + (R[j * SZ_LAM + i] - R1[j * SZ_LAM + i - 1]
                                + q * (R[j * SZ_LAM + i] - R[j * SZ_LAM + i + 1]) + FF1[j * SZ_LAM + i]) / (1. + p);
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


    void minusX1(int i1, int i2, int j, double[] R, double[] CC,
                 double[] MMu, double[] FF1, double[] R1) {
        double p, q, z = dt / (2. * dh * dh * ERH_RDS_2 * cos(j * dh) * cos(j * dh)),
                h = dh * ERH_RDS * cos(j * dh);

        for (int i = i2; i >= i1; i--) {
            p = z * (-h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[j * SZ_LAM + i + 1]);
            q = z * (h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[j * SZ_LAM + i - 1]);

            R1[j * SZ_LAM + i] = R1[j * SZ_LAM + i + 1] + (R[j * SZ_LAM + i] - R1[j * SZ_LAM + i + 1]
                    + q * (R[j * SZ_LAM + i - 1] - R[j * SZ_LAM + i]) + FF1[j * SZ_LAM + i]) / (1. + p);
        }
    }


    private void splitFi(double[] R, double[] CC, double[] MMu, double[] FF2, double[] R2) {
        boolean flg1, flg2;
        int i, j, beg = 0;
        double p, q, pm, qm, g, z = dt / (2 * dh * dh * ERH_RDS_2),
                h = dh * ERH_RDS;

        for (i = 1; i <= (SZ_LAM - 2); i++) {
            flg1 = true;
            for (j = 1; j <= (SZ_FI - 2); j++) {
                if (CC[j * SZ_LAM + i] > 0.) {
                    flg2 = true;
                } else if (CC[j * SZ_LAM + i] == 0.) flg2 = flg1;
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
                        qm = R[(j - 1) * SZ_LAM + i] + FF2[(j - 1) * SZ_LAM + i] + qm * (R[(j - 2) * SZ_LAM + i] - R[(j - 1) * SZ_LAM + i]);
                        p = pm * q + p * qm;

                        R2[(j - 1) * SZ_LAM + i] = (qm + p) / g;
                        R2[j * SZ_LAM + i] = (q + p) / g;

                        minusX2(beg, j - 2, i, R, CC, MMu, FF2, R2);
                    } else // "+", "+" (pCC[j]>0)
                    {
                        R2[j * SZ_LAM + i] = R2[(j - 1) * SZ_LAM + i] + (R[j * SZ_LAM + i] - R2[(j - 1) * SZ_LAM + i]
                                + q * (R[j * SZ_LAM + i] - R[(j + 1) * SZ_LAM + i]) + FF2[j * SZ_LAM + i]) / (1. + p);
                    }
                } else if (flg1) beg = j; // from "+" (j-1) to "-" (j)

                flg1 = flg2;
            }

            if (!flg1) {
                minusX2(beg, SZ_FI - 2, i, R, CC, MMu, FF2, R2);
            }
        }
    }


    private void minusX2(int j1, int j2, int i, double[] R, double[] CC,
                         double[] MMu, double[] FF2, double[] R2) {
        int j;
        double p, q, z = dt / (2. * dh * dh * ERH_RDS_2),
                h = dh * ERH_RDS;

        for (j = j2; j >= j1; j--) {
            p = z * (-h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[(j + 1) * SZ_LAM + i]);
            q = z * (h * CC[j * SZ_LAM + i] + MMu[j * SZ_LAM + i] + MMu[(j - 1) * SZ_LAM + i]);

            R2[j * SZ_LAM + i] = R2[(j + 1) * SZ_LAM + i] + (R[j * SZ_LAM + i] - R2[(j + 1) * SZ_LAM + i]
                    + q * (R[(j - 1) * SZ_LAM + i] - R[j * SZ_LAM + i]) + FF2[j * SZ_LAM + i]) / (1. + p);
        }
    }


}
