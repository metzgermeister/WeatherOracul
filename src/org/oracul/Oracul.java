package org.oracul;

import static java.lang.Math.*;
import static org.oracul.OraculConstants.*;
import static org.oracul.Utils.fRo;
import static org.oracul.Utils.quadrInterpol;

import org.oracul.data.ForecastInfo;
import org.oracul.data.OraculData;

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
            // SetFU(data.U, data.V, data.H, data.F1, data.F2);
            System.arraycopy(data.U, 0, data.C, 0, SZ_LAM * SZ_FI);
            // SplitLam(data.U, data.C, data.Mu1, data.F1, data.U1);

            System.arraycopy(data.V, 0, data.C, 0, SZ_LAM * SZ_FI);
            // SplitFi(data.U, data.C, data.Mu2, data.F2, data.U2);

            // Calc v
            // SetFV(data.U, data.V, data.H, data.F1, data.F2);

            System.arraycopy(data.U, 0, data.C, 0, SZ_LAM * SZ_FI);
            // SplitLam(data.V, data.C, data.Mu1, data.F1, data.V1);

            System.arraycopy(data.V, 0, data.C, 0, SZ_LAM * SZ_FI);
            // SplitFi(data.V, data.C, data.Mu2, data.F2, data.V2);

            // Final
            // SplitSum(data.U, data.U1, data.U2);
            // SplitSum(data.V, data.V1, data.V2);
        }

        // CalcNorma(2. * (i - 1) * dt / T);
        // SaveToFile('u', data.U);
        // SaveToFile('v', data.V);
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
}
