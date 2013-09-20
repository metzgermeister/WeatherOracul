package org.oracul;

import static org.oracul.OraculConstants.SZ_FI;
import static org.oracul.OraculConstants.SZ_LAM;

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

            // SetMu(data.U, data.V, data.Mu1, data.Mu2);
            // SetBoundaryCnd(2. * (i + 1) * dt / T, data.U1, data.U2, data.V1, data.V2);
            // SetActualH(2. * i * dt / T, data.H);

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
}
