package org.oracul;

import static org.oracul.util.OraculConstants.*;

import org.oracul.data.ForecastInfo;
import org.oracul.data.OraculData;

/**
 * User: Pavlo_Ivanenko Date: 9/20/13 Time: 3:01 PM
 */
public class SequentialOracul extends AbstractOracul {

    public SequentialOracul(OraculData data) {
        super(data);
    }

    @Override
    public void predict(ForecastInfo forecastInfo) {
        int i;

        for (i = 0; (i * forecastInfo.getTimeStep()) <= (forecastInfo.getTimeLimit() - forecastInfo.getTimeStep()); i++) {
            setMu(getData().U, getData().V, getData().Mu1, getData().Mu2);
            setBoundaryCnd(2. * (i + 1) * dt / T, getData().U1, getData().U2, getData().V1, getData().V2);
            setActualH(2. * i * dt / T, getData().H);

            // Calc u
            setFU(getData().U, getData().V, getData().H, getData().F1, getData().F2);
            System.arraycopy(getData().U, 0, getData().C, 0, SZ_LAM * SZ_FI);
            splitLam(getData().U, getData().C, getData().Mu1, getData().F1, getData().U1);

            System.arraycopy(getData().V, 0, getData().C, 0, SZ_LAM * SZ_FI);
            splitFi(getData().U, getData().C, getData().Mu2, getData().F2, getData().U2);

            // Calc v
            setFV(getData().U, getData().V, getData().H, getData().F1, getData().F2);

            System.arraycopy(getData().U, 0, getData().C, 0, SZ_LAM * SZ_FI);
            splitLam(getData().V, getData().C, getData().Mu1, getData().F1, getData().V1);

            System.arraycopy(getData().V, 0, getData().C, 0, SZ_LAM * SZ_FI);
            splitFi(getData().V, getData().C, getData().Mu2, getData().F2, getData().V2);

            // Final
            splitSum(getData().U, getData().U1, getData().U2);
            splitSum(getData().V, getData().V1, getData().V2);
        }

        сalcNorma(2. * (i - 1) * dt / T, getData().U2, getData().V2, getData().U, getData().V);

    }

}
