package org.oracul;

import org.oracul.data.ForecastInfo;

/**
 * User: Pavlo_Ivanenko
 * Date: 9/23/13
 * Time: 1:51 PM
 */
public interface Oracul {
    void predict(ForecastInfo forecastInfo);
}
