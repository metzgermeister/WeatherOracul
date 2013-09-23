package org.oracul.util;

import static java.lang.Math.exp;
import static org.oracul.util.OraculConstants.RO_RATE;

/**
 * User: Pavlo_Ivanenko Date: 9/20/13 Time: 6:19 PM
 */
public final class Utils {
    private Utils() {
    }

    public static double quadrInterpol(double r, double f0, double f1, double f2) {
        // f(0)=f0, f(1)=f1, f(2)=f2
        return f0 - r * (0.5 * f2 - 2. * f1 + 1.5 * f0 - r * (0.5 * f2 - f1 + 0.5 * f0));
    }

    public static double fRo(double z) {
        return 1.22 * exp(-RO_RATE * z);
    }


    public static double fG(double z, double fi) {
        double cs = Math.cos(2. * fi);
        return OraculConstants.ERH_G - 0.025928 * cs + (6.9e-5) * cs * cs - (3.086e-6) * z;
    }
}
