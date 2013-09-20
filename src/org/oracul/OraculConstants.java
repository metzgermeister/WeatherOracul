package org.oracul;

public class OraculConstants {
    public static final int SZ_LAM = 61;
    public static final int SZ_FI = 60;

    public static final double T = 24. * 3600.; // max time interval (in seconds) -- do not change
    public static final double T1 = 1. * 3600.; // calculate time interval (in seconds)
    public static final double dt = 100.; // seconds
    public static final double dh = 1.5 * Math.PI / 180.; // radians
    public static final double S = 3.0; // parameter of model dissipation
    public static final double ERH_RDS = 6.373e06; // Earth's radius in meters
    public static final double ERH_RDS_2 = 4.062e13; // Square of Earth's radius in meters
    public static final double RO_RATE = 1.07e-4;
    public static final double ERH_OMEGA = 7.292e-05; // Earth's rotation in sec^-1
    public static final double ERH_G = 9.80616; // acceleration of gravity in m/sec^2
    public static final double SRF_ALT = 11800.; // isobaric surface altitude

}
