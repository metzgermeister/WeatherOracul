package org.oracul.data;

import static org.oracul.OraculConstants.SZ_FI;
import static org.oracul.OraculConstants.SZ_LAM;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * User: Pavlo_Ivanenko Date: 9/20/13 Time: 2:47 PM
 */
public class OraculData {
    public static final int NODES_COUNT = SZ_FI * SZ_LAM;

    public final double U[] = new double[NODES_COUNT];
    public final double V[] = new double[NODES_COUNT];

    public final double H1[] = new double[NODES_COUNT];
    public final double H2[] = new double[NODES_COUNT];
    public final double H3[] = new double[NODES_COUNT];

    public final double fiBndCnd[][][][] = new double[2][3][2][SZ_LAM];
    public final double lamBndCnd[][][][] = new double[2][3][2][SZ_FI];

    public final double H[] = new double[NODES_COUNT];
    public final double U1[] = new double[NODES_COUNT];
    public final double V1[] = new double[NODES_COUNT];
    public final double U2[] = new double[NODES_COUNT];
    public final double V2[] = new double[NODES_COUNT];
    public final double F1[] = new double[NODES_COUNT];
    public final double F2[] = new double[NODES_COUNT];
    public final double Mu1[] = new double[NODES_COUNT];
    public final double Mu2[] = new double[NODES_COUNT];
    public final double C[] = new double[NODES_COUNT];

    public void importData(String inputDirectory) throws FileNotFoundException {
        int i, j, k;
        int m[] = new int[2];
        double[] pHUV[] = { U, V, H3 };

        String flNames[] = { "01Hs.inp", "01Us.inp", "01Vs.inp", "02Hs.inp", "02Us.inp", "02Vs.inp", "03Hs.inp",
                "03Us.inp", "03Vs.inp" };

        // Set initial condition for u, v
        for (i = 0; i <= 1; i++) {
            String flName = flNames[1 + i];
            System.out.println("reading file " + flName);

            File file = new File(inputDirectory + flName);
            Scanner scanner = new Scanner(file);

            for (j = 0; j <= (NODES_COUNT - 1); j++) {
                pHUV[i][j] = Double.parseDouble(scanner.next());
            }
            scanner.close();
        }

        // Set H
        pHUV[0] = H1;
        pHUV[1] = H2;
        for (i = 0; i <= 2; i++) {
            String flName = flNames[3 * i];
            System.out.println("reading file " + flName);

            File file = new File(inputDirectory + flName);
            Scanner scanner = new Scanner(file);

            for (j = 0; j <= (NODES_COUNT - 1); j++) {
                pHUV[i][j] = Double.parseDouble(scanner.next());
            }
            scanner.close();
        }

        // For boundary condition
        // [u, v] [time: 1, 2, 3] [boundary: 0, SZ] [knots]
        for (i = 0; i <= 1; i++)
            // u, v
            for (j = 0; j <= 2; j++) // time 1, 2, 3
            {
                String flName = flNames[3 * j + i + 1];
                System.out.println("reading file " + flName);

                File file = new File(inputDirectory + flName);
                Scanner scanner = new Scanner(file);

                for (k = 0; k <= (NODES_COUNT - 1); k++) {
                    m[0] = k / SZ_LAM;
                    m[1] = k % SZ_LAM;
                    double number = Double.parseDouble(scanner.next());
                    if ((m[0]) == 0) {
                        fiBndCnd[i][j][0][m[1]] = number;
                    } else if ((m[0]) == (SZ_FI - 1)) {
                        fiBndCnd[i][j][1][m[1]] = number;
                    }
                    if ((m[1]) == 0) {
                        lamBndCnd[i][j][0][m[0]] = number;
                    } else if ((m[1]) == (SZ_LAM - 1)) {
                        lamBndCnd[i][j][1][m[0]] = number;
                    }
                }
                scanner.close();
            }

    }

}
