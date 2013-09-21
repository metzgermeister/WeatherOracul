package org.oracul;

import java.io.*;
import java.util.Date;

import org.oracul.data.ForecastInfo;
import org.oracul.data.OraculData;

import static org.oracul.OraculConstants.SZ_FI;
import static org.oracul.OraculConstants.SZ_LAM;

public class Launcher {

    public static final int TIME_STEP = 1;
    public static final int TIME_LIMIT = 3600;

    public static void main(String[] args) throws IOException {
        new Launcher().doStuff();
    }

    private void doStuff() throws IOException {

        long tv1 = new Date().getTime();
        OraculData oraculData = loadDataFromFiles();
        Oracul oracul = new Oracul(oraculData);
        oracul.predict(new ForecastInfo(TIME_STEP, TIME_LIMIT));


        long tv2 = new Date().getTime();
        System.out.println("time took:" + (tv2 - tv1));

        savePrediction(oraculData);

        System.out.println("The end");

    }

    private void savePrediction(OraculData oraculData) throws IOException {
        saveToFile("OUT\\U.out", oraculData.U);
        saveToFile("OUT\\V.out", oraculData.V);
    }

    public OraculData loadDataFromFiles() throws FileNotFoundException {
        OraculData data = new OraculData();
        String inputDirectory = System.getProperty("user.dir") + "\\resources\\IN\\";
        data.importData(inputDirectory);
        return data;
    }


    private void saveToFile(String pathToFile, double[] data) throws IOException {
        System.out.println("writing data to file " + pathToFile);
        FileWriter writer = new FileWriter(pathToFile);

        for (double d : data) {
            writer.write(String.valueOf(d) + "\n");
        }
        writer.close();
    }

}
