package org.oracul;

import java.io.FileNotFoundException;
import java.util.Date;

import org.oracul.data.ForecastInfo;
import org.oracul.data.OraculData;

public class Launcher {

    public static final int TIME_STEP = 1;
    public static final int TIME_LIMIT = 3600;

    public static void main(String[] args) throws FileNotFoundException {
        new Launcher().doStuff();
    }

    private void doStuff() throws FileNotFoundException {

        long tv1 = new Date().getTime();
        OraculData oraculData = loadDataFromFiles();
        Oracul oracul = new Oracul(oraculData);
        oracul.predict(new ForecastInfo(TIME_STEP, TIME_LIMIT));

        long tv2 = new Date().getTime();
        System.out.println("time took:" + (tv2 - tv1));

        System.out.println("The end");

    }

    public OraculData loadDataFromFiles() throws FileNotFoundException {
        OraculData data = new OraculData();
        String inputDirectory = System.getProperty("user.dir") + "\\resources\\IN\\";
        data.importData(inputDirectory);
        return data;
    }

}
