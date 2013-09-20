package org.oracul.data;

/**
 * User: Pavlo_Ivanenko Date: 9/20/13 Time: 3:18 PM
 */
public class ForecastInfo {
    private final int timeStep;
    private final int timeLimit;

    public ForecastInfo(int timeStep, int timeLimit) {
        this.timeStep = timeStep;
        this.timeLimit = timeLimit;
    }

    public int getTimeStep() {
        return timeStep;
    }

    public int getTimeLimit() {
        return timeLimit;
    }

}
