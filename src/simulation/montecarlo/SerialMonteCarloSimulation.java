package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;


public class SerialMonteCarloSimulation extends MultipleTMonteCarloSimulation {

    public void run() {
        // TODO
    }

    public SerialMonteCarloSimulation(double[] T, SingleTMonteCarloSimulation[] subSimulations, long maxSweeps, long seed, MersenneTwister rnd, boolean continueFromSave, boolean realTimeEqTest, boolean parallelTempetingOff){
        super(T, subSimulations, maxSweeps, seed, rnd, continueFromSave, realTimeEqTest, parallelTempetingOff);
    }



}
