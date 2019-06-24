package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;


public class ParallelMonteCarloSimulation extends MultipleTMonteCarloSimulation {

    public void run() {
        // TODO
    }

    public ParallelMonteCarloSimulation(double[] T, SingleTMonteCarloSimulation[] subSimulations, long maxSweeps, long seed, MersenneTwister rnd, boolean continueFromSave, boolean realTimeEqTest, boolean parallelTempetingOff){
        super(T, subSimulations, maxSweeps, seed, rnd, continueFromSave, realTimeEqTest, parallelTempetingOff);
    }

}
