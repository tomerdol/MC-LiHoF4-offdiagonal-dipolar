package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;


public class SerialMonteCarloSimulation extends MultipleTMonteCarloSimulation {

    public void run() {
        // TODO
    }

    public SerialMonteCarloSimulation(final double[] T, final SingleTMonteCarloSimulation[] subSimulations, final long maxSweeps, final long seed, final MersenneTwister rnd, final boolean continueFromSave, final boolean realTimeEqTest, final boolean parallelTempetingOff, final boolean checkpoint){
        super(T, subSimulations, maxSweeps, seed, rnd, continueFromSave, realTimeEqTest, parallelTempetingOff, checkpoint);
    }



}
