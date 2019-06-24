package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.IOException;
import java.util.Arrays;
import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

public abstract class MultipleTMonteCarloSimulation extends MonteCarloSimulation implements Closeable {
    protected final boolean parallelTempetingOff;
    protected int[] acceptanceRateCount;
    protected int[] acceptanceRateSum;
    protected final double[] T;
    protected SingleTMonteCarloSimulation[] simulations;
    protected final MersenneTwister rnd;

    public void close() throws IOException {
        for (SingleTMonteCarloSimulation subSimulation : simulations){
            subSimulation.close();
        }
    }


    public MultipleTMonteCarloSimulation(double[] T, SingleTMonteCarloSimulation[] subSimulations, long maxSweeps, long seed, MersenneTwister rnd, boolean continueFromSave, boolean realTimeEqTest, boolean parallelTempetingOff){
        this.parallelTempetingOff=parallelTempetingOff;
        this.acceptanceRateCount=new int[T.length];
        this.acceptanceRateSum=new int[T.length];
        this.maxSweeps=maxSweeps;
        this.seed=seed;
        this.rnd=rnd;
        this.continueFromSave=continueFromSave;
        this.realTimeEqTest=realTimeEqTest;
        this.T = Arrays.copyOf(T, T.length);
        this.simulations = subSimulations;
    }




}
