package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.IOException;
import java.time.LocalDateTime;
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


    public MultipleTMonteCarloSimulation(final double[] T, final SingleTMonteCarloSimulation[] subSimulations, final long maxSweeps, final long seed, final MersenneTwister rnd, final boolean continueFromSave, final boolean realTimeEqTest, final boolean parallelTempetingOff, final boolean checkpoint){
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
        this.checkpoint=checkpoint;
    }

    public void printRunParameters(){

    }

    public boolean isParallelTempetingOff() {
        return parallelTempetingOff;
    }

    public double[] getT() {
        return Arrays.copyOf(T,T.length);
    }
}
