package simulation.montecarlo;

import org.apache.commons.collections4.iterators.ArrayIterator;
import org.apache.commons.math3.random.MersenneTwister;
import simulation.mmsolve.MagneticMomentsSolveIter;

import java.io.Closeable;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.Iterator;

public class MultipleTMonteCarloSimulation extends MonteCarloSimulation implements Closeable {
    // TODO check if these can be private now that I don't use serial & parallel classes
    protected final boolean parallelTempetingOff;
    protected int[] acceptanceRateCount;
    protected int[] acceptanceRateSum;
    protected final double[] T;
    protected SingleTMonteCarloSimulation[] simulations;
    protected final MersenneTwister rnd;

    public SingleTMonteCarloSimulation getIthSubSimulation(int i){ return this.simulations[i]; }

    public void initSimulation(){
        for (int i=0;i<simulations.length;i++){
            simulations[i].initSimulation();
        }
    }


    public void run(){
        run('s');
    }

    public void run(char mode){
        if (mode=='s'){
            //serial mode

        }else if (mode=='p'){
            //parallel mode

        }
    }

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
