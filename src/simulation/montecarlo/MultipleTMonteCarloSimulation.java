package simulation.montecarlo;

import javafx.concurrent.Task;
import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.FileOutputStream;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class MultipleTMonteCarloSimulation extends MonteCarloSimulation implements Closeable, Runnable {
    // TODO check if these can be private now that I don't use serial & parallel classes
    private final boolean parallelTempetingOff;
    private int[] acceptanceRateCount;
    private int[] acceptanceRateSum;
    protected final double[] T;
    private SingleTMonteCarloSimulation[] simulations;
    private final MersenneTwister rnd;
    private long sweeps;
    private transient SimulationCheckpointer checkpointer;

    public SingleTMonteCarloSimulation getIthSubSimulation(int i){ return this.simulations[i]; }

    public void initSimulation(){
        for (int i=0;i<simulations.length;i++){
            simulations[i].initSimulation();
        }
    }


    public void run() {
        run('s');
    }


    // try and swap (t)th and (t+1)th simulations.
    // return if the boolean indicating whether the swap was performed
    public void trySwitch(int t){
        double thisEnergy = simulations[t].getCurrentEnergy();
        double nextEnergy = simulations[t+1].getCurrentEnergy();
        double delta = (1/T[t] - 1/T[t + 1])*(thisEnergy - nextEnergy);
        if (!parallelTempetingOff && rnd.nextDouble() < Math.exp(delta)) {
            simulations[t].swap(simulations[t+1]);
            simulations[t].setLastSwapAcceptance(true);
        }else{
            simulations[t].setLastSwapAcceptance(false);
        }
    }

    public void run(final char mode) {
        ExecutorService executor=null;
        if (mode=='p') {
            executor = Executors.newFixedThreadPool(simulations.length);
        }
        while (sweeps<maxSweeps) {
            if (mode == 's') {
                //serial mode

                for (int i = 0; i < simulations.length; i++) {
                    simulations[i].run();
                }
                for (int i=0; i < simulations.length-1; i++){
                    trySwitch(i);
                }
            } else if (mode == 'p') {
                //parallel mode
                List<Callable<Object>> jobs = Arrays.stream(simulations).map(Executors::callable).collect( Collectors.toList() );
                try {
                    executor.invokeAll(jobs);
                }catch (InterruptedException e){
                    throw new RuntimeException("one of the monte carlo sweep threads encountered an error: " + e.getMessage());
                }
                for (int i = 0; i < simulations.length - 1; i++) {
                    trySwitch(i);
                }
            }
            sweeps++;
            if (sweeps%simulations[0].getOutWriter().getObsPrintSweepNum()==0 || sweeps==maxSweeps) {    // every obsPrintSweepNum sweeps or at the last one
                if (checkpoint) {
                    checkpointer.writeCheckpoint(this);
                }
                System.out.println(LocalDateTime.now());
            }

            if (simulations[0].getOutWriter().isPrintProgress()) System.out.println(String.format("%.2f",100.0*sweeps/maxSweeps) + "% complete                ");
        }
        if (mode=='p') {
            executor.shutdown();
            try {
                if (!executor.awaitTermination(10, TimeUnit.MINUTES)) {
                    System.err.println("threads did not finish for 10 minutes. shutting down thread pool. ");
                    executor.shutdownNow();
                }
            } catch (InterruptedException e) {
                System.err.println("error closing threads. shutting down thread pool. ");
                executor.shutdownNow();
            }
        }
    }



    public void close() throws IOException {
        for (SingleTMonteCarloSimulation subSimulation : simulations){
            subSimulation.close();
        }
    }

    public void setCheckpointer(final SimulationCheckpointer checkpointer){
        this.checkpointer=checkpointer;
    }

    public MultipleTMonteCarloSimulation(final double[] T, final SingleTMonteCarloSimulation[] subSimulations, final long maxSweeps, final long seed, final MersenneTwister rnd, final boolean continueFromSave, final boolean realTimeEqTest, final boolean parallelTempetingOff, final boolean checkpoint, final SimulationCheckpointer checkpointer){
        this.parallelTempetingOff=parallelTempetingOff;
        this.maxSweeps=maxSweeps;
        this.seed=seed;
        this.rnd=rnd;
        this.continueFromSave=continueFromSave;
        this.realTimeEqTest=realTimeEqTest;
        this.T = Arrays.copyOf(T, T.length);
        this.simulations = subSimulations;
        this.checkpoint=checkpoint;
        this.sweeps=0;
        this.checkpointer=checkpointer;
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
