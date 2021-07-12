package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

/**
 * Monte Carlo simulation with multiple temperatures (holds multiple objects of {@link SingleTMonteCarloSimulation})
 */
public class MultipleTMonteCarloSimulation extends MonteCarloSimulation implements Closeable, Runnable {
    // for serialization. should not be changed
    private static final long serialVersionUID = -7500236380863421871L;
    /** whether to turn off parallel tempering (never swap temperatures) */
    private final boolean parallelTemperingOff;
    /** Array of temperatures that are being simulated */
    protected final double[] T;
    /** Array of single temperature simulations */
    private SingleTMonteCarloSimulation[] simulations;
    /** Random number generator used for swapping temperatures in parallel tempering */
    private final MersenneTwister rnd;
    /** Current number of sweeps */
    private long sweeps;
    /** Checkpointer object used for saving and reading the simulation checkpoints */
    private transient SimulationCheckpointer checkpointer;
    /** Typical magnetic moment based on the external transverse field */
    public final double spinSize;
    /** tolerance for convergence of self-consistent calculation */
    public final double tol;
    /** Exchange interaction parameter */
    public final double J_ex;

    /**
     * Get one of the different temperature simulations
     * @param i the number of the simulation to return
     * @return A single temperature simulation object reference
     */
    public SingleTMonteCarloSimulation getIthSubSimulation(int i){ return this.simulations[i]; }

    /**
     * Initiates this simulation (by initiating all sub-simulations) and creates a checkpoint if directed
     */
    public void initSimulation(){
        for (int i=0;i<simulations.length;i++){
            simulations[i].initSimulation();
        }
        if (checkpoint) checkpointer.writeCheckpoint(this);
    }

    public void printSimulationState(){
        // Prints the current state of the lattice (spin configuration, magnetic moments and local fields)
        // for all sub-simulations
        for (int i=0;i<simulations.length;i++){
            simulations[i].printSimulationState();
        }
    }

    public void run() {
        run('s');
    }

    /**
     * try and swap (t)th and (t+1)th simulations. sets {@code lastSwapAccepted} within that sub-simulation
     * @param t simulation index to try to swap
     */
    public void trySwitch(int t){
        // try and swap adjacent temperatures.
        // see Hukushima et al. (1996) JPSJ: https://doi.org/10.1143/JPSJ.65.1604
        if (!parallelTemperingOff) {
            double thisEnergy = simulations[t].getCurrentEnergy();
            double nextEnergy = simulations[t+1].getCurrentEnergy();
            double delta = (1/T[t] - 1/T[t + 1])*(thisEnergy - nextEnergy);
            if (rnd.nextDouble() < Math.exp(delta)) {
                simulations[t].swap(simulations[t + 1]);
                simulations[t].setLastSwapAcceptance(true);
            }else{
                simulations[t].setLastSwapAcceptance(false);
            }
        }else{
            simulations[t].setLastSwapAcceptance(false);
        }
    }

    /**
     * Run the simulation
     * @param mode run mode: 's' for serial and 'p' for parallel
     */
    public void run(final char mode) {
        ExecutorService executor=null;
        if (mode=='p') {
            // create a thread pool
            executor = Executors.newFixedThreadPool(simulations.length);
        }
        while (sweeps<maxSweeps) {
            if (mode == 's') {
                //serial mode; run a MC sweeps for each of the temperatures one after the other
                for (int i = 0; i < simulations.length; i++) {
                    simulations[i].run();
                }
                // try and switch each pair of adjacent temperatures
                for (int i=0; i < simulations.length-1; i++){
                    trySwitch(i);
                }
            } else if (mode == 'p') {
                //parallel mode
                List<Callable<Object>> jobs = Arrays.stream(simulations).map(Executors::callable).collect( Collectors.toList() );
                // run one MC sweep for each of the temperatures in parallel
                try {
                    executor.invokeAll(jobs);
                }catch (InterruptedException e){
                    throw new RuntimeException("one of the monte carlo sweep threads encountered an error: " + e.getMessage());
                }
                // when all threads finish the MC sweep, try and switch each pair of adjacent temperatures
                for (int i = 0; i < simulations.length - 1; i++) {
                    trySwitch(i);
                }
            }
            sweeps++;
            if (sweeps%simulations[0].getOutWriter().getNumOfBufferedRows()==0 || sweeps==maxSweeps || sweeps==1) {    // every obsPrintSweepNum sweeps or at the last one or first one
                // write/rewrite the checkpoint
                if (checkpoint) {
                    checkpointer.writeCheckpoint(this);
                }
                // at the last sweep, also write "Done" with the time to the console
                System.out.println((sweeps==maxSweeps ? "Done: " : "") + LocalDateTime.now());
            }
            // print the progress if directed to
            if (simulations[0].getOutWriter().isPrintProgress()) System.out.println(String.format("%.2f",100.0*sweeps/maxSweeps) + "% complete                ");
        }
        if (mode=='p') {
            // close the thread pool
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

        // print out how many times each of the methods for solving the self-consistent calculation was used
        int[] methodsUsed = new int[20];
        for (int j=0; j < methodsUsed.length; j++) {
            for (int i = 0; i < simulations.length; i++) {
                methodsUsed[j] += simulations[i].getNumOfTimesMethodWasUsed(j);
            }
        }
        System.out.println("# Used methods statistics: " + Arrays.toString(methodsUsed));
    }

    /**
     * Close the simulation
     * @throws IOException
     */
    public void close() throws IOException {
        for (SingleTMonteCarloSimulation subSimulation : simulations){
            subSimulation.close();
        }
    }

    public void setCheckpointer(final SimulationCheckpointer checkpointer){
        this.checkpointer=checkpointer;
    }

    /**
     * Constructs a new {@code MultipleTMonteCarloSimulation} object with the given parameters
     * @param T array of temperatures to simulate
     * @param subSimulations array of single-temperature simulations to run under this multiple-temperature simulation
     * @param maxSweeps total number of MC sweeps
     * @param seed random number generator seed
     * @param rnd random number generator object
     * @param continueFromSave whether to continue from a previously saved simulation
     * @param parallelTemperingOff whether to turn off parallel tempering (never swap temperatures)
     * @param checkpoint whether to periodically create a checkpoint from which the simulation can be restarted
     * @param checkpointer Checkpointer object used for saving and reading the simulation checkpoints
     * @param spinSize Typical magnetic moment based on the external transverse field
     * @param tol tolerance for convergence of self-consistent calculation
     * @param J_ex Exchange interaction parameter
     */
    public MultipleTMonteCarloSimulation(final double[] T, final SingleTMonteCarloSimulation[] subSimulations, final long maxSweeps, final long seed, final MersenneTwister rnd,
                                         final boolean continueFromSave, final boolean parallelTemperingOff, final boolean checkpoint,
                                         final SimulationCheckpointer checkpointer, final double spinSize, final double tol, final double J_ex){
        this.parallelTemperingOff=parallelTemperingOff;
        this.maxSweeps=maxSweeps;
        this.seed=seed;
        this.rnd=rnd;
        this.continueFromSave=continueFromSave;
        this.T = Arrays.copyOf(T, T.length);
        this.simulations = subSimulations;
        this.checkpoint=checkpoint;
        this.sweeps=0;
        this.checkpointer=checkpointer;
        this.spinSize=spinSize;
        this.tol=tol;
        this.J_ex=J_ex;
    }

    public boolean isParallelTemperingOff() {
        return parallelTemperingOff;
    }

    /**
     * Returns the temperature schedule
     * @return a copy of the temperature schedule array
     */
    public double[] getT() {
        return Arrays.copyOf(T,T.length);
    }
}
