package simulation.montecarlo;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class MultipleTMonteCarloSimulation extends MonteCarloSimulation implements Closeable {
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

    private void updateAcceptanceRates(int t, boolean swapAcceptance){
        simulations[t].getOutWriter().writeSwapAcceptance(swapAcceptance);
        acceptanceRateCount[t]++;
        simulations[t].incAcceptanceRateCount();
        if (swapAcceptance){
            acceptanceRateSum[t]++;
            simulations[t].incAcceptanceRateSum();
        }

    }

    // try and swap (t)th and (t+1)th simulations
    public void trySwitch(int t){
        double thisEnergy = simulations[t].getCurrentEnergy();
        double nextEnergy = simulations[t+1].getCurrentEnergy();
        double delta = (1/T[t] - 1/T[t + 1])*(thisEnergy - nextEnergy);
        if (!parallelTempetingOff && rnd.nextDouble() < Math.exp(delta)) {
            simulations[t].swap(simulations[t+1]);
            updateAcceptanceRates(t, true);
        }else{
            updateAcceptanceRates(t, false);
        }


    }

    public void run(char mode) {
        while (sweeps<maxSweeps) {
            if (mode == 's') {
                //serial mode

                for (int i = 0; i < simulations.length; i++) {
                    simulations[i].run();
                }
                int i;
                for (i = 0; i < simulations.length - 1; i++) {
                    trySwitch(i);
                }
                updateAcceptanceRates(i, true); // last temperature has no acceptance rate so it is always saved as true

            } else if (mode == 'p') {
                //parallel mode

            }
            sweeps++;
            if (checkpoint) {
                checkpointer.writeCheckpoint(this);
            }

            if (simulations[0].getOutWriter().isPrintProgress()) System.out.println(String.format("%.2f",100.0*sweeps/maxSweeps) + "% complete                ");
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
