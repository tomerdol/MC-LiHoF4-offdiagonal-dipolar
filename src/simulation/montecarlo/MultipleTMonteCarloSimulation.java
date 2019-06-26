package simulation.montecarlo;

import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class MultipleTMonteCarloSimulation extends MonteCarloSimulation implements Closeable {
    // TODO check if these can be private now that I don't use serial & parallel classes
    protected final boolean parallelTempetingOff;
    protected int[] acceptanceRateCount;
    protected int[] acceptanceRateSum;
    protected final double[] T;
    protected SingleTMonteCarloSimulation[] simulations;
    protected final MersenneTwister rnd;
    protected long sweeps;
    protected transient SimulationCheckpointer checkpointer;

    public SingleTMonteCarloSimulation getIthSubSimulation(int i){ return this.simulations[i]; }

    public void initSimulation(){
        for (int i=0;i<simulations.length;i++){
            simulations[i].initSimulation();
        }
    }


    public void run() {
        run('s');
    }

    // try and swap (t)th and (t+1)th simulations
    public void trySwitch(int t){
        double thisEnergy = simulations[t].getCurrentEnergy();
        double nextEnergy = simulations[t+1].getCurrentEnergy();
        double delta = (1/T[t] - 1/T[t + 1])*(thisEnergy - nextEnergy);
        if (!parallelTempetingOff && rnd.nextDouble() < Math.exp(delta)) {
            simulations[t].swap(simulations[t+1]);
            simulations[t].getOutWriter().writeSwapAcceptance(true);
            acceptanceRateCount[t]++;
            acceptanceRateSum[t]++;
            simulations[t].incAcceptanceRateCount();
            simulations[t].incAcceptanceRateSum();
        }else{
            simulations[t].getOutWriter().writeSwapAcceptance(false);
            acceptanceRateCount[t]++;
            simulations[t].incAcceptanceRateCount();
        }


    }

    public void run(char mode) {
        while (sweeps<maxSweeps) {
            if (mode == 's') {
                //serial mode

                for (int i = 0; i < simulations.length; i++) {
                    simulations[i].run();
                }
                for (int i = 0; i < simulations.length - 1; i++) {
                    trySwitch(i);
                }

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
