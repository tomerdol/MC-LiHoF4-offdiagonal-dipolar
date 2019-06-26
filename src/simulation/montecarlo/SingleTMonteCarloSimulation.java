package simulation.montecarlo;

import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.MersenneTwister;
import simulation.mmsolve.MagneticMomentsSolveIter;

import java.io.*;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;

public class SingleTMonteCarloSimulation extends MonteCarloSimulation implements Serializable, Closeable {
    private final double T;
    private final int temperatureIndex;  // should be -1 if not part of multiple T simulation
    private long sweeps;
    private long currentBinCount;
    private final double[] binAvg;
    private final ArrayList<CircularFifoQueue<Pair<Double,Double>>> equilibratingObs;
    private Lattice lattice;
    private final MersenneTwister rnd;
    private transient OutputWriter outWriter;
    private double currentEnergy;
    // parameters for the iterative solvers
    private transient int maxIter;
    private transient double alpha;

    public SingleTMonteCarloSimulation(final double T, final int temperatureIndex, final Lattice lattice, final int numOfObservables, final long maxSweeps,
                                       final long seed, final MersenneTwister rnd, final boolean continueFromSave, final boolean realTimeEqTest,
                                       final OutputWriter out, final boolean checkpoint, final int maxIter, final double alpha) {
        this.T = T;
        this.temperatureIndex=temperatureIndex;
        this.sweeps = 0;
        this.currentBinCount = 0;
        binAvg = new double[numOfObservables];
        equilibratingObs = generateEquilibrationQueues(3,4);
        this.lattice=lattice;
        this.maxSweeps=maxSweeps;
        this.seed=seed;
        this.continueFromSave=continueFromSave;
        this.realTimeEqTest=realTimeEqTest;
        this.rnd=rnd;
        this.outWriter=out;
        this.checkpoint=checkpoint;
        this.alpha=alpha;
        this.maxIter=maxIter;

        if (lattice.energyTable==null || lattice.momentTable==null || lattice.nnArray==null || lattice.exchangeIntTable==null || lattice.intTable==null) {
            throw new NullPointerException("The Lattice given to the SingleTMonteCarloSimulation has null some pointers. ");
        }
    }

    public void close() throws IOException{
        outWriter.close();
    }

    public Lattice getLattice(){ return this.lattice; }

    public void setOutWriter(final OutputWriter outWriter){
        this.outWriter=outWriter;
    }
    public void setMaxIter(final int maxIter){
        this.maxIter=maxIter;
    }
    public void setAlpha(final double alpha){
        this.alpha=alpha;
    }

    public static void addToAvg(double[] avgArr, double[] currentValueArray){
        if (avgArr.length>>1 != currentValueArray.length){
            throw new RuntimeException("error adding observable values to bin average. array lenth mismatch. ");
        }else{
            for (int i=0;i<currentValueArray.length;i++){
                avgArr[2*i]+=currentValueArray[i];
                avgArr[2*i+1]+=currentValueArray[i]*currentValueArray[i];
            }
        }
    }

    /**
     * Calculate the standard error given the sum, sum^2 and the number of samples, N.
     * @param sum sum of samples
     * @param sum2 sum of squared samples
     * @param N number of samples
     * @return the standard error (error bar) of the given samples
     */
    public static double getStandardError(double sum, double sum2, long N){
        return Math.sqrt((sum2/N - (sum/N)*(sum/N))/N);
    }

    // check whether the given array of averages and standard deviations agree within error bars
    // TODO: change this so that it only checks whether measurements have overlapping sd errors
    public static boolean checkEquilibration(ArrayList<CircularFifoQueue<Pair<Double,Double>>> arr){
        if (arr.size()==0)
            return false;

        boolean equilibrated=true;
        for (CircularFifoQueue<Pair<Double, Double>> binsOfObservable : arr) {
            if (!binsOfObservable.isAtFullCapacity() || !equilibrated){
                // queue has less than 3 bins
                equilibrated=false;
            } else {
                double upperBound = binsOfObservable.get(0).getLeft() + binsOfObservable.get(0).getRight(), lowerBound = binsOfObservable.get(0).getLeft() - binsOfObservable.get(0).getRight();
                for (Pair<Double, Double> bin : binsOfObservable) {
                    upperBound = Math.min(upperBound, bin.getLeft() + bin.getRight());
                    lowerBound = Math.max(lowerBound, bin.getLeft() - bin.getRight());
                }
                equilibrated = upperBound > lowerBound;
                for (Pair<Double, Double> bin : binsOfObservable) {
                    if (bin.getLeft() > upperBound || bin.getLeft() < lowerBound) {
                        equilibrated = false;
                    }
                }
            }
        }

        return equilibrated;
    }

    public static boolean isAllTrue(boolean... array)
    {
        for(boolean b : array) if(!b) return false;
        return true;
    }

    public void run(){

    }

    public void printRunParameters(double[] T, String extraMessage) throws IOException{
        // print some information to the begining of the file:
        outWriter.print("#" + LocalDateTime.now(), true);
        outWriter.print("#temperature_schedule: "+ Arrays.toString(T), true);
        outWriter.print("#T="+T + ":" + temperatureIndex, true);
        //print constants:
        outWriter.print("#" + Constants.constantsToString(), true);
        outWriter.print(String.format("# Lx=%s, Ly=%s, Lz=%s, extBx=%s, maxSweeps=%s, suppressInternalTransFields=%s, " +
                        "continueFromSave=%s, maxIter=%s, bufferSize=%s, tempScheduleFileName=%s, parallelTemperingOff=%s, " +
                        "checkpoint=%s, folderName=%s, alpha=%s, verboseOutput=%s ",lattice.getLx(),lattice.getLx(),lattice.getLz(),lattice.getExtBx(), maxSweeps,lattice.isSuppressInternalTransFields(),
                continueFromSave, maxIter, outWriter.getBufferSize(),
                checkpoint, outWriter.getFolderName(),
                alpha, outWriter.isVerboseOutput()),true);
        outWriter.print("#" + Constants.locationsToString(), true);
        outWriter.print("#seed=" + seed, true);
        outWriter.print(extraMessage, true);

    }

    public void initSimulation(){
        boolean convergedConfig=false;
        int index=0;
        while (!convergedConfig) {
            lattice.randomizeConfig(rnd); // randomizes the spins and sets initial spin sizes as spinSize in the corresponding direction

            lattice.updateAllLocalFields();
            lattice.updateAllMagneticMoments(maxIter, Constants.tol, alpha);

            if (lattice.magneticMomentConvergence() > Constants.tol){
                convergedConfig=false;
            }else{
                convergedConfig=true;
            }

            if (!convergedConfig) System.out.println("initial setup convergence index (T [num]="+T+" ["+temperatureIndex+']'+"): " + ++index);
        }

        currentEnergy = lattice.getEnergy();
    }

    public static ArrayList<CircularFifoQueue<Pair<Double,Double>>> generateEquilibrationQueues(int numOfBins, int numOfObservables){
        ArrayList<CircularFifoQueue<Pair<Double,Double>>> obsArr = new ArrayList<>(numOfObservables);
        for (int observable=0;observable<numOfObservables;observable++){
            obsArr.add(new CircularFifoQueue<>(numOfBins));
        }
        return obsArr;
    }

    public double getT() {
        return T;
    }

    public long getSweeps() {
        return sweeps;
    }

    public long getCurrentBinCount() {
        return currentBinCount;
    }
}
