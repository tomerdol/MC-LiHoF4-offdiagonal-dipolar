package simulation.montecarlo;

import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.random.MersenneTwister;

import java.io.Closeable;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

public class SingleTMonteCarloSimulation extends MonteCarloSimulation implements Serializable, Closeable {
    private final double T;
    private long sweeps;
    private long currentBinCount;
    private final double[] binAvg;
    private final ArrayList<CircularFifoQueue<Pair<Double,Double>>> equilibratingObs;
    private Lattice lattice;
    private final MersenneTwister rnd;
    private final OutputWriter outWriter;

    public SingleTMonteCarloSimulation(double T, Lattice lattice, int numOfObservables, long maxSweeps, long seed, MersenneTwister rnd, boolean continueFromSave, boolean realTimeEqTest, OutputWriter out) {
        this.T = T;
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

        if (lattice.energyTable==null || lattice.momentTable==null || lattice.nnArray==null || lattice.exchangeIntTable==null || lattice.intTable==null) {
            throw new NullPointerException("The Lattice given to the SingleTMonteCarloSimulation has null some pointers. ");
        }
    }

    public void close() throws IOException{
        outWriter.close();
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

    public static ArrayList<CircularFifoQueue<Pair<Double,Double>>> generateEquilibrationQueues(int numOfBins, int numOfObservables){
        ArrayList<CircularFifoQueue<Pair<Double,Double>>> obsArr = new ArrayList<>(numOfObservables);
        for (int observable=0;observable<numOfObservables;observable++){
            obsArr.add(new CircularFifoQueue<>(numOfBins));
        }
        return obsArr;
    }
}
