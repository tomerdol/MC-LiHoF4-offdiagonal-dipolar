package simulation.montecarlo;

import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.MersenneTwister;
import simulation.mmsolve.ConvergenceException;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.io.Serializable;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;


public class SingleTMonteCarloSimulation extends MonteCarloSimulation implements Serializable, Closeable, Runnable {
    private static final long serialVersionUID = -7085068052197341667L;
    private final double T;
    private final int temperatureIndex;  // should be -1 if not part of multiple T simulation
    private final int totalNumOfTemperatures;
    private long sweeps;
    private long currentBinCount;
    private double[] binAvg;
    private final ArrayList<CircularFifoQueue<Pair<Double,Double>>> equilibratingObs;
    private Lattice lattice;
    private final MersenneTwister rnd;
    private transient OutputWriter outWriter;
    private double currentEnergy;
    private int acceptanceRateCount;
    private int acceptanceRateSum;
    private boolean[] equilibratedObs;
    private boolean equilibrated;
    private boolean lastSwapAccepted;

    public final double spinSize, tol, J_ex;
    // parameters for the iterative solvers
    private transient int maxIter;
    private transient double alpha;
    private transient BufferedWriter outProblematicConfigs;

    public SingleTMonteCarloSimulation(final double T, final int temperatureIndex, final int totalNumOfTemperatures, final Lattice lattice, final int numOfObservables, final long maxSweeps,
                                       final long seed, final MersenneTwister rnd, final boolean continueFromSave, final boolean realTimeEqTest,
                                       final OutputWriter out, final boolean checkpoint, final int maxIter, final double alpha, final BufferedWriter outProblematicConfigs,
                                       final double spinSize, final double tol, final double J_ex) {
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
        this.totalNumOfTemperatures=totalNumOfTemperatures;
        this.outProblematicConfigs=outProblematicConfigs;
        this.acceptanceRateCount=0;
        this.acceptanceRateSum=0;
        this.equilibrated=false;
        this.equilibratedObs = new boolean[4];
        this.lastSwapAccepted=false;
        this.spinSize=spinSize;
        this.tol=tol;
        this.J_ex=J_ex;

        if (lattice.energyTable==null || lattice.momentTable==null || lattice.nnArray==null || lattice.exchangeIntTable==null || lattice.intTable==null) {
            throw new NullPointerException("The Lattice given to the SingleTMonteCarloSimulation has null some pointers. ");
        }
    }

    public void addObservableToBinAvg(int numOfObservables){
        if (binAvg.length < numOfObservables) {
            double[] newBinAvg = new double[numOfObservables]; // Add column (and its sd) to the end
            int i;
            for (i = 0; i < binAvg.length; i++) {
                newBinAvg[i] = binAvg[i];
            }
            for (; i < numOfObservables; i++) {
                newBinAvg[i] = Double.NaN;
            }
            binAvg = newBinAvg;
        }
        // else, do nothing
    }

    public void swap(SingleTMonteCarloSimulation other){
        // ****************swap energies****************
        double tempEnergy = this.currentEnergy;
        this.currentEnergy = other.currentEnergy;
        other.currentEnergy = tempEnergy;
        // *********************************************

        // ****************swap lattices****************
        Lattice tempLattice = this.lattice;
        this.lattice = other.lattice;
        other.lattice = tempLattice;
        // *********************************************

    }

    public void incAcceptanceRateSum() {
        this.acceptanceRateSum++;
    }

    public void incAcceptanceRateCount() {
        this.acceptanceRateCount++;
    }

    public boolean isEquilibrated() {
        return equilibrated;
    }

    public double getCurrentEnergy() {
        return currentEnergy;
    }

    public boolean wasLastSwapAccepted() { return lastSwapAccepted; }

    public void setLastSwapAcceptance(final boolean swapAccepted) {
        this.lastSwapAccepted=swapAccepted;
        if (swapAccepted){
            this.incAcceptanceRateSum();
        }
        this.incAcceptanceRateCount();
    }

    public void setOutProblematicConfigs(final BufferedWriter outProblematicConfigs) {
        this.outProblematicConfigs = outProblematicConfigs;
    }

    /**
     * Performs a Monte Carlo metropolis step
     * @param lattice - spin lattice
     * @param T - Temperature
     * @param rnd - PRNG
     * @return Array containing the (new) energy due to local transversal fields and the energy difference in the longitudinal energy compared to the previous configuration
     */
    public static deltaEnergyAndLattice metropolisStep(Lattice lattice, double T, MersenneTwister rnd, int maxIter, double alpha, final double tol) throws ConvergenceException
    {

        double initLongEnergy = lattice.getEnergy();

        // choose random spin
        int flippedSpin = rnd.nextInt(lattice.getN());

        lattice.flipSpin(maxIter, tol, flippedSpin, alpha, rnd);

        double deltaEnergy = lattice.getEnergy() - initLongEnergy;

        double prob;
        if (deltaEnergy <= 0) {
            prob = 1;
        }
        else {
            prob = Math.exp(-(deltaEnergy )/T);
        }
        // accept new state with calculated probability
        double r = rnd.nextDouble();

        if (r > prob) {   // change is denied
			/*
            flippedSpin.flipSpin();  // the change is reversed. no need to update other longitudinal fields

            updateLocalFieldsAfterFlip(lattice, flippedSpin, intTable[2]);
            updateLocalTransFieldsAfterFlip(lattice, flippedSpin, intTable);

            updateAllMagneticMoments(lattice, intTable, momentTable, maxIter, extBx);
			*/

			/*
			for (int i=0;i<lattice.length;i++) {
				lattice[i]=tempLattice[i];
			}
			*/

            return null;

        } else {    // change is accepted
            // updateLocalFieldsAfterFlip(lattice, flippedSpin, intTable[2]);    // update other *longitudinal* (that's the reason for the [2]) fields

        }

        return new deltaEnergyAndLattice(lattice, deltaEnergy);
        //return deltaEnergy;
    }

    public static void printProblematicConfig(Lattice lattice, int flippedSpin, BufferedWriter out){
        synchronized (out) {
            try {
                singleSpin[] arr = lattice.getArray();
                out.write("#" + flippedSpin + System.lineSeparator());
                for (int i = 0; i < arr.length; i++) {
                    out.write(arr[i].getSpin() + "," + arr[i].getSpinSize() + System.lineSeparator());
                }

            } catch (IOException e) {
                System.err.println("error writing problematic config to file" + e.toString());
            }
        }
    }


    public void monteCarloSweep() {
        // a monte carlo sweep consists of ~N steps
        for (int step = 0; step < lattice.getN(); step++) {

            if (outWriter.isPrintProgress()) {
                int perc = (int) 10.0 * step /  lattice.getN();
                char[] full = new char[perc];
                char[] empty = new char[10 - perc];
                Arrays.fill(full, '#');
                Arrays.fill(empty, ' ');
                String strFull = new String(full);
                String strEmpty = new String(empty);
                System.out.print((totalNumOfTemperatures>0 ? "[" + (temperatureIndex + 1) + "/" + totalNumOfTemperatures + "] " : "") +
                        "[" + strFull + strEmpty + "]" + step + "/" +  lattice.getN() + "     \r");
            }

            Lattice tempLattice = new Lattice(lattice);
            Double deltaEnergy;
            deltaEnergyAndLattice energyAndLattice;
            try {
                energyAndLattice = metropolisStep(lattice, T, rnd, maxIter, alpha, tol);

                if (energyAndLattice!=null){
                    deltaEnergy=energyAndLattice.getDeltaEnergy();
                    lattice=energyAndLattice.getLattice();
                }else{
                    deltaEnergy=null;
                }
//                successfulSteps++;

                //System.out.println("step " + step + " successful");
            } catch (ConvergenceException e){
                printProblematicConfig(tempLattice, e.getFlippedSpin(), outProblematicConfigs);
                deltaEnergy=null;
//                failedSteps++;
                System.err.println("There was an error converging, had to abort metropolis step ("+step+ ", sweep="+sweeps+" T="+T+") and try another one. \n" + e.getMessage());
            } catch (IndexOutOfBoundsException e){
                deltaEnergy=null;
//                failedSteps++;
                System.err.println("There was an error while calculating the derivative, had to abort metropolis step ("+step+ ", sweep="+sweeps+" T="+T+") and try another one.\n" + e.toString());
            }

            if (deltaEnergy == null) {
                // reverse change
                lattice=tempLattice;
            }else{
                tempLattice=null;
                currentEnergy += deltaEnergy.doubleValue();
            }

        }
    }

    public void printSimulationState(){

        if (outWriter.getOutType() == OutputType.SPIN) {
            singleSpin[] arr;
            if (lattice.isSuppressInternalTransFields()) {
                // copy lattice to a new object that does not suppress internal fields
                Lattice tempLatticeWithAllFields = new Lattice(lattice, false);
                tempLatticeWithAllFields.updateAllLocalFields();
                arr = tempLatticeWithAllFields.getArray();  // this is not very efficient since all the spins are copied twice (once on this line and once two lines back), but it is not part of the regular run so not that terrible
            } else {
                arr = lattice.getArray(); // deep copy
            }
            for (int i = 0; i < arr.length; i++) {
                outWriter.writeObservablesPerSpin(arr[i].getN(), arr[i].getSpin(), arr[i].getSpinSize(), arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz());
            }
            outWriter.flush();
        }
    }

    public void writeObservables(){
        double m = lattice.getMagnetization();
        double[] temp = lattice.getMagneticFields();
        double[] tempSpinSizes = lattice.getSpinSizes();
        double[] mk2 = lattice.getMK2();    // m(k)^2, used later for correlation length calculation
        double transFieldMaxConfig = lattice.getTransverseFieldMaximizingNNConfigsFrac();


        outWriter.writeObservablesVerbose(sweeps, m ,currentEnergy ,temp[0] ,temp[1] ,temp[2] ,temp[3] ,temp[4] ,temp[5] ,temp[6] ,temp[7] ,temp[8], tempSpinSizes[0] ,tempSpinSizes[1] ,mk2[0], mk2[1], mk2[2], transFieldMaxConfig, lastSwapAccepted);

        if (sweeps>0) currentBinCount++;
        //System.out.println(sweeps+ " , " + sweeps/2 + " , " + currentBinCount + " , " + (sweeps&1));

        if (currentBinCount==sweeps/2 + 1 && (sweeps&1)==1){	// last step of current bin
            // restart sums and counters for next bin
            double acceptanceRateForBin;	// acceptance rate for this bin & temperature

            if (acceptanceRateCount==0) {
                acceptanceRateForBin = -1;
            }else {
                acceptanceRateForBin = ((double) acceptanceRateSum) / acceptanceRateCount;
            }

            // print data from previous bin and
            outWriter.writeObservablesBin(currentBinCount, binAvg, acceptanceRateForBin);

            // add bin averages to queue. only 3 bins are kept at any time.
            equilibratingObs.get(0).add(Pair.of(binAvg[4]/currentBinCount,getStandardError(binAvg[4],binAvg[5],currentBinCount)));		// m^2
            equilibratingObs.get(1).add(Pair.of(binAvg[6]/currentBinCount,getStandardError(binAvg[6],binAvg[7],currentBinCount)));		// E
            equilibratingObs.get(2).add(Pair.of(binAvg[28]/currentBinCount,getStandardError(binAvg[28],binAvg[29],currentBinCount)));	// mk2
            equilibratingObs.get(3).add(Pair.of(binAvg[0]/currentBinCount,getStandardError(binAvg[0],binAvg[1],currentBinCount)));		// |m|

            for (int i=0;i<4;i++) {
                equilibratedObs[i] = checkEquilibration(equilibratingObs.get(i));
            }
            equilibrated = isAllTrue(equilibratedObs);

            //System.out.println(sweeps + ") printed average of bin size "+currentBinCount + " starting from " + binStart);
            for (int avgIndex=0;avgIndex<binAvg.length;avgIndex++) binAvg[avgIndex]=0;

            acceptanceRateCount=0;
            acceptanceRateSum=0;

            //binStart=sweeps;

            currentBinCount=0;

            if (equilibrated){
                // declare that equilibration has been reached
                // TODO test the equilibration against the python script. It seems to happen a little faster than expected
                System.out.println("System equilibrated! T="+T+", sweeps="+sweeps+". Date&Time: "+ LocalDateTime.now());

            }
        }

        addToAvg(binAvg,new double[]{Math.abs(m), m , m*m, currentEnergy , temp[0] , temp[1] , temp[2] , temp[3] , temp[4] , temp[5] , temp[6] , temp[7] , tempSpinSizes[0] , tempSpinSizes[1] , mk2[0], mk2[1], mk2[2], transFieldMaxConfig});


    }

    public void close() throws IOException{
        this.outWriter.close();
    }

    public Lattice getLattice(){ return this.lattice; }

    public OutputWriter getOutWriter() {
        return outWriter;
    }

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
    public static boolean checkEquilibration(ArrayList<CircularFifoQueue<Pair<Double,Double>>> arr){
        if (arr.size()==0)
            return false;

        boolean equilibrated=true;
        for (CircularFifoQueue<Pair<Double, Double>> binsOfObservable : arr) {
            equilibrated &= checkEquilibration(binsOfObservable);
        }

        return equilibrated;
    }

    // check whether the given queue of a single observable's averages and standard deviations agree within error bars
    public static boolean checkEquilibration(CircularFifoQueue<Pair<Double, Double>> binsOfObservable){
        boolean equilibrated=true;
        if (!binsOfObservable.isAtFullCapacity()){
            // queue has less than 3 bins
            equilibrated=false;
        } else {
            double upperBound = binsOfObservable.get(0).getLeft() + binsOfObservable.get(0).getRight(), lowerBound = binsOfObservable.get(0).getLeft() - binsOfObservable.get(0).getRight();
            for (Pair<Double, Double> bin : binsOfObservable) {
                upperBound = Math.min(upperBound, bin.getLeft() + bin.getRight());
                lowerBound = Math.max(lowerBound, bin.getLeft() - bin.getRight());
            }
            equilibrated = upperBound > lowerBound;
            if (equilibrated && binsOfObservable.maxSize()>2) {
                boolean monotonicallyInc = binsOfObservable.get(0).getLeft() < binsOfObservable.get(1).getLeft();
                boolean monotonicallyDec = binsOfObservable.get(0).getLeft() > binsOfObservable.get(1).getLeft();
                for (int i=2;i<binsOfObservable.size();i++) {
                    monotonicallyInc &= binsOfObservable.get(i-1).getLeft() < binsOfObservable.get(i).getLeft();
                    monotonicallyDec &= binsOfObservable.get(i-1).getLeft() > binsOfObservable.get(i).getLeft();
                }

                // condition is that the bins not be monotonic
                equilibrated = !(monotonicallyDec || monotonicallyInc);
            }
        }
        return equilibrated;
    }

    public static boolean isAllTrue(boolean... array)
    {
        for(boolean b : array) if(!b) return false;
        return true;
    }

    public void run() {
        monteCarloSweep();
        writeObservables();
        sweeps++;

        if (sweeps%outWriter.getNumOfBufferedRows()==0 || sweeps==maxSweeps) {	// every obsPrintSweepNum sweeps or at the last one
            if (outWriter.isPrintOutputToConsole()) System.out.println("T="+T);
            outWriter.flush();
        }
    }

    public String runParameters(String version, double[] T, String extraMessage, long mutualSeed, String tempScheduleFileName, boolean parallelTemperingOff){
        return  "# VERSION: " + version + System.lineSeparator() +
                "#" + LocalDateTime.now() + System.lineSeparator() +
                "#temperature_schedule: " + Arrays.toString(T) + System.lineSeparator() +
                "#T=" + this.T + ":" + temperatureIndex + System.lineSeparator() +
                "#" + Constants.constantsToString() + System.lineSeparator() +
                String.format("# Lx=%s, Ly=%s, Lz=%s, J_ex=%f, spinSize=%.8f, tol=%4.1e, extBx=%s, extBy=%s, maxSweeps=%s, suppressInternalTransFields=%s, " +
                        "continueFromSave=%s, maxIter=%s, bufferSize=%s, tempScheduleFileName=%s, parallelTemperingOff=%s, " +
                        "checkpoint=%s, folderName=%s, alpha=%s, output=%s ",lattice.getLx(),lattice.getLx(),lattice.getLz(),J_ex,spinSize,tol,lattice.getExtBx(), lattice.getExtBy(), maxSweeps,lattice.isSuppressInternalTransFields(),
                    continueFromSave, maxIter, outWriter.getBufferSize(),
                    tempScheduleFileName, parallelTemperingOff, checkpoint, outWriter.getFolderName(),
                    alpha, outWriter.getOutType()) + System.lineSeparator() +
                "#" + Constants.locationsToString() + "," + Constants.latticeVectorsToString() + System.lineSeparator() +
                "#seed=" + mutualSeed + " (" + seed + ")" + System.lineSeparator() +
                extraMessage;
    }

    public void printRunParameters(String version, double[] T, String extraMessage, long mutualSeed, String tempScheduleFileName, boolean parallelTemperingOff) throws IOException{
        // print some information to the beginning of the results file (or console):
        outWriter.print(runParameters(version, T, extraMessage, mutualSeed, tempScheduleFileName, parallelTemperingOff), true);
        outWriter.flush();
    }

    public void initSimulation(){
        boolean convergedConfig=false;
        int index=0;
        while (!convergedConfig) {
            lattice.randomizeConfig(rnd); // randomizes the spins and sets initial spin sizes as spinSize in the corresponding direction

            lattice.updateAllLocalFields();
            lattice.updateAllMagneticMoments(2*maxIter, tol, alpha);

            if (lattice.magneticMomentConvergence() > tol){
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
