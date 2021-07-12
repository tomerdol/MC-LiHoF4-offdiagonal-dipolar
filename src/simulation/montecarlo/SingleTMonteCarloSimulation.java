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
    // for serialization. should not be changed
    private static final long serialVersionUID = -7085068052197341667L;
    /** temperature */
    private final double T;
    /** index of this simulation within a {@link MultipleTMonteCarloSimulation} */
    private final int temperatureIndex;  // should be -1 if not part of multiple T simulation
    /** total number of simulated temperatures within the containing {@link MultipleTMonteCarloSimulation} */
    private final int totalNumOfTemperatures;
    /** current number of sweeps (index) */
    private long sweeps;
    /** current number of samples in current bin */
    private long currentBinCount;
    /** accumulated samples of all observables in the current bin (at even indices) and their squares (at odd indices) */
    private double[] binAvg;
    /** the {@link Lattice} object being simulated */
    private Lattice lattice;
    /** random number generator */
    private final MersenneTwister rnd;
    /** writer object used to output the results */
    private transient OutputWriter outWriter;
    /** current energy of the system */
    private double currentEnergy;
    /** total number of attempted parallel tempering swaps within the current bin */
    private int acceptanceRateCount;
    /** number of accepted parallel tempering swaps within the current bin */
    private int acceptanceRateSum;
    /** whether the last PT swap was accepted */
    private boolean lastSwapAccepted;
    /** array that counts how many times each method for solving the self-consistent calculation was successfully used */
    private int[] methodsUsed;
    /** the typical magnetic moment of the spins */
    public final double spinSize;
    /** tolerance for convergence of the self-consistent calculation */
    public final double tol;
    /** exchange interaction parameter */
    public final double J_ex;

    // parameters for the iterative solvers
    /** Maximum iterations for iterative solver (Gauss-Seidel) */
    private transient int maxIter;
    /** Relaxation parameter for iterative solver */
    private transient double alpha;
    /** Writer to output configurations that were not successfully solved self-consistently by any of the available methods */
    private transient BufferedWriter outProblematicConfigs;

    /**
     * Constructs a new single-temperature Monte Carlo simulation
     * @param T temperature
     * @param temperatureIndex index of this simulation within a {@link MultipleTMonteCarloSimulation}
     * @param totalNumOfTemperatures total number of simulated temperatures within the containing {@link MultipleTMonteCarloSimulation}
     * @param lattice the {@link Lattice} object to be simulated
     * @param numOfObservables total number of observables to be measured and saved
     * @param maxSweeps total number of MC sweeps to perform
     * @param seed random number generator seed
     * @param rnd random number generator object
     * @param continueFromSave whether to continue the simulation from a saved checkpoint if exists
     * @param out writer object used to output the results
     * @param checkpoint whether to periodically create a checkpoint from which the simulation can be restarted
     * @param maxIter Maximum iterations for iterative solvers
     * @param alpha Relaxation parameter for iterative solver
     * @param outProblematicConfigs Writer to output configurations that were not successfully solved self-consistently by any of the available methods
     * @param spinSize the typical magnetic moment of the spins
     * @param tol tolerance for convergence of the self-consistent calculation
     * @param J_ex exchange interaction parameter
     * @throws NullPointerException if one of the tables in given lattice (intTable, nnArray, exchangeIntTable, momentTable or energyTable) are {@code null}
     */
    public SingleTMonteCarloSimulation(final double T, final int temperatureIndex, final int totalNumOfTemperatures, final Lattice lattice, final int numOfObservables, final long maxSweeps,
                                       final long seed, final MersenneTwister rnd, final boolean continueFromSave, final OutputWriter out, final boolean checkpoint,
                                       final int maxIter, final double alpha, final BufferedWriter outProblematicConfigs,
                                       final double spinSize, final double tol, final double J_ex) {
        this.T = T;
        this.temperatureIndex=temperatureIndex;
        this.sweeps = 0;
        this.currentBinCount = 0;
        binAvg = new double[numOfObservables];
        this.lattice=lattice;
        this.maxSweeps=maxSweeps;
        this.seed=seed;
        this.continueFromSave=continueFromSave;
        this.rnd=rnd;
        this.outWriter=out;
        this.checkpoint=checkpoint;
        this.alpha=alpha;
        this.maxIter=maxIter;
        this.totalNumOfTemperatures=totalNumOfTemperatures;
        this.outProblematicConfigs=outProblematicConfigs;
        this.acceptanceRateCount=0;
        this.acceptanceRateSum=0;
        this.lastSwapAccepted=false;
        this.spinSize=spinSize;
        this.tol=tol;
        this.J_ex=J_ex;
        this.methodsUsed=new int[20];

        if (lattice.energyTable==null || lattice.momentTable==null || lattice.nnArray==null || lattice.exchangeIntTable==null || lattice.intTable==null) {
            throw new NullPointerException("The Lattice given to the SingleTMonteCarloSimulation has null some pointers. ");
        }
    }

    /**
     * Adds new observables to the ones that are currently being measured (extend the binAvg array)
     * put {@code Double.NaN} as the values that will remain there until the end of the current bin
     * @param numOfObservables new number of observables
     */
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

    /**
     * Swap the lattice and energy of this simulation with that of another
     * @param other simulation with which to replace the systems
     */
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

    /**
     * increments the number of accepted swaps for the current bin
     */
    public void incAcceptanceRateSum() {
        this.acceptanceRateSum++;
    }

    /**
     * increments the total number of attempted swaps for the current bin
     */
    public void incAcceptanceRateCount() {
        this.acceptanceRateCount++;
    }

    /**
     * Returns the current energy of the system
     * @return current energy
     */
    public double getCurrentEnergy() {
        return currentEnergy;
    }

    /**
     * Indicates whether the last PT swap attempt was accepted
     * @return
     */
    public boolean wasLastSwapAccepted() { return lastSwapAccepted; }

    /**
     * Sets whether the last PT swap attempt was accepted
     * @param swapAccepted acceptance status of the last swap attempt
     */
    public void setLastSwapAcceptance(final boolean swapAccepted) {
        this.lastSwapAccepted=swapAccepted;
        if (swapAccepted){
            this.incAcceptanceRateSum();
        }
        this.incAcceptanceRateCount();
    }

    /**
     * initializes the {@code usedMethods} array. must be done after deserializing the simulation from a checkpoint
     */
    public void initMethodsUsedArr(){
        if (this.methodsUsed == null){
            this.methodsUsed=new int[20];
        }
    }

    /**
     * Sets the {@code BufferedWriter} that is used to output configurations not solved self-consistently
     * @param outProblematicConfigs the writer that is to be used to output problematic configurations
     */
    public void setOutProblematicConfigs(final BufferedWriter outProblematicConfigs) {
        this.outProblematicConfigs = outProblematicConfigs;
    }

    /**
     * Gets the number of time a given method was successfully used to solve the self-consistent calcualtion
     * @param methodIndex index of the method to check
     * @return the total number of times the given method was successfully used.
     *          can be -1 in case previous data is missing due to starting an old simulation from
     *          checkpoint from before this data was collected.
     */
    public int getNumOfTimesMethodWasUsed(final int methodIndex) {
        if (methodsUsed != null) {
            return methodsUsed[methodIndex];
        } else {
            // can happen when starting from an old checkpoint
            return -1;
        }
    }

    /**
     * Performs a Monte Carlo metropolis step
     * @param lattice spin lattice
     * @param T Temperature
     * @param rnd random number generator
     * @return object that holds the energy difference due to the MC step, the new lattice object after the step and the method that was used to solve the
     * self-consistent calculation during the step. the returned {@code Lattice} will be {@code null} if the change is rejected, in which case the simulation
     * should return to the {@code Lattice} state that was saved before the call to this method.
     * @throws ConvergenceException when the self-consistent calculation does not converge with any of the available methods
     */
    public static deltaEnergyAndLattice metropolisStep(Lattice lattice, double T, MersenneTwister rnd, int maxIter, double alpha, final double tol) throws ConvergenceException
    {

        double initLongEnergy = lattice.getEnergy();

        // choose random spin
        int flippedSpin = rnd.nextInt(lattice.getN());

        int methodUsed = lattice.flipSpin(maxIter, tol, flippedSpin, alpha, rnd);

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
            return new deltaEnergyAndLattice(null, deltaEnergy, methodUsed);
        } else {    // change is accepted
            return new deltaEnergyAndLattice(lattice, deltaEnergy, methodUsed);
        }
    }

    /**
     * Print out a problematic configuration (one that was not solved self-consistently by any of the available methods)
     * using the given writer.
     * This method is synchronized on the given {@code out} since all parallel simulation use the same writer.
     * This method should not be called very often anyway, if it is being called often there are probably other issues
     * the need to be fixed besides possible concurrency related inefficiencies
     * @param lattice the lattice from before the spin flip
     * @param flippedSpin the spin that was flipped which caused the configuration to become 'problematic'
     * @param out the output writer that will be used to write this configuration
     */
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

    /**
     * Performs a Monte Carlo sweep (N MC steps with N being the total number of spins)
     */
    public void monteCarloSweep() {
        // a monte carlo sweep consists of N steps
        for (int step = 0; step < lattice.getN(); step++) {

            if (outWriter.isPrintProgress()) {
                // create a progress bar that indicates progress within the current MC sweep
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

            // save the current lattice
            Lattice tempLattice = new Lattice(lattice);
            Double deltaEnergy;
            deltaEnergyAndLattice energyAndLattice=null;
            try {
                // perform a MC step
                energyAndLattice = metropolisStep(lattice, T, rnd, maxIter, alpha, tol);

                if (energyAndLattice.getLattice()!=null){   // the change is accepted
                    deltaEnergy=energyAndLattice.getDeltaEnergy();
                    lattice=energyAndLattice.getLattice();
                    methodsUsed[energyAndLattice.getMethodUsed()]++;    // count methods used
                }else{                                      // the change is rejected
                    deltaEnergy=null;
                    methodsUsed[energyAndLattice.getMethodUsed()]++;    // count methods used
                }
            } catch (ConvergenceException e){   // the change was unsuccessful due to a ConvergenceException
                printProblematicConfig(tempLattice, e.getFlippedSpin(), outProblematicConfigs);
                deltaEnergy=null;
                methodsUsed[0]++;    // this signifies all available methods failed
                System.err.println("There was an error converging, had to abort metropolis step ("+step+ ", sweep="+sweeps+" T="+T+") and try another one. \n" + e.getMessage());
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
                arr = tempLatticeWithAllFields.getArray();  // this is not very efficient since all the spins are copied twice (once on this line and once two lines back),
                                                            // but it is not part of the regular run so not that terrible
            } else {
                arr = lattice.getArray(); // deep copy
            }
            for (int i = 0; i < arr.length; i++) {
                outWriter.writeObservablesPerSpin(arr[i].getN(), arr[i].getSpin(), arr[i].getSpinSize(), arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz());
            }
            outWriter.flush();
        }
    }

    /**
     * Performs measurements and saves them, either to the current bin or to the output buffer
     */
    public void writeObservables(){
        // get all measurements
        double m = lattice.getMagnetization();
        double[] temp = lattice.getMagneticFields();
        double[] tempSpinSizes = lattice.getSpinSizes();
        double[] mk2 = lattice.getMK2();    // m(k)^2, used later for correlation length calculation

        // Write all of the observables at this point (following a MC sweep). this is performed only if the output type is VERBOSE in outWriter
        outWriter.writeObservablesVerbose(sweeps, m ,currentEnergy ,temp[0] ,temp[1] ,temp[2] ,temp[3] ,temp[4] ,temp[5] ,temp[6] ,temp[7] , tempSpinSizes[0] ,tempSpinSizes[1] ,mk2[0], mk2[1], mk2[2], lastSwapAccepted);

        // count how many samples are in the current bin
        if (sweeps>0) currentBinCount++;

        if (currentBinCount==sweeps/2 + 1 && (sweeps&1)==1){	// last step of current bin
            // restart sums and counters for next bin
            double acceptanceRateForBin;	// acceptance rate for this bin & temperature

            if (acceptanceRateCount==0) {   // probably some error
                acceptanceRateForBin = -1;
            }else {
                // average acceptance rate for this bin
                acceptanceRateForBin = ((double) acceptanceRateSum) / acceptanceRateCount;
            }

            // print data from previous bin. this is performed only if the output type is BIN in outWriter
            outWriter.writeObservablesBin(currentBinCount, binAvg, acceptanceRateForBin);

            // restart sums and counters
            Arrays.fill(binAvg, 0);

            acceptanceRateCount=0;
            acceptanceRateSum=0;

            currentBinCount=0;
        }

        // add the current measurement of all of the observables to the bins
        addToAvg(binAvg,new double[]{Math.abs(m), m , m*m, currentEnergy , temp[0] , temp[1] , temp[2] , temp[3] , temp[4] , temp[5] , temp[6] , temp[7] , tempSpinSizes[0] , tempSpinSizes[1] , mk2[0], mk2[1], mk2[2]});


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

    /**
     * Adds an array of current measurement to the binned measurements array with their squares
     * @param avgArr array of binned measurements
     * @param currentValueArray array of current measurements
     */
    public static void addToAvg(double[] avgArr, double[] currentValueArray){
        // the length of avgArr should be twice that of currentValueArray since it holds all
        // of the measurements and also their squares
        if (avgArr.length>>1 != currentValueArray.length){
            throw new RuntimeException("error adding observable values to bin average. array length mismatch. ");
        }else{
            for (int i=0;i<currentValueArray.length;i++){
                // sum the measurements in the even cells
                avgArr[2*i]+=currentValueArray[i];
                // and their squares in the odd cells
                avgArr[2*i+1]+=currentValueArray[i]*currentValueArray[i];
            }
        }
    }

    /**
     * Calculates the standard error given the sum, sum^2 and the number of samples, N.
     * @param sum sum of samples
     * @param sum2 sum of squared samples
     * @param N number of samples
     * @return the standard error (error bar) of the given samples
     */
    public static double getStandardError(double sum, double sum2, long N){
        return Math.sqrt((sum2/N - (sum/N)*(sum/N))/N);
    }

    /**
     * Checks whether the given array of averages and standard errors of different observables each agree within error bars
     * @param arr Array of bin measurement averages and standard errors of different observables
     * @return whether the bin averages of each observable all agree (have overlapping error ranges)
     */
    public static boolean checkEquilibration(ArrayList<CircularFifoQueue<Pair<Double,Double>>> arr){
        if (arr.size()==0)
            return false;

        boolean equilibrated=true;
        for (CircularFifoQueue<Pair<Double, Double>> binsOfObservable : arr) {
            equilibrated &= checkEquilibration(binsOfObservable);
        }
        return equilibrated;
    }

    /**
     * Checks whether the given queue of a single observable's averages and standard deviations agree within error bars
     * @param binsOfObservable queue of the last x (usually 3) bin averages and standard errors
     * @return whether all bins agree (have overlapping error ranges)
     */
    public static boolean checkEquilibration(CircularFifoQueue<Pair<Double, Double>> binsOfObservable){
        boolean equilibrated=true;
        if (!binsOfObservable.isAtFullCapacity()){
            // queue has less than 3 bins
            equilibrated=false;
        } else {
            // get the highest lower bound and the lowest upper bound
            double upperBound = binsOfObservable.get(0).getLeft() + binsOfObservable.get(0).getRight();
            double lowerBound = binsOfObservable.get(0).getLeft() - binsOfObservable.get(0).getRight();
            for (Pair<Double, Double> bin : binsOfObservable) {
                upperBound = Math.min(upperBound, bin.getLeft() + bin.getRight());
                lowerBound = Math.max(lowerBound, bin.getLeft() - bin.getRight());
            }
            // if the lowest upper bound is higher that the highest lower bound, then all bins have some overlap
            equilibrated = upperBound > lowerBound;
            // also verify that there is no monotonicity, either increasing or decreasing
            // this is no longer a used condition, but this method is not used anyhow
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

    /**
     * Checks whether a given {@code boolean} array contains all {@code true} entries.
     * @param array input array to check if all entries are {@code true}
     * @return whether all entries in {@code array} are {@code true}
     */
    public static boolean isAllTrue(boolean... array)
    {
        for(boolean b : array) if(!b) return false;
        return true;
    }

    /**
     * Run the simulation for one MC sweep and perform all necessary outputs
     */
    public void run() {
        monteCarloSweep();
        writeObservables();
        sweeps++;

        if (sweeps%outWriter.getNumOfBufferedRows()==0 || sweeps==maxSweeps) {	// every obsPrintSweepNum sweeps or at the last one
            if (outWriter.isPrintOutputToConsole()) System.out.println("T="+T);
            outWriter.flush();
        }
    }

    /**
     * Generates a string with all of the parameters and their values used in this run
     * @param version software version
     * @param T temperature array for the containing {@link MultipleTMonteCarloSimulation}
     * @param extraMessage an extra text message to append to the list of parameters
     * @param mutualSeed the seed used by the containing {@link MultipleTMonteCarloSimulation}
     * @param tempScheduleFileName name of the temperature schedule file used to set the temperatures in the containing {@link MultipleTMonteCarloSimulation}
     * @param parallelTemperingOff whether to turn off parallel tempering (never swap temperatures)
     * @return a string containing all info about the parameters used in this run
     */
    public String runParameters(String version, double[] T, String extraMessage, long mutualSeed, String tempScheduleFileName, boolean parallelTemperingOff){
        return  "# VERSION: " + version + System.lineSeparator() +
                "#" + LocalDateTime.now() + System.lineSeparator() +
                "#temperature_schedule: " + Arrays.toString(T) + System.lineSeparator() +
                "#T=" + this.T + ":" + temperatureIndex + System.lineSeparator() +
                "#" + Constants.constantsToString() + System.lineSeparator() +
                String.format("# Lx=%s, Ly=%s, Lz=%s, x=%f, J_ex=%f, spinSize=%.8f, tol=%4.1e, extBx=%s, extBy=%s, maxSweeps=%s, suppressInternalTransFields=%s, " +
                        "continueFromSave=%s, maxIter=%s, bufferSize=%s, tempScheduleFileName=%s, parallelTemperingOff=%s, " +
                        "checkpoint=%s, folderName=%s, alpha=%s, output=%s ",lattice.getLx(),lattice.getLx(),lattice.getLz(),lattice.getConcentration(),J_ex,spinSize,tol,lattice.getExtBx(), lattice.getExtBy(), maxSweeps,lattice.isSuppressInternalTransFields(),
                    continueFromSave, maxIter, outWriter.getBufferSize(),
                    tempScheduleFileName, parallelTemperingOff, checkpoint, outWriter.getFolderName(),
                    alpha, outWriter.getOutType()) + System.lineSeparator() +
                "#" + Constants.locationsToString() + "," + Constants.latticeVectorsToString() + System.lineSeparator() +
                "#seed=" + mutualSeed + " (" + seed + ")" + System.lineSeparator() +
                extraMessage;
    }

    /**
     * Prints the parameters used in this run
     * @param version software version
     * @param T temperature array for the containing {@link MultipleTMonteCarloSimulation}
     * @param extraMessage an extra text message to append to the list of parameters
     * @param mutualSeed the seed used by the containing {@link MultipleTMonteCarloSimulation}
     * @param tempScheduleFileName name of the temperature schedule file used to set the temperatures in the containing {@link MultipleTMonteCarloSimulation}
     * @param parallelTemperingOff whether to turn off parallel tempering (never swap temperatures)
     * @throws IOException
     * @see #runParameters(String, double[], String, long, String, boolean)
     */
    public void printRunParameters(String version, double[] T, String extraMessage, long mutualSeed, String tempScheduleFileName, boolean parallelTemperingOff) throws IOException{
        // print some information to the beginning of the results file (or console):
        outWriter.print(runParameters(version, T, extraMessage, mutualSeed, tempScheduleFileName, parallelTemperingOff), true);
        outWriter.flush();
    }

    /**
     * Initializes this simulation
     */
    public void initSimulation(){
        boolean convergedConfig=false;
        int index=0;
        while (!convergedConfig) {
            lattice.randomizeConfig(rnd); // randomizes the spins and sets initial spin sizes as spinSize in the corresponding direction

            // initialize all the fields and then find all of the magnetic moment self-consistently
            lattice.updateAllLocalFields();
            lattice.updateAllMagneticMoments(2*maxIter, tol, alpha);

            // if the self-consistent calculation does not succeed, try again with a new random configuration until it does
            // if this happens too often there might be a problem
            if (lattice.magneticMomentConvergence() > tol){
                convergedConfig=false;
            }else{
                convergedConfig=true;
            }

            if (!convergedConfig) System.out.println("initial setup convergence index (T [num]="+T+" ["+temperatureIndex+']'+"): " + ++index);
        }

        currentEnergy = lattice.getEnergy();
    }

    /**
     * Creates queues used to save the last x (usually 3) bin averages to track equilibration
     * @param numOfBins number of bins to keep averages from
     * @param numOfObservables number of observables to track the equilibration of
     * @return
     */
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
