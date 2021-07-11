package simulation.montecarlo;

import org.apache.commons.math3.random.MersenneTwister;
import simulation.mmsolve.ConvergenceException;
import simulation.mmsolve.fi_xi;
import simulation.mmsolve.func;

import java.io.*;
import java.util.Arrays;
import java.util.Properties;

public class Lattice implements Serializable {
    // for serialization. should not be changed
    private static final long serialVersionUID = -9119463380760410942L;
    /** number of spins in the lattice */
    private final int N;
    /** linear system size in the x and z directions */
    private final int Lx, Lz;
    /** external field in the transverse (x and y) directions  */
    private final double extBx, extBy;
    /** concentration of Ho ions */
    private final double x;
    /** whether internal transverse fields are suppressed (offdiagonal dipolar terms excluded) */
    private final boolean suppressInternalTransFields;
    /** the typical magnetic moment of the spins */
    private final double spinSize;

    // after deserialization these must be set:
    /** table of pairwise interactions (dipolar and exchange) between spins */
    transient double[][][] intTable=null;
    /** table of pairwise interactions (only exchange) between spins. not in use */
    transient double[][] exchangeIntTable=null;
    /** Object used to obtain the magnetic moment of a single spin under a given applied magnetic field */
    transient FieldTable momentTable=null;
    /** Object used to obtain the energy of a single spin under a given applied magnetic field */
    transient FieldTable energyTable=null;
    /** table of the nearest neighbors of each spin */
    transient int[][] nnArray=null;
    /** Object that performs measurements on the lattice */
    private transient ObservableExtractor measure;
    /** Object that is used to solve the self-consistent calculation */
    private transient MagneticMomentsSolveIter iterativeSolver;

    /** The array of {@link singleSpin} */
    private singleSpin[] lattice;

    /**
     * Constructs a new Lattice object
     * @param Lx linear system size along the x direction
     * @param Lz linear system size along the z direction
     * @param x concentration of magnetic ions (x=1 means no dilution)
     * @param extBx,extBy external transverse field (x and y directions)
     * @param suppressInternalTransFields whether to suppress internal transverse fields (exclude offdiagonal dipolar terms)
     * @param spinSize the typical magnetic moment of the spins
     * @param dilution boolean array that indicates which spins exist (true) and which are diluted (false)
     * @param intTable table of pairwise interactions (dipolar and exchange) between spins
     * @param exchangeIntTable Object used to obtain the energy of a single spin under a given applied magnetic field
     * @param nnArray table of the nearest neighbors of each spin
     * @param energyTable Object used to obtain the energy of a single spin under a given applied magnetic field
     * @param momentTable Object used to obtain the magnetic moment of a single spin under a given applied magnetic field
     * @param measure Object that performs measurements on the lattice
     */
    @CreatesInconsistency("If intTable, exchangeIntTable, energyTable, momentTable or measure are null")
    public Lattice(int Lx, int Lz, double x, double extBx, double extBy, boolean suppressInternalTransFields, double spinSize, boolean[] dilution, double[][][] intTable, double[][] exchangeIntTable, int[][] nnArray, FieldTable energyTable, FieldTable momentTable, final ObservableExtractor measure){
        int numOfSpins=0;    // total number of spins in the (possibly diluted) system
        for (int i=0;i<dilution.length;i++){
                if (dilution[i]) numOfSpins++;
        }
        this.N=numOfSpins;
        this.Lx=Lx;
        this.Lz=Lz;
        this.x=x;
        this.extBx=extBx;
        this.extBy=extBy;
        this.suppressInternalTransFields=suppressInternalTransFields;
        this.spinSize = spinSize;
        this.intTable=intTable;
        this.exchangeIntTable=exchangeIntTable;
        this.energyTable=energyTable;
        this.momentTable=momentTable;
        this.nnArray=nnArray;
        this.measure=measure;
        this.iterativeSolver = new MagneticMomentsSolveIter();
        // generate a new array of singleSpins, initially with all magnetic moments set to spinSize
        this.lattice=generateIsingLattice(Lx,Lz,spinSize, dilution);
        if (intTable!=null && exchangeIntTable!=null && energyTable!=null && momentTable!=null && measure!=null)
            this.updateAllLocalFields();
    }

    /**
     * Constructs a new {@code Lattice} object with no dilution (x=1)
     * @see Lattice#Lattice(int, int, double, double, boolean, double, double[][][], double[][], int[][], FieldTable, FieldTable, ObservableExtractor)
     */
    @CreatesInconsistency("If intTable, exchangeIntTable, energyTable, momentTable or measure are null")
    public Lattice(int Lx, int Lz, double extBx, double extBy, boolean suppressInternalTransFields, double spinSize, double[][][] intTable, double[][] exchangeIntTable, int[][] nnArray, FieldTable energyTable, FieldTable momentTable, final ObservableExtractor measure){
        // if no dilution array is given, create full lattice
        this(Lx, Lz, 1.0, extBx, extBy, suppressInternalTransFields, spinSize, trueArray(Lx*Lx*Lz*Constants.num_in_cell), intTable, exchangeIntTable, nnArray, energyTable, momentTable, measure);
    }

    /**
     * Deep copy {@code Lattice} object
     * @param other lattice object to copy
     */
    public Lattice(Lattice other){
        this.N=other.N;
        this.Lx=other.Lx;
        this.Lz=other.Lz;
        this.x=other.x;
        this.extBx=other.extBx;
        this.extBy=other.extBy;
        this.suppressInternalTransFields=other.suppressInternalTransFields;
        this.spinSize=other.spinSize;
        this.intTable=other.intTable;
        this.exchangeIntTable=other.exchangeIntTable;
        this.energyTable=other.energyTable;
        this.momentTable=other.momentTable;
        this.nnArray=other.nnArray;
        this.measure=other.measure;
        this.iterativeSolver = new MagneticMomentsSolveIter();


        lattice = new singleSpin[other.lattice.length];
        for (int i=0;i<other.lattice.length;i++){
            lattice[i] = new singleSpin(other.lattice[i]);
        }
    }

    /**
     * Constructor that (deep) copies a given Lattice, but sets a new value for suppressInternalFields.
     * Used for measuring what the internal transverse fields would be had they not been suppressed.
     * @param other Lattice to copy
     * @param newSuppressInternalFields New value for suppressInternalFields
     */
    @CreatesInconsistency("By design. This constructor should be used only for measurement and then discarded.")
    public Lattice(Lattice other, boolean newSuppressInternalFields){
        this.N=other.N;
        this.Lx=other.Lx;
        this.Lz=other.Lz;
        this.x=other.x;
        this.extBx=other.extBx;
        this.extBy=other.extBy;
        this.suppressInternalTransFields=newSuppressInternalFields;
        this.spinSize=other.spinSize;
        this.intTable=other.intTable;
        this.exchangeIntTable=other.exchangeIntTable;
        this.energyTable=other.energyTable;
        this.momentTable=other.momentTable;
        this.nnArray=other.nnArray;
        this.measure=other.measure;
        this.iterativeSolver = new MagneticMomentsSolveIter();

        lattice = new singleSpin[other.lattice.length];
        for (int i=0;i<other.lattice.length;i++){
            lattice[i] = new singleSpin(other.lattice[i]);
        }
    }

    /**
     * Creates a {@code boolean} array of size N filled with {@code true}
     * @param N size of array to create
     * @return an array filled with {@code true}
     */
    public static boolean[] trueArray(int N){
        boolean[] arr = new boolean[N];
        for (int i=0; i<arr.length; i++) arr[i]=true;
        return arr;
    }

    public double getConcentration() {
        return x;
    }

    public void setIntTable(double[][][] intTable) {
        this.intTable = intTable;
    }

    public void setExchangeIntTable(double[][] exchangeIntTable) {
        this.exchangeIntTable = exchangeIntTable;
    }

    public void setMomentTable(FieldTable momentTable) {
        this.momentTable = momentTable;
    }

    public void setEnergyTable(FieldTable energyTable) {
        this.energyTable = energyTable;
    }

    public void setNnArray(int[][] nnArray) {
        this.nnArray = nnArray;
    }

    public void setMeasure(ObservableExtractor measure) {
        this.measure = measure;
    }

    public void initIterativeSolver() {
        this.iterativeSolver = new MagneticMomentsSolveIter();
    }

    public ObservableExtractor getMeasure() {
        return measure;
    }

    public MagneticMomentsSolveIter getIterativeSolver() {
        return iterativeSolver;
    }

    public int getN() {
        return N;
    }

    public int getLx() {
        return Lx;
    }

    public int getLz() {
        return Lz;
    }

    public double getExtBx() {
        return extBx;
    }

    public double getExtBy() {
        return extBy;
    }

    public boolean isSuppressInternalTransFields() {
        return suppressInternalTransFields;
    }

    /**
     * Get a deep copy of the spin array
     * @return a new {@code singleSpin} array
     */
    public singleSpin[] getArray() {
        return copyLattice(this.lattice);
    }

    /**
     * Checks whether a given spin exists (or is diluted)
     * @param n spin ordinal number to check
     * @return whether the given spin exists in the lattice
     */
    public boolean spinExists(int n){
        // not very efficient, so should be used only for compatibility validation
        boolean found=false;
        for (int i=0; !found && i < lattice.length; i++){
            found = found || lattice[i].getN() == n;
        }
        return found;
    }

    /**
     * Randomize the spin configuration
     * @param rnd a random number generator
     */
    @CreatesInconsistency
    public void randomizeConfig(MersenneTwister rnd){
        for (int i=0; i<lattice.length;i++){
            if (lattice[i].getSpin()!=0){
                if (rnd.nextBoolean()) lattice[i].setSpin(1);
                else lattice[i].setSpin(-1);
                lattice[i].setSpinSize(spinSize*lattice[i].getSpin());
            }
        }
    }

    /**
     * Creates "checker board" configuration, i.e. spins are set "up" and "down" alternately
     */
    @CreatesInconsistency
    public void checkerBoard(){
        for (int i=0; i<lattice.length;i++){
            lattice[i].setSpin(i%2==0 ? 1 : -1);
            lattice[i].setSpinSize(spinSize*lattice[i].getSpin());
        }
    }

    /**
     * Flip a given spin and find new magnetic moments self-consistently
     * @param maxIter - Maximum iterations for iterative solver (Gauss-Seidel)
     * @param tol - Tolerance for self-consistent solution. Given in terms of the average deviation per spin
     * @param flipSpin - Index of the spin to flip
     * @param alpha - Relaxation parameter for iterative solver.
     * @param rnd - Random number generator used for some solvers that need random initial states
     * @return the index of the last method used to solve the self-consistent calculation
     * @throws ConvergenceException when no convergence is achieved after trying all available methods.
     */
    public int flipSpin(int maxIter, double tol, int flipSpin, double alpha, MersenneTwister rnd) throws ConvergenceException {
        singleSpin[] tempLattice = Lattice.copyLattice(lattice);    // save original lattice. This is before the spin is actually flipped.

        boolean success=false;
        // order of methods for self-consistent calculation to try
        int[] methodsToTry = new int[]{1, 8, 2, 18, 3, 5, 15, 4, 13, 7, 16, 17, 10, 14, 6, 11, 12, 9, 19};

        // build error message to which we will append the errors as the accumulate
        String errorMessage="There was an error flipping spin " + flipSpin + ". The methods tried and the resulting errors are as follows:\n" +
                Arrays.toString(methodsToTry) + "\n";

        int methodIndex;
        // go over the list of methods and try them one by one in the given order
        for (methodIndex=0;!success && methodIndex<methodsToTry.length;methodIndex++){
            try {
                lattice = solveSelfConsistentCalc(maxIter, tol, flipSpin, methodsToTry[methodIndex], nnArray, alpha, rnd);
                success=true;
            } catch (ConvergenceException e){
                errorMessage += e.toString() + "\n";
                lattice=Lattice.copyLattice(tempLattice);
            }
        }

        // possible success
        if (!success) {
            throw new ConvergenceException(errorMessage, flipSpin);
        }
        // return the last method, which is the one that succeeded
        return methodsToTry[methodIndex-1];
    }

    /**
     * Calls the simple iterative solver
     * @param maxIter maximum number of iterations allowed
     * @param tol tolerance for convergence
     * @param alpha relaxation parameter of the update step
     * @see MagneticMomentsSolveIter#updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
     */
    public void updateAllMagneticMoments(int maxIter, double tol, double alpha){
        iterativeSolver.updateAllMagneticMoments(maxIter, tol, alpha, true);
    }

    /**
     * Wrapper method that solves the self-consistent calculation following a spin-flip
     * @param maxIter maximum number of iterations allowed
     * @param tol tolerance for convergence testing
     * @param flipSpin the flipped spin
     * @param method method number to use
     * @param nnArray array of nearest neighbors (used for determining the update order)
     * @param alpha relaxation parameter of the update step
     * @param rnd random number generator to create random initial guesses
     * @return the lattice object with the solved self-consistent calculation
     * @throws ConvergenceException if the given method failed to converge to a solution within the given tolerance
     */
    private singleSpin[] solveSelfConsistentCalc(int maxIter, double tol, int flipSpin, int method, int[][] nnArray, double alpha, MersenneTwister rnd) throws ConvergenceException {
        // method descriptions:
        // 1 - regular iterative solver
        // 2 - iterative solver with gradual flipping of the spin
        // 3 - iterative solver with gradual flipping of the spin that decouples the spin from the environment
        // for 4-8 we use the previous configuration as an initial guess
        // 4 - Newton's method
        // 5 - Broyden's method, Jacobian is identity and not changing
        // 6 - Broyden's method, Jacobian is initialized and not changing
        // 7 - Broyden's method, Jacobian is initialized and changed to identity upon failure
        // 8 - Broyden's method, Jacobian is identity and initialized upon failure
        // 9-13: same as 4-8 but with a random initial guess
        // 14-18: same as 4-8 but with an initial guess that is the typical spinSize with the sign based on the spin orientation
        // 19 - iterative solver with gradual flipping of the spin that makes sure the magnetic moment is always
        // well-defined as a function of the applied field

        if ((method>=1 && method<=3) || method==19){
            int[] bfsOrder=null;
            if (nnArray!=null) {
                bfsOrder = orderBFS(lattice.length, flipSpin);
            }

            if (method==1){
                // regular iterative method
                lattice[flipSpin].flipSpin();
                iterativeSolver.updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, false);

                double convergenceDistance = magneticMomentConvergence();
                if (convergenceDistance > tol){
                    throw new ConvergenceException.Builder("max number of iterations performed without reaching convergence.", "regular")
                            .setConvergenceDistance(convergenceDistance)
                            .build();
                }else{
                    return lattice;
                }

            }else if (method==2){
                // homotopic (iterative) method

                iterativeSolver.homotopicSolve(bfsOrder, maxIter, tol, flipSpin, alpha);

            }else if (method==3){	// method==3
                // homotopic (iterative) method

                iterativeSolver.homotopicSolve2(bfsOrder, maxIter, tol, flipSpin, alpha);
            }else{  //method==19
                iterativeSolver.homotopicSolve3(bfsOrder, maxIter, tol, flipSpin, alpha);
            }
        }else if (method>=4 && method<=18) {
            // broyden's or newton's method

            lattice[flipSpin].flipSpin();

            double[] x = new double[lattice.length];

            // initialize x:
            if (method <= 8) {
                // method numbers in the range 4<=method<=8 use the given initial guess
                for (int j = 0; j < lattice.length; j++) x[j] = lattice[j].getSpinSize();
            } else if (method >= 9 && method <= 13) {
                // method numbers in the range 9<=method<=13 use a random initial guess
                for (int j = 0; j < x.length; j++) x[j] = spinSize * rnd.nextDouble() * (rnd.nextBoolean() ? 1 : -1);
                method -= 5;
            } else if (method >= 14 && method <= 18) {
                // method numbers in the range 14<=method<=18 use a "good" initial guess
                for (int j = 0; j < x.length; j++) x[j] = spinSize * lattice[j].getSpin();
                method -= 10;
            }

            fi_xi funcToSolve = new func(intTable, momentTable, lattice, extBx, extBy, suppressInternalTransFields);

            // newton's method
            if (method == 4) {
                x = simulation.mmsolve.MagneticMomentsSolve.newt(funcToSolve, x, tol);
            }
            // broyden's method, Jacobian is identity and not changing
            else if (method == 5) {
                x = simulation.mmsolve.MagneticMomentsSolve.broyden(funcToSolve, x, tol, true, false);
            }
            // broyden's method, Jacobian is initialized and not changing
            else if (method == 6) {
                x = simulation.mmsolve.MagneticMomentsSolve.broyden(funcToSolve, x, tol, false, false);
            }
            // broyden's method, Jacobian is initialized and changed to identity upon fail
            else if (method == 7) {
                x = simulation.mmsolve.MagneticMomentsSolve.broyden(funcToSolve, x, tol, false, true);
            }
            // broyden's method, Jacobian is identity and initialized upon fail
            else if (method == 8) {
                x = simulation.mmsolve.MagneticMomentsSolve.broyden(funcToSolve, x, tol, true, true);
            }

            for (int j = 0; j < lattice.length; j++) lattice[j].setSpinSize(x[j]);

            this.updateAllLocalFields();

            if (magneticMomentConvergence() > tol) { // not converged
                ConvergenceException e = new ConvergenceException.Builder("solver exited before reaching convergence. ", method == 3 ? "Newton" : "Broyden")
                        .setConvergenceDistance(magneticMomentConvergence())
                        .setFlippedSpin(flipSpin)
                        .build();
                if (funcToSolve instanceof func) {
                    e.setNumManualCalc(((func) funcToSolve).getNumManualCalc());
                }
                throw e;
            } else {
                return lattice;
            }
        }else{
            System.err.println("Method number does not exist.");
            throw new ConvergenceException.Builder("Method number does not exist.", "solveSelfConsistentCalc")
                    .build();
        }
        return lattice;
    }

    /**
     * Returns new int array which has the breadth-first order with respect to the given spin index
     * @param n is the maximum (arr.length)
     * @param root root of the tree, i.e. the flipped spin from which we start the breadth-first search
     * @return array that contains the spin numbers in the BFS order
     */
    protected int[] orderBFS(int n, int root) {
        int[] bfs = new int[n];
        int head = 0, tail = 0;             // head and tail of the queue
                                            // everything before head is already ordered
        boolean[] visited = new boolean[n]; // holds which nodes were visited and queued

        // queue the root
        bfs[tail++] = root;
        visited[root] = true;

        while (head < tail) {
            // dequeue v
            int v = bfs[head++];

            //get v's neighbors
            for (int neighbor = 0; neighbor < nnArray[v].length; neighbor++) {
                // if neighbor exists and was not already visited and not all spins were added
                if (nnArray[v][neighbor] >= 0 && !visited[nnArray[v][neighbor]] && tail < n) {
                    bfs[tail++] = nnArray[v][neighbor];
                    visited[nnArray[v][neighbor]] = true;
                }
            }
        }

        if (head < n){  // not all spins were added to the queue, which could happen if, due to dilution,
                        // the graph is not fully connected. then we just add them at the end.
            for (int i=0; i<n; i++){
                if (!visited[i]){
                    bfs[head++] = i;
                }
            }
        }
        /*
        // check all spins were queued
        boolean verify=true;
        for (int i=0;i<visited.length && verify;i++)  verify = verify && visited[i];

        // check all spins in arr are present in the new bfs
        for (int i=0;i<n && verify;i++){
            verify = verify && ArrayUtils.contains(bfs, i);
        }
        if (!verify){
            System.err.println("error with orderBFS. not all spins were inserted into the ordered array.");
            System.exit(1);
        }
        */
        return bfs;
    }

    /**
     * Deep copies the given {@code singleSpin} array
     * @param arr array to deep copy
     * @return copy of the given array
     */
    private static singleSpin[] copyLattice(singleSpin[] arr){
        singleSpin[] newLatticeArr = new singleSpin[arr.length];
        for (int i=0;i<arr.length;i++){
            newLatticeArr[i] = new singleSpin(arr[i]);
        }
        return newLatticeArr;
    }

    /**
     * Deep copies of this {@code Lattice} object
     * @return a deep copy of this {@code Lattice} object
     */
    public Lattice getCopy(){
        return new Lattice(this);
    }

    /**
     * Updates all of the local fields according to the current magnetic moments
     */
    public void updateAllLocalFields(){
        int i;

        for (i=0;i<lattice.length;i++) {
            int j;
            if (lattice[i].getSpin()!=0) {
                double Bz=0, Bx = extBx, By = extBy;
                for(j=0;j<lattice.length;j++){
                    if (lattice[j].getSpin()!=0){
                        Bz += lattice[j].getSpinSize()*intTable[2][i][j];
                        if (!suppressInternalTransFields){
                            Bx += (lattice[j].getSpinSize()) * intTable[0][i][j];
                            By += (lattice[j].getSpinSize()) * intTable[1][i][j];
                        }
                    }
                }
                lattice[i].setLocalBz(Bz);
                // if suppressInternalTransFields==true then Bx==extBx and By==extBy.
                lattice[i].setLocalBx(Bx);
                lattice[i].setLocalBy(By);
            }
        }
    }

    /**
     * Calculates the distance from convergence. Parameters are used when using a method which mixes
     * the flipped spin and therefore has a different target magnetic moment for that spin.
     * @param flipSpin the spin that is being flipped
     * @param frac at what fraction to mix the current spin state with the opposite spin state
     * @param smoothHomotopy whether to mix the spin states in a smooth way, such that the
     *                       magnetic moment of "up" becomes that of "down" with no singularities,
     *                       or naively as a weighted average
     * @return distance from convergence
     */
    private double magneticMomentConvergence(int flipSpin, double frac, boolean smoothHomotopy){
        double converged = 0;
        double max = 0;

        for (int i = 0; i < N; i++) {
            if (i!=flipSpin) {
                double shouldBe = momentTable.getValue(lattice[i].getLocalBx(), lattice[i].getLocalBy(), lattice[i].getLocalBz(), lattice[i].getSpin(), lattice[i].getPrevBIndices());
                double is = lattice[i].getSpinSize();

                converged += Math.abs(is - shouldBe);
            }else{
                double shouldBe = momentTable.getValue(lattice[i].getLocalBx(), lattice[i].getLocalBy(), lattice[i].getLocalBz(), lattice[i].getSpin(), lattice[i].getPrevBIndices(), frac, smoothHomotopy);
                double is = lattice[i].getSpinSize();

                converged += Math.abs(is - shouldBe);
            }

        }

        return Math.abs(converged/lattice.length);
    }

    /**
     * Calculates the distance from convergence, with smoothHomotopy=false.
     * @see #magneticMomentConvergence(int, double, boolean)
     */
    private double magneticMomentConvergence(int flipSpin, double frac){ return magneticMomentConvergence(flipSpin, frac, false); }

    /**
     * Calculate the magnetic moment convergence without any fractional (mixed) spins.
     * @see #magneticMomentConvergence(int, double, boolean)
     */
    public double magneticMomentConvergence(){
        return magneticMomentConvergence(-1, 0);
    }

    /**
     * Calculate the magnetic moment convergence without any fractional (mixed) spins.
     * flipSpin is ignored (which is required for homotopicSolve2), i.e. stays in the starting orientation
     * @see #magneticMomentConvergence(int, double, boolean)
     */
    private double magneticMomentConvergence(int flipSpin){
        return magneticMomentConvergence(flipSpin, 0);
    }

    /**
     * Update all magnetic fields after a specific change of one magnetic moment
     * @param i the index of the spin whose magnetic moment was just changed
     * @param prevSpinSize the previous magnetic moment (before the change)
     */
    private void updateFieldsAfterSpinSizeChange(int i, double prevSpinSize){
        int j;
        double deltaSpinSize = lattice[i].getSpinSize()-prevSpinSize;
        if (lattice[i].getSpin()!=0) {
            // change local fields (longitudinal & transverse) at all other spins by
            // removing the contribution of i there before the change and
            // adding the updated contribution
            for(j=0;j<lattice.length;j++){
                if (lattice[j].getSpin()!=0){
                    // remove interaction with prevSpinSize and add interaction with current spinSize
                    if (!suppressInternalTransFields) {
                        lattice[j].setLocalBx(lattice[j].getLocalBx() + (deltaSpinSize) * intTable[0][i][j]);
                        lattice[j].setLocalBy(lattice[j].getLocalBy() + (deltaSpinSize) * intTable[1][i][j]);
                    }
                    lattice[j].setLocalBz(lattice[j].getLocalBz() + (deltaSpinSize)*intTable[2][i][j]);

                }
            }
        }
    }

    /**
     * Verifies that all local fields are in agreement with the current mgnetic moments (within 1.0e-10)
     * @return whether or not the current local fields are consistent w/ the current magnetic moments
     */
    public boolean verifyAllLocalFields(){
        boolean ret=true;
        double tol = 1.0e-10;
        int i;
        for (i=0;i<lattice.length && ret;i++) {
            int j;

            double Bx=extBx, By=extBy, Bz=0;
            // add up all contributions to local field from other spins at spin i
            for(j=0;j<lattice.length;j++){
                if (lattice[j].getSpin()!=0){
                    // calculate interaction
                    if (!suppressInternalTransFields) {
                        Bx += (lattice[j].getSpinSize() * intTable[0][i][j]);
                        By += (lattice[j].getSpinSize() * intTable[1][i][j]);
                    }
                    Bz += (lattice[j].getSpinSize()*intTable[2][i][j]);
                }
            }

            // find max difference
            /*
            if (Math.abs(lattice[i].getLocalBx()-Bx)>max) {
                max = Math.abs(lattice[i].getLocalBx() - Bx);
            }
            if (Math.abs(lattice[i].getLocalBy()-By)>max) {
                max = Math.abs(lattice[i].getLocalBy() - By);
            }
            if (Math.abs(lattice[i].getLocalBz()-Bz)>max) {
                max = Math.abs(lattice[i].getLocalBz() - Bz);
            }
            //System.out.println("current Bz: " + arr[i].getLocalBz());
            //System.out.println("calc'd Bz: " + Bz);
            */
            // check if Bx,By or Bz are not equal to the values stored in the i-th site up to a tolerance constant
            ret &= (Math.abs(lattice[i].getLocalBx()-Bx)<tol) && (Math.abs(lattice[i].getLocalBy()-By)<tol) && (Math.abs(lattice[i].getLocalBz()-Bz)<tol);
        }

        if (!ret){
            System.err.println("bad field at spin " + i);
        }

        //System.out.println("max difference from real local field is: " + max);
        return ret;
    }

    /**
     * Reads a lattice from file and creates initial random configuration
     * Only required when there is disorder, i.e. dilution<1.0.
     * @param fileNumber - as in "config_%Lx%_%Lz%_%x%_%h%_%fileNumber%.txt
     * @param dilution - dilution of Ising lattice
     * @param Lx - number of unit cells in x direction of the lattice
     * @param Lz  - number of unit cells in z direction of the lattice
     * @param h - standard deviation of the random local field
     * @param rnd - MersenneTwister PRNG. Used to make the random starting configuration
     * @return Received Ising lattice
     * @throws IOException
     * @deprecated since the dilution profile is created at the beginning of the simulation and is received instead of created here.
     *              Use {@link #generateIsingLattice(int, int, double, boolean[])} instead.
     */
    @Deprecated
    static singleSpin[] receiveIsingLattice(int fileNumber, double dilution, int Lx, int Lz, double h, MersenneTwister rnd) throws IOException {
        final double spinSize = CrystalField.getMagneticMoment(0.0, 0.0, 0.05);

        singleSpin[] arr = null;
        try (BufferedReader in = new BufferedReader(new FileReader(System.getProperty("system") + File.separator + "data" + File.separator + "configurations" + File.separator + "config_"+Lx+"_"+Lz+"_"+dilution+"_"+h+"_"+fileNumber+".txt"))){;
            String str;
            String[] params;
            int i=0;


            if ((str = in.readLine()) != null)
                Lx=Integer.parseInt(str.split("=")[1]);
            if ((str = in.readLine()) != null)
                Lz=Integer.parseInt(str.split("=")[1]);

            arr = new singleSpin[Constants.num_in_cell*Lx*Lx*Lz];

            in.readLine();	// skip the line with the seed

            while ((str = in.readLine()) != null){
                params = str.split(",");
                try{
                    if (Integer.parseInt(params[3])!=0){	// meaning that there is a spin there
                        if (rnd.nextBoolean())	// spin +1
                            arr[i]=new singleSpin(1,0, spinSize);
                        else	// spin -1
                            arr[i]=new singleSpin(-1,0, spinSize);
                    }else{ // initialize spin but with s=0
                        arr[i]=new singleSpin(0, 0, spinSize);
                    }
                }catch(NumberFormatException nfe){
                    System.err.println("error reading line "+i+" in file "+ fileNumber+". numeric value expected.");
                    throw new RuntimeException("Error reading config file (see errors)");
                }
                i++;
            }

        } catch (IOException e) {
            System.err.println("bad input file!");
            throw new RuntimeException();
        }

        return arr;
    }

    /**
     * Generates a singleSpin array of the system
     * @param Lx Number of unit cells in x and in y directions
     * @param Lz Number of unit cells in z direction
     * @param spinSize the typical magnetic moment to assign to each spin in the returned array
     * @param dilution a {@code boolean} array that represents the dilution profile, i.e. which spins exist and which do not
     * @return An array of spins that represent a Lx*Lx*Lz crystal
     *          with all local fields set to zero and all spin sizes set to the given value.
     */
    public static singleSpin[] generateIsingLattice(int Lx, int Lz, final double spinSize, boolean[] dilution) {
        int i, j, k, l;

        // count real number of spins in lattice
        int N=0;
        for (int s=0; s< dilution.length; s++){
            if (dilution[s]) N++;
        }

        // create the array that will hold the lattice
        singleSpin[] arr = new singleSpin[N];

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Properties params = GetParamValues.getParams();
        double[][] location = new double[Constants.num_in_cell][3];	// 3D coordinate location for each of the atoms in the basis

        // fill location:
        for (l=0;l<location.length;l++){
            location[l]=GetParamValues.getLocation(params, l);
        }
        params=null;
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (i = 0; i < Lx; i++)
        {
            for (j = 0; j < Lx; j++)
            {
                for (k = 0; k < Lz; k++)
                {
                    // the spins in each unit cell are designated 0-3. see documentation for further info
                    // a lattice is created with all spins +1
                    for (l = 0; l < Constants.num_in_cell; l++)
                    {
                        int fullArrayIndex = i*Lx*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+l;
                        if (dilution[fullArrayIndex]) {
                            arr[getCompactArrayIndex(dilution, fullArrayIndex)] = new singleSpin(1, fullArrayIndex, spinSize);
                        }
                    }
                }
            }
        }

        return arr;
    }

    /**
     * Get the index in the compact spin array that corresponds to the given index in the full spin array.
     * @param dilution - array that holds the dilution config
     * @param index - index in the full array
     * @return corresponding index in the compact array
     * @throws IndexOutOfBoundsException when {@code index} is larger than the size of the given dilution array
     */
    public static int getCompactArrayIndex(boolean[] dilution, int index){
        if (index >= dilution.length) throw new IndexOutOfBoundsException("Index given is outside the bounds of the dilution array.");
        int count=0;
        for (int i=0; i<index; i++){
            if (dilution[i]){
                count++;
            }
        }
        return count;
    }


    // ********************************************   measuremets   ****************************************************
    public double[] getMagneticFields(){
        return measure.meanField(lattice);
    }

    public double[] getMK2(){ return measure.calc_mk2(lattice); }

    public double getMagnetization(){
        return measure.calcMagnetization(lattice);
    }

    public double getEnergy(){
        return measure.calcEnergy(this);
    }

    public double[] getSpinSizes(){
        return measure.calcSpinSizes(lattice);
    }

    public double getTransverseFieldMaximizingNNConfigsFrac(){
        return measure.countTransverseFieldMaximizingNNConfigs(this);
    }

    // *****************************************************************************************************************

    /**
     * Class of methods used to find a self-consistent solution for the magnetic moments
     * All methods here should be considered as ruining the singleSpin[] lattice in case an exception is thrown,
     * hence a copy should be made before calling them.
     */
    private class MagneticMomentsSolveIter {
        /**
         * Basic iterative solver, w/o homotopy features
         * @see #updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
         */
        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, boolean retryWithLowerAlpha) {
            updateAllMagneticMoments(sweepOrder, maxIter, tol, alpha, -1, false, 0, retryWithLowerAlpha, false);
        }

        /**
         * Basic iterative solver, w/o recursive call and w/ the flipped spin decoupled from the environment (stays at its initial state regardless of the field)
         * @see #updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
         */
        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, int flipSpin) {
            updateAllMagneticMoments(sweepOrder, maxIter, tol, alpha, flipSpin, false, 0, false, false);
        }

        /**
         * Basic iterative solver, w/o special ordering, w/o homotopy and w/o recursive call
         * @see #updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
         */
        public void updateAllMagneticMoments(int maxIter, double tol, double alpha) {
            updateAllMagneticMoments(null, maxIter, tol, alpha, -1, false, 0, false, false);
        }

        /**
         * Basic iterative solver, w/o special ordering
         * @see #updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
         */
        public void updateAllMagneticMoments(int maxIter, double tol, double alpha, int flipSpin, double frac) {
            updateAllMagneticMoments(null, maxIter, tol, alpha, flipSpin, false, frac, false, false);
        }

        /**
         * Basic iterative solver, w/o special ordering and w/o homotopy features
         * @see #updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
         */
        public void updateAllMagneticMoments(int maxIter, double tol, double alpha, boolean retryWithLowerAlpha) {
            updateAllMagneticMoments(null, maxIter, tol, alpha, -1, false, 0, retryWithLowerAlpha, false);
        }

        /**
         * Basic iterative solver, w/o smooth homotopy
         * @see #updateAllMagneticMoments(int[], int, double, double, int, boolean, double, boolean, boolean)
         */
        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, int flipSpin, boolean print, double frac) {
            updateAllMagneticMoments(sweepOrder, maxIter, tol, alpha, flipSpin, print, frac, false, false);
        }

        /**
         * Basic iterative solver, performs a (nonlinear) Gauss-Seidel scheme, i.e.,
         * 1. update a spin's magnetic moment based on the applied local field
         * 2. update the fields at all other sites based on the previous change
         * 3. go to the next spin
         * 4. loop through all spins in this way
         * 5. repeat until convergence is achieved
         * *** Important ***
         * This method does not verify convergence in the same way as other methods.
         * It stops when successive changes become small, which might indicate convergece,
         * but convergence must be checked in the traditional way after it is called.
         * @param sweepOrder the order in which we loop through the spins
         * @param maxIter maximum number of iterations over the entire system
         * @param tol tolerance for convergence
         * @param alpha relaxation parameter for the update in #1 (for values less than 1, we update the magnetic moment
         *              to a weighted average of the previous magnetic moment and the new magnetic moment
         * @param flipSpin the spin that is being flipped (used in the various homotopic solvers: {@link #homotopicSolve(int[], int, double, int, double)}, {@link #homotopicSolve2(int[], int, double, int, double)}, {@link #homotopicSolve3(int[], int, double, int, double)} )
         * @param print whether to print the progress. useful for debugging
         * @param frac the fraction by which to mix the previous spin with the new spin (used in the various homotopic methods)
         * @param retryWithLowerAlpha whether to try again with a lower (half) value of {@code alpha} and higher (double) number of {@code maxIter} in case convergence is not achieved after the given number of iterations
         * @param smoothHomotopy whether to change the magnetic moment's target function in a manner that does not make it ill-defined at any stage
         */
        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, int flipSpin, boolean print, double frac, boolean retryWithLowerAlpha, boolean smoothHomotopy) {
            boolean converged = false;
            int iter;
            int index;
            double sum;    // for real time convergence check
            double prevSpinSize;
            double newSpinSize;
            for (iter = 0; iter < maxIter && !converged; iter++) {

                sum = 0;    // for real-time convergence check
                for (int i = 0; i < lattice.length; i++) {
                    if (sweepOrder == null) {
                        index = i;
                    } else {
                        index = sweepOrder[i];
                    }

                    if (index != flipSpin) {
                        prevSpinSize = lattice[index].getSpinSize();
                        newSpinSize = momentTable.getValue(lattice[index].getLocalBx(), lattice[index].getLocalBy(), lattice[index].getLocalBz(), lattice[index].getSpin(), lattice[index].getPrevBIndices());
                        lattice[index].setSpinSize(prevSpinSize * (1 - alpha) + alpha * newSpinSize);
                        updateFieldsAfterSpinSizeChange(index, prevSpinSize);
                        sum += Math.abs(prevSpinSize - lattice[index].getSpinSize());
                    } else {
                        prevSpinSize = lattice[index].getSpinSize();
                        newSpinSize = momentTable.getValue(lattice[index].getLocalBx(), lattice[index].getLocalBy(), lattice[index].getLocalBz(), lattice[index].getSpin(), lattice[i].getPrevBIndices(), frac, smoothHomotopy);
                        lattice[index].setSpinSize(prevSpinSize * (1 - alpha) + alpha * newSpinSize);
                        updateFieldsAfterSpinSizeChange(index, prevSpinSize);
                        sum += Math.abs(prevSpinSize - lattice[index].getSpinSize());
                    }

                }

                if (print) {
                    System.out.print((sum / lattice.length) + " ");
                }

                if (sum / lattice.length < tol) {
                    converged = true;
                }


            }
            if (print) {
                System.out.println();
            }
            if (!converged && retryWithLowerAlpha) {
                updateAllMagneticMoments(sweepOrder, 2 * maxIter, tol, 0.5 * alpha, flipSpin, print, frac, false, smoothHomotopy);
            }
        }

        /**
         * Removes the contribution of a specific spin to the fields at all other sites
         * @param spinToRemove the index of the spin whose contribution is to be removed
         * @param magneticMomentToRemove the magnetic moment to remove from the other spins (that of {@code spinToRemove})
         */
        public void removeSpecificSpinConstibution(int spinToRemove, double magneticMomentToRemove) {
            int j;
            if (lattice[spinToRemove].getSpin() != 0) {
                // change local fields (longitudinal & transverse) at all other spins by
                // removing the contribution of i there
                for (j = 0; j < lattice.length; j++) {
                    if (lattice[j].getSpin() != 0) {
                        // remove interaction with magneticMomentToRemove
                        if (!suppressInternalTransFields) {
                            lattice[j].setLocalBx(lattice[j].getLocalBx() - magneticMomentToRemove * intTable[0][spinToRemove][j]);
                            lattice[j].setLocalBy(lattice[j].getLocalBy() - magneticMomentToRemove * intTable[1][spinToRemove][j]);
                        }
                        lattice[j].setLocalBz(lattice[j].getLocalBz() - magneticMomentToRemove * intTable[2][spinToRemove][j]);
                    }
                }
            }
        }

        /**
         * Add to the field at all sites a contribution of some spin with a given magnetic moment
         * @param flipSpin the spin whose contribution is to be added to all other sites
         * @param magneticMoment the assumed magnetic moment of {@code flipSpin}. usually a some mixture of the previous and new spin orientation
         */
        public void addMixedSpinField(int flipSpin, double magneticMoment) {
            int i;

            for (i = 0; i < lattice.length; i++) {
                if (lattice[i].getSpin() != 0) {
                    lattice[i].setLocalBz(lattice[i].getLocalBz() + magneticMoment * intTable[2][i][flipSpin]);
                    if (!suppressInternalTransFields) {
                        lattice[i].setLocalBx(lattice[i].getLocalBx() + magneticMoment * intTable[0][i][flipSpin]);
                        lattice[i].setLocalBy(lattice[i].getLocalBy() + magneticMoment * intTable[1][i][flipSpin]);
                    }
                }
            }
        }

        /**
         * Solves the self-consistent calculation by gradually flipping the given spin
         * @param bfsOrder the order in which we loop through the spins
         * @param maxIter maximum number of iterations over the entire system
         * @param tol tolerance for convergence
         * @param flipSpin the spin that is being flipped
         * @param alpha relaxation parameter for the update in #1 (for values less than 1, we update the magnetic moment
         *              to a weighted average of the previous magnetic moment and the new magnetic moment
         * @throws ConvergenceException
         */
        public void homotopicSolve(int[] bfsOrder, int maxIter, double tol, int flipSpin, double alpha) throws ConvergenceException {
            // this method tries to gradually change the target magnetic moment of the spin given as flipSpin.
            // we start by trying the full change, meaning going directly to the opposite spin orientation,
            // and only if the self-consistent calculation does not converge, then we try going half-way.
            // if there is still no convergence, we try half of that, etc.

            int numOfStepsLeft = 0; // number of steps left until the given spin is in its required state
            double perc = 1.0;      // percentage in the process of flipping (1.0 = final state, 0.0 = initial state)
            double prevPerc = 0.0;  // we keep it so we can fall back to it and try a smaller step in case of failure
            int i = 0;

            while (numOfStepsLeft >= 0) {
                i++;
                singleSpin[] tempLattice = Lattice.copyLattice(lattice);    // save the lattice in case we need to fall back to this step
                removeSpecificSpinConstibution(flipSpin, lattice[flipSpin].getSpinSize());  // remove the contribution of flipSpin from all other sites
                // create a magnetic moment value that is a mixture of the current and next orientations of the spin
                double magneticMoment = (1 - perc) * momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(),
                        lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices()) + perc * momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), -1 * lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());
                // add this mixed value to other sites
                addMixedSpinField(flipSpin, magneticMoment);
                // set this value as the magnetic moment of flipSpin
                lattice[flipSpin].setSpinSize(magneticMoment);
                // at this point the field should be consistent with the magnetic moments (though the magnetic moment of flipSpin is not a consistent with the applied field)
                updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, flipSpin, false, perc);

                if (magneticMomentConvergence(flipSpin, perc) > tol) {
                    // if this step did not converge, backtrack to the previous state and try half of the previous step
                    numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                    maxIter = (int) (1.1 * maxIter);
                    lattice = tempLattice;
                    perc = prevPerc;
                }
                prevPerc = perc;
                if (numOfStepsLeft > 0) perc = perc + (1.0 - perc) / numOfStepsLeft;
                numOfStepsLeft--;

                if (numOfStepsLeft < 0) { // last step
                    lattice[flipSpin].flipSpin();
                }

                // too many iterations, so this is not going anywhere (unless this was the last step):
                if ((i > 100 || numOfStepsLeft > Math.pow(2, 5)) && numOfStepsLeft >= 0) {
                    lattice[flipSpin].flipSpin();   // this does nothing. it is only so that magneticMomentConvergence that is called
                                                    // for the error information is run correctly.
                    ConvergenceException e = new ConvergenceException.Builder("Too many homotopic step taken without convergence. numOfStepsLeft=" + numOfStepsLeft +
                    ", i=" + i, "homotopic")
                            .setIndex(i)
                            .setConvergenceDistance(magneticMomentConvergence())
                            .build();

                    numOfStepsLeft = -1;  // this ends the loop (just in case the throw is removed)
                    throw e;
                }
            }
        }

        /**
         * Solves the self-consistent calculation by gradually flipping the given spin in a more sophisticated way than {@link #homotopicSolve(int[], int, double, int, double)}
         * @param bfsOrder the order in which we loop through the spins
         * @param maxIter maximum number of iterations over the entire system
         * @param tol tolerance for convergence
         * @param flipSpin the spin that is being flipped
         * @param alpha relaxation parameter for the update in #1 (for values less than 1, we update the magnetic moment
         *              to a weighted average of the previous magnetic moment and the new magnetic moment
         * @throws ConvergenceException
         */
        public void homotopicSolve2(int[] bfsOrder, int maxIter, double tol, int flipSpin, double alpha) throws ConvergenceException {
            // this method gradually flips the given spin but also "decouples" it from its environment in the process

            // remove the contribution of flipSpin from all other sites
            removeSpecificSpinConstibution(flipSpin, lattice[flipSpin].getSpinSize());

            int numOfStepsLeft = 0; // number of steps left until the given spin is in its required state
            double perc = 1.0;      // percentage in the process of flipping (1.0 = final state, 0.0 = initial state)
            double prevPerc = 0.0;  // we keep it so we can fall back to it and try a smaller step in case of failure
            int i = 0;

            while (numOfStepsLeft >= 0) {
                i++;
                double prevStateMoment = momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());
                double nextStateMoment = momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), (-1) * lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());
                // create a magnetic moment value that is a mixture of the current and next orientations of the spin
                double magneticMoment = perc * nextStateMoment + (1.0 - perc) * prevStateMoment;

                // save the lattice in case we need to fall back to this step
                singleSpin[] tempLattice = Lattice.copyLattice(lattice);
                // add this mixed value to other sites
                addMixedSpinField(flipSpin, magneticMoment);

                if (numOfStepsLeft > 0) {
                    // this is the basic solver, but the flipped spin is kept at its initial orientation.
                    // this way initially its magnetic moment is consistent with its local field, and should
                    // not change.
                    // as for the others, they are updated to conform to the local fields that now include a
                    // contribution from a mixture at flipSpin.
                    // one should not assume that the self-consistent calculation is solved after this call
                    updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, flipSpin);
                    // the following checks if the current state is valid, with the flipped spin at its
                    // initial state
                    if (magneticMomentConvergence(flipSpin) > tol) {
                        // if this step did not converge, backtrack to the previous state and try half of the previous step
                        maxIter = (int) (1.1 * maxIter);
                        numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                        lattice = tempLattice;
                        perc = prevPerc;
                    } else {
                        removeSpecificSpinConstibution(flipSpin, magneticMoment);
                    }

                } else {    // last step
                    lattice[flipSpin].flipSpin();
                    lattice[flipSpin].setSpinSize(magneticMoment);
                    updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, false);

                    if (magneticMomentConvergence() > tol) {
                        // backtrack
                        lattice[flipSpin].flipSpin();   // flip back
                        maxIter = (int) (1.1 * maxIter);
                        numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                        lattice = tempLattice;
                        perc = prevPerc;
                    }
                }

                prevPerc = perc;
                if (numOfStepsLeft > 0) perc = perc + (1.0 - perc) / numOfStepsLeft;
                numOfStepsLeft--;

                // too many iterations, so this is not going anywhere (unless this was the last step):
                if ((i > 150 || numOfStepsLeft > 100) && numOfStepsLeft >= 0) {
                    numOfStepsLeft = -1;  // this ends the loop
                    addMixedSpinField(flipSpin, lattice[flipSpin].getSpinSize());   // this make the fields match the moments
                    lattice[flipSpin].flipSpin();

                    ConvergenceException e = new ConvergenceException.Builder("Too many homotopic step taken without convergence. numOfStepsLeft=" + numOfStepsLeft +
                            ", i=" + i, "homotopic")
                            .setIndex(i)
                            .setConvergenceDistance(magneticMomentConvergence())
                            .build();
                    throw e;
                }
            }
        }

        /**
         * Solves the self-consistent calculation by gradually flipping the given spin in a manner that does not make the
         * target function for the magnetic moment ill-defined at any step.
         * @param bfsOrder the order in which we loop through the spins
         * @param maxIter maximum number of iterations over the entire system
         * @param tol tolerance for convergence
         * @param flipSpin the spin that is being flipped
         * @param alpha relaxation parameter for the update in #1 (for values less than 1, we update the magnetic moment
         *              to a weighted average of the previous magnetic moment and the new magnetic moment
         * @throws ConvergenceException
         * @see FieldTable#getValue(double, double, double, double, int[][], double, boolean)
         */
        public void homotopicSolve3(int[] bfsOrder, int maxIter, double tol, int flipSpin, double alpha) throws ConvergenceException {
            // This function smoothly changes the magnetic moment target function of the flipped spin while keeping it injective
            int numOfStepsLeft = 0;
            double perc = 1.0, prevPerc = 0.0;
            int i = 0;

            while (numOfStepsLeft >= 0) {
                i++;
                singleSpin[] tempLattice = Lattice.copyLattice(lattice);
                // the mixing of the previous and next spin orientations is built into the target function, which is called by passing smoothHOmotopy=true below
                updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, flipSpin, false, perc, false, true);

                if (magneticMomentConvergence(flipSpin, perc, true) > tol) {
                    // if this step did not converge, backtrack to the previous state and try half of the previous step
                    numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                    maxIter = (int) (1.1 * maxIter);
                    lattice = tempLattice;
                    perc = prevPerc;
                }
                prevPerc = perc;
                if (numOfStepsLeft > 0) perc = perc + (1.0 - perc) / numOfStepsLeft;
                numOfStepsLeft--;

                if (numOfStepsLeft < 0) { // last step
                    lattice[flipSpin].flipSpin();
                }

                // too many iterations, so this is not going anywhere (unless this was the last step):
                if ((i > 300 || numOfStepsLeft > Math.pow(2, 8)) && numOfStepsLeft >= 0) {
                    lattice[flipSpin].flipSpin();   // this does nothing. it is only so that magneticMomentConvergence that is called
                                                    // for the error information is run correctly.
                    ConvergenceException e = new ConvergenceException.Builder("Too many homotopic step taken without convergence. numOfStepsLeft=" + numOfStepsLeft +
                            ", i=" + i, "homotopic")
                            .setIndex(i)
                            .setConvergenceDistance(magneticMomentConvergence())
                            .build();

                    numOfStepsLeft = -1;  // this ends the loop (just in case the throw is removed)
                    throw e;
                }
            }
        }

    }

}
