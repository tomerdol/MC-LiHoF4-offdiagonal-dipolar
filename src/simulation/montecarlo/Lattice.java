package simulation.montecarlo;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.random.MersenneTwister;
import simulation.mmsolve.ConvergenceException;
import simulation.mmsolve.fi_xi;
import simulation.mmsolve.func;

import java.io.*;
import java.util.Arrays;
import java.util.Properties;

public class Lattice implements Serializable {
    private static final long serialVersionUID = -9119463380760410942L;
    private final int N, Lx, Lz;
    private final double extBx;
    private final boolean suppressInternalTransFields;
    private final double spinSize;
    // after deserialization these must be set:
    transient double[][][] intTable=null;
    transient double[][] exchangeIntTable=null;
    transient FieldTable momentTable=null;
    transient FieldTable energyTable=null;
    transient int[][] nnArray=null;
    private transient ObservableExtractor measure;
    private transient MagneticMomentsSolveIter iterativeSolver;

    private singleSpin[] lattice;

    @CreatesInconsistency("If intTable, exchangeIntTable, energyTable, momentTable or measure are null")
    public Lattice(int Lx, int Lz, double extBx, boolean suppressInternalTransFields, double spinSize, double[][][] intTable, double[][] exchangeIntTable, int[][] nnArray, FieldTable energyTable, FieldTable momentTable, final ObservableExtractor measure){
        this.N=Lx*Lx*Lz*4;
        this.Lx=Lx;
        this.Lz=Lz;
        this.extBx=extBx;
        this.suppressInternalTransFields=suppressInternalTransFields;
        this.spinSize = spinSize;
        this.intTable=intTable;
        this.exchangeIntTable=exchangeIntTable;
        this.energyTable=energyTable;
        this.momentTable=momentTable;
        this.nnArray=nnArray;
        this.measure=measure;
        this.iterativeSolver = new MagneticMomentsSolveIter();
        this.lattice=generateIsingLattice(Lx,Lz,spinSize);
        if (intTable!=null && exchangeIntTable!=null && energyTable!=null && momentTable!=null && measure!=null)
            this.updateAllLocalFields();
    }

    public Lattice(Lattice other){
        this.N=other.N;
        this.Lx=other.Lx;
        this.Lz=other.Lz;
        this.extBx=other.extBx;
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
        this.extBx=other.extBx;
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

    public boolean isSuppressInternalTransFields() {
        return suppressInternalTransFields;
    }

    public singleSpin[] getArray() {
        return copyLattice(this.lattice);
    }

    @CreatesInconsistency
    // randomize the spin configurations
    public void randomizeConfig(MersenneTwister rnd){
        for (int i=0; i<lattice.length;i++){
            if (lattice[i].getSpin()!=0){
                if (rnd.nextBoolean()) lattice[i].setSpin(1);
                else lattice[i].setSpin(-1);
                lattice[i].setSpinSize(spinSize*lattice[i].getSpin());
            }
        }
    }

    @CreatesInconsistency
    public void checkerBoard(){
        for (int i=0; i<lattice.length;i++){
            lattice[i].setSpin(i%2==0 ? 1 : -1);
            lattice[i].setSpinSize(spinSize*lattice[i].getSpin());
        }
    }


    public void flipSpin(int maxIter, double tol, int flipSpin, double alpha, MersenneTwister rnd) throws ConvergenceException {
        singleSpin[] tempLattice = Lattice.copyLattice(lattice);    // save original lattice. This is before the spin is actually flipped.

        boolean success=false;
        int[] methodsToTry = new int[]{1, 8, 2, 18, 3, 5, 15, 4, 13, 7, 16, 17, 10, 14, 6, 11, 12, 9, 19};

        String errorMessage="There was an error flipping spin " + flipSpin + ". The methods tried and the resulting errors are as follows:\n" +
                Arrays.toString(methodsToTry) + "\n";

        int methodIndex;
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
    }

    public void updateAllMagneticMoments(int maxIter, double tol, double alpha){
        iterativeSolver.updateAllMagneticMoments(maxIter, tol, alpha, true);
    }
    private singleSpin[] solveSelfConsistentCalc(int maxIter, double tol, int flipSpin, int method, int[][] nnArray, double alpha, MersenneTwister rnd) throws ConvergenceException {
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
                // method numbers in the range 3<=method<=7 use the given initial guess
                for (int j = 0; j < lattice.length; j++) x[j] = lattice[j].getSpinSize();
            } else if (method >= 9 && method <= 13) {
                // method numbers in the range 8<=method<=12 use a random initial guess
                for (int j = 0; j < x.length; j++) x[j] = spinSize * rnd.nextDouble() * (rnd.nextBoolean() ? 1 : -1);
                method -= 5;
            } else if (method >= 14 && method <= 18) {
                // method numbers in the range 13<=method<=17 use a "good" initial guess
                for (int j = 0; j < x.length; j++) x[j] = spinSize * lattice[j].getSpin();
                method -= 10;
            }

            fi_xi funcToSolve = new func(intTable, momentTable, lattice, extBx, suppressInternalTransFields);

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



    // returns new int array which has the BFS order with respect to the given spin index
    // n is the maximum (arr.length)
    private int[] orderBFS(int n, int root) {
        int[] bfs = new int[n];
        int head = 0, tail = 0; // head and tail of the queue
        boolean[] visited = new boolean[n];

        // queue the root
        bfs[tail++] = root;
        visited[root] = true;

        while (head < tail) {
            // dequeue v
            int v = bfs[head++];

            //get v's neighbors
            for (int neighbor = 0; neighbor < nnArray[v].length; neighbor++) {
                if (!visited[nnArray[v][neighbor]] && tail < n) {
                    bfs[tail++] = nnArray[v][neighbor];
                    visited[nnArray[v][neighbor]] = true;
                }
            }

        }


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

        return bfs;
    }

    // copies the contents of the given iterator to an array of size N
    // used for saving a configuration
    private static singleSpin[] copyLattice(singleSpin[] arr){
        singleSpin[] newLatticeArr = new singleSpin[arr.length];
        for (int i=0;i<arr.length;i++){
            newLatticeArr[i] = new singleSpin(arr[i]);
        }
        return newLatticeArr;
    }

    public Lattice getCopy(){
        return new Lattice(this);
    }


    public void updateAllLocalFields(){
        int i;

        for (i=0;i<lattice.length;i++) {
            int j;
            if (lattice[i].getSpin()!=0){
                double Bz=0, Bx = extBx, By = 0;
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
                // if suppressInternalTransFields==true then Bx==extBx and By==0.
                lattice[i].setLocalBx(Bx);
                lattice[i].setLocalBy(By);
            }

        }


    }

    private double magneticMomentConvergence(int flipSpin, double frac, boolean smoothHomotopy){
        double converged = 0;
        double max = 0;

        for (int i = 0; i < N; i++) {
            if (i!=flipSpin) {
                double shouldBe = momentTable.getValue(lattice[i].getLocalBx(), lattice[i].getLocalBy(), lattice[i].getLocalBz(), lattice[i].getSpin(), lattice[i].getPrevBIndices());
                double is = lattice[i].getSpinSize();

                converged += Math.abs(is - shouldBe);
//                if (Math.abs(is - shouldBe)>max) max=Math.abs(is - shouldBe);
            }else{
//                double shouldBe = (1-frac)*momentTable.getValue(lattice[i].getLocalBx(), lattice[i].getLocalBy(), lattice[i].getLocalBz(), lattice[i].getSpin(), lattice[i].getPrevBIndices())+frac*momentTable.getValue(lattice[i].getLocalBx(), lattice[i].getLocalBy(), lattice[i].getLocalBz(), -1*lattice[i].getSpin(), lattice[i].getPrevBIndices());
                double shouldBe = momentTable.getValue(lattice[i].getLocalBx(), lattice[i].getLocalBy(), lattice[i].getLocalBz(), lattice[i].getSpin(), lattice[i].getPrevBIndices(), frac, smoothHomotopy);
                double is = lattice[i].getSpinSize();

                converged += Math.abs(is - shouldBe);
//                if (Math.abs(is - shouldBe)>max) max=Math.abs(is - shouldBe);
            }

        }

        return Math.abs(converged/lattice.length);
    }

    // standard call
    private double magneticMomentConvergence(int flipSpin, double frac){ return magneticMomentConvergence(flipSpin, frac, false); }

    // calculate the magnetic moment convergence without any frac spins.
    // received array should have all spins in the correct state (flipSpin should have already been flipped)
    public double magneticMomentConvergence(){
        return magneticMomentConvergence(-1, 0);
    }

    // calculate the magnetic moment convergence without any frac spins.
    // flipSpin is ignored (required for homotopicSolve2)
    private double magneticMomentConvergence(int flipSpin){
        return magneticMomentConvergence(flipSpin, 0);
    }

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

    public boolean verifyAllLocalFields(){
        boolean ret=true;
        double tol = 1.0e-10;
        int i;
        //double max = 0;
        //int maxIndex=0;
        for (i=0;i<lattice.length && ret;i++) {
            int j;

            double Bx=extBx, By=0, Bz=0;
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
            // check is Bx,By or Bz are not equal to the values stored in the i-th site up to a tolerance constant
            ret &= (Math.abs(lattice[i].getLocalBx()-Bx)<tol) && (Math.abs(lattice[i].getLocalBy()-By)<tol) && (Math.abs(lattice[i].getLocalBz()-Bz)<tol);


        }

        if (!ret){
            System.err.println("bad field at spin " + i);
        }

        //System.out.println("max difference from real local field is: " + max);
        return ret;
    }

    /**
     * Reads an Ising lattice from file and creates initial random configuration
     * Only required when there is disorder, i.e. dilution<1.0.
     * @param fileNumber - as in "config_%Lx%_%Lz%_%x%_%h%_%fileNumber%.txt
     * @param dilution - dilution of Ising lattice
     * @param Lx - number of unit cells in x direction of the lattice
     * @param Lz  - number of unit cells in z direction of the lattice
     * @param h - standard deviation of the random local field
     * @param rnd - MersenneTwister PRNG. Used to make the random starting configuration
     * @return Received Ising lattice
     * @throws IOException
     */
    static singleSpin[] receiveIsingLattice(int fileNumber, double dilution, int Lx, int Lz, double h, MersenneTwister rnd) throws IOException {
        final double spinSize = CrystalField.getMagneticMoment(0.0, 0.0, 0.05);

        singleSpin[] arr = null;
        try (BufferedReader in = new BufferedReader(new FileReader("data" + File.separator + "configurations" + File.separator + "config_"+Lx+"_"+Lz+"_"+dilution+"_"+h+"_"+fileNumber+".txt"))){;
            String str;
            String[] params;
            int i=0;


            if ((str = in.readLine()) != null)
                Lx=Integer.parseInt(str.split("=")[1]);
            if ((str = in.readLine()) != null)
                Lz=Integer.parseInt(str.split("=")[1]);

            arr = new singleSpin[4*Lx*Lx*Lz];

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
     * Generates a singleSpin array of LiHo_{x}F_4
     * @param Lx - Number of unit cells in x and in y directions
     * @param Lz - Number of unit cells in z direction
     * @return An array of spins that represent a Lx*Lx*Lz lattice of LiHo_{x}F_4
     * with all local fields set to zero and all spin sizes set to the (positive) default spin size (from parameter file).
     */
    public static singleSpin[] generateIsingLattice(int Lx, int Lz, final double spinSize) {
        int i, j, k, l;
        // create the array that will hold the lattice. the array's cells correspond
        // to the unit cells of the LiHo{x}Y{x-1}F4.
        singleSpin[] arr = new singleSpin[4*Lx*Lx*Lz];

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // get location matrix (3D coordinates for each of the 4 atoms)
        Properties params = GetParamValues.getParams();
        double[][] location = new double[4][3];	// 3D coordinate location for each of the 4 atoms

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
                    for (l = 0; l < 4; l++)
                    {
                        arr[i*Lx*Lz*4+j*Lz*4+k*4+l]=new singleSpin(1,i*Lx*Lz*4+j*Lz*4+k*4+l, spinSize);
                    }
                }
            }
        }

        return arr;
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

    /*
    All methods here should be considered as ruining the singleSpin[] lattice in case an exception is thrown,
    hence a copy should be made before calling them.
     */
    private class MagneticMomentsSolveIter {

        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, boolean retryWithLowerAlpha) {
            updateAllMagneticMoments(sweepOrder, maxIter, tol, alpha, -1, false, 0, retryWithLowerAlpha, false);
        }

        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, int flipSpin) {
            updateAllMagneticMoments(sweepOrder, maxIter, tol, alpha, flipSpin, false, 0, false, false);
        }

        public void updateAllMagneticMoments(int maxIter, double tol, double alpha) {
            updateAllMagneticMoments(null, maxIter, tol, alpha, -1, false, 0, false, false);
        }

        public void updateAllMagneticMoments(int maxIter, double tol, double alpha, int flipSpin, double frac) {
            updateAllMagneticMoments(null, maxIter, tol, alpha, flipSpin, false, frac, false, false);
        }

        public void updateAllMagneticMoments(int maxIter, double tol, double alpha, boolean retryWithLowerAlpha) {
            updateAllMagneticMoments(null, maxIter, tol, alpha, -1, false, 0, retryWithLowerAlpha, false);
        }

        public void updateAllMagneticMoments(int[] sweepOrder, int maxIter, double tol, double alpha, int flipSpin, boolean print, double frac) {
            updateAllMagneticMoments(sweepOrder, maxIter, tol, alpha, flipSpin, print, frac, false, false);
        }

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
                        //if (iter==0 || Math.abs(prevSpinSize-newSpinSize)>1.0e-7) { // some significant change was actually made
                        lattice[index].setSpinSize(prevSpinSize * (1 - alpha) + alpha * newSpinSize);
                        updateFieldsAfterSpinSizeChange(index, prevSpinSize);
                        sum += Math.abs(prevSpinSize - lattice[index].getSpinSize());
                        //}
                    } else {
                        prevSpinSize = lattice[index].getSpinSize();
//                        newSpinSize = (1 - frac) * momentTable.getValue(lattice[index].getLocalBx(), lattice[index].getLocalBy(), lattice[index].getLocalBz(), lattice[index].getSpin(), lattice[i].getPrevBIndices()) + frac * momentTable.getValue(lattice[index].getLocalBx(), lattice[index].getLocalBy(), lattice[index].getLocalBz(), -1 * lattice[index].getSpin(), lattice[i].getPrevBIndices());
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
                //for (int i = 0; i < maxIter - iter; i++) System.out.print(tol+" ");
                System.out.println();
            }
            if (!converged && retryWithLowerAlpha) {
                //System.out.println(sum/lattice.length);
                //System.out.println("initiating underrelaxed gauss-seidel");
                //System.out.println();
                updateAllMagneticMoments(sweepOrder, 2 * maxIter, tol, 0.5 * alpha, flipSpin, print, frac, false, smoothHomotopy);
            } else {
                //System.out.print(iter + "/" + maxIter + "," + momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), lattice[flipSpin].getSpin()) + "," + momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), -1 * lattice[flipSpin].getSpin()) + ",");
            }
        }


        public void removeSpecificSpinConstibution(int spinToRemove, double magneticMomentToRemove) {
            int j;
            if (lattice[spinToRemove].getSpin() != 0) {
                // change local fields (longitudinal & transverse) at all other spins by
                // removing the contribution of i there
                for (j = 0; j < lattice.length; j++) {
                    if (lattice[j].getSpin() != 0) {
                        // remove interaction with prevSpinSize and add interaction with current spinSize
                        if (!suppressInternalTransFields) {
                            lattice[j].setLocalBx(lattice[j].getLocalBx() - magneticMomentToRemove * intTable[0][lattice[spinToRemove].getN()][lattice[j].getN()]);
                            lattice[j].setLocalBy(lattice[j].getLocalBy() - magneticMomentToRemove * intTable[1][lattice[spinToRemove].getN()][lattice[j].getN()]);
                        }
                        lattice[j].setLocalBz(lattice[j].getLocalBz() - magneticMomentToRemove * intTable[2][lattice[spinToRemove].getN()][lattice[j].getN()]);

                    }
                }
            }
        }

        public void addMixedSpinField(int flipSpin, double magneticMoment) {
            int i;

            for (i = 0; i < lattice.length; i++) {
                if (lattice[i].getSpin() != 0) {
                    lattice[i].setLocalBz(lattice[i].getLocalBz() + magneticMoment * intTable[2][lattice[i].getN()][lattice[flipSpin].getN()]);
                    if (!suppressInternalTransFields) {
                        lattice[i].setLocalBx(lattice[i].getLocalBx() + magneticMoment * intTable[0][lattice[i].getN()][lattice[flipSpin].getN()]);
                        lattice[i].setLocalBy(lattice[i].getLocalBy() + magneticMoment * intTable[1][lattice[i].getN()][lattice[flipSpin].getN()]);
                    }
                }
            }
        }

        public void homotopicSolve(Lattice latticeObj, int maxIter, double tol, int flipSpin, int[][] nnArray, double alpha) throws ConvergenceException {
            int[] bfsOrder = null;
            if (nnArray != null) {
                bfsOrder = orderBFS(latticeObj.getN(), flipSpin);
            }
            homotopicSolve(bfsOrder, maxIter, tol, flipSpin, alpha);
        }

        public void homotopicSolve(int[] bfsOrder, int maxIter, double tol, int flipSpin, double alpha) throws ConvergenceException {
            int numOfStepsLeft = 0;
            double perc = 1.0, prevPerc = 0.0;
            int i = 0;

            while (numOfStepsLeft >= 0) {
                i++;
                //System.out.print(" homotopy step: " + i++ + "(" + perc+"%) ("+numOfStepsLeft+")");

                singleSpin[] tempLattice = Lattice.copyLattice(lattice);
                removeSpecificSpinConstibution(flipSpin, lattice[flipSpin].getSpinSize());
                double magneticMoment = (1 - perc) * momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(),
                        lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices()) + perc * momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), -1 * lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());
                addMixedSpinField(flipSpin, magneticMoment);
                lattice[flipSpin].setSpinSize(magneticMoment);
                updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, flipSpin, false, perc);

                //System.out.print(" - "+(magneticMomentConvergence(lattice, momentTable, extBx, flipSpin,perc)<tol ? "success ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin,perc)+"), " : "failure ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin,perc)+"), "));
                if (magneticMomentConvergence(flipSpin, perc) > tol) {
                    //System.err.print("backtracking: " + perc);
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

        public void homotopicSolve2(int maxIter, double tol, int flipSpin, int[][] nnArray, double alpha) throws ConvergenceException {
            int[] bfsOrder = null;
            if (nnArray != null) {
                bfsOrder = orderBFS(lattice.length, flipSpin);
            }
            homotopicSolve2(bfsOrder, maxIter, tol, flipSpin, alpha);
        }

        public void homotopicSolve2(int[] bfsOrder, int maxIter, double tol, int flipSpin, double alpha) throws ConvergenceException {

            removeSpecificSpinConstibution(flipSpin, lattice[flipSpin].getSpinSize());

            int numOfStepsLeft = 0;
            double perc = 1.0, prevPerc = 0.0;
            int i = 0;

            while (numOfStepsLeft >= 0) {
                i++;
                //System.out.print(" homotopy step: " + i++ + "(" + perc+"%) ");

                double prevStateMoment = momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());
                double nextStateMoment = momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), (-1) * lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());
                double magneticMoment = perc * nextStateMoment + (1.0 - perc) * prevStateMoment;

                singleSpin[] tempLattice = Lattice.copyLattice(lattice);

                //updateAllLocalFields(lattice, intTable[2], flipSpin);
                //updateAllLocalTransFields(lattice, intTable, extBx, suppressInternalTransFields, flipSpin);

                addMixedSpinField(flipSpin, magneticMoment);

                if (numOfStepsLeft > 0) {
                    updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, flipSpin);
                    //System.out.print(" - "+(magneticMomentConvergence(lattice, momentTable, extBx, flipSpin)<tol ? "success ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin)+"), " : "failure ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin)+"), "));
                    if (magneticMomentConvergence(flipSpin) > tol) {
                        // backtrack
                        //maxIter *= 2;
                        maxIter = (int) (1.1 * maxIter);
                        numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                        lattice = tempLattice;
                        //magneticMoment = lattice[flipSpin].getSpinSize();
                        perc = prevPerc;
                        //System.err.print("backtracking: " + perc);

                    } else {
                        removeSpecificSpinConstibution(flipSpin, magneticMoment);
                    }

                } else {    // last step
                    lattice[flipSpin].flipSpin();
                    lattice[flipSpin].setSpinSize(magneticMoment);
                    updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, false);

                    //System.out.print(" (last step) - "+(MonteCarloMetropolis.updateFieldsAndCheckMagneticMomentConvergence(lattice, intTable, momentTable, extBx, suppressInternalTransFields)<tol ? "success ("+MonteCarloMetropolis.updateFieldsAndCheckMagneticMomentConvergence(lattice, intTable, momentTable, extBx, suppressInternalTransFields)+"), " : "failure ("+MonteCarloMetropolis.updateFieldsAndCheckMagneticMomentConvergence(lattice, intTable, momentTable, extBx, suppressInternalTransFields)+"), "));

                    if (magneticMomentConvergence() > tol) {
                        // backtrack
                        lattice[flipSpin].flipSpin();   // flip back

                        //maxIter *= 2;
                        maxIter = (int) (1.1 * maxIter);
                        numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                        lattice = tempLattice;
                        //magneticMoment = prevMagneticMoment;
                        perc = prevPerc;
                        //System.err.print("backtracking: " + perc);

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


        public void homotopicSolve3(int[] bfsOrder, int maxIter, double tol, int flipSpin, double alpha) throws ConvergenceException {
            // This function smoothly changes the magnetic moment target function of the flipped spin while keeping it injective
            int numOfStepsLeft = 0;
            double perc = 1.0, prevPerc = 0.0;
            int i = 0;
            double BzRange=5.0; // the smooth change of the target function goes from -BzRange to +BzRange

            while (numOfStepsLeft >= 0) {
                i++;
                //System.out.print(" homotopy step: " + i++ + "(" + perc+"%) ("+numOfStepsLeft+")");

                singleSpin[] tempLattice = Lattice.copyLattice(lattice);
//                removeSpecificSpinConstibution(flipSpin, lattice[flipSpin].getSpinSize());
//                double magneticMoment = (1 - perc) * momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(),
//                        lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices()) + perc * momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), -1 * lattice[flipSpin].getSpin(), lattice[flipSpin].getPrevBIndices());

//                addMixedSpinField(flipSpin, magneticMoment);
//                lattice[flipSpin].setSpinSize(magneticMoment);
                updateAllMagneticMoments(bfsOrder, maxIter, tol, alpha, flipSpin, false, perc, false, true);

                //System.out.print(" - "+(magneticMomentConvergence(lattice, momentTable, extBx, flipSpin,perc)<tol ? "success ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin,perc)+"), " : "failure ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin,perc)+"), "));
                if (magneticMomentConvergence(flipSpin, perc, true) > tol) {
                    //System.err.print("backtracking: " + perc);
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
