package simulation.mmsolve;

import simulation.montecarlo.*;

/**
 * The function whose roots will  be found using Broyden's or Newton's method.
 * extends fi_xi as it's a function from R^n to R^n.
 */
public class func extends fi_xi {
    /** The interactions table for the combined dipolar and exchange interactions between each pair of spins */
    final double[][][] intTable;
    /** The table which gives the magnetic moment of an ion for a given local 3-component applied magnetic field */
    final FieldTable momentTable;
    /** The array of all spins in the system for which a self-consistent solution is to be found */
    final singleSpin[] arr;
    /** External magnetic field B_x */
    final double extBx;
    /** External magnetic field B_y */
    final double extBy;
    /** Number of times an exact diagonalization is performed during the current self-consistent calculation.
     * This happens when we encounter a magnetic field outside the bounds of the given {@code momentTable} */
    int numManualCalc;
    /** Whether internal transverse fields are suppressed or not (offdiagonal dipolar terms included or excluded) */
    final boolean suppressInternalTransFields;


    public func(final double[][][] intTable0, final FieldTable momentTable0, final singleSpin[] arr0, double extBx0, double extBy0, boolean suppressInternalTransFields0){
        intTable = intTable0;
        momentTable=momentTable0;
        arr=arr0;
        extBx=extBx0;
        extBy=extBy0;
        numManualCalc = 0;  // we start with zero times that a "manual" exact diagonalization occurred
        suppressInternalTransFields=suppressInternalTransFields0;
    }

    public int getNumManualCalc() {
        return numManualCalc;
    }

    public double[][][] getIntTable() {
        return intTable;
    }

    public FieldTable getMomentTable() {
        return momentTable;
    }

    public singleSpin[] getArr() {
        return arr;
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
     * Evaluates the function func at point x[0..n-1].
     * @param x - The vectors of variables at which the function is to be evaluated
     * @return The vector function at point x[0..n-1], \vec f ( \vec x )
     */
    public double[] func(double[] x) throws ConvergenceException {
        double[] retVec=new double[x.length];
        for (int i = 0; i<retVec.length && numManualCalc < 20; i++) {
            // the function whose zeros we want is the difference between the given
            // magnetic moment and the magnetic moment dictated by the applied magnetic field.
            retVec[i] = x[i] - g(x, i, intTable, momentTable, arr, extBx, extBy, suppressInternalTransFields);
        }

        if (numManualCalc >= 20) {
            ConvergenceException e = new ConvergenceException.Builder("the function surpassed the permitted threshold (20) for manual calculations. ",
                    "Function evaluation")
                    .setNumManualCalc(numManualCalc)
                    .build();
            throw e;
        }

        return retVec;
    }

    /**
     * Auxiliary function which calculates the field at site i and returns the appropriate magnetic moment for that field
     * @param x - The current configuration of magnetic moments
     * @param i - The site i at which to calculate the moment
     * @param int_config_Matrix  - Interaction table
     * @param momentTable - Magnetic moment table
     * @param arr - The full system array. Used for the spin direction
     * @param extBx - External Bx magnetic field
     * @param extBy - External By magnetic field
     * @return the objective moment, i.e. what the magnetic moment should be, considering the current configuration
     */
    public double g(final double[] x, int i, double[][][] int_config_Matrix, FieldTable momentTable, singleSpin[] arr, double extBx, double extBy, boolean suppressInternalTransFields){
        // local fields start from the external values which are the same everywhere
        double[] B = new double[]{extBx, extBy, 0};

        // get the local 3-component magnetic field, {Bx,By,Bz} at site i
        for (int dim=0; dim<B.length; dim++){
            if (!suppressInternalTransFields || dim==2) {
                for (int j = 0; j < int_config_Matrix[dim][i].length; j++) {
                    B[dim] += int_config_Matrix[dim][i][j] * x[j];
                }
            }
        }

        // from the local field, get the required magnetic moment
        returnDoubleAndStatus ret = momentTable.getValue(B[0], B[1], B[2], arr[i].getSpin(), true);


        if (!ret.isSuccessful()) numManualCalc++; // to keep track of how many times exact diagonalization was performed
        return ret.getValue();
    }


}

