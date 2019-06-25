package simulation.mmsolve;

import simulation.montecarlo.*;

/**
 * The function whose roots will  be found using Broyden's or Newton's method.
 * extends fi_xi as it's a function from R^n to R^n.
 */
public class func extends fi_xi {

    final double[][][] intTable;
    final FieldTable momentTable;
    final singleSpin[] arr;
    final double extBx;
    int numManualCalc;
    final boolean suppressInternalTransFields;

    public func(final double[][][] intTable0, final FieldTable momentTable0, final singleSpin[] arr0, double extBx0, boolean suppressInternalTransFields0){
        intTable = intTable0;
        momentTable=momentTable0;
        arr=arr0;
        extBx=extBx0;
        numManualCalc =0;
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

    public boolean isSuppressInternalTransFields() {
        return suppressInternalTransFields;
    }

    /**
     * Evaluates the function func at point x[1..n]
     * @param x - The vectors of variables at which the function is to be evaluated
     * @return The vector function at point x[1..n], \vec f ( \vec x )
     */
    public double[] func(double[] x) throws ConvergenceException {
        double[] retVec=new double[x.length];
        for (int i = 0; i<retVec.length && numManualCalc <20; i++) {
            retVec[i] = x[i] - g(x, i, intTable, momentTable, arr, extBx, suppressInternalTransFields);
        }

        if (numManualCalc >=20) {
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
     * @return the objective moment, i.e. what the magnetic moment should be, considering the current configuration
     */
    public double g(final double[] x, int i, double[][][] int_config_Matrix, FieldTable momentTable, singleSpin[] arr, double extBx, boolean suppressInternalTransFields){
        double[] B = new double[]{extBx, 0, 0};

        for (int dim=0; dim<B.length; dim++){
            if (!suppressInternalTransFields || dim==2) {
                for (int j = 0; j < int_config_Matrix[dim][i].length; j++) {
                    B[dim] += int_config_Matrix[dim][i][j] * x[j];
                }
            }
        }
        returnDoubleAndStatus ret = momentTable.getValue(B[0], B[1], B[2], arr[i].getSpin(), true);


        if (!ret.isSuccessful()) numManualCalc++; // to keep track of how many times python has been called
        return ret.getValue();
    }


}

