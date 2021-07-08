package simulation.mmsolve;

import simulation.montecarlo.FieldTable;
import simulation.montecarlo.returnDoubleAndStatus;
import simulation.montecarlo.singleSpin;

/**
 * The function whose roots will be found using Broyden's Newton's method.
 * This one also receives information about the spin that's currently being flipped
 * so that it treats it accordingly (varies its required magnetic moment gradually).
 * extends fi_xi as it's a function from R^n to R^n.
 */
class funcForHomotopy extends func {
    /** The spin that is being flipped (gradually), and hence its magnetic moment is a mix of up and down. */
    final int flippedSpin;
    /** Percentage of how much of the magnetic moment of the flipped spin is in the previous orientation
     * and how much is in the next orientation. */
    final double perc;

    public funcForHomotopy(final double[][][] intTable0, final FieldTable momentTable0, final singleSpin[] arr0, double extBx0, double extBy0, final int flippedSpin0, final double perc0, final boolean suppressInternalTransFields0){
        super(intTable0, momentTable0, arr0, extBx0, extBy0, suppressInternalTransFields0);
        flippedSpin=flippedSpin0;
        perc=perc0;
    }

    public int getFlippedSpin() {
        return flippedSpin;
    }

    public double getPerc() {
        return perc;
    }

    /**
     * Evaluates the function func at point x[0..n-1]
     * @param x - The vectors of variables at which the function is to be evaluated
     * @return The vector function at point x[0..n-1], \vec f ( \vec x )
     */
    public double[] func(double[] x) throws ConvergenceException {
        double[] retVec=new double[x.length];
        for (int i = 0; i<retVec.length && numManualCalc <20; i++) {
            if (i!=flippedSpin)  retVec[i] = x[i] - g(x, i, intTable, momentTable, arr, extBx, extBy, suppressInternalTransFields);
            else retVec[i] = x[i] - g(x, i, intTable, momentTable, arr, extBx, extBy, flippedSpin);
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
     * Auxiliary function which calculates the field at site i and returns the appropriate magnetic moment for that field for the flipped spin.
     * Overloads {@link func#g(double[], int, double[][][], FieldTable, singleSpin[], double, double, boolean)} and should be
     * used only for getting the target magnetic moment of the spin that is currently being flipped.
     * @param x - The current configuration of magnetic moments
     * @param i - The site i at which to calculate the moment
     * @param int_config_Matrix  - Interaction table
     * @param momentTable - Magnetic moment table
     * @param arr - The full system array. Used for the spin direction
     * @param extBx - External Bx magnetic field
     * @param extBy - External By magnetic field
     * @param flipSpin - the spin that is being flipped. Not used, just for overloading.
     * @return the objective moment, i.e. what the magnetic moment should be, considering the current configuration
     */
    private double g(final double[] x, int i, double[][][] int_config_Matrix, FieldTable momentTable, singleSpin[] arr, double extBx, double extBy, int flipSpin){
        // this function should only be used to calculate the target magnetic moment of the flipped spin
        if (i!=flippedSpin){
            System.err.println("Wrong function called to calculate the magnetic moment. For the flipped spin g should be called with the index of the spin to flip.");
            System.exit(1);
        }

        double[] B = new double[]{extBx, extBy, 0};
        // get the local 3-component magnetic field at site i
        for (int dim=0; dim<B.length; dim++){
            if (!suppressInternalTransFields || dim==2) {
                for (int j = 0; j < int_config_Matrix[dim][i].length; j++) {
                    B[dim] += int_config_Matrix[dim][i][j] * x[j];
                }
            }
        }

        // given the local magnetic field, mix the previous and next magnetic moments with the given percentage

        returnDoubleAndStatus retNext = momentTable.getValue(B[0], B[1], B[2], -1 * arr[i].getSpin(), true);
        if (!retNext.isSuccessful()) numManualCalc++; // to keep track of how many times python has been called

        returnDoubleAndStatus retPrev = momentTable.getValue(B[0], B[1], B[2], arr[i].getSpin(), true);
        if (!retPrev.isSuccessful()) numManualCalc++; // to keep track of how many times python has been called

        return retNext.getValue()*perc + retPrev.getValue()*(1-perc);

    }

}

