package simulation.mmsolve;

import simulation.montecarlo.FieldTable;
import simulation.montecarlo.returnDoubleAndStatus;
import simulation.montecarlo.singleSpin;

/**
 * The function whose roots will  be found using Broyden's Newton's method.
 * This one also receives information about the spin that's currently being flipped
 * so that it treats it accordingly.
 * extends fi_xi as it's a function from R^n to R^n.
 */
class funcForHomotopy extends func {
    final int flippedSpin;
    final double perc;

    public funcForHomotopy(final double[][][] intTable0, final FieldTable momentTable0, final singleSpin[] arr0, double extBx0, final int flippedSpin0, final double perc0, final boolean suppressInternalTransFields0){
        super(intTable0, momentTable0, arr0, extBx0, suppressInternalTransFields0);
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
     * Evaluates the function func at point x[1..n]
     * @param x - The vectors of variables at which the function is to be evaluated
     * @return The vector function at point x[1..n], \vec f ( \vec x )
     */
    public double[] func(double[] x) throws ConvergenceException {
        double[] retVec=new double[x.length];
        for (int i=0;i<retVec.length && numCalledPython<20;i++) {
            if (i!=flippedSpin)  retVec[i] = x[i] - g(x, i, intTable, momentTable, arr, extBx, suppressInternalTransFields);
            else retVec[i] = x[i] - g(x, i, intTable, momentTable, arr, extBx, flippedSpin);
        }

        if (numCalledPython>=20) {
            ConvergenceException e = new ConvergenceException("the function surpassed the permitted threshold (20) for manual calculations. ");
            e.setErrorLocation("Function evaluation");
            e.setNumCalledPython(numCalledPython);
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
    private double g(final double[] x, int i, double[][][] int_config_Matrix, FieldTable momentTable, singleSpin[] arr, double extBx, int flipSpin){

        if (i!=flippedSpin){
            System.err.println("Wrong function called to calculate the magnetic moment. For the flipped spin g should be called with the index of the spin to flip.");
            System.exit(1);
        }

        double[] B = new double[]{extBx, 0, 0};

        for (int dim=0; dim<B.length; dim++){
            if (!suppressInternalTransFields || dim==2) {
                for (int j = 0; j < int_config_Matrix[dim][i].length; j++) {
                    B[dim] += int_config_Matrix[dim][i][j] * x[j];
                }
            }
        }

        returnDoubleAndStatus retNext = momentTable.getValue(B[0], B[1], B[2], -1 * arr[i].getSpin(), true);
        if (!retNext.isSuccessful()) numCalledPython++; // to keep track of how many times python has been called


        returnDoubleAndStatus retPrev = momentTable.getValue(B[0], B[1], B[2], arr[i].getSpin(), true);
        if (!retPrev.isSuccessful()) numCalledPython++; // to keep track of how many times python has been called
        return retNext.getValue()*perc + retPrev.getValue()*(1-perc);

    }

}

