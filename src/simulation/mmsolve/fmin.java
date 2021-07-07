package simulation.mmsolve;

/**
 * Returns 0.5*func(dot)func.
 * Also stores value of func in fvec.
 * Can be used as a minimization target to find the roots of a function from Rn to Rn.
 */
public class fmin extends f_xi {
    /** Function object from Rn to Rn */
    private final fi_xi vecFunc;
    /** Stores the last evaluation of vecFunc */
    private double[] fvec;

    /**
     * Constructs an fmin object that is (half) the sum of squares of the elements of a given
     * function from Rn to Rn. Can be used as a minimization target to find the roots of a function from Rn to Rn.
     * @param func0 - the function whose element should be squared and summed
     */
    fmin(fi_xi func0){
        vecFunc=func0;
    }

    /**
     * Evaluates fmin at the point x[0..n-1] and stores the value of vecFunc in fvec.
     * @param x - the point at which to evalue fmin
     * @return the value of fmin at the given point
     * @throws ConvergenceException when function evaluation exceeds the number of allowed exact diagonalizations
     */
    public double calc(double[] x) throws ConvergenceException {
        double sum=0;
        fvec=vecFunc.func(x);
        for (int i=0;i<x.length;i++)    sum+= fvec[i]*fvec[i];
        return 0.5*sum;
    }

    /**
     * @return the value of vecFunc from the last evaluation of fmin
     */
    double[] getFvec() {
        return fvec;
    }
}

