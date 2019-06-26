package simulation.mmsolve;

/**
 * Returns 0.5*func(dot)func . Also stores value of func in fvec.
 */
public class fmin extends f_xi {
    //Returns 0.5*func(dot)func . Also stores value of func in fvec.
    private final fi_xi vecFunc;
    private double[] fvec;

    fmin(fi_xi func0){
        vecFunc=func0;
    }

    public double calc(double[] x) throws ConvergenceException {
        double sum=0;
        fvec=vecFunc.func(x);
        for (int i=0;i<x.length;i++)    sum+= fvec[i]*fvec[i];
        return 0.5*sum;
    }

    double[] getFvec() {
        return fvec;
    }
}

