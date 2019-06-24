package simulation.mmsolve;

import simulation.montecarlo.MagneticMomentsSolve;
import simulation.montecarlo.f_xi;
import simulation.montecarlo.fi_xi;

/**
 * Returns 0.5*func(dot)func . Also stores value of func in fvec.
 */
public class fmin extends f_xi {
    //Returns 0.5*func(dot)func . Also stores value of func in fvec.
    final fi_xi vecFunc;
    double[] fvec;

    public fmin(fi_xi func0){
        vecFunc=func0;
    }

    public double calc(double[] x) throws ConvergenceException {
        double sum=0;
        fvec=vecFunc.func(x);
        for (int i=0;i<x.length;i++)    sum+= fvec[i]*fvec[i];
        return 0.5*sum;
    }

    public double[] getFvec() {
        return fvec;
    }
}

