package simulation.mmsolve;

/**
 * function from Rn to Rn
 */
public abstract class fi_xi{
    abstract double[] func(double x[]) throws ConvergenceException;
}
