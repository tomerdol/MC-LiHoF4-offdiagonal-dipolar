package simulation.mmsolve;

/**
 * function from Rn to R
 */
public abstract class f_xi{
    abstract double calc(double x[]) throws ConvergenceException;
}
