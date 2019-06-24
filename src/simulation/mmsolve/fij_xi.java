package simulation.mmsolve;

/**
 * function that receives a vector and returns a (NxN) matrix
 */
public abstract class fij_xi{
    abstract double[][] func(double x[]) throws ConvergenceException;
}
