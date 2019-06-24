package simulation.mmsolve;

import simulation.montecarlo.MonteCarloMetropolis;

/**
 * Exception for failure to converge when performing the self consistent calculation after a spin-flip.
 * @see MonteCarloMetropolis#solveSelfConsistentCalc
 */
public class ConvergenceException extends Exception {
    private int numCalledPython=0, index=-1, flippedSpin;
    String errorLocation;
    private double convergenceDistance=-1;

    // this is to be used by flip spin.
    // is should receive an error message indicating which methods were tried and failed.
    // the only important field is flippedSpin and errorMessage
    public ConvergenceException(String errorMessage, int flippedSpin0){
        super(errorMessage);
        this.flippedSpin=flippedSpin0;
    }

    // this is to be used by the solver wrapper in case the solver does not
    // return any exception but the wrapper's validation fail.
    // in this case there is no cause, index
    public ConvergenceException(String errorMessage, String errorLocation0, double convergenceDistance0){
        super(errorMessage);
        this.errorLocation=errorLocation0;
        this.convergenceDistance=convergenceDistance0;
    }

    // this is for an exception in the function evaluation
    public ConvergenceException(String errorMessage){
        super(errorMessage);
    }

    // ########################################################################
    // the following are to be used by the solving methods themselves.
    // must have fields: index, errorLocation, errorMessage
    // possible fields: numCalledPython, convergenceDistance, cause
    public ConvergenceException(String errorMessage, int index0, String errorLocation0){
        super(errorMessage);
        this.index=index0;
        this.errorLocation=errorLocation0;
    }

    public  ConvergenceException(String errorMessage, int index0, String errorLocation0, double convergenceDistance0){
        this(errorMessage, index0, errorLocation0);
        this.convergenceDistance=convergenceDistance0;
    }

    public  ConvergenceException(String errorMessage, int index0, String errorLocation0, int numCalledPython){
        this(errorMessage, index0, errorLocation0);
        this.numCalledPython=numCalledPython;
    }

    public  ConvergenceException(String errorMessage, int index0, String errorLocation0, int numCalledPython, double convergenceDistance0){
        this(errorMessage, index0, errorLocation0, numCalledPython);
        this.convergenceDistance=convergenceDistance0;
    }

    public ConvergenceException(String errorMessage, Throwable e, int index0, String errorLocation0){
        super(errorMessage, e);
        this.index=index0;
        this.errorLocation=errorLocation0;
    }

    public  ConvergenceException(String errorMessage, Throwable e, int index0, String errorLocation0, int numCalledPython){
        this(errorMessage, e, index0, errorLocation0);
        this.numCalledPython=numCalledPython;
    }

    public  ConvergenceException(String errorMessage, Throwable e, int index0, String errorLocation0, double convergenceDistance0){
        this(errorMessage, e, index0, errorLocation0);
        this.convergenceDistance=convergenceDistance0;
    }
    // ########################################################################


    public int getNumCalledPython() {
        return numCalledPython;
    }

    public void setNumCalledPython(int numCalledPython) {
        this.numCalledPython = numCalledPython;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public int getFlippedSpin() {
        return flippedSpin;
    }

    public void setFlippedSpin(int flippedSpin) {
        this.flippedSpin = flippedSpin;
    }

    public String getErrorLocation() {
        return errorLocation;
    }

    public void setErrorLocation(String errorLocation) {
        this.errorLocation = errorLocation;
    }

    public double getConvergenceDistance() {
        return convergenceDistance;
    }

    public void setConvergenceDistance(double convergenceDistance) {
        this.convergenceDistance = convergenceDistance;
    }

    @Override
    public String toString() {
        String causeError;
        if (this.getCause()!=null) {
            causeError = this.getCause().getMessage();
        }else{
            causeError="-";
        }
        return "ConvergenceException{" +
                "errorLocation=\'" + errorLocation + '\'' +
                ", index=" + index +
                ", numCalledPython=" + numCalledPython +
                ", convergenceDistance=" + convergenceDistance +
                "} " + super.getMessage() + " Caused by: " + causeError;
    }
}