package simulation.mmsolve;

import org.apache.commons.math3.random.MersenneTwister;

/**
 * Exception for failure to converge when performing the self consistent calculation after a spin-flip.
 * @see simulation.montecarlo.Lattice#solveSelfConsistentCalc(int, double, int, int, int[][], double, MersenneTwister)
 */
public class ConvergenceException extends Exception {
    private int numManualCalc=0, index=-1, flippedSpin;
    String errorLocation;
    private double convergenceDistance=-1;

    /**
     * Constructs a new ConvergenceException. Used by flipSpin
     * @param errorMessage an error message indicating which methods were tried and failed.
     * @param flippedSpin0 ordinal number of the spin attempted to be flipped.
     */
    public ConvergenceException(String errorMessage, int flippedSpin0){
        super(errorMessage);
        this.flippedSpin=flippedSpin0;
    }

    // builder pattern

    /**
     * Constructor to be used by the <i>Builder</i>
     * @param builder - builder with all required and optional configuration options that have been specified by calling methods on this builder.
     */
    private ConvergenceException(Builder builder){
        super(builder.errorMessage);
        if (builder.cause!=null)
            super.initCause(builder.cause);
        this.numManualCalc=builder.numManualCalc;
        this.convergenceDistance=builder.convergenceDistance;
        this.errorLocation=builder.errorLocation;
        this.flippedSpin=builder.flippedSpin;
        this.index=builder.index;

    }

    /**
     * A <i>builder</i> class for creating instances of {@code ConvergenceException}.
     */
    public static class Builder{
        //required
        private final String errorMessage, errorLocation;
        // optional
        private int index;
        private int numManualCalc=0;
        private double convergenceDistance=-1;
        private Throwable cause=null;
        private int flippedSpin;

        /**
         * Creates a new <i>builder</i> object with the mandatory input fields.
         * @param errorMessage - an error message indicating which methods were tried and failed.
         * @param errorLocation - where (in which method and in which stage) the convergence error occurred.
         */
        public Builder(String errorMessage, String errorLocation){
            this.errorLocation=errorLocation;
            this.errorMessage=errorMessage;
        }

        /**
         * Setss the iteration index at the time of the failed convergence.
         * @param index - the iteration index
         * @return a reference to this Builder
         */
        public Builder setIndex(int index){
            this.index=index;
            return this;
        }

        /**
         * Setss number of times a manual diagonalization (as opposed to an interpolation) was performed until the convergence error.
         * @param numManualCalc - number of manual diagonalization calculations
         * @return a reference to this Builder
         */
        public Builder setNumManualCalc(int numManualCalc){
            this.numManualCalc=numManualCalc;
            return this;
        }

        /**
         * Setss the distance (sum of absolute value of difference between the magnetic moment value and its required value, divided by the number of spins) from convergence at the time of the error.
         * @param convergenceDistance - sum of absolute value of difference between the magnetic moment value and its required value, divided by the number of spins
         * @return a reference to this Builder
         */
        public Builder setConvergenceDistance(double convergenceDistance){
            this.convergenceDistance=convergenceDistance;
            return this;
        }

        /**
         * Setss the cause of the convergence error.
         * @param cause - cause for convergence failure, given by the Broyden of Newton methods (singular Jacobian or something of that sort)
         * @return a reference to this Builder
         */
        public Builder setCause(Throwable cause){
            this.cause=cause;
            return this;
        }

        /**
         * Setss the spin that was flipped leading to the convergence error.
         * @param flippedSpin - ordinal number of the flipped spin
         * @return a reference to this Builder
         */
        public Builder setFlippedSpin(int flippedSpin){
            this.flippedSpin=flippedSpin;
            return this;
        }

        /**
         * Creates a new {@code ConvergenceException} with all configuration options that have been specified by calling methods on this builder.
         * @return the new {@code ConvergenceException}
         */
        public ConvergenceException build(){
            return new ConvergenceException(this);
        }

    }

    /**
     * Gets the number of times a manual diagonalization (as opposed to an interpolation) was performed until the convergence error.
     * @return number of manual diagonalization calculations
     */
    public int getNumManualCalcn() {
        return numManualCalc;
    }

    /**
     * Sets number of times a manual diagonalization (as opposed to an interpolation) was performed until the convergence error.
     * @param numManualCalc - number of manual diagonalization calculations
     */
    public void setNumManualCalc(int numManualCalc) {
        this.numManualCalc= numManualCalc;
    }

    /**
     * Gets the iteration index at the time of the failed convergence.
     * @return the iteration index
     */
    public int getIndex() {
        return index;
    }

    /**
     * Sets the iteration index at the time of the failed convergence.
     * @param index - the iteration index
     */
    public void setIndex(int index) {
        this.index = index;
    }

    /**
     * Gets the spin that was flipped leading to the convergence error.
     * @return ordinal number of the flipped spin
     */
    public int getFlippedSpin() {
        return flippedSpin;
    }

    /**
     * Sets the spin that was flipped leading to the convergence error.
     * @param flippedSpin - ordinal number of the flipped spin
     */
    public void setFlippedSpin(int flippedSpin) {
        this.flippedSpin = flippedSpin;
    }

    /**
     * Gets where (in which method and in which stage) the convergence error occurred.
     * @return the location of the convergence exception
     */
    public String getErrorLocation() {
        return errorLocation;
    }

    /**
     * Sets where (in which method and in which stage) the convergence error occurred.
     * @param errorLocation - the location of the convergence exception
     */
    public void setErrorLocation(String errorLocation) {
        this.errorLocation = errorLocation;
    }

    /**
     * Gets the distance (sum of absolute value of difference between the magnetic moment value and its required value, divided by the number of spins) from convergence at the time of the error.
     * @return sum of absolute value of difference between the magnetic moment value and its required value, divided by the number of spins
     */
    public double getConvergenceDistance() {
        return convergenceDistance;
    }

    /**
     * Sets the distance (sum of absolute value of difference between the magnetic moment value and its required value, divided by the number of spins) from convergence at the time of the error.
     * @param convergenceDistance - sum of absolute value of difference between the magnetic moment value and its required value, divided by the number of spins
     */
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
                ", numManualCalc=" + numManualCalc +
                ", convergenceDistance=" + convergenceDistance +
                "} " + super.getMessage() + " Caused by: " + causeError;
    }
}