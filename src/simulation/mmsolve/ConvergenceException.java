package simulation.mmsolve;

/**
 * Exception for failure to converge when performing the self consistent calculation after a spin-flip.
 * @see simulation.montecarlo.Lattice#solveSelfConsistentCalc(int, double, int, int, int[][], double)
 */
public class ConvergenceException extends Exception {
    private int numManualCalc=0, index=-1, flippedSpin;
    String errorLocation;
    private double convergenceDistance=-1;

    // this is to be used by flip spin.
    // is should receive an error message indicating which methods were tried and failed.
    // the only important field is flippedSpin and errorMessage
    public ConvergenceException(String errorMessage, int flippedSpin0){
        super(errorMessage);
        this.flippedSpin=flippedSpin0;
    }

    public ConvergenceException(Builder builder){
        super(builder.errorMessage);
        if (builder.cause!=null)
            super.initCause(builder.cause);
        this.numManualCalc=builder.numManualCalc;
        this.convergenceDistance=builder.convergenceDistance;
        this.errorLocation=builder.errorLocation;
        this.flippedSpin=builder.flippedSpin;
        this.index=builder.index;

    }

    public static class Builder{
        //required
        private final String errorMessage, errorLocation;
        // optional
        private int index;
        private int numManualCalc=0;
        private double convergenceDistance=-1;
        private Throwable cause=null;
        private int flippedSpin;

        public Builder(String errorMessage, String errorLocation){
            this.errorLocation=errorLocation;
            this.errorMessage=errorMessage;
        }

        public Builder setIndex(int index){
            this.index=index;
            return this;
        }

        public Builder setNumManualCalc(int numManualCalc){
            this.numManualCalc=numManualCalc;
            return this;
        }

        public Builder setConvergenceDistance(double convergenceDistance){
            this.convergenceDistance=convergenceDistance;
            return this;
        }

        public Builder setCause(Throwable cause){
            this.cause=cause;
            return this;
        }

        public Builder setFlippedSpin(int flippedSpin){
            this.flippedSpin=flippedSpin;
            return this;
        }

        public ConvergenceException build(){
            return new ConvergenceException(this);
        }

    }

    public int getNumManualCalcn() {
        return numManualCalc;
    }

    public void setNumManualCalc(int numManualCalc) {
        this.numManualCalc= numManualCalc;
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
                ", numManualCalc=" + numManualCalc +
                ", convergenceDistance=" + convergenceDistance +
                "} " + super.getMessage() + " Caused by: " + causeError;
    }
}