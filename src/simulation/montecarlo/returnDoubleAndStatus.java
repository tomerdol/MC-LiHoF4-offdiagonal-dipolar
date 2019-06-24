package simulation.montecarlo;

/**
 * Class used to return a value from this table in addition to a boolean used to indicate whether a manual calculation was performed.
 */
public final class returnDoubleAndStatus{
    private final boolean successful;
    private final double value;
    public returnDoubleAndStatus(final double value0, final boolean status0){
        successful=status0;
        value=value0;
    }

    public boolean isSuccessful() {
        return successful;
    }

    public double getValue() {
        return value;
    }
}
