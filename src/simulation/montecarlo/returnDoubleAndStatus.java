package simulation.montecarlo;

/**
 * Class used to return a value from a {@link FieldTable} object in addition to a
 * boolean used to indicate whether a manual calculation (exact diagonalization) was performed.
 */
public final class returnDoubleAndStatus{
    /** whether an interpolation was successfully performed or
     * an exact diagonalization had to be performed (happens when
     * the applied field is out of the bounds of the table) */
    private final boolean successful;
    /** Value returned from {@code FieldTable} */
    private final double value;

    /**
     * Constructs a return object
     * @param value0 value obtained from {@link FieldTable}
     * @param status0 indicates whether the value was obtained by interpolation ({@code true})
     *                or by exact diagonalization ({@code false}), which is much slower.
     */
    public returnDoubleAndStatus(final double value0, final boolean status0){
        successful=status0;
        value=value0;
    }

    /**
     * Whether a successful interpolation was performed
     * @return indicator of successful interpolation
     */
    public boolean isSuccessful() {
        return successful;
    }

    /**
     * Gets the value obtained from the {@link FieldTable}
     * @return value (energy or magnetic moment) obtained from the table
     * for the given applied field.
     */
    public double getValue() {
        return value;
    }
}
