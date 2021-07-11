package simulation.montecarlo;

/**
 * Type of output the simulation gives.
 */
public enum OutputType {
    /** Output all observables each MCS. */
    VERBOSE,
    /** Output observable averages over logarithmically increasing bins */
    BIN,
    /** Output spin-specific observables */
    SPIN
}
