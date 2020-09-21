package simulation.montecarlo;

/**
 * Type of output the simulation gives.
 *  VERBOSE - Output all observables each MCS.
 *  BIN - Output observable averages over logarithmically increasing bins
 *  SPIN - Output spin-specific observables
 */
public enum OutputType {
    VERBOSE,    // Output all observables each MCS
    BIN,        // Output observable averages over logarithmically increasing bins
    SPIN;       // Output spin-specific observables
}
