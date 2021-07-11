package simulation.montecarlo;

import java.io.Closeable;
import java.io.Serializable;

/**
 * Abstract class for a Monte Carlo simulation (either with multiple temperatures or with a single temperature)
 */
public abstract class MonteCarloSimulation implements Serializable, Runnable, Closeable {
    // for serialization. should not be changed
    private static final long serialVersionUID = 7948934465122188289L;
    /** Total number of MC sweeps to perform*/
    protected long maxSweeps;
    /** whether to continue from a previously saved simulation */
    protected boolean continueFromSave;
    /** Random number generator seed */
    protected long seed;
    /** whether to periodically create a checkpoint from which the simulation can be restarted */
    protected boolean checkpoint;

    /**
     * Adds more sweeps MC to the simulation
     * @param maxSweeps total number of sweeps
     */
    public void addSweeps(long maxSweeps) {
        this.maxSweeps = Math.max(maxSweeps, this.maxSweeps);
    }

    public void setContinueFromSave(boolean continueFromSave) {
        this.continueFromSave = continueFromSave;
    }

    public void setCheckpoint(boolean checkpoint) {
        this.checkpoint = checkpoint;
    }

    /**
     * Prints the current state of the lattice (spin configuration, magnetic moments and local fields)
     */
    public abstract void printSimulationState();

    /**
     * Run the simulation
     */
    public abstract void run();

    public long getSeed() {
        return seed;
    }

    public boolean isCheckpoint() {
        return checkpoint;
    }
}
