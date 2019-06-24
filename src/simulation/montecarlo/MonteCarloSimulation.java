package simulation.montecarlo;

import java.io.Closeable;
import java.io.Serializable;

public abstract class MonteCarloSimulation implements Serializable, Runnable, Closeable {
    protected long maxSweeps;
    protected boolean continueFromSave;
    protected long seed;
    protected boolean realTimeEqTest;

    public abstract void run();
}