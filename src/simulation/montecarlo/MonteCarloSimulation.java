package simulation.montecarlo;

import java.io.Closeable;
import java.io.IOException;
import java.io.Serializable;

public abstract class MonteCarloSimulation implements Serializable, Runnable, Closeable {
    protected long maxSweeps;
    protected boolean continueFromSave;
    protected long seed;
    protected boolean realTimeEqTest;
    protected boolean checkpoint;

    public abstract void run();
}
