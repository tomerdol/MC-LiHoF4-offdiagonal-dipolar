package simulation.montecarlo;

import java.io.Closeable;
import java.io.Serializable;

public abstract class MonteCarloSimulation implements Serializable, Runnable, Closeable {
    protected long maxSweeps;
    protected boolean continueFromSave;
    protected long seed;
    protected boolean realTimeEqTest;
    protected boolean checkpoint;

    public void setMaxSweeps(long maxSweeps) {
        this.maxSweeps = maxSweeps;
    }

    public void setContinueFromSave(boolean continueFromSave) {
        this.continueFromSave = continueFromSave;
    }

    public void setCheckpoint(boolean checkpoint) {
        this.checkpoint = checkpoint;
    }

    public void setRealTimeEqTest(boolean realTimeEqTest) {
        this.realTimeEqTest = realTimeEqTest;
    }

    public abstract void run();

    public long getSeed() {
        return seed;
    }
}
