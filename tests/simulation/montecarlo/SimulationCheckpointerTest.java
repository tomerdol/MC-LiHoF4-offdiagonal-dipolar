package simulation.montecarlo;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class SimulationCheckpointerTest {
    SimulationCheckpointer checkpointer = new SimulationCheckpointer("test", 4, 4, 0.0, true);

    @Before
    public void setUp() throws Exception {
        MonteCarloSimulation simulation=new MultipleTMonteCarloSimulation(new double[]{1.5,1.6,1.7}, null, 1, 1, null, false, false, false, true, checkpointer, 1,1,1);
        checkpointer.writeCheckpoint((MultipleTMonteCarloSimulation) simulation);
    }

    @Test
    public void readCheckpoint() {
        MultipleTMonteCarloSimulation readSim = checkpointer.readCheckpoint();
        assertEquals(5,readSim.tol,0.0001);
    }
}