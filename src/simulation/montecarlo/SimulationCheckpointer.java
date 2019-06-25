package simulation.montecarlo;

import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.SerializationUtils;

import java.io.File;
import java.io.FileInputStream;
import java.util.Arrays;

public class SimulationCheckpointer {
    private final String folderName;
    private final File fSaveState;
    private boolean successReadFromFile;


    public SimulationCheckpointer(final String folderName, final int Lx, final int Lz, final double extBx, final boolean suppressInternalTransFields){
        Main.makeDir("states" + File.separator, folderName);
        this.folderName = folderName;
        this.fSaveState = new File("states" + File.separator + folderName + File.separator + "save_state_" + Lx + "_" + Lz + "_" + extBx + "_" + suppressInternalTransFields + ".txt");
    }


    /**
     * Reads the checkpoint file of a simulation
     * @return monte carlo simulation object read from the file or null if reading was unsuccessful
     */
    public MonteCarloSimulation readCheckpoint(){
        MonteCarloSimulation simulation = null;
        boolean successReadFromFile = false;

        if (fSaveState.exists()){
            try (FileInputStream fis = new FileInputStream(fSaveState)) {
                // Open FileInputStream to the file
                simulation = SerializationUtils.deserialize(fis);
                successReadFromFile=true;
            }catch(Exception e) {
                // for any problem reading previous state, just continue from start
                successReadFromFile = false;
            }
        }
        return simulation;
    }

    /**
     * Verify that a given simulation is compatible with the given parameters.
     * Should be used to verify compatibility of command line arguments with simulation retrieved from checkpoint.
     * @param T Array of temperatures. for single temperature simulation a size 1 array should be given.
     * @param parallelTemperingOff
     * @param parallelMode
     * @param simulation
     * @return string containing the incompatibilities between checkpoint parameters and current parameters.
     */
    public static String verifyCheckpointCompatibility(final double[] T, final boolean parallelTemperingOff, final char parallelMode, MonteCarloSimulation simulation){
        String ret="";
        if (!verifyNumOfTemperaturesCompatibility(parallelMode, simulation)) ret+="single vs. multiple temperatures, ";

        if (!ret.isEmpty()){
            if (parallelMode=='t'){
                if (!verifySingleTemperatureCompatibility(T[0], (SingleTMonteCarloSimulation) simulation))
                    ret += "given vs. read temperature, ";
            }
            else {
                if (!verifyParallelTemperingCompatibility(parallelTemperingOff, (MultipleTMonteCarloSimulation) simulation))
                    ret += "parallel tempering, ";
                if (!verifyTemperaturesCompatibility(T, (MultipleTMonteCarloSimulation) simulation))
                    ret += "temperature schedule, ";
            }
        }

        return ret;

    }

    public static boolean verifyTemperaturesCompatibility(final double[] T, final MultipleTMonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");

        return Arrays.equals(T, simulation.getT());
    }


    public static boolean verifySingleTemperatureCompatibility(final double T, final SingleTMonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");

        return simulation.getT()==T;
    }


    public static boolean verifyParallelTemperingCompatibility(final boolean parallelTemperingOff, final MultipleTMonteCarloSimulation simulation ) {
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify parallel tempering compatibility.");

        return simulation.isParallelTempetingOff() == parallelTemperingOff;

    }

    public static boolean verifyNumOfTemperaturesCompatibility(final char parallelMode, final MonteCarloSimulation simulation) {
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify T compatibility.");

        return ((simulation instanceof SingleTMonteCarloSimulation && parallelMode!='t') ||
                (simulation instanceof MultipleTMonteCarloSimulation && parallelMode=='t'));
    }
}
