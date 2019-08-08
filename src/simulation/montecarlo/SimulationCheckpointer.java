package simulation.montecarlo;

import org.apache.commons.lang3.SerializationUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
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
    public MultipleTMonteCarloSimulation readCheckpoint(){
        MultipleTMonteCarloSimulation simulation = null;
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

    public void writeCheckpoint(MultipleTMonteCarloSimulation simulation){
        try (FileOutputStream fos = new FileOutputStream(fSaveState, false)){
            SerializationUtils.serialize(simulation, fos);
        }catch (IOException e){
            System.err.println("Error writing checkpoint file. simulation is not being saved!");
            e.printStackTrace();
        }
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
    public static String verifyCheckpointCompatibility(final double[] T, final boolean parallelTemperingOff, final char parallelMode, final double spinSize,
                                                       final double tol, final double J_ex, MonteCarloSimulation simulation){
        String ret="";
        if (!verifyNumOfTemperaturesCompatibility(parallelMode, simulation)) ret+="single vs. multiple temperatures, ";

        if (!verifyJexCompatibility(J_ex, simulation)) ret+="J_ex, ";
        if (!verifySpinSizeCompatibility(spinSize, simulation)) ret+="spinSize, ";
        if (!verifyTolCompatibility(tol, simulation)) ret+="tol, ";

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

    public static boolean verifyJexCompatibility(final double J_ex, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");
        if (simulation instanceof MultipleTMonteCarloSimulation)
            return J_ex==((MultipleTMonteCarloSimulation)simulation).J_ex;
        if (simulation instanceof SingleTMonteCarloSimulation)
            return J_ex==((SingleTMonteCarloSimulation)simulation).J_ex;
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    public static boolean verifySpinSizeCompatibility(final double spinSize, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");

        if (simulation instanceof MultipleTMonteCarloSimulation)
            return spinSize==((MultipleTMonteCarloSimulation)simulation).spinSize;
        if (simulation instanceof SingleTMonteCarloSimulation)
            return spinSize==((SingleTMonteCarloSimulation)simulation).spinSize;
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    public static boolean verifyTolCompatibility(final double tol, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");

        if (simulation instanceof MultipleTMonteCarloSimulation)
            return tol==((MultipleTMonteCarloSimulation)simulation).tol;
        if (simulation instanceof SingleTMonteCarloSimulation)
            return tol==((SingleTMonteCarloSimulation)simulation).tol;
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
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

        return ((simulation instanceof SingleTMonteCarloSimulation && parallelMode=='t') ||
                (simulation instanceof MultipleTMonteCarloSimulation && parallelMode!='t'));
    }
}
