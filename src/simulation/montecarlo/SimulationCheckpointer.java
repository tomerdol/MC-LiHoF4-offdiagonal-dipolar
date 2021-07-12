package simulation.montecarlo;

import org.apache.commons.lang3.SerializationUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

/**
 * Saves and loads a checkpoint that can be used to restart the simulation if stopped or crashed
 */
public class SimulationCheckpointer {
    /** Name of the directory where the checkpoint file will be saved */
    private final String folderName;
    /** {@link File} object which accesses the checkpoint */
    private final File fSaveState;

    /**
     * Constructs a new checkpointer object used to save and load simulation checkpoints
     * @param folderName name of the directory where the checkpoint file is saved
     * @param Lx number of unit cells in the x and y directions
     * @param Lz number of unit cells in the z direction
     * @param extBx external magnetic field in the x direction
     * @param suppressInternalTransFields whether to suppress internal transverse fields (exclude offdiagonal dipolar terms)
     * @param seed random number generator seed
     */
    public SimulationCheckpointer(final String folderName, final int Lx, final int Lz, final double extBx, final boolean suppressInternalTransFields, final long seed){
        Main.makeDir(System.getProperty("system") + File.separator + "checkpoints" + File.separator, folderName);
        this.folderName = folderName;
        this.fSaveState = new File(System.getProperty("system") + File.separator + "checkpoints" + File.separator + folderName + File.separator + "save_state_" + Lx + "_" + Lz + "_" + extBx + "_" + suppressInternalTransFields + "_" + seed + ".txt");
    }

    /**
     * Reads the checkpoint file of a simulation
     * @return {@link MultipleTMonteCarloSimulation} object read from the file or {@code null} if reading was unsuccessful
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
                // for any problem reading previous state, just start from the beginning
                successReadFromFile = false;
                System.err.println("Error reading checkpoint file: " + e.toString());
            }
        }
        return simulation;
    }

    /**
     * Writes a checkpoint file of a simulation, which can be used to start
     * the simulation again from this point
     * @param simulation the simulation object whose state is to be saved to the checkpoint file
     */
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
     * @param parallelTemperingOff whether to turn off parallel tempering (never swap temperatures)
     * @param parallelMode run mode: 's' for serial and 'p' for parallel
     * @param spinSize Typical magnetic moment based on the external transverse field
     * @param tol tolerance for convergence of self-consistent calculation
     * @param J_ex Exchange interaction parameter
     * @param seed random number generator seed
     * @param dilution boolean array that indicates which spins exist (true) and which are diluted (false)
     * @param simulation simulation whose compatibility with the other given arguments is to be checked
     * @return a string containing the incompatibilities between checkpoint parameters and current parameters.
     *          if empty, that means it's compatible.
     */
    public static String verifyCheckpointCompatibility(final double[] T, final boolean parallelTemperingOff, final char parallelMode, final double spinSize,
                                                       final double tol, final double J_ex, final long seed, final boolean[] dilution, MonteCarloSimulation simulation){
        String ret="";
        if (!verifyNumOfTemperaturesCompatibility(parallelMode, simulation)) ret+="single vs. multiple temperatures, ";

        if (!verifyJexCompatibility(J_ex, simulation)) ret+="J_ex, ";
        if (!verifySpinSizeCompatibility(spinSize, simulation)) ret+="spinSize, ";
        if (!verifyTolCompatibility(tol, simulation)) ret+="tol, ";
        if (!verifySeedCompatibility(seed, simulation)) ret+="seed, ";
        if (!verifyDilutionCompatibility(dilution, simulation)) ret+="dilution, ";

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

    /**
     * Verify the compatibility of a dilution structure with that of the given simulation.
     * @param dilution boolean array that indicates which spins exist (true) and which are diluted (false)
     * @param simulation simulation object to compare with the given dilution array for incompatibilities
     * @return whether the simulation's dilution is compatible with the given dilution array ({@code true}=compatible)
     */
    public static boolean verifyDilutionCompatibility(final boolean[] dilution, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify dilution compatibility.");
        if (simulation instanceof MultipleTMonteCarloSimulation) {
            for (int i=0;i< ((MultipleTMonteCarloSimulation) simulation).T.length;i++) {
                for (int j = 0; j < dilution.length; j++) {
                    // if a spin does not exists in the dilution array but exists in the lattice
                    if (dilution[j] == false && ((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).getLattice().spinExists(j))
                        return false;
                    // if a spin exists in the dilution array but not in the lattice
                    if (dilution[j] == true && !((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).getLattice().spinExists(j))
                        return false;
                }
            }
            // no incompatibility was found
            return true;
        }
        if (simulation instanceof SingleTMonteCarloSimulation) {
            for (int j = 0; j < dilution.length; j++) {
                // if a spin does not exists in the dilution array but exists in the lattice
                if (dilution[j] == false && ((SingleTMonteCarloSimulation) simulation).getLattice().spinExists(j))
                    return false;
                // if a spin exists in the dilution array but not in the lattice
                if (dilution[j] == true && !((SingleTMonteCarloSimulation) simulation).getLattice().spinExists(j))
                    return false;
            }
            return true;
        }
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    /**
     * Verify the compatibility of a seed with that of the given simulation.
     * @param seed random number generator seed
     * @param simulation simulation object to compare with the given seed for incompatibility
     * @return whether the simulation's seed is compatible with the given seed ({@code true}=compatible)
     */
    public static boolean verifySeedCompatibility(final long seed, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify seed compatibility.");
        if (simulation instanceof MultipleTMonteCarloSimulation)
            return seed==simulation.getSeed();
        if (simulation instanceof SingleTMonteCarloSimulation)
            return seed==simulation.getSeed();
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    /**
     * Verify the compatibility of an exchange parameter with that of the given simulation.
     * @param J_ex Exchange interaction parameter
     * @param simulation simulation object to compare with the given Exchange parameter for incompatibility
     * @return whether the simulation's exchange parameter is compatible with the given J_ex ({@code true}=compatible)
     */
    public static boolean verifyJexCompatibility(final double J_ex, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify Jex compatibility.");
        if (simulation instanceof MultipleTMonteCarloSimulation)
            return J_ex==((MultipleTMonteCarloSimulation)simulation).J_ex;
        if (simulation instanceof SingleTMonteCarloSimulation)
            return J_ex==((SingleTMonteCarloSimulation)simulation).J_ex;
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    /**
     * Verify the compatibility of a given spinSize with that of the given simulation.
     * @param spinSize Typical magnetic moment based on the external transverse field
     * @param simulation simulation object to compare with the given spinSize for incompatibility
     * @return whether the simulation's spinSize is compatible with the given spinSize ({@code true}=compatible)
     */
    public static boolean verifySpinSizeCompatibility(final double spinSize, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify spinSize compatibility.");

        if (simulation instanceof MultipleTMonteCarloSimulation)
            return spinSize==((MultipleTMonteCarloSimulation)simulation).spinSize;
        if (simulation instanceof SingleTMonteCarloSimulation)
            return spinSize==((SingleTMonteCarloSimulation)simulation).spinSize;
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    /**
     * Verify the compatibility of a given tolerance with that of the given simulation.
     * @param tol tolerance for convergence of self-consistent calculation
     * @param simulation simulation object to compare with the given tolerance for incompatibility
     * @return whether the simulation's tolerance is compatible with the given tol ({@code true}=compatible)
     */
    public static boolean verifyTolCompatibility(final double tol, final MonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify tol compatibility.");

        if (simulation instanceof MultipleTMonteCarloSimulation)
            return tol==((MultipleTMonteCarloSimulation)simulation).tol;
        if (simulation instanceof SingleTMonteCarloSimulation)
            return tol==((SingleTMonteCarloSimulation)simulation).tol;
        throw new IllegalArgumentException("simulation arg is neither an instance of SingleTMonteCarloSimulation nor of MultipleTMonteCarloSimulation");
    }

    /**
     * Verify the compatibility of a given temperature array with that of the given simulation.
     * @param T array of temperatures
     * @param simulation simulation object to compare with the given temperature array for incompatibilities
     * @return whether the simulation's temperature array is compatible with the given T[] ({@code true}=compatible)
     */
    public static boolean verifyTemperaturesCompatibility(final double[] T, final MultipleTMonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");

        return Arrays.equals(T, simulation.getT());
    }

    /**
     * Verify the compatibility of a given temperature with that of the given simulation.
     * @param T temperature
     * @param simulation simulation object to compare with the given temperature for incompatibility
     * @return whether the simulation's temperature is compatible with the given T ({@code true}=compatible)
     */
    public static boolean verifySingleTemperatureCompatibility(final double T, final SingleTMonteCarloSimulation simulation){
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify temperatures' compatibility.");

        return simulation.getT()==T;
    }

    /**
     * Verify the compatibility of a given parallel tempering switch with that of the given simulation.
     * @param parallelTemperingOff whether to turn off parallel tempering (never swap temperatures)
     * @param simulation simulation object to compare with the given parallel tempering option for incompatibility
     * @return whether the simulation's parallel tempering is compatible with the given PT option ({@code true}=compatible)
     */
    public static boolean verifyParallelTemperingCompatibility(final boolean parallelTemperingOff, final MultipleTMonteCarloSimulation simulation ) {
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify parallel tempering compatibility.");

        return simulation.isParallelTemperingOff() == parallelTemperingOff;

    }

    /**
     * Verify the compatibility of a given parallel mode with that of the given simulation.
     * @param parallelMode run mode: 't' for single temperature, 's' or 'p' for multiple temperatures.
     * @param simulation simulation object to compare with the given parallel mode option for incompatibility
     * @return whether the simulation is of the correct class, compatible with the given parallel mode ({@code true}=compatible)
     */
    public static boolean verifyNumOfTemperaturesCompatibility(final char parallelMode, final MonteCarloSimulation simulation) {
        if (simulation==null) throw new NullPointerException("Received Null simulation. cannot verify T compatibility.");

        return ((simulation instanceof SingleTMonteCarloSimulation && parallelMode=='t') ||
                (simulation instanceof MultipleTMonteCarloSimulation && parallelMode!='t'));
    }
}
