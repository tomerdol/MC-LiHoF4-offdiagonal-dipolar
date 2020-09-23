package utilities;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import simulation.montecarlo.MultipleTMonteCarloSimulation;
import simulation.montecarlo.ParseCommandLine;
import simulation.montecarlo.SimulationCheckpointer;

public class fixBackwardsCompatibility {

    public static void main(String[] args){
        int Lx=0;	// lattice x-y size
        int Lz=0;	// lattice z size
        double extBx=-1;   // external Bx
        boolean suppressInternalTransFields=false;
        long seed=0;	// should never happen;
        String folderName="default";
        int newColNum=0;

        // create Options object
        Options options = ParseCommandLine.generateOptions();
        options.addOption(Option.builder("new_col_num")
                .hasArg(true)
                .required(true)
                .type(Number.class)
                .desc("Updated number of columns the simulation prints out.")
                .build());
        // first check if help was called
        for (String s : args) {
            if (s.equals("-h") || s.equals("--help")) {
                ParseCommandLine.printHelp(options);
                System.exit(0);
            }
        }
        // then, parse command line arguments
        CommandLine commandLine = ParseCommandLine.generateCommandLine(options, args);

        // then parse the interrogate the commandLine object
        try {
            Lx = Integer.parseInt(commandLine.getOptionValues("L")[0]);
            Lz = Integer.parseInt(commandLine.getOptionValues("L")[1]);
            extBx = ((Number) commandLine.getParsedOptionValue("extBx")).doubleValue();
            if (commandLine.hasOption("suppress")) suppressInternalTransFields = commandLine.hasOption("suppress");
            if (commandLine.hasOption("seed")){
                seed = ((Number) commandLine.getParsedOptionValue("seed")).longValue();
            }else{
                throw new RuntimeException("No seed was given. A seed is mandatory for fixing backwards compatibility issues.");
            }
            if (commandLine.hasOption("name")) folderName = commandLine.getOptionValue("name");
            newColNum = ((Number) commandLine.getParsedOptionValue("new_col_num")).intValue();

        }
        catch (ArrayIndexOutOfBoundsException e){
            System.err.println("ArrayIndexOutOfBoundsException caught");
            e.printStackTrace();
            System.exit(1);
        }
        catch (NumberFormatException e){
            System.err.println("Non numeric arguments caught.");
            e.printStackTrace();
            System.exit(1);
        }
        catch (ParseException p){
            ParseCommandLine.err(args, p);
            System.exit(1);
        }

        reSaveCheckpoint(folderName, Lx, Lz, extBx, suppressInternalTransFields, seed, newColNum);
    }

    /**
     * Add a column to an old saved simulation to make it compatible with newer simulation object.
     * Reads the checkpoint, modifies the object and saves it.
     * @param folderName - project name
     * @param Lx - Number of unit cells in x and in y directions
     * @param Lz - Number of unit cells in z direction
     * @param extBx - External Bx magnetic field
     * @param suppressInternalTransFields - Suppress internal magnetic fields
     * @param seed - RNG seed
     */
    public static void reSaveCheckpoint(final String folderName, final int Lx, final int Lz, final double extBx, final boolean suppressInternalTransFields, final long seed, final int newColNum){
        SimulationCheckpointer checkpointer = new SimulationCheckpointer(folderName, Lx, Lz, extBx, suppressInternalTransFields, seed);

        MultipleTMonteCarloSimulation simulation = (MultipleTMonteCarloSimulation)checkpointer.readCheckpoint();
        for (int t=0;t<simulation.getT().length;t++){
            simulation.getIthSubSimulation(t).addObservableToBinAvg(newColNum);
        }
        checkpointer.writeCheckpoint(simulation);
    }
}
