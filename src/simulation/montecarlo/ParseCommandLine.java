package simulation.montecarlo;

import org.apache.commons.cli.*;

import java.io.PrintWriter;
import java.util.Arrays;


public class ParseCommandLine {

    /**
     * "Definition" stage of command-line parsing with Apache Commons CLI.
     * @return Definition of command-line options.
     */
    public static Options generateOptions()
    {
        final Option magneticMomentCalcMaxIter = Option.builder("max_iter")
                .required(true)
                .hasArg()
                .longOpt("magnetic_moment_calc_max_iter")
                .desc("Maximum number of iterations for the magnetic moment calculation.")
                .type(Number.class)
                .build();
        final Option systemSize = Option.builder("L")
                .required(true)
                .longOpt("system_size")
                .numberOfArgs(2)
                .desc("Number of unit cells in the x and z directions, separated by a comma.")
                .valueSeparator(',')
                .type(Number.class)
                .build();
        final Option dilution = Option.builder("d")
                .required(false)
                .longOpt("dilution")
                .hasArg()
                .desc("Dilution of Ho ions. Value between 0 and 1.")
                .type(Number.class)
                .build();
        final Option h = Option.builder("random_h")
                .required(false)
                .longOpt("simulation/montecarlo")
                .hasArg()
                .desc("Standard deviation of local random longitudinal fields.")
                .type(Number.class)
                .build();
        final Option extBx = Option.builder("extBx")
                .required(true)
                .longOpt("external_Bx")
                .hasArg()
                .desc("External magnetic field in the x direction, Bx.")
                .type(Number.class)
                .build();
        final Option maxSweeps = Option.builder("max_sweeps")
                .required(true)
                .hasArg()
                .desc("Number of Monte Carlo sweeps to run.")
                .type(Number.class)
                .build();
        final Option fileNumber = Option.builder("config_file")
                .required(false)
                .longOpt("config_file_num")
                .hasArg()
                .desc("Serial number of the configuration file to be used.")
                .type(Number.class)
                .build();
        final Option taskID = Option.builder("id")
                .required(false)
                .longOpt("task_id")
                .hasArg()
                .desc("task ID from sge system. Used in generation of random seed.")
                .type(Number.class)
                .build();
        final Option suppressInternalTransFields = Option.builder("suppress")
                .required(false)
                .longOpt("suppress_internal_trans_fields")
                .desc("Suppress internal transverse fields so that only Bex should affect the behavior.")
                .build();
        final Option continueFromSave = Option.builder("s")
                .required(true)
                .hasArg(true)
                .longOpt("continue_from_save")
                .desc("Continue from saved state or start new simulation (OVERWRITES previous results w/ same name). yes/no.")
                .build();
        final Option bufferSize = Option.builder("buffer_size")
                .required(false)
                .hasArg()
                .desc("Buffer size for output (character count).")
                .type(Number.class)
                .build();
        final Option tempSchedule = Option.builder("temp_schedule")
                .required(false)
                .longOpt("temperature_schedule")
                .hasArg()
                .desc("Location of file containing the temperature schedule.")
                .build();
        final Option seed = Option.builder("seed")
                .required(false)
                .longOpt("PRNG_seed")
                .hasArg()
                .desc("Seed for random number generator. Do not set unless for debugging.")
                .type(Number.class)
                .build();
        final Option parallelTemperingOff = Option.builder("pt_off")
                .required(false)
                .longOpt("parallel_tempering_off")
                .desc("Turn off parallel tempering.")
                .build();
        final Option printProgress = Option.builder("p")
                .required(false)
                .longOpt("print_progress")
                .desc("Print progress to terminal.")
                .build();
        final Option obsPrintSweepNum = Option.builder("print_sweep_num")
                .required(false)
                .longOpt("obs_print_sweep_num")
                .hasArg()
                .desc("Number of sweeps between output print.")
                .type(Number.class)
                .build();
        final Option printOutput = Option.builder("output")
                .required(false)
                .longOpt("print_output")
                .desc("Print the output of the simulation to the terminal. When this option is enabled the state is not saved.")
                .build();
        final Option folderName = Option.builder("name")
                .required(false)
                .hasArg()
                .longOpt("folder_name")
                .desc("Name of folder in which to save the output and save files.")
                .build();
        final Option interpolationTableName = Option.builder("interpolation_table_name")
                .required(false)
                .hasArg()
                .desc("Extention of the name of the interpolation table file.")
                .build();
        final Option alpha = Option.builder("alpha")
                .required(false)
                .hasArg()
                .type(Number.class)
                .desc("Alpha parameter for over/under relaxation of iterative moment calculation. Default is 1.0")
                .build();
        final Option verboseOutput = Option.builder("verbose")
                .required(false)
                .desc("print out the observables after each Monte Carlo sweep, as opposed to logarithmic binnig.")
                .build();
        final Option realTimeEqTest = Option.builder("test_eq")
                .desc("test equilibration in real-time and stop when it is reached.")
                .build();
        final Option tol = Option.builder("tol")
                .longOpt("tolerance")
                .required(false)
                .desc("tolerance for self consistent calculation.")
                .hasArg()
                .type(Number.class)
                .build();
        final Option parallelMode = Option.builder("mode")
                .longOpt("parallel_mode")
                .hasArg(true)
                .required(true)
                .desc("mode of operation when simulating multiple temperatures. serial (s) or parallel (p)")
                .build();
        final Option temperature = Option.builder("t")
                .longOpt("temperature")
                .hasArg(true)
                .required(false)
                .type(Number.class)
                .desc("temperature for monte carlo simulation. must be provided for single temperature simulation (if \"mode\" is not specified)")
                .build();
        final Option Jex = Option.builder("Jex")
                .longOpt("J_ex")
                .hasArg(true)
                .required(false)
                .type(Number.class)
                .desc("Exchange interaction. If not provided the value will be read from parameters.properties")
                .build();
        final Option spinSize = Option.builder("spinSize")
                .longOpt("spin_size")
                .hasArg(true)
                .required(false)
                .type(Number.class)
                .desc("Initial spin size. If not provided the value will be read from parameters.properties")
                .build();

        // taskID and seed are mutually exclusive
        //final OptionGroup seed_task = new OptionGroup();
        //seed_task.addOption(taskID);
        //seed_task.addOption(seed);

        // output and progress are mutually exclusive
        final OptionGroup prog_output = new OptionGroup();
        prog_output.addOption(printOutput);
        prog_output.addOption(printProgress);

        final Options options = new Options();
        options.addOption(magneticMomentCalcMaxIter);
        options.addOption(systemSize);
        options.addOption(dilution);
        options.addOption(h);
        options.addOption(extBx);
        options.addOption(maxSweeps);
        options.addOption(fileNumber);
        options.addOption(suppressInternalTransFields);
        options.addOption(continueFromSave);
        options.addOption(bufferSize);
        options.addOption(tempSchedule);
        options.addOption(parallelTemperingOff);
        options.addOption(obsPrintSweepNum);
        //options.addOptionGroup(seed_task);
        options.addOption(seed);
        options.addOptionGroup(prog_output);
        options.addOption(folderName);
        options.addOption(interpolationTableName);
        options.addOption(alpha);
        options.addOption(verboseOutput);
        options.addOption(realTimeEqTest);
        options.addOption(tol);
        options.addOption(parallelMode);
        options.addOption(temperature);
        options.addOption(spinSize);
        options.addOption(Jex);

        return options;
    }

    /**
     * Parsing stage of command-line processing
     *
     * @param options Options from "definition" stage.
     * @param commandLineArguments Command-line arguments provided to application.
     * @return Instance of CommandLine as parsed from the provided Options and
     *    command line arguments; may be {@code null} if there is an exception
     *    encountered while attempting to parse the command line options.
     */
    public static CommandLine generateCommandLine(
            final Options options, final String[] commandLineArguments)
    {
        final CommandLineParser cmdLineParser = new DefaultParser();
        CommandLine commandLine = null;
        try
        {
            commandLine = cmdLineParser.parse(options, commandLineArguments);
        }
        catch (ParseException parseException)
        {
            System.out.println(
                    "ERROR: Unable to parse command-line arguments "
                            + Arrays.toString(commandLineArguments) + " due to: "
                            + parseException + "\n" + "For help run with -h");
        }
        return commandLine;
    }

    /**
     * Generate usage information with Apache Commons CLI.
     *
     * @param options Instance of Options to be used to prepare
     *    usage formatter.
     */
    protected static void printUsage(final Options options)
    {
        final HelpFormatter formatter = new HelpFormatter();
        final String syntax = "simulation/montecarlo";
        System.out.println("\n=====");
        System.out.println("USAGE");
        System.out.println("=====");
        final PrintWriter pw  = new PrintWriter(System.out);
        formatter.printUsage(pw, 80, syntax, options);
        pw.flush();
    }

    /**
     * Generate help information with Apache Commons CLI.
     *
     * @param options Instance of Options to be used to prepare
     *    help formatter.
     */
    public static void printHelp(final Options options)
    {
        final HelpFormatter formatter = new HelpFormatter();
        final String syntax = "simulation/montecarlo";
        final String usageHeader = "Monte Carlo simulation of LiHoF_4";
        final String usageFooter = "";
        System.out.println("\n====");
        System.out.println("HELP");
        System.out.println("====");
        formatter.printHelp(syntax, usageHeader, options, usageFooter);
    }

    public static void err(String[] args, ParseException p){
        System.err.println(
                "ERROR: Unable to parse command-line arguments "
                        + Arrays.toString(args) + " due to: "
                        + p + System.lineSeparator() + "For help run with -h");
        System.exit(1);
    }

}
