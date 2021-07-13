package simulation.montecarlo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import utilities.GenerateSeeds;

import java.io.*;
import java.util.ArrayList;
import java.util.Properties;

public class Main {

    /**
     * Fills the sin and cos tables for mk^2
     * @param arr the array of single spins of the system
     * @param Lz # of unit cells in the z direction
     * @param Lx # of unit cells in the x and y directions
     * @param k_cos_table table to be filled with cos(k.R) for each site R in each direction k (k_x,k_y,k_z)
     * @param k_sin_table table to be filled with sin(k.R) for each site R in each direction k (k_x,k_y,k_z)
     */
    public static void create_cos_sin_tables(singleSpin[] arr, int Lz, int Lx, double[][] k_cos_table, double[][] k_sin_table){
        // for correlation length along x
        double actual_length = Lx*Constants.a;
        double k = 2*Math.PI/actual_length;
        for (int i=0; i<arr.length;i++){
            k_cos_table[0][i] = Math.cos(k*arr[i].getX(Lz,Lx));
            k_sin_table[0][i] = Math.sin(k*arr[i].getX(Lz,Lx));
        }

        // for correlation length along y
        actual_length = Lx*Constants.a;
        k = 2*Math.PI/actual_length;
        for (int i=0; i<arr.length;i++){
            k_cos_table[1][i] = Math.cos(k*arr[i].getY(Lz,Lx));
            k_sin_table[1][i] = Math.sin(k*arr[i].getY(Lz,Lx));
        }

        // for correlation length along z
        actual_length = Lz*Constants.c;
        k = 2*Math.PI/actual_length;
        for (int i=0; i<arr.length;i++){
            k_cos_table[2][i] = Math.cos(k*arr[i].getZ(Lz,Lx));
            k_sin_table[2][i] = Math.sin(k*arr[i].getZ(Lz,Lx));
        }
    }

    /**
     * Initializes an array of {@code MersenneTwister} random number generators
     * @param rndArr array of RNGs to initialize
     * @param taskID to use as part of the seed for RNGs
     * @return the seeds used to initialize the RNGs
     * @deprecated RNGs are initialized in {@link #main(String[])}
     * @see GenerateSeeds#generateSeeds(RandomGenerator, int)
     */
    @Deprecated
    public static long[] initializeMultipleRNG(MersenneTwister[] rndArr, int taskID){
        long[] seeds = new long[rndArr.length];
        for (int i=0;i<rndArr.length;i++){
            seeds[i] = initializeRNG(rndArr[i], taskID+i);
        }
        return seeds;
    }

    /**
     * Initializes a single {@code MersenneTwister} random number generator
     * @param rnd RNG to initialize
     * @param taskID to use as part of the seed
     * @return the used seeds
     * @deprecated RNGs are initialized in {@link #main(String[])}
     * @see GenerateSeeds#generateSeeds(RandomGenerator, int)
     */
    @Deprecated
    public static long initializeRNG(MersenneTwister rnd, int taskID){
        long seed;

        long time = System.currentTimeMillis();
        RandomGenerator tmpRnd = new JDKRandomGenerator();
        tmpRnd.setSeed(new int[]{(int)time, taskID});

        RandomGenerator rndMT19937 = new MersenneTwister(tmpRnd.nextLong());

        seed = Math.abs(rnd.nextLong());    // save seed
        rnd.setSeed(seed);    // reset the seed

        return seed;
    }

    /**
     * Read a temperature schedule from the given file
     * @param fname name of the file from which the temperature schedule should be read
     * @return array of temperatures to simulate
     */
    public static double[] receiveTemperatureSchedule(String fname){
        // the temperature schedule file has a single line, which lists all of the temperatures separated with commas
        double[] temperatureSchedule = null;
        try (BufferedReader br = new BufferedReader(new FileReader(fname));){
            String line;

            if ((line = br.readLine()) != null){
                String[] temps = line.split(",");
                temperatureSchedule = new double[temps.length];
                for (int i=0;i<temps.length;i++){
                    temperatureSchedule[i] = Double.parseDouble(temps[i]);
                }
            }
        }catch(Exception e){
            e.printStackTrace();
        }

        if (temperatureSchedule!=null) return temperatureSchedule;
        else throw new RuntimeException("could not read temperature schedule");
    }

    /**
     * Makes a new directory if one does not already exist
     * @param path in which to create the new directory with trailing separator char
     * @param dirName name of the new directory
     */
    public static void makeDir(String path, String dirName){
        dirName = path.concat(dirName);
        File directory = new File(dirName);
        if (! directory.exists()){
            directory.mkdir();
        }
    }

    /**
     * Creates a queue that can hold the last x bin averages and standard errors and be used to determine equilibration
     * ("last 3 bins for all observables agree within error bars")
     * @param numOfBins number of the last bins to keep in the queue
     * @param numOfTemperatures number of temperatures to track
     * @param numOfObservables number of observables to track
     * @return a queue that can hold the averages and STEs from the last X bins for all temperatures and observables
     * @deprecated since equilibration is no longer tracker in real time
     */
    @Deprecated
    public static ArrayList<ArrayList<CircularFifoQueue<Pair<Double,Double>>>> generateEquilibrationQueues(int numOfBins, int numOfTemperatures, int numOfObservables){
        ArrayList<ArrayList<CircularFifoQueue<Pair<Double,Double>>>> ret = new ArrayList<>(numOfTemperatures);
        for (int temperature=0;temperature<numOfTemperatures;temperature++){
            ArrayList<CircularFifoQueue<Pair<Double,Double>>> obsArr = new ArrayList<>(numOfObservables);
            for (int observable=0;observable<numOfObservables;observable++){
                obsArr.add(new CircularFifoQueue<>(numOfBins));
            }
            ret.add(obsArr);
        }
        return ret;
    }

    /**
     * Runs the Monte Carlo simulation
     * @param args arguments parsed using {@link ParseCommandLine}
     */
    public static void main(String[] args){
        final String VERSION = "1.4.0b";
        int Lx=0;	// lattice x-y size
        int Lz=0;	// lattice z size
        double extBx=-1;   // external Bx
        double extBy=0;   // external By
        long maxSweeps=0;	// maximum steps for the metropolis algorithm
        int taskID=123;	// sge task ID (used for random number seed)
        boolean suppressInternalTransFields=false; // whether to suppress internal transverse fields (exclude offdiagonal dipolar terms)
        boolean continueFromSave=false; // whether to continue from a saved state or start a new simulation
        int maxIter=80; // maximum number of iterations for the iterative solvers of the self-consistent calculation
        int bufferSize=0;   // size of buffer for the output (results) of the simulation
        String tempScheduleFileName = "temperature_schedule_t.txt"; // name of the temperature schedule file
        long seed=0;	// seed for the random number generator, should not stay 0
        boolean parallelTemperingOff = false;  // turn off parallel tempering
        boolean printProgress = false, printOutput=true, saveState; // for interactive runs: print the progess, print the output to
                                                                    // console and periodically save the simulation's state
        long obsPrintSweepNum=15;   // number of MC sweeps between observables are printed (and the simulation's state saved if so defined)
        String folderName="default";    // name of the project (the name of both the results directory and the checkpoints directory)
        double alpha=1.0;               // relaxation parameter for the self-consistent iterative solver
        boolean verboseOutput=true;    // whether results are printed verbosely or in bins
        char tempParallelMode='s';      // mode of running multiple temperatures, whether serial or parallel
        double tempJ_ex;                // value of exchange parameter (if not given, it is read from the parameters' file)
        double spinSize;                // initial magnetic moment value for all spins
        double tol;                     // tolerance for convergence of self-consistent calculation
        double[] T=null;                // temperature array
        String interpolationTableFileNameExtension = "";    // extension of the file name from which the magnetic moment and energy interpolation tables are read
        double x=1.0;   // concentration of magnetic ions, between 0 and 1


        // get Properties object that reads parameters from file
        Properties params = GetParamValues.getParams();
        // get default values from file
        // these values may be overwritten by command line arguments
        tempJ_ex = GetParamValues.getDoubleParam(params, "J_ex");
        tol = GetParamValues.getDoubleParam(params, "tol");


        // create Options object
        Options options = ParseCommandLine.generateOptions();
        // first check if help was called
        for (String s : args) {
            if (s.equals("-h") || s.equals("--help")) {
                ParseCommandLine.printHelp(options);
                System.exit(0);
            }
        }
        // if no command line arguments were given at all, print USAGE
        if (args.length == 0){
            ParseCommandLine.printUsage(options);
            System.exit(1);
        }
        // then, parse command line arguments
        CommandLine commandLine = ParseCommandLine.generateCommandLine(options, args);

        // then parse the interrogate the commandLine object
        try {
            if (commandLine.hasOption("mode")) tempParallelMode = commandLine.getOptionValue("mode").charAt(0);
            Lx = Integer.parseInt(commandLine.getOptionValues("L")[0]);
            Lz = Integer.parseInt(commandLine.getOptionValues("L")[1]);

            verboseOutput = commandLine.hasOption("verbose");

            // get approximate values based on L (lower L means faster progress, so higher write frequency)
            // this is just the default values. if another is given as command-line argument that is the value used
            obsPrintSweepNum = GetParamValues.getLongParam(params,"obsPrintSweepNum"+Lz);

            extBx = ((Number) commandLine.getParsedOptionValue("extBx")).doubleValue();
            if (commandLine.hasOption("extBy")) extBy = ((Number) commandLine.getParsedOptionValue("extBy")).doubleValue();
            maxSweeps = ((Number) commandLine.getParsedOptionValue("max_sweeps")).longValue();
            if (commandLine.hasOption("id")) taskID = ((Number) commandLine.getParsedOptionValue("id")).intValue();
            if (commandLine.hasOption("suppress")) suppressInternalTransFields = commandLine.hasOption("suppress");
            continueFromSave = (commandLine.getOptionValue("s").equals("y") || commandLine.getOptionValue("s").equals("yes"));
            maxIter = ((Number) commandLine.getParsedOptionValue("max_iter")).intValue();
            if (commandLine.hasOption("buffer_size")) bufferSize = ((Number) commandLine.getParsedOptionValue("buffer_size")).intValue();
            if (commandLine.hasOption("temp_schedule")) tempScheduleFileName = commandLine.getOptionValue("temp_schedule");

            // Read temperature schedule. throws runtime exception if unsuccessful
            T=receiveTemperatureSchedule(tempScheduleFileName);
            if (tempParallelMode=='s'){
                // if we are running in serial mode the simulation is slowed down by a factor of the number of temperatures.
                // thus we print observables (T.length) times more often.
                // notice 'obsPrintSweepNum' given as command line argument overrides this behaviour
                obsPrintSweepNum = obsPrintSweepNum/T.length;
            }

            seed = ((Number) commandLine.getParsedOptionValue("seed")).longValue();
            parallelTemperingOff = commandLine.hasOption("pt_off");
            printProgress = commandLine.hasOption("p");
            printOutput = commandLine.hasOption("output");
            if (commandLine.hasOption("print_sweep_num")) obsPrintSweepNum = ((Number) commandLine.getParsedOptionValue("print_sweep_num")).longValue();
            if (commandLine.hasOption("name")) folderName = commandLine.getOptionValue("name");
            if (commandLine.hasOption("alpha")) alpha = ((Number) commandLine.getParsedOptionValue("alpha")).doubleValue();

            if (commandLine.hasOption("tol")) tol = ((Number) commandLine.getParsedOptionValue("tol")).doubleValue();
            if (commandLine.hasOption("Jex")) tempJ_ex = ((Number) commandLine.getParsedOptionValue("Jex")).doubleValue();
            if (commandLine.hasOption("interpolation_table_name")) interpolationTableFileNameExtension = commandLine.getOptionValue("interpolation_table_name");
            if (commandLine.hasOption("d")) x = ((Number) commandLine.getParsedOptionValue("d")).doubleValue();

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


        saveState = !printOutput;	// when output is printed to the console the state should not be saved.
        final char parallelMode = tempParallelMode;
        final double J_ex=tempJ_ex;

        // initialize RNGs
        MersenneTwister mutualRnd = new MersenneTwister(seed);
        long[] seeds = GenerateSeeds.generateSeeds(mutualRnd, T.length);
        MersenneTwister[] rnd = new MersenneTwister[T.length];

        // Create dilution structure (based on given seed)
        boolean[] dilution = new boolean[Lx*Lx*Lz*Constants.num_in_cell];
        int N=0;    // total number of spins in the (diluted) system
        for (int i=0;i<dilution.length;i++){
            if (mutualRnd.nextDouble() < x){
                dilution[i] = true;
                N++;
            }else{
                dilution[i] = false;
            }
        }

        // first try and get spin size (initial guess) from manual calculation that diagonalizes the crystal-field hamiltonian
        // using the external Bx, By.
        // pass parameters Bx=extBx, By=extBy, Bz=0.05, spin=1, and calc for "up"
        spinSize = CrystalField.getMagneticMoment(extBx, extBy, 0.05);

        final double[][][] intTable = new double[3][N][N];      // interaction table that holds all the dipolar interactions. will be full even though it's symmetric. 1st array is x,y,z term
        final double[][] exchangeIntTable = new double[N][N];   // exchange interactions table (only exists between nearest neighbors)

        ReadInteractionsTable interactionsTableReceiver;
        if (System.getProperty("system").equals("LiHoF4")){
            interactionsTableReceiver = new ReadInteractionsTableLiHoF4(dilution);
        } else if (System.getProperty("system").equals("Fe8")){
            interactionsTableReceiver = new ReadInteractionsTableFe8();
        } else {
            throw new RuntimeException("Could not read interactions table. Illegal system name given: " + System.getProperty("system"));
        }
        ReadInteractionsTable.receiveIntTable(intTable, Lx, Lz, dilution);	// get interaction table from file

        if (suppressInternalTransFields){
            // if we suppress internal transverse fields we can just set the interaction table's x,y all to 0.0
            // as a precaution we also always pass suppressInternalTransFields to ignore these interactions
            for (int dim=0;dim<intTable.length-1;dim++){    // go over x and y only (0 and 1)
                for (int i=0;i<intTable[dim].length;i++){
                    for (int j=0;j<intTable[dim][i].length;j++){
                        intTable[dim][i][j]=0;
                    }
                }
            }
        }
        int[][] nnArray = interactionsTableReceiver.exchangeInt(exchangeIntTable, Lx, Lx, Lz, J_ex);	// receive the nearest neighbor array and fill exchangeIntTable with the exchange interaction values

        // add the exchange interaction to intTable which now holds both dipolar and exchange interactions
        for (int i=0;i<N;i++){
            for (int j=0;j<N;j++){
                intTable[2][i][j] += -0.5*exchangeIntTable[i][j]*Constants.k_B/(Constants.mu_B*Constants.g_L);
            }
        }

        // fill cos and sin tables for faster calculation of the k-space magnetization
        final double[][] k_cos_table, k_sin_table;
        {   // code block: tempLattice is discarded at the end

            // temporary lattice object used to create k_tables
            Lattice tempLattice = new Lattice(Lx, Lz, x, extBx, extBy, suppressInternalTransFields, spinSize, dilution, null, null, null, null, null, null);

            // initialize sin, cos tables for mk^2 calculation (correlation length)
            k_cos_table = new double[3][tempLattice.getN()];
            k_sin_table = new double[3][tempLattice.getN()];
            create_cos_sin_tables(tempLattice.getArray(), Lz, Lx, k_cos_table, k_sin_table);
        }

        // read the energy and magnetic moment tables from the respective files.
        // the naming convention is that they include the extrnal Bx, around which the
        // values in the tables are centered
        FieldTable energyTable = FieldTable.of(String.format("energy_up_arr_%1.2f"+interpolationTableFileNameExtension,extBx), false);
        FieldTable momentTable = FieldTable.of(String.format("magnetic_moment_up_arr_%1.2f"+interpolationTableFileNameExtension,extBx), true);

        // initialize the measurement object that will be given to the Lattice object
        ObservableExtractor measure = new ObservableExtractor(k_cos_table, k_sin_table);

        // initialize the checkpointer in charge of saving the simulation state periodically and reading the saved state
        SimulationCheckpointer checkpointer = new SimulationCheckpointer(folderName, Lx, Lz, extBx, suppressInternalTransFields, seed);
        boolean successReadFromFile = false;
        MonteCarloSimulation simulation = null;

        if (continueFromSave) {
            simulation = checkpointer.readCheckpoint();
            if (simulation!=null) {
                successReadFromFile=true;
                // verify that the read state agrees with the important given parameters
                String inconsistencies = SimulationCheckpointer.verifyCheckpointCompatibility(T, parallelTemperingOff, parallelMode, spinSize, tol, J_ex, seed, dilution, simulation);

                if (!inconsistencies.isEmpty()) {
                    System.err.println("There were some inconsistencies between the checkpoint parameters and the current parameters: " + inconsistencies);
                    System.err.println("Exiting.");
                    System.exit(1);
                }
            }
        }

        // create the results directory
        makeDir(System.getProperty("system") + File.separator + "data" + File.separator + "results" + File.separator, folderName);
        // create the directory that will hold the configuration at the end of the simulation
        makeDir(System.getProperty("system") + File.separator + "data" + File.separator + "lattice_output" + File.separator, folderName);

        SingleTMonteCarloSimulation[] subSimulations = new SingleTMonteCarloSimulation[T.length];
        // output stream that prints the configurations for which no self-consistent solution was found for later analysis
        BufferedWriter outProblematicConfigs=null;
        try {
            makeDir(System.getProperty("system") + File.separator + "data" + File.separator,"p_configs");
            // this is one BufferedWriter for all temperatures (threads). The write method is synchronized on outProblematicConfigs.
            // If there are many problematic configurations if might cause slow down (but many such configurations is problematic regardless)
            outProblematicConfigs = new BufferedWriter(new FileWriter(System.getProperty("system") + File.separator + "data" + File.separator + "p_configs" + File.separator + "problematic_"+(Lx*Lx*Lz*Constants.num_in_cell)+"_"+extBx,true));

            for (int i=0;i<T.length;i++){
                // Create file to write output into
                FileWriter out = new FileWriter(System.getProperty("system") + File.separator + "data" + File.separator + "results" + File.separator + folderName + File.separator + "table_" + Lx + "_" + Lz + "_" + extBx + "_" + T[i] + "_" + suppressInternalTransFields + "_" + seed + ".txt",
                        successReadFromFile);
                OutputWriter outputWriter = new OutputWriter.Builder(verboseOutput ? OutputType.VERBOSE : OutputType.BIN, folderName, obsPrintSweepNum, out)
                        .setPrintOutput(printOutput)
                        .setPrintProgress(printProgress)
                        .setBufferSize(bufferSize)
                        .build();


                if (successReadFromFile){
                    // initialize simulation read from checkpoint
                    // this is for fields of the sub-simulations.
                    // fields that (also) exist in the wrapper MultipleTMonteCarloSimulation are initialized later
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setOutWriter(outputWriter);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setAlpha(alpha);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setMaxIter(maxIter);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().setEnergyTable(energyTable);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().setExchangeIntTable(exchangeIntTable);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().setIntTable(intTable);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().setMomentTable(momentTable);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().setNnArray(nnArray);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().setMeasure(measure);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).getLattice().initIterativeSolver();
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).addSweeps(maxSweeps); // this potentially updates max sweeps (in case more runs are needed than initially planned)
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setOutProblematicConfigs(outProblematicConfigs);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setContinueFromSave(continueFromSave);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setCheckpoint(saveState);
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).initMethodsUsedArr();

                    // print parameters and table headers (with preceding '#') to results output file
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).printRunParameters(VERSION, T, "# successfully read saved state"+System.lineSeparator()+'#'+outputWriter.makeTableHeader().substring(1), simulation.getSeed(), tempScheduleFileName, parallelTemperingOff);
                }else{
                    // initialize a new simulation
                    Lattice lattice = new Lattice(Lx, Lz, x, extBx, extBy, suppressInternalTransFields, spinSize, dilution, intTable, exchangeIntTable, nnArray, energyTable, momentTable, measure);
                    rnd[i] = new MersenneTwister(seeds[i]);
                    subSimulations[i] = new SingleTMonteCarloSimulation(T[i], i, T.length, lattice, 34, maxSweeps, seeds[i], rnd[i], continueFromSave,
                            outputWriter, saveState, maxIter, alpha, outProblematicConfigs, spinSize, tol, J_ex);
                    // print parameters and table headers
                    subSimulations[i].printRunParameters(VERSION, T, "# unsuccessful reading checkpoint... Starting new state."+System.lineSeparator()+outputWriter.makeTableHeader(), seed, tempScheduleFileName, parallelTemperingOff);
                }
            }

            // initialization of the fields of MultipleTMonteCarloSimulation
            if (!successReadFromFile){
                simulation=new MultipleTMonteCarloSimulation(T, subSimulations, maxSweeps, seed, mutualRnd, continueFromSave, parallelTemperingOff, saveState, checkpointer, spinSize, tol, J_ex);
                ((MultipleTMonteCarloSimulation) simulation).initSimulation();
            }else{
                ((MultipleTMonteCarloSimulation)simulation).setCheckpointer(checkpointer);
                ((MultipleTMonteCarloSimulation)simulation).addSweeps(maxSweeps);
                ((MultipleTMonteCarloSimulation)simulation).setContinueFromSave(continueFromSave);
                ((MultipleTMonteCarloSimulation)simulation).setCheckpoint(saveState);
            }

            // Run simulation
            ((MultipleTMonteCarloSimulation) simulation).run(parallelMode);

            // Write final lattice states
            if (suppressInternalTransFields) ReadInteractionsTable.receiveIntTable(intTable, Lx, Lz);	// get interaction table from file AGAIN, for off-diagonal interactions that where previously set to zero
            for (int i=0;i<T.length;i++){
                // Create file to write full lattice configurations into
                // lattices are written in full only at the end of the simulation.
                // If one wishes to only write out the lattice of a finished simulation, it should be run with maxSweeps that equals the number of sweeps already done
                try (FileWriter latticeOut = new FileWriter(System.getProperty("system") + File.separator + "data" + File.separator + "lattice_output" + File.separator + folderName + File.separator + "table_" + Lx + "_" + Lz + "_" + extBx + "_" + T[i] + "_" + suppressInternalTransFields + "_" + seed + ".txt", false)) {
                    OutputWriter latticeOutputWriter = new OutputWriter.Builder(OutputType.SPIN, folderName, Constants.num_in_cell * Lx * Lx * Lz, latticeOut)
                            .build();
                    ((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).setOutWriter(latticeOutputWriter);
                    ((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).getLattice().setIntTable(intTable); // set intTable to new one with off-diagonal interactions
                    ((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).printRunParameters(VERSION, T, "# NOTICE: The local transverse (x,y) fields are \"hypothetical\", meaning that if suppressInternalTransFields is true (see above), they are effectively always (localBx=extBx,localBy=localBy)." + System.lineSeparator() + ((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).getOutWriter().makeTableHeader(), simulation.getSeed(), tempScheduleFileName, parallelTemperingOff);
                    ((MultipleTMonteCarloSimulation) simulation).getIthSubSimulation(i).printSimulationState();

                } catch (IOException e){
                    System.err.println("error writing lattice state to file: " + e.toString());
                    e.printStackTrace();
                }
            }

        }catch (IOException e) {
            System.err.println("error writing to file: " + e.toString());
            e.printStackTrace();
        }
        finally {
            // close all outputs
            if (simulation!=null){
                try{
                    ((MultipleTMonteCarloSimulation)simulation).close();
                }catch (IOException e){
                    System.err.println("error closing simulation and files: " + e.toString());
                }
            }

            if (outProblematicConfigs!=null){
                try{
                    outProblematicConfigs.close();
                }catch (IOException problematicConfigErr){
                    System.err.println("error closing problematic config file: " + problematicConfigErr.toString());
                }
            }
        }


    }

}


