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

    // fills the sin and cos tables for mk^2
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

    public static long[] initializeMultipleRNG(MersenneTwister[] rndArr, int taskID){
        long[] seeds = new long[rndArr.length];
        for (int i=0;i<rndArr.length;i++){
            seeds[i] = initializeRNG(rndArr[i], taskID+i);
        }
        return seeds;
    }

    public static long initializeRNG(MersenneTwister rnd, int taskID){
        long seed;

        long time = System.currentTimeMillis();
        RandomGenerator tmpRnd = new JDKRandomGenerator();
        tmpRnd.setSeed(new int[]{(int)time, taskID});

        RandomGenerator rndMT19937 = new MersenneTwister(tmpRnd.nextLong());

        seed = Math.abs(rnd.nextLong());    // save seed
        //seed=1402684116000210553L;
        rnd.setSeed(seed);    // reset the seed

        return seed;
    }


    public static double[] receiveTemperatureSchedule(String fname){
        double[] temperatureSchedule = null;
        try {
            @SuppressWarnings("resource")
            BufferedReader br = new BufferedReader(new FileReader(fname));
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

    public static void makeDir(String path, String dirName){
        dirName = path.concat(dirName);
        File directory = new File(dirName);
        if (! directory.exists()){
            directory.mkdir();
        }
    }


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





    public static void main(String[] args){
        final String VERSION = "1.2.0b";
        int Lx=0;	// lattice x-y size
        int Lz=0;	// lattice z size
        double extBx=-1;   // external Bx
        double extBy=0;   // external By
        long maxSweeps=0;	// maximum steps for the metropolis algorithm
        int taskID=123;	// sge task ID (used for random number seed)
        boolean suppressInternalTransFields=false;
        boolean continueFromSave=false;
        int maxIter=80;
        int bufferSize=0;
        String tempScheduleFileName = "temperature_schedule_t.txt";
        boolean receivedSeed = false;
        long seed=0;	// should never happen;
        boolean parallelTemperingOff = false;;
        boolean printProgress = false, printOutput=true, saveState;
        long obsPrintSweepNum=15;
        String folderName="default";
        double alpha=1.0;
        boolean realTimeEqTest=false;
        boolean verboseOutput=false;
        char tempParallelMode='t';
        double tempJ_ex;
        double spinSize;
        double tol;
        double[] T=null;    // temperature array
        String interpolationTableFileNameExtension = "";


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
                // thus we print observables T.length times more often.
                // notice 'obsPrintSweepNum' given as command line argument overrides this behaviour
                obsPrintSweepNum = obsPrintSweepNum/T.length;
            }

            if (commandLine.hasOption("seed")){
                receivedSeed=true;
                seed = ((Number) commandLine.getParsedOptionValue("seed")).longValue();
            }else{
                // if no seed is given then this must be a continuing simulation
                if (!continueFromSave){
                    throw new RuntimeException("No seed was given. This is only possible if the simulation is being continued from a checkpoint, in " +
                            "which case the \"continue_from_save\" option should be set to \"yes\". ");
                }
            }
            parallelTemperingOff = commandLine.hasOption("pt_off");
            printProgress = commandLine.hasOption("p");
            printOutput = commandLine.hasOption("output");
            if (commandLine.hasOption("print_sweep_num")) obsPrintSweepNum = ((Number) commandLine.getParsedOptionValue("print_sweep_num")).longValue();
            if (commandLine.hasOption("name")) folderName = commandLine.getOptionValue("name");
            if (commandLine.hasOption("alpha")) alpha = ((Number) commandLine.getParsedOptionValue("alpha")).doubleValue();
            realTimeEqTest = commandLine.hasOption("test_eq");

            if (commandLine.hasOption("tol")) tol = ((Number) commandLine.getParsedOptionValue("tol")).doubleValue();
            if (commandLine.hasOption("Jex")) tempJ_ex = ((Number) commandLine.getParsedOptionValue("Jex")).doubleValue();
            if (commandLine.hasOption("interpolation_table_name")) interpolationTableFileNameExtension = commandLine.getOptionValue("interpolation_table_name");


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


        // first try and get spin size (initial guess) from manual calculation that diagonalizes the Ho C-F hamiltonian
        // using the external Bx, By.
        // pass parameters Bx=extBx, By=extBy, Bz=0.05, spin=1, and calc for "up"
        spinSize = CrystalField.getMagneticMoment(extBx, extBy, 0.05);

        final double[][][] intTable = new double[3][Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz]; // create interaction table that holds all the dipolar interactions. will be full even though it's symmetric. 1st array is x,y,z term
        final double[][] exchangeIntTable = new double[Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz];

        ReadInteractionsTable interactionsTableReceiver;
        if (System.getProperty("system").equals("LiHoF4")){
            interactionsTableReceiver = new ReadInteractionsTableLiHoF4();
        } else if (System.getProperty("system").equals("Fe8")){
            interactionsTableReceiver = new ReadInteractionsTableFe8();
        } else {
            throw new RuntimeException("Could not read interactions table. Illegal system name given.");
        }
        ReadInteractionsTable.receiveIntTable(intTable, Lx, Lz);	// get interaction table from file

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

        // add exchange to intTable

        final int N=Lx*Lx*Lz*Constants.num_in_cell;
        for (int i=0;i<N;i++){
            for (int j=0;j<N;j++){
                intTable[2][i][j] += -0.5*exchangeIntTable[i][j]*Constants.k_B/(Constants.mu_B*Constants.g_L);
            }
        }


        final double[][] k_cos_table, k_sin_table;
        {   // code block: tempLattice is discarded at the end

            // temporary lattice objecct used to create k_tables
            Lattice tempLattice = new Lattice(Lx, Lz, extBx, extBy, suppressInternalTransFields, spinSize,null, null, null, null, null, null);

            // initialize sin, cos tables for mk^2 calculation (correlation length)
            k_cos_table = new double[3][tempLattice.getN()];
            k_sin_table = new double[3][tempLattice.getN()];
            create_cos_sin_tables(tempLattice.getArray(), Lz, Lx, k_cos_table, k_sin_table);
        }

        FieldTable energyTable = FieldTable.of(String.format("energy_up_arr_%1.2f"+interpolationTableFileNameExtension,extBx), false);
        FieldTable momentTable = FieldTable.of(String.format("magnetic_moment_up_arr_%1.2f"+interpolationTableFileNameExtension,extBx), true);

        ObservableExtractor measure = new ObservableExtractor(k_cos_table, k_sin_table);

        MersenneTwister mutualRnd = null;
        long[] seeds=null;
        MersenneTwister[] rnd = null;
        if (!receivedSeed) {
            seed = 0;   // this is a flag that should be checked later
        } else {
            if (seed!=0) {
                mutualRnd = new MersenneTwister(seed);
                seeds = GenerateSeeds.generateSeeds(mutualRnd, T.length);
                rnd = new MersenneTwister[T.length];
            }
            else {
                System.err.println("PRNG seed was not initialized for some reason!");
                System.exit(1);
            }
        }

        SimulationCheckpointer checkpointer = new SimulationCheckpointer(folderName, Lx, Lz, extBx, suppressInternalTransFields, seed);
        boolean successReadFromFile = false;
        MonteCarloSimulation simulation = null;

        if (continueFromSave) {
            simulation = checkpointer.readCheckpoint();
            if (simulation!=null) {
                successReadFromFile=true;
                String inconsistencies = SimulationCheckpointer.verifyCheckpointCompatibility(T, parallelTemperingOff, parallelMode, spinSize, tol, J_ex, seed, simulation);

                if (!inconsistencies.isEmpty()) {
                    System.err.println("There were some inconsistencies between the checkpoint parameters and the current parameters: " + inconsistencies);
                    System.err.println("Exiting.");
                    System.exit(1);

                    successReadFromFile = false;  // maybe later the simulation can be run with the new parameters.
                }
            }

        }

        makeDir(System.getProperty("system") + File.separator + "data" + File.separator + "results" + File.separator, folderName);

        makeDir(System.getProperty("system") + File.separator + "data" + File.separator + "lattice_output" + File.separator, folderName);

        SingleTMonteCarloSimulation[] subSimulations = new SingleTMonteCarloSimulation[T.length];
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
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).setRealTimeEqTest(realTimeEqTest);

                    // print parameters and table headers (with preceding '#') to results output file
                    ((MultipleTMonteCarloSimulation)simulation).getIthSubSimulation(i).printRunParameters(VERSION, T, "# successfully read saved state"+System.lineSeparator()+'#'+outputWriter.makeTableHeader().substring(1), simulation.getSeed(), tempScheduleFileName, parallelTemperingOff);
                }else{
                    // initialize new simulation
                    Lattice lattice = new Lattice(Lx, Lz, extBx, extBy, suppressInternalTransFields, spinSize, intTable, exchangeIntTable, nnArray, energyTable, momentTable, measure);
                    rnd[i] = new MersenneTwister(seeds[i]);
                    subSimulations[i] = new SingleTMonteCarloSimulation(T[i], i, T.length, lattice, 36, maxSweeps, seeds[i], rnd[i], continueFromSave,
                            realTimeEqTest, outputWriter, saveState, maxIter, alpha, outProblematicConfigs, spinSize, tol, J_ex);
                    // print parameters and table headers
                    subSimulations[i].printRunParameters(VERSION, T, "# unsuccessful reading checkpoint... Starting new state."+System.lineSeparator()+outputWriter.makeTableHeader(), seed, tempScheduleFileName, parallelTemperingOff);
                }

                //outputWriter.print(outputWriter.makeTableHeader(), true);


            }

            // initialization of the fields of MultipleTMonteCarloSimulation
            if (!successReadFromFile){
                simulation=new MultipleTMonteCarloSimulation(T, subSimulations, maxSweeps, seed, mutualRnd, continueFromSave, realTimeEqTest, parallelTemperingOff, saveState, checkpointer, spinSize, tol, J_ex);
                ((MultipleTMonteCarloSimulation) simulation).initSimulation();
            }else{
                ((MultipleTMonteCarloSimulation)simulation).setCheckpointer(checkpointer);
                ((MultipleTMonteCarloSimulation)simulation).addSweeps(maxSweeps);
                ((MultipleTMonteCarloSimulation)simulation).setContinueFromSave(continueFromSave);
                ((MultipleTMonteCarloSimulation)simulation).setCheckpoint(saveState);
                ((MultipleTMonteCarloSimulation)simulation).setRealTimeEqTest(realTimeEqTest);
            }

            // Run simulation
            ((MultipleTMonteCarloSimulation) simulation).run(parallelMode);

            // Write lattice states
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


