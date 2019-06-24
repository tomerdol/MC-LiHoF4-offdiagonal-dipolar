package simulation.montecarlo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.MersenneTwister;

import java.io.*;
import java.util.ArrayList;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

public class Main {

    // fills the sin and cos tables for mk^2
    public static void create_cos_sin_tables(singleSpin[] arr, int Lz, int Lx, double[] k_cos_table, double[] k_sin_table){
        double actual_length = Lx*Constants.a;
        double k = 2*Math.PI/actual_length;

        for (int i=0; i<arr.length;i++){
            k_cos_table[i] = Math.cos(k*arr[i].getX(Lz,Lx));
            k_sin_table[i] = Math.sin(k*arr[i].getX(Lz,Lx));
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
        // time (long) is divided into 2 integers, fileNumber and  sge taskID are all used to seed the MT so that there is no correlation between tasks
        rnd.setSeed(new int[]{(int) (time >> 32), (int) time, taskID});
        // generator warm up:
        for (int index = 0; index < 20; index++) rnd.nextInt();

        seed = Math.abs(rnd.nextLong());    // save seed
        //seed=1402684116000210553L;
        rnd.setSeed(seed);    // reset the seed

        return seed;
    }

    public static void receiveIntTable(double[][][] intTable, int Lx, int Lz){
        final double c = -Constants.mu_0*Constants.mu_B*Constants.g_L*0.25/Math.PI;	// coefficient dipolar spin-spin interaction. The minus sign is because
        // we use Ewald to calc the effective field and not the interaction

        int fileLx=0, fileLz=0;
        final int N=Lx*Lx*Lz*4;
        try (BufferedReader in = new BufferedReader(new FileReader("interactions" + File.separator + "intTable_"+Lx+"_"+Lz+".txt"))){
            String str;
            // verify Lx and Lz make sense
            if ((str = in.readLine()) != null)
                fileLx=Integer.parseInt(str.split("=")[1]);
            if ((str = in.readLine()) != null)
                fileLz=Integer.parseInt(str.split("=")[1]);
            // total number of spins does not match between file and program
            if (fileLx!=Lx || fileLz!=Lz){
                System.err.println("the file from which we read the interactions has different Lx or Lz than the ones provided.");
                System.exit(1);
            }
            // skip parameter lines:
            in.readLine();
            in.readLine();
            in.readLine();

            for (int i=0;i<N;i++){
                for (int j=i;j<N;j++){
                    // notice that we set self interaction to 0. this is because the self interaction contribution to the total energy is always the same (so it's not important)
                    if ((str = in.readLine()) != null){
                        String[] xyzInteractions = str.split(",");
                        for (int k=0;k<xyzInteractions.length;k++) {
                            intTable[k][i][j] = c*Double.parseDouble(xyzInteractions[k]);
                            intTable[k][j][i] = c*Double.parseDouble(xyzInteractions[k]);
                        }

                    } else {
                        System.err.println("reached end of interactions file before int table is full!");
                        System.exit(1);
                    }
                }
            }


            if (in.readLine()!=null) {
                System.err.println("did not finish reading interaction file and table is already full");
                System.exit(1);
            }
        } catch (IOException e) {
            System.err.println("bad interactions input file!");
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void addExchangeToNeighbor(int focusSpin, int neighbor1, int neighbor2, int neighbor3, int neighbor4, int[][] nnArray, boolean[][] nnArray_test, double[][] intTable, double J_ex){
        intTable[focusSpin][neighbor1]+=J_ex;
        intTable[neighbor1][focusSpin]+=J_ex;

        nnArray[focusSpin][0]=neighbor1;
        nnArray[neighbor1][0]=focusSpin;

        nnArray_test[focusSpin][neighbor1]=true;
        nnArray_test[neighbor1][focusSpin]=true;


        intTable[focusSpin][neighbor2]+=J_ex;
        intTable[neighbor2][focusSpin]+=J_ex;

        nnArray[focusSpin][1]=neighbor2;
        nnArray[neighbor2][1]=focusSpin;

        nnArray_test[focusSpin][neighbor2]=true;
        nnArray_test[neighbor2][focusSpin]=true;

        intTable[focusSpin][neighbor3]+=J_ex;
        intTable[neighbor3][focusSpin]+=J_ex;

        nnArray[focusSpin][2]=neighbor3;
        nnArray[neighbor3][2]=focusSpin;

        nnArray_test[focusSpin][neighbor3]=true;
        nnArray_test[neighbor3][focusSpin]=true;


        intTable[focusSpin][neighbor4]+=J_ex;
        intTable[neighbor4][focusSpin]+=J_ex;

        nnArray[focusSpin][3]=neighbor4;
        nnArray[neighbor4][3]=focusSpin;

        nnArray_test[focusSpin][neighbor4]=true;
        nnArray_test[neighbor4][focusSpin]=true;
    }

    //calculates exchange interaction with nearest neighbors
    // also returns an array of nearest neighbors
    public static int[][] exchangeInt(double[][] intTable, int Lx, int Lz, double J_ex){
        final int N = Lx*Lx*Lz*4;

        int[][] nnArray = new int[N][4];	// nearest neighbor array
        boolean[][] nnArray_test = new boolean[N][N];	// for testing

        // neighbor numbers are as follows (with respect to ion positions_1.pdf):
        // neighbor1: up-right
        // neighbor2: up-left
        // neighbor3: down-outward
        // neighbor4: down-inward

        for (int i=0;i<Lx;i++){
            for (int j=0;j<Lx;j++){
                for (int k=0;k<Lz;k++){
                    // notice we are only going through the 0th and 2nd atoms in the base since they participate in all exchange interactions

                    // nearest neighbors to 0th base atom, including periodic boundary conditions
                    int focusSpin = i*Lx*Lz*4+j*Lz*4+k*4+0;
                    int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;

                        neighbor1=i*Lx*Lz*4+j*Lz*4+k*4+1;
                        if (i==0)
                            neighbor2=(Lx-1)*Lx*Lz*4+j*Lz*4+k*4+1;
                        else
                            neighbor2=(i-1)*Lx*Lz*4+j*Lz*4+k*4+1;
                        if (k==0)
                            neighbor3=i*Lx*Lz*4+j*Lz*4+(Lz-1)*4+3;
                        else
                            neighbor3=i*Lx*Lz*4+j*Lz*4+(k-1)*4+3;

                        if (j==0 && k==0)
                            neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(Lz-1)*4+3;
                        else if(j==0){
                            neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(k-1)*4+3;
                        }
                        else if(k==0){
                            neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(Lz-1)*4+3;
                        }else{
                            neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(k-1)*4+3;
                        }

                        // put interactions in intTable and nearest neighbor indices in nnArray
                        addExchangeToNeighbor(focusSpin, neighbor1, neighbor2, neighbor3, neighbor4, nnArray, nnArray_test, intTable, J_ex);



                    // now nearest neighbors to 2nd base atom, including periodic boundary conditions
                    focusSpin = i*Lx*Lz*4+j*Lz*4+k*4+2;

                    neighbor3=i*Lx*Lz*4+j*Lz*4+k*4+1;
                    neighbor2=i*Lx*Lz*4+j*Lz*4+k*4+3;
                    if (i==Lx-1)
                        neighbor1=(0)*Lx*Lz*4+j*Lz*4+k*4+3;
                    else
                        neighbor1=(i+1)*Lx*Lz*4+j*Lz*4+k*4+3;
                    if (j==Lx-1)
                        neighbor4=i*Lx*Lz*4+(0)*Lz*4+k*4+1;
                    else
                        neighbor4=i*Lx*Lz*4+(j+1)*Lz*4+k*4+1;

                        // put interactions in intTable and nearest neighbor indices in nnArray
                    addExchangeToNeighbor(focusSpin, neighbor1, neighbor2, neighbor3, neighbor4, nnArray, nnArray_test, intTable, J_ex);

                }
            }
        }

        boolean validNNArray=true;
        for (int i=0;i<nnArray_test.length && validNNArray;i++){
            int countNearestNeighbors=0;
            for (int j=0;j<nnArray_test[i].length && validNNArray;j++){
                if (nnArray_test[i][j]) countNearestNeighbors++;
            }
            if (countNearestNeighbors!=4) validNNArray=false;
        }
        if (validNNArray){ return nnArray; }
        else {
            System.err.println("There was an error creating the nearest neighbor array. at least one of the spins has more or less than 4 neighbors.");
            System.exit(1);
            return null;
        }
    }


    /*
    @Deprecated

     * @deprecated {@see #exchangeInt(double[][], int, int, double)}

    public static int[][] exchangeInt2(double[][] intTable, int Lx, int Lz, double J_ex){
        final int N = Lx*Lx*Lz*4;

        int[][] nnArray = new int[N][4];	// nearest neighbor array
        boolean[][] nnArray_test = new boolean[N][N];	// for testing

        // neighbor numbers are as follows (with respect to ion positions_1.pdf):
        // neighbor1: up-right
        // neighbor2: up-left
        // neighbor3: down-outward
        // neighbor4: down-inward

        for (int i=0;i<Lx;i++){
            for (int j=0;j<Lx;j++){
                for (int k=0;k<Lz;k++){
                    // notice we are only going through the 0th and 2nd atoms in the base since they participate in all exchange interactions

                    // nearest neighbors to 0th base atom, including periodic boundary conditions
                    int focusSpin = i*Lx*Lz*4+j*Lz*4+k*4+0;
                    int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;
                    if (arr[focusSpin].getSpin()!=0){
                        neighbor1=i*Lx*Lz*4+j*Lz*4+k*4+1;
                        if (i==0)
                            neighbor2=(Lx-1)*Lx*Lz*4+j*Lz*4+k*4+1;
                        else
                            neighbor2=(i-1)*Lx*Lz*4+j*Lz*4+k*4+1;
                        if (k==0)
                            neighbor3=i*Lx*Lz*4+j*Lz*4+(Lz-1)*4+3;
                        else
                            neighbor3=i*Lx*Lz*4+j*Lz*4+(k-1)*4+3;

                        if (j==0 && k==0)
                            neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(Lz-1)*4+3;
                        else if(j==0){
                            neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(k-1)*4+3;
                        }
                        else if(k==0){
                            neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(Lz-1)*4+3;
                        }else{
                            neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(k-1)*4+3;
                        }

                        // put interactions in intTable and nearest neighbor indices in nnArray
                        if (arr[neighbor1].getSpin()!=0){
                            intTable[focusSpin][neighbor1]+=J_ex;
                            intTable[neighbor1][focusSpin]+=J_ex;

                            nnArray[focusSpin][0]=neighbor1;
                            nnArray[neighbor1][0]=focusSpin;

                            nnArray_test[focusSpin][neighbor1]=true;
                            nnArray_test[neighbor1][focusSpin]=true;
                        }else{
                            nnArray[focusSpin][0]=-1;	// in case neighbor 1 does not exist
                            nnArray[neighbor1][0]=-1;
                        }
                        if (arr[neighbor2].getSpin()!=0){
                            intTable[focusSpin][neighbor2]+=J_ex;
                            intTable[neighbor2][focusSpin]+=J_ex;

                            nnArray[focusSpin][1]=neighbor2;
                            nnArray[neighbor2][1]=focusSpin;

                            nnArray_test[focusSpin][neighbor2]=true;
                            nnArray_test[neighbor2][focusSpin]=true;

                        }else{
                            nnArray[focusSpin][1]=-1;	// in case neighbor 2 does not exist
                            nnArray[neighbor2][1]=-1;
                        }
                        if (arr[neighbor3].getSpin()!=0){
                            intTable[focusSpin][neighbor3]+=J_ex;
                            intTable[neighbor3][focusSpin]+=J_ex;

                            nnArray[focusSpin][2]=neighbor3;
                            nnArray[neighbor3][2]=focusSpin;

                            nnArray_test[focusSpin][neighbor3]=true;
                            nnArray_test[neighbor3][focusSpin]=true;
                        }else{
                            nnArray[focusSpin][2]=-1;	// in case neighbor 3 does not exist
                            nnArray[neighbor3][2]=-1;
                        }
                        if (arr[neighbor4].getSpin()!=0){
                            intTable[focusSpin][neighbor4]+=J_ex;
                            intTable[neighbor4][focusSpin]+=J_ex;

                            nnArray[focusSpin][3]=neighbor4;
                            nnArray[neighbor4][3]=focusSpin;

                            nnArray_test[focusSpin][neighbor4]=true;
                            nnArray_test[neighbor4][focusSpin]=true;
                        }else{
                            nnArray[focusSpin][3]=-1;	// in case neighbor 3 does not exist
                            nnArray[neighbor4][3]=-1;
                        }
                    }
                    // now nearest neighbors to 2nd base atom, including periodic boundary conditions
                    focusSpin = i*Lx*Lz*4+j*Lz*4+k*4+2;
                    if (arr[focusSpin].getSpin()!=0){
                        neighbor3=i*Lx*Lz*4+j*Lz*4+k*4+1;
                        neighbor2=i*Lx*Lz*4+j*Lz*4+k*4+3;
                        if (i==Lx-1)
                            neighbor1=(0)*Lx*Lz*4+j*Lz*4+k*4+3;
                        else
                            neighbor1=(i+1)*Lx*Lz*4+j*Lz*4+k*4+3;
                        if (j==Lx-1)
                            neighbor4=i*Lx*Lz*4+(0)*Lz*4+k*4+1;
                        else
                            neighbor4=i*Lx*Lz*4+(j+1)*Lz*4+k*4+1;

                        // put interactions in intTable and nearest neighbor indices in nnArray
                        if (arr[neighbor1].getSpin()!=0){
                            intTable[focusSpin][neighbor1]+=J_ex;
                            intTable[neighbor1][focusSpin]+=J_ex;

                            nnArray[focusSpin][0]=neighbor1;
                            nnArray[neighbor1][0]=focusSpin;

                            nnArray_test[focusSpin][neighbor1]=true;
                            nnArray_test[neighbor1][focusSpin]=true;
                        }else{
                            nnArray[focusSpin][0]=-1;	// in case neighbor 1 does not exist
                            nnArray[neighbor1][0]=-1;
                        }
                        if (arr[neighbor2].getSpin()!=0){
                            intTable[focusSpin][neighbor2]+=J_ex;
                            intTable[neighbor2][focusSpin]+=J_ex;

                            nnArray[focusSpin][1]=neighbor2;
                            nnArray[neighbor2][1]=focusSpin;

                            nnArray_test[focusSpin][neighbor2]=true;
                            nnArray_test[neighbor2][focusSpin]=true;

                        }else{
                            nnArray[focusSpin][1]=-1;	// in case neighbor 2 does not exist
                            nnArray[neighbor2][1]=-1;
                        }
                        if (arr[neighbor3].getSpin()!=0){
                            intTable[focusSpin][neighbor3]+=J_ex;
                            intTable[neighbor3][focusSpin]+=J_ex;

                            nnArray[focusSpin][2]=neighbor3;
                            nnArray[neighbor3][2]=focusSpin;

                            nnArray_test[focusSpin][neighbor3]=true;
                            nnArray_test[neighbor3][focusSpin]=true;
                        }else{
                            nnArray[focusSpin][2]=-1;	// in case neighbor 3 does not exist
                            nnArray[neighbor3][2]=-1;
                        }
                        if (arr[neighbor4].getSpin()!=0){
                            intTable[focusSpin][neighbor4]+=J_ex;
                            intTable[neighbor4][focusSpin]+=J_ex;

                            nnArray[focusSpin][3]=neighbor4;
                            nnArray[neighbor4][3]=focusSpin;

                            nnArray_test[focusSpin][neighbor4]=true;
                            nnArray_test[neighbor4][focusSpin]=true;
                        }else{
                            nnArray[focusSpin][3]=-1;	// in case neighbor 3 does not exist
                            nnArray[neighbor4][3]=-1;
                        }
                    }
                }
            }
        }

        return nnArray;
    }
    */


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

        int Lx=0;	// lattice x-y size
        int Lz=0;	// lattice z size
        double extBx=-1;   // external Bx
        long maxSweeps=0;	// maximum steps for the metropolis algorithm
        int taskID=123;	// sge task ID (used for random number seed)
        boolean suppressInternalTransFields=false;
        boolean continueFromSave=false;
        int maxIter=80;
        int bufferSize=4500;	// buffer size for output (char count) - should be enough for ~15 sweeps. might need update if more observables are printed after each MCS.
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
        char parallelMode='t';
        double temperature=-1;

        // create Options object
        Options options = ParseCommandLine.generateOptions();
        CommandLine commandLine = ParseCommandLine.generateCommandLine(options, args);

        // first check if help was called
        for (String s : args) {
            if (s.equals("-h") || s.equals("--help")) {  // or use help.getOpt() || help.getLongOpt()
                ParseCommandLine.printHelp(options);
                System.exit(0);
            }
        }

        // then parse the interrogate the commandLine object
        try {
            if (!commandLine.hasOption("mode") && (commandLine.hasOption("pt_off") || commandLine.hasOption("temp_schedule"))){
                throw new ParseException("\"mode\" tag not provided, hence only a single temperature would be run and flags \"pt_off\" and \"temp_schedule\"" +
                        "should not be given.");
            }
            if (!commandLine.hasOption("mode") && !commandLine.hasOption("temperature")){
                throw new ParseException("\"mode\" tag not provided, hence only a single temperature would be run and flags \"-t\"" +
                        "must be provided.");
            }

            if (commandLine.hasOption("mode")) parallelMode = ((Character) commandLine.getParsedOptionValue("mode")).charValue();
            if (commandLine.hasOption("t")) temperature = ((Number) commandLine.getParsedOptionValue("t")).doubleValue();

            Lx = Integer.parseInt(commandLine.getOptionValues("L")[0]);
            Lz = Integer.parseInt(commandLine.getOptionValues("L")[1]);

            verboseOutput = commandLine.hasOption("verbose");

            // get approximate values based on L (lower L means faster progress, so higher write frequency)
            // this is just the default values. if another is given as command-line argument that is the value used
            obsPrintSweepNum = GetParamValues.getLongParam(GetParamValues.getParams(),"obsPrintSweepNum"+Lz);

            extBx = ((Number) commandLine.getParsedOptionValue("extBx")).doubleValue();
            maxSweeps = ((Number) commandLine.getParsedOptionValue("max_sweeps")).longValue();
            if (commandLine.hasOption("id")) taskID = ((Number) commandLine.getParsedOptionValue("id")).intValue();
            if (commandLine.hasOption("suppress")) suppressInternalTransFields = commandLine.hasOption("suppress");
            continueFromSave = (commandLine.getOptionValue("s").equals("y") || commandLine.getOptionValue("s").equals("yes"));
            maxIter = ((Number) commandLine.getParsedOptionValue("max_iter")).intValue();
            if (commandLine.hasOption("buffer_size")) bufferSize = ((Number) commandLine.getParsedOptionValue("buffer_size")).intValue();
            if (commandLine.hasOption("temp_schedule")) tempScheduleFileName = commandLine.getOptionValue("temp_schedule");
            if (commandLine.hasOption("seed")){
                receivedSeed=true;
                seed = ((Number) commandLine.getParsedOptionValue("seed")).longValue();
            }
            parallelTemperingOff = commandLine.hasOption("pt_off");
            printProgress = commandLine.hasOption("p");
            printOutput = commandLine.hasOption("output");
            if (commandLine.hasOption("print_sweep_num")) obsPrintSweepNum = ((Number) commandLine.getParsedOptionValue("print_sweep_num")).longValue();
            if (commandLine.hasOption("name")) folderName = commandLine.getOptionValue("name");
            if (commandLine.hasOption("alpha")) alpha = ((Number) commandLine.getParsedOptionValue("alpha")).doubleValue();
            realTimeEqTest = commandLine.hasOption("test_eq");

            if (commandLine.hasOption("tol")) Constants.tol= ((Number) commandLine.getParsedOptionValue("tol")).doubleValue();

        }
        catch (ArrayIndexOutOfBoundsException e){
            System.err.println("ArrayIndexOutOfBoundsException caught");
            e.printStackTrace();
        }
        catch (NumberFormatException e){
            System.err.println("Non numeric arguments caught.");
            e.printStackTrace();
        }
        catch (ParseException p){
            ParseCommandLine.err(args, p);
        }

        saveState = !printOutput;	// when output is printed to the console the state should not be saved.

        // first try and get spin size (initial guess) from manual calculation that diagonalizes the Ho C-F hamiltonian
        // using the external Bx.
        // If that fails just use value from fit. It's not that important as it's just an initial guess
        try{
            // pass parameters Bx=extBx, By=0, Bz=0.05, spin=1, and calc for "up"
            Constants.spinSize = CrystalField.getMagneticMoment(extBx, 0.0, 0.05);
        }catch(Exception e){
            System.err.println("could not run manual calculation to get initial magnetic moment guess");
            Constants.spinSize = 5.44802407 + 1.11982863*extBx + (-2.00747245)*Math.pow(extBx,2) + (-4.20363871)*Math.pow(extBx,3)+
                    (5.20857114)*Math.pow(extBx,4) + (-2.22871652)*Math.pow(extBx,5) + (0.42596452)*Math.pow(extBx,6) + (-0.0307749)*Math.pow(extBx,7);
        }

        final double[][][] intTable = new double[3][4*Lx*Lx*Lz][4*Lx*Lx*Lz]; // create interaction table that holds all the dipolar interactions. will be full even though it's symmetric. 1st array is x,y,z term
        final double[][] exchangeIntTable = new double[4*Lx*Lx*Lz][4*Lx*Lx*Lz];

        receiveIntTable(intTable, Lx, Lz);	// get interaction table from file

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

        int[][] nnArray = exchangeInt(exchangeIntTable, Lx, Lz, Constants.J_ex);	// receive the nearest neighbor array and fill exchangeIntTable with the exchange interaction values



        // TODO:
        //  v 1. create lattices
        //  v 2. create and initialize RNG
        //  v 3. create k_tables
        //  v 4. FieldTables
        //  5. OutputWriters
        //  6. save states (try read etc.)
        //  7. problematic configs

        final double[] k_cos_table, k_sin_table;
        {   // code block: tempLattice is discarded at the end

            // temporary lattice objecct used to create k_tables
            Lattice tempLattice = new Lattice(Lx, Lz, extBx, suppressInternalTransFields, null, null, null, null, null, null);

            // initialize sin, cos tables for mk^2 calculation (correlation length)
            k_cos_table = new double[tempLattice.getN()];
            k_sin_table = new double[tempLattice.getN()];
            create_cos_sin_tables(tempLattice.getArray(), Lz, Lx, k_cos_table, k_sin_table);
        }

        FieldTable energyTable = FieldTable.of(String.format("energy_up_broyden_arr_%1.2f",extBx), false);
        FieldTable momentTable = FieldTable.of(String.format("magnetic_moment_up_broyden_arr_%1.2f",extBx), true);

        MersenneTwister mutualRnd = new MersenneTwister();
        if (!receivedSeed) {
            seed = initializeRNG(mutualRnd, taskID);
        } else {
            if (seed!=0) {
                mutualRnd.setSeed(seed);
            }
            else {
                System.err.println("PRNG seed was not initialized for some reason!");
                System.exit(1);
            }
        }

        boolean successReadFromFile=false;  // TODO: set


        MonteCarloSimulation simulation=null;
        if (parallelMode=='t'){
            // single temperature mode
            Lattice lattice = new Lattice(Lx, Lz, extBx, suppressInternalTransFields, intTable, exchangeIntTable, nnArray, energyTable, momentTable, new ObservableExtractor(k_cos_table, k_sin_table));

            try (FileWriter out = new FileWriter("analysis" + File.separator + folderName + File.separator + "table_" + Lx + "_" + Lz + "_" + extBx + "_" + temperature + "_" + suppressInternalTransFields + ".txt",successReadFromFile)) {
                OutputWriter outputWriter = new OutputWriter.Builder(verboseOutput, folderName, obsPrintSweepNum, out)
                        .setPrintOutput(printOutput)
                        .setPrintProgress(printProgress)
                        .setBufferSize(bufferSize)
                        .build();
                simulation = new SingleTMonteCarloSimulation(temperature, lattice, 30, maxSweeps, seed, mutualRnd, continueFromSave, realTimeEqTest, outputWriter);

            } catch (IOException e) {
                System.err.println("Error writing to file. ");
                e.printStackTrace();

            }
        }else{
            // multiple temperatures
            double[] T = receiveTemperatureSchedule(tempScheduleFileName);

            long[] seeds = LongStream.rangeClosed(seed+1, seed+T.length).toArray();
            MersenneTwister[] rnd = new MersenneTwister[T.length];
            SingleTMonteCarloSimulation[] subSimulations = new SingleTMonteCarloSimulation[T.length];
            File fSaveState;

            try {
                for (int i=0;i<T.length;i++){
                    Lattice lattice = new Lattice(Lx, Lz, extBx, suppressInternalTransFields, intTable, exchangeIntTable, nnArray, energyTable, momentTable, new ObservableExtractor(k_cos_table, k_sin_table));
                    rnd[i] = new MersenneTwister(seeds[i]);
                    FileWriter out = new FileWriter("analysis" + File.separator + folderName + File.separator + "table_" + Lx + "_" + Lz + "_" + extBx + "_" + temperature + "_" + suppressInternalTransFields + ".txt", successReadFromFile)
                    OutputWriter outputWriter = new OutputWriter.Builder(verboseOutput, folderName, obsPrintSweepNum, out)
                            .setPrintOutput(printOutput)
                            .setPrintProgress(printProgress)
                            .setBufferSize(bufferSize)
                            .build();

                    makeDir("states" + File.separator, folderName);
                    fSaveState = new File("states" + File.separator + folderName + File.separator + "save_state_" + Lx + "_" + Lz + "_" + dilution + "_" + h + "_" + extBx + "_" + suppressInternalTransFields + "_t.txt");
                    FileInputStream fis = null;
                    boolean successReadFromFile=false;

                    if (fSaveState.exists() && continueFromSave){
                        try{
                            // Open FileInputStream to the file
                            fis = new FileInputStream(fSaveState);
                            // Deserialize and cast into String
                            state = (ProgramState) SerializationUtils.deserialize(fis);
                            successReadFromFile=true;
                        }catch(Exception e) {
                            // for any problem reading previous state, just continue from start
                            successReadFromFile=false;
                        }finally{
                            if (fis!=null) fis.close();
                        }
                    }

                    subSimulations[i] = new SingleTMonteCarloSimulation(temperature, lattice, 30, maxSweeps, seeds[i], rnd[i], continueFromSave, realTimeEqTest, outputWriter);

                }





                if (parallelMode=='s'){
                    // serial mode
                    simulation = new ParallelMonteCarloSimulation(T, subSimulations, maxSweeps, seed, mutualRnd, continueFromSave, realTimeEqTest, parallelTemperingOff);
                }else if (parallelMode=='p'){
                    // serial mode
                    simulation = new SerialMonteCarloSimulation(T, subSimulations, maxSweeps, seed, mutualRnd, continueFromSave, realTimeEqTest, parallelTemperingOff);
                }

            }catch (IOException e) {System.err.println("error writing to file"); }
            finally {
                // close all outputs
                if (simulation!=null){
                    try{
                        simulation.close();
                    }catch (IOException e){
                        System.err.println("error closing simulation and files");
                    }
                }
            }
        }

    }

}