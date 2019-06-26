package simulation.montecarlo;

import java.io.*;
import java.time.LocalDateTime;
import java.util.*;
import java.util.stream.IntStream;

import simulation.mmsolve.ConvergenceException;
import simulation.mmsolve.MagneticMomentsSolveIter;
import simulation.mmsolve.fi_xi;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.SerializationUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.collections4.queue.CircularFifoQueue;

public class MonteCarloMetropolis {
	static int failedSteps=0, successfulSteps=0;
	static double avgSuccessTime=0, avgFailTime=0;
	static long successCount=0, failCount=0;



	public static void printArr(singleSpin[] arr, int Lz, int Lx){
		for (int i=0; i<arr.length; i++){
			double B = Math.sqrt(arr[i].getLocalBx()*arr[i].getLocalBx() + arr[i].getLocalBy()*arr[i].getLocalBy());
			System.out.print(arr[i].getX(Lz,Lx) + "," + arr[i].getY(Lz,Lx) + "," + arr[i].getZ(Lz,Lx) + "," + arr[i].getN() + "," + arr[i].getSpin() + "," + B + "\n");

		}
	}





	public static double updateFieldsAndCheckMagneticMomentConvergence(singleSpin[] arr, double[][][] intTable, FieldTable momentTable, double extBx, boolean suppressInternalTransFields, ConvergenceException e){
		singleSpin[] tempLattice = Lattice.copyLattice(arr);
		double beforeUpdate = magneticMomentConvergence(tempLattice, momentTable, extBx);
		updateAllLocalTransFields(tempLattice, intTable, extBx, suppressInternalTransFields);
		updateAllLocalFields(tempLattice, intTable[2]);
		double afterUpdate=magneticMomentConvergence(tempLattice, momentTable, extBx);
		//System.out.println(Math.abs(beforeUpdate-afterUpdate) +"  "+ verifyAllLocalFields(arr, intTable, extBx));
		if (Math.abs(beforeUpdate-afterUpdate)>1.0e-10){
			System.err.println("fields are not up to date! difference is " + (afterUpdate-beforeUpdate));
			System.err.println("fields state after (in arr): " + verifyAllLocalFields(arr, intTable, extBx));
			System.err.println("fields state after (in tempLattice): " + verifyAllLocalFields(tempLattice, intTable, extBx));
			if (e!=null) {
				System.err.println("e is not null:");
				e.printStackTrace();
			}
			System.exit(1);
		}
		return magneticMomentConvergence(tempLattice, momentTable, extBx);
	}


	public static void printProblematicConfig(singleSpin[] arr, int flippedSpin, BufferedWriter out){
		synchronized (out) {
			try {
				out.write("#" + flippedSpin + System.lineSeparator());
				for (int i = 0; i < arr.length; i++) {
					out.write(arr[i].getSpin() + "," + arr[i].getSpinSize() + System.lineSeparator());
				}

			} catch (IOException e) {
				System.err.println("error writing problematic config to file" + e.toString());
			}
		}
	}





	public static void main(String[] args) throws RuntimeException, IOException {
}



		// main program

//		MersenneTwister rnd = new MersenneTwister();

//		singleSpin[][] lattice = new singleSpin[T.length][];
		double[] currentEnergy = new double[T.length];
//		long sweeps;


//		if (!receivedSeed) {
//			seed = initializeRNG(rnd, fileNumber, taskID);
//		} else {
//			if (seed!=0) {
//				rnd.setSeed(seed);
//			}
//			else {
//				System.err.println("PRNG seed was not initialized for some reason!");
//				System.exit(1);
//			}
//		}



		// write data to file for later additional processing
//		FileWriter[] out = new FileWriter[T.length];
		// write states for checkpointing:
		ProgramState state;
		File fSaveState;


		try{
			// ********************* initialize output of checkpoint ************************
			// initialize program state
//			state = new ProgramState(null,null, 0, 0, null, null, null, null, null ,null, null);
//			makeDir("states" + File.separator, folderName);
//			fSaveState = new File("states" + File.separator + folderName + File.separator + "save_state_" + Lx + "_" + Lz + "_" + dilution + "_" + h + "_" + extBx + "_" + suppressInternalTransFields + "_t.txt");
//			FileInputStream fis = null;
//			boolean successReadFromFile=false;
//
//			if (fSaveState.exists() && continueFromSave){
//				try{
//					// Open FileInputStream to the file
//					fis = new FileInputStream(fSaveState);
//					// Deserialize and cast into String
//					state = (ProgramState) SerializationUtils.deserialize(fis);
//					successReadFromFile=true;
//				}catch(Exception e) {
//					// for any problem reading previous state, just continue from start
//					successReadFromFile=false;
//				}finally{
//					if (fis!=null) fis.close();
//				}
//			}

			// *******************************************************************************


			long[] currentBinCount;
			double[][] binAvg;
			int[] acceptanceRateCount, acceptanceRateSum;
			// array that hold the averages and sd's of the last 3 bins of 4 observables: m^2, E, mk2, m:
			ArrayList<ArrayList<CircularFifoQueue<Pair<Double,Double>>>> equilibratingObservables;

//			for (int t=0;t<T.length;t++) {
//				double temperature = T[t];
//				if (t<T.length-1) {	// last temperature is only for parallel tempering and not written to file
//
//					// ******************** initialize output of results ****************************
//					makeDir("analysis" + File.separator, folderName);
//					out[t] = new FileWriter("analysis" + File.separator + folderName + File.separator + "table_" + Lx + "_" + Lz + "_" + dilution + "_" + h + "_" + extBx + "_" + temperature + "_" + suppressInternalTransFields + "_t.txt",successReadFromFile);
//					// print some information to the begining of the file:
//					print("#" + LocalDateTime.now(), out[t], true, printOutput);
//					print("#temperature_schedule: "+Arrays.toString(T), out[t], true, printOutput);
//					print("#T="+t + ":" + temperature, out[t], true, printOutput);
//
//					//print constants:
//					print("#" + Constants.constantsToString(), out[t], true, printOutput);
//					print(String.format("# Lx=%s, Ly=%s, Lz=%s, extBx=%s, maxSweeps=%s, suppressInternalTransFields=%s, " +
//							"continueFromSave=%s, maxIter=%s, bufferSize=%s, tempScheduleFileName=%s, parallelTemperingOff=%s, " +
//							"saveState=%s, folderName=%s, alpha=%s, verboseOutput=%s ",Lx,Lx,Lz,extBx, maxSweeps,suppressInternalTransFields,
//							continueFromSave, maxIter, bufferSize, tempScheduleFileName, parallelTemperingOff, saveState, folderName,
//							alpha, verboseOutput),out[t],true,printOutput);
//					print("#" + Constants.locationsToString(), out[t], true, printOutput);

					// ******************************************************************************
				}
			}

//			if (successReadFromFile){
//				currentEnergy=state.getCurrentEnergy();
//				sweeps=state.getSweep();
//				lattice = state.getArr();
//				rnd=state.getRnd();
//				seed=state.getSeed();
//				currentBinCount=state.getCurrentBinCount();
//				binAvg=state.getBinAvg();
//				acceptanceRateCount=state.getAcceptanceRateCount();
//				acceptanceRateSum=state.getAcceptanceRateSum();
//				equilibratingObservables=state.getEquilibratingObservables();
//
//				for (int t=0;t<T.length-1;t++) {
//					print("#seed=" + seed, out[t], true, printOutput);
//					print("# successfully read saved state", out[t], true, printOutput);
//				}

//			}else{
//				sweeps=0;
/*
				currentBinCount=new long[T.length-1];
				binAvg = new double[T.length][30];
				acceptanceRateCount = new int[T.length-1];
				acceptanceRateSum = new int[T.length-1];
*/
				//equilibratingObservables = generateEquilibrationQueues(3, T.length-1, 4);


//				for (int t=0;t<T.length;t++) {
//					if (t < T.length - 1) {
////						print("#seed=" + seed, out[t], true, printOutput);
////						print("# unsuccessful reading saved state... Starting new state.", out[t], true, printOutput);
//						if (verboseOutput)
//							print(makeTableHeader(colWidths,colNames), out[t], true, printOutput);
//						else
//							//print("      binN            <|M|>           <|M|2>              <M>             <M2>             <M2>            <M22>              <E>             <E2>         <meanBx>        <meanBx2>          <stdBx>         <stdBx2>         <meanBy>        <meanBy2>          <stdBy>" +
//							//		"         <stdBy2>         <meanBz>        <meanBz2>          <stdBz>         <stdBz2>      <maxBtrans>     <maxBtrans2>       <maxBlong>      <maxBlong2>   <meanSpinSize>  <meanSpinSize2>    <stdSpinSize>   <stdSpinSize2>            <mk2>           <mk22>      swap", out[t], true, printOutput);
//							print(makeTableHeader(colWidths,colNames), out[t], true, printOutput);

//					}

//					boolean convergedConfig=false;
//					int index=0;
//					while (!convergedConfig) {
//
//						lattice[t] = randomizeConfig(lattice[t], rnd);
//						// set initial spin sizes:
//						for (int i=0;i<lattice[t].length;i++) lattice[t][i].setSpinSize(Constants.spinSize*lattice[t][i].getSpin());
//
//						updateAllLocalFields(lattice[t], intTable[2]);    // only longitudinal fields
//						updateAllLocalTransFields(lattice[t], intTable, extBx, suppressInternalTransFields);
//						MagneticMomentsSolveIter.updateAllMagneticMoments(lattice[t], intTable, momentTable, 2*maxIter, extBx, Constants.tol, alpha, suppressInternalTransFields, true);
//
//						if (MagneticMomentsSolveIter.magneticMomentConvergence(lattice[t], momentTable, extBx) > Constants.tol){
//							convergedConfig=false;
//						}else{
//							convergedConfig=true;
//						}
//
//						if (!convergedConfig) System.out.println("initial setup convergence index (T [num]="+T[t]+" ["+t+']'+"): " + ++index);
//					}
					/*
					updateAllLocalFields(lattice[t], intTable[2]);    // only longitudinal fields
					updateAllLocalTransFields(lattice[t], intTable, extBx, suppressInternalTransFields);
					updateAllMagneticMoments(lattice[t], intTable, momentTable, maxIter, extBx, Constants.tol, 1, suppressInternalTransFields);
					*/
//					currentEnergy[t] = calcEnergy(lattice[t], energyTable, exchangeIntTable);

//				}

//				if (saveState) {
//					state.updateAll(currentEnergy, seed, sweeps, rnd, lattice, currentBinCount, binAvg, acceptanceRateCount, acceptanceRateSum, equilibratingObservables);
//					FileOutputStream fos = new FileOutputStream(fSaveState, false);    // false to overwrite previous state
//					SerializationUtils.serialize(state, fos);
//					fos.close();
//				}
//			}


			// initialize string output buffers
			StringBuilder[] outputBuffer = new StringBuilder[T.length - 1];
			for (int t = 0; t < T.length - 1; t++)
				outputBuffer[t] = new StringBuilder(bufferSize);

			// initialize output buffer for problematic configurations
			BufferedWriter outProblematicConfigs = new BufferedWriter(new FileWriter("p_configs" + File.separator + "problematic_"+lattice[0].length+"_"+extBx,true));

			// initialize equilibration array
            boolean[] equilibratedT = new boolean[T.length-1];
            for (int t=0;t<equilibratedT.length;t++) equilibratedT[t] = checkEquilibration(equilibratingObservables.get(t));
            boolean equilibrationDeclared = false;

			long startTime=System.currentTimeMillis();

			long numOfBatches = maxSweeps;	// try parallel tempering replace after each monte carlo sweep
			//long numOfBatches = 1;	// don't perform parallel tempering
			for (long sweepBatch = 0; sweepBatch<numOfBatches; sweepBatch++) {
				for (long sweep = 0; sweep < maxSweeps/numOfBatches; sweep++) {

					for (int t = 0; t < T.length; t++) {

						double temperature = T[t];
						// a monte carlo sweep consists of ~N steps
						for (int step = 0; step < lattice[t].length; step++) {


							if (printProgress) {
								int perc = (int) 10.0 * step / lattice[t].length;
								char[] full = new char[perc];
								char[] empty = new char[10 - perc];
								Arrays.fill(full, '#');
								Arrays.fill(empty, ' ');
								String strFull = new String(full);
								String strEmpty = new String(empty);
								System.out.print("[" + (t + 1) + "/" + T.length + "] [" + strFull + strEmpty + "]" + step + "/" + lattice[t].length + "     \r");
							}

							singleSpin[] tempLattice = Lattice.copyLattice(lattice[t]);
							Double deltaEnergy;
							deltaEnergyAndLattice energyAndLattice;
							try {
								energyAndLattice = MonteCarloSweep.metropolisStep(lattice, intTable, exchangeIntTable, temperature, rnd, energyTable, momentTable, suppressInternalTransFields, extBx, maxIter, method, nnArray, alpha, tolCheckMult);

                                if (energyAndLattice!=null){
                                    deltaEnergy=energyAndLattice.getDeltaEnergy();
                                    lattice[t]=energyAndLattice.getLattice();
                                }else{
                                    deltaEnergy=null;
                                }
                                successfulSteps++;

								//System.out.println("step " + step + " successful");
							} catch (ConvergenceException e){
								printProblematicConfig(tempLattice, e.getFlippedSpin(), outProblematicConfigs);
								deltaEnergy=null;
								failedSteps++;
								System.err.println("There was an error converging, had to abort metropolis step ("+step+ ", sweep="+sweepBatch+" T="+temperature+") and try another one. \n" + e.getMessage());
							} catch (IndexOutOfBoundsException e){
								deltaEnergy=null;
								failedSteps++;
								System.err.println("There was an error while calculating the derivative, had to abort metropolis step ("+step+ ", sweep="+sweepBatch+" T="+temperature+") and try another one.\n" + e.toString());
							}

							if (deltaEnergy == null) {
								// reverse change
								lattice[t]=tempLattice;
							}else{
								tempLattice=null;
								currentEnergy[t] += deltaEnergy.doubleValue();
							}

						}

						double m = calcMagnetization(lattice[t]);
						double[] temp = meanField(lattice[t]);
						double[] tempSpinSizes = calcSpinSizes(lattice[t]);
						double mk2 = calc_mk2(lattice[t], k_cos_table, k_sin_table);    // m(k)^2, used later for correlation length calculation

						if (t < T.length - 1) {
							Formatter formatter = new Formatter(outputBuffer[t]);
							if (verboseOutput) {
								//outputBuffer[t].append(sweeps + "\t" + m + "\t" + currentEnergy[t] + "\t" + temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\t" + temp[4] + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + tempSpinSizes[0] + "\t" + tempSpinSizes[1] + "\t" + mk2);
								//formatter.format("% 10d % 16.9g % 18.13g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g % 16.9g",sweeps, m ,currentEnergy[t] ,temp[0] ,temp[1] ,temp[2] ,temp[3] ,temp[4] ,temp[5] ,temp[6] ,temp[7] ,tempSpinSizes[0] ,tempSpinSizes[1] ,mk2);
								formatter.format(makeTableRowFormat(colWidths, new char[]{'d','g','g','g','g','g','g','g','g','g','g','g','g','g','g'}),sweeps, m ,currentEnergy[t] ,temp[0] ,temp[1] ,temp[2] ,temp[3] ,temp[4] ,temp[5] ,temp[6] ,temp[7] ,tempSpinSizes[0] ,tempSpinSizes[1] ,mk2);
                                if (sweep < maxSweeps / numOfBatches - 1)
									outputBuffer[t].append(System.lineSeparator());
							}
							if (sweeps>0) currentBinCount[t]++;
							//System.out.println(sweeps+ " , " + sweeps/2 + " , " + currentBinCount + " , " + (sweeps&1));

							if (currentBinCount[t]==sweeps/2 + 1 && (sweeps&1)==1){	// last step of current bin
								// restart sums and counters for next bin
								double acceptanceRateForBin;	// acceptance rate for this bin & temperature

								if (acceptanceRateCount[t]==0) {
									acceptanceRateForBin = -1;
								}else {
									acceptanceRateForBin = ((double) acceptanceRateSum[t]) / acceptanceRateCount[t];
								}

								// print data from previous bin and
								if (!verboseOutput){
									formatter.format("% "+colWidths[0]+"d %s % "+(colWidths[colWidths.length-1]-2)+".2f%n", currentBinCount[t], avgArrToString(binAvg[t], currentBinCount[t],colWidths,1) , acceptanceRateForBin);
								}

								// add bin averages to queue. only 3 bins are kept at any time.
								equilibratingObservables.get(t).get(0).add(Pair.of(binAvg[t][4]/currentBinCount[t],getStandardError(binAvg[t][4],binAvg[t][5],currentBinCount[t])));		// m^2
								equilibratingObservables.get(t).get(1).add(Pair.of(binAvg[t][6]/currentBinCount[t],getStandardError(binAvg[t][6],binAvg[t][7],currentBinCount[t])));		// E
								equilibratingObservables.get(t).get(2).add(Pair.of(binAvg[t][28]/currentBinCount[t],getStandardError(binAvg[t][28],binAvg[t][29],currentBinCount[t])));	// mk2
								equilibratingObservables.get(t).get(3).add(Pair.of(binAvg[t][0]/currentBinCount[t],getStandardError(binAvg[t][0],binAvg[t][1],currentBinCount[t])));		// |m|
                                equilibratedT[t] = checkEquilibration(equilibratingObservables.get(t));

								//System.out.println(sweeps + ") printed average of bin size "+currentBinCount + " starting from " + binStart);
								for (int avgIndex=0;avgIndex<binAvg[t].length;avgIndex++) binAvg[t][avgIndex]=0;

								acceptanceRateCount[t]=0;
								acceptanceRateSum[t]=0;

								//binStart=sweeps;

								currentBinCount[t]=0;

                                if (t==T.length-1 && !equilibrationDeclared && isAllTrue(equilibratedT)){
                                    // declare that equilibration has been reached
                                    System.out.println("System equilibrated! sweeps="+sweeps+". Date&Time: "+LocalDateTime.now());
                                    equilibrationDeclared=true;

                                    // TODO change output type from bin to verbose

                                } else if (t==T.length-1 && equilibrationDeclared && !isAllTrue(equilibratedT)){
                                    // equilibration has previously been declared and now it turns out to have been a mistake
                                    System.out.println("System not in equilibrium! sweeps="+sweeps+". Date&Time: "+LocalDateTime.now());
                                    equilibrationDeclared=false;
                                }
							}

							addToAvg(binAvg[t],new double[]{Math.abs(m), m , m*m, currentEnergy[t] , temp[0] , temp[1] , temp[2] , temp[3] , temp[4] , temp[5] , temp[6] , temp[7] , tempSpinSizes[0] , tempSpinSizes[1] , mk2});

						}
					}   // end of temperature loop

					sweeps++;

				}
				// parallel tempering scheme
				for (int t = 0; t < T.length - 1; t++) {
					double thisEnergy = currentEnergy[t];
					double nextEnergy = currentEnergy[t + 1];
					double delta = (1/T[t] - 1/T[t + 1])*(thisEnergy - nextEnergy);
					if (!parallelTemperingOff && rnd.nextDouble() < Math.exp(delta)) {
						// ****************swap energies****************
						double tempEnergy = currentEnergy[t];
						currentEnergy[t] = currentEnergy[t + 1];
						currentEnergy[t + 1] = tempEnergy;

						// *********************************************

						// ****************swap lattices****************
						singleSpin[] tempLattice = lattice[t];
						lattice[t] = lattice[t + 1];
						lattice[t + 1] = tempLattice;
						// *********************************************

						if (verboseOutput) {
							// write swap acceptance rate
							outputBuffer[t].append(new String(new char[colWidths[colWidths.length-1]-2]).replace('\0', ' ') + "1");
							outputBuffer[t].append(System.lineSeparator());
						}else{
							// sum swap acceptance for bin average
							acceptanceRateCount[t]++;
							acceptanceRateSum[t]++;
						}
					} else {
						if (verboseOutput) {
							// swap acceptance rate
							outputBuffer[t].append(new String(new char[colWidths[colWidths.length-1]-2]).replace('\0', ' ') + "0");
							outputBuffer[t].append(System.lineSeparator());
						}else{
							// sum swap acceptance for bin average
							acceptanceRateCount[t]++;
						}
					}
				}


				if (sweeps%obsPrintSweepNum==0 || sweepBatch==numOfBatches-1) {	// every obsPrintSweepNum sweeps or at the last one
					// save state & flush output
					for (int t = 0; t < T.length - 1; t++) {
						print(outputBuffer[t].toString(), out[t], false, printOutput);
						out[t].flush();
						if (outputBuffer[t].length()>bufferSize) System.err.println("temperature #"+t+" exceeded its buffer with char count of "+outputBuffer[t].length()+'/'+bufferSize);
						outputBuffer[t]=null;
						outputBuffer[t]=new StringBuilder(bufferSize);
					}

					if (saveState) {
						state.updateAll(currentEnergy, seed, sweeps, rnd, lattice, currentBinCount, binAvg, acceptanceRateCount, acceptanceRateSum, equilibratingObservables);
						FileOutputStream fos = new FileOutputStream(fSaveState, false);    // false to overwrite previous state
						SerializationUtils.serialize(state, fos);
						fos.close();
					}

//					System.out.println(LocalDateTime.now());

				}
				if (printProgress) System.out.println(String.format("%.2f",100*(double)sweepBatch/numOfBatches) + "% complete                ");

			}

			System.out.println("Success %: "+successfulSteps+"/"+(failedSteps+successfulSteps)+" - "+100*((double)successfulSteps/(successfulSteps+failedSteps))+"%");
			System.out.println("total time: " + (System.currentTimeMillis()-startTime));
			System.out.println();
			long totalCount=failCount+successCount;
			System.out.println("average success time: "+(avgSuccessTime/successCount));
			System.out.println("average fail time: "+(avgFailTime/failCount));
			System.out.println("success rate: " + successCount + "/" + totalCount + "=" + (1.0*successCount)/totalCount);
			System.out.println("fail rate: " + failCount + "/" + totalCount + "=" + (1.0*failCount)/totalCount);
			System.out.println((avgSuccessTime/successCount) + "\t" + (1.0*successCount)/totalCount + "\t" + (avgFailTime/failCount) + "\t" + (1.0*failCount)/totalCount);


			// close outProblematicConfigs
			try {
				outProblematicConfigs.close();
			} catch (IOException e) {
				System.out.println("error closing problematic configs file");
			}
		}


//		catch (IOException e) {System.err.println("error writing to file"); }
//		finally {
//			// close all outputs
//			for (int t=0;t<T.length-1;t++) {
//				if (out[t] != null) {
//					try {
//						out[t].close();
//					} catch (IOException e) {
//						System.out.println("error closing file");
//					}
//
//				}
//			}
//		}






	}

	/*
	@Deprecated
	public static void solveSelfConsistentCalc(singleSpin[] lattice, double[][][] intTable, FieldTable momentTable,
											   int maxIter, double extBx, double tol, boolean suppressInternalTransFields) throws ConvergenceException
	{
		//updateAllLocalFields(lattice, intTable[2]);    // only longitudinal fields
		//updateAllLocalTransFields(lattice, intTable, extBx, suppressInternalTransFields);

		// first try the regular algorithm (gauss-seidel like)
		updateAllMagneticMoments(lattice, intTable, momentTable, maxIter, extBx, tol, 0.97, suppressInternalTransFields);

		//System.out.println(magneticMomentConvergence(lattice, momentTable, extBx));

		// if the regular method does not converge, try the broyden method:
		if (magneticMomentConvergence(lattice, momentTable, extBx) > tol) {
			String broydenVersionUsed="Broyden1";
			//System.out.println("convergence was " + magneticMomentConvergence(lattice, momentTable, extBx) + ". Running broyden");
			// make x vector
			double[] x = new double[lattice.length];
			for (int j=0;j<lattice.length;j++) x[j]=lattice[j].getSpinSize();// + (rnd.nextDouble()-0.5)*0.1;

			fi_xi funcBroyden = new func(intTable, momentTable, lattice, extBx, suppressInternalTransFields);
			try {
				x = simulation.mmsolve.broyden(funcBroyden, x);
			} catch(ConvergenceException e){
				broydenVersionUsed="Broyden2";
				funcBroyden = new func(intTable, momentTable, lattice, extBx, suppressInternalTransFields);
				//System.err.println(e.toString());
				//System.err.println("running broyden2");
				for (int j=0;j<lattice.length;j++) x[j]=lattice[j].getSpinSize();   // reinitialize x
				x = simulation.mmsolve.broyden2(new func(intTable, momentTable, lattice, extBx, suppressInternalTransFields), x);
				//System.out.println("broyden2 successful!");
			}
			// If we reach this point, simulation.mmsolve did not throw ConvergenceException so we assume Broyden has converged successfully


			for (int j=0;j<lattice.length;j++) lattice[j].setSpinSize(x[j]);

			updateAllLocalFields(lattice, intTable[2]);    // only longitudinal fields
			updateAllLocalTransFields(lattice, intTable, extBx, suppressInternalTransFields);


			if (magneticMomentConvergence(lattice, momentTable, extBx) > tol) { // not converged
				throw new ConvergenceException("Error: Self consistent calculation did not converge! ", broydenVersionUsed, funcBroyden, 0, magneticMomentConvergence(lattice, momentTable, extBx));
			}
		}

	}

	 */

}

class deltaEnergyAndLattice {
	private Lattice lattice;
	private double deltaEnergy;

	public deltaEnergyAndLattice(Lattice lattice, double deltaEnergy) {
		this.lattice = lattice;
		this.deltaEnergy = deltaEnergy;
	}

	public Lattice getLattice() {
		return lattice;
	}

	public void setLattice(Lattice lattice) {
		this.lattice = lattice;
	}

	public double getDeltaEnergy() {
		return deltaEnergy;
	}

	public void setDeltaEnergy(double deltaEnergy) {
		this.deltaEnergy = deltaEnergy;
	}
}


/**
 * Saves important variables and states so that the program can continue from where it left off in case of a crash.
 * @author Tomer
 *
 */
class ProgramState implements Serializable {	// TODO: add some sort of deserialize verification (Lx, Lz etc. match current values)
	private double[] currentEnergy;
	private long seed, sweep;
	private MersenneTwister rnd;
	private singleSpin[][] arr;
	private long[] currentBinCount;
	private double[][] binAvg;
	private int[] acceptanceRateCount;
	private int[] acceptanceRateSum;
	private ArrayList<ArrayList<CircularFifoQueue<Pair<Double,Double>>>> equilibratingObservables;

	/**
	 * Constructor for Program state
	 * @param currentTransEnergy - current transverse energy array
	 * @param currentLongEnergy - current longitudinal energy array
	 * @param seed - RNG seed
	 * @param sweep - sweep number
	 * @param rnd - random number generator (MersenneTwister)
	 * @param arr - spin lattice arrays
	 */
	public ProgramState(double[] currentTransEnergy, double[] currentLongEnergy, long seed, long sweep, MersenneTwister rnd, singleSpin[][] arr, long[] currentBinCount, double[][] binAvg,
						int[] acceptanceRateCount, int[] acceptanceRateSum, ArrayList<ArrayList<CircularFifoQueue<Pair<Double,Double>>>> equilibratingObservables) {
		this.currentEnergy = currentEnergy;
		this.seed = seed;
		this.sweep = sweep;
		this.rnd = rnd;
		this.arr = arr;
		this.acceptanceRateCount = acceptanceRateCount;
		this.acceptanceRateSum = acceptanceRateSum;
		this.binAvg = binAvg;
		this.currentBinCount = currentBinCount;
		this.equilibratingObservables = equilibratingObservables;
	}


	/**
	 * update program state (to be used before saving)
	 * @param currentEnergy - current energy array
	 * @param seed - RNG seed
	 * @param sweep - sweep number
	 * @param rnd - random number generator (MersenneTwister)
	 * @param arr - spin lattice arrays
	 */
	public void updateAll(double[] currentEnergy, long seed, long sweep, MersenneTwister rnd, singleSpin[][] arr, long[] currentBinCount, double[][] binAvg,
						  int[] acceptanceRateCount, int[] acceptanceRateSum, ArrayList<ArrayList<CircularFifoQueue<Pair<Double,Double>>>> equilibratingObservables) {
		this.currentEnergy = currentEnergy;
		this.seed = seed;
		this.sweep = sweep;
		this.rnd = rnd;
		this.arr = arr;
		this.acceptanceRateCount = acceptanceRateCount;
		this.acceptanceRateSum = acceptanceRateSum;
		this.binAvg = binAvg;
		this.currentBinCount = currentBinCount;
		this.equilibratingObservables = equilibratingObservables;
	}

	public double[] getCurrentEnergy() {
		return currentEnergy;
	}

	public void setCurrentEnergy(double[] currentTransEnergy) {
		this.currentEnergy = currentEnergy;
	}


	public long getSeed() {
		return seed;
	}

	public void setSeed(long seed) {
		this.seed = seed;
	}

	public long getSweep() {
		return sweep;
	}

	public void setSweep(long sweep) {
		this.sweep = sweep;
	}

	public MersenneTwister getRnd() {
		return rnd;
	}

	public void setRnd(MersenneTwister rnd) {
		this.rnd = rnd;
	}

	public singleSpin[][] getArr() {
		return arr;
	}

	public void setArr(singleSpin[][] arr) {
		this.arr = arr;
	}

	public long[] getCurrentBinCount() {
		return currentBinCount;
	}

	public void setCurrentBinCount(long[] currentBinCount) {
		this.currentBinCount = currentBinCount;
	}

	public double[][] getBinAvg() {
		return binAvg;
	}

	public void setBinAvg(double[][] binAvg) {
		this.binAvg = binAvg;
	}

	public int[] getAcceptanceRateCount() {
		return acceptanceRateCount;
	}

	public void setAcceptanceRateCount(int[] acceptanceRateCount) {
		this.acceptanceRateCount = acceptanceRateCount;
	}

	public int[] getAcceptanceRateSum() {
		return acceptanceRateSum;
	}

	public void setAcceptanceRateSum(int[] acceptanceRateSum) {
		this.acceptanceRateSum = acceptanceRateSum;
	}

    public ArrayList<ArrayList<CircularFifoQueue<Pair<Double, Double>>>> getEquilibratingObservables() {
        return equilibratingObservables;
    }

    public void setEquilibratingObservables(ArrayList<ArrayList<CircularFifoQueue<Pair<Double, Double>>>> equilibratingObservables) {
        this.equilibratingObservables = equilibratingObservables;
    }
}
