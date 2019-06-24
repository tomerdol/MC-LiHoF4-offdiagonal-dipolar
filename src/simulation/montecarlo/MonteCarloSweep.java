package simulation.montecarlo;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.MersenneTwister;
import simulation.mmsolve.ConvergenceException;

import java.io.BufferedWriter;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.Formatter;

import static simulation.montecarlo.MonteCarloMetropolis.*;

public class MonteCarloSweep implements Runnable{
    int t;                  // temperature index
    double temperature;
    Lattice lattice// lattice of singleSpin
    boolean printProgress;  // flag indicating whether to print the progress to the console
    int totalNumOfTemperatures;
    MersenneTwister rnd;
    int maxIter;
    double alpha;
    BufferedWriter outProblematicConfigs;
    long sweepBatch;
    double currentEnergy;
    boolean writeToFile;

    public void run(){
        // a monte carlo sweep consists of ~N steps
        for (int step = 0; step < lattice.getN(); step++) {

            if (printProgress) {
                int perc = (int) 10.0 * step /  lattice.getN();
                char[] full = new char[perc];
                char[] empty = new char[10 - perc];
                Arrays.fill(full, '#');
                Arrays.fill(empty, ' ');
                String strFull = new String(full);
                String strEmpty = new String(empty);
                System.out.print("[" + (t + 1) + "/" + totalNumOfTemperatures + "] [" + strFull + strEmpty + "]" + step + "/" +  lattice.length + "     \r");
            }

            Lattice tempLattice = new Lattice(lattice);
            Double deltaEnergy;
            deltaEnergyAndLattice energyAndLattice;
            try {
                energyAndLattice = metropolisStep(lattice, temperature, rnd, maxIter, alpha);

                if (energyAndLattice!=null){
                    deltaEnergy=energyAndLattice.getDeltaEnergy();
                     lattice=energyAndLattice.getLattice();
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
                 lattice=tempLattice;
            }else{
                tempLattice=null;
                currentEnergy += deltaEnergy.doubleValue();
            }

        }

        double m = lattice.getMagnetization();
        double[] temp = lattice.getMagneticFields();
        double[] tempSpinSizes = lattice.getSpinSizes();
        double mk2 = lattice.getMK2();    // m(k)^2, used later for correlation length calculation

        if (writeToFile) {
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
                    System.out.println("System equilibrated! sweeps="+sweeps+". Date&Time: "+ LocalDateTime.now());
                    equilibrationDeclared=true;
                } else if (t==T.length-1 && equilibrationDeclared && !isAllTrue(equilibratedT)){
                    // equilibration has previously been declared and now it turns out to have been a mistake
                    System.out.println("System not in equilibrium! sweeps="+sweeps+". Date&Time: "+LocalDateTime.now());
                    equilibrationDeclared=false;
                }
            }

            addToAvg(binAvg[t],new double[]{Math.abs(m), m , m*m, currentEnergy[t] , temp[0] , temp[1] , temp[2] , temp[3] , temp[4] , temp[5] , temp[6] , temp[7] , tempSpinSizes[0] , tempSpinSizes[1] , mk2});

        }
    }


    /**
     * Performs a Monte Carlo metropolis step
     * @param lattice - spin lattice
     * @param T - Temperature
     * @param rnd - PRNG
     * @return Array containing the (new) energy due to local transversal fields and the energy difference in the longitudinal energy compared to the previous configuration
     */
    public static deltaEnergyAndLattice metropolisStep(Lattice lattice, double T, MersenneTwister rnd, int maxIter, double alpha) throws ConvergenceException
    {

        double initLongEnergy = lattice.getEnergy();

        // choose random spin
        int flippedSpin = rnd.nextInt(lattice.getN());

        lattice.flipSpin(maxIter, Constants.tol, flippedSpin, alpha);

        double deltaEnergy = lattice.getEnergy() - initLongEnergy;

        double prob;
        if (deltaEnergy <= 0) {
            prob = 1;
        }
        else {
            prob = Math.exp(-(deltaEnergy )/T);
        }
        // accept new state with calculated probability
        double r = rnd.nextDouble();

        if (r > prob) {   // change is denied
			/*
            flippedSpin.flipSpin();  // the change is reversed. no need to update other longitudinal fields

            updateLocalFieldsAfterFlip(lattice, flippedSpin, intTable[2]);
            updateLocalTransFieldsAfterFlip(lattice, flippedSpin, intTable);

            updateAllMagneticMoments(lattice, intTable, momentTable, maxIter, extBx);
			*/

			/*
			for (int i=0;i<lattice.length;i++) {
				lattice[i]=tempLattice[i];
			}
			*/

            return null;

        } else {    // change is accepted
            // updateLocalFieldsAfterFlip(lattice, flippedSpin, intTable[2]);    // update other *longitudinal* (that's the reason for the [2]) fields

        }

        return new deltaEnergyAndLattice(lattice, deltaEnergy);
        //return deltaEnergy;
    }

}
