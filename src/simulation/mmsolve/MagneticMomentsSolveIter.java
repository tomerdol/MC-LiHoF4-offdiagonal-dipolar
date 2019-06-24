package simulation.mmsolve;

import simulation.montecarlo.*;

public class MagneticMomentsSolveIter {
/*
    public static void updateAllLocalTransFields(singleSpin[] arr, double[][][] intTable, double extBx, boolean suppressInternalTransFields, int flipSpin){
        int i;
        for (i=0;i<arr.length;i++) {
            int j;
            if (!suppressInternalTransFields) {
                if (arr[i].getSpin() != 0) {
                    double Bx = 0, By = 0;
                    // add up all contributions to transversal field from other spins at spin i
                    for (j = 0; j < arr.length; j++) {
                        if (arr[j].getSpin() != 0 && i != j && j!=flipSpin) {
                            // calculate interaction
                            Bx += (arr[j].getSpinSize()) * intTable[0][arr[i].getN()][arr[j].getN()];
                            By += (arr[j].getSpinSize()) * intTable[1][arr[i].getN()][arr[j].getN()];
                        }
                    }

                    arr[i].setLocalBx(Bx + extBx);
                    arr[i].setLocalBy(By);
                }
            } else {
                arr[i].setLocalBx(extBx);
                arr[i].setLocalBy(0.0);
            }

        }
    }

    public static void updateAllLocalTransFields(singleSpin[] arr, double[][][] intTable, double extBx, boolean suppressInternalTransFields){
        updateAllLocalTransFields(arr, intTable, extBx, suppressInternalTransFields, -1);
    }
    */

    /*
    @Deprecated
    public static singleSpin[] homotopicSolve2(singleSpin[] lattice, double[][][] intTable, FieldTable momentTable,
                                               int maxIter, double extBx, double tol, boolean suppressInternalTransFields, int flipSpin, int[][] nnArray, boolean reverse, double alpha, boolean alphaEscalation) throws ConvergenceException
    {
        int[] bfsOrder=null;
        if (nnArray!=null) {
            bfsOrder = MagneticMomentsSolveIter.orderBFS(lattice.length, flipSpin, nnArray);
        }

        if (nnArray!=null && reverse) {
            // reverse bfsOrder
            for (int i = 0; i < bfsOrder.length / 2; i++) {
                int temp = bfsOrder[i];
                bfsOrder[i] = bfsOrder[bfsOrder.length - i - 1];
                bfsOrder[bfsOrder.length - i - 1] = temp;
            }
        }

        removeSpecificSpinConstibution(lattice, flipSpin, lattice[flipSpin].getSpinSize(), intTable, suppressInternalTransFields);

        int numOfStepsLeft=0;
        double perc = 1.0, prevPerc=0.0;
        int i=0;

        while (numOfStepsLeft>=0) {
            i++;
            //System.out.print(" homotopy step: " + i++ + "(" + perc+"%) ");

            double prevStateMoment = momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), lattice[flipSpin].getSpin());
            double nextStateMoment = momentTable.getValue(lattice[flipSpin].getLocalBx(), lattice[flipSpin].getLocalBy(), lattice[flipSpin].getLocalBz(), (-1) * lattice[flipSpin].getSpin());
            double magneticMoment = perc * nextStateMoment + (1.0 - perc) * prevStateMoment;

            singleSpin[] tempLattice = Temp.copyLattice(lattice);

            //updateAllLocalFields(lattice, intTable[2], flipSpin);
            //updateAllLocalTransFields(lattice, intTable, extBx, suppressInternalTransFields, flipSpin);

            addMixedSpinField(lattice, intTable, flipSpin, magneticMoment, suppressInternalTransFields);

            if (numOfStepsLeft>0) {
                updateAllMagneticMoments(lattice, bfsOrder, intTable, momentTable, maxIter, extBx, tol, alpha, suppressInternalTransFields, flipSpin);
                //System.out.print(" - "+(magneticMomentConvergence(lattice, momentTable, extBx, flipSpin)<tol ? "success ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin)+"), " : "failure ("+magneticMomentConvergence(lattice, momentTable, extBx, flipSpin)+"), "));
                if (magneticMomentConvergence(lattice, momentTable, extBx, flipSpin) > tol) {
                    // backtrack
                    //maxIter *= 2;
                    maxIter = (int)(1.1*maxIter);
                    numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                    lattice = tempLattice;
                    //magneticMoment = lattice[flipSpin].getSpinSize();
                    perc = prevPerc;
                    //System.err.print("backtracking: " + perc);

                }else {
                    removeSpecificSpinConstibution(lattice, flipSpin, magneticMoment, intTable, suppressInternalTransFields);
                }

            } else {    // last step
                lattice[flipSpin].flipSpin();
                lattice[flipSpin].setSpinSize(magneticMoment);

                updateAllMagneticMoments(lattice, bfsOrder, intTable, momentTable, maxIter, extBx, tol, alpha, suppressInternalTransFields, -1);

                //System.out.print(" (last step) - "+(MonteCarloMetropolis.updateFieldsAndCheckMagneticMomentConvergence(lattice, intTable, momentTable, extBx, suppressInternalTransFields)<tol ? "success ("+MonteCarloMetropolis.updateFieldsAndCheckMagneticMomentConvergence(lattice, intTable, momentTable, extBx, suppressInternalTransFields)+"), " : "failure ("+MonteCarloMetropolis.updateFieldsAndCheckMagneticMomentConvergence(lattice, intTable, momentTable, extBx, suppressInternalTransFields)+"), "));

                if (magneticMomentConvergence(lattice, momentTable, extBx) > tol) {
                    // backtrack
                    lattice[flipSpin].flipSpin();   // flip back

                    //maxIter *= 2;
                    maxIter = (int)(1.1*maxIter);
                    numOfStepsLeft = (numOfStepsLeft + 1) * 2;
                    lattice = tempLattice;
                    //magneticMoment = prevMagneticMoment;
                    perc = prevPerc;
                    //System.err.print("backtracking: " + perc);

                }
            }

            prevPerc=perc;
            if (numOfStepsLeft > 0) perc = perc + (1.0 - perc)/numOfStepsLeft;
            numOfStepsLeft--;


            //if (i>30) System.out.println(i);
            // too many iterations, so this is not going anywhere (unless this was the last step):
            if ((i>150 || numOfStepsLeft>100) && numOfStepsLeft>=0){
                numOfStepsLeft=-1;  // this ends the loop
                addMixedSpinField(lattice, intTable, flipSpin, lattice[flipSpin].getSpinSize(), suppressInternalTransFields);   // this make the fields match the moments
                lattice[flipSpin].flipSpin();
                if (!alphaEscalation) {
                    //System.err.println("started broyden escalation");
                    try {
                        solveSelfConsistentCalc(lattice, intTable, momentTable, maxIter, extBx, tol, suppressInternalTransFields);
                    } catch (ConvergenceException e) {
                        // possibly add some info about convergence before Broyden was called
                        e.setFlippedSpin(flipSpin);
                        throw e;
                    }
                } else {
                    //System.err.println("started alpha escalation");
                    boolean success=false;
                    while (alpha>0.3 && !success) { // for initial alpha=1.0 this means 12 tries
                        updateAllMagneticMoments(lattice, bfsOrder, intTable, momentTable, maxIter, extBx, tol, alpha, suppressInternalTransFields, -1);
                        alpha=alpha/0.9;
                        maxIter=(int)(maxIter*1.2);
                        if (magneticMomentConvergence(lattice, momentTable, extBx) < tol) {
                            success=true;
                        }
                    }
                    if (!success){
                        ConvergenceException e = new ConvergenceException("alpha escalation failed. numOfStepsLeft="+numOfStepsLeft+" before escalation.", "Regular", null, i, magneticMomentConvergence(lattice, momentTable, extBx));
                        e.setFlippedSpin(flipSpin);
                        throw e;
                    }

                }
            }
        }

        totalStep++;
        if (i>maxHomotopicStep) maxHomotopicStep=i;
        homotopicStep+=i;

        return lattice;
    }

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

    @Deprecated
        public static void addMixedSpinField(singleSpin[] arr, double[][][] intTable, int flipSpin, double magneticMoment, boolean suppressInternalTransFields){
        int i;

        for (i=0;i<arr.length;i++) {
            if (arr[i].getSpin()!=0){
                if (i!=flipSpin){
                    arr[i].setLocalBz(arr[i].getLocalBz() + magneticMoment*intTable[2][arr[i].getN()][arr[flipSpin].getN()]);
                    if (!suppressInternalTransFields) {
                        arr[i].setLocalBx(arr[i].getLocalBx() + magneticMoment*intTable[0][arr[i].getN()][arr[flipSpin].getN()]);
                        arr[i].setLocalBy(arr[i].getLocalBy() + magneticMoment*intTable[1][arr[i].getN()][arr[flipSpin].getN()]);
                    }
                }
            }
        }
    }


    @Deprecated
    public static void removeSpecificSpinConstibution(singleSpin[] arr, int spinToRemove, double magneticMomentToRemove, double[][][] intTable, boolean suppressInternalTransFields){
        int j;
        if (arr[spinToRemove].getSpin()!=0) {
            // change local fields (longitudinal & transverse) at all other spins by
            // removing the contribution of i there
            for(j=0;j<arr.length;j++){
                if (arr[j].getSpin()!=0 && spinToRemove!=j){
                    // remove interaction with prevSpinSize and add interaction with current spinSize
                    if (!suppressInternalTransFields) {
                        arr[j].setLocalBx(arr[j].getLocalBx() - magneticMomentToRemove * intTable[0][arr[spinToRemove].getN()][arr[j].getN()]);
                        arr[j].setLocalBy(arr[j].getLocalBy() - magneticMomentToRemove * intTable[1][arr[spinToRemove].getN()][arr[j].getN()]);
                    }
                    arr[j].setLocalBz(arr[j].getLocalBz() - magneticMomentToRemove * intTable[2][arr[spinToRemove].getN()][arr[j].getN()]);

                }
            }
        }
    }


    public static void updateAllMagneticMoments(singleSpin[] arr, int[] sweepOrder, double[][][] intTable, FieldTable momentTable, int maxIter, double extBx, double tol, double alpha, boolean suppressInternalTransFields, int flipSpin) {

        boolean converged = false;
        int iter;
        int index;
        double sum = 0;    // for real time convergence check
        //System.out.println("start:");
        for (iter = 0; iter < maxIter && !converged; iter++) {

            sum = 0;    // for real-time convergence check
            for (int i = 0; i < arr.length; i++) {
                if (sweepOrder==null){
                    index=i;
                }else{
                    index=sweepOrder[i];
                }

                if (index != flipSpin) {
                    double prevSpinSize = arr[index].getSpinSize();
                    arr[index].setSpinSize(prevSpinSize * (1 - alpha) + alpha * momentTable.getValue(arr[index].getLocalBx(), arr[index].getLocalBy(), arr[index].getLocalBz(), arr[index].getSpin()));
                    MonteCarloMetropolis.updateFieldsAfterSpinSizeChange(arr, arr[index], intTable, prevSpinSize, suppressInternalTransFields);

                    sum += Math.abs(prevSpinSize - arr[index].getSpinSize());
                }
            }

            if (sum / arr.length < tol*0.1) {

                converged = true;
            }


            //System.out.println(iter + "\t" + sum / arr.length + "\t" + MonteCarloMetropolis.magneticMomentConvergence(arr, momentTable, extBx));

        }

        if (!converged && alpha > 0.85) {
            //System.out.println(sum/arr.length);
            //System.out.println("initiating underrelaxed gauss-seidel");
            //System.out.println();
            updateAllMagneticMoments(arr, sweepOrder, intTable, momentTable, 2 * maxIter, extBx, tol, 0.5, suppressInternalTransFields, flipSpin);
        }

    }

    public static void updateAllMagneticMoments(singleSpin[] arr, double[][][] intTable, FieldTable momentTable, int maxIter, double extBx, double tol, double alpha, boolean suppressInternalTransFields){
        updateAllMagneticMoments(arr, null, intTable, momentTable, maxIter, extBx, tol, alpha, suppressInternalTransFields, -1);
    }

    @Deprecated
    public static double magneticMomentConvergence(singleSpin[] arr, FieldTable momentTable, double extBx, int flipSpin){
        double converged = 0;

        for (int i = 0; i < arr.length; i++) {
            if (i!=flipSpin) {
                double shouldBe = momentTable.getValue(arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz(), arr[i].getSpin());
                double is = arr[i].getSpinSize();

                converged += Math.abs(is - shouldBe);
            }
        }

        return Math.abs(converged/arr.length);
    }


    */

}
