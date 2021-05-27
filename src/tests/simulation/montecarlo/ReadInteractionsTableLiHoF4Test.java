package simulation.montecarlo;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

@RunWith(Parameterized.class)
public class ReadInteractionsTableLiHoF4Test{
    double[][][] intTable;
    int Lx;
    int Lz;
    double J_ex;
    boolean[] dilution;

    public ReadInteractionsTableLiHoF4Test(double[][][] intTable, int Lx, int Lz, double J_ex, boolean[] dilution) {
        this.intTable = intTable;
        this.Lx = Lx;
        this.Lz = Lz;
        this.J_ex = J_ex;
        this.dilution = dilution;
    }

    @Parameterized.Parameters(name = "{index}: ({0},{1},{2},{3})")
    public static Collection<Object[]> input(){
        int Lx=7, Lz=7;
        int tests=10;
        int n=4*Lx*Lx*Lz;
        double x=0.3;
        Object[][] arr = new Object[tests][];
        // random search
        Random rnd = new Random();
        for (int i=0;i<arr.length;i++){
            boolean[] inputArray = new boolean[n];
            int N=0;
            for (int j=0;j<inputArray.length;j++){
                if (rnd.nextDouble() < x) {
                    inputArray[j] = true;
                    N++;
                }else{
                    inputArray[j] = false;
                }
            }

            final double[][][] intTable = new double[3][N][N];
            arr[i] = new Object[]{intTable, Lx, Lz, 1.16e-3, inputArray};
        }

        return Arrays.asList(arr);
    }



    @Test
    public void exchangeInt() throws Exception {
        ReadInteractionsTable.receiveIntTable(intTable, Lx, Lz, dilution);	// get interaction table from file
        boolean[] noDilution = new boolean[4*Lx*Lx*Lz];
        for (int i=0;i<noDilution.length;i++) noDilution[i] = true;
        double[][][] intTable2 = new double[3][4*Lx*Lx*Lz][4*Lx*Lx*Lz];
        ReadInteractionsTable.receiveIntTable(intTable2, Lx, Lz, noDilution);	// get interaction table from file
        int n = 0;
        for (int i=0;i<dilution.length;i++){
            if (dilution[i]) n++;
//            System.out.print(dilution[i] + " ");
        }
//        System.out.println();

        final double[][] exchangeIntTable = new double[n][n];

        ReadInteractionsTable interactionsTableReceiver;
        if (System.getProperty("system").equals("LiHoF4")){
            interactionsTableReceiver = new ReadInteractionsTableLiHoF4(dilution);
        } else if (System.getProperty("system").equals("Fe8")){
            interactionsTableReceiver = new ReadInteractionsTableFe8();
        } else {
            throw new RuntimeException("Could not read interactions table. Illegal system name given.");
        }

        int[][] nnArray = interactionsTableReceiver.exchangeInt(exchangeIntTable, Lx, Lx, Lz, J_ex);	// receive the nearest neighbor array and fill exchangeIntTable with the exchange interaction values

        final double[][] exchangeIntTable2 = new double[4*Lx*Lx*Lz][4*Lx*Lx*Lz];
        ReadInteractionsTable interactionsTableReceiver2;
        if (System.getProperty("system").equals("LiHoF4")){
            interactionsTableReceiver2 = new ReadInteractionsTableLiHoF4(noDilution);
        } else if (System.getProperty("system").equals("Fe8")){
            interactionsTableReceiver2 = new ReadInteractionsTableFe8();
        } else {
            throw new RuntimeException("Could not read interactions table. Illegal system name given.");
        }
        int[][] nnArray2 = interactionsTableReceiver2.exchangeInt(exchangeIntTable2, Lx, Lx, Lz, J_ex);	// receive the nearest neighbor array and fill exchangeIntTable with the exchange interaction values

        // add exchange to intTable
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                intTable[2][i][j] += -0.5*exchangeIntTable[i][j]*Constants.k_B/(Constants.mu_B*Constants.g_L);
            }
        }

        // add exchange to intTable2
        for (int i=0;i<4*Lx*Lx*Lz;i++){
            for (int j=0;j<4*Lx*Lx*Lz;j++){
                intTable2[2][i][j] += -0.5*exchangeIntTable2[i][j]*Constants.k_B/(Constants.mu_B*Constants.g_L);
            }
        }

        int[] correspondenceArray = new int[n];
        int count=0;
        for (int i=0;i<dilution.length;i++){
            if (dilution[i]){
                correspondenceArray[count++] = i;
            }
        }

        // testing intTable
        for (int i = 0; i < intTable[0].length; i++) {
            for (int j = 0; j < intTable[0][i].length; j++) {
                for (int k=0; k<intTable.length; k++) {
//                    System.out.print(intTable[k][i][j] + ",");
                    assertEquals(intTable[k][i][j], intTable2[k][correspondenceArray[i]][correspondenceArray[j]], 0.0000001);
                }
//                System.out.println();
            }
//            System.out.println();
        }

        // testing nnArray
        for (int i = 0; i < nnArray.length; i++) {
            System.out.print(i + "-" + correspondenceArray[i] + " : ");
            for (int j = 0; j < nnArray[i].length; j++) {
                if (nnArray[i][j]>=0){
                    System.out.print(correspondenceArray[nnArray[i][j]] + ",");
                    System.out.print(nnArray2[correspondenceArray[i]][j] + ",");
                    assertEquals(correspondenceArray[nnArray[i][j]], nnArray2[correspondenceArray[i]][j]);
                }
            }
            System.out.println();
        }

        assertTrue(nnArray.length == n);

    }
}