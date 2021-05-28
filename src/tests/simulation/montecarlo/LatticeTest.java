package simulation.montecarlo;

import org.apache.commons.lang3.ArrayUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

@RunWith(Parameterized.class)
public class LatticeTest {
    double[][][] intTable;
    int Lx;
    int Lz;
    double J_ex;
    boolean[] dilution;

    public LatticeTest(double[][][] intTable, int Lx, int Lz, double J_ex, boolean[] dilution) {
        this.intTable = intTable;
        this.Lx = Lx;
        this.Lz = Lz;
        this.J_ex = J_ex;
        this.dilution = dilution;
    }

    @Parameterized.Parameters(name = "{index}: ({0},{1},{2},{3})")
    public static Collection<Object[]> input(){
        int Lx=7, Lz=7;
        int tests=100;
        int n=4*Lx*Lx*Lz;
        double x=0.6;
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
    public void orderBFS() {
        double extBx=0.0;   // external Bx
        double extBy=0;   // external By
        boolean suppressInternalTransFields=false;
        double spinSize=5.5;

        ReadInteractionsTable.receiveIntTable(intTable, Lx, Lz, dilution);	// get interaction table from file
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

        // add exchange to intTable
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                intTable[2][i][j] += -0.5*exchangeIntTable[i][j]*Constants.k_B/(Constants.mu_B*Constants.g_L);
            }
        }

        Lattice lattice = new Lattice(Lx, Lz, extBx, extBy, suppressInternalTransFields, spinSize,dilution, intTable, exchangeIntTable, nnArray, null, null, null);

        for (int test=0; test<n; test++){
            System.out.println("*** testing spins " + test + " ***");
            int[] orderBFS = lattice.orderBFS(n, test);

            boolean verify = orderBFS.length == n; // if orderBFS is of different length than n there is
                                                            // obviously a problem
            for (int i=0;i<n && verify;i++){
                // check that each number between 0 and the number of spins appears in orderBFS
                verify = verify && ArrayUtils.contains(orderBFS, i);
            }
            assertTrue(verify);
        }

    }

    @Test
    public void getCompactArrayIndex() {
        int n = 0;
        for (int i=0;i<dilution.length;i++){
            if (dilution[i]) n++;
//            System.out.print(dilution[i] + " ");
        }

        ReadInteractionsTable interactionsTableReceiver;
        if (System.getProperty("system").equals("LiHoF4")){
            interactionsTableReceiver = new ReadInteractionsTableLiHoF4(dilution);
        } else if (System.getProperty("system").equals("Fe8")){
            interactionsTableReceiver = new ReadInteractionsTableFe8();
        } else {
            throw new RuntimeException("Could not read interactions table. Illegal system name given.");
        }

        for (int i=0; i<4*Lx*Lx*Lz; i++){
            assertEquals(interactionsTableReceiver.correspondenceArray[i], Lattice.getCompactArrayIndex(dilution, i));
        }
    }
}