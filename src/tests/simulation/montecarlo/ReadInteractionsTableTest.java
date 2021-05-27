package simulation.montecarlo;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import static org.junit.Assert.*;
@RunWith(Parameterized.class)
public class ReadInteractionsTableTest {
    double[][][] intTable;
    int Lx;
    int Lz;
    boolean[] dilution;

    public ReadInteractionsTableTest(double[][][] intTable, int lx, int lz, boolean[] dilution) {
        this.intTable = intTable;
        Lx = lx;
        Lz = lz;
        this.dilution = dilution;
    }

    @Parameterized.Parameters(name = "{index}: ({0},{1},{2},{3})")
    public static Collection<Object[]> input(){
        int Lx=7, Lz=7;
        int tests=10;
        int n=4*Lx*Lx*Lz;
        double x=0.5;
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
            arr[i] = new Object[]{intTable, Lx, Lz, inputArray};
        }

        return Arrays.asList(arr);
    }


    @Test
    public void receiveIntTable() {
        ReadInteractionsTable.receiveIntTable(intTable, Lx, Lz, dilution);	// get interaction table from file
        boolean[] noDilution = new boolean[4*Lx*Lx*Lz];
        for (int i=0;i<noDilution.length;i++) noDilution[i] = true;
        double[][][] intTable2 = new double[3][4*Lx*Lx*Lz][4*Lx*Lx*Lz];
        ReadInteractionsTable.receiveIntTable(intTable2, Lx, Lz, noDilution);	// get interaction table from file
        int n = 0;
        for (int i=0;i<dilution.length;i++){
            if (dilution[i]) n++;
            System.out.print(dilution[i] + " ");
        }
        System.out.println();

        int[] correspondenceArray = new int[n];
        int count=0;
        for (int i=0;i<dilution.length;i++){
            if (dilution[i]){
                correspondenceArray[count++] = i;
            }
        }

        // testing
        for (int i = 0; i < intTable[0].length; i++) {
            for (int j = 0; j < intTable[0][i].length; j++) {
                for (int k=0; k<intTable.length; k++) {
                    System.out.print(intTable[k][i][j] + ",");
                    assertEquals(intTable[k][i][j], intTable2[k][correspondenceArray[i]][correspondenceArray[j]], 0.0000001);

                }
                System.out.println();
            }
            System.out.println();
        }
    }


}