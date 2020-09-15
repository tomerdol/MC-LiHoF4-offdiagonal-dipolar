package simulation.montecarlo;


import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;

import static org.junit.Assert.assertTrue;
@RunWith(Parameterized.class)
public class FieldTableTest {
    FieldTable momentTable;
    double bx, by, bz;
    boolean s;

    public FieldTableTest(double bx, double by, double bz, boolean s){
        this.bx=bx;
        this.by=by;
        this.bz=bz;
        this.s=s;
    }

    @Parameterized.Parameters(name = "{index}: ({0},{1},{2})")
    public static Collection<Object[]> input(){
//        double minB=-2.02857006;
//        double maxB=2.02857006;
        double minB=-1.2;
        double maxB=1.2;
        int numOfBs=11;
        Object[][] arr = new Object[numOfBs*numOfBs*numOfBs][];

        // methodological search
        /*
        for (int i=0;i<numOfBs;i++){
            for (int j=0;j<numOfBs;j++){
                for (int k=0;k<numOfBs;k++){
                    // for spin down change 2nd parameter to false
                    arr[i*numOfBs*numOfBs+j*numOfBs+k+0] = new Object[]{minB+i*(maxB-minB)/numOfBs,minB+j*(maxB-minB)/numOfBs,minB+k*(maxB-minB)/numOfBs,true};
                }
            }
        }

         */

        // random search
        Random rnd = new Random();
        for (int i=0;i<arr.length;i++){
            arr[i] = new Object[]{minB+(maxB-minB)*rnd.nextDouble(),minB+(maxB-minB)*rnd.nextDouble(),minB+(maxB-minB)*rnd.nextDouble(),true};
        }

        return Arrays.asList(arr);
    }

    @Before
    public void setUp() throws Exception {
        double extBx = 0.0;
        momentTable = FieldTable.of(String.format("magnetic_moment_up_arr_%1.2f",extBx), true);
    }

    @Test
    public void getValue() {
        double spin = s ? 1 : -1;
        double tableValue = momentTable.getValue(bx, by, bz, spin);
        double exactValue = momentTable.manualCalcValue(bx, by, bz*(spin), true)*Math.signum(spin);
        // test whether relative error is <5%
        System.out.println(Math.abs(tableValue-exactValue)/Math.abs(exactValue));
        assertTrue(Math.abs(tableValue-exactValue)/Math.abs(exactValue)<.05);
    }
}