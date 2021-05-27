package simulation.montecarlo;

import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

public class ObservableExtractorTest {
    Lattice lattice;
    ObservableExtractor measure = new ObservableExtractor(null, null);
    int Lx=10;	// lattice x-y size
    int Lz=10;	// lattice z size
    @Before
    public void setUp() throws Exception {
        double J_ex = 0.0;
        String interpolationTableFileNameExtension = "_0.014";
        double extBx = 0.0;
        double extBy = 0.0;
        boolean suppressInternalTransFields = true;
        double spinSize = 10;
        final double[][] exchangeIntTable = new double[Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz];
        final double[][][] intTable = new double[3][Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz]; // create interaction table that holds all the dipolar interactions. will be full even though it's symmetric. 1st array is x,y,z term
        boolean[] dilution = new boolean[Constants.num_in_cell*Lx*Lx*Lz];
        for (int i=0; i< dilution.length; i++) dilution[i]=true;

        ReadInteractionsTable interactionsTableReceiver;
        if (System.getProperty("system").equals("LiHoF4")){
            interactionsTableReceiver = new ReadInteractionsTableLiHoF4(dilution);
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

        final int N=Lx*Lx*Lz*Constants.num_in_cell;
        for (int i=0;i<N;i++){
            for (int j=0;j<N;j++){
                intTable[2][i][j] += -0.5*exchangeIntTable[i][j]*Constants.k_B/(Constants.mu_B*Constants.g_L);
            }
        }
        FieldTable energyTable = FieldTable.of(String.format("energy_up_arr_%1.2f"+interpolationTableFileNameExtension,extBx), false);
        FieldTable momentTable = FieldTable.of(String.format("magnetic_moment_up_arr_%1.2f"+interpolationTableFileNameExtension,extBx), true);

        lattice = new Lattice(Lx, Lz, extBx, extBy, suppressInternalTransFields, spinSize, intTable, exchangeIntTable, nnArray, energyTable, momentTable, measure);
    }

    @Test
    public void singleEnergy() {
        boolean convergedConfig=false;
        int index=0;
        while (!convergedConfig) {
            lattice.updateAllLocalFields();
            lattice.updateAllMagneticMoments(160, 1.0e-4, 0.95);

            if (lattice.magneticMomentConvergence() > 1.0e-4){
                convergedConfig=false;
            }else{
                convergedConfig=true;
            }

            if (!convergedConfig) System.out.println("setup convergence index: " + ++index);
        }

        singleSpin[] arr = lattice.getArray();
        int i=0;
        System.out.println("Single molecule energy: " + lattice.energyTable.getValue(arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz(), arr[i].getSpin(),arr[i].getPrevBIndices()));
    }

    @Test
    public void calcEnergy() {
        System.out.println(ObservableExtractor.calcEnergy(lattice) / (10*10*10));
//        System.out.println(lattice.getArray()[0].);
    }

    @Test
    public void calcSpinSizes() {
        System.out.println(ObservableExtractor.calcSpinSizes(lattice.getArray())[0]);
    }

    @Test
    public void meanField() {
        System.out.println(Arrays.toString(ObservableExtractor.meanField(lattice.getArray())));
    }
}