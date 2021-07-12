package simulation.montecarlo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Reads the interaction between pairs of spins from a file
 */
public abstract class ReadInteractionsTable {
    /** Array which for each spin in an undiluted system, hold its index in the compact array
     * that represents the diluted system */
    int[] correspondenceArray;
    /** total number of spins in the (possibly diluted) system */
    int N;

    /**
     * Fills the {@link #correspondenceArray correspondence array}
     * @param dilution {@code boolean} array which indicates which spins exist and which do not
     */
    protected void setCorrespondenceArray(boolean[] dilution){
        correspondenceArray = new int[dilution.length];
        int count=0;
        for (int i=0;i<dilution.length;i++){
            if (dilution[i]){
                correspondenceArray[i] = count++;
            } else {
                correspondenceArray[i] = -1;
            }
        }
        N=count;
    }

    /**
     * Receive the interactions table for an undiluted system
     * @param intTable interaction table array to fill
     * @param Lx number of unit cells in the x and y directions
     * @param Lz number of unit cells in the z direction
     */
    public static void receiveIntTable(double[][][] intTable, int Lx, int Lz){
        boolean[] noDilution = new boolean[Constants.num_in_cell*Lx*Lx*Lz];
        for (int i=0;i<noDilution.length;i++) noDilution[i] = true;
        receiveIntTable(intTable, Lx, Lz, noDilution);
    }

    /**
     * Receive the interactions table for a diluted system.
     * The interactions table is a 2D symmetric array which holds
     * in cell [i,j] the interaction between spins i and j.
     * There are 3 such arrays (first dimension of intTable) which hold the
     * xz, yz and zz interactions in that order.
     * The file only holds half of the symmetric array (incl. the diagonal) but the
     * values read fill the full symmetric array.
     * @param intTable interaction table array to fill
     * @param Lx number of unit cells in the x and y directions
     * @param Lz number of unit cells in the z direction
     * @param dilution {@code boolean} array which indicates which spins exist and which do not
     */
    public static void receiveIntTable(double[][][] intTable, int Lx, int Lz, boolean[] dilution){
        final double c = -Constants.mu_0*Constants.mu_B*Constants.g_L*0.25/Math.PI;	// coefficient dipolar spin-spin interaction. The minus sign is because
                                                                                    // we use Ewald to calc the effective field and not the interaction

        int fileLx=0, fileLz=0;
        final int N=Lx*Lx*Lz*Constants.num_in_cell;
        // the file starts with the lines:
        // Lx=?
        // Lz=?
        // followed by three lines or parameters used in its creation (not needed here)
        try (BufferedReader in = new BufferedReader(new FileReader(System.getProperty("system") + File.separator + "data" + File.separator + "interactions" + File.separator + "intTable_"+Lx+"_"+Lz+".txt"))){
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

            // indices i & j go over the full interaction table (without regard for dilution)
            // indices m & n go over the compact interaction table which only holds interactions
            // between existing spins
            for (int i=0, m=0;i<N;i++){
                for (int j=i, n=m;j<N;j++){
                    if ((str = in.readLine()) != null){
                        if (dilution[i] && dilution[j]){
                            // each line hold 3 interactions: xz, yz, zz
                            // separated by commas
                            String[] xyzInteractions = str.split(",");
                            for (int k = 0; k < xyzInteractions.length; k++) {

                                intTable[k][m][n] = c * Double.parseDouble(xyzInteractions[k]);
                                intTable[k][n][m] = c * Double.parseDouble(xyzInteractions[k]);

                                if (k == 2) {  // for the zz interactions we multiply by half to avoid double counting.
                                    intTable[k][m][n] *= 0.5;
                                    if (m!=n) intTable[k][n][m] *= 0.5;
                                }
                            }
                            n++;
                        }
                    } else {
                        System.err.println("reached end of interactions file before int table is full!");
                        System.exit(1);
                    }
                }
                if (dilution[i]){
                    m++;
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

    /**
     * Calculates the zz exchange interaction with nearest neighbors
     * and also returns an array of nearest neighbors
     * @param intTable exchange interaction array to fill
     * @param Lx number of unit cells in the x directions
     * @param Ly number of unit cells in the y directions
     * @param Lz number of unit cells in the z directions
     * @param J_ex strength of the exchange interaction
     * @return array of nearest neighbors: cell [i,j] hold the index (in the compact array) of the j-th neighbor of spin i
     */
    public abstract int[][] exchangeInt(double[][] intTable, int Lx, int Ly, int Lz, double J_ex);

    /**
     * Adds the exchange interaction to the neighbors of the given spin and adds them as its neighbors in nnArray
     * @param focusSpin spin in focus
     * @param neighbors array of the neighbors of {@code focusSpin}, by their compact indices
     * @param nnArray array of nearest neighbors to fill
     * @param nnArray_test full 2D array used to indicate the spins i and j (compact indices) are neighbors
     *                     by having the cells [i,j] and [j,i] set to {@code true}.
     * @param intTable array of exchange interaction to fill
     * @param J_ex strength of the exchange interaction
     */
    static void addExchangeToNeighbor(int focusSpin, int[] neighbors, int[][] nnArray, boolean[][] nnArray_test, double[][] intTable, double J_ex){
        // loop through the given neighbors of focusSpin and add the exchange interaction between them and focusSpin to the intTable
        for (int i=0;i<neighbors.length;i++){
            if (focusSpin >=0 && neighbors[i] >=0) {
                intTable[focusSpin][neighbors[i]] += J_ex;
                intTable[neighbors[i]][focusSpin] += J_ex;

                // also set them as neighbors
                nnArray_test[focusSpin][neighbors[i]] = true;
                nnArray_test[neighbors[i]][focusSpin] = true;
            }

            // only fill nnArray rows that correspond to existing spins,
            // but do allow -1 values for non-existing neighbors
            // (this should be checked when using nnArray)
            if (focusSpin >=0) nnArray[focusSpin][i] = neighbors[i];
            if (neighbors[i] >=0) nnArray[neighbors[i]][i] = focusSpin;

        }
    }
}
