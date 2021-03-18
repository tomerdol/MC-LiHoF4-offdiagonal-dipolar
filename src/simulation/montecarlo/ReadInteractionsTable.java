package simulation.montecarlo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public abstract class ReadInteractionsTable {
    public static void receiveIntTable(double[][][] intTable, int Lx, int Lz){
        final double c = -Constants.mu_0*Constants.mu_B*Constants.g_L*0.25/Math.PI;	// coefficient dipolar spin-spin interaction. The minus sign is because
        // we use Ewald to calc the effective field and not the interaction

        int fileLx=0, fileLz=0;
        final int N=Lx*Lx*Lz*Constants.num_in_cell;
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

            for (int i=0;i<N;i++){
                for (int j=i;j<N;j++){
                    if ((str = in.readLine()) != null){
                        String[] xyzInteractions = str.split(",");
                        for (int k=0;k<xyzInteractions.length;k++) {

                            intTable[k][i][j] = c * Double.parseDouble(xyzInteractions[k]);
                            intTable[k][j][i] = c * Double.parseDouble(xyzInteractions[k]);
                            if (k==2){  // for the zz interactions we multiply by half to avoid double counting.
                                intTable[k][i][j] *= 0.5;
                                intTable[k][j][i] *= 0.5;
                            }
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

    public abstract int[][] exchangeInt(double[][] intTable, int Lx, int Ly, int Lz, double J_ex);

    static void addExchangeToNeighbor(int focusSpin, int[] neighbors, int[][] nnArray, boolean[][] nnArray_test, double[][] intTable, double J_ex){
        for (int i=0;i<neighbors.length;i++){
            intTable[focusSpin][neighbors[i]]+=J_ex;
            intTable[neighbors[i]][focusSpin]+=J_ex;

            nnArray[focusSpin][i]=neighbors[i];
            nnArray[neighbors[i]][i]=focusSpin;

            nnArray_test[focusSpin][neighbors[i]]=true;
            nnArray_test[neighbors[i]][focusSpin]=true;
        }
    }
}
