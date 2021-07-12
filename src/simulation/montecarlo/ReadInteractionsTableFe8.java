package simulation.montecarlo;

/**
 * Reads the interaction between pairs of spins from a file for the Fe8 system.
 * Does not support dilution at this stage.
 */
public class ReadInteractionsTableFe8 extends ReadInteractionsTable{

    /**
     * Calculates zz exchange interaction with nearest neighbors for Fe8
     * and also returns an array of nearest neighbors.
     * Since the exchange interactions in the Fe8 system are negligible, {@code J_ex} should be
     * zero, and the main use of this is to get the returned nearest-neighbor array
     * @param intTable exchange interaction array to fill
     * @param Lx number of unit cells in the x directions
     * @param Ly number of unit cells in the y directions
     * @param Lz number of unit cells in the z directions
     * @param J_ex strength of the exchange interaction
     * @return array of nearest neighbors: cell [i,j] hold the index (in the compact array) of the j-th neighbor of spin i
     */
    public int[][] exchangeInt(double[][] intTable, int Lx, int Ly, int Lz, double J_ex){
        final int N = Lx*Ly*Lz*Constants.num_in_cell;
        final int numOfNeighbors;
        // for a single unit cell system, there are no neighbors (a spin cannot be its own neighbor)
        if (Lx==1 && Ly==1 && Lz==1) numOfNeighbors = 0;
        // similarly, if the system has only one unit cell in any direction, there are no neighbors in that direction
        else if ((Lx==1 && Ly==1 && Lz==2) || (Lx==1 && Ly==2 && Lz==1) || (Lx==2 && Ly==1 && Lz==1)) numOfNeighbors = 1;
        else if ((Lx==1 && Ly==2 && Lz==2) || (Lx==2 && Ly==2 && Lz==1) || (Lx==2 && Ly==1 && Lz==2)) numOfNeighbors = 2;
        else numOfNeighbors = 6;    // number of nearest neighbors for each spin.
        int[][] nnArray = new int[N][numOfNeighbors];;	// nearest neighbor array

        boolean[][] nnArray_test = new boolean[N][N];	// for testing the validity of the returned nnArray

        // neighbor numbers are as follows
        // neighbor1: along +x
        // neighbor2: along -x
        // neighbor3: along +y
        // neighbor4: along -y
        // neighbor5: along +z
        // neighbor6: along -z

        for (int i=0;i<Lx;i++){
            for (int j=0;j<Ly;j++){
                for (int k=0;k<Lz;k++){
                    // nearest neighbors, including periodic boundary conditions
                    int focusSpin = i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell;

                    int neighbor1=((i+1)%Lx)*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell;
                    int neighbor2=((Lx+i-1)%Lx)*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell;
                    int neighbor3=i*Ly*Lz*Constants.num_in_cell+((j+1)%Ly)*Lz*Constants.num_in_cell+k*Constants.num_in_cell;
                    int neighbor4=i*Ly*Lz*Constants.num_in_cell+((Ly+j-1)%Ly)*Lz*Constants.num_in_cell+k*Constants.num_in_cell;
                    int neighbor5=i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+((k+1)%Lz)*Constants.num_in_cell;
                    int neighbor6=i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+((Lz+k-1)%Lz)*Constants.num_in_cell;

                    // put interactions in intTable and nearest neighbor indices in nnArray
                    addExchangeToNeighbor(focusSpin, new int[]{neighbor1, neighbor2, neighbor3, neighbor4, neighbor5, neighbor6}, nnArray, nnArray_test, intTable, J_ex);
                }
            }
        }

        // verify that all nearest neighbor pairs are accounted for in nnArray
        boolean validNNArray=true;
        for (int i=0;i<nnArray_test.length && validNNArray;i++){
            int countNearestNeighbors=0;
            for (int j=0;j<nnArray_test[i].length && validNNArray;j++){
                if (nnArray_test[i][j]) countNearestNeighbors++;
            }
            if (countNearestNeighbors!=numOfNeighbors) validNNArray=false;
        }
        if (validNNArray){ return nnArray; }
        else {
            System.err.println("There was an error creating the nearest neighbor array. at least one of the spins has more or less than 4 neighbors.");
            System.exit(1);
            return null;
        }
    }

}
