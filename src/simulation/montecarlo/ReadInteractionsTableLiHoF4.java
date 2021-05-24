package simulation.montecarlo;

public class ReadInteractionsTableLiHoF4 extends ReadInteractionsTable{

    //calculates exchange interaction with nearest neighbors
    // also returns an array of nearest neighbors
    public int[][] exchangeInt(double[][] intTable, int Lx, int Ly, int Lz, double J_ex){
        final int N = Lx*Ly*Lz*Constants.num_in_cell;
        final int numOfNeighbors = (Lx==1 && Ly==1 && Lz==1) ? 2 : 4;    // number of nearest neighbors for each spin.
        // Lx=Ly=Lz=1 is a unique case where each spin has only 2
        // distinct neighbors
        int[][] nnArray = new int[N][4];;	// nearest neighbor array

        boolean[][] nnArray_test = new boolean[N][N];	// for testing

        // neighbor numbers are as follows, for the 0th and 2nd atoms in the base (with respect to ion positions_1.pdf):
        // neighbor1: up-right
        // neighbor2: up-left
        // neighbor3: down-outward
        // neighbor4: down-inward

        // for the 1st and 3rd atoms, it is the same, up to parity inversion, i.e.,
        // neighbor1: down-left
        // neighbor2: down-right
        // neighbor3: up-inward
        // neighbor4: up-outward
        // notice this is relevant mostly for how the nearest-neighbors are ordered in nnArray, which is filled in addExchangeToNeighbor()

        for (int i=0;i<Lx;i++){
            for (int j=0;j<Ly;j++){
                for (int k=0;k<Lz;k++){
                    // notice we are only going through the 0th and 2nd atoms in the base since they participate in all exchange interactions

                    // nearest neighbors to 0th base atom, including periodic boundary conditions
                    int focusSpin = i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+0;

                    int neighbor1=i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+1;
                    int neighbor2=((Lx+i-1)%Lx)*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+1;
                    int neighbor3=i*Ly*Lz*Constants.num_in_cell+((Ly+j-1)%Ly)*Lz*Constants.num_in_cell+((Lz+k-1)%Lz)*Constants.num_in_cell+3;
                    int neighbor4=i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+((Lz+k-1)%Lz)*Constants.num_in_cell+3;

                    // put interactions in intTable and nearest neighbor indices in nnArray
                    addExchangeToNeighbor(focusSpin, new int[] {neighbor1, neighbor2, neighbor3, neighbor4}, nnArray, nnArray_test, intTable, J_ex);



                    // now nearest neighbors to 2nd base atom, including periodic boundary conditions
                    focusSpin = i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+2;

                    neighbor3=i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+1;
                    neighbor2=i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+3;
                    neighbor1=((i+1)%Lx)*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+3;
                    neighbor4=i*Ly*Lz*Constants.num_in_cell+((j+1)%Ly)*Lz*Constants.num_in_cell+k*Constants.num_in_cell+1;
                    // put interactions in intTable and nearest neighbor indices in nnArray
                    addExchangeToNeighbor(focusSpin, new int[]{neighbor1, neighbor2, neighbor3, neighbor4}, nnArray, nnArray_test, intTable, J_ex);

                }
            }
        }

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