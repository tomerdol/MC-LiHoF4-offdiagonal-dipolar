package simulation.montecarlo;

/**
 * Extracts measurement from a {@link Lattice} object
 */
public class ObservableExtractor {
    private double[][] k_cos_table, k_sin_table;

    /**
     * Constructs a new {@code ObservableExtractor} object
     * @param k_cos_table table with cos(k.R) for each site R in each direction k (k_x,k_y,k_z)
     * @param k_sin_table table with sin(k.R) for each site R in each direction k (k_x,k_y,k_z)
     */
    public ObservableExtractor(double[][] k_cos_table, double[][] k_sin_table){
        this.k_cos_table=k_cos_table;
        this.k_sin_table=k_sin_table;
    }

    /**
     * Calculates the mean (over the lattice) magnetic fields and their standard deviations
     * @param arr the lattice array
     * @return an array of means and standard deviations of Bx, By and Bz in the format: [<Bx>,std(Bx),<By>,std(By),<Bz>,std(Bz),maxB_transverse,maxB_long]
     */
    public static double[] meanField(singleSpin[] arr){

        double meanBx = 0, meanBy=0, meanBz=0, stdBx=0, stdBy=0, stdBz=0;
        double maxBtrans = 0, maxBlong = 0;
        for (int i=0;i<arr.length;i++){
            meanBx += arr[i].getLocalBx();
            stdBx += arr[i].getLocalBx()*arr[i].getLocalBx();
            meanBy += arr[i].getLocalBy();
            stdBy += arr[i].getLocalBy()*arr[i].getLocalBy();
            meanBz += arr[i].getLocalBz();
            stdBz += arr[i].getLocalBz()*arr[i].getLocalBz();

            if (Math.sqrt(arr[i].getLocalBx()*arr[i].getLocalBx() + arr[i].getLocalBy()*arr[i].getLocalBy()) > maxBtrans){
                maxBtrans = Math.sqrt(arr[i].getLocalBx()*arr[i].getLocalBx() + arr[i].getLocalBy()*arr[i].getLocalBy());
            }
            if (Math.abs(arr[i].getLocalBz()) > maxBlong){
                maxBlong = Math.abs(arr[i].getLocalBz());
            }
        }

        meanBx = meanBx / arr.length;
        meanBy = meanBy / arr.length;
        meanBz = meanBz / arr.length;
        if ((stdBx/arr.length) - meanBx*meanBx > 0) stdBx = Math.sqrt((stdBx/arr.length) - meanBx*meanBx);
        else stdBx=0;
        if ((stdBy/arr.length) - meanBy*meanBy > 0) stdBy = Math.sqrt((stdBy/arr.length) - meanBy*meanBy);
        else stdBy=0;
        if ((stdBz/arr.length) - meanBz*meanBz > 0) stdBz = Math.sqrt((stdBz/arr.length) - meanBz*meanBz);
        else stdBz=0;
        return new double[]{meanBx, stdBx, meanBy, stdBy, meanBz, stdBz, maxBtrans, maxBlong};	// format [<Bx>,std(Bx),<By>,std(By),<Bz>,std(Bz),max(B_transverse),max(Bz)]
    }

    /**
     * Calculates the k-space magnetization squared for k_min used for the finite-size correlation length
     * See, e.g. Andresen et al. (2013) PRL: https://link.aps.org/doi/10.1103/PhysRevLett.111.177202
     * @param arr array of spins
     * @return length 3 array with |m(k_min)|^2 for k_min in each of the spatial directions x,y,z
     */
    public double[] calc_mk2(singleSpin[] arr){
        double[] mk2 = new double[3];
        for (int dim=0;dim<3;dim++) {
            double mkCos=0, mkSin=0;
            for (int i = 0; i < arr.length; i++) {
                mkCos += arr[i].getSpinSize() * k_cos_table[dim][i];
                mkSin += arr[i].getSpinSize() * k_sin_table[dim][i];
            }
            mk2[dim] = (mkCos*mkCos + mkSin*mkSin)/(arr.length*arr.length);
        }
        return mk2;
    }

    /**
     * Calculates the total magnetization of the system
     * @param arr the spin system
     * @return the (normalized) magnetization
     */
    public static double calcMagnetization(singleSpin[] arr){
        double magnetization=0;
        int total=arr.length;
        for (int i=0; i<arr.length; i++){
            magnetization += arr[i].getSpinSize();
        }
        return magnetization/total;
    }


    /**
     * Calculates the total energy of the system
     * @param lattice the lattice whose energy should be calculated
     */
    public static double calcEnergy(Lattice lattice){
        singleSpin[] arr = lattice.getArray();
        FieldTable energyTable = lattice.energyTable;

        // calculate the total energy
        double energy=0;
        for(int i=0;i<arr.length;i++){
            energy += energyTable.getValue(arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz(), arr[i].getSpin(),arr[i].getPrevBIndices());
        }

        return energy;
    }

    /**
     * Calculates only the energy due to exchange interactions.
     * Not used since exchange interactions are already embedded in the standard interactions table
     * @param lattice the lattice whose energy should be calculated
     * @return the energy of the system due only to exchange interactions
     */
    public static double calcExchangeEnergy(Lattice lattice) {
        singleSpin[] arr = lattice.getArray();
        double[][] exchangeIntTable = lattice.exchangeIntTable;

        double exchangeEnergy = 0;
        for (int i=0; i<arr.length;i++){
            for (int j=i+1;j<arr.length;j++){
                exchangeEnergy += arr[i].getSpinSize()*arr[j].getSpinSize()*exchangeIntTable[i][j];
            }
        }
        return exchangeEnergy;
    }

    /**
     * Calculate mean and standard deviation of the magnetic moment size over the lattice
     * @param arr the lattice array
     * @return array of length 2 which holds the spin size mean and the spin size std.
     */
    public static double[] calcSpinSizes(singleSpin[] arr){
        double mean=0,std=0;
        for (int i=0;i<arr.length;i++){
            mean += Math.abs(arr[i].getSpinSize());
            std += arr[i].getSpinSize()*arr[i].getSpinSize();
        }
        mean = mean / arr.length;
        std = std / arr.length;
        std = std - mean*mean;

        return new double[]{mean, Math.sqrt(std)};
    }

    /**
     * Calculates the fraction of spins whose 4 nearest neighbors are oriented in one of the configurations that maximize the local transverse field on those spins.
     * These configurations are ones where for each of the two pairs of spins that share a common plane (which also includes the z axis), both spins are oriented opposite to one another.
     * The identification of such configurations relies on the order of neighbors in lattice.nnArray
     * @param lattice the lattice object
     * @return the fraction of such spins out of the total number of spins.
     */
    public static double countTransverseFieldMaximizingNNConfigs(Lattice lattice){
        int count=0;
        singleSpin[] arr = lattice.getArray();

        for (int i=0;i<arr.length;i++){
            // identify transverse field maximizing configurations:
            // the pairs of spins that share a common plane with the z axis are neighbors 0&1 and 2&3
            if (lattice.nnArray[i][0] >=0 && lattice.nnArray[i][1] >=0 &&
                    arr[lattice.nnArray[i][0]].getSpin()*arr[lattice.nnArray[i][1]].getSpin()==-1 &&
                    lattice.nnArray[i][2] >=0 && lattice.nnArray[i][3] >=0 &&
                    arr[lattice.nnArray[i][2]].getSpin()*arr[lattice.nnArray[i][3]].getSpin()==-1){
                boolean allNeighborSpinSizesLarge=true;
                for (int j=0; j<lattice.nnArray[i].length && allNeighborSpinSizesLarge; j++) {
                    allNeighborSpinSizesLarge &= arr[lattice.nnArray[i][j]].getSpin()*arr[lattice.nnArray[i][j]].getSpinSize() > 1.0;   // should (almost) always be positive
                }
                if (allNeighborSpinSizesLarge) count++;
            }
        }
        return 1.0*count/arr.length;
    }



}
