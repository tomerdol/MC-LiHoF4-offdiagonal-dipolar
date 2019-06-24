package simulation.montecarlo;

public class ObservableExtractor {
    private double[] k_cos_table, k_sin_table;

    public ObservableExtractor(double[] k_cos_table, double[] k_sin_table){
        this.k_cos_table=k_cos_table;
        this.k_sin_table=k_sin_table;
    }

    /**
     * calculates the mean (over the lattice) magnetic fields and their standard deviations
     * @param arr - the lattice array
     * @return an array of means and standard deviations of Bx, By and Bz in the format: [<Bx>,std(Bx),<By>,std(By),<Bz>,std(Bz),maxB_transverse]
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
        return new double[]{meanBx, stdBx, meanBy, stdBy, meanBz, stdBz, maxBtrans, maxBlong};	// format [<Bx>,std(Bx),<By>,std(By)]
    }


    public double calc_mk2(singleSpin[] arr){
        double mkx=0, mky=0;
        for (int i=0;i<arr.length;i++){
            //mkx += arr[i].getSpin()*k_cos_table[i];
            //mky += arr[i].getSpin()*k_sin_table[i];
            //mkx += arr[i].getSpin()*arr[i].getSpinSize()*k_cos_table[i];
            //mky += arr[i].getSpin()*arr[i].getSpinSize()*k_sin_table[i];
            mkx += arr[i].getSpinSize()*k_cos_table[i];
            mky += arr[i].getSpinSize()*k_sin_table[i];
        }
        return (mkx*mkx + mky*mky)/(arr.length*arr.length);
    }


    public static double calcMagnetization(singleSpin[] arr){
        double magnetization=0;
        int total=0;
        for (int i=0; i<arr.length; i++){
            //magnetization += arr[i].getSpin();
            magnetization += arr[i].getSpinSize();
            //magnetization += arr[i].getSpin()*arr[i].getSpinSize();
            total+=Math.abs(arr[i].getSpin());
        }
        return magnetization/total;
    }


    // calculates the total energy of the system
    public static double calcEnergy(Lattice lattice){
        singleSpin[] arr = lattice.getArray();
        FieldTable energyTable = lattice.energyTable;

        //double c = Constants.g_L*Constants.mu_B/Constants.k_B;
        // calculate the total energy
        double energy=0;
        for(int i=0;i<arr.length;i++){
            energy += energyTable.getValue(arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz(), arr[i].getSpin(),arr[i].getPrevBIndices());
            //energy += arr[i].getSpin()*arr[i].getSpinSize()*c*arr[i].getLocalBz();
        }

        return 0.5*energy + calcExchangeEnergy(lattice);
    }

    public static double calcExchangeEnergy(Lattice lattice) {
        singleSpin[] arr = lattice.getArray();
        double[][] exchangeIntTable = lattice.exchangeIntTable;

        double exchangeEnergy = 0;
        for (int i=0; i<arr.length;i++){
            for (int j=i+1;j<arr.length;j++){
                //exchangeEnergy += arr[i].getSpin()*arr[i].getSpinSize()*arr[j].getSpin()*arr[j].getSpinSize()*exchangeIntTable[arr[i].getN()][arr[j].getN()];
                exchangeEnergy += arr[i].getSpinSize()*arr[j].getSpinSize()*exchangeIntTable[arr[i].getN()][arr[j].getN()];
            }
        }
        return exchangeEnergy;
    }


    /**
     * Calculate mean and std of the magnetic moment size over the lattice
     * @param arr - the lattice array
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




}
