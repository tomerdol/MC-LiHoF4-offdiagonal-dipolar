package simulation.montecarlo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Random;

/**
 * Holds a table of precalculated values and performs interpolation to get requested values.
 * Used to hold energy & magnetic moment of a Ho ion under different (Bx,By,Bz)
 */
final public class FieldTable {
    private final double[][][] table;
    private final double[][] values;
    private final boolean transposedRequiresMinus;      // when calculating the magnetic moment of a down-spin, one can calculate the magnetic moment of an up-spin in with Bz->-Bz, but the
                                                        // result requires a minus sign. hence for FieldTable that contains magnetic moments, transposedRequiresMinus should be true.
    private final String name;
    private final int[] dj;   // min distance to determine whether consecutive calls are correlated
    //private final boolean[] cor;    // tracks whether consecutive calls are correlated (separately for x,y,z)

    /**
     * Reads the energy_arr/magnetic_moment_arr file into a table of double and return a FieldTable object
     * @param fileName - name of table file to be read (without .txt extension)
     * @return FieldTable object that has a table of average energy drop or magnetic moment due to different combinations of local fields (Bx,By,Bz) and also the Bx,By,Bz values themselves
     */
    public static FieldTable of(String fileName, boolean transposedRequiresMinus) {
        double[][][] table = null;
        double[][] Bvalues = new double[3][];
        try{
            String line=null;
            String[] values=null;
            @SuppressWarnings("resource")
            BufferedReader br = new BufferedReader(new FileReader("interactions" + File.separator + fileName + ".txt"));

            // first line is a header
            line = br.readLine();
            if (line!=null){
                // parse array shape line
                String[] shape = line.replaceAll("\\(|\\)|\\#|\\s|\\t","").split(",");
                for (int i=0;i<Bvalues.length;i++){
                    Bvalues[i]=new double[Integer.parseInt(shape[i])];
                }

                // parse B values (Bz,By,Bx)
                for (int i=0;i<Bvalues.length;i++){
                    try{
                        line = br.readLine();
                        String[] singleBvaleus = line.split("\\[|\\]")[1].replaceAll("\\(|\\)|\\#|\\s|\\t","").split(",");
                        Bvalues[i]=toDouble(singleBvaleus,0,singleBvaleus.length);
                    } catch (ArrayIndexOutOfBoundsException e){
                        throw new RuntimeException("Mismatch between shape at head of table file and values at next 3 lines");
                    }
                }

                // start reading table
                table = new double[Bvalues[0].length][Bvalues[1].length][Bvalues[2].length];
                int Yindex=0;
                int Zindex=0;
                while ((line = br.readLine()) != null) {
                    if (!line.equals("# New Bz slice")) {
                        values = line.split("\\t|\\s");
                        table[Zindex][Yindex++] = toDouble(values, 0, values.length);
                    }else{
                        Zindex++;
                        Yindex=0;
                    }
                }
            }else{
                System.err.println("transverse field table file is empty");
            }
        }catch(Exception e){
            e.printStackTrace();
        }

        if (table!=null) return new FieldTable(table,Bvalues,transposedRequiresMinus, fileName);
        else throw new RuntimeException("trans_field table not read");
    }


    /**
     * Turns an array of strings into an array of doubles
     * @param inArray - input array (strings)
     * @return an array of doubles represented by the strings
     */
    public static double[] toDouble(String[] inArray, int start, int end){
        double[] retArray = new double[end-start];
        int index=0;
        for (int i=start;i<end;i++){
            retArray[index++] = Double.parseDouble(inArray[i]);
        }

        return retArray;
    }

    /**
     * FieldTable constructor
     * @param table 3D double array of the values (usually read from a file by {@link FieldTable#of(String, boolean)}
     * @param values 2D double array (must be 3 by N) representing the B values (Bx,By,Bz)
     * @param transposedRequiresMinus0 when calculating the magnetic moment of a down-spin, one can calculate the magnetic moment of an up-spin in with Bz->-Bz, but the
     *                                 result requires a minus sign. Hence for a FieldTable that contains magnetic moments, transposedRequiresMinus should be true, and
 *                                     for a FieldTable that contains energies it should be false.
     * @param name Name of this FieldTable. Options are "magnetic_moment_up_broyden" and "energy_up_broyden" to be consistent with the python script previously used.
     */
    private FieldTable(double[][][] table, double[][] values, boolean transposedRequiresMinus0, String name) {
        boolean validInput=true;
        validInput = validInput && values[0].length==table.length;
        for (int z=0;z<table.length && validInput;z++){
            validInput = validInput && values[1].length==table[z].length;
            for (int y=0;y<table[z].length && validInput;y++){
                validInput = validInput && values[2].length==table[z][y].length;
            }
        }
        if (!validInput)    throw new IllegalArgumentException("Invalid FieldTable arguments! mismatch between B values and table values");

        this.table = table;
        this.values = values;
        this.transposedRequiresMinus = transposedRequiresMinus0;
        this.name=name;
        this.dj=new int[]{Math.max(1,(int)Math.pow(values[0].length,0.25)),
                Math.max(1,(int)Math.pow(values[1].length,0.25)),
                Math.max(1,(int)Math.pow(values[2].length,0.25))};

        //this.cor = new boolean[3];
    }

    public boolean isTransposedRequiresMinus() {
        return transposedRequiresMinus;
    }

    public String getName() {
        return name;
    }

    @Override
    public String toString() {
        String ret="";
        ret = ret + "# Bz: [" + Arrays.toString(values[0]) + "]\n";
        ret = ret + "# By: [" + Arrays.toString(values[1]) + "]\n";
        ret = ret + "# Bx: [" + Arrays.toString(values[2]) + "]\n";

        for (int z=0;z<this.table.length;z++){
            ret = ret + "# New Bz slice\n";
            for (int y=0;y<this.table[z].length;y++){
                ret = ret + Arrays.toString(table[z][y]) +"\n";
            }
        }
        return ret;
    }


    /**
     * Bilinear interpolation. see picture at http://supercomputingblog.com/graphics/coding-bilinear-interpolation/
     * @return Intermediate value withing square based on the values at the corners
     */
    public static double bilinearInterpolation(double x1, double x2, double x, double y1, double y2, double y, double Q11, double Q12, double Q21, double Q22) {
        if (Q11==Q22){
            return Q11;
        }else{

            double R1 = ((x2 - x)/(x2 - x1))*Q11 + ((x - x1)/(x2 - x1))*Q21;
            double R2 = ((x2 - x)/(x2 - x1))*Q12 + ((x - x1)/(x2 - x1))*Q22;

            double ret = ((y2 - y)/(y2 - y1))*R1 + ((y - y1)/(y2 - y1))*R2;


            // since bilinear interpolation is now only used by trilinearInterpolation the verifications are done there.
            /*
            if ((ret > Q11 && ret > Q22 && ret > Q21 && ret > Q12)){
                System.err.println("error in bilinear interpolation: ret="+ret+" Q11="+Q11+" Q22="+Q22+" Q21="+Q21+" Q12="+Q12);
                return Math.max(Math.max(Math.max(Q11, Q12), Q21),Q22);	// return the maximum value
                //throw new RuntimeException("error in bilinear interpolation: ret="+ret+" Q11="+Q11+" Q22="+Q22);
            }else if((ret < Q11 && ret < Q22 && ret < Q21 && ret < Q12)){
                System.err.println("error in bilinear interpolation: ret="+ret+" Q11="+Q11+" Q22="+Q22+" Q21="+Q21+" Q12="+Q12);
                return Math.min(Math.min(Math.min(Q11, Q12), Q21),Q22);	// return the minimum value
            }else{
                return ret;
            }
            */
            return ret;
        }
    }

    /**
     * Trilinear interpolation. Simply performs linear interpolation between 2 values obtained with {@link FieldTable#bilinearInterpolation}
     * @return Intermediate value withing cube based on the values at the corners
     */
    public static double trilinearInterpolation(double x1, double x2, double x, double y1, double y2, double y, double z1, double z2, double z, double Q111, double Q121, double Q211, double Q221, double Q112, double Q122, double Q212, double Q222){
        if (Q111==Q222){
            return Q111;
        }else {
            double triy1 = bilinearInterpolation(x1, x2, x, y1, y2, y, Q111, Q121, Q211, Q221);
            double triy2 = bilinearInterpolation(x1, x2, x, y1, y2, y, Q112, Q122, Q212, Q222);

            double ret = triy1 + (z - z1) * (triy2 - triy1) / (z2 - z1);


            // verification is turned off for better performance
            // make sure return result is within the received cube
            /*
            double[] arr = new double[]{Q111, Q121, Q211, Q221, Q112, Q122, Q212, Q222};
            double max=arr[0], min=arr[0];
            for (int i=1;i<arr.length;i++){
                if (arr[i]>max) max=arr[i];
                if (arr[i]<min) min=arr[i];
            }

            if (ret>max){
                System.err.println("error in trilinear interpolation: ret="+ret+" , "+Arrays.toString(arr));
                return max;	// return the maximum value
                //throw new RuntimeException("error in bilinear interpolation: ret="+ret+" Q11="+Q11+" Q22="+Q22);
            }else if(ret<min){
                System.err.println("error in trilinear interpolation: ret="+ret+" , " + Arrays.toString(arr));
                return min;	// return the minimum value
            }else{
                return ret;
            }
            */
            return ret;
        }
    }

    // binary search. returns indices of the 2 closes values in a length-2 array: {lower,higher}
    // if the value is outside the bounds of the table a, returns null

    /**
     * Binary search. Returns indices of the 2 closest values in a length-2 array: {lower,higher}
     * If the value is outside the bounds of the table a, returns null
     * @param value value to find
     * @param a array in which to find the value
     * @return Indices of the 2 closest values in a length-2 array: {lower,higher}, or <code>null</code> if the value is outside the bounds of the received table.
     * @see FieldTable#hunt
     */
    public static int[] binarySearch(double value, double[] a){
        if(value < a[0] || value > a[a.length-1]) {
            // value is not in the table, so return null
            return null;
        }

        int lo = 0;
        int hi = a.length - 1;

        while (lo <= hi) {
            int mid = (hi + lo) / 2;

            if (value < a[mid]) {
                hi = mid - 1;
            } else if (value > a[mid]) {
                lo = mid + 1;
            } else {
                return new int[]{mid,mid};
            }
        }

        // at this point, lo == hi + 1
        return new int[]{hi,lo};
    }

    /**
     * Given a value, and an array, a, return an array int[]{high, low} such that value is between a[low] and a[high].
     * Altervative to {@link #binarySearch(double, double[]) binarySearch} that looks near the received startIndex.
     * Useful when subsequent calls should have similar values.
     * @param value value to find
     * @param a array in which to find the value
     * @param startIndex index from which the search sould start
     * @return Indices of the 2 closest values in a length-2 array: {lower,higher}, or <code>null</code> if the value is outside the bounds of the received table.
     * @see #binarySearch(double, double[])
     */
    public static int[] hunt(double value, double[] a, int startIndex){
        // Given a value, value, and an array, a, return an array int[]{high, low} such that value is
        // between a[low] and a[high]. The values in a must be monotonically increasing.
        // The returned value is not less than 0, nor greater than n-1. If the value is above a[n-1] or below a[0],
        // null is returned.
        // This function, in contrast with binarySearch, assumes subsequent calls have similar values, so the search
        // starts by looking near the previous found int[]{high, low}.

        if(value < a[0] || value > a[a.length-1]) {
            // value is not in the table, so return null
            return null;
        }

        int inc = 1;
        int lo = startIndex, hi;

        if (lo < 0 || lo > a.length-1) {
            lo=0;
            hi=a.length-1;
        } else {
            if (value >= a[lo]){    // hunt up
                for(;;){
                    hi = lo+inc;
                    if (hi >= a.length) { hi = a.length-1; break; }     // off end of table
                    else if (value < a[hi]) break;                      // found bracket
                    else {
                        lo = hi;
                        inc += inc;
                    }
                }
            } else {    // hunt down
                hi = lo;
                for (;;){
                    lo = lo - inc;
                    if (lo<=0) {lo=0; break;}                           // off end of table
                    else if (value >= a[lo]) break;                     // found bracket
                    else{
                        hi=lo;
                        inc += inc;
                    }
                }
            }
        }

        while (hi-lo > 1) {
            int mid = (hi + lo) >> 1;

            if (value >= a[mid]) {
                lo = mid;
            } else {
                hi = mid;
            }
        }


        return new int[]{lo,lo+1};
    }

    @Deprecated
    /**
     * @deprecated replaced by {@link #manualCalcValue} that uses ojAlgo to calc eigenenergies/magnetic-moments when necessary.
     */
    public static double pythonCalcValue(double Bx, double By, double Bz, String name){
        double ret = 0;
        try{
            ProcessBuilder pb = new ProcessBuilder("python3.6","single_ion_local_field_broyden.py",""+Bx,""+By,""+Bz,name);
            Process p = pb.start();

            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            ret = Double.parseDouble(in.readLine());
        }catch(Exception e){
            System.err.println("error running python: " + e);
        }

        return ret;
    }

    /** manually (build and diagonalize the hamiltonian) calculate the magnetic moment or the eigenenergy (of a spin up)
     * @param Bx Bx field
     * @param By By field
     * @param Bz Bz field
     * @param magneticMoment if true, calculate the magnetic moment. if false, calculate the eigenenergy.
     * @return Eigenenergy or magnetic moment of a Ho ion under the received magnetic field and crystal field potential.
     */
    public static double manualCalcValue(double Bx, double By, double Bz, boolean magneticMoment){
        if (!magneticMoment)
            return CrystalField.getEnergy(Bx, By, Bz);
        else
            return CrystalField.getMagneticMoment(Bx, By, Bz);
    }

    /**
     * Looks for the indices of the received Bx, By & Bz. Uses either {@link #hunt} or {@link #binarySearch(double, double[])} depending on the
     * received indicator.
     * @param bx Bx field
     * @param by By field
     * @param bz Bz field
     * @param prevBIndices
     * @return 3 by 2 int array of closest indices to the received point. The formatting is:
     * {{Bz_below, Bz_above},
     *  {By_below, By_above},
     *  {Bx_below, Bx_above}}
     */
    public int[][] searchForIndices(double bx, double by, double bz, int[][] prevBIndices) {
        int[] closestBxIndices, closestByIndices, closestBzIndices;
        if (prevBIndices==null) {
            // Bx
            closestBxIndices = binarySearch(bx, values[2]);    // look for bx in the Bx(2) values array
            // By
            closestByIndices = binarySearch(by, values[1]);    // look for by in the By(1) values array
            // Bz
            closestBzIndices = binarySearch(bz, values[0]);    // look for bz in the Bz(0) values array
        }else{
            // Bx
            if (prevBIndices[1][2]!=0) closestBxIndices = hunt(bx, values[2], prevBIndices[0][2]);    // look for bx in the Bx(2) values array
            else closestBxIndices = binarySearch(bx, values[2]);

            // By
            if (prevBIndices[1][1]!=0) closestByIndices = hunt(by, values[1], prevBIndices[0][1]);    // look for by in the By(1) values array
            else closestByIndices = binarySearch(by, values[1]);

            // Bz
            if (prevBIndices[1][0]!=0) closestBzIndices = hunt(bz, values[0], prevBIndices[0][0]);    // look for bz in the Bz(0) values array
            else closestBzIndices = binarySearch(bz, values[0]);
        }
        return new int[][]{closestBzIndices, closestByIndices, closestBxIndices};
    }

    /**
     * The main function that retrieves values from the table.
     * Given Bx, By & Bz, returns the interpolated value.
     * Also if prevBIndices is received, updates its elements to the closest indices of this search.
     * @param bx Bx value
     * @param by By value
     * @param bz Bz value
     * @param prevBIndices Array of the (lower) indices where the previous field was found. Could be useful if the
     *                    current search is expected to be close to that. Could be <code>null</code>, and then full binary-search is performed.
     * @return Interpolated value
     * @throws IndexOutOfBoundsException if the received bx, by or bz are outside the bounds of the table.
     */
    public double findInTable(double bx, double by, double bz, int[][] prevBIndices) throws IndexOutOfBoundsException {
        double ret;
        int [][] closestBIndices = searchForIndices(bx,by,bz,prevBIndices);
        int[] closestBxIndices=closestBIndices[2], closestByIndices=closestBIndices[1], closestBzIndices=closestBIndices[0];

        if (closestBxIndices != null && closestByIndices != null && closestBzIndices != null) {
            //table contains the energy drop due to the local transverse field
            // we interpolate according to the 8 nearest pre-calculated values
            ret = trilinearInterpolation(values[2][closestBxIndices[0]], values[2][closestBxIndices[1]], bx, values[1][closestByIndices[0]], values[1][closestByIndices[1]], by, values[0][closestBzIndices[0]], values[0][closestBzIndices[1]], bz,
                    table[closestBzIndices[0]][closestByIndices[0]][closestBxIndices[0]], table[closestBzIndices[0]][closestByIndices[1]][closestBxIndices[0]], table[closestBzIndices[0]][closestByIndices[0]][closestBxIndices[1]], table[closestBzIndices[0]][closestByIndices[1]][closestBxIndices[1]],
                    table[closestBzIndices[1]][closestByIndices[0]][closestBxIndices[0]], table[closestBzIndices[1]][closestByIndices[1]][closestBxIndices[0]], table[closestBzIndices[1]][closestByIndices[0]][closestBxIndices[1]], table[closestBzIndices[1]][closestByIndices[1]][closestBxIndices[1]]);
            if (prevBIndices!=null) {
                prevBIndices[1][0] = Math.abs(closestBzIndices[0] - prevBIndices[0][0]) < dj[0] ? 1 : 0;
                prevBIndices[1][1] = Math.abs(closestByIndices[0] - prevBIndices[0][1]) < dj[1] ? 1 : 0;
                prevBIndices[1][2] = Math.abs(closestBxIndices[0] - prevBIndices[0][2]) < dj[2] ? 1 : 0;
                prevBIndices[0][0]=closestBzIndices[0];
                prevBIndices[0][1]=closestByIndices[0];
                prevBIndices[0][2]=closestBxIndices[0];
            }
        } else {
            throw new IndexOutOfBoundsException("value out of bounds of the received table");
        }

        return ret;
    }

    /**
     * Get the partial derivative of the (linear) interpolated function. diff_by = 0,1,2 which means partial differentiation by z,y,x.
     * @param diff_by Partial differentiation by 0,1,2 which means by z,y,x. (in this order)
     * @param b_I Array of the magnetic field at which the partial derivative should be performed. format is {bz,by,bx}
     * @param s spin
     * @param prevBIndices same as in {@link #findInTable(double, double, double, int[][])}
     * @return partial derivative of the (linear) interpolated function at point b_I.
     */
    public double getDerivative(int diff_by, double[] b_I, double s, int[][] prevBIndices) {

        // copy b_I[] into b[] so that it is not ruined
        double[] b = new double[b_I.length];
        b[0]=b_I[0]*s;  // if s<0 it's the same as if bz -> -bz
        for (int i=1;i<b_I.length;i++) b[i]=b_I[i];

        double ret;

        int[][] closestBIndices = searchForIndices(b[2],b[1],b[0],prevBIndices);
        if (closestBIndices[0] != null && closestBIndices[1] != null && closestBIndices[2] != null) {

            int const1, const2;
            if (diff_by == 0) {
                const1 = 1;
                const2 = 2;
            } else if (diff_by == 1) {
                const1 = 2;
                const2 = 0;
            } else if (diff_by == 2) {
                const1 = 0;
                const2 = 1;
            } else {
                const1 = -1;
                const2 = -1;
                System.err.println("diff_by only accepts integer between 0 and 2");
                System.exit(1);
            }
            double triy1 = bilinearInterpolation(values[const2][closestBIndices[const2][0]], values[const2][closestBIndices[const2][1]], b[const2], values[const1][closestBIndices[const1][0]], values[const1][closestBIndices[const1][1]], b[const1],
                    table[closestBIndices[0][0]][closestBIndices[1][0]][closestBIndices[2][0]],table[closestBIndices[0][booleanToInt(const1==0)]][closestBIndices[1][booleanToInt(const1==1)]][closestBIndices[2][booleanToInt(const1==2)]], table[closestBIndices[0][booleanToInt(const2==0)]][closestBIndices[1][booleanToInt(const2==1)]][closestBIndices[2][booleanToInt(const2==2)]], table[closestBIndices[0][booleanToInt(const1==0 || const2==0)]][closestBIndices[1][booleanToInt(const1==1 || const2==1)]][closestBIndices[2][booleanToInt(const1==2 || const2==2)]]);
            double triy2 = bilinearInterpolation(values[const2][closestBIndices[const2][0]], values[const2][closestBIndices[const2][1]], b[const2], values[const1][closestBIndices[const1][0]], values[const1][closestBIndices[const1][1]], b[const1],
                    table[closestBIndices[0][booleanToInt(diff_by==0)]][closestBIndices[1][booleanToInt(diff_by==1)]][closestBIndices[2][booleanToInt(diff_by==2)]],table[closestBIndices[0][booleanToInt(const1==0 || diff_by==0)]][closestBIndices[1][booleanToInt(const1==1 || diff_by==1)]][closestBIndices[2][booleanToInt(const1==2 || diff_by==2)]], table[closestBIndices[0][booleanToInt(const2==0 || diff_by==0)]][closestBIndices[1][booleanToInt(const2==1 || diff_by==1)]][closestBIndices[2][booleanToInt(const2==2 || diff_by==2)]], table[closestBIndices[0][1]][closestBIndices[1][1]][closestBIndices[2][1]]);

            ret = (triy2 - triy1) / (values[diff_by][closestBIndices[diff_by][1]] - values[diff_by][closestBIndices[diff_by][0]]);
        } else {
            throw new IndexOutOfBoundsException("value out of bounds of the received table. Cannot calculate derivative.");
        }

        if (transposedRequiresMinus){
            if (diff_by!=0) ret = Math.signum(s)*ret;
        }


        return ret;
    }

    /**
     * Convert true to 1 and false to 0.
     * @param value value to be converted
     * @return 1 if value is true and 0 if value is false
     */
    public static int booleanToInt(boolean value) {
        // Convert true to 1 and false to 0.
        return value ? 1 : 0;
    }

    /**
     * get value from this table (by interpolation)
     * @param bx Bx field
     * @param by By field
     * @param bz Bz field
     * @param s spin
     * @param returnSuccessStatus return the success status (success is when manual search is not performed). Used to monitor how many times a manual calculation is performed.
     *                            The value is not actually used, this is only to overload the function.
     * @param prevBIndices see {@link #findInTable(double, double, double, int[][])}
     * @return Interpolated value from this table
     */
    public returnDoubleAndStatus getValue(double bx, double by, double bz, double s, boolean returnSuccessStatus, int[][] prevBIndices) {
        // if s<0 it's the same as if bz -> -bz
        bz=bz*s;
        double ret;
        boolean status=true;
        try{
            ret = findInTable(bx,by,bz,prevBIndices);
        } catch(IndexOutOfBoundsException e){
            //bx, by or bz are outside the range of the table
            // so the the value calculated directly
            // this takes a lot of time and is not recommended as a frequent solution

            // the call and name "magnetic_moment_up_broyden" is kept to preserve compatibility with the python script
            ret = manualCalcValue(bx, by, bz, this.name.substring(0, name.length() - 9).equals("magnetic_moment_up_broyden"));
            status=false;
            //System.err.println("field was not in table. " + name.substring(0, name.length() - 20) + " had to be calculated manually using python. Bx=" + bx + " By=" + by + " Bz=" + bz + "|| the value calculated is: " + ret);
        }

        if (transposedRequiresMinus)    ret = Math.signum(s)*ret;

        return new returnDoubleAndStatus(ret, status);

    }

    /**
     * @see #getValue(double, double, double, double, boolean, int[][])
     */
    public double getValue(double bx, double by, double bz, double s, int[][] prevBIndices) {
        return getValue(bx,by,bz,s,false, prevBIndices).getValue();
    }

    /**
     * @see #getValue(double, double, double, double, boolean, int[][])
     */
    public double getValue(double bx, double by, double bz, double s) {
        return getValue(bx,by,bz,s,false, null).getValue();
    }

    /**
     * @see #getValue(double, double, double, double, boolean, int[][])
     */
    public returnDoubleAndStatus getValue(double bx, double by, double bz, double s, boolean returnSuccessStatus) {
        return getValue(bx,by,bz,s,returnSuccessStatus, null);
    }

    // just for testing getDerivative
    public static double[] getNumericalDerivative(double[] b, double s, FieldTable momentTable){

        int numOfSteps=0;

        double bx=b[2], by=b[1], bz=b[0];
        double[] derivative = new double[3];
        double[] prev_derivative = new double[3];
        double delta=0.002;
        double tol=1.0e-10;

        prev_derivative[2] = (momentTable.getValue(bx+delta, by, bz, s)-momentTable.getValue(bx-delta, by, bz, s))/(2*delta);
        prev_derivative[1] = (momentTable.getValue(bx, by+delta, bz, s)-momentTable.getValue(bx, by-delta, bz, s))/(2*delta);
        prev_derivative[0] = (momentTable.getValue(bx, by, bz+delta, s)-momentTable.getValue(bx, by, bz-delta, s))/(2*delta);
        while (numOfSteps<3){
            delta*=0.5;
            derivative[2] = (momentTable.getValue(bx+delta, by, bz, s)-momentTable.getValue(bx-delta, by, bz, s))/(2*delta);
            derivative[1] = (momentTable.getValue(bx, by+delta, bz, s)-momentTable.getValue(bx, by-delta, bz, s))/(2*delta);
            derivative[0] = (momentTable.getValue(bx, by, bz+delta, s)-momentTable.getValue(bx, by, bz-delta, s))/(2*delta);
            if (Math.abs(derivative[0]-prev_derivative[0])<tol && Math.abs(derivative[1]-prev_derivative[1])<tol && Math.abs(derivative[2]-prev_derivative[2])<tol){
                numOfSteps++;
            }else{
                System.arraycopy(derivative, 0, prev_derivative, 0, derivative.length);
            }
            //System.out.println(derivative[2]+","+derivative[1]+","+derivative[0]+" delta="+delta);
        }

        return derivative;
    }

    // testing getDerivative
    public static void testDerivative(FieldTable momentTable, double extBx){
        Random rnd = new Random();
        boolean allSuccess=true;
        boolean[] fails = new boolean[5000];
        for (int i=0;i<5000; i++) {
            double bx = extBx + 4 * rnd.nextDouble() - 2, by = 4 * rnd.nextDouble() - 2, bz = 4 * rnd.nextDouble() - 2;
            double[] b = new double[]{bz,by,bx};
            int s = rnd.nextBoolean() ? 1 : -1;

            double[] realDerivative = getNumericalDerivative(b,s, momentTable);
            double[] testedDerivative = new double[]{momentTable.getDerivative(0, b, s, null),
                    momentTable.getDerivative(1, b, s, null),
                    momentTable.getDerivative(2, b, s, null)};

            if (realDerivative[0]==0 && realDerivative[1]==0 && realDerivative[2]==0){
                i--;
            }else {
                for (int j = 0; j < realDerivative.length; j++) {

                    if ((Math.abs(realDerivative[j] - testedDerivative[j]) > 1e-8)) {
                        System.out.println("("+i+")");
                        System.out.println("real: " + Arrays.toString(realDerivative));
                        System.out.println("tested: " + Arrays.toString(testedDerivative));
                        System.out.println(s + ", df(" + bx + "," + by + "," + bz + ")/d" + j + " : " + (Math.abs(realDerivative[j] - testedDerivative[j]) < 1e-8));

                        allSuccess = false;
                        fails[i]=true;
                        System.out.println("error!");
                        //System.exit(1);
                    }
                }
            }

        }
        System.out.println("Success: " + (allSuccess ? "Yes" : "No"));
        System.out.println("fails at: ");
        int count=0;
        for (int i=0;i<fails.length;i++){
            if (fails[i]==true){
                System.out.print(i+", ");
                count++;
            }
        }
        System.out.println("success rate: "+(100.0-100.0*count/fails.length)+"%");
    }

}

