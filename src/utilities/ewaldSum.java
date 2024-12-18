package utilities;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.special.Erf;
import simulation.montecarlo.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Properties;

import static simulation.montecarlo.Main.makeDir;

public class ewaldSum {

	/**
	 * Convergence tests for initial verification of ewald summation in 3D. this compares ewald results with direct calculation
	 * @param arr array of spins
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y directions
	 * @param alpha Ewald parameter
	 * @param N total number of spins
	 * @param real_cutoff cutoff in real space
	 * @param k_cutoff cutoff in reciprocal space
	 * @param direct_cutoff cutoff for direct calculation
	 * @param ellipsoidRatio initial ratio of the ellipsoid which contains all unit cells included in the direct calculation
	 */
	public static void convergenceTests(singleSpin[] arr, int Lz, int Lx, double alpha, int N, int real_cutoff, int k_cutoff, int direct_cutoff, double ellipsoidRatio){
		// get 2 close spins
		int i=(int)N/2;
		int j=i+1;

		// calculate their interaction using the Ewald method
		double ewald_result = calcSum3D(arr[i],arr[j],Lz,Lx,alpha,real_cutoff,k_cutoff)[0];
		System.out.println("ratio\tewald\tdirect(1)\tdirect(2)\tdirect(3)\tdirect(4)\tdirect(5)\tdirect(6)\tdirect(7)\tdirect(8)\tdirect(9)\tdirect(10)");

		// calculate their interaction by summing directly all copies within an ellipsoid which becomes narrower (closer to a needle shape)
		double[][] direct_result;
		for (double eR=1;eR<=ellipsoidRatio;eR=eR+(ellipsoidRatio/100)) {
			direct_result = realCalcSum3DE(arr[i], arr[j], Lz, Lx, Lx, Lz*direct_cutoff, eR);
			System.out.print(eR + "\t" + ewald_result);
			for (int index=0;index<direct_result.length;index++){
				System.out.print("\t" + direct_result[index][0]);
			}
			System.out.println();
		}
	}

	/**
	 * More advanced convergence tests than {@link #convergenceTests(singleSpin[], int, int, double, int, int, int, int, double)}.
	 * This is meant to find optimal values for alpha, real_cutoff, k_cutoff
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y directions
	 * @param real_cutoff cutoff in real space
	 * @param max_k_cutoff maximal cutoff in reciprocal space
	 */
	public static void convergenceTests2(int Lz, int Lx, int real_cutoff, int max_k_cutoff){
		// get 2 close spins.
		// spins i,j,k is somewhere in the middle of the system
		int i=(int)(Lx/2);
		int j=(int)(Lx/2);
		int k=(int)(Lz/2);
		int testSpin=(i)*Lx*Lz*Constants.num_in_cell+(j)*Lz*Constants.num_in_cell+(k)*Constants.num_in_cell+0;
		int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;

		// get the neighbors of testSpin
		neighbor1=i*Lx*Lz*4+j*Lz*4+k*4+1;
		if (i==0)
			neighbor2=(Lx-1)*Lx*Lz*4+j*Lz*4+k*4+1;
		else
			neighbor2=(i-1)*Lx*Lz*4+j*Lz*4+k*4+1;
		if (k==0)
			neighbor3=i*Lx*Lz*4+j*Lz*4+(Lz-1)*4+3;
		else
			neighbor3=i*Lx*Lz*4+j*Lz*4+(k-1)*4+3;

		if (j==0 && k==0)
			neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(Lz-1)*4+3;
		else if(j==0){
			neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(k-1)*4+3;
		}
		else if(k==0){
			neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(Lz-1)*4+3;
		}else{
			neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(k-1)*4+3;
		}
		double alpha;

		final double[][][] intTable = new double[3][Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz];
		final double[][] exchangeIntTable = new double[Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz];	// all zeros


		System.out.println("alpha\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12");
		// sweep over values of alpha, the Ewald parameter
		// look for convergence of the system's energy in a
		// checkerboard configuration
		for (alpha = 0.1 / Lz; alpha < 5.0 / Lz; alpha += 0.1 / Lz) {
			System.out.print(alpha*Lz);
			for(int k_cutoff=1;k_cutoff<=max_k_cutoff;k_cutoff++) {
				Lattice lattice = new Lattice(Lx, Lz, 0.0, 0.0, false, 5.51,intTable, exchangeIntTable, null, null, null, null);
				fillIntTable(lattice.getArray(), Lz, Lx, alpha, real_cutoff, k_cutoff, intTable);
				lattice.checkerBoard();
				lattice.updateAllLocalFields();
				System.out.print("\t" + (calcEnergy(lattice)/(Lx*Lx*Lz*Constants.num_in_cell)));
			}
			System.out.println();
		}
	}

	/**
	 * Most advanced convergence tests.
	 * The output is used by the jupyter notebooks in /ewald_convergence when called by
	 * {@link #multipleSizeTestsWithvariableKcutoff(int[], int, int[])} and {@link #multipleSizeTestsWithvariableRealcutoff(int[], int[], int)}.
	 * Meant to find optimal values for alpha, real_cutoff, k_cutoff
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y directions
	 * @param real_cutoff cutoff in real space
	 * @param k_cutoff cutoff in reciprocal space
	 */
	public static void convergenceTests3(int Lz, int Lx, int real_cutoff, int k_cutoff){
		// get a test-spin in the middle of the system
		int i=(int)(Lx/2);
		int j=(int)(Lx/2);
		int k=(int)(Lz/2);
		int testSpin=(i)*Lx*Lz*Constants.num_in_cell+(j)*Lz*Constants.num_in_cell+(k)*Constants.num_in_cell+0;
		int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;

		// get the neighbors of the test spin
		if (System.getProperty("system").equals("LiHoF4")) {
			// Nearest neighbors for LiHoF4
			neighbor1 = i * Lx * Lz * Constants.num_in_cell + j * Lz * Constants.num_in_cell + k * Constants.num_in_cell + 1;
			neighbor2 = ((Lx + i - 1) % Lx) * Lx * Lz * Constants.num_in_cell + j * Lz * Constants.num_in_cell + k * Constants.num_in_cell + 1;
			neighbor3 = i * Lx * Lz * Constants.num_in_cell + ((Lx + j - 1) % Lx) * Lz * Constants.num_in_cell + ((Lz + k - 1) % Lz) * Constants.num_in_cell + 3;
			neighbor4 = i * Lx * Lz * Constants.num_in_cell + j * Lz * Constants.num_in_cell + ((Lz + k - 1) % Lz) * Constants.num_in_cell + 3;
		}else if (System.getProperty("system").equals("Fe8")){
			// Nearest neighbors for Fe8
			neighbor1 = i * Lx * Lz * Constants.num_in_cell + j * Lz * Constants.num_in_cell + k * Constants.num_in_cell;
			neighbor2 = ((Lx + i - 1) % Lx) * Lx * Lz * Constants.num_in_cell + j * Lz * Constants.num_in_cell + k * Constants.num_in_cell;
			neighbor3 = i * Lx * Lz * Constants.num_in_cell + ((Lx + j - 1) % Lx) * Lz * Constants.num_in_cell + ((Lz + k - 1) % Lz) * Constants.num_in_cell;
			neighbor4 = i * Lx * Lz * Constants.num_in_cell + j * Lz * Constants.num_in_cell + ((Lz + k - 1) % Lz) * Constants.num_in_cell;
		}
		double alpha;

		final double[][][] intTable = new double[3][Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz];
		final double[][] exchangeIntTable = new double[Constants.num_in_cell*Lx*Lx*Lz][Constants.num_in_cell*Lx*Lx*Lz];	// all zeros

		Lattice lattice = new Lattice(Lx, Lz, 0.0, 0.0, false, 5.51, intTable, exchangeIntTable, null, null, null, null);
		singleSpin[] arr = lattice.getArray();
		System.out.println("#real_cutoff="+real_cutoff);
		System.out.println("#reciprocal_cutoff="+k_cutoff);

		// print out the neighbors and their locations relative to the test spin
		DecimalFormat df = new DecimalFormat("0.0000000000");
		System.out.println("#r0 - r0 : " + Arrays.toString(arr[testSpin].distance(arr[testSpin], Lz, Lx)));
		System.out.println("#r0 - r1 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor1], Lz, Lx)));
		System.out.println("#r0 - r2 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor2], Lz, Lx)));
		System.out.println("#r0 - r3 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor3], Lz, Lx)));
		System.out.println("#r0 - r4 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor4], Lz, Lx)));

		// print out all interactions with each of the neighbors and the self-interaction
		System.out.println("alpha\t0x\t0y\t0z\t1x\t1y\t1z\t2x\t2y\t2z\t3x\t3y\t3z\t4x\t4y\t4z");
		for (alpha = 0.1 / Lz; alpha < 5.0 / Lz; alpha += 0.1 / Lz) {
			System.out.print(df.format(alpha*Lz)+"\t");

			double[] interaction0=ewaldSum.calcSum3D(arr[testSpin], arr[testSpin], Lz, Lx, alpha, real_cutoff, k_cutoff);
			System.out.print(df.format(interaction0[0])+"\t"+df.format(interaction0[1])+"\t"+df.format(interaction0[2])+"\t");
			double[] interaction1=ewaldSum.calcSum3D(arr[testSpin], arr[neighbor1], Lz, Lx, alpha, real_cutoff, k_cutoff);
			System.out.print(df.format(interaction1[0])+"\t"+df.format(interaction1[1])+"\t"+df.format(interaction1[2])+"\t");
			double[] interaction2=ewaldSum.calcSum3D(arr[testSpin], arr[neighbor2], Lz, Lx, alpha, real_cutoff, k_cutoff);
			System.out.print(df.format(interaction2[0])+"\t"+df.format(interaction2[1])+"\t"+df.format(interaction2[2])+"\t");
			double[] interaction3=ewaldSum.calcSum3D(arr[testSpin], arr[neighbor3], Lz, Lx, alpha, real_cutoff, k_cutoff);
			System.out.print(df.format(interaction3[0])+"\t"+df.format(interaction3[1])+"\t"+df.format(interaction3[2])+"\t");
			double[] interaction4=ewaldSum.calcSum3D(arr[testSpin], arr[neighbor4], Lz, Lx, alpha, real_cutoff, k_cutoff);
			System.out.print(df.format(interaction4[0])+"\t"+df.format(interaction4[1])+"\t"+df.format(interaction4[2]));

			System.out.println();
		}
	}

	/**
	 * Print ewald convergence tests for multiple Lz's for analysis by the jupyter notebooks in /ewald_convergence
	 * @param Lz number of unit cells in the z direction
	 * @param real_cutoff constant real space cutoff
	 * @param k_cutoffs array of reciprocal space cutoffs
	 */
	public static void multipleSizeTestsWithvariableKcutoff(int[] Lz, int real_cutoff, int[] k_cutoffs){
		for (int i=0;i<Lz.length;i++) {
			try {
				makeDir("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\Documents-Tomer\\thesis\\ewald_convergence\\Fe8\\", "L=" + Lz[i]);
				for (int j=0;j<k_cutoffs.length;j++) {
					PrintStream fileOut = new PrintStream("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\Documents-Tomer\\thesis\\ewald_convergence\\Fe8\\L=" + Lz[i] + "\\" +k_cutoffs[j]+"_"+real_cutoff+".txt");
					System.setOut(fileOut);
					convergenceTests3(Lz[i],Lz[i],real_cutoff,k_cutoffs[j]);
				}
			} catch (FileNotFoundException e){
				e.printStackTrace();
			}
		}
	}

	/**
	 * Print ewald convergence tests for multiple Lz's for analysis by the jupyter notebooks in /ewald_convergence
	 * @param Lz number of unit cells in the z direction
	 * @param real_cutoffs array of real space cutoffs
	 * @param k_cutoff constant reciprocal space cutoff
	 */
	public static void multipleSizeTestsWithvariableRealcutoff(int[] Lz, int[] real_cutoffs, int k_cutoff){

		for (int i=0;i<Lz.length;i++) {
			try {
				makeDir("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\Documents-Tomer\\thesis\\ewald_convergence\\Fe8\\", "L=" + Lz[i]);
				for (int j=0;j<real_cutoffs.length;j++) {
					PrintStream fileOut = new PrintStream("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\Documents-Tomer\\thesis\\ewald_convergence\\Fe8\\L=" + Lz[i] + "\\" +k_cutoff+"_"+real_cutoffs[j]+".txt");
					System.setOut(fileOut);
					convergenceTests3(Lz[i],Lz[i],real_cutoffs[j],k_cutoff);
				}
			} catch (FileNotFoundException e){
				e.printStackTrace();
			}
		}
	}

	/**
	 * Calculates the energy of the system
	 * @param lattice the lattice
	 * @return the energy of the lattice
	 */
	public static double calcEnergy(Lattice lattice){
		double energy=0;
		singleSpin[] arr = lattice.getArray();
		for(int i=0;i<arr.length;i++){
			energy += CrystalField.getEnergy(arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz(), arr[i].getSpin());
		}
		return energy;
	}

	/**
	 * Create the tables of interaction using Ewald's method.
	 * Also used for testing the convergence with the given parameters when the
	 * appropriate lines are uncommented.
	 * @param args
	 */
	public static void main(String[] args){
		// read parameters from file:
		Properties params = GetParamValues.getParams();
		int real_cutoff=GetParamValues.getIntParam(params, "real_cutoff");
		int k_cutoff=GetParamValues.getIntParam(params, "k_cutoff");
		double alpha = GetParamValues.getDoubleParam(params, "alpha");
		params=null;
		
		int Lx=0;	// lattice x-y size
		int Lz=0;	// lattice z size
		k_cutoff=0;
		real_cutoff=0;
        // for convergence tests
        int direct_cutoff = 100;
		double ellipsoidRatio = 50;

        // get lattice parameters as command line arguments
        try {
        	Lx = Integer.parseInt(args[0]);
        	Lz = Integer.parseInt(args[1]);
			k_cutoff = Integer.parseInt(args[2]);
			real_cutoff = Integer.parseInt(args[3]);
        	// for convergence tests using convergenceTests()
        	// real_cutoff = Integer.parseInt(args[2]);
			// k_cutoff = Integer.parseInt(args[3]);
			// direct_cutoff = Integer.parseInt(args[4]);
			// ellipsoidRatio = Double.parseDouble(args[5]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.err.println("ArrayIndexOutOfBoundsException caught " + e.toString());
        }
        catch (NumberFormatException e){
			System.err.println("Argument not a valid number " + e.toString());
		}

//      uncomment the next 4 lines to run convergence tests to find optimal real- and reciprocal-space values
//		multipleSizeTestsWithvariableKcutoff(new int[]{3,4,5,6,7,8,9,10},5,new int[]{2,3,4,5,6,7,8,9,10,11,12});
//		multipleSizeTestsWithvariableRealcutoff(new int[]{3,4,5,6,7,8,9,10},new int[]{2,3,4,5,6,7,8,9,10,11,12},2);
//		convergenceTests3(Lz, Lx, real_cutoff, k_cutoff);
//      System.exit(0);

		// if convergence tests were not run, do the Ewald calculation
        try (BufferedWriter out = new BufferedWriter(new FileWriter(System.getProperty("system") + File.separator + "data" + File.separator + "interactions" + File.separator + "intTable_"+Lx+"_"+Lz+".txt"))){
            
            //print lattice sizes
            out.write("Lx="+Lx);
            out.newLine();
            out.write("Lz="+Lz);
            out.newLine();
            out.write("real_cutoff="+real_cutoff);
            out.newLine();
            out.write("k_cutoff="+k_cutoff);
            out.newLine();
            out.write("alpha="+alpha);
            out.newLine();

            // generate a lattice object, just to get the singleSpin array from it
            singleSpin[] arr = (new Lattice(Lx, Lz, 0.0, 0.0, false, 5, null, null, null, null, null, null)).getArray();

            // fill interactions table and print it to the file 'out'
			if (System.getProperty("system").equals("LiHoF4")) {
				// in LiHoF4 the c is parallel to the z axis
				fillIntTable(arr, Lz, Lx, alpha / (Constants.c*Lz), real_cutoff, k_cutoff, out);
			}else if (System.getProperty("system").equals("Fe8")){
				// in Fe8 the a is parallel to the z axis
				fillIntTable(arr, Lz, Lx, alpha / (Constants.a*Lz), real_cutoff, k_cutoff, out);
			} else {
				throw new RuntimeException("Cannot create Ewald table. Illegal system name given.");
			}

        }
        catch (IOException e) { System.out.println("bad file"); }

        // other convergence tests
//		int N = Lx*Lx*Lz*Constants.num_in_cell;
//		singleSpin[] arr = GenerateLattice.generate_ising_lattice(Lx, Lz, 1, 0, null);
//		convergenceTests(arr,Lz,Lx,1.0/Lz,N,real_cutoff,k_cutoff,direct_cutoff, ellipsoidRatio);
//		convergenceTests2(arr,Lz,Lx,real_cutoff,k_cutoff);
    }

	/**
	 * Fills a given interactions table with dipolar interactions calculated using Ewald's method.
	 * This fills the interaction with the appropriate units and prefactor (hence it is not used for the simulation,
	 * where these are added separately)
	 * @param arr spin array
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y directions
	 * @param alpha Ewald parameter
	 * @param real_cutoff real-space cutoff
	 * @param k_cutoff reciprocal-space cutoff
	 * @param intTable interactions tables to fill
	 */
	public static void fillIntTable(singleSpin[] arr, int Lz, int Lx, double alpha, int real_cutoff, int k_cutoff, final double[][][] intTable){
		final double c = -Constants.mu_0*Constants.mu_B*Constants.g_L*0.25/Math.PI;	// coefficient dipolar spin-spin interaction. The minus sign is because
		for (int i=0;i<arr.length;i++){
			for (int j=i;j<arr.length;j++){
				double[] interaction = new double[3];
				interaction = ewaldSum.calcSum3D(arr[i], arr[j], Lz, Lx, alpha, real_cutoff, k_cutoff);

				for (int k=0;k<interaction.length;k++) {
					double value = c * interaction[k];
					value = Math.floor(value * 1.0e7) / 1.0e7;	// truncate at 7 digits
					intTable[k][i][j] = value;
					intTable[k][j][i] = value;
					if (k==2){  // for the zz interactions we multiply by half to avoid double counting.
						intTable[k][i][j] *= 0.5;
						intTable[k][j][i] *= 0.5;
					}
				}
			}
		}
	}

	/**
	 * Prints out the dipolar interactions calculated using Ewald's method.
	 * This fills the interaction with the bare dipolar interaction (it is the one used in the simulation)
	 * @param arr spin array
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y directions
	 * @param alpha Ewald parameter
	 * @param real_cutoff real-space cutoff
	 * @param k_cutoff reciprocal-space cutoff
	 * @param out the output stream where the calc'd interactions are saved
	 * @throws IOException
	 */
	public static void fillIntTable(singleSpin[] arr, int Lz, int Lx, double alpha, int real_cutoff, int k_cutoff, BufferedWriter out) throws IOException{
		DecimalFormat df = new DecimalFormat("0.0000000000");
		for (int i=0;i<arr.length;i++){
			// we only save half of the table since it is symmetric
			for (int j=i;j<arr.length;j++){
				double[] interaction = new double[3];
				interaction = ewaldSum.calcSum3D(arr[i], arr[j], Lz, Lx, alpha, real_cutoff, k_cutoff);
				// also print out the results to the console
				System.out.println("("+i+","+j+")" + " : " + interaction[0] + "," + interaction[1] + "," + interaction[2]);
				// print out to the dile
				out.write(df.format(interaction[0])+","+df.format(interaction[1])+","+df.format(interaction[2]));	// print 10 significant digits
				out.newLine();
			}	
		}
	}


	// call calcSum3D with Lx=Ly as is often the case

	/**
	 * Calls {@link #calcSum3D(singleSpin, singleSpin, int, int, double, int, int)} with Lx=Ly as is often the case
	 */
	public static double[] calcSum3D(singleSpin i, singleSpin j, int Lz, int Lx, double alpha, int realN_cutoff, int reciprocalN_cutoff){
		return calcSum3D(i, j, Lz, Lx, Lx, alpha, realN_cutoff, reciprocalN_cutoff);
	}

	/**
	 * The correct implementation of the Ewald summation method
	 * @param i spin object {@code i}
	 * @param j spin object {@code j}
	 * @param Lz number of unit cells in the z direction
	 * @param Ly number of unit cells in the y direction
	 * @param Lx number of unit cells in the x direction
	 * @param alpha Ewald parameter
	 * @param realN_cutoff real space cutoff
	 * @param reciprocalN_cutoff reciprocal space cutoff
	 * @return the zx, zy, zz interactions in a size 3 array: {zx, zy, zz}
	 * @see <a href="http://arxiv.org/abs/0912.3469" target="_top">arXiv:0912.3469 [cond-mat] - Appendix C</a>
	 */
	public static double[] calcSum3D(singleSpin i, singleSpin j, int Lz, int Ly, int Lx, double alpha, int realN_cutoff, int reciprocalN_cutoff){
		/*
		Calculates the dipolar interaction between two spin, i and j, using the Ewald summation method.
		See Appendix C in arXiv:0912.3469 [cond-mat] for a concise explanation.
		 */

		// init real and reciprocal sums
		double realSumZ=0, realSumX=0, realSumY=0;
		double reciprocalSumZ=0;
		double reciprocalSumY=0;
		double reciprocalSumX=0;

		// generate displacement vector between spins i and j
		Vector3D rij = i.getLocation(Lz,Ly).subtract(j.getLocation(Lz,Ly));

		// generate displacement vectors the size of the lattice in each of the primitive directions
		Vector3D latticeLz = Constants.primitiveLatticeVectors[2].scalarMultiply(Lz);
		Vector3D latticeLy = Constants.primitiveLatticeVectors[1].scalarMultiply(Ly);
		Vector3D latticeLx = Constants.primitiveLatticeVectors[0].scalarMultiply(Lx);

		double r, rx, ry, rz, B, C;
		int nx,ny,nz;

		// sum vectors in real space
		for (nx=-realN_cutoff;nx<=realN_cutoff;nx++){
			for (ny=-realN_cutoff;ny<=realN_cutoff;ny++){
				for (nz=-realN_cutoff;nz<=realN_cutoff;nz++){
					if (nx!=0 || ny!=0 || nz!=0 || i!=j){	// exclude self interaction (the prime in eq. (C1) in the referenced preprint)

						Vector3D rijCopy = new Vector3D(1, rij);	// distance between i and the copy of j that is located
																	// nx unit cells along primitive direction 0
																	// ny unit cells along primitive direction 1 and
																	// nz unit cells along primitive direction 2

						// displacement vector between spin i and the copy of spin j
						rijCopy=rijCopy.add(nz,latticeLz);
						rijCopy=rijCopy.add(ny,latticeLy);
						rijCopy=rijCopy.add(nx,latticeLx);

						// components of said vector
						rz=rijCopy.getZ();
						ry=rijCopy.getY();
						rx=rijCopy.getX();
						r = rijCopy.getNorm();

						// eq. (C12)
						B = (Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*Math.exp(-alpha*alpha*r*r))/(r*r*r);

						// eq. (C13)
						C = (3*Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*(3+2*alpha*alpha*r*r)*Math.exp(-alpha*alpha*r*r))/(Math.pow(r, 5));

						// eq. (C15)
						realSumZ += B - rz*rz*C;
						realSumX -= rx*rz*C;
						realSumY -= ry*rz*C;
					}
				}
			}
		}

		// reciprocal lattice primitive vectors
		Vector3D b1 =  Vector3D.crossProduct(latticeLy,latticeLz).scalarMultiply(2 * Math.PI / latticeLx.dotProduct(Vector3D.crossProduct(latticeLy,latticeLz)));
		Vector3D b2 =  Vector3D.crossProduct(latticeLz,latticeLx).scalarMultiply(2 * Math.PI / latticeLx.dotProduct(Vector3D.crossProduct(latticeLy,latticeLz)));
		Vector3D b3 =  Vector3D.crossProduct(latticeLx,latticeLy).scalarMultiply(2 * Math.PI / latticeLx.dotProduct(Vector3D.crossProduct(latticeLy,latticeLz)));

		for (int mx=-reciprocalN_cutoff;mx<=reciprocalN_cutoff;mx++){
			for (int my=-reciprocalN_cutoff;my<=reciprocalN_cutoff;my++) {
				for (int mz = 0; mz <= reciprocalN_cutoff; mz++) {
					// the following condition makes sure we only count half the (reciprocal) space which is what turns the exp to cos
					// also implements the condition G!=0
					if (mz > 0 || (mz==0 && mx>0) || (mz==0 && mx==0 && my>0)) {
						Vector3D G = Vector3D.ZERO;	// reciprocal lattice vector
						G=G.add(mz,b3);
						G=G.add(my,b2);
						G=G.add(mx,b1);

						double Gz=G.getZ();
						double Gy=G.getY();
						double Gx=G.getX();
						double G2 = G.getNormSq();

						double Gdotrij = G.dotProduct(rij);
						double cosGdotrij = Math.cos(Gdotrij);

						// eq. (C16)
						reciprocalSumX += 2*cosGdotrij*(((Gx * Gz) / G2) * Math.exp(-G2 / (4 * alpha * alpha)));
						reciprocalSumY += 2*cosGdotrij*(((Gy * Gz) / G2) * Math.exp(-G2 / (4 * alpha * alpha)));
						reciprocalSumZ += 2*cosGdotrij*(((Gz * Gz) / G2) * Math.exp(-G2 / (4 * alpha * alpha)));
					}
				}
			}
		}

		double[] ret = new double[3];

		// self interaction correction term to be omitted from the zz interaction for i==j, eq. (C17):
		double selfInteraction = i.getN()==j.getN() ? 4*Math.pow(alpha,3)/(3*Math.sqrt(Math.PI)) : 0;

		ret[0] = realSumX + (4*Math.PI/(latticeLx.dotProduct(Vector3D.crossProduct(latticeLy,latticeLz))))*reciprocalSumX;	//zx interaction
		ret[1] = realSumY + (4*Math.PI/(latticeLx.dotProduct(Vector3D.crossProduct(latticeLy,latticeLz))))*reciprocalSumY;	//zy interaction
		ret[2] = realSumZ + (4*Math.PI/(latticeLx.dotProduct(Vector3D.crossProduct(latticeLy,latticeLz))))*reciprocalSumZ - selfInteraction;	//zz interaction
		return ret;

	}


	// bottom line: this is the right implementation of ewald sum for this project
	// returns the zx, zy, zz interations in a size 3 array: {zx, zy, zz}

	/**
	 * Calculates the dipolar interaction between two spins using the Ewald summation method.
	 * Works only for Orthorhombic, Tetragonal and Cubic crystal systems.
	 * @deprecated use {@link #calcSum3D(singleSpin, singleSpin, int, int, int, double, int, int)}
	 */
	@Deprecated
	public static double[] calcSum3D2(singleSpin i, singleSpin j, int Lz, int Ly, int Lx, double alpha, int realN_cutoff, int reciprocalN_cutoff){
		double realSumZ=0, realSumX=0, realSumY=0, reciprocalSumZ=0, reciprocalSumX=0, reciprocalSumY=0;
		double actual_height = Lz* Constants.c;
		double actual_length = Lx*Constants.a;
		double actual_width = Ly*Constants.b;

		double r, rx, ry, rz, B, C;
		int x,y,z;
		
		
		for (x=-realN_cutoff;x<=realN_cutoff;x++){
			for (y=-realN_cutoff;y<=realN_cutoff;y++){
				for (z=-realN_cutoff;z<=realN_cutoff;z++){
					if (x!=0 || y!=0 || z!=0 || i!=j){	// exclude self interaction
						rz=(j.getZ(Lz,Lx)+z*actual_height)-i.getZ(Lz,Lx);
						ry=(j.getY(Lz,Lx)+y*actual_length)-i.getY(Lz,Lx);
						rx=(j.getX(Lz,Lx)+x*actual_length)-i.getX(Lz,Lx);
						r = Math.sqrt(rx*rx + ry*ry + rz*rz);
						
						B = (Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*Math.exp(-alpha*alpha*r*r))/(r*r*r);
						C = (3*Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*(3+2*alpha*alpha*r*r)*Math.exp(-alpha*alpha*r*r))/(Math.pow(r, 5));
						
						realSumZ += B - rz*rz*C;
						realSumX -= rx*rz*C;
						realSumY -= ry*rz*C;
					}
				}
			}
		}
		
		// distances in the lattice cell
		rx=i.getX(Lz,Lx)-j.getX(Lz,Lx);
		ry=i.getY(Lz,Lx)-j.getY(Lz,Lx);
		rz=i.getZ(Lz,Lx)-j.getZ(Lz,Lx);
		
		for (int mx=-reciprocalN_cutoff;mx<=reciprocalN_cutoff;mx++){
			for (int my=-reciprocalN_cutoff;my<=reciprocalN_cutoff;my++) {
				for (int mz = 0; mz <= reciprocalN_cutoff; mz++) {
					// this condition makes sure we only count half the (reciprocal) space which is what turns the exp to cos
					// also implements the condition G!=0
					if (mz > 0 || (mz==0 && mx>0) || (mz==0 && mx==0 && my>0)) {
						double Gx = 2 * Math.PI * mx / actual_length;
						double Gy = 2 * Math.PI * my / actual_length;
						double Gz = 2 * Math.PI * mz / actual_height;
						double G2 = Gx * Gx + Gy * Gy + Gz * Gz;

						reciprocalSumX += ((Gx * Gz) / G2) * Math.exp(-G2 / (4 * alpha * alpha)) * 2 * Math.cos(Gx * rx + Gy * ry + Gz * rz);
						reciprocalSumY += ((Gy * Gz) / G2) * Math.exp(-G2 / (4 * alpha * alpha)) * 2 * Math.cos(Gx * rx + Gy * ry + Gz * rz);
						reciprocalSumZ += ((Gz * Gz) / G2) * Math.exp(-G2 / (4 * alpha * alpha)) * 2 * Math.cos(Gx * rx + Gy * ry + Gz * rz);
					}
				}
			}
		}

		double[] ret = new double[3];

		// self interaction correction term to be omitted from the zz interaction for i==j:
		double selfInteraction = i.getN()==j.getN() ? 4*Math.pow(alpha,3)/(3*Math.sqrt(Math.PI)) : 0;

		ret[0] = realSumX + (4*Math.PI/(actual_height*actual_length*actual_length))*reciprocalSumX;	//zx interaction
		ret[1] = realSumY + (4*Math.PI/(actual_height*actual_length*actual_length))*reciprocalSumY;	//zy interaction
		ret[2] = realSumZ + (4*Math.PI/(actual_height*actual_length*actual_length))*reciprocalSumZ - selfInteraction;	//zz interaction
		return ret;
	
	}


	/**
	 * Calculates the dipolar interaction between two spins naively implementing periodic boundary conditions
	 * @param i spin object {@code i}
	 * @param j spin object {@code j}
	 * @param Lz number of unit cells in the z direction
	 * @param Ly number of unit cells in the y direction
	 * @param Lx number of unit cells in the x direction
	 * @return array containing 3 components of the dipolar interaction {xz,yz,zz}
	 */
	public static double[] dipolarInteraction(singleSpin i, singleSpin j, int Lz, int Ly, int Lx) {

		// get displacement with PBC
		Vector3D displacement = i.getLocation(Lz, Ly).subtract(j.getLocation(Lz, Ly));
		if (displacement.dotProduct(Constants.primitiveLatticeVectors[0].normalize()) > 0.5*Lx*Constants.primitiveLatticeVectors[0].getNorm()){
			displacement = displacement.subtract(Constants.primitiveLatticeVectors[0]);
		}
		if (displacement.dotProduct(Constants.primitiveLatticeVectors[0].normalize()) <= -0.5*Lx*Constants.primitiveLatticeVectors[0].getNorm()){
			displacement = displacement.add(Constants.primitiveLatticeVectors[0]);
		}
		if (displacement.dotProduct(Constants.primitiveLatticeVectors[1].normalize()) > 0.5*Ly*Constants.primitiveLatticeVectors[1].getNorm()){
			displacement = displacement.subtract(Constants.primitiveLatticeVectors[1]);
		}
		if (displacement.dotProduct(Constants.primitiveLatticeVectors[1].normalize()) <= -0.5*Ly*Constants.primitiveLatticeVectors[1].getNorm()){
			displacement = displacement.add(Constants.primitiveLatticeVectors[1]);
		}
		if (displacement.dotProduct(Constants.primitiveLatticeVectors[2].normalize()) > 0.5*Lz*Constants.primitiveLatticeVectors[2].getNorm()){
			displacement = displacement.subtract(Constants.primitiveLatticeVectors[2]);
		}
		if (displacement.dotProduct(Constants.primitiveLatticeVectors[2].normalize()) <= -0.5*Lz*Constants.primitiveLatticeVectors[2].getNorm()){
			displacement = displacement.add(Constants.primitiveLatticeVectors[2]);
		}

		// calculate dipolar interaction with PBC
		if (i==j)	return new double[]{0,0,0};
		double r=displacement.getNorm();
		double rz=displacement.getZ();
		double rx=displacement.getX();
		double ry=displacement.getY();
		return new double[]{(r*r-3*rx*rz)/Math.pow(r, 5),
							(r*r-3*ry*rz)/Math.pow(r, 5),
							(r*r-3*rz*rz)/Math.pow(r, 5)};
	}

	/**
	 * Calculates the dipolar interaction between two spins using periodic boundary conditions
	 * and summing periodic copies of the lattice in cubes of size {@code N_cutoff}.
	 * This is not the correct limiting shape without a demagnetization term.
	 * @param i spin object {@code i}
	 * @param j spin object {@code j}
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y direction
	 * @param N_cutoff limiting size of the cube that includes copies of the lattice that are summed
	 * @return dipolar interaction between two spins with PBC
	 */
	public static double realCalcSum3D(singleSpin i, singleSpin j, int Lz, int Lx, int N_cutoff){
		double sum = 0;
		double actual_height = Lz*Constants.c;
		double actual_length = Lx*Constants.a;
		for (int l=-N_cutoff;l<=N_cutoff;l++){
			//System.out.println(l);
			for (int m=-N_cutoff;m<=N_cutoff;m++){
				//System.out.println("0");
				for (int n=-N_cutoff;n<=N_cutoff;n++){
					if (l!=0 || m!=0 || n!=0 || i!=j){
						double rz=(j.getZ(Lz,Lx)+n*actual_height)-i.getZ(Lz,Lx);
						double r = Math.sqrt((j.getX(Lz,Lx)+l*actual_length-i.getX(Lz,Lx))*(j.getX(Lz,Lx)+l*actual_length-i.getX(Lz,Lx)) + (j.getY(Lz,Lx)+m*actual_length-i.getY(Lz,Lx))*(j.getY(Lz,Lx)+m*actual_length-i.getY(Lz,Lx)) + rz*rz);
						
						sum += (r*r - 3*rz*rz)/Math.pow(r, 5);
						/*
						if (n!=0){
							rz=(j.getZ()-n*Lz)-i.getZ();
							r = Math.sqrt((j.getX()-i.getX())*(j.getX()-i.getX()) + (j.getY()-i.getY())*(j.getY()-i.getY()) + rz*rz);
							
							sum += (r*r - 3*rz*rz)/Math.pow(r, 5);
						}
						*/
					}
				}
			}
		}
		return sum;
	}

	/**
	 * Calculates the dipolar interaction between two spins using periodic boundary conditions
	 * and summing periodic copies of the lattice in ellipsoids of a given ratio and a given length
	 * of its c axis.
	 * The limiting behavior of this calculation (with small ratios should mimic the Ewald method
	 * without a demagnetization term).
	 * @param i spin object {@code i}
	 * @param j spin object {@code j}
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x direction
	 * @param Ly number of unit cells in the y direction
	 * @param c_max length of the c-axis of the ellipsoid (the long axis, parallel to the magnetic easy axis of the system)
	 * @param ellipsoidRatio the ratio between the ellipsoid's c and a and b axes.
	 * @return
	 */
	public static double[][] realCalcSum3DE(singleSpin i, singleSpin j, int Lz, int Ly, int Lx, double c_max, double ellipsoidRatio){
		double sumX = 0, sumY = 0, sumZ = 0;
		Vector3D rij = i.getLocation(Lz,Ly).subtract(j.getLocation(Lz,Ly));

		Vector3D latticeLz = Constants.primitiveLatticeVectors[2].scalarMultiply(Lz);
		Vector3D latticeLy = Constants.primitiveLatticeVectors[1].scalarMultiply(Ly);
		Vector3D latticeLx = Constants.primitiveLatticeVectors[0].scalarMultiply(Lx);
		// n ~ z ; m ~ y ; l ~ x
		double a_max = c_max / ellipsoidRatio;
		double b_max = c_max / ellipsoidRatio;

		double[][] ellipsoidRadii = new double[10][3];
		ellipsoidRadiiFill(ellipsoidRadii,c_max,ellipsoidRatio);

		double[][] res = new double[10][3];	// maximum ellipsoid axis is divided up to 10 parts
		zeroArr(res);

		// there is no good way to know which cells will be inside the ellipsoid so we take double the length
		for (int n=-2*((int)(c_max/Lz)-1);n<=2*((int)(c_max/Lz)+1);n++){
			for (int l=-2*((int)(c_max/Lx)-1);l<=2*((int)(c_max/Lx)+1);l++){
				for (int m=-2*((int)(c_max/Ly)-1);m<=2*((int)(c_max/Ly)+1);m++){
					if (i!=j) {

						Vector3D rijCopy = new Vector3D(1, rij);	// distance between i and the copy of j that is located
						// nx unit cells along primitive direction 0
						// ny unit cells along primitive direction 1 and
						// nz unit cells along primitive direction 2

						rijCopy=rijCopy.add(n,latticeLz);
						rijCopy=rijCopy.add(m,latticeLy);
						rijCopy=rijCopy.add(l,latticeLx);
						double rz=rijCopy.getZ();
						double ry=rijCopy.getY();
						double rx=rijCopy.getX();
						double r = rijCopy.getNorm();

						// ellipsoid:

						for (int radiusIndex = ellipsoidRadii.length-1; radiusIndex >= 0; radiusIndex--) {
							if ((rx * rx) / (ellipsoidRadii[radiusIndex][0] * ellipsoidRadii[radiusIndex][0]) + (ry * ry) / (ellipsoidRadii[radiusIndex][1] * ellipsoidRadii[radiusIndex][1]) + (rz * rz) / (ellipsoidRadii[radiusIndex][2] * ellipsoidRadii[radiusIndex][2]) <= 1) {
								res[radiusIndex][2] += (r * r - 3 * rz * rz) / Math.pow(r, 5);
								res[radiusIndex][0] += (-3 * rz * rx) / Math.pow(r, 5);
								res[radiusIndex][1] += (-3 * rz * ry) / Math.pow(r, 5);
							}

						}
					}
				}
			}
		}


		return res;
	}

	/**
	 * Fills the given array with axes of 10 ellipsoids smaller than the one defined by {@code c_max}
	 * @param arr array of ellipsoid radii to fill
	 * @param c_max maximal c axis radius
	 * @param ellipsoidRatio ratio between the c axis and the a and b axes.
	 */
	private static void ellipsoidRadiiFill(double[][] arr, double c_max, double ellipsoidRatio){
		double step = c_max / 10;
		for (int i=arr.length-1;i>=0;i--){
			arr[i][2] = c_max;					//c
			arr[i][1] = c_max / ellipsoidRatio;	//b
			arr[i][0] = c_max / ellipsoidRatio;	//a
			c_max = c_max - step;
		}
	}

	/**
	 * Fills a {@code double} 2D array with zeros.
	 * @param arr array to fill w/ zeros
	 */
	private static void zeroArr(double[][] arr){
		for (double[] a : arr){
			for (double b : a){
				b=0;
			}
		}
	}



}
