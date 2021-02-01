package utilities;

import org.apache.commons.math3.special.Erf;
import simulation.montecarlo.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Properties;

import static simulation.montecarlo.Main.makeDir;

public class ewaldSum {

	// convergence tests for initial verification of ewald summation in 3D. this compares ewald results with direct calculation
	public static void convergenceTests(singleSpin[] arr, int Lz, int Lx, double a, double alpha, int N, int real_cutoff, int k_cutoff, int direct_cutoff, double ellipsoidRatio){
		// get 2 close spins
		int i=(int)(arr.length*0.5);
		int j=i+1;

		double ewald_result = calcSum3D(arr[i],arr[j],Lz,Lx,alpha,real_cutoff,k_cutoff)[0];
		System.out.println("ratio\tewald\tdirect(1)\tdirect(2)\tdirect(3)\tdirect(4)\tdirect(5)\tdirect(6)\tdirect(7)\tdirect(8)\tdirect(9)\tdirect(10)");

		double[][] direct_result;
		for (double eR=1;eR<ellipsoidRatio;eR=eR+(ellipsoidRatio/1000)) {
			direct_result = realCalcSum3DE(arr[i], arr[j], Lz, Lx, Lz*direct_cutoff, eR);
			System.out.print(eR + "\t" + ewald_result);
			for (int index=0;index<direct_result.length;index++){
				System.out.print("\t" + direct_result[index][0]);
			}
			System.out.println();
		}

/*
		int real_cutoff=4;
		int k_cutoff=2;
		for (alpha=1.0/Lz;alpha<5/Lz;alpha+=0.1/Lz)	{
			double[][][] intTable = new double[3][arr.length][arr.length]; // create interaction table that holds all the dipolar interactions. will be full even though it's symmetric. 1st array is x,y,z terms
			for (int i=0;i<arr.length;i++){
				for (int j=i;j<arr.length;j++){
					if (i!=j)
						intTable[i][j] = D*ewaldSum.calcSum3D(arr[i], arr[j], c*Lz, a*Lx, alpha, real_cutoff, k_cutoff);
				}
			}
			MonteCarloMetropolis.updateAllLocalFields(arr, intTable);

		}
*/

	}

	// more advanced convergence tests. this is meant to find optimal values for alpha, real_cutoff, k_cutoff
	public static void convergenceTests2(int Lz, int Lx, int real_cutoff, int max_k_cutoff){
		// get 2 close spins
		int i=(int)(Lx/2);
		int j=(int)(Lx/2);
		int k=(int)(Lz/2);
		int testSpin=(i)*Lx*Lz*4+(j)*Lz*4+(k)*4+0;
		int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;

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

		final double[][][] intTable = new double[3][4*Lx*Lx*Lz][4*Lx*Lx*Lz];
		final double[][] exchangeIntTable = new double[4*Lx*Lx*Lz][4*Lx*Lx*Lz];	// all zeros


		//System.out.println("alpha\txz\tyz\tzz");
		System.out.println("alpha\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12");
//		System.out.println(i +" " + j);
		for (alpha = 0.1 / Lz; alpha < 5.0 / Lz; alpha += 0.1 / Lz) {
			System.out.print(alpha*Lz);
			//for(real_cutoff=1;real_cutoff<=12;real_cutoff++) {
			for(int k_cutoff=1;k_cutoff<=max_k_cutoff;k_cutoff++) {
				Lattice lattice = new Lattice(Lx, Lz, 0.0, false, 5.51,intTable, exchangeIntTable, null, null, null, null);
				fillIntTable(lattice.getArray(), Lz, Lx, alpha, real_cutoff, k_cutoff, intTable);
				lattice.checkerBoard();
				lattice.updateAllLocalFields();
				System.out.print("\t" + (calcEnergy(lattice)/(Lx*Lx*Lz*4)));
				//System.out.println((alpha * Lz) + "\t" + ewald_result[0] + "\t" + ewald_result[1] + "\t" + ewald_result[2]);
			}
			System.out.println();
		}

	}

	public static void convergenceTests3(int Lz, int Lx, int real_cutoff, int k_cutoff){
		// get 2 close spins
		int i=(int)(Lx/2);
		int j=(int)(Lx/2);
		int k=(int)(Lz/2);
		int testSpin=(i)*Lx*Lz*4+(j)*Lz*4+(k)*4+0;
		int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;

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

		final double[][][] intTable = new double[3][4*Lx*Lx*Lz][4*Lx*Lx*Lz];
		final double[][] exchangeIntTable = new double[4*Lx*Lx*Lz][4*Lx*Lx*Lz];	// all zeros

		Lattice lattice = new Lattice(Lx, Lz, 0.0, false, 5.51, intTable, exchangeIntTable, null, null, null, null);
		singleSpin[] arr = lattice.getArray();
		System.out.println("#real_cutoff="+real_cutoff);
		System.out.println("#reciprocal_cutoff="+k_cutoff);

		DecimalFormat df = new DecimalFormat("0.0000000000");
		System.out.println("r0 - r0 : " + Arrays.toString(arr[testSpin].distance(arr[testSpin], Lz, Lx)));
		System.out.println("r0 - r1 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor1], Lz, Lx)));
		System.out.println("r0 - r2 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor2], Lz, Lx)));
		System.out.println("r0 - r3 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor3], Lz, Lx)));
		System.out.println("r0 - r4 : " + Arrays.toString(arr[testSpin].distance(arr[neighbor4], Lz, Lx)));

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

	public static void multipleSizeTestsWithvariableKcutoff(int[] Lz, int real_cutoff, int[] k_cutoffs){

		for (int i=0;i<Lz.length;i++) {
			try {
				makeDir("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\thesis\\ewald_convergence\\", "L=" + Lz[i]);
				for (int j=0;j<k_cutoffs.length;j++) {
					PrintStream fileOut = new PrintStream("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\thesis\\ewald_convergence\\L=" + Lz[i] + "\\" +k_cutoffs[j]+"_"+real_cutoff+".txt");
					System.setOut(fileOut);
					convergenceTests3(Lz[i],Lz[i],real_cutoff,k_cutoffs[j]);
				}
			} catch (FileNotFoundException e){
				e.printStackTrace();
			}
		}
	}

	public static void multipleSizeTestsWithvariableRealcutoff(int[] Lz, int[] real_cutoffs, int k_cutoff){

		for (int i=0;i<Lz.length;i++) {
			try {
				makeDir("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\thesis\\ewald_convergence\\", "L=" + Lz[i]);
				for (int j=0;j<real_cutoffs.length;j++) {
					PrintStream fileOut = new PrintStream("C:\\Users\\Tomer\\OneDrive - post.bgu.ac.il\\thesis\\ewald_convergence\\L=" + Lz[i] + "\\" +k_cutoff+"_"+real_cutoffs[j]+".txt");
					System.setOut(fileOut);
					convergenceTests3(Lz[i],Lz[i],real_cutoffs[j],k_cutoff);
				}
			} catch (FileNotFoundException e){
				e.printStackTrace();
			}
		}
	}

	public static double calcEnergy(Lattice lattice){
		double energy=0;
		singleSpin[] arr = lattice.getArray();
		for(int i=0;i<arr.length;i++){
			energy += CrystalField.getEnergy(arr[i].getLocalBx(), arr[i].getLocalBy(), arr[i].getLocalBz(), arr[i].getSpin());
		}
		return energy;
	}


	public static void main(String[] args){
		// read parameters from file:
		Properties params = GetParamValues.getParams();
		final double a=GetParamValues.getDoubleParam(params, "a");
		final double c=GetParamValues.getDoubleParam(params, "c");
		int real_cutoff=GetParamValues.getIntParam(params, "real_cutoff");
		int k_cutoff=GetParamValues.getIntParam(params, "k_cutoff");
		double alpha = GetParamValues.getDoubleParam(params, "alpha");
		params=null;
		
		int Lx=0;	// lattice x-y size
		int Lz=0;	// lattice z size
		k_cutoff=0;
		real_cutoff=0;
        // for convergence tests
        //int direct_cutoff = 1;
		//double ellipsoidRatio = 10;

        // get lattice parameters as command line arguments
        try {
        	Lx = Integer.parseInt(args[0]);
        	Lz = Integer.parseInt(args[1]);
			k_cutoff = Integer.parseInt(args[2]);
			real_cutoff = Integer.parseInt(args[3]);
        	// for convergence tests
        	//real_cutoff = Integer.parseInt(args[2]);
			//k_cutoff = Integer.parseInt(args[3]);
			//direct_cutoff = Integer.parseInt(args[4]);
			//ellipsoidRatio = Double.parseDouble(args[5]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.err.println("ArrayIndexOutOfBoundsException caught " + e.toString());
        }
        catch (NumberFormatException e){
			System.err.println("Argument not a valid number " + e.toString());
		}
		//multipleSizeTestsWithvariableKcutoff(new int[]{3,4,5,6,7,9,10},2,new int[]{2,3,4,5,6,7,8,9,10,11,12});
		//multipleSizeTestsWithvariableRealcutoff(new int[]{3,4,5,6,7,9,10},new int[]{2,3,4,5,6,7,8,9,10,11,12},6);
		//convergenceTests3(Lz, Lx, real_cutoff, k_cutoff);

        //System.exit(0);

        try (BufferedWriter out = new BufferedWriter(new FileWriter("data" + File.separator + "interactions" + File.separator + "intTable_"+Lx+"_"+Lz+".txt"))){
            
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
            
            int N = Lx*Lx*Lz*4;

            singleSpin[] arr = GenerateLattice.generate_ising_lattice(Lx, Lz, a, c, 1, 0, null);
            /*
            for (int i=0;i<arr.length;i++){
            	if (i%2==0){
            		arr[i].setSpin(-1);
            	}else{
            		arr[i].setSpin(1);
            	}
            }
            */

            // fill interactions table and print it to the file 'out'
            fillIntTable(arr,Lz,Lx,alpha/(c*Lz),real_cutoff,k_cutoff,out);

        }
        catch (IOException e) { System.out.println("bad file"); }

        // for convergence tests
		//int N = Lx*Lx*Lz*4;
		//singleSpin[] arr = GenerateLattice.generate_ising_lattice(Lx, Lz, a, c, 1, 0, null);
		//convergenceTests(arr,Lz,Lx,a,2/(c*Lz),N,real_cutoff,k_cutoff,direct_cutoff, ellipsoidRatio);
		//convergenceTests2(arr,Lz,Lx,real_cutoff,k_cutoff);





    }

	public static void fillIntTable(singleSpin[] arr, int Lz, int Lx, double alpha, int real_cutoff, int k_cutoff, final double[][][] intTable){
		final double c = -Constants.mu_0*Constants.mu_B*Constants.g_L*0.25/Math.PI;	// coefficient dipolar spin-spin interaction. The minus sign is because
		for (int i=0;i<arr.length;i++){
			for (int j=i;j<arr.length;j++){
				double[] interaction = new double[3];
				interaction = ewaldSum.calcSum3D(arr[i], arr[j], Lz, Lx, alpha, real_cutoff, k_cutoff);

				//System.out.println(Double.toString(-D*interaction));
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
	
	public static void fillIntTable(singleSpin[] arr, int Lz, int Lx, double alpha, int real_cutoff, int k_cutoff, BufferedWriter out) throws IOException{
		DecimalFormat df = new DecimalFormat("0.0000000000");
		for (int i=0;i<arr.length;i++){
			for (int j=i;j<arr.length;j++){
				double[] interaction = new double[3];
				interaction = ewaldSum.calcSum3D(arr[i], arr[j], Lz, Lx, alpha, real_cutoff, k_cutoff);
				
				//System.out.println(Double.toString(-D*interaction));
				System.out.println("("+i+","+j+")" + " : " + interaction[0] + "," + interaction[1] + "," + interaction[2]);
				out.write(df.format(interaction[0])+","+df.format(interaction[1])+","+df.format(interaction[2]));	// print 10 significant digits
				out.newLine();
			}	
		}
	}
	

	// bottom line: this is the right implementation of ewald sum for this project
	// returns the zx, zy, zz interations in a size 3 array: {zx, zy, zz}
	public static double[] calcSum3D(singleSpin i, singleSpin j, int Lz, int Lx, double alpha, int realN_cutoff, int reciprocalN_cutoff){
		
		double realSumZ=0, realSumX=0, realSumY=0, reciprocalSumZ=0, reciprocalSumX=0, reciprocalSumY=0;
		double actual_height = Lz* Constants.c;
		double actual_length = Lx*Constants.a;

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


	public static double[][] realCalcSum3DE(singleSpin i, singleSpin j, int Lz, int Lx, double c_max, double ellipsoidRatio){
		double sumX = 0, sumY = 0, sumZ = 0;
		double actual_height = Lz*Constants.c;
		double actual_length = Lx*Constants.a;
		// n ~ z ; m ~ y ; l ~ x
		double a_max = c_max / ellipsoidRatio;
		double b_max = c_max / ellipsoidRatio;

		double[][] ellipsoidRadii = new double[10][3];
		ellipsoidRadiiFill(ellipsoidRadii,c_max,ellipsoidRatio);

		double[][] res = new double[10][3];	// maximum ellipsoid axis is divided up to 10 parts
		zeroArr(res);
		for (int n=-(int)(c_max/actual_height)-1;n<=(int)(c_max/actual_height)+1;n++){
			for (int l=-(int)(a_max/actual_length)-1;l<=(int)(a_max/actual_length)+1;l++){
				for (int m=-(int)(b_max/actual_length)-1;m<=(int)(b_max/actual_length)+1;m++){
					if (i!=j) {
						double rz = (j.getZ(Lz,Lx) + n * actual_height) - i.getZ(Lz,Lx);
						double rx = (j.getX(Lz,Lx) + l * actual_length) - i.getX(Lz,Lx);
						double ry = (j.getY(Lz,Lx) + m * actual_length) - i.getY(Lz,Lx);

						// ellipsoid:

						for (int radiusIndex = ellipsoidRadii.length-1; radiusIndex >= 0; radiusIndex--) {
							if ((rx * rx) / (ellipsoidRadii[radiusIndex][0] * ellipsoidRadii[radiusIndex][0]) + (ry * ry) / (ellipsoidRadii[radiusIndex][1] * ellipsoidRadii[radiusIndex][1]) + (rz * rz) / (ellipsoidRadii[radiusIndex][2] * ellipsoidRadii[radiusIndex][2]) <= 1) {
								double r = Math.sqrt(rx * rx + ry * ry + rz * rz);

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

	private static void ellipsoidRadiiFill(double[][] arr, double c_max, double ellipsoidRatio){
		double step = c_max / 10;
		for (int i=arr.length-1;i>=0;i--){
			arr[i][2] = c_max;					//c
			arr[i][1] = c_max / ellipsoidRatio;	//b
			arr[i][0] = c_max / ellipsoidRatio;	//a
			c_max = c_max - step;
		}
	}

	private static void zeroArr(double[][] arr){
		for (double[] a : arr){
			for (double b : a){
				b=0;
			}
		}
	}



}
