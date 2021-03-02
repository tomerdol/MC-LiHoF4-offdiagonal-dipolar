package utilities;

import org.apache.commons.math3.random.MersenneTwister;
import simulation.montecarlo.Constants;
import simulation.montecarlo.CrystalField;
import simulation.montecarlo.GetParamValues;
import simulation.montecarlo.singleSpin;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;


public class GenerateLattice {

	/**
	 * Generates a singleSpin array of LiHo_{x}Y_{1-x}F_4
	 * @param Lx - Number of unit cells in x and in y directions
	 * @param Lz - Number of unit cells in z direction
	 * @param dilution - Dilution of Holmium atoms (x)
	 * @param h - Random field parameter
	 * @param rnd_spin - PRNG for random dilution. If null is given then an empty lattice will be returned (all spin 0) which still can be useful for Ewald summation.
	 * @return An array of spins that represent a Lx*Lx*Lz lattice of LiHo_{x}Y_{1-x}F_4
	 */
	public static singleSpin[] generate_ising_lattice(int Lx, int Lz, double dilution, double h, MersenneTwister rnd_spin) {
		return generate_ising_lattice(Lx, Lx, Lz, dilution, h, rnd_spin);
	}

	/**
	 * Generates a singleSpin array of LiHo_{x}Y_{1-x}F_4
	 * @param Lx - Number of unit cells in x direction
	 * @param Ly - Number of unit cells in y direction
	 * @param Lz - Number of unit cells in z direction
	 * @param dilution - Dilution of Holmium atoms (x)
	 * @param h - Random field parameter
	 * @param rnd_spin - PRNG for random dilution. If null is given then an empty lattice will be returned (all spin 0) which still can be useful for Ewald summation.
	 * @return An array of spins that represent a Lx*Lx*Lz lattice of LiHo_{x}Y_{1-x}F_4
	 */
	public static singleSpin[] generate_ising_lattice(int Lx, int Ly, int Lz, double dilution, double h, MersenneTwister rnd_spin) {
		final double spinSize = CrystalField.getMagneticMoment(0.0, 0.0, 0.05);

        int i, j, k, l;
        // create the array that will hold the lattice. the array's cells correspond
        // to the unit cells of the LiHo{x}Y{x-1}F4.
        singleSpin[] arr = new singleSpin[Constants.num_in_cell*Lx*Ly*Lz];
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // get location matrix (3D coordinates for each of the 4 atoms)
        Properties params = GetParamValues.getParams();
        double[][] location = new double[Constants.num_in_cell][3];	// 3D coordinate location for each of the atoms in the basis
        
        // fill location:
        for (l=0;l<location.length;l++){
        	location[l]=GetParamValues.getLocation(params, l);
        }
        params=null;
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        for (i = 0; i < Lx; i++)
        {
            for (j = 0; j < Ly; j++)
            {
            	for (k = 0; k < Lz; k++)
            	{
            		// the spins in each unit cell are designated 0-3. see documentation for further info
            		// a lattice is created with all spins +1
            		for (l = 0; l < Constants.num_in_cell; l++)
            		{
	            		if (rnd_spin!=null && rnd_spin.nextDouble()<dilution){
	            			int s=1;
	            			//if (rnd_spin.nextBoolean()) s=-1;
	            			arr[i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+l]=new singleSpin(1,i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+l, spinSize);
	            		}else{
		                	arr[i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+l]=new singleSpin(0,i*Ly*Lz*Constants.num_in_cell+j*Lz*Constants.num_in_cell+k*Constants.num_in_cell+l, spinSize);
	            		}
            		}
            	}
            }
        }
        
        return arr;
    }

	public static void printArr(singleSpin[] arr, int Lz, int Lx, BufferedWriter out){
		for (int i=0; i<arr.length; i++){
			try{

//				out.write(arr[i].getX(Lz,Lx)+","+arr[i].getY(Lz,Lx)+","+arr[i].getZ(Lz,Lx)+","+arr[i].getSpin()+","+arr[i].getH()+","+arr[i].getN());
				out.newLine();

			}
			catch (IOException e) {}
		}
	}
	
	public static void main(String[] args) {
		// read a,c from parameters file:
		Properties params = GetParamValues.getParams();
		final double a=GetParamValues.getDoubleParam(params, "a");
		final double c=GetParamValues.getDoubleParam(params, "c");
		params=null;
		
		int Lx;	// lattice x-y size
		int Lz;	// lattice z size
        double dilution;	// starting dilution
        double h;
        int numOfConfigurations = 2;
        Lx=6;
        Lz=16;
        dilution=0.4;
        h=0;
        int startConfigNum = 1;
        
        // get lattice parameters as command line arguments
        try {
        	Lx = Integer.parseInt(args[0]);
        	Lz = Integer.parseInt(args[1]);
            dilution = Double.parseDouble(args[2]);
            h = Double.parseDouble(args[3]);
            startConfigNum = Integer.parseInt(args[4]);
            numOfConfigurations = Integer.parseInt(args[5]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("ArrayIndexOutOfBoundsException caught");
        }
        catch (NumberFormatException e){}
        
        for (int i=startConfigNum;i<=numOfConfigurations;i++){
        	BufferedWriter out = null;
        	try {
	            out = new BufferedWriter(new FileWriter("data" + File.separator + "configurations" + File.separator + "config_"+Lx+"_"+Lz+"_"+dilution+"_"+h+"_"+i+".txt"));
	            
	            //initialize PRNG
	            long seed = System.currentTimeMillis();
	            MersenneTwister rnd = new MersenneTwister(seed);
	            
	            //print lattice sizes
	            out.write("Lx="+Lx);
	            out.newLine();
	            out.write("Lz="+Lz);
	            out.newLine();
	            out.write("seed="+seed);
	            out.newLine();
	            
	            // print the lattice itself
	            printArr(generate_ising_lattice(Lx,Lz,dilution,h, rnd), Lz, Lx, out);
	            
	            out.close();
	        }
	        catch (IOException e) { System.out.println("bad file"); }
	        finally { 
	        	if (out!=null) {
	        		try {out.close(); }
	        		catch (IOException e) {System.out.println("error closing file"); };
	        	}
	        }
        }

	}

}
