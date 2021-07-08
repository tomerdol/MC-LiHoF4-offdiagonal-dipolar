package simulation.montecarlo;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.Arrays;
import java.util.Properties;

/**
 * Constants to be used throughout the program. They are read from a file (parameters_{system}.properties) by the class {@link GetParamValues}
 */
public class Constants {
    /** sizes of unit cell a,b,c */
    public static final double a, b, c;
    /** Vacuum permeability */
    public static final double mu_0;
    /** Bohr magneton */
    public static final double mu_B;
    /** Lande g-factor of the ions */
    public static final double g_L;
    /** Boltzmann constant */
    public static final double k_B;
    /** Number of spins in the basis */
    public static final int num_in_cell;
    /** An array containing the three primitive lattice vectors of the Bravais lattice */
    public static final Vector3D[] primitiveLatticeVectors;
    /** An array containing the fractional positions of spins in the unit cell (for a lattice w/ a basis) */
    public static final double[][] basis;
    /** The name of the system; e.g. LiHoF4 or Fe8 */
    public static final String systemName;

    static{
        // This is received as a -D variable from the command line
        systemName=System.getProperty("system");

        // read a,c from parameters file:
        Properties params = GetParamValues.getParams();
        // lattice unit cell sizes:
        primitiveLatticeVectors = new Vector3D[3];    // 3 dimensions
        primitiveLatticeVectors[0] = new Vector3D(GetParamValues.getLatticeVector(params, "a"));
        primitiveLatticeVectors[1] = new Vector3D(GetParamValues.getLatticeVector(params, "b"));
        primitiveLatticeVectors[2] = new Vector3D(GetParamValues.getLatticeVector(params, "c"));
        a=primitiveLatticeVectors[0].getNorm();    // a is x
        b=primitiveLatticeVectors[1].getNorm();    // b is y
        c=primitiveLatticeVectors[2].getNorm();    // c is z

        mu_0=GetParamValues.getDoubleParam(params, "mu_0");	// Vacuum permeability
        mu_B=GetParamValues.getDoubleParam(params, "mu_B");	// Bohr Magneton
        k_B=GetParamValues.getDoubleParam(params,"k_B");    // Boltzmann constant
        g_L=GetParamValues.getDoubleParam(params, "g_L");	    // g-factor

        num_in_cell = GetParamValues.getIntParam(params, "num_in_cell");    // number of spins in unit cell
        basis = new double[num_in_cell][3];	// 3D coordinate location for each of the atoms in the basis

        // fill location:
        for (int l=0;l<basis.length;l++){
            basis[l]=GetParamValues.getLocation(params, l);
        }

        params=null;
    }

    /**
     * Constructs a string containing all physical constant and their values
     * @return a string of physical constants separated by "|"
     */
    public static String constantsToString(){
        return String.format("a=%6.3e|c=%6.3e|k_B=%8.5e|mu_0=%8.5e|mu_B=%8.5e|g_L=%f|num_in_cell=%d|",a,c,k_B,mu_0,mu_B,g_L,num_in_cell);
    }

    /**
     * Constructs a string containing locations of spins in the basis
     * @return the basis of the crystal, separated by "|"
     */
    public static String locationsToString(){
        String ret="locations:";
        for (int i=0;i<basis.length;i++){
            ret = ret + Arrays.toString(basis[i]) + "|";
        }
        return ret;
    }

    /**
     * Constructs a string containing the primitive lattice vectors
     * @return the primitive lattice vectors, separated by "|"
     */
    public static String latticeVectorsToString(){
        String ret="Lattice Vectors:";
        for (int i=0;i<primitiveLatticeVectors.length;i++){
            ret = ret + primitiveLatticeVectors[i].toString() + "|";
        }
        return ret;
    }

}
