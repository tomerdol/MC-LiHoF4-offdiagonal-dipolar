package simulation.montecarlo;

import java.util.Arrays;
import java.util.Properties;

/**
 * Constants to be used throughout the program. They are read from a file (parameters.properties) by the class {@link GetParamValues GetParamValues}
 */
public class Constants {
    public static final double a ,c, mu_0, mu_B, g_L, k_B, J_ex;
    public static double spinSize, tol; // TODO: move these to MonteCarloSimulation (and J_ex)
    public static final int num_in_cell;
    public static final double[][] location;

    static{
        // read a,c from parameters file:
        Properties params = GetParamValues.getParams();
        // lattice unit cell sizes:
        a=GetParamValues.getDoubleParam(params, "a");
        c=GetParamValues.getDoubleParam(params, "c");
        mu_0=GetParamValues.getDoubleParam(params, "mu_0");	// Vacuum permeability
        mu_B=GetParamValues.getDoubleParam(params, "mu_B");	// Bohr Magneton
        k_B=GetParamValues.getDoubleParam(params,"k_B");    // Boltzmann constant
        g_L=GetParamValues.getDoubleParam(params, "g_L");	    // g-factor
        spinSize=GetParamValues.getDoubleParam(params, "s");	// s_z size
        // minimum percentage to verify ground state:
        tol = GetParamValues.getDoubleParam(params, "tol");

        num_in_cell = GetParamValues.getIntParam(params, "num_in_cell");    // number of spins in unit cell
        J_ex = GetParamValues.getDoubleParam(params, "J_ex");
        location = new double[num_in_cell][3];	// 3D coordinate location for each of the 4 atoms

        // fill location:
        for (int l=0;l<location.length;l++){
            location[l]=GetParamValues.getLocation(params, l);
        }

        params=null;
    }

    public static String constantsToString(){
        return String.format("a=%6.3e|c=%6.3e|J_ex=%f|k_B=%8.5e|mu_0=%8.5e|mu_B=%8.5e|g_L=%f|spinSize=%.8f|tol=%4.1e|num_in_cell=%d|",a,c,J_ex,k_B,mu_0,mu_B,g_L,spinSize,tol,num_in_cell);
    }

    public static String locationsToString(){
        String ret="locations:";
        for (int i=0;i<location.length;i++){
            ret = ret + Arrays.toString(location[i]) + "|";
        }
        return ret;
    }


}
