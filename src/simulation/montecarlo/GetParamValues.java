package simulation.montecarlo;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * Get parameters and constants from the parameter file "parameters.properties".
 * @see Constants
 */
public class GetParamValues {

	/**
	 * Creates a {@link Properties} object that enables reads from the "parameters.properties" file
	 * @return Properties object to read from the "parameters.properties" file
	 */
	public static Properties getParams() {

		Properties params = new Properties();
		InputStream input = null;
		
		try {
		
			input = new FileInputStream("parameters.properties");
		
			// load a properties file
			params.load(input);
		
		} catch (IOException ex) {
			ex.printStackTrace();
		} finally {
			if (input != null) {
				try {
					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return params;
	}

	/**
	 * Receive an integer parameter from the parameter file
	 * @param params Properties object
	 * @param param Name of the parameter to receive
	 * @return Value of the (int) parameter
	 * @throws RuntimeException when the given name does not match any of the possible int parameters
	 */
	public static int getIntParam(Properties params, String param){
		// parameters are either doubles or integers
		int p=0;
		try{
			if (param.equals("real_cutoff") || param.equals("k_cutoff") || param.equals("num_in_cell")){	// make sure param is expected to be integer
				if (params.getProperty(param)!=null)
					p = Integer.parseInt(params.getProperty(param));
				else
					throw new RuntimeException("value not found in parameter file");
			}else{
				throw new RuntimeException("error reading parameter (integer expected)");
			}

		}catch(Exception e){
			throw new RuntimeException("error reading parameter");
		}
		
		return p;
	}

	/**
	 * Receive an integer parameter from the parameter file
	 * @param params Properties object
	 * @param param Name of the parameter to receive
	 * @return Value of the (int) parameter
	 * @throws RuntimeException when the given name does not match any of the possible int parameters
	 */
	public static long getLongParam(Properties params, String param){
		// parameters are either doubles or integers
		long p=0;
		try{
			if (param.substring(0,param.length()-1).equals("obsPrintSweepNum")){	// make sure param is expected to be integer
				if (params.getProperty(param)!=null)
					p = Long.parseLong(params.getProperty(param));
				else
					throw new RuntimeException("Value not found in parameter file");
			}else{
				throw new RuntimeException("error reading parameter (expected long)");
			}

		}catch(Exception e){
			throw new RuntimeException("error reading parameter");
		}

		return p;
	}

	/**
	 * Receive the location within the unit cell of one of the ions, by number (see "ion_locations.pdf" in the documentation)
	 * @param params Properties object
	 * @param num index of the ion (in the basis) whose location should be given
	 * @return location on the ion within the unit cell
	 */
	public static double[] getLocation(Properties params, int num){
		double[] location = new double[3];	// 3 coordinates: x,y,z
		String str = params.getProperty(Integer.toString(num) + "loc");
		if (str==null)
			throw new RuntimeException("value not found in parameter file");
		try{
			for (int i=0;i<location.length;i++){
				location[i] = Double.parseDouble(str.split(",")[i]);
			}
		}catch(Exception e){
			throw new RuntimeException("error reading location parameter");
		}
		
		return location;
	}

	/**
	 * Receive a double parameter from the parameter file
	 * @param params Properties object
	 * @param param Name of the parameter to receive
	 * @return Value of the (double) parameter
	 * @throws RuntimeException when the given name does not match any of the possible double parameters
	 */
	public static double getDoubleParam(Properties params, String param){
		// parameters are either doubles or integers
		double p = 0;
		try{
			// param is expected to be double
			if (param.equals("a") || param.equals("c") || param.equals("D") || param.equals("J_ex")
					|| param.equals("precision") || param.equals("alpha") || param.equals("tol") || param.equals("mu_0")
					|| param.equals("mu_B") || param.equals("g_L") || param.equals("s") || param.equals("k_B")){
				if (params.getProperty(param)!=null)
					p = Double.parseDouble(params.getProperty(param));
				else
					throw new RuntimeException("value not found in parameter file");
			}else{
				throw new RuntimeException("error reading parameter (double expected)");
			}

		}catch(Exception e){
			throw new RuntimeException("error reading parameter");
		}
		
		return p;
	}

}
