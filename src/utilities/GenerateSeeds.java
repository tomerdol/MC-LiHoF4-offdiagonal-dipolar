package utilities;

import org.apache.commons.math3.random.*;

/**
 * Generates seed for the random number generators used in the simulation
 */
public class GenerateSeeds {

    /**
     * Generates a desired number of random seeds
     * @param args the task ID and the number of desired number of seeds
     */
    public static void main(String[] args){
        int taskID=0, numOfSeeds=0;
        // get the desired number of seeds and the task ID from the command line
        try {
            taskID = Integer.parseInt(args[0]);
            numOfSeeds = Integer.parseInt(args[1]);
        } catch (NumberFormatException e) {
            System.err.println("input should be numeric");
            e.printStackTrace();
            System.exit(1);
        }
        // also use the system's time
        long time = System.currentTimeMillis();
        // using Apache's math3 library, create a long seed
        RandomGenerator tmpRnd = new JDKRandomGenerator();
        tmpRnd.setSeed(new int[]{(int)time, taskID});

        // feed the seed to a MersenneTwister RNG object
        RandomGenerator rndMT19937 = new MersenneTwister(tmpRnd.nextLong());

        // use the MersenneTwister to generate how many seeds are needed
        long[] seeds = generateSeeds(rndMT19937, numOfSeeds);
        for (int i=0;i<numOfSeeds;i++) System.out.println(seeds[i]);

    }

    /**
     * Generates a desired number of random seeds while ensuring there are no accidental duplicates
     * @param rnd random number generator
     * @param numOfSeeds desired number of seeds to create
     * @return array of generated seeds
     */
    public static long[] generateSeeds(RandomGenerator rnd, int numOfSeeds){
        long[] seeds = new long[numOfSeeds];
        // generate the (always positive) seeds
        for (int i=0;i<numOfSeeds;i++)  seeds[i]=Math.abs(rnd.nextLong());

        // check that there are no duplicates
        boolean duplicates=false;
        for (int j=0;j<seeds.length && !duplicates;j++)
            for (int k=j+1;k<seeds.length && !duplicates;k++)
                if (k!=j && seeds[k] == seeds[j])
                    duplicates=true;
        if (duplicates){
            throw new RuntimeException("found duplicates in created seeds. this is highly unlikely and needs to be checked. If everything is ok just run again. ");
        }else{
            return seeds;
        }
    }
}
