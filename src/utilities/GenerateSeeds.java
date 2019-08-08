package utilities;

import org.apache.commons.math3.random.*;

public class GenerateSeeds {
    public static void main(String[] args){
        int taskID=0, numOfSeeds=0;
        try {
            taskID = Integer.parseInt(args[0]);
            numOfSeeds = Integer.parseInt(args[1]);
        } catch (NumberFormatException e) {
            System.err.println("input should be numeric");
            e.printStackTrace();
            System.exit(1);
        }

        long time = System.currentTimeMillis();
        RandomGenerator tmpRnd = new JDKRandomGenerator();
        tmpRnd.setSeed(new int[]{(int)time, taskID});

        RandomGenerator rndMT19937 = new MersenneTwister(tmpRnd.nextLong());

        long[] seeds = generateSeeds(rndMT19937, numOfSeeds);
        for (int i=0;i<numOfSeeds;i++) System.out.println(seeds[i]);


    }

    public static long[] generateSeeds(RandomGenerator rnd, int numOfSeeds){
        long[] seeds = new long[numOfSeeds];
        for (int i=0;i<numOfSeeds;i++)  seeds[i]=Math.abs(rnd.nextLong());

        // check no duplicates
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
