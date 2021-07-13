package simulation.montecarlo;

import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Formatter;
import java.util.InputMismatchException;
import java.util.stream.IntStream;

/**
 * Manages the output of a simulation
 */
public class OutputWriter implements Closeable {
    /** Whether to print the progress of the simulation using a terminal progress bar and percentage */
    private final boolean printProgress;
    /** Whether to print the output to the console (as opposed to saving to a file) */
    private final boolean printOutputToConsole;
    /** Number of rows to keep buffered before flushing */
    private final long numOfBufferedRows;
    /** Name of the project and the directory where the results will be saved */
    private final String folderName;
    /** The type of output - verbose, binned or per spin info */
    private OutputType outType;
    /** Array containing the column widths of the output (# of characters) */
    private final int[] colWidths;
    /** Array containing the column headers */
    private final String[] colNames;
    /** The size (in # of chars) of the output buffer */
    private final int bufferSize;
    /** {@code FileWriter} that performs the output */
    private final FileWriter out;
    /** The buffer that holds buffered results until they are outputted */
    private StringBuilder outputBuffer;

    public long getNumOfBufferedRows() {
        return numOfBufferedRows;
    }

    public boolean isPrintProgress() {
        return printProgress;
    }

    public boolean isPrintOutputToConsole() {
        return printOutputToConsole;
    }

    /**
     * Writes out the result lines held in the buffer
     */
    public void flush(){
        try {
            print(outputBuffer.toString(), false);
            if (!printOutputToConsole) out.flush();
            if (outputBuffer.length() > bufferSize) // warn that the buffer was not large enough to hold the data
                                                    // (had to be extended which is inefficient)
                System.err.println("output buffer exceeded its buffer with char count of " + outputBuffer.length() + '/' + bufferSize);
            outputBuffer = null;
            outputBuffer = new StringBuilder(bufferSize);
        }catch (IOException e){
            throw new RuntimeException("Error saving results: " + e.getMessage());
        }
    }

    /**
     * Uses the {@code FileWriter} or the Console to write out the given string.
     * @param toPrint String to print (to console or to a file)
     * @param lineBreak whether a line break should be appended at the end of the string
     * @throws IOException
     */
    public void print(String toPrint, boolean lineBreak) throws IOException {
        if (printOutputToConsole && lineBreak) {
            System.out.println(toPrint);
        }else if (printOutputToConsole){
            System.out.print(toPrint);
        }else if (!printOutputToConsole && lineBreak){
            out.write(toPrint + System.lineSeparator());
        }else{
            out.write(toPrint);
        }
    }

    public void close() throws IOException{
        this.out.close();
    }

    /**
     * Constructs table headers using the {@code colNames} array
     * @return a row of column names separated by whitespace so that they constitute table headers
     */
    public String makeTableHeader(){
        if (colNames.length!=colNames.length)
            throw new InputMismatchException("array of column names should be the same length as array of column widths.");
        int totalChrsForRow = IntStream.of(colWidths).sum()+colWidths.length+10;    //+10 just to be safe
        StringBuilder headerRow = new StringBuilder(totalChrsForRow);
        // construct the column headers
        for (int i=0;i<colNames.length;i++){
            if (colNames[i].length()>colWidths[i]){
                throw new InputMismatchException("the column header "+colNames[i]+" is longer than its designated" +
                        "column width of "+colWidths[i]);
            } else {
                // pad the column header with the required amount of whitespace so that its width,
                // which is the length of its name + some whitespace, agrees with the requirement in colWidths
                for (int pad=0;pad<=colWidths[i]-colNames[i].length();pad++){
                    headerRow.append(" ");
                }
                // write the column name
                headerRow.append(colNames[i]);
            }
        }
        if (headerRow.length()>totalChrsForRow)
            System.err.println("header row builder not allocated enough space: " + headerRow.length() + "/" + (totalChrsForRow+10));
        return headerRow.toString();
    }

    /**
     * Writes a row of results into {@code outputBuffer}
     * @param sweeps row number (sweep index)
     * @param m magnetization
     * @param currentEnergy energy
     * @param magField0 mean of local Bx
     * @param magField1 std of local Bx
     * @param magField2 mean of local By
     * @param magField3 std of local By
     * @param magField4 mean of local Bz
     * @param magField5 std of local Bz
     * @param magField6 maximal transverse magnetic field
     * @param magField7 maximal longitudinal magnetic field
     * @param spinSizes0 mean of magnetic moments
     * @param spinSizes1 std of magnetic moments
     * @param mk2x |m(k_x)|^2
     * @param mk2y |m(k_y)|^2
     * @param mk2z |m(k_z)|^2
     * @param lastSwapAccepted whether the last swap attempt between temperature (in parallel tempering) was accepted
     */
    public void writeObservablesVerbose(final long sweeps, final double m, final double currentEnergy, final double magField0, final double magField1, final double magField2,
                                 final double magField3, final double magField4, final double magField5, final double magField6, final double magField7,
                                 final double spinSizes0, final double spinSizes1, final double mk2x, final double mk2y, final double mk2z, final boolean lastSwapAccepted){
        if (outType == OutputType.VERBOSE) {
            Formatter formatter = new Formatter(outputBuffer);
            formatter.format(makeTableRowFormat(new char[]{'d', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'c'}), sweeps, m, currentEnergy, magField0, magField1, magField2, magField3, magField4, magField5, magField6, magField7, spinSizes0, spinSizes1, mk2x, mk2y, mk2z, (lastSwapAccepted ? '1' : '0'));
            outputBuffer.append(System.lineSeparator());
        }
    }

    /**
     * Writes a row of the lattice output to {@code outputBuffer}
     * @param n row number (spin index)
     * @param spin the spin orientation
     * @param spinSize the magnetic moment of the spin
     * @param localBx local magnetic field in the x direction
     * @param localBy local magnetic field in the y direction
     * @param localBz local magnetic field in the z direction
     */
    public void writeObservablesPerSpin(final int n, final int spin, final double spinSize, final double localBx, final double localBy, final double localBz){
        if (outType == OutputType.SPIN) {
            Formatter formatter = new Formatter(outputBuffer);
            formatter.format(makeTableRowFormat(new char[]{'d', 'd', 'g', 'g', 'g', 'g'}), n, spin, spinSize, localBx, localBy, localBz);
            outputBuffer.append(System.lineSeparator());
        }
    }

    /**
     * Writes a row of results into {@code outputBuffer} in the binned format
     * @param currentBinCount number of samples in the current bin
     * @param binAvg an array of observable averages over the bin
     * @param acceptanceRateForBin average acceptance rate for the bin
     */
    public void writeObservablesBin(final long currentBinCount, final double[] binAvg, final double acceptanceRateForBin){
        if (outType == OutputType.BIN){
            Formatter formatter = new Formatter(outputBuffer);
            // the first and last columns are treated separately, all the rest is just pasted in the middle by the "%s"
            formatter.format(" % "+colWidths[0]+"d%s % "+(colWidths[colWidths.length-1])+".2f%n", currentBinCount, avgArrToString(binAvg, currentBinCount, 1) , acceptanceRateForBin);
        }
    }

    /**
     * Appends the last acceptance state of the parallel tempering attempt to the end of a row
     * @param swapAccepted whether the last attempt was accepted
     */
    public void writeSwapAcceptance(final boolean swapAccepted){
        if (outType == OutputType.VERBOSE){
            outputBuffer.append(new String(new char[colWidths[colWidths.length-1]-2]).replace('\0', ' ') + (swapAccepted ? '1' : '0'));
            outputBuffer.append(System.lineSeparator());
        }
    }

    /**
     * Formats the data to a row with the required widths and numeric types
     * @param colTypes array containing the numeric types of the columns
     * @return String used by a java Formatter object to construct a row with the given widths in the given format
     * @see <a href="https://docs.oracle.com/javase/7/docs/api/java/util/Formatter.html" target="_top">Class Formatter</a>
     */
    public String makeTableRowFormat(char[] colTypes){
        if (colWidths.length!=colTypes.length)
            throw new InputMismatchException("array of column types should be the same length as array of column widths. colWidths.length=" + colWidths.length + ", colTypes.length="+colTypes.length);
        int totalChrsForRow = IntStream.of(colWidths).sum()+colWidths.length+10;    // +10 just to be safe
        StringBuilder rowFormat = new StringBuilder(totalChrsForRow);
        for (int i=0;i<colTypes.length;i++){
            // % 10d
            if (colTypes[i]=='g')
                rowFormat.append(" % "+colWidths[i]+'.'+(colWidths[i]-8)+colTypes[i]);
            else if (colTypes[i]=='d')
                rowFormat.append(" % "+colWidths[i]+colTypes[i]);
            else if (colTypes[i]=='c')
                rowFormat.append(" %"+colWidths[i]+colTypes[i]);
            else
                throw new IllegalArgumentException("column types can be only 'd' or 'g'");
        }
        if (rowFormat.length()>totalChrsForRow)
            System.err.println("string format builder not allocated enough space: "+ rowFormat.length()+'/'+(colTypes.length*9));
        return rowFormat.toString();
    }

    /**
     * Given an array of (un-averaged) bin data to write as a row in the results, automatically averages over the bin
     * and constructs a string for that row in the given format (correct column widths)
     * @param avgArr array of bin data to average and format in columns
     * @param N number of samples in bin
     * @param startIndex where to start in the colWidths array (currently the first column is treated separately, so should be 1)
     * @return a string of the averaged bin data
     * @see <a href="https://docs.oracle.com/javase/7/docs/api/java/util/Formatter.html" target="_top">Class Formatter</a>
     */
    public String avgArrToString(double[] avgArr, long N, int startIndex) {
        if (avgArr.length+startIndex>colWidths.length)
            throw new InputMismatchException("start index too high such that startIndex+avgArr.length > colWidths.length.");
        int totalChrsForRow = IntStream.of(colWidths).sum()+ colWidths.length + 10; // +10 just to be safe
        StringBuilder str = new StringBuilder(totalChrsForRow);
        Formatter formatter = new Formatter(str);
        int i;
        for (i=0;i<avgArr.length;i++){
            formatter.format(" % "+(colWidths[startIndex+i])+'.'+(colWidths[startIndex+i]-8)+"g",(avgArr[i]/N));
        }
        if (str.length()>totalChrsForRow)
            System.err.println("avg arr row builder not allocated enough space: " + str.length() + "/" + (totalChrsForRow));
        return str.toString();
    }

    public int getBufferSize(){ return this.bufferSize; }

    public String getFolderName() { return this.folderName; }

    public OutputType getOutType() { return outType;
    }

    /**
     * Creates a new {@code OutputWriter} with all configuration options that have been specified by calling methods on the given builder.
     * @param builder builder with all required and optional configuration options that have been specified by calling methods on it.
     */
    public OutputWriter(Builder builder) {
        this.printProgress = builder.printProgress;
        this.printOutputToConsole = builder.printOutputToConsole;
        this.numOfBufferedRows = builder.numOfBufferedRows;
        this.folderName = builder.folderName;
        this.outType = builder.outType;
        this.colNames = builder.colNames;
        this.colWidths = builder.colWidths;
        this.bufferSize = builder.bufferSize;
        this.out = builder.out;
        this.outputBuffer = new StringBuilder(bufferSize);
    }

    /**
     * A <i>builder</i> class for creating instances of {@code OutputWriter}
     */
    public static class Builder {
        // required
        private OutputType outType;
        private final String folderName;
        private final long numOfBufferedRows;
        private final FileWriter out;

        // optional - default values below
        private boolean printProgress=false;
        private boolean printOutputToConsole=false;
        // optional - default values to be determined later based on other input
        private int bufferSize=0;

        // created without outside input
        private final int[] colWidths;
        private final String[] colNames;

        /**
         * Creates a new <i>builder</i> object with the mandatory input fields
         * @param outType type of output (verbose, binned, per spin)
         * @param folderName name of the project which is also the name of the directory where results will be saved
         * @param numOfBufferedRows number of rows to keep in the buffer between flushes
         * @param out {@code FileWriter} object used for outputting the results
         */
        public Builder(OutputType outType, String folderName, long numOfBufferedRows, FileWriter out) {
            this.outType = outType;
            this.folderName = folderName;
            this.numOfBufferedRows = numOfBufferedRows;
            this.out = out;
            this.colWidths=makeColHeaderWidths();
            this.colNames=makeColHeaderNames();
        }

        /**
         * Sets whether a progress bar and progress percentage should be printed to the console
         * @param printProgress print the simulation progress to the console
         * @return a reference to this Builder
         */
        public Builder setPrintProgress(boolean printProgress) {
            this.printProgress = printProgress;
            return this;
        }

        /**
         * Sets whether the results should be printed to the console (as opposed to saved to a file)
         * @param printOutputToConsole print the simulation results to the console
         * @return a reference to this Builder
         */
        public Builder setPrintOutput(boolean printOutputToConsole) {
            this.printOutputToConsole = printOutputToConsole;
            return this;
        }

        /**
         * Sets the size (in # of chars) of the buffer used to keep results between writes to the file
         * If not called, the buffer size will be determined by the number of observables and specified column widths
         * @param bufferSize number of characters to keep in the buffer between flushes
         * @return a reference to this Builder
         */
        public Builder setBufferSize(int bufferSize) {
            this.bufferSize = bufferSize;
            return this;
        }

        /**
         * Creates a new {@code OutputWriter} with all configuration options that have been specified by calling methods on this builder.
         * @return the new {@code OutputWriter}
         */
        public OutputWriter build() {
            if (bufferSize==0){
                // includes for each buffered row: column widths, 1 extra space for each column as separator, line separator, 2 extra for safety.
                this.bufferSize = (int)(numOfBufferedRows)*(IntStream.of(colWidths).sum() + colWidths.length + System.lineSeparator().length() + 2);
            }
            return new OutputWriter(this);
        }

        /**
         * Sets the column header of the output table, based on the output type
         * @return an array of column names
         */
        private String[] makeColHeaderNames(){
            // set output table parameters
            String[] colNames;
            switch (outType){
                case VERBOSE:
                    colNames = new String[]{"index", "Magnetization", "Energy", "meanBx", "stdBx", "meanBy", "stdBy", "meanBz", "stdBz", "maxBtrans", "maxBlong", "meanSpinSize", "stdSpinSize", "mk2x", "mk2y", "mk2z", "swap"};
                    break;
                case BIN:
                    colNames = new String[]{"binN", "<|M|>", "<|M|2>", "<M>", "<M2>", "<M2>", "<M22>", "<E>", "<E2>", "<meanBx>", "<meanBx2>",
                                            "<stdBx>", "<stdBx2>", "<meanBy>", "<meanBy2>","<stdBy>", "<stdBy2>", "<meanBz>", "<meanBz2>", "<stdBz>", "<stdBz2>", "<maxBtrans>",
                                            "<maxBtrans2>", "<maxBlong>", "<maxBlong2>", "<meanSpinSize>",  "<meanSpinSize2>", "<stdSpinSize>", "<stdSpinSize2>", "<mk2x>", "<mk2x2>", "<mk2y>", "<mk2y2>","<mk2z>", "<mk2z2>", "swap"};
                    break;
                case SPIN:
                    colNames = new String[]{"n", "spin", "spinSize", "localBx", "localBy", "localBz"};
                    break;
                default:
                    throw new IllegalStateException("Unexpected value: " + outType);
            }
            return colNames;
        }

        /**
         * Sets the column widths (# of characters per column) of the output table, based on the output type
         * @return an array of column widths
         */
        private int[] makeColHeaderWidths(){
            // set output table parameters
            int[] colWidths;
            switch (outType){
                case VERBOSE:
                    colWidths = new int[]{10, 17, 19, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 5};
                    break;
                case BIN:
                    colWidths = new int[]{10, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 10};
                    break;
                case SPIN:
                    colWidths = new int[]{6, 5, 17, 17, 17, 17};
                    break;
                default:
                    throw new IllegalStateException("Unexpected value: " + outType);
            }
            return colWidths;
        }
    }

}
