package simulation.montecarlo;

import java.io.Closeable;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Formatter;
import java.util.InputMismatchException;
import java.util.stream.IntStream;

public class OutputWriter implements Closeable {
    private final boolean printProgress, printOutputToConsole;
    private final long numOfBufferedRows;
    private final String folderName;
    private OutputType outType;
    private final int[] colWidths;
    private final String[] colNames;
    private final int bufferSize;
    private final FileWriter out;
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

    public void flush(){
        try {
            print(outputBuffer.toString(), false);
            if (!printOutputToConsole) out.flush();
            if (outputBuffer.length() > bufferSize)
                System.err.println("output buffer exceeded its buffer with char count of " + outputBuffer.length() + '/' + bufferSize);
            outputBuffer = null;
            outputBuffer = new StringBuilder(bufferSize);
        }catch (IOException e){
            throw new RuntimeException("Error saving results: " + e.getMessage());
        }
    }

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


    public String makeTableHeader(){
        if (colNames.length!=colNames.length)
            throw new InputMismatchException("array of column names should be the same length as array of column widths.");
        int totalChrsForRow = IntStream.of(colWidths).sum();
        StringBuilder headerRow = new StringBuilder(totalChrsForRow+10);
        for (int i=0;i<colNames.length;i++){
            if (colNames[i].length()>colWidths[i]){
                throw new InputMismatchException("the column header "+colNames[i]+" is longer than its designated" +
                        "column width of "+colWidths[i]);
            } else {
                for (int pad=0;pad<colWidths[i]-colNames[i].length();pad++){
                    headerRow.append(" ");
                }
                headerRow.append(colNames[i]);
            }
        }
        if (headerRow.length()>totalChrsForRow+10)
            System.err.println("header row builder not allocated enough space: " + headerRow.length() + "/" + (totalChrsForRow+10));
        return headerRow.toString();
    }

    public void writeObservablesVerbose(final long sweeps, final double m, final double currentEnergy, final double magField0, final double magField1, final double magField2,
                                 final double magField3, final double magField4, final double magField5, final double magField6, final double magField7, final double magField8,
                                 final double spinSizes0, final double spinSizes1, final double mk2, final boolean lastSwapAccepted){
        if (outType == OutputType.VERBOSE) {
            Formatter formatter = new Formatter(outputBuffer);
            formatter.format(makeTableRowFormat(new char[]{'d', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'c'}), sweeps, m, currentEnergy, magField0, magField1, magField2, magField3, magField4, magField5, magField6, magField7, magField8, spinSizes0, spinSizes1, mk2, (lastSwapAccepted ? '1' : '0'));
            outputBuffer.append(System.lineSeparator());
        }
    }

    public void writeObservablesPerSpin(final int n, final int spin, final double spinSize, final double localBx, final double localBy, final double localBz){
        if (outType == OutputType.SPIN) {
            Formatter formatter = new Formatter(outputBuffer);
            formatter.format(makeTableRowFormat(new char[]{'d', 'd', 'g', 'g', 'g'}), n, spin, spinSize, localBx, localBy, localBz);
            outputBuffer.append(System.lineSeparator());
        }
    }

    public void writeObservablesBin(final long currentBinCount, final double[] binAvg, final double acceptanceRateForBin){
        if (outType == OutputType.BIN){
            Formatter formatter = new Formatter(outputBuffer);
            formatter.format("% "+colWidths[0]+"d %s % "+(colWidths[colWidths.length-1]-2)+".2f%n", currentBinCount, avgArrToString(binAvg, currentBinCount, 1) , acceptanceRateForBin);
        }
    }

    public void writeSwapAcceptance(final boolean swapAccepted){
        if (outType == OutputType.VERBOSE){
            outputBuffer.append(new String(new char[colWidths[colWidths.length-1]-2]).replace('\0', ' ') + (swapAccepted ? '1' : '0'));
            outputBuffer.append(System.lineSeparator());
        }
    }

    public String makeTableRowFormat(char[] colTypes){
        if (colWidths.length!=colTypes.length)
            throw new InputMismatchException("array of column types should be the same length as array of column widths. colWidths.length=" + colWidths.length + ", colTypes.length="+colTypes.length);
        StringBuilder rowFormat = new StringBuilder(colTypes.length*10);
        for (int i=0;i<colTypes.length;i++){	// the last column is printed separately
            // % 10d
            if (colTypes[i]=='g')
                rowFormat.append("% "+(colWidths[i]-1)+'.'+(colWidths[i]-8)+colTypes[i]+' ');
            else if (colTypes[i]=='d')
                rowFormat.append("%"+(colWidths[i])+colTypes[i]+' ');
            else if (colTypes[i]=='c')
                rowFormat.append("%"+(colWidths[i]-1)+colTypes[i]+' ');
            else
                throw new IllegalArgumentException("column types can be only 'd' or 'g'");
        }
//		System.out.println(rowFormat.toString());
        if (rowFormat.length()>colTypes.length*9)
            System.err.println("string format builder not allocated enough space: "+ rowFormat.length()+'/'+(colTypes.length*9));
        return rowFormat.toString();
    }


    public String avgArrToString(double[] avgArr, long N, int startIndex) {
        if (avgArr.length+startIndex>colWidths.length)
            throw new InputMismatchException("start index too high such that startIndex+avgArr.length > colWidths.length.");
        int sum = IntStream.of(colWidths).sum();
        StringBuilder str = new StringBuilder(sum+10);
        Formatter formatter = new Formatter(str);
        int i;
        for (i=0;i<avgArr.length-1;i++){
            formatter.format("% "+(colWidths[startIndex+i]-1)+'.'+(colWidths[startIndex+i]-8)+"g ",(avgArr[i]/N));
        }
        formatter.format("% "+(colWidths[startIndex+i]-1)+'.'+(colWidths[startIndex+i]-8)+"g ",(avgArr[i]/N));
        //System.out.println("String length: " + str.length());
        if (str.length()>sum+10)
            System.err.println("avg arr row builder not allocated enough space: " + str.length() + "/" + (sum+10));
        return str.toString();
    }

    public int getBufferSize(){ return this.bufferSize; }

    public String getFolderName() { return this.folderName; }

    public OutputType getOutType() { return outType;
    }

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

        public Builder(OutputType outType, String folderName, long numOfBufferedRows, FileWriter out) {
            this.outType = outType;
            this.folderName = folderName;
            this.numOfBufferedRows = numOfBufferedRows;
            this.out = out;
            this.colWidths=makeColHeaderWidths();
            this.colNames=makeColHeaderNames();
        }

        public Builder setPrintProgress(boolean printProgress) {
            this.printProgress = printProgress;
            return this;
        }

        public Builder setPrintOutput(boolean printOutputToConsole) {
            this.printOutputToConsole = printOutputToConsole;
            return this;
        }

        public Builder setBufferSize(int bufferSize) {
            this.bufferSize = bufferSize;
            return this;
        }

        public OutputWriter build() {
            if (bufferSize==0){
                this.bufferSize = (int)(numOfBufferedRows)*(IntStream.of(colWidths).sum() + 4);	// +4 includes 2 for the line separator and 2 as extra buffer just in case
            }
            return new OutputWriter(this);
        }

        private String[] makeColHeaderNames(){
            // set output table parameters
            String[] colNames;
            switch (outType){
                case VERBOSE:
                    colNames = new String[]{"index", "Magnetization", "Energy", "meanBx", "stdBx", "meanBy", "stdBy", "meanBz", "stdBz", "maxBtrans", "maxBlong", "percBz", "meanSpinSize", "stdSpinSize", "mk2", "swap"};
                    break;
                case BIN:
                    colNames = new String[]{"binN", "<|M|>", "<|M|2>", "<M>", "<M2>", "<M2>", "<M22>", "<E>", "<E2>", "<meanBx>", "<meanBx2>",
                                            "<stdBx>", "<stdBx2>", "<meanBy>", "<meanBy2>","<stdBy>", "<stdBy2>", "<meanBz>", "<meanBz2>", "<stdBz>", "<stdBz2>", "<maxBtrans>",
                                            "<maxBtrans2>", "<maxBlong>", "<maxBlong2>", "<percBz>", "<percBz2>", "<meanSpinSize>",  "<meanSpinSize2>", "<stdSpinSize>", "<stdSpinSize2>", "<mk2>", "<mk22>", "swap"};
                    break;
                case SPIN:
                    colNames = new String[]{"n", "spin", "spinSize", "localBx", "localBy", "localBz"};
                    break;
                default:
                    throw new IllegalStateException("Unexpected value: " + outType);
            }
            return colNames;
        }

        private int[] makeColHeaderWidths(){
            // set output table parameters
            int[] colWidths;
            switch (outType){
                case VERBOSE:
                    colWidths = new int[]{10, 17, 19, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 5};
                    break;
                case BIN:
                    colWidths = new int[]{10, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 10};
                    break;
                case SPIN:
                    colWidths = new int[]{10, 5, 17, 17, 17, 17};
                    break;
                default:
                    throw new IllegalStateException("Unexpected value: " + outType);
            }
            return colWidths;
        }
    }

}
