package simulation.montecarlo;

import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Formatter;
import java.util.InputMismatchException;
import java.util.stream.IntStream;

public class OutputWriter implements Closeable {
    private final boolean printProgress, printOutputToConsole;
    private final long obsPrintSweepNum;
    private final String folderName;
    private final boolean verboseOutput;
    private final int[] colWidths;
    private final String[] colNames;
    private final int bufferSize;
    private final FileWriter out;

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
        out.close();
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

    public String makeTableRowFormat(char[] colTypes){
        if (colWidths.length!=colTypes.length)
            throw new InputMismatchException("array of column types should be the same length as array of column widths.");
        StringBuilder rowFormat = new StringBuilder(colTypes.length*9);
        for (int i=0;i<colTypes.length-1;i++){	// the last column is printed separately
            // % 10d
            if (colTypes[i]=='g')
                rowFormat.append("% "+(colWidths[i]-1)+'.'+(colWidths[i]-8)+colTypes[i]+' ');
            else if (colTypes[i]=='d')
                rowFormat.append("% "+(colWidths[i])+colTypes[i]+' ');
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


    public static void makeDir(String path, String dirName){
        dirName = path.concat(dirName);
        File directory = new File(dirName);
        if (! directory.exists()){
            directory.mkdir();
        }
    }

    public OutputWriter(Builder builder) {
        this.printProgress = builder.printProgress;
        this.printOutputToConsole = builder.printOutputToConsole;
        this.obsPrintSweepNum = builder.obsPrintSweepNum;
        this.folderName = builder.folderName;
        this.verboseOutput = builder.verboseOutput;
        this.colNames = builder.colNames;
        this.colWidths = builder.colWidths;
        this.bufferSize = builder.bufferSize;
        this.out = builder.out;
    }

    public static class Builder {
        // required
        private final boolean verboseOutput;
        private final String folderName;
        private final long obsPrintSweepNum;
        private final FileWriter out;

        // optional - default values below
        private boolean printProgress=false;
        private boolean printOutputToConsole=false;
        // optional - default values to be determined later based on other input
        private int bufferSize=0;

        // created without outside input
        private final int[] colWidths=makeColHeaderWidths();
        private final String[] colNames=makeColHeaderNames();

        public Builder(boolean verboseOutput, String folderName, long obsPrintSweepNum, FileWriter out) {
            this.verboseOutput = verboseOutput;
            this.folderName = folderName;
            this.obsPrintSweepNum=obsPrintSweepNum;
            this.out = out;
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
                this.bufferSize = (int)(obsPrintSweepNum)*(IntStream.of(colWidths).sum() + 4);	// +4 includes 2 for the line separator and 2 as extra buffer just in case
            }
            return new OutputWriter(this);
        }

        private String[] makeColHeaderNames(){
            // set output table parameters
            String[] colNames;
            if (verboseOutput) {
                colNames = new String[]{"index", "Magnetization", "Energy", "meanBx", "stdBx", "meanBy", "stdBy", "meanBz", "stdBz", "maxBtrans", "maxBlong", "meanSpinSize", "stdSpinSize", "mk2", "swap"};
            }else{
                colNames = new String[]{"binN", "<|M|>", "<|M|2>", "<M>", "<M2>", "<M2>", "<M22>", "<E>", "<E2>", "<meanBx>", "<meanBx2>",
                        "<stdBx>", "<stdBx2>", "<meanBy>", "<meanBy2>","<stdBy>", "<stdBy2>", "<meanBz>", "<meanBz2>", "<stdBz>", "<stdBz2>", "<maxBtrans>",
                        "<maxBtrans2>", "<maxBlong>", "<maxBlong2>", "<meanSpinSize>",  "<meanSpinSize2>", "<stdSpinSize>", "<stdSpinSize2>", "<mk2>", "<mk22>", "swap"};
            }
            return colNames;
        }

        private int[] makeColHeaderWidths(){
            // set output table parameters
            int[] colWidths;
            if (verboseOutput) {
                colWidths = new int[]{10, 17, 19, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 5};
            }else{
                colWidths = new int[]{10, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 10};
            }
            return colWidths;
        }
    }

}
