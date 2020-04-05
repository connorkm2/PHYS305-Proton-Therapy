package protontherapy;

import java.io.*;
import java.util.Arrays;


/// this import is needed for the file input/output functionality
import java.io.*;

class Voxel
{ 
    private boolean plotPP;
    private double [] binlow, binhigh;
    private double [] binwidth;
    private int nbins;
    private double[][] binCentre;
    private String histname;

    // double array to store the actual histogram data
    private double[][][] sumBinEnergy;

    private long [] underflow, overflow;
    private long [] nfilled;
    
    //3D array of Voxels for generating z slice matrices
    private double[][][] zSlices;

    // constructor for the class Histogram
    public Voxel(int numberOfBins, double [] start, double [] end, String name, boolean plotBraggPeaks)
    {
        // store the parameters and setup the histogram
        // note that parameters need to have different names than class variables
        nbins = numberOfBins;
        
        // switches some functions to plot individual bragg peaks
        plotPP = plotBraggPeaks;
        
        binlow = start; // beginning coordinates
        binhigh = end; // end coordinates
        
        histname = name;
                
        binwidth = new double[3];
        
        // variables 
        for(int i = 0; i < binlow.length; i++){
            binwidth[i] = (binhigh[i] - binlow[i]) / (double) nbins;
        }
        sumBinEnergy = new double[ProtonTherapy.energies.length][3][nbins];
        underflow = new long[3];
        overflow = new long[3];
        nfilled = new long[3];
        
        zSlices = new double[nbins][nbins][nbins];
                
        // calculate the centre of each bin for all dimensions
        binCentre = new double[3][nbins];
        for(int a = 0; a < 3; a++){
            for (int i = 0; i < nbins; i++) {
                binCentre[a][i] = binlow[a] + (i+0.5)*binwidth[a];
            }
        }
    }
    // returns number of bins
    public int getNbins()
    {
        return nbins;
    }
    // returns underflow variable
    public long getUnderflow(int i)
    {
        return underflow[i];
    }
    // returns overflow variable
    public long getOverflow(int i)
    {
        return overflow[i];
    }
    // returns number variable
    public long getNfilled(int i)
    {
        return nfilled[i];
    }
    
//    double energy: the value of energy depositited. 
//    int ke: position value of energy of current itteration in main of ProtonTherapy.
//    If output is normal combined bragg peak (ke = 0)
    public void fill(double energy, Particle p, int ke){
        fillzSlices(energy, p);
        
        if(!plotPP){ ke = 0; }
        double [] position = {p.x, p.y, p.z};
        
        for(int i = 0; i < binlow.length; i++){
            if (position[i] < binlow[i]) {
                underflow[i]++;
                // increases overflow
            } else if (position[i] >= binhigh[i]) {
                overflow[i]++;
                // filling bins 
            } else {

                int ibin = (int) ((position[i] - binlow[i])/binwidth[i]);
                sumBinEnergy[ke][i][ibin] = sumBinEnergy[ke][i][ibin] + energy;//*(ProtonTherapy.energies[ke][1]);
            }
                // increases filled bin count by one
                nfilled[i]++;
            }
        
    }
    
//    Will fill an array of slices of z with the sum of energy deposited 
//    in those voxels. 
//    Essentially uses the histogram bins to genertate a 3D matrix of values at 
//    a coordinate of x,y,z
    public void fillzSlices(double energy, Particle p){
        int xBin = (int) ((p.x - binlow[0])/binwidth[0]);
        int yBin = (int) ((p.y - binlow[1])/binwidth[1]);
        int zBin = (int) ((p.z - binlow[2])/binwidth[2]);
                
        zSlices[zBin][xBin][yBin] = zSlices[zBin][xBin][yBin] + energy;      
    }
    
//    i : selection for x,y,z
//    int ke: position value of energy of current itteration in main of ProtonTherapy.
    public double getBinEnergy(int ke, int i, int nbin)
    {
        // returns the contents on bin 'nbin' to the user
        return sumBinEnergy[ke][i][nbin];
    }
        
    //-------------------------------------
    public void print()
    {
        for(int a = 0; a < 3; a++){
            for (int bin = 0; bin < getNbins(); bin++) {
                System.out.println("Bin " + bin + " = " +getBinEnergy(0, a, bin));
            }
            System.out.println("The number of fills = " + getNfilled(a));
            System.out.println("Underflow = " + getUnderflow(a)
                + ", Overflow = " + getOverflow(a));
        }

    }
    
    //-------------------------------------
    public void writeToDisk(String filename)
    {
        String [] plane = {"x","y","z"};
        for(int i = 0; i < 3; i++){
            // this sends the output to a file with name "filename"
            // the block with try { ... } catch (IOException e) { ... } is needed to handle the case,
            // where opening the file fails, e.g. disk is full or similar
            PrintWriter outputFile;
            try {
                outputFile = new PrintWriter(plane[i]+"_"+filename);
            } catch (IOException e) {
                System.err.println("Failed to open file " + plane[i]+"_"+filename + ". Histogram data was not saved.");
                return;
            }

            // Write the file as a comma seperated file (.csv) so it can be read it into EXCEL
            // first some general information about the histogram
            outputFile.println("histname, " + histname);
            outputFile.println("binlow, " + binlow[i]);
            outputFile.println("binwidth, " + binwidth[i]);
            outputFile.println("nbins, " + nbins);
            outputFile.println("underflow, " + underflow[i]);
            outputFile.println("overflow, " + overflow[i]);

            // now make a loop to write the contents of each bin to disk, one number at a time
            // together with the x-coordinate of the centre of each bin.
            for (int n = 0; n < nbins; n++) {
                // comma separated values
                outputFile.println(n + "," + binCentre[i][n] + "," + getBinEnergy(0, i, n));
            }
            outputFile.close(); // close the output file
            System.out.println(plane[i]+"_"+filename+" written!");
        }
    }
    
//    This function is simply an alternative to the original writeToDisk above. It 
//    combines x,y,z to one file to make plotting easier when debugging.
    public void writeToDiskCombined(String filename){
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file "+filename + ". Histogram data was not saved.");
            return;
        }
        
        //outputFile.println("x, Energy, y, Energy, z, Energy");
        for(int n = 0; n < nbins; n++){
            for(int i = 0; i < 3; i++){
                outputFile.print(binCentre[i][n] + "," + getBinEnergy(0, i, n)+",");
            }
            outputFile.println();
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");

    }
    
//    This fucntion outputs the individual pristine peaks of the simulation.
//    Can be useful for investigating weightings and also for nice plots.
//    This fucntion is only called if plotPP is true. This is set when creating 
//    an instance on Voxel class.
    public void writeToDiskPP(String filename){
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file "+filename + ". Histogram data was not saved.");
            return;
        }
        outputFile.print("z,");
        for(int i = 0; i < ProtonTherapy.energies.length; i++){
            outputFile.print("PP #"+(i+1)+" ("+ProtonTherapy.energies[i][0]+"MeV),");
        }
        outputFile.println();
        for(int n = 0; n < nbins; n++){
            outputFile.print(binCentre[2][n]+",");
            for(int i = 0; i < ProtonTherapy.energies.length; i++){
                outputFile.print(getBinEnergy(i, 2, n)+",");
            }
            outputFile.println();
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
        
//    Outputs a csv file which consists of enrgy deposition data in a grid of x,y 
//    for a particular slice of z (specific bin value). 
//    double depth: the value of z in meters of the desired slice.
    public void writeZSlice(double depth){
        int zSliceNum = (int) ((depth - binlow[2])/binwidth[2]);
        
        String filename = zSliceNum+"_zSlice.csv";
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file "+filename + ". Histogram data was not saved.");
            return;
        }
        for(int i = 0; i < nbins; i++){
            for(int a = 0; a < nbins; a++){
                outputFile.print(zSlices[zSliceNum][a][i] +",");
            }
            outputFile.println();
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
    
    public void writeProfile(){
        String filename = "dose_profile.csv";
        PrintWriter outputFile;
        try{
            outputFile = new PrintWriter(filename);
        } catch(IOException e){
            System.err.println("Failed to open file."+filename+" Histogram data was not saved");
            return;
        }
        for(int i = 0; i < nbins; i++){
            double value = 0;
            for(int a = 0; a < nbins; a++){
                value = value + zSlices[a][i][nbins/2];
            }
            outputFile.println(binCentre[0][i]+","+ value);            
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
    
//    This function is called so that it makes use of the output fucntions above 
//    makes it easier for dunping data.
    public void writeData(double depth, String filename){
        if(plotPP == true){
            writeToDiskPP("bragg_peaks.csv");
        }else{
            writeToDiskCombined(filename);
        }
        writeProfile();
        writeZSlice(depth);

    }

}
