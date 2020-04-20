package protontherapy;

import java.io.*;
import java.util.Arrays;


/// this import is needed for the file input/output functionality
import java.io.*;

class Voxel extends Parameters 
{ 
    private boolean plotPP;
    private double [] binlow, binhigh;
    private double [] binwidth;
    // changed from private 
    private int nbins;
    private double[][] binCentre;
    private String histname;

    // double array to store the actual histogram data
    private double[][] zSliceEnergyPP;

    private long underflow, overflow;
    private long nfilled;
    
    //3D array of Voxels for generating z slice matrices
    private double[][][] voxels;
    
    //DVH info
    private int DVHnbins;
    private double[][] DVHvalues;
    
    // voxel coords arrays
    private double [] voxelx;
    private double [] voxely;
    private double [] voxelz;

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
        zSliceEnergyPP = new double[ProtonTherapy.energies.length][nbins];
        underflow = 0;
        overflow = 0;
        nfilled = 0;
        
        voxels = new double[nbins][nbins][nbins];
                
        // calculate the centre of each bin for all dimensions
        binCentre = new double[3][nbins];
        for(int a = 0; a < 3; a++){
            for (int i = 0; i < nbins; i++) {
                binCentre[a][i] = binlow[a] + (i+0.5)*binwidth[a];
            }
        }
        DVHnbins = 150;
        DVHvalues = new double[2][DVHnbins];
        
        System.out.println(Arrays.deepToString(binCentre));
    }

    // returns number of bins
    public int getNbins()
    {
        return nbins;
    }
    // returns underflow variable
    public long getUnderflow()
    {
        return underflow;
    }
    // returns overflow variable
    public long getOverflow()
    {
        return overflow;
    }
    // returns number variable
    public long getNfilled()
    {
        return nfilled;
    }
    public int getDVHnbins(){
        return DVHnbins;
    }
    public double[][] getDVHvalues(){
        return DVHvalues;
    }
    public double getVoxelVolume(){
        return (binwidth[0]*binwidth[1]*binwidth[2]);
    }
    

    

    
    
    
//    double energy: the value of energy depositited. 
//    int ke: position value of energy of current itteration in main of ProtonTherapy.
//    If output is normal combined bragg peak (ke = 0)
    public void fill(double energy, Particle p, int ke){
        //fillVoxels(energy, p);
        
        int xBin = (int) ((p.x - binlow[0])/binwidth[0]);
        int yBin = (int) ((p.y - binlow[1])/binwidth[1]);
        int zBin = (int) ((p.z - binlow[2])/binwidth[2]);
        
        if(p.x < binlow[0] || p.y < binlow[1] || p.z < binlow[2]){
            underflow++;
        }else if(p.x > binhigh[0] || p.y > binhigh[1] || p.z > binhigh[2]){
            overflow++;
        }else{
            voxels[zBin][xBin][yBin] = voxels[zBin][xBin][yBin] + energy;
            if(plotPP){
                zSliceEnergyPP[ke][zBin] += energy;
            }
        }
        
    }
    
    // calculates absorbed dose for each Zslice
    public double [] getAbsorbedDose(double [][][] voxels, 
                                    double x0, double y0, double z0,
                                    double x1, double y1, double z1, 
                                    double rhoin) {
        // initialising energy slice array
        double [][] total_energy = new double [2][nbins];
        // calculates volume of each voxel
        double voxelVolume = (x1 - x0)/nbins * (y1 - y0)/nbins * (z1 - z0)/nbins;
        // calculates voxel mass
        double voxelMass = rhoin*voxelVolume;
        // initialises absorbed dose array
        double [] AbsorbedDose = new double [nbins];
        
        // for z slices
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
            // for voxels in y
                for (int k = 0; k < nbins; k++) {
                // summing x and y energy deposited in each z slice
                total_energy[0][i] += voxels[i][j][k];
                }
                }
            // total energy absorbed in slice / total slice mass
            AbsorbedDose[i] = total_energy[0][i]/(voxelMass*Math.pow(nbins, 2));
            
            }
        return AbsorbedDose;
        }   
    
    
//    i : selection for x,y,z
//    int ke: position value of energy of current itteration in main of ProtonTherapy.
    public double getSliceEnergyPP(int ke, int nbin)
    {
        // returns the contents on bin 'nbin' to the user
        return zSliceEnergyPP[ke][nbin];
    }
    
    public void convertVoxelsToIsodose(){
        int centerbin = (int) ((0 - binlow[0])/binwidth[0]);
        int centerZBin = (int) ((this.getTumourCenterPos() - binlow[2])/binwidth[2]);
        double maxValue = 0;
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    if(maxValue<voxels[zi][xi][yi]){
                        maxValue = voxels[zi][xi][yi];
                    }
                }
            }
        }
        //System.out.println("max value: "+maxValue);
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    voxels[zi][xi][yi] = (voxels[zi][xi][yi]/voxels[centerZBin][50][50])*100;                    
                }
            }
        }
    }
    
    public boolean isInSphere(int xBin, int yBin, int zBin){
        int centerBin = (int) ((this.getTumourCenterPos() - binlow[2])/binwidth[2]);
        double rad = 0.018;
        double radVoxel = Math.sqrt(Math.pow(binCentre[0][xBin], 2)
                                +Math.pow(binCentre[1][yBin], 2)
                                +Math.pow(binCentre[2][zBin]-binCentre[2][centerBin], 2));
        return (radVoxel<rad);
    }
    
    public double[][][] generateOtherDVH(){
        Histogram[] DVH = generateDVH();
        double[] totalVoxels = new double[2];
        double[][][] DVHs = new double[2][2][DVHnbins];
        //double[][] DVH = new double[2][DVHnbins];
        
        for(int a = 0; a < DVH.length; a++){
            for(int i = 0; i < DVHnbins; i++){
                DVHvalues[a][i] = DVH[a].getContent(i);
                totalVoxels[a] += DVH[a].getContent(i); 
            }
        }
        for(int a = 0; a < DVH.length; a++){
            for(int i = 0; i < DVHnbins; i++){
                DVHs[0][a][i] = DVHvalues[a][i]/totalVoxels[a];
                if(i==0){
                    DVHs[1][a][i] = 1-DVHs[0][a][i];
                }else{
                    DVHs[1][a][i] = DVHs[1][a][i-1]-DVHs[0][a][i]; 
                }
            }
        }
        return DVHs;
    }
    
    public Histogram[] generateDVH(){
        Histogram DVHTumour = new Histogram(DVHnbins, 0, 150, "DVHTumour");
        Histogram DVHBody = new Histogram(DVHnbins, 0, 150, "DVHBody");

        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    if(isInSphere(xi,yi,zi)){
                        DVHTumour.fill((voxels[zi][xi][yi]/voxels[zi][50][50])*100);                      
                    }else{
                        DVHBody.fill((voxels[zi][xi][yi]/voxels[zi][50][50])*100);
                    }
                }
            }
        }
        Histogram[] hists = {DVHTumour, DVHBody};        
        return hists;
    }
        
    //-------------------------------------
    public void print()
    {
        for(int a = 0; a < 3; a++){
            for (int bin = 0; bin < getNbins(); bin++) {
                System.out.println("Bin " + bin + " = " +getSliceEnergyPP(0,bin));
            }
            System.out.println("The number of fills = " + getNfilled());
            System.out.println("Underflow = " + getUnderflow()
                + ", Overflow = " + getOverflow());
        }

    }
    
    //-------------------------------------
//    This function is called so that it makes use of the output fucntions below 
//    makes it easier for dumping data.
    public void writeData(double depth, String filename){
        //convertVoxelsToIsodose();
        writeSOBP();
        writeProfile();
        writeZSlice(depth);
        writeDVH();
        writeToDiskPP("bragg_peaks.csv");
        writeInfo();
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
//        outputFile.print("z,");
//        for(int i = 0; i < ProtonTherapy.energies.length; i++){
//            outputFile.print("PP #"+(i+1)+" ("+ProtonTherapy.energies[i][0]+"MeV),");
//        }
//        outputFile.println();
        for(int n = 0; n < nbins; n++){
            outputFile.print(binCentre[2][n]+",");
            for(int i = 0; i < ProtonTherapy.energies.length; i++){
                outputFile.print(getSliceEnergyPP(i, n)+",");
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
                outputFile.print(voxels[zSliceNum][a][i] +",");
            }
            outputFile.println();
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
    
    // write file with the absorbed dose for each z slice
    public void writeDose() {
        double [] absorbed_dose = new double [nbins];
        String filename = "AbsorbedDose.csv";
        
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file "+filename + ". Histogram data was not saved.");
            return;
            
        }
        // calculating absorbed dose and writing absorbed dose file
        absorbed_dose = getAbsorbedDose(voxels, -0.20, -0.20, 0.22,            // start x, y, z
                             0.20, 0.20, 0.72, 1);
        System.out.println("horse");
        System.out.println(Arrays.toString(absorbed_dose));
        for(int i = 0; i < nbins; i++) {
            outputFile.print(absorbed_dose[i] +",");
        }
        
        outputFile.println();
        outputFile.close();
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
                value = value + voxels[a][i][nbins/2];
            }
            outputFile.println(binCentre[0][i]+","+ value);            
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
    
    public void writeSOBP(){
        String filename = "SOBP.csv";
        PrintWriter outputFile;
        try{
            outputFile = new PrintWriter(filename);
        } catch(IOException e){
            System.err.println("Failed to open file."+filename+" Histogram data was not saved");
            return;
        }
        for(int zi = 0; zi < nbins; zi++){
            double value = 0;
            for(int xi = 0; xi < nbins; xi++){
                for(int yi = 0; yi < nbins; yi++){
                    value += voxels[zi][xi][yi];
                } 
            }
            outputFile.println(binCentre[2][zi]+","+value);
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");        
    }
    
    public void writeDVH(){
        double[][][] DVHs = generateOtherDVH();
        String filename = "DVH.csv";
        PrintWriter outputFile;
        try{
            outputFile = new PrintWriter(filename);
        } catch(IOException e){
            System.err.println("Failed to open file."+filename+" Histogram data was not saved");
            return;
        }
        
        for(int i = 0; i < DVHnbins; i++){
            outputFile.println(i+","+DVHs[0][0][i]+","+DVHs[1][0][i]+","+DVHs[0][1][i]+","+DVHs[1][1][i]);
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!"); 
    }
    
    public void writeInfo(){
        //double TCP = this.getTCP(this.getDVHvalues(), this.getVoxelVolume());
        
        System.out.println("*** START INFO OUTPUT ***");
        //System.out.println("TCP: "+TCP);
    }

}
