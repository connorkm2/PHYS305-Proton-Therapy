package protontherapy;

import java.io.*;
import java.util.Arrays;
import java.util.stream.*;

public class RBE extends Parameters {
    
    double rho = 1000;
   
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
    private double[][][] LET;
    
    // initialisng arrays to store dE and dz for each voxel.
    double [][][] dE = new double [nbins][nbins][nbins];
    double [][][] dEdz = new double [nbins][nbins][nbins];
    
    
    //DVH info
    private int DVHnbins;
    private double[][] DVHvalues;    
    
    public RBE(int numberOfBins, double [] start, double [] end, String name) {
        // store the parameters and setup the histogram
        // note that parameters need to have different names than class variables
        nbins = numberOfBins;
        
        binlow = start; // beginning coordinates
        binhigh = end; // end coordinates
        
        histname = name;
                
        binwidth = new double[3];
        
        // variables 
        for(int i = 0; i < binlow.length; i++){
            binwidth[i] = (binhigh[i] - binlow[i]) / (double) nbins;
        }
        //zSliceEnergyPP = new double[ProtonTherapy.energies.length][nbins];
        underflow = 0;
        overflow = 0;
        nfilled = 0;
        
        voxels = new double [nbins][nbins][nbins];
        LET = new double [nbins][nbins][nbins];
        dE = new double [nbins][nbins][nbins];
        dEdz = new double [nbins][nbins][nbins];
                
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
    
    // returns data from LET hist
    public double [] getLET(int ke, int nbin)
    {
        // returns the contents on bin 'nbin' to the user
        return LET[ke][nbin];
    }

    // need to modify this for calculating LET 
    public void fill(double LETParticle, Particle p, int ke){
    //fillVoxels(energy, p);

    int xBin = (int) ((p.x - binlow[0])/binwidth[0]);
    int yBin = (int) ((p.y - binlow[1])/binwidth[1]);
    int zBin = (int) ((p.z - binlow[2])/binwidth[2]);

    if(p.x < binlow[0] || p.y < binlow[1] || p.z < binlow[2]){
        underflow++;
    }else if(p.x > binhigh[0] || p.y > binhigh[1] || p.z > binhigh[2]){
        overflow++;
    }else{
        LET[zBin][xBin][yBin] = LET[zBin][xBin][yBin] + LETParticle;
    }
    }
    
    // ***how do I reference this directly from voxel class?
        // checks if particle is in tumour
    public boolean isInSphere(int xBin, int yBin, int zBin){
    int centerBin = (int) ((this.getTumourCenterPos() - binlow[2])/binwidth[2]);
    double rad = 0.018;
    double radVoxel = Math.sqrt(Math.pow(binCentre[0][xBin], 2)
                            +Math.pow(binCentre[1][yBin], 2)
                            +Math.pow(binCentre[2][zBin]-binCentre[2][centerBin], 2));
    return (radVoxel<rad);
    }
    

    // creates histogram of LET energies
    public Histogram [] LET(double [] EnergyLossArray, double [] distance, String filename){
        Histogram LET = new Histogram(DVHnbins, 0, 150, "LET");
        
        // calculating LET for each particle 
        double totalEloss = finddE(EnergyLossArray);
        double dist = findTrackLength(distance);
        double dE_dz = totalEloss/dist;
        double LETParticle = (totalEloss*dE_dz*(1/rho))/(totalEloss);
    
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    // if particle is in tumour volume, fill LET hist
                    if(isInSphere(xi,yi,zi)){
                        LET.fill(LETParticle);                      
                    }
                }
            }
        }
        Histogram[] hists = {LET};        
        return hists;
    }
    
    
    public void writeLET(int ke, String var)
    {
        String filename = var+".csv";
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file " + filename + ". Histogram data was not saved.");
            return;
        }

        // now make a loop to write the contents of each bin to disk, one number at a time
        // together with the x-coordinate of the centre of each bin.
        for (int n = 0; n < nbins; n++) {
            // comma separated values
            outputFile.println(n + "," + binCentre[n] + "," + getLET(ke, n));
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
 

    
    
    
//// Calculating RBE for Different Models ///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
// CASE 1 - RBE = 1.1
    
    public double [][][] simpleRBE(double RBE, double [][][] dE, int nbins) {
       
        // initialising RBE weighted dose for simple case
        double [][][] simpleRBEWeight = new double [nbins][nbins][nbins];
        
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                    simpleRBEWeight[i][j][k] = RBE*dE[i][j][k];
                }
            }
        }
        
        return simpleRBEWeight;
    }
    
// CURRENTLY ALL METHODS BELOW ARE VOID - NEED TO WORK OUT HOW TO INCLUDE FILL FUNCTION
// CASE 2 - CARABE-FERNANDEZ MODEL
    
    // calculates RBE max and associated weighted dose
    public void CarFerRBEmax(double [][][] dE, double [][][] LET, double alpha_beta, int nbins) {
        // initialising arrays
        double [][][] RBEmax = new double [nbins][nbins][nbins];
        double [][][] maxCarFerRBEWeight = new double [nbins][nbins][nbins];
        
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                RBEmax[i][j][k] = 0.834 + 0.154*(2.686/(alpha_beta))*LET[i][j][k];
                maxCarFerRBEWeight[i][j][k] = RBEmax[i][j][k]*dE[i][j][k];
    }
            }
        }
        
    }
    
    // calculates RBE min and associated weighted dose
    public void CarFerRBEmin(double [][][] dE, double [][][] LET, double alpha_beta, int nbins) {
        // initialising arrays
        double [][][] RBEmin = new double [nbins][nbins][nbins];
        double [][][] minCarFerRBEWeight = new double [nbins][nbins][nbins];
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                RBEmin[i][j][k] = 0.834 + 0.154*(2.686/(alpha_beta))*LET[i][j][k];
                minCarFerRBEWeight[i][j][k] = RBEmin[i][j][k]*dE[i][j][k];
        
                }
            }
        }
        
    }
   

        
    // CASE 3 - WEDENBERG MODEL
    
    public void WendelRBEmin(double [][][] dE, double [][][] LET, double alpha_beta, int nbins, double RBEmin) {
        // initialising arrays
        double [][][] minWendelRBEWeight = new double [nbins][nbins][nbins];
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                minWendelRBEWeight[i][j][k] = RBEmin*dE[i][j][k];
        
                }
            }
        }
        
    }

        public void WendelRBEmax(double [][][] dE, double [][][] LET, double alpha_beta, int nbins) {
        // initialising arrays
        double [][][] maxWendelRBEWeight = new double [nbins][nbins][nbins];
        double [][][] RBEmax = new double [nbins][nbins][nbins];
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                RBEmax[i][j][k] = 1.00 + (0.434/(alpha_beta))*LET[i][j][k];
                maxWendelRBEWeight[i][j][k] = RBEmax[i][j][k]*dE[i][j][k];
        
                }
            }
        }
        
        }
        
    
    public double [] getPhantomStart(double [] phantomPosition) {
           double [] start = Arrays.copyOfRange(phantomPosition, 0, 2); 
           return start;
    }

    public double [] getPhantomEnd(double [] phantomPosition) {
            double [] end = Arrays.copyOfRange(phantomPosition, 3, 5);
           return end;
    }

    // returns energy lost by single particle in said voxel
    public double finddE (double [] EnergyLossArray) {
        double totalEloss = 0;
        // summing all elements in array
        for (int i = 0; i < EnergyLossArray.length; i++) {
            totalEloss += EnergyLossArray[i];
        }
        return totalEloss;
    }

    // returns total distance travelled by particle
    public double findTrackLength(double [] stepsize) {
        double totalTrackLength = 0;
        // summing all elements in array
        for (int i = 0; i < stepsize.length; i++) {
            totalTrackLength += stepsize[i];
        }
        return totalTrackLength;  
    }    

        // returns voxel x coordinates
    public double [] getVoxelx(double [] phantomPosition, int nbins){
        // initialise voxel x array
        double [] voxelx = new double [nbins+1];
        // returns start and end x values
        double startx = phantomPosition[0];
        double endx = phantomPosition[3];
        // calculates step length in x
        double stepdist = (phantomPosition[3] - phantomPosition[0])/nbins; 
        // fills voxelx array with x coords
        for (int i = 0; i < nbins; i++){
            voxelx[i] = startx + i*stepdist;
        }
        return voxelx;
    }
    
    // returns voxel y coordinates
    public double [] getVoxely(double [] phantomPosition, int nbins){
        // initialise voxel y array
        double [] voxely = new double [nbins+1];
        // returns start and end y values
        double starty = phantomPosition[1];
        double endy = phantomPosition[4];
        // calculates step length in y
        double stepdist = (phantomPosition[4] - phantomPosition[1])/nbins; 
        // fills voxely array with y coords
        for (int i = 0; i < nbins; i++){
            voxely[i] = starty + i*stepdist;
        }
        return voxely;
    }
    
    // returns voxel z coordinates
    public double [] getVoxelz(double [] phantomPosition, int nbins){
        // initialise voxel z array
        double [] voxelz = new double [nbins+1];
        // returns start and end z values
        double startz = phantomPosition[2];
        double endz = phantomPosition[5];
        // calculates step length in z
        double stepdist = (phantomPosition[5] - phantomPosition[2])/nbins; 
        // fills voxelz array with z coords
        for (int i = 0; i < nbins; i++){
            voxelz[i] = startz + i*stepdist;
        }
        return voxelz;
    }
}
