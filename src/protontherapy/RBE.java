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
    
    // returns data from LET hist
    public double getContent(int zBin, int xBin, int yBin)
    {
        // returns the contents on bin 'nbin' to the user
        return voxels[zBin][xBin][yBin];
    }
    

    // fills voxels with values
    public void fill(double LETParticle, Particle p, int ke){
    //fillVoxels(energy, p);
    voxels = new double[nbins][nbins][nbins];

    int xBin = (int) ((p.x - binlow[0])/binwidth[0]);
    int yBin = (int) ((p.y - binlow[1])/binwidth[1]);
    int zBin = (int) ((p.z - binlow[2])/binwidth[2]);

    if(p.x < binlow[0] || p.y < binlow[1] || p.z < binlow[2]){
        underflow++;
    }else if(p.x > binhigh[0] || p.y > binhigh[1] || p.z > binhigh[2]){
        overflow++;
    }else{
        voxels[zBin][xBin][yBin] = voxels[zBin][xBin][yBin] + LETParticle;
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
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
            // comma separated values
            outputFile.println(zi + "," + binCentre[zi] + "," + getContent(zi, xi, yi));
        }
        }
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
 

//// Calculating RBE for Different Models ///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
// CASE 1 - RBE = 1.1
    
    // returns RBE weighted dose
    public double getSimpleRBEweight(double totalEloss) {
                double RBE = 1.1;
                double simpleRBEweight = RBE*totalEloss;      
        return simpleRBEweight;
    }
    
    // returns RBE weighted dose histogram for case 1
    public Histogram [] simpleRBEHist(int ke, int nbin, double [] EnergyLossArray) {
        Histogram simpleRBE = new Histogram(DVHnbins, 0, 150, "Simple RBE");
        
        double totalEloss = finddE(EnergyLossArray);
        double simpleRBEweight = getSimpleRBEweight(totalEloss);
        
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    // if particle is in tumour volume, fill LET hist
                    if(isInSphere(xi,yi,zi)){
                        simpleRBE.fill(simpleRBEweight);
                    }
                }
            }
        }
        Histogram[] hists = {simpleRBE};        
        return hists;

    }
    
        public void writeSimpleRBE(int ke, String filename)
    {
        filename = filename+".csv";
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file " + filename + ". Histogram data was not saved.");
            return;
        }

        // now make a loop to write the contents of each bin to disk, one number at a time
        // together with the x-coordinate of the centre of each bin.
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){ 
            // comma separated values
            outputFile.println(zi + "," + binCentre[zi] + "," + getContent(zi, xi, yi));
        }
        }
        }
        outputFile.close(); // close the output file
        System.out.println(filename+" written!");
    }
    

// CASE 2 - CARABE-FERNANDEZ MODEL
    
    // returns CF RBE max
    public double [] getCarFerRBEmax(double alpha_beta, int ke, int nbin) {
        // get LET information from LET histogram
        double [] LET = new double [nbins];
        LET = getLET(ke, nbin);
        
        // calculate RBEmax for each bin using LET
        double [] CarFerRBEmax = new double [nbins];
        for (int n = 0; n < nbins; n++)
            CarFerRBEmax[n] = 0.834 + (0.154*(2.686/(alpha_beta))*LET[n]);

        return CarFerRBEmax;
    }
    
    // returns CF RBE min
    public double [] getCarFerRBEmin(double alpha_beta, int ke, int nbin) {
        // get LET information from LET histogram
        double [] LET = new double [nbins];
        LET = getLET(ke, nbin);

        // calculate RBEmax for each bin using LET
        double [] CarFerRBEmin = new double [nbins];
        for (int n = 0; n < nbins; n++)
            CarFerRBEmin[n] = 1.09 + 0.006*(2.686/(alpha_beta))*LET[n];

        return CarFerRBEmin;
    }
    
        // returns RBE weighted dose histogram for case 1 max
    public Histogram [] CarFerMinHist(int ke, int nbin, double [] EnergyLossArray) {
        Histogram CarFerMinHist = new Histogram(DVHnbins, 0, 150, "Carabe-Fernandez Min");
        
        double totalEloss = finddE(EnergyLossArray);
        double [] CarFerRBE = getCarFerRBEmin(alpha_beta, ke, nbins);
        
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    // if particle is in tumour volume, fill LET hist
                    if(isInSphere(xi,yi,zi)){
                        // weight dose based on RBE of z bin and fill array
                        double CarFerRBEweight = CarFerRBE[zi]*totalEloss;
                        CarFerMinHist.fill(CarFerRBEweight);
                    }
                }
            }
        }
        Histogram[] hists = {CarFerMinHist};        
        return hists;

    }
    
        // returns RBE weighted dose histogram for case 1 min
        public Histogram [] CarFerMaxHist(int ke, int nbin, double [] EnergyLossArray) {
        Histogram CarFerMaxHist = new Histogram(DVHnbins, 0, 150, "Carabe-Fernandez Max");
        
        double totalEloss = finddE(EnergyLossArray);
        double [] CarFerRBE = getCarFerRBEmax(alpha_beta, ke, nbins);
        
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    // if particle is in tumour volume, fill LET hist
                    if(isInSphere(xi,yi,zi)){
                        // weight dose based on RBE of z bin and fill array
                        double CarFerRBEweight = CarFerRBE[zi]*totalEloss;
                        CarFerMaxHist.fill(CarFerRBEweight);
                    }
                }
            }
        }
        Histogram[] hists = {CarFerMaxHist};        
        return hists;

    }
        

    
        
    // CASE 3 - WEDENBERG MODEL
        
    public double [] getWedenRBEmax(double alpha_beta, int ke, int nbin) {
        // get LET information from LET histogram
        double [] LET = new double [nbins];
        LET = getLET(ke, nbin);
        
        // calculate RBEmax for each bin using LET
        double [] WedenRBEmax = new double [nbins];
        for (int n = 0; n < nbins; n++)
            WedenRBEmax[n] = 1.00 + (0.434/(alpha_beta))*LET[n];

        return WedenRBEmax;
    }
    
    // returns CF RBE min
    public double getWedenRBEmin(double totalEloss) {
                double RBE = 1.0;
                double WedenRBEmin = RBE*totalEloss;      
        return WedenRBEmin;
    }
        
    // returns RBE weighted dose histogram for CASE 3 MIN
    public Histogram [] WedenMinHist(int ke, int nbin, double [] EnergyLossArray) {
        Histogram WedenMinHist = new Histogram(DVHnbins, 0, 150, "Wedenberg Min");
        
        double totalEloss = finddE(EnergyLossArray);
        // returns RBE weighted energy
        double WedenRBEmin = getWedenRBEmin(totalEloss);
        
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    // if particle is in tumour volume, fill LET hist
                    if(isInSphere(xi,yi,zi)){
                        WedenMinHist.fill(WedenRBEmin);
                    }
                }
            }
        }
        Histogram[] hists = {WedenMinHist};        
        return hists;

    }
    
        // returns RBE weighted dose histogram for CASE 3 MAX
        public Histogram [] WedenMaxHist(int ke, int nbin, double [] EnergyLossArray) {
        Histogram WedenMaxHist = new Histogram(DVHnbins, 0, 150, "Wedenberg Max");
        
        double totalEloss = finddE(EnergyLossArray);
        double [] WedenRBE = getWedenRBEmax(alpha_beta, ke, nbins);
        
        for(int zi = 0; zi<nbins; zi++){
            for(int xi = 0; xi<nbins;xi++){
                for(int yi = 0; yi<nbins;yi++){
                    // if particle is in tumour volume, fill LET hist
                    if(isInSphere(xi,yi,zi)){
                        // weight dose based on RBE of z bin and fill array
                        double WedenRBEweight = WedenRBE[zi]*totalEloss;
                        WedenMaxHist.fill(WedenRBEweight);
                    }
                }
            }
        }
        Histogram[] hists = {WedenMaxHist};        
        return hists;

    }   
        
        
    

        
    
        
}
        
    

