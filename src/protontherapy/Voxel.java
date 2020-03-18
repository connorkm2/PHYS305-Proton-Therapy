package protontherapy;

import java.io.*;

/// this import is needed for the file input/output functionality
import java.io.*;

class Voxel
{       
    private double [] binlow, binhigh;
    private double [] binwidth;
    private int nbins;
    private double[][] binCentre;
    private String histname;

    // double array to store the actual histogram data
    private double[][] sumBinEnergy;

    private long [] underflow, overflow;
    private long [] nfilled;

    // constructor for the class Histogram
    public Voxel(int numberOfBins, double [] start, double [] end, String name)
    {
        // store the parameters and setup the histogram
        // note that parameters need to have different names than class variables
        nbins = numberOfBins;
        
        // need bins in 3 dimensions
        binlow = start; // beginning voxel coordinate
        binhigh = end; // end voxel coordinate
        
        
        histname = name;
        
        //System.out.println(binwidth[0]);
        
        binwidth = new double[3];
        
        // variables 
        for(int i = 0; i < binlow.length; i++){
            binwidth[i] = (binhigh[i] - binlow[i]) / (double) nbins;
        }
        sumBinEnergy = new double[3][nbins];
        underflow = new long[3];
        overflow = new long[3];
        nfilled = new long[3];
                
        // calculate and save the z coordinate of the centre of each bin
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
    
//    // fills voxels with energy
//    public void fill_z(double energy, Particle p) 
//    {
//        // increases underflow count
//        if (p.z < binlow_z) {
//            underflow++;
//        // increases overflow
//        } else if (p.z >= binhigh_z) {
//            overflow++;
//        // filling bins 
//        } else {
//            
//            int ibin = (int) ((p.z - binlow_z)/binwidth);
////            System.out.println("Cheetah");
////            System.out.println(p.z);
////            System.out.println(binlow_z);
////            System.out.println(binwidth);
//            sumBinEnergy[ibin] = sumBinEnergy[ibin] + energy;
//            //sumWeights[ibin] = sumWeights[ibin] + 1.0;
//        }
//        // increases filled bin count by one
//        nfilled++;
//    }
    
    public void fill(double energy, Particle p){
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
                //            System.out.println("Cheetah");
                //            System.out.println(p.z);
                //            System.out.println(binlow_z);
                //            System.out.println(binwidth);
                sumBinEnergy[i][ibin] = sumBinEnergy[i][ibin] + energy;
                //sumWeights[ibin] = sumWeights[ibin] + 1.0;
            }
                // increases filled bin count by one
                nfilled[i]++;
            }
        
    }

    // 
    public double getBinEnergy(int i, int nbin)
    {
        // returns the contents on bin 'nbin' to the user
        return sumBinEnergy[i][nbin];
    }

    public double getError(int i, int nbin)
    {
        // returns the error on bin 'nbin' to the user
        return Math.sqrt(sumBinEnergy[i][nbin]);
    }
    
    //-------------------------------------
    public void print()
    {
        for(int a = 0; a < 3; a++){
            for (int bin = 0; bin < getNbins(); bin++) {
                System.out.println("Bin " + bin + " = " +getBinEnergy(a, bin)
                                   + " +- " + getError(a, bin));
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
                outputFile.println(n + "," + binCentre[i][n] + "," + getBinEnergy(i, n) + "," + getError(i, n));
            }
            outputFile.close(); // close the output file
            System.out.println("File written!");
        }
    }
}
