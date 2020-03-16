package protontherapy;

import java.io.*;

/// this import is needed for the file input/output functionality
import java.io.*;

class Voxel
{       
    private double binlow_z, binhigh_z;
    private double binwidth;
    private int nbins;
    private double[] binCentre;
    private String histname;

    // double array to store the actual histogram data
    private double[] sumBinEnergy;

    private long underflow, overflow;
    private long nfilled;

    // constructor for the class Histogram
    public Voxel(int numberOfBins, double start, double end, String name)
    {
        // store the parameters and setup the histogram
        // note that parameters need to have different names than class variables
        nbins = numberOfBins;
        
        // need bins in 3 dimensions
        binlow_z = start; // beginning voxel coordinate
        binhigh_z = end; // end voxel coordinate
        
        
        histname = name;
        
        // variables 
        binwidth = (binhigh_z - binlow_z) / (double) nbins;
        sumBinEnergy = new double[nbins];
        underflow = 0;
        overflow = 0;
        nfilled = 0;
                
        // calculate and save the z coordinate of the centre of each bin
        binCentre = new double[nbins];
        for (int i = 0; i < nbins; i++) {
            binCentre[i] = binlow_z + (i+0.5)*binwidth;
        }
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
    
    // fills voxels with energy
    public void fill_z(double z_energy, Particle p) 
    {
        // increases underflow count
        if (p.z < binlow_z) {
            underflow++;
        // increases overflow
        } else if (p.z >= binhigh_z) {
            overflow++;
        // filling bins 
        } else {
            
            int ibin = (int) ((p.z - binlow_z)/binwidth);
//            System.out.println("Cheetah");
//            System.out.println(p.z);
//            System.out.println(binlow_z);
//            System.out.println(binwidth);
            sumBinEnergy[ibin] = sumBinEnergy[ibin] + z_energy;
            //sumWeights[ibin] = sumWeights[ibin] + 1.0;
        }
        // increases filled bin count by one
        nfilled++;
    }

    // 
    public double getBinEnergy(int nbin)
    {
        // returns the contents on bin 'nbin' to the user
        return sumBinEnergy[nbin];
    }

    public double getError(int nbin)
    {
        // returns the error on bin 'nbin' to the user
        return Math.sqrt(sumBinEnergy[nbin]);
    }
    
    //-------------------------------------
    public void print()
    {
        for (int bin = 0; bin < getNbins(); bin++) {
            System.out.println("Bin " + bin + " = " +getBinEnergy(bin)
                               + " +- " + getError(bin));
        }
        System.out.println("The number of fills = " + getNfilled());
        System.out.println("Underflow = " + getUnderflow()
                           + ", Overflow = " + getOverflow());
    }
    
    //-------------------------------------
    public void writeToDisk(String filename)
    {
        // this sends the output to a file with name "filename"
        // the block with try { ... } catch (IOException e) { ... } is needed to handle the case,
        // where opening the file fails, e.g. disk is full or similar
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file " + filename + ". Histogram data was not saved.");
            return;
        }

        // Write the file as a comma seperated file (.csv) so it can be read it into EXCEL
        // first some general information about the histogram
        outputFile.println("histname, " + histname);
        outputFile.println("binlow, " + binlow_z);
        outputFile.println("binwidth, " + binwidth);
        outputFile.println("nbins, " + nbins);
        outputFile.println("underflow, " + underflow);
        outputFile.println("overflow, " + overflow);

        // now make a loop to write the contents of each bin to disk, one number at a time
        // together with the x-coordinate of the centre of each bin.
        for (int n = 0; n < nbins; n++) {
            // comma separated values
            outputFile.println(n + "," + binCentre[n] + "," + getBinEnergy(n) + "," + getError(n));
        }
        outputFile.close(); // close the output file
        System.out.println("File written!");
    }
}
