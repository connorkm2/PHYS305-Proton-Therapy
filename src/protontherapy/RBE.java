package protontherapy;

import java.io.*;
import java.util.Arrays;
import java.util.stream.*;

public class RBE {
    
public double [][][] LET(double x, double y, double z, int nbins, int numberOfEvents, double [] phantomPosition,
                            double [] voxelx, double [] voxely, double [] voxelz, double [] EnergyLossArray, double [] stepsize) {
    
    // intialising array to store LET of each voxel
    double [][][] LET = new double [nbins][nbins][nbins];
    
    // initialisng arrays to store dE and dz for each voxel.
//    double [][][] dE = new double [nbins][nbins][nbins];
//    double [][][] dz = new double [nbins][nbins][nbins];
    
    double [] start = getPhantomStart(phantomPosition);
    double [] end = getPhantomEnd(phantomPosition);
    
    voxelx = getVoxelx(phantomPosition, nbins);
    voxely = getVoxely(phantomPosition, nbins);
    voxelz = getVoxelz(phantomPosition, nbins);
        
    // for all particles
    for (int a = 0; a < numberOfEvents; a++) {
        
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                    
                    // checking if particle is in certain voxel volume
                    if (voxelx[i] <= x) && (x <= voxelx[i+1]) && 
                        (y <= voxely[j]) && (y <= voxely[j+1]) &&
                        (z <= voxelz[k]) && (z <= voxelz[k+1]) { 
                        
                        // calculate energy loss
//                        double dE[i][j][k] = finddE(EnergyLossArray);
                        
                        // calculate track length
//                        double dz[i][j][k] = findTrackLength(stepsize);
                    }
    
    }    
}
}
        return LET;
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


    // initialise voxel coordinates
//    Voxel voxel_x = new Voxel(nbins, start, end, "RBE", false);
//    voxel_x.getVoxelx(phantomPosition, nbins);
//    
//    Voxel voxel_y = new Voxel(nbins, start, end, "RBE", false);
//    voxel_y.getVoxely(phantomPosition, nbins);
//    
//    Voxel voxel_z = new Voxel(nbins, start, end, "RBE", false);
//    voxel_z.getVoxelz(phantomPosition, nbins);

    


