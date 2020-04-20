package protontherapy;

import java.io.*;
import java.util.Arrays;
import java.util.stream.*;

public class RBE {
    
public double [][][] LET(double x, double y, double z, int nbins, int numberOfEvents, double [] phantomPosition,
                            double [] voxelx, double [] voxely, double [] voxelz, double [] EnergyLossArray, double [] stepsize) {
    
    // density
    double rho = 1000; //kgm^-3

    // intialising array to store LET of each voxel
    double [][][] LET = new double [nbins][nbins][nbins];
    
    // initialisng arrays to store dE and dz for each voxel.
    double [][][] dE = new double [nbins][nbins][nbins];
    double [][][] dEdz = new double [nbins][nbins][nbins];
    
    double [] start = getPhantomStart(phantomPosition);
    double [] end = getPhantomEnd(phantomPosition);
    
    voxelx = getVoxelx(phantomPosition, nbins);
    voxely = getVoxely(phantomPosition, nbins);
    voxelz = getVoxelz(phantomPosition, nbins);
        
    // for voxels in z
    for (int i = 0; i < nbins; i++){
        // for voxels in x
        for (int j = 0; j < nbins; j++) {
            // for voxels in y
            for (int k = 0; k < nbins; k++) {
                
                    // for all particles
                    for (int a = 0; a < numberOfEvents; a++) {
                    
                    // checking if particle is in certain voxel volume
                    if ((voxelx[i] <= x) && (x <= voxelx[i+1]) && 
                        (voxely[j] <= y) && (y <= voxely[j+1]) &&
                        (voxelz[k] <=  z) && (z <= voxelz[k+1])) { 
                        
                        // calculate energy loss in specific voxel and sum to find total energy gained by voxel
                           double dEProton = finddE(EnergyLossArray);
                           dE[i][j][k] += dEProton;
                        
                        // calculate track length of particular particle
                           double dz = findTrackLength(stepsize);
                           
                        // calculate dE/dz for each voxel and sum 
                           dEdz[i][j][k] += (dEProton/dz);
                           
                    }
    
                }  
                // calculate final LET for each voxel
                LET[i][j][k] = (dE[i][j][k]*dEdz[i][j][k]*(1/rho))/dE[i][j][k];
}
}
}

        return LET;
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
    
    // CASE 2 - CARABE-FERNANDEZ MODEL
        public double [][][] CarFerRBE(double [][][] dE, int nbins) {
        double alpha_beta = 2.686;
            
        // initialising RBE weighted dose
        double [][][] CarFerRBEWeight = new double [nbins][nbins][nbins];
        
        // initialising array for variable RBE
        double [][][] RBEmax = new double [nbins][nbins][nbins];
        double [][][] RBEmin = new double [nbins][nbins][nbins];
        
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                RBEmax[i][j][k] = 0.834 + 0.154*(2.686/(alpha_beta));
                    
                CarFerRBEWeight[i][j][k] = RBE*dE[i][j][k];
                }
            }
        }
        
        return CarFerRBEWeight;
        }


    
}

