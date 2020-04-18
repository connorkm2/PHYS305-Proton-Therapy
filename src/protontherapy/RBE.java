package protontherapy;

import java.io.*;
import java.util.Arrays;
import java.util.stream.*;

public class RBE {
    
public double [] getPhantomStart(double [] phantomPosition) {
       double [] start = Arrays.copyOfRange(phantomPosition, 0, 2); 
       return start;
}

public double [] getPhantomEnd(double [] phantomPosition) {
        double [] end = Arrays.copyOfRange(phantomPosition, 3, 5);
       return end;
}
    
public double [][][] LET(double x, double y, double z, int nbins, int numberOfEvents, double [] phantomPosition) {
    // intialising array to store LET of each voxel
    double [][][] LET = new double [nbins][nbins][nbins];
    
    double [] start = getPhantomStart(phantomPosition);
    double [] end = getPhantomEnd(phantomPosition);
    
    // initialise voxel coordinates
    Voxel voxel_x = new Voxel(nbins, start, end, "RBE", false);
    voxel_x.getVoxelx(phantomPosition, nbins);
    
    Voxel voxel_y = new Voxel(nbins, start, end, "RBE", false);
    voxel_y.getVoxely(phantomPosition, nbins);
    
    Voxel voxel_z = new Voxel(nbins, start, end, "RBE", false);
    voxel_z.getVoxelz(phantomPosition, nbins);
        
    // for all particles
    for (int a = 0; a < numberOfEvents; a++) {
        
        // for voxels in z
        for (int i = 0; i < nbins; i++){
            // for voxels in x
            for (int j = 0; j < nbins; j++) {
                // for voxels in y
                for (int k = 0; k < nbins; k++) {
                    
                    // checking if particle is in certain voxel volume
                    if (x > voxel_x[i]) && (x < voxel_x[i+1]) && 
                        (y > voxel_y[j]) && (y < voxel_y[j+1]) &&
                        (z > voxel_z[k]) && (z < voxel_z[k+1]) { 
                            
                        
                    }
                

    
    /*
    for all voxels
        for all particles
            if particle is in particular voxel
                calculate dE
                store in array
    
                calculate dz
                store in array
    
    
        calculate LET(v) for each voxel
        store in array
    
    */
    
    
                            }    
}
}
        return LET;
}
}




    


