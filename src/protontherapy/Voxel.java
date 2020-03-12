package protontherapy;

import java.io.*;

// class to create voxels
class Voxel
{    
    // defining variables
    int nvoxel = 0;
    int i = 0;
    
    // tumour dimensions - hard-coded in atm, need to relate to proton therapy class
    double x_tumour = 0.4; //m
    double y_tumour = 0.4; //m
    double z_tumour = 0.3; //m
    
    // total number of voxels desired
    int total_voxel = 10;
    
    // calculating voxel dimensions
    double voxel_x = x_tumour/total_voxel;        
    double voxel_y = y_tumour/total_voxel;
    double voxel_z = z_tumour/total_voxel;
            
    // initializes empty arrays for x,y,z coords
    double [] x_coords = new double[total_voxel];
    double [] y_coords = new double[total_voxel];
    double [] z_coords = new double[total_voxel];

    // x array - *** what is going wrong here?
    for (int i = 0; i <= total_voxel ; i++) {
        x_coords[i] = i*voxel_x;
    }
    
    // y array
    for (int i = 0; i <= total_voxel ; i++) {
        y_coords[i] = i*voxel_y;
    }
    
    // z array
    for (int i = 0; i <= total_voxel ; i++) {
        z_coords[i] = i*voxel_z;           
    } 
   /*
    public [] AddVoxel(double x_coords[i], double y_coords[i], double z_coords[i],
                         double x_coords[i+1], double y_coords[i+1], double z_coords[i+1]);
    {
        if (nvoxel >= total_voxel) {
            return -1;
        }
        
        type[nvoxel] = 1;
        shapes[nshapes] = new double[6];
        shapes[nshapes][0] = x0;
        shapes[nshapes][1] = y0;
        shapes[nshapes][2] = z0;
        shapes[nshapes][3] = x1;
        shapes[nshapes][4] = y1;
        shapes[nshapes][5] = z1;

      

        nshapes++;
        return (nshapes-1);
    }*/
}