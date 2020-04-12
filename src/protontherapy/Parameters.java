/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package protontherapy;

/**
 *
 * @author Cara
 */

// class to input SOBP parameters, as well as calculate NTCP and TCP
public class Parameters {
    
    public static double tumourDepthInPatient = 7; // cm // depth to center of tumour
    public static double[] tumourDimensions = {3,3,14}; //w,h,d [cm]
    
    public static final double[] phantomPosition = {-0.2, -0.2, 0.22,    //x0,y0,z0 [m]
                                               0.2, 0.2, 0.72};     //x1,y1,z1 [m]
    public static final double[] scatererPosition = {-0.20, -0.20, 0.01,    //x0,y0,z0 [m]
                                                    0.20, 0.20, 0.02};   //x1,y1,z1 [m]
    
    
//NTCP [Normal Tissue Complications Probability] and TCP [Tumour Control Probability] calculation method

// IS ABSORBED DOSE 2D OR 1D 
/*
public double getTCP(double nbins, double voxelVolume, double [][] absorbed_dose) {
	// voxels in ith dose bin
	double m_i = Math.pow(nbins, 2);
	
	// clonogen (tumour cell) density
	double clonogenDensity = 10000000; //cm^-3
	
	// number of treatment fractions [*** is this the same as number of slices??? ***]
	double n = absorbed_dose.length;
	
	// alpha and alpha-beta factors [taken from paper connor sent over]
	double alpha = 0.26; //Gy^-1
	double alpha_beta = 10; //Gy
	
	// initialising remaining clonogens (in each slice) array
	double [] clonogensLeft = new double [nbins]; 
	
	for (int i = 0; i < absorbed_dose.length; i++) {
		// looping over each value of dose
		clonogensLeft[i] = Math.exp(-n * alpha * absorbed_dose[i] * (1 + (absorbed_dose[i]/alpha_beta))) * 
                        m_i * voxelVolume * clonogenDensity;
								
		double clonogensLeftTotal += clonogensLeft[i];
	}
	System.out.println(clonogensLeft);
	System.out.println(clonogensLeftTotal);
	
	// TCP var
	double TCP = Math.exp(-clonogensLeftTotal);
	System.out.println("The TCP is:" + TCP);
	
	return TCP;
}


public double getNTCP(double n, double [][] absorbed_dose, double m_i, double nbins) {
	// total dose
	for (int i = 0; i < absorbed_dose.length; i++) {
		double [] totalDose = n * absorbed_dose[i];
	}
	
	// fractional volume receiving total dose
	// *** SHOULD THIS BE MULTIPLED BY NBINS??? - nslices not yet initialised!
	double v_i = m_i/(m_i * nslices);
	
	// dose that leads to 50% NTCP 
	// value for lung taken from - https://iopscience.iop.org/article/10.1088/0031-9155/53/3/014/pdf
	double TD_50 = 29.9; //Gy
	
	// other variables 
	double s = 1.25;
	double gamma = 1.00;
	double e = 2.718;
	
	// Compound multiplies the term inside the brackets 
	for (int i = 0; i < absorbed_dose.length; i++) {
		double term_1 *= Math.pow(1 - Math.pow(Math.pow(- e * (e * gamma * (1 - totalDose[i]/TD_50)))), s), v_i);
	}
	
	double NTCP = 100 * Math.pow((1 - term_1), (1/s));
	
	return NTCP;
}*/

    
    public double getTumourZ0(){
        double z0 = ((phantomPosition[2]*100)+tumourDepthInPatient)-(tumourDimensions[2]/2);
        System.out.println("z0 = "+z0/100);
        return z0/100;
    }
    
    public double getTumourZ1(){
        double z1 = ((phantomPosition[2]*100)+tumourDepthInPatient)+(tumourDimensions[2]/2);
        System.out.println("z1 = "+z1/100);
        return z1/100;
    }
    
    public double getTumourX0(){
        double x0 = -(tumourDimensions[0]/2);
        System.out.println("x0 = "+x0/100);
        return x0/100;
    }
    
    public double getTumourX1(){
        double x1 = (tumourDimensions[0]/2);
        System.out.println("x1 = "+x1/100);
        return x1/100;
    }
    
    public double getTumourY0(){
        double y0 = -(tumourDimensions[1]/2);
        System.out.println("y0 = "+y0/100);
        return y0/100;
    }
    
    public double getTumourY1(){
        double y1 = (tumourDimensions[1]/2);
        System.out.println("y1 = "+y1/100);
        return y1/100;
    }
    
    public double getTumourCenterPos(){
        return (tumourDepthInPatient/100)+phantomPosition[2]; // m
    }
    
}
