package protontherapy;
// this import is needed for the file input/output functionality
import java.io.*;
import java.util.stream.*;

class Track
{
    // store momentum components and positions of the Particle
    // mass is a fixed value
    
    float mass;

    float [] px;
    float [] py;
    float [] pz;

    float [] t;
    float [] x;
    float [] y;
    float [] z;

    int lastpos;

    Track(double m, int steps) {
        mass = (float) m;
        px = new float[steps];
        py = new float[steps];
        pz = new float[steps];
        t = new float[steps];
        x = new float[steps];
        y = new float[steps];
        z = new float[steps];
        lastpos = 0;
    }

    void savePositionMomentum(Particle p) {
        if (lastpos >= px.length) {
            System.err.println("Out of space in output track storage, allow for more steps.");
            return;
        }
        px[lastpos] = (float) p.px;
        py[lastpos] = (float) p.py;
        pz[lastpos] = (float) p.pz;
        x[lastpos] = (float) p.x;
        y[lastpos] = (float) p.y;
        z[lastpos] = (float) p.z;
        t[lastpos] = (float) p.t;
        lastpos++;
    }
        
    public void writeToDisk(String filename)
    {
        // this sends the output to a file with name "filename"
        // the block with try { ... } catch (IOException e) { ... } is needed to handle the case,
        // where opening the file fails, e.g. disk is full or similar
        PrintWriter outputFile;
        try {
            outputFile = new PrintWriter(filename);
        } catch (IOException e) {
            System.err.println("Failed to open file " + filename + ". Track data was not saved.");
            return;
        }

        // now make a loop to write the contents of the arrays to disk, one number at a time
        outputFile.println("time [s], x [m], y [m], z[m], E [MeV], px [MeV], py [MeV], pz [MeV]");
        for (int n = 0; n < lastpos; n++) {
            // comma separated values
            double E = Math.sqrt(mass*mass + px[n]*px[n] + py[n]*py[n] + pz[n]*pz[n]);
            outputFile.println(t[n] + "," + x[n] + "," + y[n] + "," + z[n] + ","
                               + E + "," + px[n] + "," + py[n] + "," + pz[n]);
        }
        outputFile.close(); // close the output file
        return;
    }
    
    // returns array of step sizes
    public double [] findStepSize(int i, int steps) {
        // initialising array to store distance of each step
       double [] trackLength = new double[steps];
       
       // calculates step distance in x, y, z
       double stepDist = Math.pow(x[i]-x[i-1], 2);
       stepDist = stepDist + Math.pow(y[i] - y[i-1], 2);
       stepDist = stepDist + Math.pow(z[i] - z[i-1], 2);
       
       trackLength[i] = stepDist;
       
       // need to sum total distance, where does this happen? - do I only need to do it in Z?
       
    return trackLength;
    }
    
    /*
    public double [] getEnergyArray(int steps, output) {
            
    // calculate total energy loss
    double [] EnergyLossArray = new double [steps];
    EnergyLossArray[n] = Experiment.doEloss(output, output.distance(lastStep), ke);
    
    return EnergyLossArray;
    }
    */
    

}