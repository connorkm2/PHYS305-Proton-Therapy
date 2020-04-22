package protontherapy;

import static java.lang.Math.sin;
import java.util.Arrays;
import java.util.Random;

class ProtonTherapy extends Parameters
{
    
    /*  PARAMTERS  */
    static final double tumourDepth = 6; // cm

    // parameters used for the ParticleTracker
    // total time to track (seconds), number of time steps, use/don't use RK4
    static final double time = 1E-7;
    static final int nsteps = 100000;
    static final boolean useRungeKutta4 = true;

    // minimum size of experimental features,
    // ensure this is a factor ~10 smaller than the thinest elements of the experiment
    static final double minfeaturesize = 0.0005;

    // start values
    //static final double startMomentum = 150.; // MeV
    static final double startKineticEnergy = 245.; // MeV
    static final double startAngle = 0;      // Radians
    
    // Number of events to simulate (ie. the number of particles)
    static int numberOfEvents = 100000;
    
    // Energy ranges for when we add multiple energy ranges to flatten dose area
    static final double [][] energies = getEnergiesNEW(250);
    
    static final double delta = 0.2; // for det hists
    static final double delta_A = 0.05; //For first 4 hists gen and sim theta
    
    static final int Nb = 400; // for det hists
    static final int Nb_A = 200; //For first 4 hists gen and sim theta
    
    static Random randGen = new Random();
    
    public static void main (String [] args )
    {
        System.out.println(Arrays.deepToString(getEnergiesNEW(200)));
        // setup histograms for analysis
        Histogram hist_gen_mom = new Histogram(50, 0., 150*1.01, "initial generated Momentum");
        Histogram hist_sim_mom = new Histogram(50, 0., 150*1.01, "simulated final Momentum");

        Histogram hist_gen_theta_zx = new Histogram(Nb_A, startAngle-delta_A, startAngle+delta_A, "initial generated Theta z-x");
        Histogram hist_gen_theta_zy = new Histogram(Nb_A, -delta_A, delta_A, "initial generated Theta z-y");
        Histogram hist_sim_theta_zx = new Histogram(Nb_A, startAngle-delta_A, startAngle+delta_A, "simulated Theta z-x");
        Histogram hist_sim_theta_zy = new Histogram(Nb_A, -delta_A, delta_A, "simulated Theta z-y");
        
        // initialising detector histograms
        Histogram hist_det_theta_zx = new Histogram(Nb, -0.2, 0.2, "measured theta z-x");
        Histogram hist_det_theta_zy = new Histogram(Nb, -0.2, 0.2, "detector theta z-y");
        
        
        // Define the genotrical properties of the experiment in method SetupExperiment()
        Geometry Experiment = SetupExperiment();
                
        for(int ke = 0; ke < energies.length; ke++){
            // start of main loop: run the simulation numberOfEvents times
            
//            // Added for weighting the subsequent beams
            int np = (int) Math.rint( numberOfEvents*energies[ke][1]);
//           System.out.println(numberOfEvents);
//            System.out.println("Rabbit");
            System.out.println(np);
//            System.out.println(energies[ke][0]);
//            System.out.println(energies[ke][1]);
            
            for (int nev = 0; nev < np; nev++) {

                if (nev % 1000 == 0) {
                    //System.out.println("Simulating event " + nev);
                }

                // get the particles of the event to simulate
                Particle [] Particles_gen = GetParticles(energies[ke][0]);

                // simulate propagation of each generated particle,
                // store output in Particles_sim and Tracks_sim
                Particle [] Particles_sim = new Particle[Particles_gen.length];
                Track [] Tracks_sim = new Track[Particles_gen.length];

                for (int ip = 0; ip < Particles_gen.length; ip++) {
                     //some output (need to disable before running large numbers of events!)
    //                 System.out.println("Simulating particle " + ip + " of event " + nev);
    //                 Particles_gen[ip].print();

                    ParticleTracker tracker = new ParticleTracker(Particles_gen[ip], time, nsteps, useRungeKutta4);

                    Particles_sim[ip] = tracker.track(Experiment, ke);

//                     System.out.println("Output particle");
//                     Particles_sim[ip].print();

                    // save the full simulated track for later analysis
                    Tracks_sim[ip] = tracker.getTrack();


                    // write scatter plot for event 0, particle 0 to disk into file "output_particle.csv"
                    if (nev == 0 && ip == 0) {
                        //Tracks_sim[ip].writeToDisk("output_particle.csv");
                    }
                }
                // end of simulated particle propagation
                
                // simulate detection of each particle in each element from the simulated tracks
                // this is just for dumping the simulation to the screen
    //             for (int ip = 0; ip < Tracks_sim.length; ip++) {
    //              double [][] detection_txyz = Experiment.detectParticles(Tracks_sim[ip]);
    //
    //                 for (int idet = 1; idet < Experiment.getNshapes(); idet++) {
    //                     System.out.println("Particle " + ip + " detection in volume " + idet);
    //                     System.out.println("(t,x,y,z) = (" + detection_txyz[idet][0] + ", "
    //                                     + detection_txyz[idet][1] + ", "
    //                                     + detection_txyz[idet][2] + ", "
    //                                     + detection_txyz[idet][3] + ")");
    //                 }
    //             }

                // at this stage the simulation is done and we analyse the output
                // typically we don't want to store thousands of tracks,
                // but rather calculate some interesting quantities and make histograms of the distributions

                // the following analysis is specific to single-particle events with two "detectors"
                // it would look different in more complex cases

                // retrieve initial generated particle momentum and fill histogram
                hist_gen_mom.fill(Particles_gen[0].momentum());
                // retrieve simulated particle momentum at the end of the simulation and fill histogram
                hist_sim_mom.fill(Particles_sim[0].momentum());

                // calculate theta angles in the z-x and z-y planes
                // theta ~ atan2(x, z)

                // generated - this should be equal to given startAngle and zero
                double gen_theta_zx = Math.atan2(Particles_gen[0].px, Particles_gen[0].pz);
                hist_gen_theta_zx.fill(gen_theta_zx);
                double gen_theta_zy = Math.atan2(Particles_gen[0].py, Particles_gen[0].pz);
                hist_gen_theta_zy.fill(gen_theta_zy);

                // same after simulation - muon will have scattered around a bit
                double sim_theta_zx = Math.atan2(Particles_sim[0].px, Particles_sim[0].pz);
                hist_sim_theta_zx.fill(sim_theta_zx);
    //            System.out.println("cat");
    //            //stem.out.println(Particles_sim[0].px);
    //            System.out.println(Particles_sim[0].pz);
                //Particles_sim[0].print();

                double sim_theta_zy = Math.atan2(Particles_sim[0].py, Particles_sim[0].pz);
                hist_sim_theta_zy.fill(sim_theta_zy);

                // after detection: reconstruct the angle from the two detected positions!
                // the detectors have volume number 2+3 (see printout)
                double [][] detection_txyz = Experiment.detectParticles(Tracks_sim[0]);
                double x_det = detection_txyz[3][1]; // x-coo in detector 2
                double y_det = detection_txyz[3][2]; // y-coo in detector 2
                double z_det = detection_txyz[3][3]; // z-coo in detector 2

                // calculating the detector theta angles
                double det_theta_zx = Math.atan2(x_det, z_det);
                double det_theta_zy = Math.atan2(y_det, z_det);

                // filling histograms
                hist_det_theta_zx.fill(det_theta_zx);
                hist_det_theta_zy.fill(det_theta_zy);

                // end of analysis


            }
        }
        // end of main event loop

        // write out histograms for plotting and futher analysis
//        hist_gen_mom.writeToDisk("gen_mom.csv");
//        hist_sim_mom.writeToDisk("sim_mom.csv");
//        
//        hist_gen_theta_zx.writeToDisk("gen_theta_zx.csv");
//        hist_gen_theta_zy.writeToDisk("gen_theta_zy.csv");
//        hist_sim_theta_zx.writeToDisk("sim_theta_zx.csv");
//        hist_sim_theta_zy.writeToDisk("sim_theta_zy.csv");
        
        // writing to disk for no smearing
        hist_det_theta_zx.writeToDisk("det_theta_zx.csv");
        hist_det_theta_zy.writeToDisk("det_theta_zy.csv");
        
        Experiment.writeEnergyHist(0.3,"energy_hist.csv");

    }

    public static Geometry SetupExperiment ()
    {
        // example setup the experiment

        Geometry Experiment = new Geometry(minfeaturesize);
        
        // this line defines the size of the experiment in vacuum
        double[] pos1 =      {-0.5, -0.5, 0.,                // start x, y, z
                             0.5, 0.5, 1.5};  // end   x, y, z
        Experiment.AddCuboid(pos1,
                             0., 0., 0., "Vacuum");                     // zeros for "vacuum"

        // Block of tantalum of thickness 1cm
        Experiment.AddCuboid(scatererPosition,   // end   x, y, z
                             16.65, 73, 180.94788, "Tantalum scatterer");           // density, Z, A
        
        // water phantom
        Experiment.AddCuboid(phantomPosition,   // end   x, y, z
                             1, 7.42, 18.015, "Water phantom");           // density, Z, A
        
        // two 1mm-thin "silicon detectors" 10cm and 20cm after the iron block
        double[] pos2 =      {-0.5, -0.5, 0.219, // start x, y, z
                             0.5, 0.5, 0.22};   // end   x, y, z
        Experiment.AddCuboid(pos2,
                             2.33, 14, 28.085, "Si detector");                 // density, Z, A
        
         //Contoured Scatterer
        Experiment.AddContour(-0.2, -0.2, 0.05,
                             16.65, 73, 180.94788, "Contour scatter");
        
//        //Aperture
//        Experiment.AddAperture(0.03, 0.2, 0.2, 0.21,
//                               11.34, 82, 207.2, "Aperture");
       
        
        Experiment.Print();

        return Experiment;
    }

    
    public static Particle[] GetParticles(double energy)
    {
        //System.out.println("p");
        // example to simulate just one proton starting at (0,0,0)
        // with a total momentum startMomentum and theta=startAngle
        // we follow the particle physics "convention"
        // to have the z-axis in the (approximate) direction of the beam
        // this just sets up the array (for a case where one event has more than one particle)
        Particle [] Particles_gen = new Particle[1];
        
        // converting input kinetic energy to momentum
        double startMomentum = Math.sqrt((938+energy)*(938+energy)-(938*938));
//        System.out.println(startMomentum);
        
        // create particle and set properties
        Particles_gen[0] = new Particle();

        // initial momentum px,py,pz (MeV)
        double phi = 0;
        Particles_gen[0].px = startMomentum*Math.sin(startAngle)*Math.cos(phi);
        Particles_gen[0].py = startMomentum*Math.sin(startAngle)*Math.sin(phi);
        Particles_gen[0].pz = startMomentum*Math.cos(startAngle);

        // Set charge and mass of a positive proton
        Particles_gen[0].m = 938;
        Particles_gen[0].Q = +1;  

        double randValue = (-0.05)*(0.05+0.05) * randGen.nextGaussian();

        // initial position (x,y,z) = (0,0,0)
        Particles_gen[0].x = randValue;
        Particles_gen[0].y = randValue;
        Particles_gen[0].z = 0;//randValue;

        return Particles_gen;
    }
    
    public static double[][] getEnergies(double start, double end){
        int steps = (int) (end-start) +1;
        double [][] values = new double[steps][2];
        double fraction = 1;
        
        for(int i = 0; i < steps; i++){
            values[i][0] = end-i;
            
            // The fraction by which to reduce the nev
            values[i][1] = fraction;
            fraction = fraction - ((double)1/(10*steps));
        }
        
        return values;
    }
    
//    This new function will determine the range of energies based on a fixed 
//    distance between pristine peaks as described in the refernece material
//    it uses the relationship (R ~ alpha*E^1.8) which is taken from literature.
//    alpha is taken as 0.0022 and R is given in cm.
//    The range in energy generated here is small and may need some testing,
//    but it is based on the paramters from the book.
    
    public static double[][] getEnergiesNEW(double startE){
        int numPeaks = 16;
        
        double R_0 = (0.0022*Math.pow(startE, 1.8)); // cm
        // coluumn 1 stores energies, column 2 stores weights 
        double [][] energies = new double[numPeaks][2];
        double steps = R_0/numPeaks;
        
//        energies[0][0] = startE;
//        energies[0][1] = 1 - Math.pow(1 - 1/(2*(double)numPeaks), (1 - 1/p));
        
        for(int i = 0; i < numPeaks; i++){
            //double nextEnergy = Math.pow(((R_0-(0.6*i))/0.0022),1/1.8);
            //System.out.println(i);
            if(i==0){
            energies[i][0] = startE - i*steps;
            energies[i][1] = 1;    
            }else{
                energies[i][0] = startE - i*steps;
                energies[i][1] = 0.3*Math.pow(i, -0.80) + 0.1;
                //energies[i][1] = 0.4*Math.pow(i, -0.50) ;//+ 0.1;
            }
//            System.out.println(energies[i][1]);
//            System.out.println(energies[i][0]);
//            System.out.println("dog");
        }
        
        double tumourDepth = R_0 - (0.0022*Math.pow(energies[numPeaks-1][0], 1.8));
        System.out.println("Depth of tumour: "+tumourDepth);
        return energies;
        
    }
    
        public static double [][] getEnergiesNEW2(double startE){
        double p = 1.04;
        int numPeaks = 10;
        double R_0 = 0.0022*Math.pow(startE, p);
        double targetWidth = R_0*0.15;
//        System.out.println(R_0);
//        System.out.println(targetWidth);
//        System.out.println("dog");


        double [][] energies = new double[numPeaks+1][2];

        for(int k = 0; k < numPeaks;k++){
            //if(k == 6){ p = 1.4; }
            double range  = (1-((1-((double)k/(double)numPeaks))*targetWidth))*R_0;
            energies[k][0] = Math.pow((range/0.0022), 1/p);
            //energies[k][0] = Math.pow(((R_0-(0.06*(double)(numPeaks-k)))/0.0022),1/p);
            //System.out.println(energies[k][0]);

            if(k == 0){
                energies[k][1] = Math.pow((1-(1/(2*(double)numPeaks))),1-(1/p));
//                energies[k][1] = (Math.pow(1-((1/(double)numPeaks)*((double)k - 0.5)),1-(1/p))
//                        - Math.pow(1-((1/(double)numPeaks)*((double)k + 0.5)),1-(1/p)))*10;
            }else if(k == numPeaks){
                energies[k][1] = Math.pow((1/(2*(double)numPeaks)), 1-(1/p))*10;
            }else{
                energies[k][1] = (Math.pow(1-((1/(double)numPeaks)*((double)k - 0.5)),1-(1/p))
                        - Math.pow(1-((1/(double)numPeaks)*((double)k + 0.5)),1-(1/p)))*10;
            }
            //System.out.println(k);
        }

        energies[numPeaks][0] = startE;
        energies[numPeaks][1] = 0.85;

        return energies;

    }

}