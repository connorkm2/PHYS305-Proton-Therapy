package protontherapy;

import java.util.Random;

class ProtonTherapy
{
    // Program to run a simulation of particles in a "experiment"

    // The program makes use of the class Particle (almost identical to the one used in week 4)
    // to store the components of the position fourvector (time, x, y, z)
    // and the four-momentum (technically stored as mass and px, py, pz)

    // The particle tracking is performed by the class ParticleTracker,
    // which is almost identical to the class used in week 4, however it is extended
    // to incorporate energy loss, multiple scattering
    // and a simple adaptive algorithm to ensure that single steps end approximately
    // at boundaries betwen different experimental features
    
    // The "experimental geometry" is defined in the class Geometry:
    //    * An experiment is a collection of "volumes", that are numbered with a unique
    //      identifier from 0 to the number of experimental features.
    //    * Currently all experiemntal features are cuboids of different sizes and materials,
    //      with the sides aligned to the coordinate system axes. The example is a block of iron
    //      (user-definable length) + two "planar detectors"
    //    * Internally to the "Geometry" class are two helper classes that implement the
    //      formulas for calculation of energy loss (class "Energy loss")
    //      and multiple scattering angles (class "MCS")
    //    * The main functionality is to check, if a certain particle is inside a certain volume
    //      and apply energy loss and multiple scattering during simulation
    //    * The class also provides a simple mechanism to detect changes in the volume as
    //      the particle propagates and suggest an adapted step length to keep one step within
    //      one volume (the granularity of this scan is adjusted with "minfeaturesize")
    //    * At the end, the class is used to "detect" particle signals in certain volumes (detectors)
    //
    //
    // At the end of the simulation of each event, one may analyse the
    // results and fill histograms. Examples are provided to perform calculations using the:
    //    * Generated particles (Particles_gen)
    //    * Simulated particles (Particles_sim) - these include the effect of energy loss and
    //      multiple scattering
    //    * Detector response (Particles_det) - these provide measurement points that can be
    //      further used to reconstruct particles like in a real experiment, where
    //      Particles_gen and Particles_sim are unknown

    // parameters used for the ParticleTracker
    // total time to track (seconds), number of time steps, use/don't use RK4
    static final double time = 1E-7;
    static final int nsteps = 100000;
    static final boolean useRungeKutta4 = false;

    // minimum size of experimental features,
    // ensure this is a factor ~10 smaller than the thinest elements of the experiment
    static final double minfeaturesize = 0.0005;

    // start values
    //static final double startMomentum = 150.; // MeV
    static final double startKineticEnergy = 245.; // MeV
    static final double startAngle = 0;      // Radians
    
    // Number of events to simulate (ie. the number of particles)
    static final int numberOfEvents = 10000;
    
    // Energy ranges for when we add multiple energy ranges to flatten dose area
    static final double [] energies = getEnergies(230, 250);
    
    static final double delta = 0.2; // for det hists
    static final double delta_A = 0.05; //For first 4 hists gen and sim theta
    
    static final int Nb = 100; // for det hists
    static final int Nb_A = 200; //For first 4 hists gen and sim theta
    
    static Random randGen = new Random();
    
    public static void main (String [] args )
    {
        // setup histograms for analysis
        Histogram hist_gen_mom = new Histogram(50, 0., 150*1.01, "initial generated Momentum");
        Histogram hist_sim_mom = new Histogram(50, 0., 150*1.01, "simulated final Momentum");

        Histogram hist_gen_theta_zx = new Histogram(Nb_A, startAngle-delta_A, startAngle+delta_A, "initial generated Theta z-x");
        Histogram hist_gen_theta_zy = new Histogram(Nb_A, -delta_A, delta_A, "initial generated Theta z-y");
        Histogram hist_sim_theta_zx = new Histogram(Nb_A, startAngle-delta_A, startAngle+delta_A, "simulated Theta z-x");
        Histogram hist_sim_theta_zy = new Histogram(Nb_A, -delta_A, delta_A, "simulated Theta z-y");
        
        // initialising detector histograms
        Histogram hist_det_theta_zx2 = new Histogram(Nb, startAngle-delta, startAngle+delta, "measured theta z-x");
        Histogram hist_det_theta_zy2 = new Histogram(Nb, 0, 3, "detector theta z-y");
        Histogram hist_det_theta_zx3 = new Histogram(Nb, startAngle-delta, startAngle+delta, "measured theta z-x");
        Histogram hist_det_theta_zy3 = new Histogram(Nb, 0, 3, "detector theta z-y");
        
        // initialising smeared histograms
        Histogram hist_det_theta_zx2_smear = new Histogram(Nb, startAngle-delta, startAngle+delta, "measured theta z-x");
        Histogram hist_det_theta_zy2_smear = new Histogram(Nb, 0, 3, "detector theta z-y smear");
        Histogram hist_det_theta_zx3_smear = new Histogram(Nb, startAngle-delta, startAngle+delta, "measured theta z-x");
        Histogram hist_det_theta_zy3_smear = new Histogram(Nb, 0, 3, "detector theta z-y smear");
        
        
        // Define the genotrical properties of the experiment in method SetupExperiment()
        Geometry Experiment = SetupExperiment();
                
        for(int ke = 0; ke < energies.length; ke++){
            // start of main loop: run the simulation numberOfEvents times
            for (int nev = 0; nev < numberOfEvents; nev++) {

                if (nev % 1000 == 0) {
                    //System.out.println("Simulating event " + nev);
                }

                // get the particles of the event to simulate
                Particle [] Particles_gen = GetParticles(energies[ke]);

                // simulate propagation of each generated particle,
                // store output in Particles_sim and Tracks_sim
                Particle [] Particles_sim = new Particle[Particles_gen.length];
                Track [] Tracks_sim = new Track[Particles_gen.length];

                for (int ip = 0; ip < Particles_gen.length; ip++) {
                     //some output (need to disable before running large numbers of events!)
    //                 System.out.println("Simulating particle " + ip + " of event " + nev);
    //                 Particles_gen[ip].print();

                    ParticleTracker tracker = new ParticleTracker(Particles_gen[ip], time, nsteps, useRungeKutta4);

                    Particles_sim[ip] = tracker.track(Experiment);

                    // System.out.println("Output particle");
                    // Particles_sim[ip].print();

                    // save the full simulated track for later analysis
                    Tracks_sim[ip] = tracker.getTrack();


                    // write scatter plot for event 0, particle 0 to disk into file "output_particle.csv"
                    if (nev == 0 && ip == 0) {
                        Tracks_sim[ip].writeToDisk("output_particle.csv");
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
                double x_det2 = detection_txyz[2][1]; // x-coo in detector 2
                double x_det3 = detection_txyz[3][1]; // x-coo in detector 3
                double y_det2 = detection_txyz[2][2]; // y-coo in detector 2
                double y_det3 = detection_txyz[3][2]; // y-coo in detector 3
                double z_det2 = detection_txyz[2][3]; // z-coo in detector 2
                double z_det3 = detection_txyz[3][3]; // z-coo in detector 3

                // calculating the detector theta angles
                double det_theta_zx2 = Math.atan2(x_det2, z_det2);
                double det_theta_zy2 = Math.atan2(x_det2, y_det2);
                double det_theta_zx3 = Math.atan2(x_det3, z_det3);
                double det_theta_zy3 = Math.atan2(x_det3, y_det3);

                // filling histograms
                hist_det_theta_zx2.fill(det_theta_zx2);
                hist_det_theta_zy2.fill(det_theta_zy2);
                hist_det_theta_zx3.fill(det_theta_zx3);
                hist_det_theta_zy3.fill(det_theta_zy3);

                // smearing
                double stdev = 0.005;
                double smearing = randGen.nextGaussian()*stdev;

                // smearing histograms
                hist_det_theta_zx2_smear.fill(det_theta_zx2 + smearing);
                hist_det_theta_zy2_smear.fill(det_theta_zy2 + smearing);
                hist_det_theta_zx3_smear.fill(det_theta_zx3 + smearing);
                hist_det_theta_zy3_smear.fill(det_theta_zy3 + smearing);

                //System.out.println("antelope");
                //System.out.println(smearing);

                // end of analysis


            }
        }
        // end of main event loop

        // write out histograms for plotting and futher analysis
        hist_gen_mom.writeToDisk("gen_mom.csv");
        hist_sim_mom.writeToDisk("sim_mom.csv");
        
        hist_gen_theta_zx.writeToDisk("Agen_theta_zx.csv");
        hist_gen_theta_zy.writeToDisk("Agen_theta_zy.csv");
        hist_sim_theta_zx.writeToDisk("Asim_theta_zx.csv");
        hist_sim_theta_zy.writeToDisk("Asim_theta_zy.csv");
        
        // writing to disk for no smearing
        hist_det_theta_zx2.writeToDisk("Bdet_theta_zx2.csv");
        hist_det_theta_zy2.writeToDisk("Bdet_theta_zy2.csv");
        hist_det_theta_zx3.writeToDisk("Bdet_theta_zx3.csv");
        hist_det_theta_zy3.writeToDisk("Bdet_theta_zy3.csv");
        
        // writing to disk for smearing
        hist_det_theta_zx2_smear.writeToDisk("Cdet_theta_zx2_smear.csv");
        hist_det_theta_zy2_smear.writeToDisk("Cdet_theta_zy2_smear.csv");
        hist_det_theta_zx3_smear.writeToDisk("Cdet_theta_zx3_smear.csv");
        hist_det_theta_zy3_smear.writeToDisk("Cdet_theta_zy3_smear.csv");
        
        Experiment.writeEnergyHist("energy_hist.csv");

    }

    public static Geometry SetupExperiment ()
    {
        // example setup the experiment

        Geometry Experiment = new Geometry(minfeaturesize);
        
        // this line defines the size of the experiment in vacuum
        Experiment.AddCuboid(-0.5, -0.5, 0.,                // start x, y, z

                             0.5, 0.5, 0.84,  // end   x, y, z

                             0., 0., 0.);                     // zeros for "vacuum"

        // Block of tantalum of thickness 1cm
        Experiment.AddCuboid(-0.20, -0.20, 0.2,            // start x, y, z
                             0.20, 0.20, 0.21,   // end   x, y, z
                             16.65, 73, 180.94788);           // density, Z, A
                
        // water phantom
        Experiment.AddCuboid(-0.20, -0.20, 0.21,            // start x, y, z
                             0.20, 0.20, 0.71,   // end   x, y, z
                             1, 7.42, 18.015);           // density, Z, A
        
        // two 1mm-thin "silicon detectors" 10cm and 20cm after the iron block
        Experiment.AddCuboid(-0.5, -0.5, 0.915, // start x, y, z
                             0.5, 0.5, 0.92,   // end   x, y, z
                             2.33, 14, 28.085);                 // density, Z, A
        
        Experiment.AddCuboid(-0.5, -0.5, 0.93, // start x, y, z
                             0.45, 0.45, 0.94,   // end   x, y, z
                             2.33, 14, 28.085);                 // density, Z, A
        
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

        double randValue = (-0.2)*(0.2+0.2) * randGen.nextDouble();

        // initial position (x,y,z) = (0,0,0)
        Particles_gen[0].x = 0;
        Particles_gen[0].y = 0;
        Particles_gen[0].z = randValue;

        return Particles_gen;
    }
    
    public static double[] getEnergies(double start, double end){
        int steps = (int) (end-start);
        double [] values = new double[steps];
        
        for(int i = 0; i < steps; i++){
            values[i] = start+i;
        }
        
        return values;
    }
}