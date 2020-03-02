package protontherapy;

class EnergyLoss
{
    
    private double K;
    private double z;
    private double M;
    private double rho;
    private double Z;
    private double A;
    private double m; // MeV // electron mass
    private double W;
    private double I;
        
    public EnergyLoss(double rho_in, double Z_in, double A_in)
    {
        rho = rho_in; // g/m-3
        Z = Z_in;
        A = A_in;
        m = 0.511;
        I = 0.0000135*Z;
        
        K = 0.307075; // MeV m^2
        
        // For a muon
        z = -1; //incident particle charge
        M = 106; // MeV
    }
    
    public double getEnergyLoss(Particle p)
    {
        W = (2*m*p.beta()*p.beta()*p.gamma()*p.gamma()) / (1 + (2*p.gamma()*(m/M)) + ((m/M)*(m/M)));
//        System.out.println(W);

        
        double var1 = 0.5*Math.log((2*m*p.beta()*p.beta()*p.gamma()*p.gamma()*W)/(I*I))-(p.beta()*p.beta());
        double x = ((2*m*p.beta()*p.beta()*p.gamma()*p.gamma()*W)/(I*I))-(p.beta()*p.beta());
//        System.out.println(Z);
//        System.out.println(x);
//        System.out.println("mouse");
        double energyLoss = (K*z*z*rho*(Z/A)*(1/(p.beta()*p.beta())))*var1;
        
//        System.out.println(var1);
//        System.out.println(energyLoss);
//        System.out.println("hamster");
        
         // shall return energy loss in MeV/m
        return energyLoss*100;
    }
    
}