package protontherapy;


class MCS
{
    private double Z;
    private double rho;
    private double A;
    private double z;
    
    public MCS(double rho_in, double Z_in, double A_in)
    {
        Z = Z_in;
        A = A_in;
        rho = rho_in;
        
        //for a muon
        z = -1;
    }

    public double getX0()
    {
        
        double sqrtZ = Math.sqrt(Z);
        double top = 716.4*A;
        double bot = rho*Z*(Z+1)*Math.log(287/(sqrtZ));
        double X0 = top/bot ;   // shall return X0 in m
        
//        System.out.println("Rabbit");
//        System.out.println(Z);
//        System.out.println(bot);
        
        return X0/100;
    }

    public double getTheta0(Particle part, double x)
    {
        double X0 = getX0();
        double var1 = z*Math.sqrt(x/X0)*(1+(0.038*Math.log(x/X0)));
        double Theta0 = -(13.6/(part.beta()*part.momentum()))*var1;// shall return Theta0 for material thickness x
//        
//        System.out.println("fish");
//        System.out.println(X0);
//        System.out.println(var1);
//        System.out.println(Theta0);

        return Theta0;
    }
}