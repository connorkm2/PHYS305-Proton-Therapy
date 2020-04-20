package protontherapy;

public class Pair {
    
    private final Particle output;
    private final double [] EnergyLossArray;
    private final double [] stepsize;

    public Pair(Particle output, double[] EnergyLossArray, double [] stepsize) {
        this.output = output;
        this.EnergyLossArray = EnergyLossArray;
        this.stepsize = stepsize;
    }


    public Particle getOutput() {
        return output;
    }

    public double [] getEnergyLossArray() {
        return EnergyLossArray;
    }
    
    public double [] getStepsize() {
        return stepsize;
    }
    

public Pair Result(Particle output, double[] EnergyLossArray, double [] stepsize) {

    return new Pair(output, EnergyLossArray, stepsize);
}



}

