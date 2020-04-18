package protontherapy;

public class Pair {
    
    private final Particle output;
    private final double [] EnergyLossArray;

    public Pair(Particle output, double[] EnergyLossArray) {
        this.output = output;
        this.EnergyLossArray = EnergyLossArray;
    }

    public Particle getOutput() {
        return output;
    }

    public double [] getEnergyLossArray() {
        return EnergyLossArray;
    }

public Pair Result(Particle output, double[] EnergyLossArray) {

    return new Pair(output, EnergyLossArray);
}

}

