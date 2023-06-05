package STAM;

public class STAMMomentApproximation {


    public double [] partialVector;


    STAMMomentApproximation(int N){
        this.partialVector = new double[(N+1) * 6 + 4];
    }

}
