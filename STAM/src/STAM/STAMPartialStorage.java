package STAM;

//Dummy class for storing the partials at the node
public class STAMPartialStorage {

    public double [] partialVector;

    STAMPartialStorage(int N){
        this.partialVector = new double[N * 6];
    }

}