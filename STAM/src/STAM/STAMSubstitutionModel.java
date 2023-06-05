package STAM;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

public class STAMSubstitutionModel extends SubstitutionModel.Base {
    public Input<RealParameter> m_k = new Input<RealParameter>("Kimura_k", "Instantaneous rate of mutating from the A to the C");

    public Input<RealParameter> m_u = new Input<RealParameter>("Kimura_u", "Instantaneous rate of mutating from the A to the T");



    public Input<RealParameter> thetaInput = new Input<RealParameter>("theta", "population size parameter with one value for each node in the tree");
    public Input<RealParameter> coalesentRateInput = new Input<RealParameter>("coalescentRate", "coalescent rate parameter with one value for each node in the tree", Input.Validate.XOR, thetaInput);

    public STAMSubstitutionModel() {
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

    }

    @Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {}

    @Override
    public double[] getFrequencies() {return null;}

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {return null;}

    @Override
    public boolean canReturnComplexDiagonalization() {return false;}

    @Override
    public boolean canHandleDataType(DataType dataType) {return true;}
}
