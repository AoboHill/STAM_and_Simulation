package STAM;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;


@Description("Standard prior for SnAP analysis, consisting of a Yule prior on the tree " +
        "(parameterized by lambda) " +
        "and gamma distribution over the theta values " +
        "(with parameters alpha and beta). " +
        "Thetas are represented by the coalescenceRate parameter where values are theta=2/coalescenceRate")
public class STAMPrior extends Distribution {
    public Input<RealParameter> m_pAlpha = new Input<RealParameter>("alpha", "Alpha parameter for the gamma prior on population size (theta) values", Validate.REQUIRED);
    public Input<RealParameter> m_pBeta = new Input<RealParameter>("beta", "Beta parameter for the gamma prior on population size (theta) values", Validate.REQUIRED);
    public Input<RealParameter> m_pKappa = new Input<RealParameter>("kappa", "prior parameter -- see docs for details");
    public Input<RealParameter> m_ptheta = new Input<RealParameter>("theta", "Populations sizes for the nodes in the tree", Validate.REQUIRED);
    public Input<RealParameter> m_pLambda = new Input<RealParameter>("lambda", "Birth rate for the Yule model prior on the species tree");//, Validate.REQUIRED);
    public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations"); //, Validate.REQUIRED);


    enum Priors {
        gamma,inverseGamma,CIR,uniform
    }

    public Input<STAMPrior.Priors> m_priors = new Input<STAMPrior.Priors>("rateprior", "prior on rates. " +
            "This can be " + Arrays.toString(STAMPrior.Priors.values()) + " (default 'gamma')", STAMPrior.Priors.gamma, STAMPrior.Priors.values());

    int PRIORCHOICE = 0;

    @Override
    public void initAndValidate() {
        if (m_pKappa.get() == null) {
            System.err.println("WARNING: kappa parameter not set for SnAPPrior. using default value of 1.0");
            m_pKappa.setValue(new RealParameter(new Double[]{1.0}), this);
        }

        // determine rate prior
        switch (m_priors.get()) {
            case gamma: PRIORCHOICE =0;break;
            case inverseGamma: PRIORCHOICE =1;break;
            case CIR: PRIORCHOICE =2;break;
            case uniform: PRIORCHOICE =3;break;
        }
        System.err.println("Rate prior = " + STAMPrior.Priors.values()[PRIORCHOICE] + "");
    }


    @Override
    public double calculateLogP() {
        logP = 0.0;

        double alpha = m_pAlpha.get().getValue();
        double beta = m_pBeta.get().getValue();
        double kappa = m_pKappa.get().getValue();

        if (outsideBounds(m_pAlpha.get()) || outsideBounds(m_pBeta.get()) || outsideBounds(m_pKappa.get()) || outsideBounds(m_ptheta.get())) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        Tree tree = m_pTree.get();
        double heightsum = tree.getRoot().getHeight();
        heightsum += heightSum(tree.getRoot());
        int nspecies = (tree.getNodeCount() + 1) / 2;
        double lambda = m_pLambda.get().getValue();

        double mu = 0.0; //Death process

        if (mu==0.0) {
            //Yule Model
            logP += (nspecies-1)*Math.log(lambda) - lambda*heightsum;
        } else {
            //Birth death model.    See Thompson 1975, pg 56
            List<Double> allHeights = getSortedHeights(tree.getRoot());
            Iterator<Double> p = allHeights.iterator();
            double x1 = p.next();
            double p0n= p0(x1,lambda,mu);

            for(int n=2; n<nspecies;n++) {
                double xn = p.next();
                logP += Math.log(mu*p1(xn,lambda,mu)/p0n);
            }
        }


        //Gamma values in tree
        RealParameter theta = m_ptheta.get();



        if (PRIORCHOICE == 0) {
            //Assume independent gamma distributions for thetas.

            //We assume that 2/r has a gamma(alpha,beta) distribution. That means that r has density proportional to
            // 1/(r^2)  * GAMMA(2/r|alpha,beta)
            //which has log (alpha - 1.0)*Math.log(2.0/r) - (beta *(2.0/ r)) - 2*log(r), which in turn simplifies to the expr. below (w/ consts)

            for (int iNode = 0; iNode < theta.getDimension(); iNode++) {
                double r = theta.getValue(iNode);
                logP += (alpha-1) * Math.log(r) - r * beta;
            }
        } else if (PRIORCHOICE == 1) {

        } else if (PRIORCHOICE == 2) {



        } else {

        }






        return logP;
    } // calculateLogLikelihood

    private boolean outsideBounds(RealParameter realParameter) {
        if (realParameter == null) {
            return false;
        }
        for (int i = 0; i < realParameter.getDimension(); i++) {
            double d = realParameter.getArrayValue(i);
            if (d < realParameter.getLower() || d > realParameter.getUpper()) {
                return true;
            }
        }
        return false;
    }

    double heightSum(Node node) {
        if (node.isLeaf()) {
            return 0;
        } else {
            double h = node.getHeight();
            h += heightSum(node.getLeft());
            if (node.getRight() != null) {
                h += heightSum(node.getRight());
            }
            return h;
        }
    } // heightSum

    //Returns a list of branching times in the tree, sorted in an decreasing sequence. First one is
    //the height of the mrca of the tree.
    List<Double> getSortedHeights(Node node) {
        return null;
    }

    private double p0(double t, double lambda, double mu) {
        return mu*(1-Math.exp(-(lambda-mu)*t))/(lambda - mu*Math.exp(-(lambda-mu)*t));
    }

    private double p1(double t, double lambda, double mu) {
        double denominator = (lambda - mu*Math.exp(-(lambda-mu)*t));
        return (lambda - mu)*(lambda - mu) * Math.exp(-(lambda-mu)*t) / (denominator * denominator);
    }

    @Override public List<String> getArguments() {return null;}
    @Override public List<String> getConditions() {return null;}
    @Override public void sample(State state, Random random) {};
} // class SSSPrior

