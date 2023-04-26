package STAM;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import snap.NodeData;
import snapper.SnapperTreeLikelihood;
import STAM.STAMTreeLikeLihoodCore;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;

@Description("Tree Likelihood Function for Single Site Nucleotide Sequences on a tree.")

public class STAMTreeLikelihood extends TreeLikelihood{
    public Input<Integer> NInput = new Input<Integer>("N", "Number of bins", 10);

    public Input<Boolean> m_bInitFromTree = new Input<Boolean>("initFromTree", "whether to initialize coalescenceRate from starting tree values (if true), or vice versa (if false)", false);
    public Input<String> m_pPattern = new Input<String>("pattern", "pattern of metadata element associated with this parameter in the tree");


    // Snapper paramter but not used in STAM
    public Input<Boolean> useLogLikelihoodCorrection = new Input<Boolean>("useLogLikelihoodCorrection", "use correction of log likelihood for the purpose of calculating " +
            "Bayes factors for different species assignments. There is (almost) no computational cost involved for the MCMC chain, but the log likelihood " +
            "might be reported as positive number with this correction since the likelihood is not a proper likelihood any more.", true);

    public Input<Boolean> showPatternLikelihoodsAndQuit = new Input<Boolean>("showPatternLikelihoodsAndQuit", "print out likelihoods for all patterns for the starting state, then quit", false);
    public Input<Boolean> mutationOnlyAtRoot = new Input<Boolean>("mutationOnlyAtRoot", "Emulate the likelihood calculation of RoyChoudhury et al (2008) which assumes that mutations occur only in the ancestral (root) population", false);
    public Input<Boolean> hasDominantMarkers = new Input<Boolean>("dominant", "indicate that alleles are dominant (default false)", false);

    public Input<Boolean> m_usenNonPolymorphic = new Input<Boolean>("non-polymorphic",
            "Check box if there was no pre-filtering of sites to remove all constant sites. " +
                    "Leave unchecked if constant sites had been removed or systematically not selected (e.g. SNP data). The likelihoods will be adjusted according.",
            //"Whether to use non-polymorphic data in the sequences. " +
            //"If true, constant-sites in the data will be used as part of the likelihood calculation. " +
            //"If false (the default) constant sites will be removed from the sequence data and a normalization factor is " +
            //"calculated for the likelihood.",
            true);

    public Input<Integer> m_numUnfilteredSitesInput = new Input<Integer>("number of sites which were not filtered to remove constant sites",  "Number of sites not pre-filtered.  (default =0). This setting ignored unless non-polymorphic set to TRUE", 0);


    final private Input<List<TreeLikelihood>> likelihoodsInput = new Input<>("*","",new ArrayList<>());

    public STAMTreeLikelihood() throws Exception {
        // suppress some validation rules
        siteModelInput.setRule(Validate.OPTIONAL);
    }

    /** shadow variable of m_pData input */
    protected Alignment m_data2;

    protected STAMData m_stamdata;

    /** SampleSizes = #lineages per taxon **/
    int [] m_nSampleSizes;
    /** likelihood core, doing the actual hard work of calculating the likelihood **/
    STAMTreeLikeLihoodCore m_core;

    /** some variable for shadowing inputs **/
    boolean m_bUsenNonPolymorphic;
    int 	m_numUnfilteredSites;
    boolean m_bMutationOnlyAtRoot;
    boolean m_bHasDominantMarkers;

    double m_fP0 = 0.0, m_fP1 = 0.0;
    double m_fStoredP0 = 0.0, m_fStoredP1 = 0.0;
    double ascLogP = Double.NaN, storedAscLogP = Double.NaN;

    int N;

    STAMSubstitutionModel m_substitutionmodel;


    // Correction so that the returned value is a likelihood instead
    // of a sufficient statistic for the likelihood
    protected double m_fLogLikelihoodCorrection = 0;

    // Sampled parameter equal to the number of sites which have been removed from the data during ascertainment
    IntegerParameter ascSiteCount;



    private ExecutorService pool = null;
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();


    @Override
    public void initAndValidate() {
        // check that alignment has same taxa as tree
        //if (!(dataInput.get() instanceof Data)) {
        //	throw new IllegalArgumentException("The data input should be a snap.Data object");
        //}



        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of leaves in the tree does not match the number of sequences");
        }

        m_bUsenNonPolymorphic = m_usenNonPolymorphic.get();
        m_numUnfilteredSites = m_numUnfilteredSitesInput.get();
        m_bMutationOnlyAtRoot = mutationOnlyAtRoot.get();
        m_bHasDominantMarkers = hasDominantMarkers.get();


        m_siteModel = (SiteModel.Base) siteModelInput.get();

        N =  NInput.get();

        TreeInterface tree = treeInput.get();

        m_substitutionmodel = ((STAMSubstitutionModel)m_siteModel.substModelInput.get());


        if (m_substitutionmodel.thetaInput.get() != null) {
            // parameterised as thetas
            Input<RealParameter> thetaInput = m_substitutionmodel.thetaInput;

            Double [] values = new Double[tree.getNodeCount()];
            String sTheta = "";
            if (m_bInitFromTree.get() == true) {
                tree.getMetaData(tree.getRoot(), values, m_pPattern.get());
                for (Double d : values) {
                    sTheta += d + " ";
                }
            } else {
                List<Double> sValues = thetaInput.get().valuesInput.get();
                for (int i = 0; i < values.length; i++) {
                    values[i] = sValues.get(i % sValues.size());
                    sTheta += values[i] + " ";
                }
                tree.setMetaData(tree.getRoot(), values, m_pPattern.get());
            }
            RealParameter pTheta = thetaInput.get();
            RealParameter theta = new RealParameter();
            theta.initByName("value", sTheta, "upper", pTheta.getUpper(), "lower", pTheta.getLower(), "dimension", values.length);
            theta.setID(pTheta.getID());
            thetaInput.get().assignFrom(theta);
        } else {
            // parameterised as coalescentRates
            Input<RealParameter> coalesecntRateInput = m_substitutionmodel.coalesentRateInput;

            Double [] values = new Double[tree.getNodeCount()];
            String sRates = "";
            if (m_bInitFromTree.get() == true) {
                tree.getMetaData(tree.getRoot(), values, m_pPattern.get());
                for (Double d : values) {
                    sRates += 2.0/d + " ";
                }
            } else {
                List<Double> sValues = coalesecntRateInput.get().valuesInput.get();
                for (int i = 0; i < values.length; i++) {
                    values[i] = sValues.get(i % sValues.size());
                    sRates += values[i] + " ";
                }
                tree.setMetaData(tree.getRoot(), values, m_pPattern.get());
            }
            RealParameter pCoalescentRate = coalesecntRateInput.get();
            RealParameter coalescentRate = new RealParameter();
            coalescentRate.initByName("value", sRates, "upper", pCoalescentRate.getUpper(), "lower", pCoalescentRate.getLower(), "dimension", values.length);
            coalescentRate.setID(pCoalescentRate.getID());
            coalesecntRateInput.get().assignFrom(coalescentRate);
        }


        m_data2 = dataInput.get();

        // cast it to our own data type
        m_stamdata = (STAM.STAMData) m_data2;


        m_stamdata.initAndValidate();



        Integer [] nSampleSizes = m_data2.getStateCounts().toArray(new Integer[0]);

        m_nSampleSizes = new int[nSampleSizes.length];
        //System.out.println("yes here");
        for (int i = 0; i < nSampleSizes.length; i++) {
            m_nSampleSizes[i] = nSampleSizes[i];
            //System.out.println(nSampleSizes[i]);
        }

        //System.out.println("stop");
        if (!(treeInput.get().getRoot() instanceof NodeData)) {
            throw new IllegalArgumentException("Tree has nodes of the wrong type. NodeData expected, but found " +
                    treeInput.get().getRoot().getClass().getName());
        }

        // Get the number of pattern

        int numPatterns = m_stamdata.getPatternCount();
        //System.out.println("numPatterns");
        //System.out.println(numPatterns);

        /*
        There are potential modifications for the numPatterns, since here we are using the count data.
        * */

        if (branchRateModel == null) {
            branchRateModel = new StrictClockModel();
        }
        if (!(branchRateModel instanceof StrictClockModel)) {
            //We assume that the mutation rate (but not theta) is constant for the species tree.
            throw new IllegalArgumentException("Only strict clock model allowed for branchRateModel, not " + branchRateModel.getClass().getName());
        }

        m_core = new STAMTreeLikeLihoodCore(treeInput.get().getRoot().getNodeCount(), numPatterns,N);




        /* Here we just suppose that we do not need this likelihood correction
        // calculate Likelihood Correction.
        // When the assignment of individuals to populations/species is fixed, the allele counts in each population are sufficient
        // statistics for the species tree parameters. However when testing species assignments this is no longer the case.
        // To address this we multiply the likelihood computed from allele counts by the probability of observing
        // the given sequences given those allele counts (and the species assignments).

        // verified by DJB on 27/7/2021 that this term is correct and required when
        // estimating marginal likelihoods with different species assignments (as
        // for species delimitation using the BFD or BFD* method).
        m_fLogLikelihoodCorrection = 0;
        if (useLogLikelihoodCorrection.get()) {
            // RRB: note that increasing the number of constant sites
            // does not change the m_fLogLikelihoodCorrection since the
            // contribution of constant sites is zero. This means,
            // m_fLogLikelihoodCorrection does not need to be recalculated
            // when ascSiteCount changes.
            // DJB: This is true, but only until we start looking at non-constant sites being ascertained.
            for (int i = 0; i < numPatterns; i++) {
                int [] thisSite = m_data2.getPattern(i);  //count of red alleles for this site
                int [] lineageCounts = getPatternLineagCounts(m_data2, i); //count of total lineages for this site
                for (int j = 0; j < thisSite.length; j++) {
                    m_fLogLikelihoodCorrection -= logBinom(thisSite[j], lineageCounts[j]) * m_data2.getPatternWeight(i);
                }
            }
        }
        System.err.println("Log Likelihood Correction = " + m_fLogLikelihoodCorrection);

        * */



        int nodeCount = tree.getNodeCount();
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        patternLogLikelihoods = new double[numPatterns];
        m_fRootPartials = new double[numPatterns * N * 6];

        //System.out.println(m_siteModel.getCategoryCount());
        m_core.initialize(
                nodeCount,
                numPatterns,
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

        double [] f = new double[N];


        // Get the pattern
        int [][][] thisSite = m_stamdata.getCountsPatterns();



        for (int j = 0; j < thisSite[0].length; j++) {
                m_core.setLeafPartials(j,thisSite);

            }
    }

    @Override
    public double calculateLogP() {

        try {
            // get current tree
            NodeData root = (NodeData) treeInput.get().getRoot();
            //Double [] theta = m_substitutionmodel.thetaInput.get().getValues();
            // assing gamma values to tree
//	    	if (m_pGamma.get().somethingIsDirty()) {
//	    		// sync gammas in parameter with gammas in tree, if necessary
//	    		m_pGamma.get().prepare();
//	    	}

            //double u = m_substitutionmodel.m_pU.get().getValue();
            //double v  = m_substitutionmodel.m_pV.get().getValue();
            boolean useCache = false;
            //boolean useCache = false;
            boolean dprint = showPatternLikelihoodsAndQuit.get();
            if (dprint) {
                System.out.println("Log Likelihood Correction = " + m_fLogLikelihoodCorrection);
            }


            double [] fCategoryRates = m_siteModel.getCategoryRates(null);
            if (branchRateModel != null) {
                double branchRate = branchRateModel.getRateForBranch(null);
                for (int i = 0; i < fCategoryRates.length; i++) {
                    fCategoryRates[i] *= branchRate;
                }
            }



            traverse(root);



            // amalgamate site probabilities over patterns
            int numPatterns = m_stamdata.getPatternCount();
            // claculate log prob
            logP = 0;
            // Site_L
            //System.out.println("PATTERNS");
            for(int id = 0; id < numPatterns - (m_bUsenNonPolymorphic ? 0 : 2); id++) {
                //System.out.println("haha, you arrive here");
                double freq = m_stamdata.getPatternWeight(id);

                double siteL = patternLogLikelihoods[id];

                if (Double.isInfinite(siteL))
                {
                    logP = -10e100;
                    break;
                }
                logP += (double)freq * siteL;

                //	System.out.println(Arrays.toString(m_data2.getPattern(id)) + " " + id + " " + (siteL) + " "+ (freq));
            }
            // correction for constant sites. If we are sampling the numbers of constant sites
            // (stored in ascSiteCount) then we include these probabilities. Otherwise we
            // assume that we want conditional likelihood, in which case we divide
            // by the probability that a site is not ascertained (or more correctly,
            // subtract the log probability.
            if (!m_bUsenNonPolymorphic) {
                m_fP0 =  patternLogLikelihoods[numPatterns - 2];
                m_fP1 =  patternLogLikelihoods[numPatterns - 1];
                if (ascSiteCount != null) {
                    ascLogP = (double)ascSiteCount.getValue(0) * Math.log(m_fP0) +
                            (double)ascSiteCount.getValue(1) * Math.log(m_fP1);
                    logP += ascLogP;
                } else {
                    logP -= (double) (m_stamdata.getSiteCount() - m_numUnfilteredSites) * Math.log(1.0 - m_fP0 - m_fP1); //Correct likelihoods for those sites which were pre-filtered (removing constant sites)
                }
            }

            if (useLogLikelihoodCorrection.get()) {
                logP += m_fLogLikelihoodCorrection;
            }


//			logP = m_core.computeLogLikelihood(root, u , v,
//	    			m_nSampleSizes,
//	    			m_data2,
//	    			coalescenceRate,
//	    			fCategoryRates, fCategoryProportions,
//	    			useCache,
//	    			m_bUsenNonPolymorphic,
//	    			dprint /*= false*/);
            if(Double.isNaN(logP)){
                System.out.println("NaN occur");
                logP = -10e100;
                //logP = 0;
            }
            //System.out.println(logP);
            return logP;

        } catch (Exception e) {
            e.printStackTrace();
            return 0;
        }
    } // calculateLogLikelihood

    @Override
    protected boolean requiresRecalculation() {

        boolean isDirty = super.requiresRecalculation();
        if (ascSiteCount != null && ascSiteCount.somethingIsDirty()) {
            isDirty = true;
        }
        return isDirty;
    }

    /** CalculationNode methods **/
    @Override
    public void store() {
        storedLogP = logP;

        m_core.m_bReuseCache = true;
        m_fStoredP0 = m_fP0;
        m_fStoredP1 = m_fP1;
        storedAscLogP = ascLogP;
        super.store();

    }

    @Override
    public void restore() {
        logP = storedLogP;

        m_core.m_bReuseCache = false;
        m_fP0 = m_fStoredP0;
        m_fP1 = m_fStoredP1;
        ascLogP = storedAscLogP;
        super.restore();
        }

    @Override public List<String> getArguments() {return null;}
    @Override public List<String> getConditions() {return null;}
    @Override public void sample(State state, Random random) {};

    @Override
    public void init(PrintStream out) {
        super.init(out);
        if (!m_bUsenNonPolymorphic) {
            out.append("P0\t");
            out.append("P1\t");
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        super.log(nSample, out);
        if (!m_bUsenNonPolymorphic) {
            out.append(m_fP0 + "\t");
            out.append(m_fP1 + "\t");
        }
    }

    public double getProbVariableSites() {
        if (!m_bUsenNonPolymorphic) {
            return 1.0 - m_fP0 - m_fP1;
        } else {
            return 1.0;
        }
    }

    // return contribution of ascertained sites to log likelihood
    public double getAscSitesLogP() {
        if (ascSiteCount != null) {
            return ascLogP;
        } else {
            return 0.0;
        }
    }


    /* Assumes there IS a branch rate model as opposed to traverse() */
    protected int traverse(final Node node) {



        int update = (node.isDirty() | hasDirt);

        update = Tree.IS_FILTHY;



        final int nodeIndex = node.getNr();


        //System.out.println(branchRate);
        final double branchTime = node.getLength();
        //System.out.println("The time "+ branchTime);


        // This is for stationary distribution calculation
        double[] freq = new double[N*6];
        double[] freq1 = new double[N*6];
        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex] ||
                (m_substitutionmodel.thetaInput.get() != null && m_substitutionmodel.thetaInput.get().isDirty(nodeIndex)) ||
                (m_substitutionmodel.coalesentRateInput.get() != null && m_substitutionmodel.coalesentRateInput.get().isDirty(nodeIndex)) )) {
            //m_substitutionmodel.thetaInput.get().getStoredValue(nodeIndex) != m_substitutionmodel.thetaInput.get().getValue(nodeIndex)
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            m_core.setNodeMatrixForUpdate(nodeIndex);
            Double[] theta = getThetas();
            double k = m_substitutionmodel.m_k.get().getValue();
            double u = m_substitutionmodel.m_u.get().getValue();
            double[] fCategoryRates = m_siteModel.getCategoryRates(null);
            double [] time = m_core.time[nodeIndex];
            double[] waitingForScaled = buildQ1(k,u);
            //System.out.println("before theta" + Arrays.toString(waitingForScaled));
            //System.out.println(Arrays.toString(waitingForScaled));
            double[][][] Bins = STAMTipLikelihood.getTipBins(N);
            double scaledTheta = theta[nodeIndex] / fCategoryRates[0];
            double[] rateMatrix = scaledBytehta(waitingForScaled,scaledTheta);
            time[0] =  branchTime;
            //System.out.println("rate matrix" + Arrays.toString(rateMatrix));
            double[] it = rateToTransition(Bins,time[0],rateMatrix,N);
            double[] stationary = rateToTransition(Bins,100000,rateMatrix,N);
            System.arraycopy(stationary, 0, freq, 0, N * 6);
            System.arraycopy(normalize(freq),0,freq1,0,N*6);
            //System.out.println("transition " + Arrays.toString(it));


            // This part should be used to account for site heterogeneity but currently not used
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                //double scaledTheta = theta[nodeIndex] / fCategoryRates[i];
                //final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                //final double jointBranchRate = m_siteModel.getRateForCategory(i, node) ;
                //System.out.println("jointbranchrate" + jointBranchRate);
                //System.out.println("scale" + m_siteModel.getRateForCategory(i, node));
                //time[i] =  branchTime/scaledTheta;// * scaledTheta / 2;
                //System.out.println("After theta" + Arrays.toString(rateMatrix));
                //System.out.println(Arrays.toString(rateMatrix));
                // System.out.println(node.getNr() + " " + Arrays.toString(Q.Q));
                //System.out.println(Arrays.toString(rateMatrix));
                //System.out.println(branchRate);
                //System.out.println(branchTime);
                //System.out.println(Arrays.toString(rateToTransition(Bins, time[i], rateMatrix, N)));
                //long startTime=System.nanoTime();

                //m_core.setNodeMatrix(nodeIndex, i, rateToTransition(Bins,time[i],rateMatrix,N));
                m_core.setNodeMatrix(nodeIndex, i, it);
                //long endTime=System.nanoTime(); //获取结束时间
                //System.out.println("程序运行时间 （total）： "+(endTime-startTime)+"ns");


            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                m_core.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    m_core.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    //long startTime=System.nanoTime(); //获取结束时间

                    m_core.calculatePartials(childNum1, childNum2, nodeIndex);
                    //long endTime=System.nanoTime(); //获取结束时间
                    //System.out.println("程序运行时间 （along the branch）： "+(endTime-startTime)+"ns");
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {

                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    //System.out.println(Arrays.toString(proportions));
                    m_core.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    //m_core.calculateLogLikelihoods_1(m_fRootPartials, patternLogLikelihoods);

                    m_core.calculateLogLikelihoods_1(m_fRootPartials, freq1, patternLogLikelihoods);
                    //System.out.println(Arrays.toString(m_fRootPartials));
                }

            }
        }
        return update;
    } // traverse

    public double[] rateToTransition(double[][][] Bins, double time, double[] matrix_, int N) {

        double[] matrix = new double [N * 6 * N * 6];

        int index = 0;
        for (double[][] starts : Bins){
            for (double[] start : starts){
                //long startTime=System.nanoTime();
                STAMMomentApproximation mom_ = new STAMMomentApproximation(start,time,new Array2d(4,4,matrix_));
                mom_.Kimura();

                double[] waitToBenormalized = new double[6 * N];
                int count = 0;
               // long startTime1=System.nanoTime();
                STAMProbability prob = new STAMProbability(mom_.Ext, mom_.Varxt,N);
                //long startTime2=System.nanoTime();
                for (double[][] aims: Bins){
                    for (double[] aim : aims){

                        //System.out.println(Arrays.toString(mom_.Ext));
                        //System.out.println(mom_.Varxt);
                        //System.out.println("------------");
                        //System.out.println(Arrays.toString(start));
                        //System.out.println(Arrays.toString(aim));
                        //System.out.println(time);
                        //System.out.println(Arrays.toString(mom_.Ext));
                        //System.out.println(mom_.Varxt);
                        double density = prob.HierarchBeta(aim,N);

                        if (Double.isInfinite(density)){
                            density = 0;
                        }
                        //System.out.println(density);
                        //System.out.println("************");
                        //if (Double.isNaN(density)){
                        //    matrix[index] = 0;
                        //} else{
                        //matrix[index] = density;
                        waitToBenormalized[count] = density;
                        count++;
                        //}
                        //System.out.println("程序运行时间 （calculate beta）： "+(endTime-startTime)+"ns");
                        index ++ ;
                    }

                }

                //System.out.println("before normalize:" + Arrays.toString(waitToBenormalized));
                double[] afterNormalized = normalize(waitToBenormalized);
                //System.out.println(Arrays.toString(afterNormalized));
                System.arraycopy(afterNormalized,0,matrix,index-6*N,6*N);
                //System.out.println(Arrays.toString(matrix));
            }
        }

        return matrix;
    }

    double[] normalize(double[] vector){
        double sum = Arrays.stream(vector).sum();
        if (sum == 0){
            return new double[6 * N];
        }
        double[] v = new double[vector.length];
        for (int i = 0; i < vector.length; i++){
            v[i] = vector[i]/sum;
        }
        return v;
    }

    private Double[] getThetas() {
        Double [] thetas;
        if (m_substitutionmodel.thetaInput.get() != null) {
            //System.out.println("we now at theta");
            thetas = m_substitutionmodel.thetaInput.get().getValues();
        } else {
            //System.out.println("we now at coalesent");
            thetas = m_substitutionmodel.coalesentRateInput.get().getValues();
            for (int i = 0; i < thetas.length; i++) {
                thetas[i] = 2.0/thetas[i];
            }
        }
        return thetas;
    }


    // wait to be changed
    private double[] buildQ(double k, double u){
        double[] U_I = new double[]{


                -(k+2)*u, k * u, u, u,
                k * u, -(k+2)*u, u, u,
                u, u, -(k+2)*u, k * u,
                u, u, k * u, -(k+2)*u
        };



        double[] stationaryDist = new double[]{0.25,0.25,0.25,0.25};

        double[] diagonalOfU_I = new double[]{U_I[0],U_I[5],U_I[10],U_I[15]};



        for (int i = 0; i < stationaryDist.length; i++){
            stationaryDist[i] = stationaryDist[i] * diagonalOfU_I[i];
        }

        //System.out.println(Arrays.toString(stationaryDist));
        double sumStationary = Arrays.stream(stationaryDist).sum();
        //System.out.println(sumStationary);

        for (int i = 0; i < U_I.length; i++){
            U_I[i] = (U_I[i]/(-sumStationary));

        }

        double[] Q = U_I;

        return Q;


    }

    private double[] buildQ1(double k, double u){
        double[] U_I = new double[]{


                -(k+2)*u, k * u, u, u,
                k * u, -(k+2)*u, u, u,
                u, u, -(k+2)*u, k * u,
                u, u, k * u, -(k+2)*u
        };



        return U_I;


    }

    private double[] scaledBytehta(double[] input,double theta){
        double[] x = input;
        for (int i = 0 ; i < input.length; i++){
            input[i] = input[i] * theta;
        }
        return x;
    }



    }


