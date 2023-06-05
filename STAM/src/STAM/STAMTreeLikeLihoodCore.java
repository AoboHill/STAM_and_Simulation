package STAM;

import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.tree.Node;

import java.util.Arrays;

public class STAMTreeLikeLihoodCore extends BeerLikelihoodCore {

    public boolean m_bReuseCache = false;
    final double MIN_STEP = 0.0001;
    final static boolean debug = false;

    // 2 x #nodes x #patterns at bottom of branch
    STAMMomentApproximation[][][] mom;
    public double[][] time;
    int N;

    double[] v1;
    double[] v1cache;
    double[] v2;
    double[] v2cache;
    double[] Q1;
    double[] Q2;
    double[] delta;
    double theta;



    /**
     * maps state to list of sitepatterns that contain that state
     * stateMap[nodeIndex][leaf state value][pattern]
     **/
    int[][][] stateMap;
    int[][][][] stateMap_1;


    public STAMTreeLikeLihoodCore(Node root, Alignment data, int N, double theta) {
        this(root.getNodeCount(), data.getPatternCount(), N,theta);
    }

    public STAMTreeLikeLihoodCore(int nodeCount, int patternCount, int N, double theta) {
        super((N+1) * 6 + 4);
        this.theta = theta;
        this.N = N;
        mom = new STAMMomentApproximation[2][nodeCount][patternCount];
        for (int i = 0; i < nodeCount; i++) {
            for (int j = 0; j < patternCount; j++) {
                mom[0][i][j] = new STAMMomentApproximation(N);
                mom[1][i][j] = new STAMMomentApproximation(N);
            }
        }
        v1 = new double[(N+1) * 6 + 4];
        v1cache = new double[((N+1) * 6 + 4)* patternCount ];
        v2 = new double[(N+1) * 6 + 4];
        v2cache = new double[((N+1) * 6 + 4) * patternCount];
        Q1 = new double[((N+1) * 6 + 4) * ((N+1) * 6 + 4)];
        Q2 = new double[((N+1) * 6 + 4) * ((N+1) * 6 + 4)];


    }

    @Override
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories,
                           boolean useAmbiguities) {
        super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);
        time = new double[nodeCount][matrixCount];
    }


    public void updatePartial(double time, double[] matrix_, double[] a){

        Array2d transitionMatrix = new Array2d((N+1) * 6 + 4, (N+1) * 6 + 4, matrix_);

        double[] copy = transitionMatrix.mulrowVectorLeft(a);

        System.arraycopy(copy, 0, a, 0, (N+1) * 6 + 4);

    }

//    public void updatePartial(double time, double[] matrix_, double[] a){
//        double[][] doubleDMatrix = to2DMatrix(matrix_);
//        System.out.println(Arrays.deepToString(doubleDMatrix));
//        System.out.println(Arrays.toString(a));
//        double[] copy = LUDecomposition.solveLU(doubleDMatrix,a);
//        System.out.println(Arrays.toString(copy));
//        System.arraycopy(copy, 0, a, 0, N * 6);
//    }


    public double[][] to2DMatrix(double[] matrix){
        double[][] it = new double[(N+1) * 6 + 4][(N+1) * 6 + 4];
        for (int i = 0; i < (N+1) * 6 + 4; i++){
            System.arraycopy(matrix,0,it[i],0,it.length);
        }

        return it;
    }


    public void setLeafPartials(int nodeIndex,int[][][] countsPatterns) {

        for (int k = 0; k < nrOfPatterns; k++) {
            //System.out.println( "Info about input " + Arrays.toString(countsPatterns[k][nodeIndex]));
            mom[currentPartialsIndex[nodeIndex]][nodeIndex][k].partialVector = getTipLikelihoods(countsPatterns[k][nodeIndex]);
            //System.out.println("info about output " + Arrays.toString(getTipLikelihoods(countsPatterns[k][nodeIndex], N)));
            //System.out.println(Arrays.toString(countsPatterns[k][nodeIndex]));
            //System.out.println(Arrays.toString(mom[currentPartialsIndex[nodeIndex]][nodeIndex][k].partialVector));

        }

    }

    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {
        calculateParentsPartials(nodeIndex1, mom[currentPartialsIndex[nodeIndex1]][nodeIndex1],
                matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1], nodeIndex2,
                mom[currentPartialsIndex[nodeIndex2]][nodeIndex2],
                matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                mom[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
    }

    public void calculateParentsPartials(int nodeIndex1, STAMMomentApproximation[] partials1, double[] matrices1,
                                         int nodeIndex2, STAMMomentApproximation[] partials2, double[] matrices2, STAMMomentApproximation[] partials3
    ){

        for (int l = 0; l < nrOfMatrices; l++) {
            int w = l * matrixSize;
            System.arraycopy(matrices1, w, Q1, 0, matrixSize);
            System.arraycopy(matrices2, w, Q2, 0, matrixSize);
            for (int k = 0; k < nrOfPatterns; k++) {

                System.arraycopy(partials1[k].partialVector, 0, v1, 0, N * 6);
                updatePartial(time[nodeIndex1][l], Q1, v1);

                System.arraycopy(partials2[k].partialVector, 0, v2, 0, N * 6);
                updatePartial(time[nodeIndex2][l], Q2, v2);

                for (int i = 0; i < (N+1) * 6 + 4; i++) {
                    partials3[k].partialVector[i] = Math.max(v1[i] * v2[i], 0);
                    //System.out.println(Arrays.toString(partials3[k].partialVector));
                    // f[i] = (v1[i] * v2[i]);
                }
                //System.out.println(Arrays.toString(partials3[k].partialVector));
            }
        }

    }

    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
        calculateIntegratePartials(mom[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
    }


    protected void calculateIntegratePartials(STAMMomentApproximation[] inPartials, double[] proportions,
                                              double[] outPartials) {

        int u = 0;
        // int v = 0;
        double[] f;
        for (int k = 0; k < nrOfPatterns; k++) {
            f = inPartials[k].partialVector;
            for (int i = 0; i < nrOfStates; i++) {
                outPartials[u] = f[i] * proportions[0];
                u++;
            }
        }
        //System.out.println(nrOfMatrices);

        for (int l = 1; l < nrOfMatrices; l++) {
            u = 0;

            for (int k = 0; k < nrOfPatterns; k++) {
                f = inPartials[k].partialVector;

                for (int i = 0; i < nrOfStates; i++) {

                    System.out.println(proportions[l]);
                    outPartials[u] += f[i] * proportions[l];

                    u++;
                }
            }
        }
    }

    public void calculateLogLikelihoods_1(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {
            // System.out.println("here");
            // System.out.println(c.a[0]);
            double sum = 0;
            for (int i = 1; i < (N+1) * 6 + 4; i += 1) {
                // System.out.println(frequencies[i]);
                // System.out.println(i);
                // sum += Math.max(0, c.a[i] * frequencies[i]);
                sum += partials[v] / ((N+1) * 6 + 4);
                v ++;
            }
            //System.out.println("sum");
            //System.out.println(sum);
            if (sum <= 0) {
                // sum should be > 0, but due to numerical issues can be <= 0
                //System.out.println(Arrays.toString(partials));
                if (sum <= -0.1) {
                    // should never get here since this indicates very bad numerical behaviour
                    //System.out.println("sum: " + sum);
                }
                sum = 1e-200;
            }
            //System.out.println(sum);
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
        //System.out.println(Arrays.toString(outLogLikelihoods));
    }



    public static double cumprod(double[] x) {
        double result = 1;
        for (double j : x) {
            result *= j;
        }
        return result;
    }

    public static double cumprod1(double x){
        double result = 1;
        for (double i = 1; i <= x; i++){
            result = result * i;
        }
        return result;
    }
    public static double multinomialPdf(int[] sample, double[] probs){

        double[] container = new double[]{cumprod1(sample[0]),cumprod1(sample[1]),
                cumprod1(sample[2]),cumprod1(sample[3])};

        return (cumprod1(Arrays.stream(sample).sum())/cumprod(container))
                * Math.pow(probs[0],sample[0]) * Math.pow(probs[1],sample[1])
                * Math.pow(probs[2],sample[2]) * Math.pow(probs[2],sample[3]);


    }




    public static double[][] getTipBins(double nBins, double[][] CL_Points){
        double[][] tipBins = new double[(int) nBins][4];
        int[][] combinations = new int[][]{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        int[][] complementary = new int[][]{{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}};

        tipBins[0] = new double[]{1.0,0,0,0};
        tipBins[1] = new double[]{0,1.0,0,0};
        tipBins[2] = new double[]{0,0,1.0,0};
        tipBins[3] = new double[]{0,0,0,1.0};


        int index = 0;


        int i = 0;
        for (int j = 4; j < nBins; j++){
            tipBins[j][combinations[i][0]] = (CL_Points[index][0] + CL_Points[index][1])/2;
            tipBins[j][combinations[i][1]] = 1 - (CL_Points[index][0] + CL_Points[index][1])/2;
            tipBins[j][complementary[i][0]] = 0;
            tipBins[j][complementary[i][1]] = 0;
            index = index + 1;
            if(index >= CL_Points.length){
                index = 0;
            }
            if ((j - 3) % CL_Points.length ==0){
                i = i + 1;
            }

        }





        return tipBins;

    }

    public static double[][] CL_bins(int nsize, double M) {
        double[][] mymatrix = new double[nsize + 1][3];
        double[] myc = new double[nsize];

        for (int n = 0; n < nsize; n++) {
            myc[n] = Math.cos((2.0 * (n + 1) - 1.0) * Math.PI / (2.0 * nsize));
        }

        double b = 1.0 - 1.0/M;
        double a = 1.0/M;
        for (int n = 0; n < nsize; n++) {
            myc[n] = (b - a) / 2.0 * myc[n] + (a + b) / 2.0;
        }

        // Reverse the array myc
        for (int i = 0; i < myc.length / 2; i++) {
            double temp = myc[i];
            myc[i] = myc[myc.length - 1 - i];
            myc[myc.length - 1 - i] = temp;
        }

        double i = 1.0/M;
        for (int n = 0; n < nsize; n++) {
            mymatrix[n][0] = i;
            mymatrix[n][1] = myc[n];
            i = myc[n];
        }

        mymatrix[nsize][0] = i;
        mymatrix[nsize][1] = 1.0 - 1.0/M;

        for (int j = 0; j < mymatrix.length; j++) {
            mymatrix[j][2] = mymatrix[j][1] - mymatrix[j][0];
        }

        return mymatrix;
    }

    /**
     * Returns P(tip rate | beta0, beta1, epsilon, trait0, trait1) ~
     * Gaussian(mean = beta0 * trait0 + beta1 * trait1 + ... + epsilon, sd = epsilon)
     * within the tip rate interval (a,b)
     * @param aimBin number of bins and its values to calculate the tip likelihoods
     * @param numNucleotide number of alleles
     * @return tip likelihood between intervals a and b
     */
    public double getTipLikelihood(double[] aimBin, int[] numNucleotide) {
        return multinomialPdf(numNucleotide,aimBin);
    }

    public double[] getTipLikelihoods(int[] numNucletide){
        double[][] tipBins = getTipBins((N+1) * 6 + 4,CL_bins(N,theta));
        double[] tipLikelihoods = new double[(N+1) * 6 + 4];
        int index = 0;
        for (double[] eachBin : tipBins){
            tipLikelihoods[index] = multinomialPdf(numNucletide, eachBin);
            index += 1;

        }

        return tipLikelihoods;

    }


}
