package STAM;

import java.util.Arrays;

public class STAMTipLikelihood {



    private int numBins;


    public STAMTipLikelihood(int numBins) {
        this.numBins = numBins;
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



    public static double[][][] getTipBins(double nBins){
        double[][][] tipBins = new double[6][(int) nBins][4];
        int[][] combinations = new int[][]{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        int[][] complementary = new int[][]{{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}};


        for (int i = 0; i < 6; i++){
            for (int j = 0; j < nBins; j++){
                tipBins[i][j][combinations[i][0]] = 1/nBins * (j) + 1/(2*nBins);
                tipBins[i][j][combinations[i][1]] = 1 - 1/nBins * (j) - 1/(2*nBins);
                tipBins[i][j][complementary[i][0]] = 0;
                tipBins[i][j][complementary[i][1]] = 0;
            }
        }
        return tipBins;

    }


    public static double[][][] getTipBins1(double nBins){
        double[][][] tipBins = new double[6][(int) nBins][4];
        int[][] combinations = new int[][]{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        int[][] complementary = new int[][]{{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}};


        for (int i = 0; i < 6; i++){
            for (int j = 0; j < nBins; j++){
                tipBins[i][j][combinations[i][0]] = 1/nBins * (j) + 1/(2*nBins);
                tipBins[i][j][combinations[i][1]] = 1 - 1/nBins * (j) - 1/(2*nBins);
                tipBins[i][j][complementary[i][0]] = 0;
                tipBins[i][j][complementary[i][1]] = 0;
            }
        }
        return tipBins;

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

    public double[] getTipLikelihoods(int[] numNucletide, int nBins){
        double[][][] tipBins = getTipBins(numBins);
        double[] tipLikelihoods = new double[6 * nBins];
        int index = 0;
        for (double[][] eachCombination : tipBins){
            for (double[] eachBin : eachCombination){
                tipLikelihoods[index] = multinomialPdf(numNucletide, eachBin);
            }
        }

        return tipLikelihoods;

    }




}
