package STAM;

import org.apache.commons.math.special.Gamma;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import static java.lang.Math.abs;

public class STAMProbabilityVA {

    /*
    // This function can output the estimations of the parameters
    public double[] HierarchBetaApprox(double[] Ext, Array2d Varxt){
        double muOmega = Ext[0] + Ext[1];
        double muEta1 = Ext[0]/muOmega;
        double muEta2 = Ext[2]/(1-muOmega);
        double VarOmega = Varxt.at(1,1) + Varxt.at(2,2) + 2 * Varxt.at(1,2);
        double VarEta1 = (Varxt.at(1,1)- VarOmega * (muEta1 * muEta1))/(VarOmega + muOmega * muOmega);
        double VarEta2 = (Varxt.at(3,3)- VarOmega * (muEta2 * muEta2))/(VarOmega + (1-muOmega) * (1-muOmega));
        double phiOmega = muOmega * (1-muOmega)/VarOmega - 1;
        double phiEta1 = muEta1 * (1-muEta1)/VarEta1 - 1;
        double phiEta2 = muEta2 * (1-muEta2)/VarEta2 - 1;

        double[] parameterList = new double[6];
        parameterList[0] = muOmega;
        parameterList[1] = VarOmega;
        parameterList[2] = muEta1;
        parameterList[3] = VarEta1;
        parameterList[4] = muEta2;
        parameterList[5] = VarEta2;

        return parameterList;

    }

     */

    // This function gives the transition probability

    double N;

    double alphaAG;
    double betaAG;
    double alphaA;
    double betaA;
    double alphaC;
    double betaC;

    double[] Ext;

    Array2d Varxt;

    double BetaF_AG;

    double BetaF_A;

    double BetaF_C;

    double muOmega;
    double muEta1;
    double muEta2;
    double VarOmega;
    //double VarEta1 = (Varxt.at(1,1)- VarOmega * (muEta1 * muEta1))/(VarOmega + muOmega * muOmega);
    double VarEta1;
    double VarEta2;

    double phiOmega;
    double phiEta1;
    double phiEta2;

    double normailzeConstant1;

    double normailzeConstant2;

    double normailzeConstant3;

    double normailzeConstant4;

    double normailzeConstant5;

    double normailzeConstant6;

    double theta;

    double prob_spike0AG;
    double prob_spike0A;
    double prob_spike0C;

    double prob_spike1AG;
    double prob_spike1A;
    double prob_spike1C;



    STAMProbabilityVA(double[] Ext, Array2d Varxt, double N, double theta) throws MathException {
        this.N = N;
        this.Ext = Ext;
        this.Varxt = Varxt;
        this.muOmega = Ext[0] + Ext[1];
        this.muEta1 = Ext[0]/muOmega;
        this.muEta2 = Ext[2]/(1-muOmega);
        this.VarOmega = Varxt.at(1,1) + Varxt.at(2,2) + 2 * Varxt.at(1,2);
        //double VarEta1 = (Varxt.at(1,1)- VarOmega * (muEta1 * muEta1))/(VarOmega + muOmega * muOmega);
        this.VarEta1 = (Varxt.at(1,1)+Varxt.at(2,2)-
                (muEta1*muEta1+(1-muEta1)*(1-muEta1))*VarOmega)
                /(2*(muOmega*muOmega+VarOmega));
        this.VarEta2 = (Varxt.at(3,3)+Varxt.at(4,4)-
                (muEta2*muEta2+(1-muEta2)*(1-muEta2))*VarOmega)
                /(2*((1-muOmega)*(1-muOmega)+VarOmega));

        //double VarEta2 = (Varxt.at(3,3)- VarOmega * (muEta2 * muEta2))/(VarOmega + (1-muOmega) * (1-muOmega));



        //if (VarOmega > muOmega*(1-muOmega)){
        //    VarOmega = muOmega*(1-muOmega) - 0.0000000000000001;
        //}

        //if (VarEta1 > muEta1*(1-muEta1)){
        //    VarEta1 = muEta1*(1-muEta1) - 0.0000000000000001;
        //}

        //if (VarEta2 > muEta2*(1-muEta2)){
        //    VarEta2 = muEta2*(1-muEta2) - 0.0000000000000001;
        //}





        this.phiOmega = abs(muOmega * (1-muOmega)/VarOmega - 1);
        this.phiEta1 = abs(muEta1 * (1-muEta1)/VarEta1 - 1);
        this.phiEta2 = abs(muEta2 * (1-muEta2)/VarEta2 - 1);



        this.alphaAG = muOmega*phiOmega;
        this.betaAG = phiOmega*(1-muOmega);
        this.alphaA = muEta1*phiEta1;
        this.betaA = phiEta1*(1-muEta1);
        this.alphaC = muEta2*phiEta2;
        this.betaC = phiEta2*(1-muEta2);

        this.BetaF_AG = BetaF(alphaAG,betaAG);

        this.BetaF_A = BetaF(alphaA,betaA);

        this.BetaF_C = BetaF(alphaC,betaC);




        this.theta = theta;

        this.prob_spike0AG = prob_spike0(alphaAG,betaAG);

        this.prob_spike0A = prob_spike0(alphaA,betaA);

        this.prob_spike0C = prob_spike0(alphaC,betaC);


        this.prob_spike1AG = prob_spike1(alphaAG,betaAG);

        this.prob_spike1A = prob_spike1(alphaA,betaA);

        this.prob_spike1C = prob_spike1(alphaC,betaC);



    }

    public double BetaF(double alpha, double beta){
        return Math.exp(Gamma.logGamma(alpha) + Gamma.logGamma(beta) - Gamma.logGamma(alpha + beta));
    }

    public double unNormalizedDensity(double aim, double alpha, double beta){
        return Math.exp((alpha - 1) * Math.log(aim) + (beta - 1) * Math.log(1-aim));
    }

    public double prob_spike0(double alpha, double beta) throws MathException {
        if (alpha > Math.pow(10, 3) && beta > Math.pow(10, 3)) {
            return 0;
        } else if (alpha < Math.pow(10, 3) && beta < Math.pow(10, 3)){
            return Math.pow(theta, -alpha) / (alpha * BetaF(alpha, beta));
        }
        BetaDistributionImpl betadist = new BetaDistributionImpl(alpha,beta);
        return betadist.cumulativeProbability(1/theta);
    }

    public double prob_spike1(double alpha, double beta) throws MathException {
        if (alpha > Math.pow(10, 3) && beta > Math.pow(10, 3)) {
            return 0;
        } else if (alpha < Math.pow(10, 3) && beta < Math.pow(10, 3)){
            return Math.pow(theta, -beta) / (beta * BetaF(alpha, beta));
        }
        BetaDistributionImpl betadist = new BetaDistributionImpl(alpha,beta);
        return 1 - betadist.cumulativeProbability(1-1/theta);

    }




    public double HierarchBeta(double[] alleleFreq , double theta){

        if (alleleFreq[0] == 1 && alleleFreq[1] == 0 && alleleFreq[2] == 0 && alleleFreq[3] == 0){
            return this.prob_spike1AG * this.prob_spike1A;
        } else if (alleleFreq[0] == 0 && alleleFreq[1] == 1 && alleleFreq[2] == 0 && alleleFreq[3] == 0){
            return this.prob_spike1AG * this.prob_spike0A;
        } else if (alleleFreq[0] == 0 && alleleFreq[1] == 0 && alleleFreq[2] == 1 && alleleFreq[3] == 0){
            return this.prob_spike0AG * this.prob_spike1C;
        } else if(alleleFreq[0] == 0 && alleleFreq[1] == 0 && alleleFreq[2] == 0 && alleleFreq[3] == 1){
            return this.prob_spike0AG * this.prob_spike0C;
        } else if(alleleFreq[2] != 0 && alleleFreq[3] != 0){
            return unNormalizedDensity(alleleFreq[2],alphaC,betaC)/BetaF_C * prob_spike0AG;
        } else if(alleleFreq[0] != 0 && alleleFreq[1] != 0){
            return unNormalizedDensity(alleleFreq[0],alphaA,betaA)/BetaF_A * prob_spike1AG;
        } else if(alleleFreq[0] + alleleFreq[2] == 0){
            return unNormalizedDensity(alleleFreq[1],alphaAG,betaAG)/BetaF_AG *
                    prob_spike0A * prob_spike0C / (alleleFreq[1] *(1-alleleFreq[1]));
        } else if(alleleFreq[0] + alleleFreq[3] == 0){
            return unNormalizedDensity(alleleFreq[1],alphaAG,betaAG)/BetaF_AG *
                    prob_spike0A * prob_spike1C / (alleleFreq[1] *(1-alleleFreq[1]));
        }else if(alleleFreq[1] + alleleFreq[2] == 0){
            return unNormalizedDensity(alleleFreq[0],alphaAG,betaAG)/BetaF_AG *
                    prob_spike1A * prob_spike0C / (alleleFreq[0] *(1-alleleFreq[0]));
        }else if(alleleFreq[1] + alleleFreq[3] == 0){
            return unNormalizedDensity(alleleFreq[0],alphaAG,betaAG)/BetaF_AG *
                    prob_spike1A * prob_spike1C / (alleleFreq[0] *(1-alleleFreq[0]));
        }
        return 0;




    }


    /*
    public static void main(String[] args) throws Exception {
        Array2d Q = new Array2d(4,4,new double[]{-0.0006250,  0.0001875,  0.0002500,  0.0001875,
                0.0001875, -0.0006250,  0.0001875,  0.0002500,
                0.0002500,  0.0001875, -0.0006250,  0.0001875,
                0.0001875,  0.0002500,  0.0001875, -0.0006250});

        MomentApproximation mom = new MomentApproximation(new double[]{0.1,0.2,0.3,0.4},0.5,Q);
        mom.Kimura();
        Probability prob = new Probability();
        double density = prob.HierarchBeta(mom.Ext,mom.Varxt, new double[]{0.1,0.2,0.3,0.4});

    }
    * */


}
