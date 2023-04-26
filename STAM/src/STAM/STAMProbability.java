package STAM;

import org.apache.commons.math.special.Gamma;

public class STAMProbability {

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


    STAMProbability(double[] Ext, Array2d Varxt, double N){
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





        this.phiOmega = muOmega * (1-muOmega)/VarOmega - 1;
        this.phiEta1 = muEta1 * (1-muEta1)/VarEta1 - 1;
        this.phiEta2 = muEta2 * (1-muEta2)/VarEta2 - 1;



        this.alphaAG = muOmega*phiOmega;
        this.betaAG = phiOmega*(1-muOmega);
        this.alphaA = muEta1*phiEta1;
        this.betaA = phiEta1*(1-muEta1);
        this.alphaC = muEta2*phiEta2;
        this.betaC = phiEta2*(1-muEta2);

        this.BetaF_AG = BetaF(alphaAG,betaAG);

        this.BetaF_A = BetaF(alphaA,betaA);

        this.BetaF_C = BetaF(alphaC,betaC);

        this.normailzeConstant1 = Math.pow(N,-betaAG) / (betaAG * BetaF_AG);

        this.normailzeConstant2 = Math.pow(N,-betaA) / (betaA * BetaF_A)
                * Math.pow(N,-betaC) / (betaC * BetaF_C);

        this.normailzeConstant3 = Math.pow(N,-betaA) / (betaA * BetaF_A)
                * Math.pow(N,-alphaC) / (alphaC * BetaF_C);

        this.normailzeConstant4 =Math.pow(N,-alphaA) / (alphaA * BetaF_A)
                * Math.pow(N,-betaC) / (betaC * BetaF_C);

        this.normailzeConstant5 =Math.pow(N,-alphaA) / (alphaA * BetaF_A)
                * Math.pow(N,-alphaC) / (alphaC * BetaF_C);

        this.normailzeConstant6 = Math.pow(N,-betaAG) / (alphaAG * BetaF_AG);


    }

    public double BetaF(double alpha, double beta){
        return Math.exp(Gamma.logGamma(alpha) + Gamma.logGamma(beta) - Gamma.logGamma(alpha + beta));
    }

    public double unNormalizedDensity(double aim, double alpha, double beta){
        return Math.exp((alpha - 1) * Math.log(aim) + (beta - 1) * Math.log1p(-aim));
    }


    public double HierarchBeta(double[] alleleFreq , double N){

        /*
        double muOmega = Ext[0] + Ext[1];
        double muEta1 = Ext[0]/muOmega;
        double muEta2 = Ext[2]/(1-muOmega);
        double VarOmega = Varxt.at(1,1) + Varxt.at(2,2) + 2 * Varxt.at(1,2);
        //double VarEta1 = (Varxt.at(1,1)- VarOmega * (muEta1 * muEta1))/(VarOmega + muOmega * muOmega);
        double VarEta1 = (Varxt.at(1,1)+Varxt.at(2,2)-
                (muEta1*muEta1+(1-muEta1)*(1-muEta1))*VarOmega)
                /(2*(muOmega*muOmega+VarOmega));
        double VarEta2 = (Varxt.at(3,3)+Varxt.at(4,4)-
                (muEta2*muEta2+(1-muEta2)*(1-muEta2))*VarOmega)
                /(2*(muOmega*muOmega+VarOmega));

        //double VarEta2 = (Varxt.at(3,3)- VarOmega * (muEta2 * muEta2))/(VarOmega + (1-muOmega) * (1-muOmega));


        if (VarOmega > muOmega*(1-muOmega)){
            VarOmega = muOmega*(1-muOmega) - 0.0000000000000001;
        }

        if (VarEta1 > muEta1*(1-muEta1)){
            VarEta1 = muEta1*(1-muEta1) - 0.0000000000000001;
        }

        if (VarEta2 > muEta2*(1-muEta2)){
            VarEta2 = muEta2*(1-muEta2) - 0.0000000000000001;
        }





        double phiOmega = muOmega * (1-muOmega)/VarOmega - 1;
        double phiEta1 = muEta1 * (1-muEta1)/VarEta1 - 1;
        double phiEta2 = muEta2 * (1-muEta2)/VarEta2 - 1;



        double alphaAG = muOmega*phiOmega;
        double betaAG = phiOmega*(1-muOmega);
        double alphaA = muEta1*phiEta1;
        double betaA = phiEta1*(1-muEta1);
        double alphaC = muEta2*phiEta2;
        double betaC = phiEta2*(1-muEta2);
         */



        /*/
        if (muOmega*phiOmega > 0 ){
 alphaAG = muOmega*phiOmega;
 } else{
 alphaAG = 0.0000000000000001;
 }

 if (phiOmega*(1-muOmega) > 0 ){
 betaAG = phiOmega*(1-muOmega);
 } else{
 betaAG = 0.0000000000000001;
 }

 if (muEta1*phiEta1 > 0 ){
 alphaA = muEta1*phiEta1;
 } else{
 alphaA = 0.0000000000000001;
 }

 if (phiEta1*(1-muEta1) > 0 ){
 betaA = phiEta1*(1-muEta1);
 } else{
 betaA = 0.0000000000000001;
 }

 if (muEta2*phiEta2 > 0 ){
 alphaC = muEta2*phiEta2;
 } else{
 alphaC = 0.0000000000000001;
 }

 if (phiEta2*(1-muEta2) > 0 ){
 betaC = phiEta2*(1-muEta2);
 } else{
 betaC = 0.0000000000000001;
 }

         */













        double component1 = alleleFreq[0] + alleleFreq[1];
        double component2 = alleleFreq[0]/(alleleFreq[0] + alleleFreq[1]);
        double component3 = alleleFreq[2]/(alleleFreq[2] + alleleFreq[3]);

        if (alleleFreq[0] + alleleFreq[1] == 1){
            double n = unNormalizedDensity(component2,alphaA,betaA) / BetaF_A;
            double a = normailzeConstant1 * n;
            if (Double.isNaN(a)){
                if (Double.isNaN(n)){
                    if (component2 == 0.5 && !Double.isNaN(normailzeConstant1) && normailzeConstant1 != 0){
                        return 1000000;
                    } else {
                        return 0;
                    }
                }
                return 0;
            } else {
                if (Double.isNaN(a)){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!");
                }
                return a;
            }
//            if (Double.isNaN(a)){
//                System.out.println("*******************************");
//                System.out.println("1");
//                System.out.println(Arrays.toString(Ext));
//                System.out.println(Varxt);
//                System.out.println("muOmega = "+muOmega);
//                System.out.println("muEta1 = "+muEta1);
//                System.out.println("muEta2 = "+muEta2);
//                System.out.println("VarOmega = "+VarOmega);
//                System.out.println("VarEta1 = "+VarEta1);
//                System.out.println("VarEta2 = "+VarEta2);
//                System.out.println("phiOmega = "+phiOmega);
//                System.out.println("phiEta1 = "+phiEta1);
//                System.out.println("phiEta2 = "+phiEta2);
//                System.out.println("alphaAG = " + alphaAG);
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("alphaA = " + alphaA);
//                System.out.println("betaA = " + betaA);
//                System.out.println("alphaC = " + alphaC);
//                System.out.println("betaC = " + betaC);
//                //System.out.println("component1 = "+ component1);
//                System.out.println("component2 = "+ component2);
//                System.out.println("The formula:" + "Math.pow(N,-betaAG) / (betaAG * BetaF(alphaAG,betaAG)) * distEta1.density(component2)");
//                System.out.println("Math.pow(N,-betaAG) = " + Math.pow(N,-betaAG));
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("BetaF(alphaAG,betaAG) = " + BetaF(alphaAG,betaAG));
//                System.out.println("distEta1.density(component2) = " + unNormalizedDensity(component2,alphaA,betaA));
//                //System.out.println("component3 = " + component3);
//                System.out.println("*******************************");
//                return 0;
//            }
//            return a;
        } else if (alleleFreq[0] + alleleFreq[2] == 1) {
            double n = unNormalizedDensity(component1,alphaAG,betaAG)/ BetaF_AG;
            double a = 1 / (alleleFreq[0] * alleleFreq[2]) *  normailzeConstant2 * n;
            if (Double.isNaN(a)){
                if (Double.isNaN(n)){
                    if (component1 == 0.5 && !Double.isNaN(normailzeConstant2) && normailzeConstant2 != 0){
                        return 1000000;
                    } else {
                        return 0;
                    }
                }
                return 0;
            } else {
                if (Double.isNaN(a)){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!");
                }
                return a;
            }
//            if (Double.isNaN(a)){
//                System.out.println("*******************************");
//                System.out.println("2");
//                System.out.println(Arrays.toString(Ext));
//                System.out.println(Varxt);
//                System.out.println("muOmega = "+muOmega);
//                System.out.println("muEta1 = "+muEta1);
//                System.out.println("muEta2 = "+muEta2);
//                System.out.println("VarOmega = "+VarOmega);
//                System.out.println("VarEta1 = "+VarEta1);
//                System.out.println("VarEta2 = "+VarEta2);
//                System.out.println("phiOmega = "+phiOmega);
//                System.out.println("phiEta1 = "+phiEta1);
//                System.out.println("phiEta2 = "+phiEta2);
//                System.out.println("alphaAG = " + alphaAG);
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("alphaA = " + alphaA);
//                System.out.println("betaA = " + betaA);
//                System.out.println("alphaC = " + alphaC);
//                System.out.println("betaC = " + betaC);
//                System.out.println("component1 = "+ component1);
//                //System.out.println("component2 = "+ component2);
//                //System.out.println("component3 = " + component3);
//                System.out.println("The formula:"+"1 / (alleleFreq[0] * alleleFreq[2])\n" +
//                        "                    * distOmega.density(component1)\n" +
//                        "                    * Math.pow(N,-betaA) / (betaA * BetaF(alphaA,betaA))\n" +
//                        "                    * Math.pow(N,-betaC) / (betaC * BetaF(alphaC,betaC))");
//                System.out.println("distOmega.density(component1) = " + unNormalizedDensity(component1,alphaAG,betaAG) );
//                System.out.println("Math.pow(N,-betaA) = " + Math.pow(N,-betaA));
//                System.out.println("betaA * BetaF(alphaA,betaA) = " + betaA * BetaF(alphaA,betaA));
//                System.out.println("Math.pow(N,-betaC) = " + Math.pow(N,-betaC));
//                System.out.println("betaC * BetaF(alphaC,betaC) = " + betaC * BetaF(alphaC,betaC));
//                System.out.println("*******************************");
//                return 0;
//            }
//            return a;

        } else if (alleleFreq[0] + alleleFreq[3] == 1) {
            double n = unNormalizedDensity(component1,alphaAG,betaAG) / BetaF_AG;
            double a = 1 / (alleleFreq[0] * alleleFreq[3]) * normailzeConstant3 * n;
            if (Double.isNaN(a)){
                if (Double.isNaN(n)){
                    if (component1 == 0.5&& !Double.isNaN(normailzeConstant3) && normailzeConstant3 != 0){
                        return 1000000;
                    } else {
                        return 0;
                    }
                }
                return 0;
            } else {
                if (Double.isNaN(a)){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!");
                }
                return a;
            }
//            if (Double.isNaN(a)){
//                System.out.println("*******************************");
//                System.out.println("3");
//                System.out.println(Arrays.toString(Ext));
//                System.out.println(Varxt);
//                System.out.println("muOmega = "+muOmega);
//                System.out.println("muEta1 = "+muEta1);
//                System.out.println("muEta2 = "+muEta2);
//                System.out.println("VarOmega = "+VarOmega);
//                System.out.println("VarEta1 = "+VarEta1);
//                System.out.println("VarEta2 = "+VarEta2);
//                System.out.println("phiOmega = "+phiOmega);
//                System.out.println("phiEta1 = "+phiEta1);
//                System.out.println("phiEta2 = "+phiEta2);
//                System.out.println("alphaAG = " + alphaAG);
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("alphaA = " + alphaA);
//                System.out.println("betaA = " + betaA);
//                System.out.println("alphaC = " + alphaC);
//                System.out.println("betaC = " + betaC);
//                System.out.println("component1 = "+ component1);
//                System.out.println("The formula:" + "1 / (alleleFreq[0] * alleleFreq[3])\n" +
//                        "                    * distOmega.density(component1)\n" +
//                        "                    * Math.pow(N,-betaA) / (betaA * BetaF(alphaA,betaA))\n" +
//                        "                    * Math.pow(N,-alphaC) / (alphaC * BetaF(alphaC,betaC))");
//                System.out.println("distOmega.density(component1) = " + unNormalizedDensity(component1,alphaAG,betaAG));
//                System.out.println("Math.pow(N,-betaA) = " + Math.pow(N,-betaA));
//                System.out.println("betaA * BetaF(alphaA,betaA) = " + betaA * BetaF(alphaA,betaA));
//                System.out.println("Math.pow(N,-alphaC) = " + Math.pow(N,-alphaC));
//                System.out.println("alphaC * BetaF(alphaC,betaC) = " + alphaC * BetaF(alphaC,betaC));
//                //System.out.println("component2 = "+ component2);
//                //System.out.println("component3 = " + component3);
//                System.out.println("*******************************");
//                return 0;
//            }
//            return a;
        } else if (alleleFreq[1] + alleleFreq[2] == 1) {
            double n = unNormalizedDensity(component1,alphaAG,betaAG)/ BetaF_AG;
            double a = 1 / (alleleFreq[1] * alleleFreq[2]) * normailzeConstant4 * n;
            if (Double.isNaN(a)){
                if (Double.isNaN(n)){
                    if (component1 == 0.5&& !Double.isNaN(normailzeConstant4) && normailzeConstant4 != 0){
                        return 1000000;
                    } else {
                        return 0;
                    }
                }
                return 0;
            } else {
                if (Double.isNaN(a)){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!");
                }
                return a;
            }
//            if (Double.isNaN(a)){
//                System.out.println("*******************************");
//                System.out.println("4");
//                System.out.println(Arrays.toString(Ext));
//                System.out.println(Varxt);
//                System.out.println("muOmega = "+muOmega);
//                System.out.println("muEta1 = "+muEta1);
//                System.out.println("muEta2 = "+muEta2);
//                System.out.println("VarOmega = "+VarOmega);
//                System.out.println("VarEta1 = "+VarEta1);
//                System.out.println("VarEta2 = "+VarEta2);
//                System.out.println("phiOmega = "+phiOmega);
//                System.out.println("phiEta1 = "+phiEta1);
//                System.out.println("phiEta2 = "+phiEta2);
//                System.out.println("alphaAG = " + alphaAG);
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("alphaA = " + alphaA);
//                System.out.println("betaA = " + betaA);
//                System.out.println("alphaC = " + alphaC);
//                System.out.println("betaC = " + betaC);
//                System.out.println("component1 = "+ component1);
//                System.out.println("The formula:" + "1 / (alleleFreq[1] * alleleFreq[2])\n" +
//                        "                    * distOmega.density(component1)\n" +
//                        "                    * Math.pow(N,-alphaA) / (alphaA * BetaF(alphaA,betaA))\n" +
//                        "                    * Math.pow(N,-betaC) / (betaC * BetaF(alphaC,betaC))");
//                System.out.println("distOmega.density(component1) = " +unNormalizedDensity(component1,alphaAG,betaAG));
//                System.out.println("Math.pow(N,-alphaA) = " + Math.pow(N,-alphaA) );
//                System.out.println("alphaA * BetaF(alphaA,betaA) = " + alphaA * BetaF(alphaA,betaA));
//                System.out.println("Math.pow(N,-betaC) = " + Math.pow(N,-betaC));
//                System.out.println("betaC * BetaF(alphaC,betaC) = " + betaC * BetaF(alphaC,betaC));
//                //System.out.println("component2 = "+ component2);
//                //System.out.println("component3 = " + component3);
//                System.out.println("*******************************");
//                return 0;
//            }
//            return a;
        } else if (alleleFreq[1] + alleleFreq[3] == 1) {

            double n = unNormalizedDensity(component1,alphaAG,betaAG)/ BetaF_AG;
            double a = 1 / (alleleFreq[1] * alleleFreq[3]) * normailzeConstant5 * n;

            if (Double.isNaN(a)){
                if (Double.isNaN(n)){
                    if (component1 == 0.5&& !Double.isNaN(normailzeConstant5) && normailzeConstant5 != 0){
                        return 1000000;
                    } else {
                        return 0;
                    }
                }
                return 0;
            } else {
                if (Double.isNaN(a)){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!");
                }
                return a;
            }
//            if (Double.isNaN(a)){
//                System.out.println("*******************************");
//                System.out.println("5");
//                System.out.println(Arrays.toString(Ext));
//                System.out.println(Varxt);
//                System.out.println("muOmega = "+muOmega);
//                System.out.println("muEta1 = "+muEta1);
//                System.out.println("muEta2 = "+muEta2);
//                System.out.println("VarOmega = "+VarOmega);
//                System.out.println("VarEta1 = "+VarEta1);
//                System.out.println("VarEta2 = "+VarEta2);
//                System.out.println("phiOmega = "+phiOmega);
//                System.out.println("phiEta1 = "+phiEta1);
//                System.out.println("phiEta2 = "+phiEta2);
//                System.out.println("alphaAG = " + alphaAG);
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("alphaA = " + alphaA);
//                System.out.println("betaA = " + betaA);
//                System.out.println("alphaC = " + alphaC);
//                System.out.println("betaC = " + betaC);
//                System.out.println("component1 = "+ component1);
//                System.out.println("The formula" +"1 / (alleleFreq[1] * alleleFreq[3])\n" +
//                        "                    * distOmega.density(component1)\n" +
//                        "                    * Math.pow(N,-alphaA) / (alphaA * BetaF(alphaA,betaA))\n" +
//                        "                    * Math.pow(N,-alphaC) / (alphaC * BetaF(alphaC,betaC))");
//                System.out.println("distOmega.density(component1) = " + unNormalizedDensity(component1,alphaAG,betaAG));
//                System.out.println("Math.pow(N,-alphaA) = " + Math.pow(N,-alphaA));
//                System.out.println("alphaA * BetaF(alphaA,betaA) = " + alphaA * BetaF(alphaA,betaA));
//                System.out.println("Math.pow(N,-alphaC) = " + Math.pow(N,-alphaC));
//                System.out.println("alphaC * BetaF(alphaC,betaC) = " + alphaC * BetaF(alphaC,betaC));
//                //System.out.println("component2 = "+ component2);
//                //System.out.println("component3 = " + component3);
//                System.out.println("*******************************");
//                return 0;
//            }
//            return a;

        } else if (alleleFreq[2] + alleleFreq[3] == 1) {
            double n = unNormalizedDensity(component3, alphaC,betaC) / BetaF(alphaC,betaC);
            double a = normailzeConstant6 * n;
            if (Double.isNaN(a)){
                if (Double.isNaN(n)){
                    if (component3 == 0.5&& !Double.isNaN(normailzeConstant6) && normailzeConstant6 != 0){
                        return 1000000;
                    } else {
                        return 0;
                    }
                }
                return 0;
            } else {
                if (Double.isNaN(a)){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!");
                }
                return a;
            }
//            if (Double.isNaN(a)){
//                System.out.println("*******************************");
//                System.out.println("6");
//                System.out.println(Arrays.toString(Ext));
//                System.out.println(Varxt);
//                System.out.println("muOmega = "+muOmega);
//                System.out.println("muEta1 = "+muEta1);
//                System.out.println("muEta2 = "+muEta2);
//                System.out.println("VarOmega = "+VarOmega);
//                System.out.println("VarEta1 = "+VarEta1);
//                System.out.println("VarEta2 = "+VarEta2);
//                System.out.println("phiOmega = "+phiOmega);
//                System.out.println("phiEta1 = "+phiEta1);
//                System.out.println("phiEta2 = "+phiEta2);
//                System.out.println("alphaAG = " + alphaAG);
//                System.out.println("betaAG = " + betaAG);
//                System.out.println("alphaA = " + alphaA);
//                System.out.println("betaA = " + betaA);
//                System.out.println("alphaC = " + alphaC);
//                System.out.println("betaC = " + betaC);
//                System.out.println("The formula" + "Math.pow(N,-betaAG) / (alphaAG * BetaF(alphaAG,betaAG)) * distEta2.density(component3)");
//                System.out.println("Math.pow(N,-betaAG) = " + Math.pow(N,-betaAG));
//                System.out.println("alphaAG * BetaF(alphaAG,betaAG) = " + alphaAG * BetaF(alphaAG,betaAG));
//                System.out.println("distEta2.density(component3) = " +unNormalizedDensity(component3, alphaC,betaC));
//                //System.out.println("component1 = "+ component1);
//                //System.out.println("component2 = "+ component2);
//                System.out.println("component3 = " + component3);
//                System.out.println("*******************************");
//                return 0;
//            }
//            return a;

        } else{
            return 0;
        }





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
