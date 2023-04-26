package STAM;


import beast.core.Description;

@Description("Give the mean and covariance matrix required for moment approximation method")

public class STAMMomentApproximation {



    // Initial allele frequence
    double[] x0;

    // Time
    double t;


    // Rate matrix
    Array2d Q;

    // Number of the alleles
    int dim = 4;

    // Expectation
    double[] Ext;

    // Covariance matrix
    Array2d Varxt;


    // Constructor for a moment result object
    public STAMMomentApproximation(double[] x0, double t, Array2d Q){
        this.x0 = x0;
        this.t = t;
        this.Q = Q;
    }

    public void Kimura(){
        //long startTime=System.nanoTime();
        double[] matrix = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        Array2d E = new Array2d(dim,dim,matrix);
        double alpha = Q.at(1,2);
        //System.out.println(alpha);
        double beta = Q.at(1,3);
        //System.out.println(beta);
        double[] a = new double[]{0.25,0.25,0.5};
        double[] b = new double[]{0,4*beta,2*(alpha+beta)};
        Array2d A1 = E;
        Array2d A1Copy = E.copy();
        Array2d A2 = new Array2d(dim,dim,new double[]{1,1,-1,-1,1,1,-1,-1,-1,-1,1,1,-1,-1,1,1});
        Array2d A2Copy = A2.copy();
        Array2d A3 = new Array2d(dim,dim,new double[]{1,-1,0,0,-1,1,0,0,0,0,1,-1,0,0,-1,1});
        Array2d A3Copy = A3.copy();
        Array2d[] A = new Array2d[]{A1,A2,A3};




        A1Copy.scalarMul(a[0]);
        A1Copy.scalarMul(Math.exp(-b[0]*t));

        //System.out.println(A1Copy);
        A2Copy.scalarMul(a[1]);
        A2Copy.scalarMul(Math.exp(-b[1]*t));
        //System.out.println(A2Copy);
        A3Copy.scalarMul(a[2]);
        A3Copy.scalarMul(Math.exp(-b[2]*t));
        //System.out.println(A3Copy);
        A1Copy.addMatrix(A2Copy);
        A1Copy.addMatrix(A3Copy);

        Array2d eQt = A1Copy;
        //System.out.println(eQt);

        Ext = A1Copy.mulrowVectorLeft(x0);

        Array2d sum = new Array2d(4,4,new double[]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});


        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                for (int k = 0 ; k < 3; k++){
                    Array2d ACopy = A[i].copy();
                    ACopy.scalarMul(a[i] * a[j] * a[k] * gijk(i,j,k,b));
                    Array2d Diag = diag(A[j].mulrowVectorLeft(x0));
                    ACopy.multiplyRight(Diag);
                    ACopy.multiplyRight(A[k]);
                    sum.addMatrix(ACopy);
                }
            }
        }

        //long endTime=System.nanoTime(); //获取结束时间
        //System.out.println("程序运行时间 （计算Qmatrix）： "+(endTime-startTime)+"ns");
        Array2d eQtT = eQt.transpose();
        eQtT.mulcolVectorRight(x0);
        Array2d resultVec = vecMul(eQtT.mulcolVectorRight(x0),x0);
        //System.out.println(resultVec);
        resultVec.multiplyRight(eQt);
        resultVec.scalarMul(Math.exp(-t)-1);
        sum.addMatrix(resultVec);
        Varxt = sum;


        // For testing
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
                /(2*((1-muOmega)*(1-muOmega)+VarOmega));

        //double VarEta2 = (Varxt.at(3,3)- VarOmega * (muEta2 * muEta2))/(VarOmega + (1-muOmega) * (1-muOmega));






        double phiOmega = muOmega * (1-muOmega)/VarOmega - 1;
        double phiEta1 = muEta1 * (1-muEta1)/VarEta1 - 1;
        double phiEta2 = muEta2 * (1-muEta2)/VarEta2 - 1;



        double alphaAG = muOmega*phiOmega;
        double betaAG = phiOmega*(1-muOmega);
        double alphaA = muEta1*phiEta1;
        double betaA = phiEta1*(1-muEta1);
        double alphaC = muEta2*phiEta2;
        double betaC = phiEta2*(1-muEta2);


        if (alphaAG <= 0 || betaAG <= 0){
            System.out.println("Ext = "+Arrays.toString(Ext));
            System.out.println("Varxt = ");
            System.out.println(Varxt);
            System.out.println("t = " + t);
            System.out.println("x0 = " + Arrays.toString(x0));
            System.out.println("Q = " + Q);
            System.out.println("alphaAG or betaAG become negative");
            System.out.println("alphaAG = " + alphaAG);
            System.out.println("betaAG = " + betaAG);
            System.out.println("***************************************************");

        }

        if (alphaA <= 0 || betaA <= 0){
            System.out.println("Ext = "+Arrays.toString(Ext));
            System.out.println("Varxt = ");
            System.out.println(Varxt);
            System.out.println("t = " + t);
            System.out.println("x0 = " + Arrays.toString(x0));
            System.out.println("Q = " + Q);
            System.out.println("alphaA or betaA become negative");
            System.out.println("alphaA = " + alphaA);
            System.out.println("betaA = " + betaA);
            System.out.println("***************************************************");
        }

        if (alphaC <= 0 || betaC <= 0){
            System.out.println("Ext = "+Arrays.toString(Ext));
            System.out.println("Varxt = ");
            System.out.println(Varxt);
            System.out.println("t = " + t);
            System.out.println("x0 = " + Arrays.toString(x0));
            System.out.println("Q = " + Q);
            System.out.println("alphaC or betaC become negative");
            System.out.println("alphaC = " + alphaC);
            System.out.println("betaAC = " + betaC);
            System.out.println("***************************************************");
        }
         */

        // For testing

    }

    // help function for Kimura Model
    public double gijk(int i, int j, int k, double[] b){
        if (1 + b[i] + b[k] == b[j]){
            return t*Math.exp(-b[j]*t);
        } else{
            return (Math.exp(-b[j]*t) - Math.exp(-(1+b[i]+b[k])*t))/(1+b[i]+b[k]-b[j]);
        }
    }

    // help function for Kimura Model (subject to 4 dimensional)

    public Array2d diag(double[] diags){
        double[] result = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        result[0] = diags[0];
        result[5] = diags[1];
        result[10] = diags[2];
        result[15] = diags[3];
        return new Array2d(4,4,result);
    }

    // help function for Kimura (subject to 4 dimensional)
    public Array2d vecMul(double[] v1, double[] v2){
        Array2d matrix = new Array2d(4,4);
        for (int j = 1; j < 5; j++){
            for (int i = 1; i < 5; i++){
                matrix.set(i,j,v1[i-1]*v2[j-1]);
            }
        }
        return matrix;
    }




}
