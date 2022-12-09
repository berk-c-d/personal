import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.stream.IntStream;


public class cost_version3_esra1_m {

    public static void main(String[] args) throws FileNotFoundException {
        Instant inst1 = Instant.now();

        int BigM=10000; int nresult =42;
        int period =5;
        double[] B0 = new double[period];
        double[] B1 = new double[period];
        double[] B2 = new double[period];
        double[] B3 = new double[period];
        int inv = 0;
        int j = 0;
        int pp = 0;


        //Base Case
        int unitCost = 40;

        double holdingCost = 2;

        double backorderCost = 20;

        int salvageValue = 0;

        int fixedCost = 25;

        int returnCredit = 80;

        int retailPrice = 80;

        double alpha = 0.5; //return rate

        double lambda = 4; //demand rate

        double alphaSimulation = 0.5;

        double lambdaSimulation = 4;

        int m = 10; //max demand

        //Intercept
        B0[0]=0;
        B0[1]=643.800653902418
        ;
        B0[2]=440.40165902393
        ;
        B0[3]=270.931722295635
        ;
        B0[4]=145.461518383893
        ;

        //B1_Previous_sale
        B1[0]=0;
        B1[1]=42.6431417173849
        ;
        B1[2]=42.1093295446437
        ;
        B1[3]=42.1093295446437
        ;
        B1[4]=40.3002076368292
        ;

        //B2_Inv_min_brk_pos
        B2[0]=0;
        B2[1]=1.4239134152982
        ;
        B2[2]=0.472005889651649
        ;
        B2[3]=-0.260388064701157
        ;
        B2[4]=0.539082572245109
        ;


        //B3_Prev_mult_Inv_min_brk_neg
        B3[0]=0;
        B3[1]=-37.5725734746556
        ;
        B3[2]=-51.5009157588997
        ;
        B3[3]=-62.9101891728928
        ;
        B3[4]=-74.9655721635784
        ;


        //--------------Inventory level after ordering---------------------------------------------------------------------------------
        IntStream streamS = IntStream.range(-((period-1) * m), ((period)*m)+1);;
        int MM = (period-1) * m; //neutralizes negative S
        int[] S = streamS.toArray();
        int s = S.length;
        //------------------------------------------------------------------------------------------------------------------------------

        int sl =((period+2) * m)+1; //max # of sales

        double[] poissonDistProbability = new double[m+1];
        poissonDistProbability = poissonDist (lambda, m);
        double[][] binomialDistProbability = new double[sl][sl+1];
        binomialDistProbability = binomialDist(sl, alpha);

        double[] arr = new double[nresult];
        arr = SimulationMatrix_v3adapted(period, m, alphaSimulation, lambdaSimulation, unitCost, fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, salvageValue,nresult, inv , j,S,pp, sl, s, alpha, lambda, BigM, MM,  poissonDistProbability, binomialDistProbability,B0,B1,B2,B3);

        System.out.println(Arrays.toString(arr));

        Instant inst2 = Instant.now();
        System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
    }

    public static double[] SimulationMatrix_v3adapted(int period, int m, double alphaSimulation, double lambdaSimulation, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost, int salvageValue,int nresult, int inv , int j,int [] S,int pp,int sl, int s, double alpha, double lambda, int BigM, int MM,  double[] poissonDistProbability, double[][] binomialDistProbability,double[] B0,double[] B1,double[] B2,double[] B3) {

        int nn = 100000;

        double[] arr = new double[nresult];
        double[] arr_sum = new double[nresult];
        double[] arr_mean = new double[nresult];


        for (int i = 0; i < nn; i++) {
            arr = ProcessMatrix_v3adapted(period, m, alphaSimulation, lambdaSimulation, unitCost, fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, salvageValue, nresult, inv , j, S, pp,sl, s, alpha, lambda, BigM, MM, poissonDistProbability, binomialDistProbability, B0,B1,B2,B3);

            for(int jj=0;jj<arr.length;jj++){
                arr_sum[jj]+=arr[jj];
            }
        }

        for(int jj=0;jj<arr.length;jj++){
            arr_mean[jj]=arr_sum[jj]/nn;
        }
        return arr_mean;

    }

    public static double[] ProcessMatrix_v3adapted(int period, int m, double alphaSimulation, double lambdaSimulation, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost, int salvageValue, int nresult, int inv , int j,int [] S,int pp,int sl, int s, double alpha, double lambda, int BigM, int MM,  double[] poissonDistProbability, double[][] binomialDistProbability,double[] B0,double[] B1,double[] B2,double[] B3) {


        double[] result = new double[nresult];

        int[] I = new int[period + 1]; //inventory
        int[] R = new int[period + 1]; //return
        int[] O = new int[period]; //order
        int[] D = new int[period]; //demand
        int[] BO = new int[period]; //backorder
        int[] IO = new int[period + 1]; //inventory on hand
        int[] SL = new int[period + 1]; //sales
        int[] HI = new int[period]; // holding inventory


        int[] numberOfOrder = new int[period];


        int[] Z_bigSDouble = new int[period];

        for (int p = 0; p < period; p++) {
            D[p] = getPoissonRandom(lambdaSimulation);
            if (D[p]>m) {
                D[p]=m;
            }
        }

        I[0] = 0;
        R[0] = 0;

        Z_bigSDouble[0]  =  periodMatrixmodifiedv6( period, s,  m, alpha, lambda, unitCost, fixedCost,  returnCredit, holdingCost, backorderCost,  poissonDistProbability, B0, B1, B2, B3);

        O[0] = Z_bigSDouble[0];

        BO[0] = Math.max(0, D[0] - I[0] - O[0] - R[0]);
        IO[0] = Math.max(0, I[0]);
        SL[0] = Math.min(IO[0] + O[0], D[0]);
        HI[0] = Math.max(0, I[0] + O[0] + R[0] - D[0]);


        for (int i = 1; i < period-1; i++) {

            I[i] = I[i - 1] + O[i - 1] + R[i - 1] - D[i - 1];
            R[i] = getBinomialRandom(SL[i - 1], alphaSimulation);

            Z_bigSDouble[i]  = periodMatrixmodifiedv4(I[i],SL[i-1], S, period, i, sl, s, m, alpha, lambda, unitCost, fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, BigM, MM, poissonDistProbability, binomialDistProbability, B0, B1,B2,B3);


            O[i] = Math.max(0, (Z_bigSDouble[i]- I[i]));
            BO[i] = Math.max(0, D[i] - I[i] - O[i] - R[i]);
            IO[i] = Math.max(0, I[i]);
            SL[i] = Math.min(IO[i] + O[i] + R[i], D[i] + BO[i - 1]);
            HI[i] = Math.max(0, I[i] + O[i] + R[i] - D[i]);

        }
        I[period-1] = I[period-2] + O[period-2] + R[period-2] - D[period-2];
        R[period-1] = getBinomialRandom(SL[period-2], alphaSimulation);

        Z_bigSDouble[period-1]  = periodMatrixmodifiedv5(I[period-1],SL[period-2], S, period, period-1,sl, s, m, alpha, lambda, unitCost, fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, salvageValue,BigM, MM, poissonDistProbability, binomialDistProbability);


        O[period-1] = Math.max(0, (Z_bigSDouble[period-1]- I[period-1]));


        BO[period-1] = Math.max(0, D[period-1] - I[period-1] - O[period-1] - R[period-1]);
        IO[period-1] = Math.max(0, I[period-1]);
        SL[period-1] = Math.min(IO[period-1] + O[period-1] + R[period-1], D[period-1] + BO[period-2]);
        HI[period-1] = Math.max(0, I[period-1] + O[period-1] + R[period-1] - D[period-1]);


        I[period] = I[period - 1] + O[period - 1] + R[period - 1] - D[period - 1];
        R[period] = getBinomialRandom(SL[period - 1], alphaSimulation);

        double procurementCost = 0;
        double totalOrderingCost = 0;
        double totalHoldingCost = 0;
        double totalBackorderCost = 0;
        double returnCost = 0;
        double returnCost1 = 0;
        double totalEndingInventoryCost = 0;
        double salvage=0;
        int sumNumberOfOrder =0;


        for (int p = 0; p < period; p++) {

            procurementCost += O[p] * unitCost;

            if (O[p] > 0) {
                numberOfOrder[p] = 1;
            } else {
                numberOfOrder[p] = 0;
            }

            totalOrderingCost += numberOfOrder[p] * fixedCost;

            sumNumberOfOrder += numberOfOrder[p] ;

            totalHoldingCost += HI[p] * holdingCost;

            totalBackorderCost += BO[p] * backorderCost;

            if (I[period] + R[period] < 0) {
                totalEndingInventoryCost = (I[period] + R[period]) * -(unitCost)
                        - I[period] * alphaSimulation * (returnCredit - salvageValue);
            } else if (I[period] < 0) {
                totalEndingInventoryCost = (I[period] + R[period]) * -salvageValue
                        - I[period] * alphaSimulation * (returnCredit - salvageValue);
            } else {
                totalEndingInventoryCost = (I[period] + R[period]) * -salvageValue;
            }

            if (I[period] < 0) {
                SL[period] = -(I[period]);
            } else {
                SL[period] = 0;
            }

            returnCost1 += R[p] * returnCredit;
            returnCost = returnCost1 + R[period] * returnCredit;

            salvage=I[period]+R[period];


        }

////////////////////////////////////////////////////////////////////////


        double pc = procurementCost;
        double oc = totalOrderingCost;
        double hc = totalHoldingCost;
        double bc = totalBackorderCost;
        double tec = totalEndingInventoryCost;
        double rc = returnCost;
        double SLV= salvage;


        double totalCost = procurementCost + totalOrderingCost + totalHoldingCost + totalBackorderCost + totalEndingInventoryCost + returnCost;

        double tc = totalCost;


        result = new double[] {tc, pc, oc, hc, bc,tec, rc,  SLV,tc, O[0],O[1],O[2],O[3],O[4], Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3],Z_bigSDouble[4], BO[0],BO[1],BO[2],BO[3],BO[4], HI[0],HI[1],HI[2],HI[3],HI[4], R[1],R[2],R[3],R[4],R[5], SL[0],SL[1],SL[2],SL[3],SL[4],SL[5], I[5]  };


        //   result = new double[] {tc};

        return result;

    }

    public static int periodMatrixmodifiedv4(int inv , int j, int [] S,int period, int pp, int sl, int s, int m, double alpha, double lambda, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost,int BigM, int MM,   double[] poissonDistProbability, double[][] binomialDistProbability,double[] B0,double[] B1,double[] B2,double[] B3) {

        double total;
        double total2;
        int Z_bigSDouble = 0;
        int SS = 0;
        double z_line;
        //----------------------------------------

        double min = 100000;

        for (int i = 0; i < s; i++) { //S up to order
            total2 = 0;
            for (int jj = 0; jj <= m; jj++) { //demand
                total = 0;
                for (int k = 0; k <= j; k++) { //return
                    double prev_sale=Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv));
                    double inventory=S[i] + k - jj;
                    double break_point=lambda*(period-pp-1)-alpha*prev_sale;
                    double break_point_pos=Math.max(break_point,0);
                    double inv_min_brk_point=inventory-break_point_pos;
                    double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                    double inv_min_brk_point_neg=Math.min(0,inv_min_brk_point);
                    double reg_result=B0[pp+1] +B1[pp+1]*prev_sale +B2[pp+1]*inv_min_brk_point_pos +B3[pp+1]*inv_min_brk_point_neg*prev_sale;
                    //System.out.println("pp+1: "+(pp+1)+" Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);


                    if (inv < S[i]) {
                        z_line = (binomialDistProbability[j][k]) *  (poissonDistProbability[jj]) * ((S[i]- inv) * unitCost
                                + fixedCost +  Math.max(0,(S[i] + k - jj))* holdingCost
                                + Math.max(0,(jj - S[i]  - k))* backorderCost + k * returnCredit
                                +reg_result);



                        //+ B1[pp+1]*Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv))
                        //+ B2[pp+1]*Math.min(0,(S[i] + k - jj-Math.max(0,(lambda*(period-pp)-alpha*Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv))))))
                        //+ B3[pp+1]*Math.max(0,(S[i] + k - jj-Math.max(0,(lambda*(period-pp)-alpha*Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv)))))));
                    } else if (inv ==S[i]) {
                        z_line = (binomialDistProbability[j][k]) * (poissonDistProbability[jj]) *(Math.max(0,(S[i] + k - jj))* holdingCost
                                + Math.max(0,(jj - S[i]  - k))* backorderCost + k * returnCredit
                                +reg_result);

                        //+ B1[pp+1]*Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv))
                        //+ B2[pp+1]*Math.min(0,(S[i] + k - jj-Math.max(0,(lambda*(period-pp)-alpha*Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv))))))
                        //+ B3[pp+1]*Math.max(0,(S[i] + k - jj-Math.max(0,(lambda*(period-pp)-alpha*Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv)))))));
                    } else {
                        z_line = BigM;
                    }
                    total += z_line;

                }
                total2 +=total;
                //   System.out.println("Total: " +  i+"  "+jj+"  "+total);
            }

            //     System.out.println("Total: " +  i+"  "+"  "+total2);
            if (total2< min) {
                min = total2;
                SS = i;
            }
            //     System.out.println("Total: " +  i+"  "+"jj"+"  "+min);
        }

        Z_bigSDouble = SS - MM;

        //   System.out.println("Total: " +  1+"  "+"jj"+"  "+ Z_bigSDouble);


        return Z_bigSDouble;
    }

    public static int periodMatrixmodifiedv5(int inv , int j, int [] S,int period, int pp , int sl, int s, int m, double alpha, double lambda, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost,int salvageValue, int BigM, int MM,  double[] poissonDistProbability, double[][] binomialDistProbability) {

        double total;
        double total2;
        int Z_bigSDouble = 0;
        int SS = 0;
        double z_line;

        //----------------------------------------

        double min = 100000;

        for (int i = 0; i < s; i++) { //S up to order
            total2 = 0;
            for (int jj = 0; jj <= m; jj++) { //demand
                total = 0;
                for (int k = 0; k <= j; k++) { //return
                    if (inv < S[i]) {
                        z_line = (binomialDistProbability[j][k]) *  (poissonDistProbability[jj]) * ((S[i]- inv) * unitCost
                                + fixedCost +  (Math.max(0,(S[i] + k - jj)))* holdingCost
                                +(Math.max(0,(jj - S[i]  - k)))* backorderCost + k * returnCredit
                                + Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv)) *alpha * returnCredit+ ( Math.min(0, inv+k)*-unitCost )+ ( Math.min(0, inv+k)* alpha * (returnCredit - salvageValue) )
                                +( Math.max(0, inv+k)*-salvageValue));
                    } else if (inv ==S[i]) {
                        z_line = (binomialDistProbability[j][k]) * (poissonDistProbability[jj]) *((Math.max(0,(S[i] + k - jj)))* holdingCost
                                +(Math.max(0,(jj - S[i]  - k)))* backorderCost + k * returnCredit
                                + Math.min(Math.max(0, inv) + Math.max(0, (S[i] - inv)) + k, jj - Math.min(0, inv))*alpha * returnCredit+ ( Math.min(0, inv+k)*-unitCost )+ ( Math.min(0, inv+k)* alpha * (returnCredit - salvageValue) )
                                +( Math.max(0, inv+k)*-salvageValue));
                    } else {
                        z_line = BigM;
                    }
                    total += z_line;
                    //    System.out.println("Total: " +  i+"  "+jj+"  "+total);
                }
                total2 +=total;
            }
            //    System.out.println("Total: " +  i+"  "+"  "+total2);
            if (total2< min) {
                min = total2;
                SS = i;
            }

        }

        Z_bigSDouble = SS - MM;


        return Z_bigSDouble;
    }

    public static int periodMatrixmodifiedv6(int period, int s, int m, double alpha, double lambda, int unitCost, int fixedCost,  int returnCredit, double holdingCost, double backorderCost,  double[] poissonDistProbability, double[] B0,double[] B1,double[] B2,double[] B3) {

        double total;
        int Z_bigSDouble = 0;
        int SS = 0;
        double z_line;
        //----------------------------------------


        double min = 100000;
        for (int i = 0; i < s; i++) { //S up to order
            total = 0;
            for (int jj = 0; jj <= m; jj++) { //demand
                double prev_sale=Math.min(i,jj);
                double break_point=lambda*(period-1)-alpha*Math.min(i,jj);
                double break_point_pos=Math.max(break_point,0);
                double inventory=i-jj;
                double inv_min_brk_point=inventory-break_point_pos;
                double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                double inv_min_brk_point_neg=Math.min(0,inv_min_brk_point);
                double reg_result=B0[1] +B1[1]*prev_sale +B2[1]*inv_min_brk_point_pos +B3[1]*inv_min_brk_point_neg*prev_sale;
                //System.out.println("Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);

                if (i>0) {
                    z_line = (poissonDistProbability[jj]) * (i * unitCost
                            + fixedCost +  (Math.max(0,(i - jj)))* holdingCost
                            +(Math.max(0,(jj - i)))* backorderCost
                            +reg_result);       //Upcoming period's cost

                }else {
                    z_line = (poissonDistProbability[jj]) * (
                            (Math.max(0,(i - jj)))* holdingCost
                                    +(Math.max(0,(jj - i)))* backorderCost
                                    + reg_result);

                }
                //System.out.println(z_line);
                total += z_line;

            }

            // System.out.println("Total: " +  i+"  "+"  "+total);
            if (total< min) {
                min = total;
                SS = i;
            }
        }

        Z_bigSDouble = SS;

        return Z_bigSDouble;
    }



    public static int getPoissonRandom(double lambda) {

        double L = Math.exp(-lambda);

        double p = 1.0;

        int k = 0;

        do {

            k++;

            p *= Math.random();

        } while (p > L);

        return k - 1;

    }

    public static int getBinomialRandom(int n, double alpha) {

        int x = 0;

        for (int i = 0; i < n; i++) {

            if (Math.random() < alpha)

                x++;

        }

        return x;

    }

    static double nCr(int n, int k) {

        if (k > n / 2)

            k = n - k;

        double answer = 1;

        for (int i = 1; i <= k; i++) {

            answer *= (n - k + i);

            answer /= i;

        }

        return answer;

    }

    private static double[][] binomialDist(int sl, double alpha) {

        double[][] P = new double[sl][sl];



        for (int j = 0; j <sl ; j++) {
            for (int k = 0; k <=j ; k++) {
                P[j][k]=nCr(j, k) * (double) Math.pow(alpha, k) * (double) Math.pow(1 - alpha, j - k);
            }
        }


        return P;

    }

    private static double[] poissonDist(double lambda,int m) {

        double[] P = new double[m+1];

        P[0]= Math.pow(Math.E, -1 * lambda);

        double Prest=0;

        for (int i = 1; i <m ; i++) {
            P[i]=P[i-1]*(lambda/i);
            Prest +=P[i];
        }

        P[m]=1-P[0]-Prest;

        return P;

    }


}
