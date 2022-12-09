import java.io.*;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.stream.IntStream;


public class cost_version1_berk {

    public static void main(String[] args) throws IOException {
        Instant inst1 = Instant.now();

        int BigM=10000; int nresult =42;
        int period =5;

        //Base Case
        /*
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
        int sl =((period+2) * m)+1; //max # of sales


        //-----------------------------------------------------
        int prev_sale_range=sl;
        int period_range=5;
        int num_of_betas=3;
        double[][][] beta_vals= new double[period_range][sl][num_of_betas];     //beta_vals[period][prevsale][0]=Intercept
                                                                                //beta_vals[period][prevsale][1]=InvInv min brk pos---
                                                                                //beta_vals[period][prevsale][2]=InvInv min brk neg---


        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results\\Regression Results\\Version_1\\Results for all periods_version1_for base case.csv";



         */
        //Case 5
        int case_number=5;
        double alpha_all=0.4;
        double holdingCost = 2;
        double backorderCost = 50;
        int fixedCost = 40;
        int m = 25; //max demand
        double lambda_all=8;

        int unitCost = 40;
        int salvageValue = 0;

        int returnCredit = 80;

        int retailPrice = 80;

        double alpha = alpha_all; //return rate

        double lambda = lambda_all; //demand rate

        double alphaSimulation = alpha_all;

        double lambdaSimulation = lambda_all;


        int sl =((period+2) * m)+1; //max # of sales


        //-----------------------------------------------------
        int prev_sale_range=sl;
        int period_range=5;
        int num_of_betas=3;
        double[][][] beta_vals= new double[period_range][sl][num_of_betas];     //beta_vals[period][prevsale][0]=Intercept
                                                                                //beta_vals[period][prevsale][1]=InvInv min brk pos---
                                                                                //beta_vals[period][prevsale][2]=InvInv min brk neg---


        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results\\Regression Results\\Version_1\\Results for all periods_version1_for_run number_5.csv";



        //----------------------------------------------------------------------------------------------------------------
        File file = new File(csvFile);
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line = " ";
        String[] tempArr;
        line = br.readLine();

        while ((line = br.readLine()) != null) {
            tempArr = line.split(",");
            beta_vals[Integer.parseInt(tempArr[0])][Integer.parseInt(tempArr[1])][0]= Double.parseDouble(tempArr[3]); //Intercept
            beta_vals[Integer.parseInt(tempArr[0])][Integer.parseInt(tempArr[1])][1]=Double.parseDouble(tempArr[5]); //Inv min brk pos
            beta_vals[Integer.parseInt(tempArr[0])][Integer.parseInt(tempArr[1])][2]=Double.parseDouble(tempArr[4]) ;//Inv min brk neg

            //for (String tempStr: tempArr) {
              //  System.out.print(tempStr + " ");
            //}
            //System.out.println();
        }
        br.close();

        //--------------Inventory level after ordering---------------------------------------------------------------------------------
        IntStream streamS = IntStream.range(-((period-1) * m), ((period)*m)+1);;
        int MM = (period-1) * m; //neutralizes negative S
        int[] S = streamS.toArray();
        int s = S.length;
        //------------------------------------------------------------------------------------------------------------------------------

        double[] poissonDistProbability = new double[m+1];
        poissonDistProbability = poissonDist (lambda, m);
        double[][] binomialDistProbability = new double[sl][sl+1];
        binomialDistProbability = binomialDist(sl, alpha);

        double[] arr = new double[nresult];
        PrintWriter Solution= new PrintWriter("Integrated_reg_progression_ver1_for_case_"+case_number+".csv");
        arr = SimulationMatrix(period, m, alphaSimulation, lambdaSimulation, unitCost, fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, salvageValue,nresult ,S, sl, s, alpha, lambda, BigM, MM,  poissonDistProbability, binomialDistProbability,Solution,beta_vals);

        System.out.println(Arrays.toString(arr));

        Instant inst2 = Instant.now();
        System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
    }

    public static double[] SimulationMatrix(int period, int m, double alphaSimulation, double lambdaSimulation, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost, int salvageValue,int nresult ,int [] S,int sl, int s, double alpha, double lambda, int BigM, int MM,  double[] poissonDistProbability, double[][] binomialDistProbability,PrintWriter Solution,double[][][] beta_vals) {
        //Solution.println("Order_up_to,Cost");
        Solution.println("Order_up_to,Demand,Inventory, prev_sale,break_point_pos,i_min_k_neg, i_min_k_pos, reg_result,zline,total");



        //int nn = 100000;
        int nn=1;

        double[] arr = new double[nresult];
        double[] arr_sum = new double[nresult];
        double[] arr_mean = new double[nresult];


        for (int run = 0; run < nn; run++) {

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
            double[][] Z_bigSDoublee = new double[period][2];
            double[][] Orderupto_cost_pair_for_period_0 = new double[s][2];

            for (int p = 0; p < period; p++) {
                D[p] = getPoissonRandom(lambdaSimulation);
                if (D[p]>m) {
                    D[p]=m;
                }
            }

            I[0] = 0;
            R[0] = 0;


            //---------------------------------------Z_bigSDouble[0]
            double total;
            int SS = 0;
            double z_line;

            double min = 100000;
            for (int i = 0; i < s; i++) { //S up to order
                total = 0;
                for (int jj = 0; jj <= m; jj++) { //demand
                    double prev_sale=Math.min(i,jj);
                    int prev_sale_int=Math.min(i,jj);
                    double break_point=lambda*(period-1)-alpha*Math.min(i,jj);
                    double break_point_pos=Math.max(break_point,0);
                    double inventory=i-jj;
                    double inv_min_brk_point=inventory-break_point_pos;
                    double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                    double inv_min_brk_point_neg=Math.min(0,inv_min_brk_point);
                    //double reg_result=B0[1] +B1[1]*prev_sale +B2[1]*inv_min_brk_point_pos +B3[1]*inv_min_brk_point_neg*prev_sale;
                    double reg_result_for_period1=beta_vals[1][prev_sale_int][0]+beta_vals[1][prev_sale_int][1]*inv_min_brk_point_pos+beta_vals[1][prev_sale_int][2]*inv_min_brk_point_neg;
                    //System.out.println("Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);

                    if (i>0) {
                        z_line = (poissonDistProbability[jj]) * (i * unitCost
                                + fixedCost +  (Math.max(0,(i - jj)))* holdingCost
                                +(Math.max(0,(jj - i)))* backorderCost
                                +reg_result_for_period1);       //Upcoming period's cost

                    }else {
                        z_line = (poissonDistProbability[jj]) * (
                                (Math.max(0,(i - jj)))* holdingCost
                                        +(Math.max(0,(jj - i)))* backorderCost
                                        + reg_result_for_period1);

                    }
                    //System.out.println(z_line);
                    total += z_line;

                    Solution.println(i+","+jj+","+inventory+","+prev_sale+","+break_point_pos+","+inv_min_brk_point_neg+","+inv_min_brk_point_pos+","+reg_result_for_period1+","+z_line+","+total);
                }
                //Solution.println(i+","+total);

                Orderupto_cost_pair_for_period_0[i][0]=i;
                Orderupto_cost_pair_for_period_0[i][1]=total;
                //Solution.println(Arrays.toString(Orderupto_cost_pair_for_period_0[i]).replace("[","").replace("]",""));


                // System.out.println("Total: " +  i+"  "+"  "+total);
                if (total< min) {
                    min = total;
                    SS = i;

                }

            }

            Z_bigSDouble[0] = SS;
            Z_bigSDoublee[0][0]=SS;
            Z_bigSDoublee[0][1]=min;//Cost değeri
            //Solution.println(Arrays.toString(Z_bigSDoublee[0]).replace("[","").replace("]",""));
            //Burada SS DEĞERİNİ VE min değerini yazdırmak gerekli

            //----------------------------------------------------------------------------------------------------------------------------------------------

            O[0] = Z_bigSDouble[0];

            BO[0] = Math.max(0, D[0] - I[0] - O[0] - R[0]);
            IO[0] = Math.max(0, I[0]);
            SL[0] = Math.min(IO[0] + O[0], D[0]);
            HI[0] = Math.max(0, I[0] + O[0] + R[0] - D[0]);

            //---------------------------------------------------Z_bigSDouble[i]-------------------------------------------------------------------------------------------
            for (int i = 1; i < period-1; i++) {

                I[i] = I[i - 1] + O[i - 1] + R[i - 1] - D[i - 1];
                R[i] = getBinomialRandom(SL[i - 1], alphaSimulation);

                double total_;
                double total2_;
                int SS_ = 0;
                double z_line_;


                double min_ = 100000;

                for (int i_ = 0; i_ < s; i_++) { //S up to order
                    total2_ = 0;
                    for (int jj_ = 0; jj_ <= m; jj_++) { //demand
                        total_ = 0;
                        for (int k_ = 0; k_ <=SL[i-1] ; k_++) { //returniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii <=j vardı sorulack
                            double prev_sale=Math.min(Math.max(0, I[i]) + Math.max(0, (S[i_] - I[i])) + k_, jj_ - Math.min(0, I[i]));
                            int prev_sale_int=Math.min(Math.max(0, I[i]) + Math.max(0, (S[i_] - I[i])) + k_, jj_ - Math.min(0, I[i]));
                            double inventory=S[i_] + k_ - jj_;
                            double break_point=lambda*(period-i-1)-alpha*prev_sale;
                            double break_point_pos=Math.max(break_point,0);
                            double inv_min_brk_point=inventory-break_point_pos;
                            double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                            double inv_min_brk_point_neg=Math.min(0,inv_min_brk_point);
                            //double reg_result=B0[i+1] +B1[i+1]*prev_sale +B2[i+1]*inv_min_brk_point_pos +B3[i+1]*inv_min_brk_point_neg*prev_sale;
                            double reg_result_=beta_vals[i+1][prev_sale_int][0]+beta_vals[i+1][prev_sale_int][1]*inv_min_brk_point_pos+beta_vals[i+1][prev_sale_int][2]*inv_min_brk_point_neg;


                            //System.out.println("pp+1: "+(pp+1)+" Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);


                            if (I[i] < S[i_]) {
                                z_line_ = (binomialDistProbability[SL[i-1]][k_]) *  (poissonDistProbability[jj_]) * ((S[i_]- I[i]) * unitCost
                                        + fixedCost +  Math.max(0,(S[i_] + k_ - jj_))* holdingCost
                                        + Math.max(0,(jj_ - S[i_]  - k_))* backorderCost + k_ * returnCredit
                                        +reg_result_);

                            } else if (I[i] ==S[i_]) {
                                z_line_ = (binomialDistProbability[SL[i-1]][k_]) * (poissonDistProbability[jj_]) *(Math.max(0,(S[i_] + k_ - jj_))* holdingCost
                                        + Math.max(0,(jj_ - S[i_]  - k_))* backorderCost + k_ * returnCredit
                                        +reg_result_);

                            } else {
                                z_line_ = BigM;
                            }
                            total_ += z_line_;

                        }
                        total2_ +=total_;
                        //   System.out.println("Total: " +  i+"  "+jj+"  "+total);
                    }

                    //     System.out.println("Total: " +  i+"  "+"  "+total2);
                    if (total2_< min_) {
                        min_ = total2_;
                        SS_ = i_;
                    }
                    //     System.out.println("Total: " +  i+"  "+"jj"+"  "+min);
                    //Z_bigSDouble[i] = SS_ - MM;
                }
                Z_bigSDouble[i] = SS_ - MM;

                //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                O[i] = Math.max(0, (Z_bigSDouble[i]- I[i]));
                BO[i] = Math.max(0, D[i] - I[i] - O[i] - R[i]);
                IO[i] = Math.max(0, I[i]);
                SL[i] = Math.min(IO[i] + O[i] + R[i], D[i] + BO[i - 1]);
                HI[i] = Math.max(0, I[i] + O[i] + R[i] - D[i]);

            }
            I[period-1] = I[period-2] + O[period-2] + R[period-2] - D[period-2];
            R[period-1] = getBinomialRandom(SL[period-2], alphaSimulation);


            //--------------------------------------------------Z_bigSDouble[period-1] -------------------------------------------------------------------------------------------

            double total__;
            double total2__;
            int SS__ = 0;
            double z_line__;


            double min__ = 100000;

            for (int i__ = 0; i__ < s; i__++) { //S up to order
                total2__ = 0;
                for (int jj__ = 0; jj__ <= m; jj__++) { //demand
                    total__ = 0;
                    for (int k__ = 0; k__ <= SL[period-2]; k__++) { //return
                        if (I[period-1] < S[i__]) {
                            z_line__ = (binomialDistProbability[SL[period-2]][k__]) *  (poissonDistProbability[jj__]) * ((S[i__]- I[period-1]) * unitCost
                                    + fixedCost +  (Math.max(0,(S[i__] + k__ - jj__)))* holdingCost
                                    +(Math.max(0,(jj__ - S[i__]  - k__)))* backorderCost + k__ * returnCredit
                                    + Math.min(Math.max(0, I[period-1]) + Math.max(0, (S[i__] - I[period-1])) + k__, jj__ - Math.min(0, I[period-1])) *alpha * returnCredit+ ( Math.min(0, I[period-1]+k__)*-unitCost )+ ( Math.min(0, I[period-1]+k__)* alpha * (returnCredit - salvageValue) )
                                    +( Math.max(0, I[period-1]+k__)*-salvageValue));
                        } else if (I[period-1] ==S[i__]) {
                            z_line__ = (binomialDistProbability[SL[period-2]][k__]) * (poissonDistProbability[jj__]) *((Math.max(0,(S[i__] + k__ - jj__)))* holdingCost
                                    +(Math.max(0,(jj__ - S[i__]  - k__)))* backorderCost + k__ * returnCredit
                                    + Math.min(Math.max(0, I[period-1]) + Math.max(0, (S[i__] - I[period-1])) + k__, jj__ - Math.min(0, I[period-1]))*alpha * returnCredit+ ( Math.min(0, I[period-1]+k__)*-unitCost )+ ( Math.min(0, I[period-1]+k__)* alpha * (returnCredit - salvageValue) )
                                    +( Math.max(0, I[period-1]+k__)*-salvageValue));
                        } else {
                            z_line__ = BigM;
                        }
                        total__ += z_line__;
                        //    System.out.println("Total: " +  i+"  "+jj+"  "+total);
                    }
                    total2__ +=total__;
                }
                //    System.out.println("Total: " +  i+"  "+"  "+total2);
                if (total2__< min__) {
                    min__ = total2__;
                    SS__ = i__;
                }

            }

            Z_bigSDouble[period-1] = SS__ - MM;

            //--------------------------------------------------------------------------------------------------------------------------------------------


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


            arr = new double[] {tc, pc, oc, hc, bc,tec, rc,  SLV,tc, O[0],O[1],O[2],O[3],O[4], Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3],Z_bigSDouble[4], BO[0],BO[1],BO[2],BO[3],BO[4], HI[0],HI[1],HI[2],HI[3],HI[4], R[1],R[2],R[3],R[4],R[5], SL[0],SL[1],SL[2],SL[3],SL[4],SL[5], I[5]  };


            for(int jj=0;jj<arr.length;jj++){
                arr_sum[jj]+=arr[jj];
            }
        }
        Solution.close();

        for(int jj=0;jj<arr.length;jj++){
            arr_mean[jj]=arr_sum[jj]/nn;
        }
        return arr_mean;


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
