import java.io.*;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.stream.IntStream;


public class cost_version4_berk_fixed {

    public static void main(String[] args) throws IOException {
        Instant inst1 = Instant.now();

        int BigM=10000; int nresult =46;



        /*
        //Basecase 0.5	2	20	25	10	4
        double alpha_all=0.5;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost = 25;
        int m = 10; //max demand
        double lambda_all=4;
        String case_number="base";

        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASEBASE.csv";


         //Case-1 0.4	2	20	40	25	8
        double alpha_all=0.4;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost = 40;
        int m = 25; //max demand
        double lambda_all=8;
        String case_number="1";

        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE1.csv";


         //Case-2  0.4	2.4	20	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 20;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="2";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE2.csv";

            //Case-3  	0.4	2	32	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 32;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="3";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE3.csv";


            //Case-4 0.4	2.4	32	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 32;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="4";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE4.csv";


            //Case-5  0.4	2	50	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 50;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="5";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE5.csv";


            //Case-6  0.5	2	50	40	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 50;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="6";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE6.csv";


            //Case-7 0.4	2.4	50	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 50;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="7";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE7.csv";


            //Case-8 0.4	2	20	50	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 50;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="8";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE8.csv";


            //Case-9  0.4	2.4	20	50	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 20;
            int fixedCost = 50;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="9";

            String csvFile = "C:\\Users\\berkc\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE9.csv";


            //Case-10  0.5	2	20	100	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="10";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE10.csv";


            //Case-11  0.5	2	20	45	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 45;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="11";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE11.csv";


            //Case-12  0.4	2	50	80	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 50;
            int fixedCost = 80;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="12";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE12.csv";


            //Case-13  0	2	20	100	25	8
            double alpha_all=0;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="13";

            String csvFile = "C:\\Users\\berkc\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE13.csv";


            //Case-14  0.1	2	20	100	25	8
            double alpha_all=0.1;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="14";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE14.csv";

            //Case-15 0.4	2	5	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 5;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="15";

            String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE15.csv";

             //Case-16 0.5	2	5	100	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 5;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="16";

              String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_version_4_CASE16.csv";

    */

        //Case-18  0.5	2	20	100	25	8 period=2
        double alpha_all=0.5;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost = 100;
        int m = 25; //max demand
        double lambda_all=8;
        String case_number="18";

        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\All Results-2\\Regression Results\\Results for all periods_Multiperiod_version2_CASE_0.csv";

        int period =4;

        double[] B0 = new double[period];
        double[] B1 = new double[period];
        double[] B2 = new double[period];
        double[] B3 = new double[period];
        double[] B4 = new double[period];
        double[] B5 = new double[period];


        int nn=100000;   //100000
        //----------------------------------------------------------------------------------------------------------------
        File file = new File(csvFile);
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line = " ";
        String[] tempArr;
        line = br.readLine();

        while ((line = br.readLine()) != null) {
            tempArr = line.split(",");
            //Intercept
            B0[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[2]);
            B0[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[2]);
            B0[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[2]);
            B0[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[2]);
            //Previous Sale
            B1[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[3]);
            B1[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[3]);
            B1[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[3]);
            B1[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[3]);
            //ınv-t_pos
            B2[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[4]);
            B2[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[4]);
            B2[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[4]);
            B2[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[4]);
            //ınv-t_pos^2
            B3[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[5]);
            B3[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[5]);
            B3[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[5]);
            B3[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[5]);
            //ınv-t_neg
            B4[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[6]);
            B4[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[6]);
            B4[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[6]);
            B4[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[6]);
            //ınv-t_neg^2
            B5[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[7]);
            B5[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[7]);
            B5[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[7]);
            B5[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[7]);

        }
        br.close();


        int unitCost = 40;

        int salvageValue = 0;

        int returnCredit = 80;

        int retailPrice = 80;

        double alpha = alpha_all; //return rate

        double lambda = lambda_all; //demand rate

        double alphaSimulation = alpha_all;

        double lambdaSimulation = lambda_all;

        int period_for_calculation=5;

        //--------------Inventory level after ordering---------------------------------------------------------------------------------
        IntStream streamS = IntStream.range(-((period_for_calculation-1) * m), ((period_for_calculation)*m)+1);;
        int MM = (period_for_calculation-1) * m; //neutralizes negative S
        int[] S = streamS.toArray();
        int s = S.length;
        //------------------------------------------------------------------------------------------------------------------------------

        int sl =((period_for_calculation+2) * m)+1; //max # of sales

        double[] poissonDistProbability = new double[m+1];
        poissonDistProbability = poissonDist (lambda, m);
        double[][] binomialDistProbability = new double[sl][sl+1];
        binomialDistProbability = binomialDist(sl, alpha);

        double[] arr = new double[nresult];
        PrintWriter Solution;
        PrintWriter Solution1;

        if(nn==1){Solution= new PrintWriter("Integrated_reg_progression_ver4_for_case_"+case_number+"period0.csv");}
        else{ Solution= new PrintWriter("xxxx");}
        if(nn==1){Solution1= new PrintWriter("Integrated_reg_progression_ver4_for_case_"+case_number+"period1.csv");}
        else{ Solution1= new PrintWriter("xxxx");}

        arr = SimulationMatrix(period, m, alphaSimulation, lambdaSimulation, unitCost, fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, salvageValue,nresult ,S, sl, s, alpha, lambda, BigM, MM,  poissonDistProbability, binomialDistProbability,B0,B1,B2,B3,B4,B5,Solution,nn,Solution1);
        PrintWriter Simulation_Solution;
        if(nn==1){Simulation_Solution= new PrintWriter("xxxxxx");}
        else{ Simulation_Solution= new PrintWriter("Integrated_reg_solution_ver4_for_case_"+case_number+".csv");}
        //PrintWriter Simulation_Solution= new PrintWriter("Integrated_reg_solution_ver4_for_case_"+case_number+".csv");
        Simulation_Solution.println("tc, pc, oc, hc, bc,tec, rc,  SLV,tc, O[0],O[1],O[2],O[3],O[4], Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3],Z_bigSDouble[4], BO[0],BO[1],BO[2],BO[3],BO[4], HI[0],HI[1],HI[2],HI[3],HI[4], R[1],R[2],R[3],R[4],R[5], SL[0],SL[1],SL[2],SL[3],SL[4],SL[5], I[0], I[1], I[2], I[3], I[4], I[5]");
        Simulation_Solution.println(Arrays.toString(arr));
        Simulation_Solution.close();
        System.out.println(Arrays.toString(arr));
        Instant inst2 = Instant.now();
        System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
    }

    public static double[] SimulationMatrix(int period, int m, double alphaSimulation, double lambdaSimulation, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost, int salvageValue,int nresult ,int [] S,int sl, int s, double alpha, double lambda, int BigM, int MM,  double[] poissonDistProbability, double[][] binomialDistProbability,double[] B0,double[] B1,double[] B2,double[] B3,double[] B4,double[] B5,PrintWriter Solution,int nn,PrintWriter Solution1) {
        //Solution.println("Order_up_to,Cost");
        Solution.println("Order_up_to,Demand,Inventory, prev_sale,break_point_pos,i_min_k_neg, i_min_k_pos, reg_result,zline,total");
        Solution1.println("Order_up_to,Demand,Return,Inventory, prev_sale,break_point_pos,i_min_k_neg, i_min_k_pos, reg_result,zline_,total");



        //int nn = 100000;

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
                    double break_point=lambda*(period-1)-alpha*Math.min(i,jj);
                    double break_point_pos=Math.max(break_point,0);
                    double inventory=i-jj;
                    double inv_min_brk_point=inventory-break_point_pos;
                    double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                    double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                    double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                    double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                    double reg_result=B0[1] +B1[1]*prev_sale +B2[1]*inv_min_brk_point_pos +B3[1]*inv_min_brk_point_pos_sqr+B4[1]*inv_min_brk_point_neg +B5[1]*inv_min_brk_point_neg_sqr;
                    //System.out.println("Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);

                    z_line = (poissonDistProbability[jj]) * (i * unitCost
                            + Math.min(1, i)*fixedCost +  (Math.max(0,(i - jj)))* holdingCost
                            +(Math.max(0,(jj - i)))* backorderCost
                            +reg_result);       //Upcoming period's cost


                    //System.out.println(z_line);
                    total += z_line;
                    if(nn==1) {
                        Solution.println(i+","+jj+","+inventory+","+prev_sale+","+break_point_pos+","+inv_min_brk_point_neg+","+inv_min_brk_point_pos+","+reg_result+","+z_line+","+total);
                    }}

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
                            double inventory=S[i_] + k_ - jj_;
                            double break_point=lambda*(period-i-1)-alpha*prev_sale;
                            double break_point_pos=Math.max(break_point,0);
                            double inv_min_brk_point=inventory-break_point_pos;
                            double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                            double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                            double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                            double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                            double reg_result=B0[i+1] +B1[i+1]*prev_sale +B2[i+1]*inv_min_brk_point_pos +B3[i+1]*inv_min_brk_point_pos_sqr+B4[i+1]*inv_min_brk_point_neg +B5[i+1]*inv_min_brk_point_neg_sqr;
                            //System.out.println("pp+1: "+(pp+1)+" Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);

                            /*
                            if (I[i] < S[i_]) {
                                z_line_ = (binomialDistProbability[SL[i-1]][k_]) *  (poissonDistProbability[jj_]) * ((S[i_]- I[i]) * unitCost
                                        + fixedCost +  Math.max(0,(S[i_] + k_ - jj_))* holdingCost
                                        + Math.max(0,(jj_ - S[i_]  - k_))* backorderCost + k_ * returnCredit
                                        +reg_result);

                            } else if (I[i] ==S[i_]) {
                                z_line_ = (binomialDistProbability[SL[i-1]][k_]) * (poissonDistProbability[jj_]) *(Math.max(0,(S[i_] + k_ - jj_))* holdingCost
                                        + Math.max(0,(jj_ - S[i_]  - k_))* backorderCost + k_ * returnCredit
                                        +reg_result);

                            }
                             */

                            if (I[i] <= S[i_]) { z_line_ = (binomialDistProbability[SL[i-1]][k_]) *  (poissonDistProbability[jj_]) * ((S[i_]- I[i]) * unitCost
                                    + Math.min(1, S[i_]-I[i])*fixedCost +  Math.max(0,(S[i_] + k_ - jj_))* holdingCost
                                    + Math.max(0,(jj_ - S[i_]  - k_))* backorderCost + k_ * returnCredit +reg_result);}
                            else {
                                z_line_ = BigM;
                            }


                            total_ += z_line_;
                            if(nn==1 & i==1) {
                                Solution1.println(i_+","+jj_+","+k_+","+inventory+","+prev_sale+","+break_point_pos+","+inv_min_brk_point_neg+","+inv_min_brk_point_pos+","+reg_result+","+z_line_+","+total_);
                            }


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
            double total3__;


            double min__ = 100000;


            /*
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
             */
            for (int i__ = 0; i__ < s; i__++) { //S up to order
                total3__ = 0;
                for (int jj__ = 0; jj__ <= m; jj__++) { //demand
                    total2__ = 0;
                    for (int k__ = 0; k__ <= SL[period-1-1]; k__++) { //return
                        total__ =0;
                        int prev_sale=Math.min(Math.max(0,I[period-1]) + Math.max(0, (S[i__] - I[period-1])) + k__, jj__ - Math.min(0, I[period-1]));
                        int inventory=S[i__] + k__ - jj__;

                        for (int kk = 0; kk <= prev_sale; kk++) { //return
                            if (I[period-1] <= S[i__]) {
                                z_line__ = poissonDistProbability[jj__] * binomialDistProbability[SL[period-1-1]][k__] * binomialDistProbability[prev_sale][kk] *
                                        ((S[i__]- I[period-1]) * unitCost
                                                + Math.min(1, S[i__]-I[period-1])*fixedCost +  Math.max(0,(S[i__] + k__ - jj__))* holdingCost
                                                + Math.max(0,(jj__ - S[i__]  - k__))* backorderCost + k__ * returnCredit
                                                + Math.min(0,inventory + kk) * -unitCost
                                                + Math.max(0,inventory + kk) * -salvageValue
                                                + kk * returnCredit
                                                - Math.min(0,inventory) * alpha * (returnCredit - salvageValue));
                            } else {
                                z_line__ = BigM;
                            }
                            total__ += z_line__;
                            //    System.out.println("Total: " +  i+"  "+jj+"  "+total);
                        }
                        total2__ +=total__;
                    }
                    total3__ +=total2__;
                }
                //          System.out.println("Total: " +  i+" " +pp+"  "+"  "+total2);
                if (total3__< min__) {
                    min__ = total3__;
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
                            - I[period] * alphaSimulation * (returnCredit - salvageValue)
                    +R[period] * returnCredit;
                } else if (I[period] < 0) {
                    totalEndingInventoryCost = (I[period] + R[period]) * -salvageValue
                            - I[period] * alphaSimulation * (returnCredit - salvageValue)
                    + R[period] * returnCredit;
                } else {
                    totalEndingInventoryCost = (I[period] + R[period]) * -salvageValue+
                            R[period] * returnCredit;
                }

                if (I[period] < 0) {
                    SL[period] = -(I[period]);
                } else {
                    SL[period] = 0;
                }

                returnCost += R[p] * returnCredit;


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

            // result = new double[] {tc, pc, oc, hc, bc,tec, rc, SLV,tc, O[0],O[1], bigSValue[0], bigSValue[1], BO[0],BO[1], HI[0],HI[1],R[1],R[2], SL[0],SL[1],SL[2], I[0],I[1]};
            arr = new double[] {tc, pc, oc, hc, bc,tec, rc,SLV,tc,O[0],O[1],Z_bigSDouble[0],Z_bigSDouble[1], BO[0],BO[1], HI[0],HI[1], R[1],R[2], SL[0],SL[1], I[0], I[1]};
            //arr = new double[] {tc, pc, oc, hc, bc,tec, rc,SLV,tc,O[0],O[1],O[2],O[3],O[4],Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3],Z_bigSDouble[4], BO[0],BO[1],BO[2],BO[3],BO[4], HI[0],HI[1],HI[2],HI[3],HI[4], R[1],R[2],R[3],R[4],R[5], SL[0],SL[1],SL[2],SL[3],SL[4],SL[5], I[0], I[1], I[2], I[3], I[4], I[5]  };


            for(int jj=0;jj<arr.length;jj++){
                arr_sum[jj]+=arr[jj];
            }
        }
        Solution.close();
        Solution1.close();

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
