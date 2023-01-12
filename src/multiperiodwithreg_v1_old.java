import java.io.*;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.stream.IntStream;


public class multiperiodwithreg_v1_old {

    public static void main(String[] args) throws IOException {
        Instant inst1 = Instant.now();

        int BigM=10000; int nresult =46;
        int period =4;

        double[] B0 = new double[period];
        double[] B1S1 = new double[period];
        double[] B1S2 = new double[period];
        double[] B2 = new double[period];
        double[] B3 = new double[period];
        double[] B4 = new double[period];
        double[] B5 = new double[period];



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



        //Case-14  0.1	2	20	100	25	8
        double alpha_all=0.1;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost = 100;
        int m = 25; //max demand
        double lambda_all=8;
        String case_number="14";
            */

        //Basecase 0.1 0.4 2 20	25	10	4
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost = 25;
        int m = 10; //max demand
        double lambda_all=4;
        String case_number="base";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_base.csv";


        //----------------------------------------------------------------------------------------------------------------
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
            //Previous Sale 1
            B1S1[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[3]);
            //Previous Sale 2
            B1S2[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[4]);
            //ınv-t_pos
            B2[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[5]);
            //ınv-t_pos^2
            B3[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[6]);
            //ınv-t_neg
            B4[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[7]);
            //ınv-t_neg^2
            B5[Integer.parseInt(tempArr[0])]=Double.parseDouble(tempArr[8]);

        }
        br.close();



        double alpha1 = alpha1_all;

        double alpha2 = alpha2_all;

        double alpha = alpha2/(1-alpha1);

        double alpha1Simulation = alpha1_all;

        double alpha2Simulation = alpha2_all;

        double alphaSimulation  = alpha2Simulation/(1-alpha1Simulation);

        double lambda = lambda_all; //demand rate

        double lambdaSimulation = lambda_all;

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
        PrintWriter Solution;
        if(nn==1){Solution= new PrintWriter("Integrated_reg_progression_multiperiod_for_case_"+case_number+".csv");}
        else{ Solution= new PrintWriter("xxxx");}

        arr = SimulationMatrix(period, m, alpha1Simulation,alpha2Simulation,alphaSimulation, lambdaSimulation,
                unitCost, fixedCost, retailPrice, holdingCost, backorderCost, salvageValue,nresult ,S, sl, s,
                alpha1, alpha2,alpha,
                lambda, BigM, MM,  poissonDistProbability, binomialDistProbability,B0,B1S1,B1S2,B2,B3,B4,B5,Solution,nn);
        PrintWriter Simulation_Solution;
        if(nn==1){Simulation_Solution= new PrintWriter("xxxxxx");}
        else{ Simulation_Solution= new PrintWriter("Integrated_reg_solution_multiperiod_for_case_"+case_number+".csv");}
        //PrintWriter Simulation_Solution= new PrintWriter("Integrated_reg_solution_ver4_for_case_"+case_number+".csv");
        Simulation_Solution.println("tc, pc, oc, hc, bc,tec, rc,  SLV,tc, O[0],O[1],O[2],O[3],O[4], Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3],Z_bigSDouble[4], BO[0],BO[1],BO[2],BO[3],BO[4], HI[0],HI[1],HI[2],HI[3],HI[4], R[1],R[2],R[3],R[4],R[5], SL[0],SL[1],SL[2],SL[3],SL[4],SL[5], I[0], I[1], I[2], I[3], I[4], I[5]");
        Simulation_Solution.println(Arrays.toString(arr));
        Simulation_Solution.close();
        System.out.println(Arrays.toString(arr));
        Instant inst2 = Instant.now();
        System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
    }

    public static double[] SimulationMatrix(int period, int m, double alpha1Simulation,double alpha2Simulation,double alphaSimulation,
                                            double lambdaSimulation, int unitCost, int fixedCost, int [] retailPrice, double holdingCost,
                                            double backorderCost, int salvageValue,int nresult ,int [] S,int sl, int s, double alpha1,double alpha2,double alpha,
                                            double lambda, int BigM, int MM,  double[] poissonDistProbability, double[][] binomialDistProbability,
                                            double[] B0,double[] B1S1,double[] B1S2, double[] B2,double[] B3,double[] B4,double[] B5,PrintWriter Solution,int nn) {
        //Solution.println("Order_up_to,Cost");
        Solution.println("Order_up_to,Demand,Inventory, prev_sale,break_point_pos,i_min_k_neg, i_min_k_pos, reg_result,zline,total");



        //int nn = 100000;

        double[] arr = new double[nresult];
        double[] arr_sum = new double[nresult];
        double[] arr_mean = new double[nresult];


        for (int run = 0; run < nn; run++) {

            int[] I = new int[period + 1];

            int[] R1 = new int[period + 1];

            int[] R2 = new int[period + 2];

            int[] O = new int[period];

            int[] D = new int[period];

            int[] BO = new int[period];

            int[] IO = new int[period + 1];

            int[] SL = new int[period + 1];

            int[] HI = new int[period];

            int[] SNR = new int[period+2];

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
            R1[0] = 0;


            //---------------------------------------Z_bigSDouble[0]
            double total;
            int SS = 0;
            double z_line;

            double min = 100000;
            for (int i = 0; i < s; i++) { //S up to order
                total = 0;
                for (int jj = 0; jj <= m; jj++) { //demand
                    double prev_sale=Math.min(i,jj);
                    double break_point=lambda*(period-1)-alpha1*Math.min(i,jj);
                    double break_point_pos=Math.max(break_point,0);
                    double inventory=i-jj;
                    double inv_min_brk_point=inventory-break_point_pos;
                    double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                    double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                    double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                    double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                    double reg_result=B0[1] +B1S1[1]*prev_sale +B2[1]*inv_min_brk_point_pos +B3[1]*inv_min_brk_point_pos_sqr+B4[1]*inv_min_brk_point_neg +B5[1]*inv_min_brk_point_neg_sqr;
                    //System.out.println("Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);

                    if (i==0) {
                        z_line = (poissonDistProbability[jj])
                                *(jj* backorderCost
                                + reg_result);

                    } else {

                        z_line = (poissonDistProbability[jj])

                                *(i * unitCost + fixedCost
                                + Math.max(0,(i - jj))* holdingCost
                                + Math.max(0,(jj - i ))* backorderCost
                                + reg_result);
                    }

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

            BO[0] = Math.max(0, D[0] - I[0] - O[0] - R1[0]);
            IO[0] = Math.max(0, I[0]);
            SL[0] = Math.min(IO[0] + O[0], D[0]);
            HI[0] = Math.max(0, I[0] + O[0] + R1[0] - D[0]);


            //---------------------------------------------------Z_bigSDouble[1]-------------------------------------------------------------------------------------------


            I[1] = I[0] + O[0] + R1[0] - D[0];
            R1[1] = getBinomialRandom(SL[0], alpha1Simulation);

            double total___;
            double total2___;
            int SS___ = 0;
            double z_line___;


            double min___ = 100000;

            for (int i___ = 0; i___ < s; i___++) { //S up to order
                total2___ = 0;
                for (int jj___ = 0; jj___ <= m; jj___++) { //demand
                    total___ = 0;
                    for (int k___ = 0; k___ <=SL[0] ; k___++) {
                        double prev_sale=Math.min(Math.max(0, I[1]) + Math.max(0, (S[i___] - I[1])) + k___, jj___ - Math.min(0, I[1]));
                        double snr = SL[0] - k___;
                        double inventory=S[i___] + k___ - jj___;
                        double break_point=lambda*(period-1-1)-alpha1*prev_sale-alpha2*snr;
                        double break_point_pos=Math.max(break_point,0);
                        double inv_min_brk_point=inventory-break_point_pos;
                        double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                        double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                        double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                        double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                        double reg_result=B0[2] +B1S1[2]*prev_sale+B1S2[2]*snr +B2[2]*inv_min_brk_point_pos +B3[2]*inv_min_brk_point_pos_sqr+B4[2]*inv_min_brk_point_neg +B5[2]*inv_min_brk_point_neg_sqr;
                        //System.out.println("pp+1: "+(pp+1)+" Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);


                        if (I[1] < S[i___]) { z_line___ = (binomialDistProbability[SL[0]][k___]) *  (poissonDistProbability[jj___]) *
                                ((S[i___]- I[1]) * unitCost  + fixedCost
                                        +  Math.max(0,(S[i___] + k___ - jj___))* holdingCost
                                        + Math.max(0,(jj___ - S[i___]  - k___))* backorderCost
                                        + k___ * retailPrice[0]
                                        +reg_result);

                        } else if (I[1] == S[i___]) { z_line___ = (binomialDistProbability[SL[0]][k___]) *  (poissonDistProbability[jj___]) *
                                (Math.max(0,(S[i___] + k___ - jj___))* holdingCost
                                        + Math.max(0,(jj___ - S[i___]  - k___))* backorderCost
                                        + k___ * retailPrice[0]
                                        + reg_result);

                        } else {
                            z_line___ = BigM;
                        }


                        total___ += z_line___;

                    }
                    total2___ +=total___;
                    //   System.out.println("Total: " +  i+"  "+jj+"  "+total);
                }

                //     System.out.println("Total: " +  i+"  "+"  "+total2);
                if (total2___< min___) {
                    min___ = total2___;
                    SS___ = i___;
                }
                //     System.out.println("Total: " +  i+"  "+"jj"+"  "+min);
                //Z_bigSDouble[i] = SS_ - MM;
            }
            Z_bigSDouble[1] = SS___ - MM;

            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            O[1] = Math.max(0, (Z_bigSDouble[1]- I[1]));
            BO[1] = Math.max(0, D[1] - I[1] - O[1] - R1[1]);
            IO[1] = Math.max(0, I[1]);
            SL[1] = Math.min(IO[1] + O[1] + R1[1], D[1] + BO[0]);
            HI[1] = Math.max(0, I[1] + O[1] + R1[1] - D[1]);
            R2[1] = 0;



            //---------------------------------------------------Z_bigSDouble[i]-------------------------------------------------------------------------------------------
            for (int i = 2; i < period-1; i++) {

                I[i] = I[i-1] + O[i-1] + R1[i-1] + R2[i-1] - D[i-1];

                R1[i] = getBinomialRandom(SL[i-1], alpha1Simulation);

                SNR[i]=SL[i-2]-R1[i-1];

                R2[i] = getBinomialRandom(SNR[i], alphaSimulation);



                double total_;
                double total2_;
                int SS_ = 0;
                double z_line_;


                double min_ = 100000;

                for (int i_ = 0; i_ < s; i_++) { //S up to order
                    total2_ = 0;
                    for (int jj_ = 0; jj_ <= m; jj_++) { //demand
                        total_ = 0;
                        for (int k_ = 0; k_ <=SL[i-1] ; k_++) {

                            double prev_sale=Math.min(Math.max(0, I[i]) + Math.max(0, (S[i_] - I[i])) + k_, jj_ - Math.min(0, I[i]));
                            double snr = SL[i-1] - k_;
                            double inventory=S[i_] + k_ - jj_;
                            double break_point=lambda*(period-i-1)-alpha1*prev_sale-alpha2*snr;
                            double break_point_pos=Math.max(break_point,0);
                            double inv_min_brk_point=inventory-break_point_pos;
                            double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                            double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                            double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                            double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                            double reg_result=B0[i+1] +B1S1[i+1]*prev_sale +B1S2[i+1]*snr
                                    +B2[i+1]*inv_min_brk_point_pos +B3[i+1]*inv_min_brk_point_pos_sqr+B4[i+1]*inv_min_brk_point_neg
                                    +B5[i+1]*inv_min_brk_point_neg_sqr;
                            //System.out.println("pp+1: "+(pp+1)+" Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);


                            if (I[i] < S[i_]) {
                                z_line_ = (binomialDistProbability[SL[i-1]][k_]) * (poissonDistProbability[jj_])
                                        * ((S[i_]- I[i]) * unitCost + fixedCost
                                        +  Math.max(0,(S[i_] + k_ + R2[i]- jj_))* holdingCost
                                        + Math.max(0,(jj_ - S[i_]  - k_ - R2[i]))* backorderCost
                                        + k_ * retailPrice[i-1] + R2[i]* retailPrice[i-2]
                                        +reg_result);

                            } else if (I[i] ==S[i_]) {
                                z_line_ = (binomialDistProbability[SL[i-1]][k_]) * (poissonDistProbability[jj_])
                                        *(Math.max(0,(S[i_] + k_ + R2[i] - jj_))* holdingCost
                                        + Math.max(0,(jj_ - S[i_]  - k_ - R2[i]))* backorderCost
                                        + k_ * retailPrice[i-1] + R2[i]* retailPrice[i-2]
                                        +reg_result);

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


                O[i] = (int)Math.max(0, (Z_bigSDouble[i] - I[i]));

                BO[i] = Math.max(0, D[i] - I[i] - O[i] - R1[i]- R2[i]);

                IO[i] = Math.max(0, I[i]);

                SL[i] = Math.min(IO[i] + O[i] + R1[i] + R2[i], D[i] + BO[i-1]);

                HI[i] = Math.max(0, I[i] + O[i] + R1[i] + R2[i] - D[i]);

            }
            I[period-1] = I[period-2] + O[period-2] + R1[period-2] - D[period-2];

            R1[period-1] = getBinomialRandom(SL[period-2], alpha1Simulation);

            SNR[period-1]=SL[period-3]-R1[period-2];

            R2[period-1] = getBinomialRandom(SNR[period-1], alphaSimulation);


            //--------------------------------------------------Z_bigSDouble[period-1] -------------------------------------------------------------------------------------------

            double total__;
            double total2__;
            int SS__ = 0;
            double z_line__;
            double total3__;


            double min__ = 100000;


            for (int i__ = 0; i__ < s; i__++) { //S up to order
                total3__ = 0;
                for (int jj__ = 0; jj__ <= m; jj__++) { //demand
                    total2__ = 0;
                    for (int k__ = 0; k__ <= SL[period-1-1]; k__++) { //return
                        total__ =0;
                        int prev_sale=Math.min(Math.max(0,I[period-1]) + Math.max(0, (S[i__] - I[period-1])) + k__, jj__ - Math.min(0, I[period-1]));
                        int inventory=S[i__] + k__ - jj__;

                        for (int kk = 0; kk <= prev_sale; kk++) { //return
                            if (I[period-1] < S[i__]) {
                                z_line__ = poissonDistProbability[jj__] * binomialDistProbability[SL[period-1-1]][k__] * binomialDistProbability[prev_sale][kk] *
                                        ((S[i__]- I[period-1]) * unitCost + fixedCost
                                                + Math.max(0,(S[i__] + k__ + kk - jj__))* holdingCost
                                                + Math.max(0,(jj__ - S[i__]  - k__ - kk))* backorderCost
                                                + k__ * retailPrice[period - 2] + kk *retailPrice[period - 1]
                                                + Math.min(0,inventory + kk + (SL[period-1-1]-k__) ) * -unitCost
                                                + Math.max(0,inventory + kk + (SL[period-1-1]-k__) ) * -salvageValue
                                                + kk * retailPrice[period-1] + (SL[period-1-1]-k__) * retailPrice[period-2]
                                                + (SL[period-1-1]-k__) * alpha *  (retailPrice[period-1] - salvageValue)
                                                - Math.min(0,inventory) * (alpha1) * (retailPrice[period-1] - salvageValue)
                                                - Math.min(0,inventory) *(1-alpha1)  * (alpha) * (retailPrice[period-2] - salvageValue));


                            }
                            else if (I[period-1] == S[i__]) {
                                z_line__ = poissonDistProbability[jj__] * binomialDistProbability[SL[period-1-1]][k__] * binomialDistProbability[prev_sale][kk] *
                                        (Math.max(0,(S[i__] + k__ + kk - jj__))* holdingCost
                                                + Math.max(0,(jj__ - S[i__]  - k__ - kk))* backorderCost
                                                + k__ * retailPrice[period - 2] + kk *retailPrice[period - 1]
                                                + Math.min(0,inventory + kk + (SL[period-1-1]-k__) ) * -unitCost
                                                + Math.max(0,inventory + kk + (SL[period-1-1]-k__) ) * -salvageValue
                                                + kk * retailPrice[period-1] + (SL[period-1-1]-k__) * retailPrice[period-2]
                                                + (SL[period-1-1]-k__) * alpha *  (retailPrice[period-1] - salvageValue)
                                                - Math.min(0,inventory) * (alpha1) * (retailPrice[period-1] - salvageValue)
                                                - Math.min(0,inventory) *(1-alpha1)  * (alpha) * (retailPrice[period-2] - salvageValue));


                            }


                            else {
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


            BO[period-1] = Math.max(0, D[period-1] - I[period-1] - O[period-1] - R1[period-1] - R2[period-1]);
            IO[period-1] = Math.max(0, I[period-1]);
            SL[period-1] = Math.min(IO[period-1] + O[period-1] + R1[period-1] + R2[period-1], D[period-1] + BO[period-2]);
            HI[period-1] = Math.max(0, I[period-1] + O[period-1] + R1[period-1] + R2[period-1] - D[period-1]);


            I[period] = I[period - 1] + O[period - 1] + R1[period-1] + R2[period-1] - D[period - 1];

            R1[period] = getBinomialRandom(SL[period - 1], alphaSimulation);

            SNR[period]=SL[period-2]-R1[period-1];

            R2[period] = getBinomialRandom(SNR[period], alphaSimulation);

            double procurementCost = 0;

            double totalOrderingCost = 0;

            double totalHoldingCost = 0;

            double totalBackorderCost = 0;

            double returnCost = 0;

            double returnCost1 = 0;

            double totalEndingInventoryCost = 0;


            for (int p = 0; p < period; p++) {


                procurementCost += O[p] * unitCost;

                totalOrderingCost += Math.min(1,O[p]) * fixedCost;

                totalHoldingCost += HI[p] * holdingCost;

                totalBackorderCost += BO[p] * backorderCost;


            }
            for (int p = 2; p < period; p++) {
                returnCost1 += R1[p]* retailPrice[p-1] + R2[p]* retailPrice[p-2];
            }

            returnCost += returnCost1 + R1[1]* retailPrice[0] ;

            SNR[period]=SL[period-2]-R1[period-1];

            R2[period] = getBinomialRandom(SNR[period], alphaSimulation);

            SL[period]= -(Math.min(0,I[period]));

            SNR[period+1]=SL[period-1]-R1[period];

            R2[period+1] = getBinomialRandom(SNR[period+1], alphaSimulation);

            totalEndingInventoryCost =
                    + Math.min(0,(I[period] + R1[period] + R2[period])) * -unitCost
                            + Math.max(0,(I[period] + R1[period] + R2[period] )) * -salvageValue
                            + R1[period] * retailPrice[period-1] + R2[period] * retailPrice[period-2]
                            + R2[period+1] * (retailPrice[period-1] - salvageValue)
                            - Math.min(0,I[period]) * alpha1Simulation  * (retailPrice[period-1] - salvageValue)
                            - Math.min(0,I[period]) * (1-alpha1Simulation) * (alphaSimulation)  * (retailPrice[period-2] - salvageValue);

            double salvage=I[period]+R1[period]+R2[period];
            double SLV= salvage;

            double totalCost = procurementCost + totalOrderingCost + totalHoldingCost + totalBackorderCost + totalEndingInventoryCost + returnCost;


            double pc = procurementCost;

            double oc = totalOrderingCost;

            double hc = totalHoldingCost;

            double bc = totalBackorderCost;

            double tec = totalEndingInventoryCost;

            double rc = returnCost;

            double tc = totalCost;


            arr = new double[] {tc, pc, oc, hc, bc,tec, rc,tc,
                    I[0], I[1], I[2], I[3], I[4],
                    D[0],D[1],D[2],D[3],
                    Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3],
                    O[0],O[1],O[2],O[3],
                    BO[0],BO[1],BO[2],BO[3],
                    HI[0],HI[1],HI[2],HI[3],
                    SL[0],SL[1],SL[2],SL[3],SL[4],
                    R1[1],R1[2],R1[3],R1[4]
                    ,R2[2],R2[3],R2[4],R2[5]

                      };


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
