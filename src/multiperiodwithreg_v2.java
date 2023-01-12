import java.io.*;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.stream.IntStream;


public class multiperiodwithreg_v2 {

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
        //Basecase
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";



         //Case-1
         double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost =2 ;
        double backorderCost = 20;
        int fixedCost = 40;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="1";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-2
         double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost =2.4 ;
        double backorderCost =20 ;
        int fixedCost =40 ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="2";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-3
         double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = 2;
        double backorderCost = 50;
        int fixedCost = 40;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="3";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-4
         double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost =100 ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="4";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-5
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost =2 ;
        double backorderCost =5 ;
        int fixedCost = 100;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="5";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-6
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="6";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-7
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="7";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-8
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="8";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-9
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="9";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


         //Case-10
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="10";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";


        //Empty
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = ;
        double backorderCost = ;
        int fixedCost = ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";

            */

        //Basecase
        double alpha1_all=0.1;
        double alpha2_all=0.4;
        double holdingCost = 2;
        double backorderCost = 20;
        int fixedCost =25 ;
        int m = 10; //max demand
        double lambda_all=4;

        String case_number="base";
        int unitCost = 20;
        int salvageValue = 0;
        int [] retailPrice =  {80,80,80,80,80,80,80};
        String csvFile = "C:\\Users\\berkc\\IdeaProjects\\personal\\Multiperiod Results-2\\Regression Results\\Regression Coefficients for all periods_Multiperiod_CASE_"+case_number+".csv";



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

        //--------------Inventory-----------------------------------------------------------------------------------------------------
        IntStream streamII = IntStream.range(-(period * m),((period)*4*m ));
        int M =period * m; //neutralizes negative I
        int[] II = streamII.toArray();
        int n = II.length;
        //-----------------------------------------------------------------------------------------------------------------------------

        //--------------Inventory level after ordering---------------------------------------------------------------------------------
        IntStream streamS = IntStream.range(-((period-1) * m), ((period)*m));
        int MM = (period-1) * m; //neutralizes negative S
        int[] S = streamS.toArray();
        int s = S.length;
        //------------------------------------------------------------------------------------------------------------------------------

        int sl =((period+1) * m)+1;

        double[] poissonDistProbability = new double[m+1];
        poissonDistProbability = poissonDist (lambda, m);

        double[][] binomialDistProbability1 = new double[sl][sl+1];
        binomialDistProbability1 = binomialDist(sl, alpha1);

        double[][] binomialDistProbability2 = new double[sl][sl+1];
        binomialDistProbability2 = binomialDist(sl, alpha);

        double[][][][][] Z_bigSDoubleee = new double[period+1][2][n][sl][sl];

        Z_bigSDoubleee[period][0] = ZEndcalculate(II, period,n, sl, m,  alpha1, alpha2, alpha, backorderCost,
                salvageValue, retailPrice,fixedCost,unitCost,
                binomialDistProbability1,  binomialDistProbability2);



        Z_bigSDoubleee[period-1] = periodMatrixmodifiedv4_2(II , S, period,  Z_bigSDoubleee[period][0], n, sl, s, m,
                lambda, unitCost, fixedCost,
                retailPrice, holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,
                binomialDistProbability1, binomialDistProbability2, period-1);

        double[] arr = new double[nresult];
        PrintWriter Solution;
        if(nn==1){Solution= new PrintWriter("Integrated_reg_progression_multiperiod_for_case_"+case_number+".csv");}
        else{ Solution= new PrintWriter("xxxx");}

        arr = SimulationMatrix(period, m, alpha1Simulation,alpha2Simulation,alphaSimulation, lambdaSimulation,
                unitCost, fixedCost, retailPrice, holdingCost, backorderCost, salvageValue,nresult ,S, sl, s,
                alpha1, alpha2,alpha,
                lambda, BigM, MM,  poissonDistProbability, binomialDistProbability1, binomialDistProbability2,
                B0,B1S1,B1S2,B2,B3,B4,B5,Solution,nn,II,M,n, Z_bigSDoubleee[period-1]);
        PrintWriter Simulation_Solution;
        if(nn==1){Simulation_Solution= new PrintWriter("xxxxxx");}
        else{ Simulation_Solution= new PrintWriter("Integrated_reg_solution_multiperiod_for_case_"+case_number+".csv");}
        //PrintWriter Simulation_Solution= new PrintWriter("Integrated_reg_solution_ver4_for_case_"+case_number+".csv");
        Simulation_Solution.println("tc, pc, oc, hc, bc,tec, rc,tc, I[0], I[1], I[2], I[3], I[4], D[0],D[1],D[2],D[3], Z_bigSDouble[0],Z_bigSDouble[1],Z_bigSDouble[2],Z_bigSDouble[3], O[0],O[1],O[2],O[3], BO[0],BO[1],BO[2],BO[3], HI[0],HI[1],HI[2],HI[3], SL[0],SL[1],SL[2],SL[3],SL[4], R1[1],R1[2],R1[3],R1[4],R2[2],R2[3],R2[4],R2[5]");
        Simulation_Solution.println(Arrays.toString(arr));
        Simulation_Solution.close();
        System.out.println(Arrays.toString(arr));
        Instant inst2 = Instant.now();
        System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
    }



    public static double[] SimulationMatrix(int period, int m, double alpha1Simulation,double alpha2Simulation,double alphaSimulation,
                                            double lambdaSimulation, int unitCost, int fixedCost, int [] retailPrice, double holdingCost,
                                            double backorderCost, int salvageValue,int nresult ,int [] S,int sl, int s, double alpha1,
                                            double alpha2,double alpha, double lambda, int BigM, int MM,  double[] poissonDistProbability,
                                            double[][] binomialDistProbability1, double[][] binomialDistProbability2,
                                            double[] B0,double[] B1S1,double[] B1S2, double[] B2,double[] B3,double[] B4,double[] B5,
                                            PrintWriter Solution,int nn, int []II, int M, int n, double [][][][] Z_bigSDoubleee) {
        //Solution.println("Order_up_to,Cost");
        Solution.println("Order_up_to,Demand,Inventory, prev_sale,break_point_pos,i_min_k_neg, i_min_k_pos, reg_result,zline,total");




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
                    double break_point=lambda*(period-1)-alpha1*prev_sale;
                    double break_point_pos=Math.max(break_point,0);
                    double inventory=i-jj;
                    double inv_min_brk_point=inventory-break_point_pos;
                    double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                    double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                    double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                    double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                    double reg_result=B0[1] +B1S1[1]*prev_sale +B1S2[1] +B2[1]*inv_min_brk_point_pos +B3[1]*inv_min_brk_point_pos_sqr+B4[1]*inv_min_brk_point_neg +B5[1]*inv_min_brk_point_neg_sqr;
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
                    //   if(nn==1) {
                    //       Solution.println(i+","+jj+","+inventory+","+prev_sale+","+break_point_pos+","+inv_min_brk_point_neg+","+inv_min_brk_point_pos+","+reg_result+","+z_line+","+total);
                    //  }

                }

                //        Orderupto_cost_pair_for_period_0[i][0]=i;
                //      Orderupto_cost_pair_for_period_0[i][1]=total;
                //Solution.println(Arrays.toString(Orderupto_cost_pair_for_period_0[i]).replace("[","").replace("]",""));


                // System.out.println("Total: " +  i+"  "+"  "+total);
                if (total< min) {
                    min = total;
                    SS = i;

                }

            }
            Z_bigSDouble[0] = SS;
            //   Z_bigSDouble[0] = 11;
            Z_bigSDoublee[0][0]=SS;
            Z_bigSDoublee[0][1]=min;//Cost değeri
            //Solution.println(Arrays.toString(Z_bigSDoublee[0]).replace("[","").replace("]",""));
            //Burada SS DEĞERİNİ VE min değerini yazdırmak gerekli

            //----------------------------------------------------------------------------------------------------------------------------------------------

            O[0] = (int)Z_bigSDouble[0];

            BO[0] = Math.max(0, D[0] - I[0] - O[0] );
            IO[0] = Math.max(0, I[0]);
            SL[0] = Math.min(IO[0] + O[0], D[0]);
            HI[0] = Math.max(0, I[0] + O[0] - D[0]);


            //---------------------------------------------------Z_bigSDouble[1]-------------------------------------------------------------------------------------------


            I[1] = I[0] + O[0] - D[0];
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
                        double break_point=lambda*(period-2)-alpha1*prev_sale-alpha*snr;
                        double break_point_pos=Math.max(break_point,0);
                        double inv_min_brk_point=inventory-break_point_pos;
                        double inv_min_brk_point_pos=Math.max(0,inv_min_brk_point);
                        double inv_min_brk_point_neg=Math.max(0,inv_min_brk_point*-1);
                        double inv_min_brk_point_pos_sqr=inv_min_brk_point_pos*inv_min_brk_point_pos;
                        double inv_min_brk_point_neg_sqr=inv_min_brk_point_neg*inv_min_brk_point_neg;
                        double reg_result=B0[2] +B1S1[2]*prev_sale+B1S2[2]*snr +B2[2]*inv_min_brk_point_pos +B3[2]*inv_min_brk_point_pos_sqr+B4[2]*inv_min_brk_point_neg +B5[2]*inv_min_brk_point_neg_sqr;
                        //System.out.println("pp+1: "+(pp+1)+" Inventory: "+inventory+" prev_sale: "+prev_sale+" break_point_pos: "+break_point_pos+" i_min_k_neg: "+inv_min_brk_point_neg+" i_min_k_pos: "+inv_min_brk_point_pos+" reg_result: "+reg_result);


                        if (I[1] < S[i___]) { z_line___ = (binomialDistProbability1[SL[0]][k___]) *  (poissonDistProbability[jj___]) *
                                ((S[i___]- I[1]) * unitCost  + fixedCost
                                        +  Math.max(0,(S[i___] + k___ - jj___))* holdingCost
                                        + Math.max(0,(jj___ - S[i___]  - k___))* backorderCost
                                        + k___ * retailPrice[0]
                                        +reg_result);

                        } else if (I[1] == S[i___]) { z_line___ = (binomialDistProbability1[SL[0]][k___]) *  (poissonDistProbability[jj___]) *
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
            O[1] = (int)Math.max(0, (Z_bigSDouble[1]- I[1]));
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
                double total1_;
                double total2_;
                int SS_ = 0;
                double z_line_;


                double min_ = 100000;

                for (int i_ = 0; i_ < s; i_++) { //S up to order
                    total2_ = 0;
                    for (int jj_ = 0; jj_ <= m; jj_++) { //demand
                        total1_ = 0;
                        for (int k_ = 0; k_ <=SL[i-1] ; k_++) {
                            total_ = 0;

                            for (int k2_ = 0; k2_ <= SNR[i] ; k2_++) {

                                double prev_sale=Math.min(Math.max(0, I[i]) + Math.max(0, (S[i_] - I[i])) + k_+ k2_, jj_ - Math.min(0, I[i]));
                                double snr = SL[i-1] - k_;
                                double inventory=S[i_] + k_ + k2_- jj_;
                                double break_point=lambda*(period-i-1)-alpha1*prev_sale-alpha*snr;
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
                                    z_line_ = (binomialDistProbability1[SL[i-1]][k_]) * (binomialDistProbability2[SNR[i]][k2_])
                                            *(poissonDistProbability[jj_])
                                            * ((S[i_]- I[i]) * unitCost + fixedCost
                                            +  Math.max(0,(S[i_] + k_ + k2_- jj_))* holdingCost
                                            + Math.max(0,(jj_ - S[i_]  - k_ - k2_))* backorderCost
                                            + k_ * retailPrice[i-1] + k2_* retailPrice[i-2]
                                            +reg_result);

                                } else if (I[i] ==S[i_]) {
                                    z_line_ = (binomialDistProbability1[SL[i-1]][k_]) * (binomialDistProbability2[SNR[i]][k2_])
                                            * (poissonDistProbability[jj_])
                                            *(Math.max(0,(S[i_] + k_ + k2_ - jj_))* holdingCost
                                            + Math.max(0,(jj_ - S[i_]  - k_ - k2_))* backorderCost
                                            + k_ * retailPrice[i-1] + k2_* retailPrice[i-2]
                                            +reg_result);

                                } else {
                                    z_line_ = BigM;
                                }


                                total_ += z_line_;
                            }
                            total1_ +=total_;
                        }
                        total2_ +=total1_;
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
            I[period-1] = I[period-2] + O[period-2] + R1[period-2] + R2[period-2] - D[period-2];

            R1[period-1] = getBinomialRandom(SL[period-2], alpha1Simulation);

            SNR[period-1]=SL[period-3]-R1[period-2];

            R2[period-1] = getBinomialRandom(SNR[period-1], alphaSimulation);


            //--------------------------------------------------Z_bigSDouble[period-1] -------------------------------------------------------------------------------------------
/*
            double total__;
            double total1__;
            double total2__;
            int SS__ = 0;
            double z_line__;
            double total3__;
            double total4__;


            double min__ = 100000;


            for (int i__ = 0; i__ < s; i__++) { //S up to order

            	total4__ = 0;
                for (int jj__ = 0; jj__ <= m; jj__++) { //demand

                	total3__ = 0;
                    for (int k__ = 0; k__ <= SL[period-2]; k__++) { //return //R1[period-1]


                    	total2__ =0;


                        for (int k2__ = 0; k2__ <= SNR[period-1]; k2__++) { //return //R2[period-1]

                        total1__ =0;
                        int prev_sale__=Math.min(Math.max(0,I[period-1]) + Math.max(0, (S[i__] - I[period-1])) + k__ + k2__, jj__ - Math.min(0, I[period-1]));
                        int inventory__=S[i__] + k__ + k2__ - jj__;



                        for (int kk = 0; kk <= prev_sale__; kk++) { //return //R1[period]

                        	total__ =0;
                        	 int snr__ = SL[period-2] - k__;
                        	 int snr_kk = prev_sale__ - kk;

                        	for (int kk2 = 0; kk2 <= snr__; kk2++) { //return  //R2[period]

                        		if (I[period-1] < S[i__]) {
                                    z_line__ = (poissonDistProbability[jj__])* (binomialDistProbability1[SL[period-2]][k__])
                                    		 * (binomialDistProbability2[SNR[period-1]][k2__])* (binomialDistProbability1[prev_sale__][kk])
                                    		 * (binomialDistProbability2[snr__][kk2])*
                                    		   ((S[i__]- I[period-1]) * unitCost + fixedCost
                                    		   + Math.max(0,(S[i__] + k__ + k2__ - jj__))* holdingCost
                                               + Math.max(0,(jj__ - S[i__]  - k__ - k2__))* backorderCost
                                               + k__ * retailPrice[period - 2] + k2__ *retailPrice[period - 3]


                                               + Math.min(0,inventory__ + kk + kk2 ) * -unitCost
                                               + Math.max(0,inventory__ + kk + kk2 ) * -salvageValue
                                               + kk * retailPrice[period-1] + kk2 * retailPrice[period-2]
                                               + snr_kk * alpha *  (retailPrice[period-1] - salvageValue)
                                               - Math.min(0,inventory__) * (alpha1) * (retailPrice[period-1] - salvageValue)
                                               - Math.min(0,inventory__) *(1-alpha1)  * (alpha) * (retailPrice[period-2] - salvageValue));


                        		}
                        		if (I[period-1] == S[i__]) {
                        			z_line__ = (poissonDistProbability[jj__])* (binomialDistProbability1[SL[period-2]][k__])
                                   		 * (binomialDistProbability2[SNR[period-1]][k2__])* (binomialDistProbability1[prev_sale__][kk])
                                   		 * (binomialDistProbability2[snr__][kk2])*
                                   		   (Math.max(0,(S[i__] + k__ + k2__ - jj__))* holdingCost
                                              + Math.max(0,(jj__ - S[i__]  - k__ - k2__))* backorderCost
                                              + k__ * retailPrice[period - 2] + k2__ *retailPrice[period - 3]


                                              + Math.min(0,inventory__ + kk + kk2 ) * -unitCost
                                              + Math.max(0,inventory__ + kk + kk2 ) * -salvageValue
                                              + kk * retailPrice[period-1] + kk2 * retailPrice[period-2]
                                              + snr_kk * alpha *  (retailPrice[period-1] - salvageValue)
                                              - Math.min(0,inventory__) * (alpha1) * (retailPrice[period-1] - salvageValue)
                                              - Math.min(0,inventory__) *(1-alpha1)  * (alpha) * (retailPrice[period-2] - salvageValue));


                            }


                            else {
                                z_line__ = BigM;
                            }
                            total__ += z_line__;
                            //    System.out.println("Total: " +  i+"  "+jj+"  "+total);
                        }
                        total1__ += total__;
                        //    System.out.println("Total: " +  i+"  "+jj+"  "+total);
                    }
                        total2__ +=total1__;

                    }
                    total3__ +=total2__;
                }
                total4__ +=total3__;
                }
                //          System.out.println("Total: " +  i+" " +pp+"  "+"  "+total2);
                if (total4__< min__) {
                    min__ = total4__;
                    SS__ = i__;
                }

            }





          //  Z_bigSDouble[period-1] = SS__ - MM;

            Z_bigSDouble[period-1] = 2;
*/





            Z_bigSDouble[period-1] =  (int)Z_bigSDoubleee[1][I[period-1] + M][SL[period-2]][SNR[period-1]] ;


            O[period-1] = (int)Math.max(0, (Z_bigSDouble[period-1]- I[period-1]));


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

    public static double[][][][] periodMatrixmodifiedv4_2(int [] II , int [] S,int period, double[][][] ZN,
                                                          int n, int sl, int s, int m, double lambda, int unitCost, int fixedCost, int [] retailPrice,
                                                          double holdingCost, double backorderCost,int BigM, int MM, int M,
                                                          double[] poissonDistProbability,double[][] binomialDistProbability1,double[][] binomialDistProbability2,
                                                          int p) {

        double total;
        double total1;
        double total2;

        double[][][][] Z_bigSDouble = new double[2][n][sl][sl];
        int[][][] SS = new int[n][sl][sl];
        double z_line;
        //----------------------------------------


        for (int inv = 0; inv < n; inv++) { //Inventory

            for (int j1 = 0; j1 < sl; j1++) { //1 period ago sold previous sales

                for (int j2 = 0; j2 < sl; j2++) { //2 period ago sold but not returned yet
                    double min = 100000;
                    for (int i = 0; i < s; i++) { //S up to order
                        total2 = 0;
                        for (int jj = 0; jj <= m; jj++) { //demand
                            total1 = 0;
                            for (int k1 = 0; k1 <= j1; k1++) { //return1
                                total = 0;
                                for (int k2 = 0; k2 <= j2; k2++) { //return2


                                    if (II[inv] < S[i]) {
                                        z_line = (binomialDistProbability1[j1][k1])
                                                *(binomialDistProbability2[j2][k2])
                                                *(poissonDistProbability[jj])

                                                *((S[i] - II[inv]) * unitCost + fixedCost
                                                + Math.max(0,(S[i] + k1 + k2 - jj))* holdingCost
                                                + Math.max(0,(jj - S[i] - k1 - k2))* backorderCost
                                                + k1 * retailPrice[p-1] + k2 * retailPrice[p-2]
                                                + ZN[S[i] + k1 + k2 - jj + M][Math.min(

                                                (Math.max(0, II[inv]) + Math.max(0, (S[i] - II[inv])) + k1 + k2 ),

                                                (jj - (Math.min(0, II[inv]))))][j1-k1]);
                                    }  else if (II[inv] == S[i]) {
                                        z_line = (binomialDistProbability1[j1][k1])
                                                *(binomialDistProbability2[j2][k2])
                                                *(poissonDistProbability[jj])

                                                *(Math.max(0,(S[i] + k1 + k2 - jj))* holdingCost
                                                + Math.max(0,(jj - S[i] - k1 - k2))* backorderCost
                                                + k1 * retailPrice[p-1] + k2 * retailPrice[p-2]
                                                + ZN[S[i] + k1 + k2 - jj + M][Math.min(

                                                (Math.max(0, II[inv]) + Math.max(0, (S[i] - II[inv])) + k1 + k2 ),

                                                (jj - (Math.min(0, II[inv]))))][j1-k1]);
                                    } else {
                                        z_line = BigM;
                                    }



                                    total += z_line;
                                }
                                total1 +=total;
                            }

                            total2 +=total1;


                        }


                        if (total2< min) {


                            min = total2;
                            SS[inv][j1][j2] = i;
                        }
                        //result2[inv][j][i] = total2;
                    }
                    Z_bigSDouble[0][inv][j1][j2] = min;
                    Z_bigSDouble[1][inv][j1][j2] = SS[inv][j1][j2] - MM;


                }
            }
        }
        return Z_bigSDouble;
    }


    public static double[][][] ZEndcalculate(int [] II, int period, int n, int sl, int m, double alpha1, double alpha2,double alpha,
                                             double backorderCost,
                                             int salvageValue, int [] retailPrice, int fixedCost, int unitCost,
                                             double[][] binomialDistProbability1, double[][] binomialDistProbability2) {


        double[][][] ZEnd = new double[n][sl][sl];


        double total;
        double total1;


        double z_line;

        for (int inv = 0; inv < n; inv++) { //Inventory

            for (int j1 = 0; j1 < sl; j1++) { //1 period ago sold previous sales

                for (int j2 = 0; j2 < sl; j2++) { //2 period ago sold but not returned yet

                    total1 = 0;
                    for (int k1 = 0; k1 <= j1; k1++) { //return1
                        total = 0;
                        for (int k2 = 0; k2 <= j2; k2++) { //return2

                            z_line = (binomialDistProbability1[j1][k1])
                                    * (binomialDistProbability2[j2][k2])
                                    * ( Math.min(0,II[inv] + k1 + k2) * -unitCost
                                    + Math.max(0,II[inv] + k1 + k2) * -salvageValue
                                    + k1 * retailPrice[period-1] + k2 * retailPrice[period-2]
                                    + (j1-k1)*alpha * (retailPrice[period-1] - salvageValue)
                                    - Math.min(0,II[inv]) * (alpha1) * (retailPrice[period-1] - salvageValue)
                                    - Math.min(0,II[inv]) *(1-alpha1)  * (alpha) * (retailPrice[period-2] - salvageValue));
                            total += z_line;
                            //  total += 0;
                        }
                        total1 +=total;
                    }


                    ZEnd[inv][j1][j2] = total1;

                }
            }

        }




        return  ZEnd;

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
