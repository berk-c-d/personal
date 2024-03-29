import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

//Last Updated 13.12.2022
public class multiperiod_1 {

    public static void main(String[] args) throws FileNotFoundException {
        Instant inst1 = Instant.now();
        int BigM=10000; int nresult =46;
        int nExperiment = 1;


        for (int i = 0; i < nExperiment; i++) {

            System.out.println(" Experiment: " + i);
            /*
            //Basecase 0.5	2	20	25	10	4
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 25;
            int m = 10; //max demand
            double lambda_all=4;
            String case_number="base";


            //Case-1 0.4	2	20	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="1";


            //Case-2  0.4	2.4	20	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 20;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="2";


            //Case-3  	0.4	2	32	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 32;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="3";


            //Case-4 0.4	2.4	32	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 32;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="4";


            //Case-5  0.4	2	50	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 50;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="5";


            //Case-6  0.5	2	50	40	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 50;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="6";


            //Case-7 0.4	2.4	50	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 50;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="7";


            //Case-8 0.4	2	20	50	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 50;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="8";


            //Case-9  0.4	2.4	20	50	25	8
            double alpha_all=0.4;
            double holdingCost = 2.4;
            double backorderCost = 20;
            int fixedCost = 50;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="9";


            //Case-10  0.5	2	20	100	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="10";


            //Case-11  0.5	2	20	45	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 45;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="11";


            //Case-12  0.4	2	50	80	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 50;
            int fixedCost = 80;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="12";


            //Case-13  0	2	20	100	25	8
            double alpha_all=0;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="13";


            //Case-14  0.1	2	20	100	25	8
            double alpha_all=0.1;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="14";

            //Case-15 0.4	2	5	40	25	8
            double alpha_all=0.4;
            double holdingCost = 2;
            double backorderCost = 5;
            int fixedCost = 40;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="15";

             //Case-16 0.5	2	5	100	25	8
            double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 5;
            int fixedCost = 100;
            int m = 25; //max demand
            double lambda_all=8;
            String case_number="16";




 */


            //Basecase 0.5	2	20	25	10	4
            //  double alpha_all=0.5;
            double holdingCost = 2;
            double backorderCost = 20;
            int fixedCost = 25;
            int m = 6; //max demand
            double lambda_all=3;
            String case_number="base";
            int period =4; //number of periods starting from 0, period is also used for terminal period  n=5
            double alpha1_all=0.1;
            double alpha2_all=0.4;


            double alpha1 = alpha1_all;

            double alpha2 = alpha2_all;

            double alpha = alpha2_all/(1-alpha1_all);

            double alpha1Simulation = alpha1_all;

            double alpha2Simulation = alpha2_all;

            double alphaSimulation  = alpha2Simulation/(1-alpha1Simulation);

            int unitCost = 20;

            int salvageValue = 0;  ///normalde 0 dı 0 ken arrayler repeat edebilirrr!!!!!!!!!!!

            int returnCredit1 = 80; //geri dönen ürün fiyatı
            int returnCredit2 = 80; //geri dönen ürün fiyatı

            //     double alpha = alpha_all; //return rate

            double lambda = lambda_all; //demand rate

            //   double alphaSimulation = alpha_all;

            double lambdaSimulation =lambda_all;

            int period_for_calculation=4;
            //--------------Inventory-----------------------------------------------------------------------------------------------------
            IntStream streamI = IntStream.range(-(period_for_calculation * m),((period_for_calculation)*4*m ));
            int M =period_for_calculation * m; //neutralizes negative I
            int[] I = streamI.toArray();
            int n = I.length;
            //-----------------------------------------------------------------------------------------------------------------------------

            //--------------Inventory level after ordering---------------------------------------------------------------------------------
            IntStream streamS = IntStream.range(-((period_for_calculation-1) * m), ((period_for_calculation)*(m))+1); //!!!!!!!!!!!!!son sayıdaki bigm bunu 1 arttırınca düzeldi!!!!!!!!!!!!!!!!!!!!!!!!!!
            int MM = (period_for_calculation-1) * m; //neutralizes negative S
            int[] S = streamS.toArray();
            int s = S.length;
            //------------------------------------------------------------------------------------------------------------------------------

            int sl =((period_for_calculation+1) * m)+1; //max # of sales  ////!!!!!bu değer arttırıldığı zaman possible conditionların oluştuğu dosyadaki en altta bulunan bigm ler ortadan kayboldu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!çünkü binomdistdeki j indexi dolaylı yoldan etkiliyor

            double[] poissonDistProbability = new double[m+1];
            poissonDistProbability = poissonDist (lambda, m);  //Max demand gönderiliyor çünkü 0dan bu değere kadar demand gelme olasılıkları hesaplanıp 1D arraye kaydediliyor.

            double[][] binomialDistProbability1 = new double[sl][sl+1];
            binomialDistProbability1 = binomialDist(sl, alpha1);

            double[][] binomialDistProbability2 = new double[sl][sl+1];
            binomialDistProbability2 = binomialDist(sl, alpha);


            double[][][][][] Z_bigSDouble = new double[period_for_calculation+1][2][n][sl][sl];
            Z_bigSDouble[period][0] = ZEndcalculate(I, period, n, sl, m, alpha1, alpha2,alpha, backorderCost,
                    salvageValue, returnCredit1, returnCredit2, fixedCost, unitCost,
                    binomialDistProbability1, binomialDistProbability2);

            PrintWriter Solutionnn= new PrintWriter("opt_progression_for_case_"+case_number+".csv");
            Solutionnn.println("Order_up_to,Demand,Inventory, prev_sale,total,total2");

            for (int p = period - 1; p >= 2; p--) {
                Z_bigSDouble[p] = periodMatrixmodifiedv4_2(I , S, period,  Z_bigSDouble[p+1][0], n, sl, s, m, lambda, unitCost, fixedCost,
                		returnCredit1, returnCredit2, holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,
                		binomialDistProbability1, binomialDistProbability2, Solutionnn, p);
                }

                Z_bigSDouble[1] = periodMatrixmodifiedv4_2_period1(I , S, period,  Z_bigSDouble[2][0], n, sl, s, m, lambda, unitCost, fixedCost,
                		returnCredit1, returnCredit2, holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,
                		binomialDistProbability1, binomialDistProbability2, Solutionnn, 1);


                Z_bigSDouble[0] = periodMatrixmodifiedv4_2_period0(I , S, period,  Z_bigSDouble[1][0], n, sl, s, m, lambda, unitCost, fixedCost,
                        returnCredit1, returnCredit2, holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,
                        binomialDistProbability1, binomialDistProbability2, Solutionnn, 0);

            Solutionnn.close();



            write_sols_for_Z_bigSDouble_v5_multi(period,alpha1,alpha2,alpha,lambda,m, period,"Write File Name and fix method",0,Z_bigSDouble[0][0].length,0,Z_bigSDouble[0][0][0].length,0,Z_bigSDouble[0][0][0][0].length,Z_bigSDouble,M,MM);
            write_sols_for_Z_bigSDouble_v5_multi_pos(period,alpha1,alpha2,alpha,lambda,m, period,"POS_COST_hold_"+holdingCost+"backo_"+backorderCost+"fix_"+fixedCost+"case_number_"+case_number,0,Z_bigSDouble[0][0].length,0,Z_bigSDouble[0][0][0].length,0,Z_bigSDouble[0][0][0][0].length,Z_bigSDouble,M,MM);
            //-------------------------------------------------NAIVE--------------------------------------------------------------------------------

            double OP = Z_bigSDouble[0][0][M][0][0];

            double[] arr = new double[nresult];
            arr = SimulationMatrix_v3adapted(period, Z_bigSDouble, OP, m, alpha1Simulation, alpha2Simulation, alphaSimulation,
                    lambdaSimulation, unitCost, fixedCost, returnCredit1, returnCredit2, holdingCost, backorderCost,
                    salvageValue, M, nresult);

            //Experiment[i] = arr;
            //PrintWriter Solution = new PrintWriter("max_per_"+String.valueOf(period)+"__maxdem_"+String.valueOf(m)+"__alpha_"+String.valueOf(alphaSimulation)+"__lambda_"+String.valueOf(lambdaSimulation)+"_Simulation_result.txt");
            PrintWriter Solution = new PrintWriter("SIMULATION RESULT FOR CASE NUMBER "+case_number+".csv");
            Solution.println("tc, pc, oc, hc, bc,tec, rc,SLV,tc, O[0],O[1],O[2],O[3],O[4],bigSValue[0], bigSValue[1], bigSValue[2], bigSValue[3], bigSValue[4], BO[0],BO[1],BO[2],BO[3],BO[4], HI[0],HI[1],HI[2],HI[3],HI[4],R[1],R[2],R[3],R[4],R[5],SL[0],SL[1],SL[2],SL[3],SL[4],SL[5], I[0],I[1],I[2],I[3],I[4],I[5]");
            Solution.println(Arrays.toString(arr));

            System.out.println(Arrays.toString(arr));
            //System.out.println(Arrays.deepToString(Z_bigSDouble[period][0]));
            System.out.println(OP);

            Solution.close();
        }
        //Solution.close();


        Instant inst2 = Instant.now();
        System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
    }
    public static void write_sols_for_Z_bigSDouble_v2(int maxperiod,double alphaval,double lambdaval, int maxdemand,int period,String filename,int inventory_range_low,int inventory_range_high,int order_up_to_range_low,int order_up_to_range_high,double[][][][] Z_bigSDouble,int M, int MM ) throws FileNotFoundException {
        double[] resultt = new double[3];
        for (int p =0; p <=period; p++) {
            String ssss=String.valueOf(p);
            String ssss1=String.valueOf(inventory_range_low-M);
            String ssss2=String.valueOf(inventory_range_high-M);
            String ssss3=String.valueOf(order_up_to_range_low-MM);
            String ssss4=String.valueOf(order_up_to_range_high-MM);
            String maxperiods=String.valueOf(maxperiod);
            String alphavals=String.valueOf(alphaval);
            String lambdavals=String.valueOf(lambdaval);
            String maxdemands=String.valueOf(maxdemand);
            PrintWriter Solutionn= new PrintWriter("maxper_"+maxperiods+"__maxdem_"+maxdemands+"__alpha_"+alphavals+"__lambda_"+lambdavals+"-"+filename+" for period "+ssss+"__Inventory range "+ssss1+" to"+ssss2+"__Orderupto range "+ssss3+" to "+ssss4+".csv");
            Solutionn.println("Inventory,Previous_sale,Profit");
            for (int inv = inventory_range_low; inv < inventory_range_high; inv++) {
                for (int j = order_up_to_range_low; j <  order_up_to_range_high; j++) {
                    resultt=new double[]{inv-M,j,Z_bigSDouble[p][0][inv][j]};
                    Solutionn.println(Arrays.toString(resultt).replace("[","").replace("]",""));
                }
            }
            Solutionn.close();
        }


    }
    public static void write_sols_for_Z_bigSDouble_v3(int maxperiod,double alphaval,double lambdaval, int maxdemand,int period,String filename,int inventory_range_low,int inventory_range_high,int order_up_to_range_low,int order_up_to_range_high,double[][][][] Z_bigSDouble,int M, int MM) throws FileNotFoundException {
        double[] resultt = new double[3];
        for (int p =0; p <=period; p++) {
            String ssss=String.valueOf(p);
            String ssss1=String.valueOf(inventory_range_low-M);
            String ssss2=String.valueOf(inventory_range_high-M);
            String ssss3=String.valueOf(order_up_to_range_low-MM);
            String ssss4=String.valueOf(order_up_to_range_high-MM);
            String maxperiods=String.valueOf(maxperiod);
            String alphavals=String.valueOf(alphaval);
            String lambdavals=String.valueOf(lambdaval);
            String maxdemands=String.valueOf(maxdemand);
            PrintWriter Solutionn= new PrintWriter("maxper_"+maxperiods+"__maxdem_"+maxdemands+"__alpha_"+alphavals+"__lambda_"+lambdavals+"-"+filename+" for period "+ssss+"__Inventory range "+ssss1+" to"+ssss2+"__Orderupto range "+ssss3+" to "+ssss4+".csv");
            Solutionn.println("Inventory,Previous_sale,Cost");

            //  for (int inv = -maxdemand*p+M; inv < maxdemand*maxperiod+1+M; inv++) {
            //   for (int j = 0; j < (p*maxdemand)+1; j++) {
            //  if((   (inv<M)    && ((((inv-M)*-1)+j)<=p*maxdemand) )  || (p==0 && inv-M==0 && j==0)  ||  ((p!=0) &&  (inv>=M)   &&  (inv-M+j<=maxperiod*maxdemand) ) ) {
            //  resultt = new double[]{inv - M, j , Z_bigSDouble[p][0][inv][j]};
            // Solutionn.println(Arrays.toString(resultt).replace("[", "").replace("]", ""));
            //}}}
            for (int inv = inventory_range_low; inv < inventory_range_high; inv++) {
                for (int j = order_up_to_range_low; j < order_up_to_range_high; j++) {
                    if((   (inv<M)    && ((((inv-M)*-1)+j)<=p*maxdemand) )  || (p==0 && inv-M==0 && j==0)  ||  ((p!=0) &&  (inv>=M)   &&  (inv-M+j<=maxperiod*maxdemand) && j<=p*maxdemand ) ) {
                        resultt = new double[]{inv - M, j , Z_bigSDouble[p][0][inv][j]};
                        Solutionn.println(Arrays.toString(resultt).replace("[", "").replace("]", ""));
                    }}
            }
            Solutionn.close();
        }


    }
    public static void write_sols_for_Z_bigSDouble_v4(int maxperiod,double alphaval,double lambdaval, int maxdemand,int period,String filename,int inventory_range_low,int inventory_range_high,int order_up_to_range_low,int order_up_to_range_high,double[][][][] Z_bigSDouble,int M, int MM) throws FileNotFoundException {
        double[] resultt = new double[3];
        for (int p =0; p <=period; p++) {
            String ssss=String.valueOf(p);
            String ssss1=String.valueOf(inventory_range_low-M);
            String ssss2=String.valueOf(inventory_range_high-M);
            String ssss3=String.valueOf(order_up_to_range_low-MM);
            String ssss4=String.valueOf(order_up_to_range_high-MM);
            String maxperiods=String.valueOf(maxperiod);
            String alphavals=String.valueOf(alphaval);
            String lambdavals=String.valueOf(lambdaval);
            String maxdemands=String.valueOf(maxdemand);
            PrintWriter Solutionn= new PrintWriter("maxper_"+maxperiods+"_maxdem_"+maxdemands+"_alpha_"+alphavals+"_lambda_"+lambdavals+"-"+filename+"period "+ssss+"__Inventory_"+ssss1+"to"+ssss2+"_prevs_"+ssss3+" to "+ssss4+".csv");
            Solutionn.println("Inventory,Previous_sale,Cost,Break_Point,Break_Point_pos,Inv_brk_pos,Inv_brk_neg,Inv_brk_pos_sqr,Inv_brk_neg_sqr");
            double break_p=0;
            double break_p_pos;
            double inv_min_brk_pos;
            double inv_min_brk_neg;
            for (int inv = inventory_range_low; inv < inventory_range_high; inv++) {
                for (int j = order_up_to_range_low; j < order_up_to_range_high; j++) {
                    if((   (inv<M)    && ((((inv-M)*-1)+j)<=p*maxdemand) )  || (p==0 && inv-M==0 && j==0)  ||  ((p!=0) &&  (inv>=M)   &&  (inv-M+j<=maxperiod*maxdemand) && j<=p*maxdemand ) ) {
                        break_p=lambdaval*(maxperiod-p)-alphaval*j;
                        if(break_p<=0){break_p_pos=0;}
                        else{break_p_pos=break_p;}

                        if(inv-M-break_p_pos>=0){inv_min_brk_pos=inv-M-break_p_pos;
                            inv_min_brk_neg=0;}
                        else{inv_min_brk_neg=(inv-M-break_p_pos)*-1;
                            inv_min_brk_pos=0;
                        }
                        resultt = new double[]{inv - M, j , Z_bigSDouble[p][0][inv][j],break_p,break_p_pos,inv_min_brk_pos,inv_min_brk_neg,inv_min_brk_pos*inv_min_brk_pos,inv_min_brk_neg*inv_min_brk_neg};
                        Solutionn.println(Arrays.toString(resultt).replace("[", "").replace("]", ""));
                    }}
            }
            Solutionn.close();
        }


    }
    public static void write_sols_for_Z_bigSDouble_v5_multi(int maxperiod,double alphaval1,double alphaval2,double alphaval3,double lambdaval, int maxdemand,int period,String filename,int inventory_range_low,int inventory_range_high,int prevsale_1_range_low,int prevsale_1_range_high,int prevsale_2_range_low,int prevsale_2_range_high,double[][][][][] Z_bigSDouble,int M, int MM) throws FileNotFoundException {
        double[] resultt = new double[3];
        for (int p =0; p <=period; p++) {
            String ssss=String.valueOf(p);/*
            String ssss1=String.valueOf(inventory_range_low-M);
            String ssss2=String.valueOf(inventory_range_high-M);
            String ssss3=String.valueOf(prevsale_1_range_low-MM);
            String ssss4=String.valueOf(prevsale_1_range_high-MM);
            String maxperiods=String.valueOf(maxperiod);
            String alphavals=String.valueOf(alphaval1);
            String lambdavals=String.valueOf(lambdaval);
            String maxdemands=String.valueOf(maxdemand);

            */
            //PrintWriter Solutionn= new PrintWriter("maxper_"+maxperiods+"_maxdem_"+maxdemands+"_alpha_"+alphavals+"_lambda_"+lambdavals+"-"+filename+"period "+ssss+"__Inventory_"+ssss1+"to"+ssss2+"_prevs_"+ssss3+" to "+ssss4+".csv");
            PrintWriter Solutionn= new PrintWriter("Allcond_multiperiod_for_period_"+ssss+".csv");
            Solutionn.println("Inventory,Previous_sale1,Previous_sale2,Cost,Break_Point,Break_Point_pos,Inv_brk_pos,Inv_brk_neg,Inv_brk_pos_sqr,Inv_brk_neg_sqr");
            double break_p=0;
            double break_p_pos;
            double inv_min_brk_pos;
            double inv_min_brk_neg;
            for (int inv = inventory_range_low; inv < inventory_range_high; inv++) {
                for (int j1 = prevsale_1_range_low; j1 < prevsale_1_range_high; j1++) {
                    for (int j2 = prevsale_2_range_low; j2 < prevsale_2_range_high; j2++) {
                    //if((   (inv<M)    && ((((inv-M)*-1)+j)<=p*maxdemand) )  || (p==0 && inv-M==0 && j==0)  ||  ((p!=0) &&  (inv>=M)   &&  (inv-M+j<=maxperiod*maxdemand) && j<=p*maxdemand ) ) {
                        //----------------------------------------------------------------------------
                        break_p=lambdaval*(maxperiod-p)-alphaval1*j1;//not completeeeeeeeeeeeeee

                        if(break_p<=0){break_p_pos=0;}
                        else{break_p_pos=break_p;}

                        if(inv-M-break_p_pos>=0){inv_min_brk_pos=inv-M-break_p_pos;
                            inv_min_brk_neg=0;}
                        else{inv_min_brk_neg=(inv-M-break_p_pos)*-1;
                            inv_min_brk_pos=0;
                        }
                        //----------------------------------------------------------------------------
                        resultt = new double[]{inv - M, j1,j2 , Z_bigSDouble[p][0][inv][j1][j2],break_p,break_p_pos,inv_min_brk_pos,inv_min_brk_neg,inv_min_brk_pos*inv_min_brk_pos,inv_min_brk_neg*inv_min_brk_neg};
                        Solutionn.println(Arrays.toString(resultt).replace("[", "").replace("]", ""));
                    //}
                    }}}
            Solutionn.close();
        }


    }

    public static void write_sols_for_Z_bigSDouble_v5_multi_pos(int maxperiod,double alphaval1,double alphaval2,double alphaval3,double lambdaval, int maxdemand,int period,String filename,int inventory_range_low,int inventory_range_high,int prevsale_1_range_low,int prevsale_1_range_high,int prevsale_2_range_low,int prevsale_2_range_high,double[][][][][] Z_bigSDouble,int M, int MM) throws FileNotFoundException {
        double[] resultt = new double[3];
        for (int p =0; p <=period; p++) {
            String ssss=String.valueOf(p);/*
            String ssss1=String.valueOf(inventory_range_low-M);
            String ssss2=String.valueOf(inventory_range_high-M);
            String ssss3=String.valueOf(prevsale_1_range_low-MM);
            String ssss4=String.valueOf(prevsale_1_range_high-MM);
            String maxperiods=String.valueOf(maxperiod);
            String alphavals=String.valueOf(alphaval1);
            String lambdavals=String.valueOf(lambdaval);
            String maxdemands=String.valueOf(maxdemand);

            */
            //PrintWriter Solutionn= new PrintWriter("maxper_"+maxperiods+"_maxdem_"+maxdemands+"_alpha_"+alphavals+"_lambda_"+lambdavals+"-"+filename+"period "+ssss+"__Inventory_"+ssss1+"to"+ssss2+"_prevs_"+ssss3+" to "+ssss4+".csv");
            PrintWriter Solutionn= new PrintWriter("POScond_multiperiod_for_period_"+ssss+".csv");
            Solutionn.println("Inventory,Previous_sale1,Previous_sale2,Cost,Break_Point,Break_Point_pos,Inv_brk_pos,Inv_brk_neg,Inv_brk_pos_sqr,Inv_brk_neg_sqr");
            double break_p=0;
            double break_p_pos;
            double inv_min_brk_pos;
            double inv_min_brk_neg;
            for (int inv = inventory_range_low; inv < inventory_range_high; inv++) {
                for (int j1 = prevsale_1_range_low; j1 < prevsale_1_range_high; j1++) {
                    for (int j2 = prevsale_2_range_low; j2 < prevsale_2_range_high; j2++) {

                        if(     (p==0 && inv-M==0 && j1==0 && j2==0) || //For Period 0

                                (p==1  && j2==0 && ((inv<M)  && ((((inv-M)*-1)+j1)<=p*maxdemand) ) )|| //For Period 1
                                (p==1  && j2==0 && (inv>=M)   &&  (inv-M+j1<=maxperiod*maxdemand) && j1<=p*maxdemand  && ((inv-M+j1+j2)<=maxperiod*maxdemand)) || //For Period 1

                                ((p!=1) && (inv<M)  && ((((inv-M)*-1)+j1)<=p*maxdemand) && (j2<=(p-1)*maxdemand)&& ((((inv-M)*-1)+j1+j2)<=maxperiod*maxdemand) ) ||  //
                                ((p!=1) &&  (p!=0) &&  (inv>=M)   &&  (inv-M+j1<=maxperiod*maxdemand) && (j1<=p*maxdemand) && (j2<=(p-1)*maxdemand)  && ((inv-M+j1+j2)<=maxperiod*maxdemand) ) ) {
                        //----------------------------------------------------------------------------
                        break_p=lambdaval*(maxperiod-p)-alphaval1*j1-alphaval2*j2;//not completeeeeeeeeeeeeee

                        if(break_p<=0){break_p_pos=0;}
                        else{break_p_pos=break_p;}

                        if(inv-M-break_p_pos>=0){inv_min_brk_pos=inv-M-break_p_pos;
                            inv_min_brk_neg=0;}
                        else{inv_min_brk_neg=(inv-M-break_p_pos)*-1;
                            inv_min_brk_pos=0;
                        }
                        //----------------------------------------------------------------------------
                        resultt = new double[]{inv - M, j1,j2 , Z_bigSDouble[p][0][inv][j1][j2],break_p,break_p_pos,inv_min_brk_pos,inv_min_brk_neg,inv_min_brk_pos*inv_min_brk_pos,inv_min_brk_neg*inv_min_brk_neg};
                        Solutionn.println(Arrays.toString(resultt).replace("[", "").replace("]", ""));
                        }
                    }}}
            Solutionn.close();
        }


    }

    public static double[] SimulationMatrix_v3adapted(int period, double[][][][][] Z_bigSDouble, double OP, int m,
                                                      double alpha1Simulation, double alpha2Simulation, double alphaSimulation, double lambdaSimulation,
                                                      int unitCost, int fixedCost, int returnCredit1,int returnCredit2, double holdingCost,
                                                      double backorderCost, int salvageValue, int M, int nresult) {
        int n = 100000;

        double[] arr = new double[nresult];
        double[] arr_sum = new double[nresult];
        double[] arr_mean = new double[nresult];

        //double[][] arr2 = new double[n][nresult];

        for (int i = 0; i < n; i++) {
            arr = ProcessMatrix_v3adapted(period, Z_bigSDouble, OP, m, alpha1Simulation, alpha2Simulation, alphaSimulation,
                    lambdaSimulation, unitCost, fixedCost, returnCredit1, returnCredit2, holdingCost, backorderCost,
                    salvageValue, M, nresult);

            for(int j=0;j<arr.length;j++){
                arr_sum[j]+=arr[j];

            }

        }

        for(int j=0;j<arr.length;j++){
            arr_mean[j]=arr_sum[j]/n;

        }


        //arr = meanOfSecondInstance(arr2, nresult);
        return arr_mean;

    }

    public static double[] ProcessMatrix_v3adapted(int period, double[][][][][] Z_bigSDouble, double OP, int m,
                                                   double alpha1Simulation, double alpha2Simulation, double alphaSimulation, double lambdaSimulation,
                                                   int unitCost, int fixedCost, int returnCredit1,int returnCredit2, double holdingCost,
                                                   double backorderCost, int salvageValue, int M, int nresult) {

        double[] result = new double[nresult];

        int[] I = new int[period + 1];

        int[] R1 = new int[period + 1];

        int[] R2 = new int[period + 1];

        int[] O = new int[period];

        int[] D = new int[period];

        int[] BO = new int[period];

        int[] IO = new int[period + 1];

        int[] SL = new int[period + 1];

        int[] HI = new int[period];

        int[] SNR = new int[period+2];

        int[] BigS = new int[period];

        for (int p = 0; p < period; p++) {

            D[p] = getPoissonRandom(lambdaSimulation);

            if (D[p] > m) {

                D[p] = m;

            }

        }

        I[0] = 0;

        O[0] = (int)Math.max(0, Z_bigSDouble[0][1][M][0][0]);

        BigS[0]= (int)Z_bigSDouble[0][1][M][0][0];

        BO[0] = Math.max(0, D[0] - I[0] - O[0] );

        IO[0] = Math.max(0, I[0]);

        SL[0] = Math.min(IO[0] + O[0], D[0]);

        HI[0] = Math.max(0, I[0] + O[0] - D[0]);

        R1[0] = 0;


        I[1] = I[0] + O[0] + R1[0] - D[0];

        R1[1] = getBinomialRandom(SL[0], alpha1Simulation);

        R2[1] = 0;

        O[1] = (int)Math.max(0, (Z_bigSDouble[1][1][I[1] + M][SL[0]][0] - I[1]));

        BigS[1]= (int)Z_bigSDouble[1][1][I[1] + M][SL[0]][0];

        BO[1] = Math.max(0, D[1] - I[1] - O[1] - R1[1] - R2[1]);

        IO[1] = Math.max(0, I[1]);

        SL[1] = Math.min(IO[1] + O[1] + R1[1] + R2[1], D[1] + BO[0]);

        HI[1] = Math.max(0, I[1] + O[1] + R1[1] + R2[1] - D[1]);


        SNR[1]=0;


        for (int i = 2; i < period; i++) {

            I[i] = I[i-1] + O[i-1] + R1[i-1] + R2[i-1] - D[i-1];

            R1[i] = getBinomialRandom(SL[i-1], alpha1Simulation);

            SNR[i]=SL[i-2]-R1[i-1];

            R2[i] = getBinomialRandom(SNR[i], alphaSimulation);

            O[i] = (int)Math.max(0, (Z_bigSDouble[i][1][I[i] + M][SL[i-1]][SNR[i]] - I[i]));

            BigS[i]= (int)Z_bigSDouble[i][1][I[i] + M][SL[i-1]][SNR[i]];

            BO[i] = Math.max(0, D[i] - I[i] - O[i] - R1[i]- R2[i]);

            IO[i] = Math.max(0, I[i]);

            SL[i] = Math.min(IO[i] + O[i] + R1[i] + R2[i], D[i] + BO[i-1]);

            HI[i] = Math.max(0, I[i] + O[i] + R1[i] + R2[i] - D[i]);


        }

        I[period] = I[period-1] + O[period-1] + R1[period-1] + R2[period-1] - D[period-1];

        R1[period] = getBinomialRandom(SL[period-1], alpha1Simulation);

        SNR[period]=SL[period-2]-R1[period-1];

        R2[period] = getBinomialRandom(SNR[period], alphaSimulation);


        double procurementCost = 0;

        double totalOrderingCost = 0;

        double totalHoldingCost = 0;

        double totalBackorderCost = 0;

        double returnCost = 0;

        double totalEndingInventoryCost = 0;




        for (int p = 0; p < period; p++) {

            int[] numberOfOrder = new int[period];

            procurementCost += O[p] * unitCost;

            if (O[p] > 0) {

                numberOfOrder[p] = 1;

            } else {

                numberOfOrder[p] = 0;

            }

            totalOrderingCost += numberOfOrder[p] * fixedCost;

            totalHoldingCost += HI[p] * holdingCost;

            totalBackorderCost += BO[p] * backorderCost;

            SL[period]= -(Math.min(0,I[period]));
            SNR[period +1]=SL[period-1]-R1[period];

            returnCost += R1[p]* returnCredit1 + R2[p]* returnCredit2;

        }

        totalEndingInventoryCost =
                + Math.min(0,(I[period] + R1[period] + R2[period])) * -unitCost
                        + Math.max(0,(I[period] + R1[period] + R2[period] )) * -salvageValue // + SNR[period +1] *alphaSimulation
                        + R1[period] * returnCredit1 + R2[period] * returnCredit2
                        + SNR[period +1]*alphaSimulation*(returnCredit2 - salvageValue)
                - Math.min(0,I[period]) * alpha1Simulation  * (returnCredit1 - salvageValue)
                - Math.min(0,I[period]) * (1-alpha1Simulation) * (alphaSimulation)  * (returnCredit2 -salvageValue);

        double salvage=I[period]+R1[period]+R2[period];
        double SLV= salvage;

        double totalCost = procurementCost + totalOrderingCost + totalHoldingCost + totalBackorderCost

                + totalEndingInventoryCost + returnCost;


        double pc = procurementCost;

        double oc = totalOrderingCost;

        double hc = totalHoldingCost;

        double bc = totalBackorderCost;

        double tec = totalEndingInventoryCost;

        double rc = returnCost;

        double tc = totalCost;



        result = new double[] {OP,tc, pc, oc, hc, bc,tec, rc, SLV,tc,BigS[0],BigS[1],BigS[2],BigS[3]};


        return result;

    }
    public static double[][][][] periodMatrixmodifiedv4_2(int [] I , int [] S,int period, double[][][] ZN,
                                                          int n, int sl, int s, int m, double lambda, int unitCost, int fixedCost, int returnCredit1,
                                                          int returnCredit2, double holdingCost, double backorderCost,int BigM, int MM, int M,
                                                          double[] poissonDistProbability,double[][] binomialDistProbability1,double[][] binomialDistProbability2,
                                                          PrintWriter Solutionnn,int p) {

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


                                    if (I[inv] <= S[i]) {
                                        z_line = (binomialDistProbability1[j1][k1])
                                                *(binomialDistProbability2[j2][k2])
                                                *(poissonDistProbability[jj])

                                                *((S[i] - I[inv]) * unitCost + Math.min(1, S[i]-I[i]) *fixedCost
                                                + Math.max(0,(S[i] + k1 + k2 - jj))* holdingCost
                                                + Math.max(0,(jj - S[i] - k1 - k2))* backorderCost
                                                + k1 * returnCredit1 + k2 * returnCredit2
                                                + ZN[S[i] + k1 + k2 - jj + M][Math.min(

                                                (Math.max(0, I[inv]) + Math.max(0, (S[i] - I[inv])) + k1 + k2),

                                                (jj - (Math.min(0, I[inv]))))][j1-k1]);

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

    public static double[][][][] periodMatrixmodifiedv4_2_period1(int [] I , int [] S,int period, double[][][] ZN,
                                                                  int n, int sl, int s, int m, double lambda, int unitCost, int fixedCost, int returnCredit1,
                                                                  int returnCredit2, double holdingCost, double backorderCost,int BigM, int MM, int M,
                                                                  double[] poissonDistProbability,double[][] binomialDistProbability1,double[][] binomialDistProbability2,
                                                                  PrintWriter Solutionnn,int p) {

        double total;

        double total2;

        double[][][][] Z_bigSDouble = new double[2][n][sl][sl];
        int[][] SS = new int[n][sl];
        double z_line;
        //----------------------------------------


        for (int inv = 0; inv < n; inv++) { //Inventory

            for (int j1 = 0; j1 < sl; j1++) { //1 period ago sold previous sales
                double min = 100000;
                for (int i = 0; i < s; i++) { //S up to order
                    total2 = 0;
                    for (int jj = 0; jj <= m; jj++) { //demand
                        total = 0;
                        for (int k1 = 0; k1 <= j1; k1++) { //return1



                            if (I[inv] <= S[i]) {
                                z_line = (binomialDistProbability1[j1][k1])

                                        *(poissonDistProbability[jj])

                                        *((S[i] - I[inv]) * unitCost + Math.min(1, S[i]-I[i]) *fixedCost
                                        + Math.max(0,(S[i] + k1 - jj))* holdingCost
                                        + Math.max(0,(jj - S[i] - k1 ))* backorderCost
                                        + k1 * returnCredit1
                                        + ZN[S[i] + k1  - jj + M][Math.min(

                                        (Math.max(0, I[inv]) + Math.max(0, (S[i] - I[inv])) + k1 ),

                                        (jj - (Math.min(0, I[inv]))))][j1-k1]);

                            } else {
                                z_line = BigM;
                            }

                            total += z_line;


                        }

                        total2 +=total;


                    }


                    if (total2< min) {


                        min = total2;
                        SS[inv][j1] = i;
                    }

                }
                Z_bigSDouble[0][inv][j1][0] = min;
                Z_bigSDouble[1][inv][j1][0] = SS[inv][j1] - MM;


            }
        }
        return Z_bigSDouble;
    }

    public static double[][][][] periodMatrixmodifiedv4_2_period0(int [] I , int [] S,int period, double[][][] ZN,
                                                                  int n, int sl, int s, int m, double lambda, int unitCost, int fixedCost, int returnCredit1,
                                                                  int returnCredit2, double holdingCost, double backorderCost,int BigM, int MM, int M,
                                                                  double[] poissonDistProbability,double[][] binomialDistProbability1,double[][] binomialDistProbability2,
                                                                  PrintWriter Solutionnn,int p) {

        double total;

        double[][][][] Z_bigSDouble = new double[2][n][sl][sl];
        int SS =0;
        double z_line;

        //----------------------------------------


        double min = 100000;
        for (int i = 0; i < s; i++) { //S up to order
            total = 0;
            for (int jj = 0; jj <= m; jj++) { //demand


                z_line = (poissonDistProbability[jj])

                        *(i * unitCost + Math.min(1, i) *fixedCost
                        + Math.max(0,(i - jj))* holdingCost
                        + Math.max(0,(jj - i ))* backorderCost
                        + ZN[i - jj + M][Math.min(i,jj)][0]);

                total += z_line;

            }


            if (total< min) {


                min = total;
                SS = i;
            }

        }
        Z_bigSDouble[0][M][0][0] = min;
        Z_bigSDouble[1][M][0][0] = SS;




        return Z_bigSDouble;
    }

    public static double[][][] ZEndcalculate(int [] I, int period, int n, int sl, int m, double alpha1, double alpha2, double alpha,double backorderCost,
                                             int salvageValue, int returnCredit1, int returnCredit2, int fixedCost, int unitCost,
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
                                    * ( Math.min(0,I[inv] + k1 + k2) * -unitCost
                                    + Math.max(0,I[inv] + k1 + k2) * -salvageValue
                                    + k1 * returnCredit1 + k2 * returnCredit2
                                    + (j1-k1)*alpha * (returnCredit2 - salvageValue)
                                    - Math.min(0,I[inv]) * (alpha1) * (returnCredit1 - salvageValue)
                                    - Math.min(0,I[inv]) *(1-alpha1)  * (alpha) * (returnCredit2 -salvageValue));
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
    private static int getPoissonRandomwebf(double mean) {
        Random r = new Random();
        double L = Math.exp(-mean);
        int k = 0;
        double p = 1.0;
        do {
            p = p * r.nextDouble();
            k++;
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

    public static double[][] binomialDist(int sl, double alpha) {

        double[][] P = new double[sl][sl];



        for (int j = 0; j <sl ; j++) {
            for (int k = 0; k <=j ; k++) {
                P[j][k]=nCr(j, k) * (double) Math.pow(alpha, k) * (double) Math.pow(1 - alpha, j - k);
            }
        }


        return P;

    }

    public static double[] poissonDist(double lambda,int m) {

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
