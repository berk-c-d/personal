import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.stream.IntStream;


public class ProfitPoisson_BackOrder_Revisited_FINAL_04052022_original {

	public static void main(String[] args) throws FileNotFoundException {
		Instant inst1 = Instant.now();

	int BigM=10000; int nresult =70;

	PrintWriter Solution0= new PrintWriter("Solution0_Z0_01092022_Noreturn.txt");
	PrintWriter Solution1= new PrintWriter("Solution1_Z1_01092022_Noreturn.txt");
	PrintWriter Solution2= new PrintWriter("Solution2_Z2_01092022_Noreturn.txt");
	PrintWriter Solution3= new PrintWriter("Solution3_Z3_01092022_Noreturn.txt");
	PrintWriter Solution4= new PrintWriter("Solution4_Z4_01092022_Noreturn.txt");
	PrintWriter Solution5= new PrintWriter("Solution5_Z5_01092022_Noreturn.txt");
	
		PrintWriter Solution = new PrintWriter("Solution_oneway10timesallparameter_06052022_rev.txt");

			int nExperiment = 1;
				
		

		double[][] Experiment = new double[nExperiment][nresult];
/*
		double[][] data = new double[nExperiment][15];


				Scanner input = new Scanner(new File("Input_parameter_06052022.txt"));

		for (int i = 0; i < nExperiment; i++) {

		for (int j = 0; j < 15; j++) {

		data[i][j] = input.nextDouble();

		}

		}

		input.close();
		*/
		//For batch sensitivity analysis, use this

		for (int i = 0; i < nExperiment; i++) {

		System.out.println(" Experiment: " + i);
		
		/*		int unitCost = (int) data[i][0];

		double holdingCost = (double) data[i][1];

		double backorderCost = (double) data[i][2];

		int salvageValue = (int) data[i][3];

		int fixedCost = (int) data[i][4];

		int returnCredit = (int) data[i][5];

		int retailPrice = (int) data[i][6];

		int period = (int) data[i][7];

		double alpha = (double) data[i][8];

		double lambda = (double) data[i][9];
		
		double alphaC = (double) data[i][10];

		double lambdaC = (double) data[i][11];
		
		double alphaSimulation = (double) data[i][12];

		double lambdaSimulation = (double) data[i][13];
		
		int m = (int) data[i][14];

	*/
	//Base case scenario parameter set:	
		int unitCost = 40;

		double holdingCost = 2;

		double backorderCost = 20;

		int salvageValue = 30;

		int fixedCost = 25;

		int returnCredit = 80;

		int retailPrice = 80;

		int period = 5; //number of periods starting from 0, period is also used for terminal period

		double alpha = 0.5; //return rate

		double lambda = 4; //demand rate
		
		double alphaC = 0;

		double lambdaC = 2;
		
		double alphaSimulation = 0.5;

		double lambdaSimulation = 4;

		int m = 8; //max demand
	

		
		IntStream streamI = IntStream.range(-(period * m), (period *2 *m)); //inventory
		int[] I = streamI.toArray();
		IntStream streamS = IntStream.range(-((period-1) * m), ((period-1)*m)); //BigS - UptoOrder level (inventory level after ordering)
		int[] S = streamS.toArray();
		int MM = (period-1) * m; //neutralizes negative S 
		int M =period * m; //neutralizes negative I
		int sl =((period+1) * m)+1; //max # of sales
		int n = I.length; //max # of inventory
		int s = S.length; //max # of up to order level (inventory level after ordering)
	
		// Two retailers (smart, naive) naive retailer is labeled with capital C in the phrases. (C-Conventional)
		
		double[] poissonDistProbability = new double[m+1];
		poissonDistProbability = poissonDist (lambda, m);

		double[][] binomialDistProbability = new double[sl][sl+1];
		binomialDistProbability = binomialDist(sl, alpha);
		
		double[] poissonDistProbabilityC = new double[m+1];
		poissonDistProbabilityC = poissonDist (lambdaC, m);

		double[][] binomialDistProbabilityC = new double[sl][sl+1];
		binomialDistProbabilityC = binomialDist(sl, alphaC);
		
		double[][][] Z = new double[period + 1][n][sl]; //Z denotes profit to go.
	
		//periods are denoted with p starting from 0, and "period" is the terminal period. 

		Z[period] = ZEndcalculate(I, period, n,sl, m, alpha, unitCost, salvageValue, retailPrice, returnCredit,

		fixedCost,unitCost, binomialDistProbability);
		
		double[][] J = new double[s + m][sl]; //J denotes newsvendor cost (holding cost and backorder/lostsales cost in a period)

		J = Jcalculate(S, period, s,sl, m, alpha, lambda,holdingCost, backorderCost, MM,poissonDistProbability,binomialDistProbability);

		for (int p = period - 1; p >= 0; p--) {

		Z[p] = periodMatrix(I, S,period, Z[p + 1],J, n, sl, s, m, alpha, lambda, unitCost, fixedCost, retailPrice,

		returnCredit, holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,binomialDistProbability);

		}
		
		double[] result0 = new double[3];
		for (int inv = M; inv <  m+M; inv++) {
			for (int j = 0; j <  m; j++) {
				
					result0=new double[]{inv,j,Z[0][inv][j]};
				Solution0.println(Arrays.toString(result0));
				
					}
				}
		double[] result1 = new double[3];
		for (int inv = M; inv <  m+M; inv++) {
			for (int j = 0; j <  m; j++) {
				
					result1=new double[]{inv,j,Z[1][inv][j]};
				Solution1.println(Arrays.toString(result1));
				
					}
				}
			
	
			double[] result2 = new double[3];
			for (int inv = M; inv <  m+M; inv++) {
				for (int j = 0; j <  m; j++) {
					result2=new double[]{inv,j,Z[2][inv][j]};
				Solution2.println(Arrays.toString(result2));
				
					}
				}
			
			
			double[] result3 = new double[3];
			for (int inv = M; inv <  m+M; inv++) {
				for (int j = 0; j <  m; j++) {
					
					result3=new double[]{inv,j,Z[3][inv][j]};
				Solution3.println(Arrays.toString(result3));
				
					}
				}
			
			double[] result4 = new double[3];
			for (int inv = M; inv < m+M; inv++) {
				for (int j = 0; j < m; j++) {
				
					result4=new double[]{inv,j,Z[4][inv][j]};
				Solution4.println(Arrays.toString(result4));
				
					}
				}
			
			
			double[] result5 = new double[3];
			for (int inv = M; inv <  m+M; inv++) {
				for (int j = 0; j <  m; j++) {
				
					result5=new double[]{inv,j,Z[5][inv][j]};
				Solution5.println(Arrays.toString(result5));
				
					}
				}
			
		
		//BigS - UptoOrder level (inventory level after ordering)
		int[][][] bigS = new int[period][n][sl]; 

		for (int p = period - 1; p >= 0; p--) {

		bigS[p] = periodMatrixIndex(I,S,period, Z[p + 1],J,n, sl, s, m, alpha, lambda, unitCost, fixedCost,

		retailPrice, returnCredit, holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,binomialDistProbability);

		}

double[][][] ZC = new double[period + 1][n][sl];
		

		ZC[period] = ZEndcalculate(I, period, n,sl, m, alphaC, unitCost, salvageValue, retailPrice, returnCredit,

		fixedCost,unitCost, binomialDistProbabilityC);
		
		double[][] JC = new double[s + m][sl];

		JC = Jcalculate(S, period, s,sl, m, alphaC, lambdaC,holdingCost, backorderCost, MM,poissonDistProbabilityC,binomialDistProbabilityC);

		for (int p = period - 1; p >= 0; p--) {

		ZC[p] = periodMatrix(I, S,period, ZC[p + 1],JC, n, sl, s, m, alphaC, lambdaC, unitCost, fixedCost, retailPrice,

		returnCredit, holdingCost, backorderCost, BigM, MM, M, poissonDistProbabilityC,binomialDistProbabilityC);

		}
		
		

		int[][][] bigSC = new int[period][n][sl];

		for (int p = period - 1; p >= 0; p--) {

		bigSC[p] = periodMatrixIndex(I, S,period, ZC[p + 1],JC, n, sl, s, m, alphaC, lambdaC, unitCost, fixedCost, retailPrice,

				returnCredit, holdingCost, backorderCost, BigM, MM, M, poissonDistProbabilityC,binomialDistProbabilityC);
		}

		double OP = Z[0][M][0];
		double OPC = ZC[0][M][0];
	
		double[] arr = new double[nresult];

		arr = SimulationMatrix(period, bigS,bigSC,  OP, OPC, m, alphaSimulation, lambdaSimulation, unitCost,
				fixedCost, retailPrice, returnCredit, holdingCost, backorderCost, salvageValue, M, nresult);

		
		Experiment[i] = arr;

		Solution.println(Arrays.toString(arr));
		
		
		
		System.out.println(Arrays.toString(arr));
		System.out.println(OP);
		System.out.println(OPC);
		

		}
		Solution.close();
		Solution0.close(); 
		Solution1.close(); 
		Solution2.close(); 
		Solution3.close();
		Solution4.close();
		Solution5.close();


		Instant inst2 = Instant.now();
		System.out.println("Elapsed Time: "+ Duration.between(inst1, inst2).toString());
		}

////////////////Simulation for n many random demand and return/////////////////////
	
		public static double[] SimulationMatrix(int period, int[][][] bigS,int[][][] bigSC,

		double OP, double OPC,int m, double alphaSimulation, double lambdaSimulation, int unitCost, int fixedCost,
		int retailPrice, int returnCredit, double holdingCost, double backorderCost, int salvageValue,int M, int nresult) {

		int n = 100000;

		double[] arr = new double[nresult];

		double[][] arr2 = new double[n][nresult];

		for (int i = 0; i < n; i++) {

		arr2[i] = ProcessMatrix(period, bigS, bigSC,OP, OPC, m, alphaSimulation, lambdaSimulation,

		unitCost, fixedCost, retailPrice,

		returnCredit, holdingCost, backorderCost, salvageValue, M, nresult);

		}

		arr = meanOfSecondInstance(arr2, nresult);
		return arr;

		}
////////////////Simulation for 1 random demand and return/////////////////////
		public static double[] ProcessMatrix(int period, int[][][] bigS, int[][][] bigSC, double OP, double OPC,

		int m, double alphaSimulation, double lambdaSimulation, int unitCost, int fixedCost, int retailPrice,

		int returnCredit, double holdingCost, double backorderCost, int salvageValue, int M, int nresult) {

		double[] result = new double[nresult];

		int[] I = new int[period + 1]; //inventory
		int[] R = new int[period + 1]; //return
		int[] O = new int[period]; //order
		int[] D = new int[period]; //demand
		int[] BO = new int[period]; //backorder
		int[] IO = new int[period + 1]; //inventory on hand
		int[] SL = new int[period + 1]; //sales
		int[] HI = new int[period]; // holding inventory

		int[] IC = new int[period + 1];
		int[] RC = new int[period + 1];
		int[] OC = new int[period];
		int[] BOC = new int[period];
		int[] IOC = new int[period + 1];
		int[] SLC = new int[period + 1];
		int[] HIC = new int[period];
		

		int[] numberOfOrder = new int[period];
		int[] numberOfOrderC = new int[period];

		for (int p = 0; p < period; p++) {
		D[p] = getPoissonRandom(lambdaSimulation);
		if (D[p]>m) {
		D[p]=m;
		}
		}

		I[0] = 0;
		R[0] = 0;


		O[0] = Math.max(0, (bigS[0][M][0]));


		BO[0] = Math.max(0, D[0] - I[0] - O[0] - R[0]);
		IO[0] = Math.max(0, I[0]);
		SL[0] = Math.min(IO[0] + O[0], D[0]);
		HI[0] = Math.max(0, I[0] + O[0] + R[0] - D[0]);

		IC[0] = 0;
		RC[0] = 0;


		OC[0] = Math.max(0, (bigSC[0][M][0]));


		BOC[0] = Math.max(0, D[0] - IC[0] - OC[0] - RC[0]);
		IOC[0] = Math.max(0, IC[0]);
		SLC[0] = Math.min(IOC[0] + OC[0], D[0]);
		HIC[0] = Math.max(0, IC[0] + OC[0] + RC[0] - D[0]);

		for (int i = 1; i < period; i++) {

		I[i] = I[i - 1] + O[i - 1] + R[i - 1] - D[i - 1];
		R[i] = getBinomialRandom(SL[i - 1], alphaSimulation);


		O[i] = Math.max(0, (bigS[i][I[i] + M][SL[i - 1]] - I[i]));


		BO[i] = Math.max(0, D[i] - I[i] - O[i] - R[i]);
		IO[i] = Math.max(0, I[i]);
		SL[i] = Math.min(IO[i] + O[i] + R[i], D[i] + BO[i - 1]);
		HI[i] = Math.max(0, I[i] + O[i] + R[i] - D[i]);

		IC[i] = IC[i - 1] + OC[i - 1] + RC[i - 1] - D[i - 1];
		RC[i] = getBinomialRandom(SLC[i - 1], alphaSimulation);

		OC[i] = Math.max(0, (bigSC[i][IC[i] + M][SLC[i - 1]] - IC[i]));

		BOC[i] = Math.max(0, D[i] - IC[i] - OC[i] - RC[i]);
		IOC[i] = Math.max(0, IC[i]);
		SLC[i] = Math.min(IOC[i] + OC[i] + RC[i], D[i] + BOC[i - 1]);
		HIC[i] = Math.max(0, IC[i] + OC[i] + RC[i] - D[i]);


		}



		I[period] = I[period - 1] + O[period - 1] + R[period - 1] - D[period - 1];
		R[period] = getBinomialRandom(SL[period - 1], alphaSimulation);

		IC[period] = IC[period - 1] + OC[period - 1] + RC[period - 1] - D[period - 1];
		RC[period] = getBinomialRandom(SLC[period - 1], alphaSimulation);



		double procurementCost = 0;
		double totalOrderingCost = 0;
		double totalHoldingCost = 0;
		double totalBackorderCost = 0;
		double returnCost = 0;
		double returnCost1 = 0;
		double totalEndingInventoryCost = 0;
		double Revenue = 0;
		double Revenue1 = 0;


		double procurementCostC = 0;
		double totalOrderingCostC = 0;
		double totalHoldingCostC = 0;
		double totalBackorderCostC = 0;
		double returnCostC = 0;
		double returnCost1C = 0;
		double totalEndingInventoryCostC = 0;
		double RevenueC = 0;
		double Revenue1C = 0;

		int sumNumberOfOrder =0;
		double CSL = 0;
		int sumNumberOfOrderC =0;
		double CSLC = 0;
		
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

		Revenue1 += SL[p] * retailPrice;
		Revenue = Revenue1 + SL[period] * retailPrice;

		

		
		procurementCostC += OC[p] * unitCost;

		if (OC[p] > 0) {
		numberOfOrderC[p] = 1;
		} else {
		numberOfOrderC[p] = 0;
		}

		totalOrderingCostC += numberOfOrderC[p] * fixedCost;
		
		sumNumberOfOrderC += numberOfOrderC[p] ;

		totalHoldingCostC += HIC[p] * holdingCost;

		totalBackorderCostC += BOC[p] * backorderCost;

		if (IC[period] + RC[period] < 0) {
		totalEndingInventoryCostC = (IC[period] + RC[period]) * -(unitCost)
		- IC[period] * alphaSimulation * (returnCredit - salvageValue);
		} else if (IC[period] < 0) {
		totalEndingInventoryCostC = (IC[period] + RC[period]) * -salvageValue
		- IC[period] * alphaSimulation * (returnCredit - salvageValue);
		} else {
		totalEndingInventoryCostC = (IC[period] + RC[period]) * -salvageValue;
		}

		if (IC[period] < 0) {
		SLC[period] = -(IC[period]);
		} else {
		SLC[period] = 0;
		}

		returnCost1C += RC[p] * returnCredit;
		returnCostC = returnCost1C + RC[period] * returnCredit;

		Revenue1C += SLC[p] * retailPrice;
		RevenueC = Revenue1C + SLC[period] * retailPrice;
		
	
		

		double salvage=I[period]+R[period];
		double SLV= salvage;

		double salvageC=IC[period]+RC[period];
		double SLVC= salvageC;
		}
		
		
		
	
		int[] orderPeriods = new int[period+1];
		int [] cycleBO= new int [period+1];
		int [] numberOfBO= new int [period];
		
		int[] orderPeriodsC = new int[period+1];
		int [] cycleBOC= new int [period+1];
		int [] numberOfBOC= new int [period];
		

///////////////////////////////////////////////////////
// CSL : Cycle Service Level for the smart retailer (Type 1)
		
		int j=0;
		int p=0;

		while (p<period) {
			
			if (numberOfOrder[p]==1) {
				
				orderPeriods[j]=p ;      
			
				j=j+1;
				
				}  
			
		p=p+1;
		
		orderPeriods[j]=period;
		
		} 

		

		int cycleService=0; 
		cycleBO[period] =0;

		for (int i = 0; i < sumNumberOfOrder ; i++) {

			for (int t = orderPeriods[i]; t < orderPeriods[i+1]; t++) {
				
				cycleBO[i] +=BO[t];
				
			}
			
		}		
			
		
		
		for (int i = 0; i < sumNumberOfOrder; i++) {
			
			
				if(cycleBO[i]==0) {
					
					cycleService=cycleService+1		;  
			
		}
		}

		CSL =(double)cycleService/sumNumberOfOrder ;
		
		
////////////////////////////////////////////////////////////
// CSLC : Cycle Service Level for the naive retailer (Type 1)
		
		int jj=0;
		int pp=0;

		while (pp<period) {
			
			if (numberOfOrderC[pp]==1) {
				
				orderPeriodsC[jj]=pp ;      
			
				jj=jj+1;
				
				}  
			
		pp=pp+1;
		
		orderPeriodsC[jj]=period;
		
		} 

		

		int cycleServiceC=0; 
		
		cycleBOC[period] =0;

		for (int i = 0; i < sumNumberOfOrderC ; i++) {

			for (int t = orderPeriodsC[i]; t < orderPeriodsC[i+1]; t++) {
				
				cycleBOC[i] +=BOC[t];
				
			}
			
		}		
			
		
		
		for (int i = 0; i < sumNumberOfOrderC; i++) {
			
			
				if(cycleBOC[i]==0) {
					
					cycleServiceC=cycleServiceC+1		;  
			
		}
		}

		CSLC =(double)cycleServiceC/sumNumberOfOrderC ;
		
///////////////////////////////////////////////////////////////////////////////////////
// fraverage : Fill rate for the smart retailer	(Type 2)	
		double fr =0;
		double fraverage =0;
		double frC =0;
		double fraverageC =0;
		double demand = 0;
		
		
		
		double[] inventoryonhand = new double[period];

		for (int ppp = 0; ppp < period ; ppp++) {
			
			demand +=D[ppp];

		inventoryonhand[ppp] = IO[ppp] + O[ppp] + R[ppp];


		fr += (Math.min(D[ppp], inventoryonhand[ppp])) ;


		fraverage =  fr /(double) demand;
		
///////////////////////////////////////////////////////////////////////
// fraverageC : Fill rate for the naive retailer	(Type 2)	

		double[] inventoryonhandC = new double[period];

		
		

		inventoryonhandC[ppp] = IOC[ppp] + OC[ppp] + RC[ppp];

		;

		
		frC += (Math.min(D[ppp], inventoryonhandC[ppp])) ;
	


		

		fraverageC = frC /demand;
		}
////////////////////////////////////////////////////////////////////////
		
		double salvage=I[period]+R[period];
		double SLV= salvage;

		double salvageC=IC[period]+RC[period];
		double SLVC= salvageC;
		
		double pc = procurementCost;
		double oc = totalOrderingCost;
		double hc = totalHoldingCost;
		double bc = totalBackorderCost;
		double tec = totalEndingInventoryCost;
		double rc = returnCost;
		double rv = Revenue;
		double OPR = OP;


		double pcC = procurementCostC;
		double ocC = totalOrderingCostC;
		double hcC = totalHoldingCostC;
		double bcC = totalBackorderCostC;
		double tecC = totalEndingInventoryCostC;
		double rcC = returnCostC;
		double rvC = RevenueC;
		double OPRC = OPC;

		double totalCost = procurementCost + totalOrderingCost + totalHoldingCost + totalBackorderCost
				+ totalEndingInventoryCost + returnCost;

				double Profit = Revenue - totalCost;
				
		double totalCostC = procurementCostC + totalOrderingCostC + totalHoldingCostC + totalBackorderCostC
		+ totalEndingInventoryCostC + returnCostC;

		double ProfitC = RevenueC - totalCostC;


		double tc = totalCost;
		double pr = Profit;
		double tcC = totalCostC;
		double prC = ProfitC;
		//pr profit , OPR optimal profit
		result = new double[] {pr, -OPR ,prC, -OPRC ,CSL,fraverage,CSLC,fraverageC, pc, oc, hc, bc,tec, rc, tc, rv,SLV,pcC, ocC, 
				hcC, bcC, tecC, rcC, tcC, rvC, SLVC, O[0], O[1],O[2],O[3],BO[0],BO[1],BO[2],BO[3],HI[0],HI[1],HI[2],HI[3],
				R[1],R[2],R[3],R[4],SL[0],SL[1],SL[2],SL[3],SL[4],I[4], OC[0], OC[1],OC[2],OC[3],BOC[0],BOC[1],BOC[2],BOC[3],HIC[0],
				HIC[1],HIC[2],HIC[3],RC[1],RC[2],RC[3],RC[4],SLC[0],SLC[1],SLC[2],SLC[3],SLC[4],IC[4]}; 
		
				
		return result;

		}
		
//////////////////Profit for each upto order level for each period////////////////////////////////	
		
		private static double[][][] periodMatrix1(int [] I, int [] S,int period, double[][] ZN,double[][] J,int n, int sl,int s, int m, double alpha, double lambda, 

				int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost, int BigM, int MM, int M,
				double[] poissonDistProbability,double[][] binomialDistProbability) {


				double[][][][][] ZZZZ = new double[n][sl][s][m+1][sl];
				double[][][][] ZZZ = new double[n][sl][s][m+1];
				double[][][] ZZ = new double[n][sl][s];

					
				
				for (int inv = 0; inv < n; inv++) { //Inventory

				for (int j = 0; j < sl; j++) { //previous sales

				for (int i = 0; i < s; i++) { //S up to order

				for (int jj = 0; jj <= m; jj++) { //demand

				for (int k = 0; k <= j; k++) { //return

				if (I[inv] < S[i]) {

				ZZZZ[inv][j][i][jj][k] = (binomialDistProbability[j][k]) * (poissonDistProbability[jj])

				* ((S[i]- I[inv]) * unitCost + fixedCost + J[i][j] + k * returnCredit

						+ ZN[S[i] + k - jj + M][Math.min(
								(Math.max(0, I[inv]) + Math.max(0, (S[i] - I[inv])) + k),
								(jj - (Math.min(0, I[inv]))))]
								- (Math.min((Math.max(0, I[inv]) + Math.max(0, (S[i] - I[inv])) + k),
								(jj - (Math.min(0, I[inv]))))) * retailPrice);

				} else if (I[inv] ==S[i]) {

				ZZZZ[inv][j][i][jj][k] =(binomialDistProbability[j][k]) * (poissonDistProbability[jj])

				* ( J[i][j] + k * returnCredit

						+ ZN[S[i] + k - jj + M][Math.min(
								(Math.max(0, I[inv]) + Math.max(0, (S[i] - I[inv])) + k),
								(jj - (Math.min(0, I[inv]))))]
								- (Math.min((Math.max(0, I[inv]) + Math.max(0, (S[i] - I[inv])) + k),
								(jj - (Math.min(0, I[inv]))))) * retailPrice);

				} else {

				ZZZZ[inv][j][i][jj][k] = BigM;

				}

				}

				}

				}

				}

				}

				ZZZ = sumOfFifthInstance(ZZZZ, sl, s, m);

				ZZ = sumOfFourthInstance(ZZZ, sl, s, m);

				return ZZ;

				}
//////////////////Optimal profit for each period////////////////////////////////	

				private static double[][] periodMatrix(int [] I , int [] S,int period, double[][] ZN, double[][] J,int n, int sl, int s, int m, double alpha,

				double lambda, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost, double backorderCost,int BigM, int MM, int M,
				double[] poissonDistProbability,double[][] binomialDistProbability) {

				double[][][] ZZ = new double[n][sl][s];

				double[][] Z = new double[n][sl];

				ZZ = periodMatrix1(I,S, period, ZN,J, n,sl,s, m, alpha, lambda, unitCost, fixedCost, retailPrice, returnCredit,

				holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,binomialDistProbability);

				Z = minOfThirdInstance(ZZ, sl, s);

				return Z;

				}
//////////////////Optimal upto order levels for each period////////////////////////////////	
				private static int[][] periodMatrixIndex(int [] I, int [] S, int period, double[][] ZN, double[][] J,int n, int sl, int s, int m,

				double alpha,

				double lambda, int unitCost, int fixedCost, int retailPrice, int returnCredit, double holdingCost,

				double backorderCost, int BigM, int MM, int M, double[] poissonDistProbability,double[][] binomialDistProbability) {

				int[][] SS = new int[n][sl];

				double[][][] ZZ = new double[n][sl][s];

				int[][] bigS = new int[n][sl];

				ZZ = periodMatrix1(I,S, period, ZN,J, n,sl,s, m, alpha, lambda, unitCost, fixedCost, retailPrice, returnCredit,

				holdingCost, backorderCost, BigM, MM, M, poissonDistProbability,binomialDistProbability);

				SS = minIndexOfThirdInstance(ZZ, sl, s);

				for (int inv = 0; inv < n; inv++) {

				for (int j = 0; j < sl; j++) {

				bigS[inv][j] = SS[inv][j] - MM;

				}

				}

				return bigS;

				}
				
//////////////////Newsvendor cost////////////////////////////////				
				private static double[][] Jcalculate(int [] S, int period, int s, int sl, int m, double alpha, double lambda,double holdingCost,

				double backorderCost, int MM, double[] poissonDistProbability,double[][] binomialDistProbability) {

				
	
				double[][][][] JJJ = new double[s][sl][m + 1][sl];

				double[][][] JJ = new double[s][sl][m + 1];

				double[][] J = new double[s][sl];

				for (int i = 0; i < s; i++) {

				for (int j = 0; j < sl; j++) {

				for (int jj = 0; jj <= m; jj++) {

				for (int k = 0; k <= j; k++) {

				if (S[i]  + k - jj > 0) {

				JJJ[i][j][jj][k] = (binomialDistProbability[j][k]) *poissonDistProbability[jj] * (S[i] + k - jj)

				* holdingCost;

				} else {

				JJJ[i][j][jj][k] = (binomialDistProbability[j][k]) * poissonDistProbability[jj] * (jj - S[i]  - k)

						* backorderCost;
				

				}

				}

				}

				}

				}

				JJ = sumOfFourthInstanceJ(JJJ, sl, m);

				J = sumOfThirdInstanceJ(JJ, sl, m);

				return J;

				}
//////////////////Optimal profit for terminal period////////////////////////////////	
				private static double[][] ZEndcalculate(int [] I, int period, int n, int sl, int m, double alpha, int backorderCost, int salvageValue,

				int retailPrice, int returnCredit, int fixedCost, int unitCost,double[][] binomialDistProbability) {

				
					double[][][] ZZEnd = new double[n][sl][sl];
					double[][] ZEnd = new double[n][sl];


				for (int inv = 0; inv < n; inv++) {

				for (int j = 0; j < sl; j++) {

				for (int k = 0; k <= j; k++) {

				
				if ((I[inv] + k) < 0) {
					ZZEnd[inv][j][k] = (binomialDistProbability[j][k]) * (((I[inv] + k) * -unitCost)  + k * returnCredit
					- I[inv] * alpha * (returnCredit - salvageValue) + retailPrice * I[inv]);
					} else if (I[inv] < 0) {
					ZZEnd[inv][j][k] = (binomialDistProbability[j][k]) * (((I[inv] + k) * -salvageValue) + k * returnCredit
					- I[inv] * alpha * (returnCredit - salvageValue) + retailPrice * I[inv]);
					} else {
					ZZEnd[inv][j][k] = (binomialDistProbability[j][k]) * (((I[inv] + k) * -salvageValue) + k * returnCredit);
					}
				}

				}

				}

				ZEnd = sumOfThirdInstance(ZZEnd, sl);

				return ZEnd;

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

		static int nCr(int n, int k) {

		if (k > n / 2)

		k = n - k;

		int answer = 1;

		for (int i = 1; i <= k; i++) {

		answer *= (n - k + i);

		answer /= i;

		}

		return answer;

		}

		static double binoPdf(int n, int k, double alpha) {

		return nCr(n, k) * (double) Math.pow(alpha, k) * (double) Math.pow(1 - alpha, n - k);

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

		public static double[][] sumOfThirdInstance(double[][][] array, int sl) {

		double[][] result = new double[array.length][sl];

		double total;

		for (int inv = 0; inv < array.length; inv++) {

		for (int j = 0; j < sl; j++) {

		total = 0;

		for (int k = 0; k <= j; k++) {

		total += array[inv][j][k];

		}

		result[inv][j] = total;

		}

		}

		return result;

		}

		public static double[][][] sumOfFourthInstanceJ(double[][][][] array, int sl, int m) {

		double[][][] result = new double[array.length][sl][m + 1];

		double total;

		for (int i = 0; i < array.length; i++) {

		for (int j = 0; j < sl; j++) {

		for (int jj = 0; jj <= m; jj++) {

		total = 0;

		for (int k = 0; k <= j; k++) {

		total += array[i][j][jj][k];

		}

		result[i][j][jj] = total;

		}

		}

		}

		return result;

		}

		public static double[][][] sumOfFourthInstance(double[][][][] array, int sl, int s, int m) {

		double[][][] result = new double[array.length][sl][s];

		double total;

		for (int inv = 0; inv < array.length; inv++) {

		for (int j = 0; j < sl; j++) {

		for (int i = 0; i < s; i++) {

		total = 0;

		for (int jj = 0; jj <= m; jj++) {

		total += array[inv][j][i][jj];

		}

		result[inv][j][i] = total;

		}

		}

		}

		return result;

		}

		public static double[][][][] sumOfFifthInstance(double[][][][][] array, int sl, int s, int m) {

		double[][][][] result = new double[array.length][sl][s][m + 1];

		double total;

		for (int inv = 0; inv < array.length; inv++) {

		for (int j = 0; j < sl; j++) {

		for (int i = 0; i < s; i++) {

		for (int jj = 0; jj <= m; jj++) {

		total = 0;

		for (int k = 0; k <= j; k++) {

		total += array[inv][j][i][jj][k];

		}

		result[inv][j][i][jj] = total;

		}

		}

		}

		}

		return result;

		}

		public static double[][] meanOfThirdInstance(double[][][] array, int array2size, int array3size) {

		double[][] result = new double[array.length][array2size];

		double mean;

		for (int i = 0; i < array.length; i++) {

		for (int j = 0; j < array2size; j++) {

		mean = 0;

		for (int jj = 0; jj <= array3size; jj++) {

		mean += array[i][j][jj] / (double) array3size;

		}

		result[i][j] = mean;

		}

		}

		return result;

		}

		public static double[][] sumOfThirdInstanceJ(double[][][] array, int sl, int m) {

		double[][] result = new double[array.length][sl];

		double sum;

		for (int i = 0; i < array.length; i++) {

		for (int j = 0; j < sl; j++) {

		sum = 0;

		for (int jj = 0; jj <= m; jj++) {

		sum += array[i][j][jj];

		}

		result[i][j] = sum;

		}

		}

		return result;

		}

		public static double[][][] meanOfFourthInstance(double[][][][] array, int array2size, int array3size,

		int array4size) {

		double[][][] result = new double[array.length][array2size][array3size];

		double mean;

		for (int inv = 0; inv < array.length; inv++) {

		for (int j = 0; j < array2size; j++) {

		for (int i = 0; i < array3size; i++) {

		mean = 0;

		for (int jj = 0; jj <= array4size; jj++) {

		mean += array[inv][j][i][jj] / (double) array4size;

		}

		result[inv][j][i] = mean;

		}

		}

		}

		return result;

		}

		public static double[][] minOfThirdInstance(double[][][] array, int array2size, int array3size) {

		double[][] result = new double[array.length][array2size];

		for (int inv = 0; inv < array.length; inv++) {

		for (int j = 0; j < array2size; j++) {

		double min = 100000;

		for (int i = 0; i < array3size; i++) {

		if (array[inv][j][i] < min) {

		min = array[inv][j][i];

		}

		}

		result[inv][j] = min;

		}

		}

		return result;

		}

		public static int[][] minIndexOfThirdInstance(double[][][] array, int array2size, int array3size) {

		int[][] result = new int[array.length][array2size];

		for (int inv = 0; inv < array.length; inv++) {

		for (int j = 0; j < array2size; j++) {

		double min = 10000;

		for (int i = 0; i < array3size; i++) {

		if (array[inv][j][i] < min) {

		min = array[inv][j][i];

		result[inv][j] = i;

		}

		}

		}

		}

		return result;

		}

		public static int[] minIndexOfSecondInstance(int[][] array, int array2size) {

		int[] result = new int[array.length];

		for (int j = 0; j < array2size; j++) {

		double min = 100;

		for (int inv = 0; inv < array2size && array[j][inv] != 0; inv++) {

		if (array[j][inv] < min) {

		min = array[j][inv];

		result[j] = inv;

		}

		}

		}

		return result;

		}

		

		

		public static double Mean(double[] array) {

		double mean = 0;

		for (int i = 0; i < array.length; i++) {

		mean += array[i] / (double) array.length;

		}

		return mean;

		}

		public static double[] meanOfSecondInstance(double[][] array, int array2size) {

		double[] result = new double[array2size];

		double mean;

		for (int i = 0; i < array2size; i++) {

		mean = 0;

		for (int j = 0; j < array.length; j++) {

		mean += array[j][i] / (double) array.length;

		}

		result[i] = mean;

		}

		return result;

		}

		
		}
