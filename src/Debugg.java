import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Random;

public class Debugg {
    public static void main(String[] args) throws FileNotFoundException {
        int period=5;
        int lambdaSimulation=14;
        int[] D= new int[period];
        int sum=0;
        double avg=0;


        for (int p = 0; p < period; p++) {
            D[p] = getPoissonRandomweb(lambdaSimulation);
            //if (D[p]>m) {
            //D[p]=m;
            //}
            sum+=D[p];
        }
        avg=sum/period;
        System.out.println(Arrays.toString(D));
        System.out.println(avg);

    }
    private static int getPoissonRandomweb(double mean) {
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
    public static double[][] binomialDist(int sl, double alpha) {

        double[][] P = new double[sl][sl];



        for (int j = 0; j <sl ; j++) {
            for (int k = 0; k <=j ; k++) {
                P[j][k]=nCr(j, k) * (double) Math.pow(alpha, k) * (double) Math.pow(1 - alpha, j - k);
            }
        }


        return P;

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
    public static int getBinomialRandom(int n, double alpha) {

        int x = 0;

        for (int i = 0; i < n; i++) {

            if (Math.random() < alpha)

                x++;

        }

        return x;

    }


}
