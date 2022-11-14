
class HelloWorld {
    public static void main(String[] args) {
        int n =30;
        int r =15;
        if (r > n / 2)

            r = n - r;

        double answer = 1;

        for (int i = 1; i <= r; i++) {

            answer *= (n - r + i);

            answer = answer/i;

        }
        if(answer<0){System.out.println(n);
            System.out.println("  ");
            System.out.println(r);
            System.out.println("  ");
            System.out.println(answer);

        }


    }
}