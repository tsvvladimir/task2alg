public class Main {
    public static void testJacobi() {
        System.out.println("hello Jacobi!");
        //test Jacobi
        int n = 5;
        int maxsteps = 25;
        //Matrix A = Matrix.randomsym(n, n);
        Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {5.0, 7.0}} );
        Matrix help1 = new Matrix(n, 1);
        for(int i = 0; i < n; i++)
            help1.data[i][0] = 1;
        //Matrix B = new Matrix(A.times(help1));
        Matrix B = new Matrix(new double[][] {{11}, {13}});
        Jacobi j = new Jacobi(A, B);
        j.dsplit();
        Matrix result = j.execute(maxsteps);
        result.show();
        System.out.println("bye bye Jacobi!");

    }

    public static void testSOR() {
        //test SOR
        System.out.println("hello SOR!");
        int n = 5;
        int maxsteps = 25;
        double w = 1.9;
        //Matrix A = Matrix.randomsym(n, n);
        Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {5.0, 7.0}} );
        Matrix help1 = new Matrix(n, 1);
        for(int i = 0; i < n; i++)
            help1.data[i][0] = 1;
        //B = new Matrix(A.times(help1));
        Matrix B = new Matrix(new double[][] {{11}, {13}});
        SOR sor = new SOR(A, B);
        sor.dsplit();
        //Matrix result = sor.execute1(maxsteps, w);
        Matrix result = sor.execute2(maxsteps, w);
        result.show();
        System.out.println("bye bye SOR!");
    }
    public static void main(String[] args) {
	// write your code here
     //testJacobi();
     testSOR();
    }
}
