public class Main {
    public static void testJacobi() {
        System.out.println("hello Jacobi!");
        //test Jacobi
        int n = 5;
        int maxsteps = 25;
        //Matrix A = Matrix.randomsym(n, n);
        Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {5.0, 7.0}} );
        Matrix help1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++)
            help1.data[i][0] = 1;
        Matrix B = new Matrix(A.times(help1));
        //Matrix B = new Matrix(new double[][] {{11}, {13}});
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
        double w = 1.3;
        //Matrix A = Matrix.randomsym(n, n);
        Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {1.0, 7.0}} );
        Matrix help1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++)
            help1.data[i][0] = 1;
        Matrix B = new Matrix(A.times(help1));
        System.out.println("matrix B:");
        B.show();
        //Matrix B = new Matrix(new double[][] {{11}, {13}});
        SOR sor = new SOR(A, B);
        sor.dsplit();
        Matrix result = sor.execute1(maxsteps, w);
        //Matrix result = sor.execute2(maxsteps, w);
        result.show();
        System.out.println("bye bye SOR!");
    }
    public static void  testoptimumiterJacobi() {
        //test optimum Jacobi
        System.out.println("hello optimJacobi!");
        int n = 5;
        int maxsteps = 25;
        double w = 1.3;
        //Matrix A = Matrix.randomsym(n, n);
        double m = -0.3;
        double M = 1.3;
        Matrix A = new Matrix(new double[][] {{0.1, 1, 0, 0, 0}, {0, 0.1, 0, 0, 0}, {0, 0, 0.5, 1, 0}, {0, 0, 0, 0.5, 1}, {0, 0, 0, 0, 0.5}} );
        A = A.plus(A.transpose()).muldig(0.5);
        Matrix help1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++)
            help1.data[i][0] = 1;
        Matrix B = new Matrix(A.times(help1));
        System.out.println("matrix B:");
        B.show();
        //Matrix B = new Matrix(new double[][] {{11}, {13}});
        optimJacobi opJ = new optimJacobi(A, B, m, M);
        opJ.dsplit();
        //Matrix result = opJ.execute1(maxsteps, w);
        Matrix result = opJ.execute(maxsteps);
        result.show();
        System.out.println("bye bye optimJacobi!");
    }
    public static void  testcheba() {
        //test optimum Jacobi
        System.out.println("hello Chebyshev!");
        int n = 5;
        int maxsteps = 25;
        double w = 1.3;
        //Matrix A = Matrix.randomsym(n, n);
        double m = 0.1;
        double M = 0.5;
        Matrix A = new Matrix(new double[][] {{0.1, 1, 0, 0, 0}, {0, 0.1, 0, 0, 0}, {0, 0, 0.5, 1, 0}, {0, 0, 0, 0.5, 1}, {0, 0, 0, 0, 0.5}} );
        Matrix help1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++)
            help1.data[i][0] = 1;
        Matrix B = new Matrix(A.times(help1));
        System.out.println("matrix B:");
        B.show();
        //Matrix B = new Matrix(new double[][] {{11}, {13}});
        Cheb ch = new Cheb(A, B, m, M);
        ch.dsplit();
        //Matrix result = opJ.execute1(maxsteps, w);
        Matrix result = ch.execute(maxsteps);
        result.show();
        System.out.println("bye bye Chebyshev!");
    }

    public static void testCG() {
        //test SOR
        System.out.println("hello CG!");
        int n = 5;
        int maxsteps = 10;
        double w = 1.3;
        //Matrix A = Matrix.randomsym(n, n);
        Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {1.0, 7.0}} );
        Matrix help1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++)
            help1.data[i][0] = 1;
        Matrix B = new Matrix(A.times(help1));
        System.out.println("matrix B:");
        B.show();
        //Matrix B = new Matrix(new double[][] {{11}, {13}});
        CG cg = new CG(A, B);
        //cg.dsplit();
        Matrix result = cg.execute1(maxsteps);
        //Matrix result = sor.execute2(maxsteps, w);
        result.show();
        System.out.println("bye bye CG!");
    }
    public static void main(String[] args) {
	// write your code here
     //testJacobi();
     //testSOR();
       testoptimumiterJacobi();
        //testcheba();
        //testCG();
    }
}
