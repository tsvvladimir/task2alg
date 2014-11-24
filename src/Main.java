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
        double accur = 0.001;
        Matrix result = j.execute(maxsteps, accur);
        result.show();
        System.out.println("bye bye Jacobi!");

    }

    public static void testSOR() {
        //test SOR
        System.out.println("hello SOR!");
        int n = 5;
        int maxsteps = 25;
        double w = 1.0;
        //Matrix A = Matrix.randomsym(n, n);
        Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {5.0, 7.0}} );
        Matrix help1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++)
            help1.data[i][0] = 1;
        Matrix B = new Matrix(A.times(help1));
        System.out.println("matrix B:");
        B.show();
        //Matrix B = new Matrix(new double[][] {{11}, {13}});
        SOR sor = new SOR(A, B);
        sor.dsplit();
        double accur = 0.001;
        Matrix result = sor.execute1(maxsteps, w, accur);
        //Matrix result = sor.execute2(maxsteps, w);
        result.show();
        System.out.println("bye bye SOR!");
    }
    public static void  testoptimumiterJacobi() {
        //test optimum Jacobi
        System.out.println("hello optimJacobi!");
        //int n = 5;
        int maxsteps = 100;
        //double w = 1.3;
        //Matrix A = Matrix.randomsym(n, n);
        //double m = 0.1;
        //double M = 0.5;
        double m = -0.02;
        double M = 0.7;
        //Matrix A = new Matrix(new double[][] {{0.1, 1, 0, 0, 0}, {0, 0.1, 0, 0, 0}, {0, 0, 0.5, 1, 0}, {0, 0, 0, 0.5, 1}, {0, 0, 0, 0, 0.5}} );
        //Matrix A = new Matrix(new double[][] {{0.01,0.1,0.1}, {0.1,0.7,0}, {0.1,0,0.5}});
        Matrix A = new Matrix(new double[][] {{100.0,0.1,0.1}, {0.1,10.0/7.0,0}, {0.1,0,2.0}});
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
        double m = -0.02;
        double M = 0.7;
        //Matrix A = new Matrix(new double[][] {{0.1, 1, 0, 0, 0}, {0, 0.1, 0, 0, 0}, {0, 0, 0.5, 1, 0}, {0, 0, 0, 0.5, 1}, {0, 0, 0, 0, 0.5}} );
        Matrix A = new Matrix(new double[][] {{100.0,0.1,0.1}, {0.1,10.0/7.0,0}, {0.1,0,2.0}});
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
        //Matrix A = new Matrix(new double[][] {{2.0, 1.0}, {1.0, 7.0}} );
        Matrix A = new Matrix(new double[][] {{100.0,0.1,0.1}, {0.1,10.0/7.0,0}, {0.1,0,2.0}});
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

    public static void testSSORblock() {
        //SSOR block
        System.out.println("try it");

        int n = 9;

        //Matrix A = new Matrix(new double[][]{{1,0, 0}, {0,1,0}, {0,0,1}});//
        //Matrix b = new Matrix(new double[][]{{1}, {2}, {3}});//

        Matrix A = Matrix.randomsym(n, n);
        Matrix b = Matrix.random(n,1);

        //Matrix approxim = new Matrix(new double[][]{{10}, {10}, {10}});//
        Matrix approxim =  Matrix.random(n, 1);


        SSORblock sb = new SSORblock(A, b, approxim, 1, 3);
        sb.printResultChangingW(null);
    }

    public static void main(String[] args) {

     //testSSORblock();

     //testJacobi();
     testSOR();

     //testoptimumiterJacobi();
     //testcheba();
     //testCG();

    }
}
