/**
 * Created by vladimirtsvetkov on 21/10/14.
 */
public class CG {
    Matrix A;
    Matrix B;
    Matrix D;
    Matrix L;
    Matrix U;
    public CG(Matrix A, Matrix B) {
        this.A = A;
        this.B = B;
    }
    public  CG(Matrix D, Matrix L, Matrix U, Matrix B) {
        this.D = D;
        this.L = L;
        this.U = U;
        this.B = B;
    }
    public void dsplit() {
        Matrix D = new Matrix(A.M, A.M);
        Matrix L = new Matrix(A.M, A.M);
        Matrix U = new Matrix(A.M, A.M);
        for(int i = 0; i < A.M; i++) {
            D.setElement(i, i, A.GetElement(i, i));
        }
        for(int i = 0; i < (A.M - 1); i++)
            for(int j = i + 1; j < A.N; j++)
                U.setElement(i, j, A.GetElement(i, j));
        for(int i = 0; i < (A.M - 1); i++)
            for(int j = i + 1; j < A.N; j++)
                L.setElement(j, i, A.GetElement(j, i));
        this.D = D;
        this.L = L.muldig(-1);
        this.U = U.muldig(-1);
        System.out.println("A");
        this.A.show();
        System.out.println("D");
        this.D.show();
        System.out.println("L");
        this.L.show();
        System.out.println("U");
        this.U.show();
        System.out.println("Solution:");
    }
    public Matrix execute1(int stepmax) {
        Matrix xk = new Matrix(A.M, 1);
        Matrix xk1 = new Matrix(A.M, 1);
        for(int i = 0; i < A.M; i++) {
            xk1.setElement(i, 0, 1.0);
        }
        Matrix dk = B.minus(A.times(xk));
        Matrix rk = dk;

        for(int step = 0; step < stepmax; step++) {

            double alphaknom = rk.transpose().times(rk).GetElement(0, 0);
            double alphakdenom = dk.transpose().times(A).times(dk).GetElement(0, 0);
            double alphak = alphaknom/alphakdenom;

            xk1 = xk.plus(dk.muldig(alphak));
            xk = xk1;
            Matrix rkold = rk;
            rk = rk.minus(A.times(dk).muldig(alphak));
            double bknom = rk.transpose().times(rk).GetElement(0, 0);
            double bkdenom = rkold.transpose().times(rkold).GetElement(0, 0);
            double bk = bknom/bkdenom;
            dk = rk.plus(dk.muldig(bk));
            xk.show();
            System.out.println("\n");
            Matrix nevyz = (B.minus(A.times(xk)));
            double norm = nevyz.norm2();
            if (norm < 0.0001) {
                System.out.println("Congrats!!!!!! steps CG:" + step + "\n");
                break;
            }
        }
        return xk1;
    }

    public Matrix execute2(int stepmax, double w) {
        Matrix xk = new Matrix(D.M, 1);
        Matrix xk1 = new Matrix(D.M, 1);
        for(int i = 0; i < D.M; i++) {
            xk1.setElement(i, 0, 1.0);
        }
        for(int step = 0; step < stepmax; step++) {
            Matrix help1 = L.plus(D.muldig(1/w));
            Matrix help2 = help1.degMin1();
            Matrix help3 = D.muldig(1/w).minus(D).minus(U);
            xk1 = help2.times(help3.times(xk).plus(B));
            xk = xk1;
            xk.show();
            System.out.println("\n");
        }
        return xk1;
    }
}
