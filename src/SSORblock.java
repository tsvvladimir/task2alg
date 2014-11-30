/**
 * Created by natalia on 09.11.14.
 */
public class SSORblock {

    Matrix A;
    Matrix b;
    Matrix x;   //start approx

    int sizeofBlock;

    double w;

    public SSORblock(Matrix A, Matrix b, Matrix x, double w, int sizeofBlock) {
        this.A = A;
        this.b = b;
        this.x = x;
        this.w = w;
        this.sizeofBlock = sizeofBlock;
    }

    //if order=1 - count forward x_k/2, if order=-1 count back x_k
    private Matrix countBackOrForward(Matrix prevX, int order) {
        prevX.show();
        int numOfBlocks = A.data.length / sizeofBlock;
        Matrix xk = new Matrix(prevX);

        int from, to;
        if(order > 0) {
            from = 0; to = numOfBlocks;
        } else {
            from = numOfBlocks-1; to = -1;
        }

        //every iteration counts block i of xk
        for (int i = from; i != to; i += order) {
            //get block i of x from previos step (k-1)
            Matrix x_kmin1_i = prevX.getBlock(sizeofBlock*i, sizeofBlock*(i+1), 0, 1);
            //System.out.print("\nx^k-1:\n"); x_kmin1_i.show();
            //get cur block i of b
            Matrix b_i = b.getBlock(sizeofBlock*i, sizeofBlock*(i+1), 0, 1);
            //System.out.print("\nbi:\n"); b_i.show();
            //get block ii of A
            Matrix a_ii = A.getBlock(sizeofBlock*i, sizeofBlock*(i+1), sizeofBlock*i, sizeofBlock*(i+1));
            //System.out.print("\naii:\n"); a_ii.show();

            //count sum on all j < i of a_ij * x_j^(k)
            Matrix sum1 = Matrix.Null(sizeofBlock, 1);
            for (int j = from; j != i; j += order) {
                Matrix a_ij = A.getBlock(i*sizeofBlock, (i+1)*sizeofBlock, sizeofBlock*j, sizeofBlock*(j+1));
                //System.out.print("\na ij:\n"); a_ij.show();
                Matrix x_j = prevX.getBlock(sizeofBlock*j, sizeofBlock*(j+1), 0, 1);    //!!!in lectures - xk
                //System.out.print("\nx j:\n"); x_j.show();
                //a_ij.times(x_j).show();
                sum1 = sum1.minus(a_ij.times(x_j));
                //System.out.print("\nsum1:\n"); sum1.show();

            }
            //System.out.print("\nsum1:\n"); sum1.show();

            //count sum on all j > i of a_ij * x_j^(k-1)
            Matrix sum2 = Matrix.Null(sizeofBlock, 1);
            for (int j = i + order; j != to; j += order) {
                Matrix a_ij = A.getBlock(i*sizeofBlock, (i+1)*sizeofBlock, sizeofBlock*j, sizeofBlock*(j+1));
                Matrix x_j = xk.getBlock(sizeofBlock*j, sizeofBlock*(j+1), 0, 1);   //!!!in lectures - x_kmin1
                sum2 = sum2.minus(a_ij.times(x_j));
            }
            //System.out.print("\nsum2:\n"); sum2.show();

            //count block i of xk
            Matrix firstElem = x_kmin1_i.mulOnConstant(1 - w);
            Matrix secondElem = sum1.plus(sum2).plus(b_i);

            Matrix xk_i = firstElem.plus(a_ii.degMin1().times(secondElem).mulOnConstant(w));
            //System.out.print("\nx^k i:\n"); xk_i.show();

            //copy block xk_i to xk
            for (int z = 0; z < sizeofBlock; z++)
                xk.setElement(i*sizeofBlock+z, 0, xk_i.GetElement(z, 0));
            //System.out.print("\nxk:\n"); xk.show();
            //System.out.print("\n");
        }
        xk.show();
        return xk;
    }

    public int count(double accuracy) {
        if (A == null || A.data.length % sizeofBlock != 0) {
            System.out.println("error: trying to run SSOR to null matrix or illigal sizof block");
            return 0;
        }

        int numOfSteps = 0;
        boolean goOn = true;

        Matrix prevX = new Matrix(x);   //vector x got on step k-1

        while (goOn) {
            numOfSteps++;

            //count forward
            Matrix xHalfK = countBackOrForward(prevX, 1);
            //now we have x k/2, so need to count xk
            //count back
            Matrix xk = countBackOrForward(xHalfK, -1);
            //Matrix xk = xHalfK;

            //!!!must be first if, but it doesn't work
            //if (xk.minus(prevX).norm2() < accuracy.norm2()) {//!!!may be it's better to count |xk| and |precX| first
            //if (numOfSteps > 5){//


            if (b.minus(A.times(xk)).norm2() < accuracy) {//accuracy.norm2()) {
                goOn = false;
                System.out.println("check: b-Ax=");
                b.show();
                System.out.println();
                Matrix must_be_B = A.times(xk);
                must_be_B.show();
                //b.minus(A.times(xk)).show();
            }
            prevX = new Matrix(xk);

            //b.minus(A.times(xk)).show();
            System.out.println();
            //if algorithm doesn't work we finish manually - in case of infinity
            if (numOfSteps > 100) {

                System.out.println("check: b-Ax=");
                b.show();
                System.out.println();
                Matrix must_be_B = A.times(xk);
                must_be_B.show();
                //b.minus(A.times(xk)).show();
                //xk.show();
                System.out.println("check: number of steps is more then 1000, it's not good;)");
                goOn = false;
            }
        }

        return numOfSteps;
    }

    public void printResultChangingW(Matrix accuracy) {
         /*if (accuracy == null) {
             accuracy = new Matrix(x);
             for (int i = 0; i < x.data[0].length; i++)
                 accuracy.setElement(i, 0, 0.00000000000000000001);
         }*/
        double accur = 0.00001;

        for (double w = 1.1; w < 1.2; w += 0.1) {
            this.w = w;

            System.out.print("steps: (" + w + ","+ count(accur)+");\n ");

        }
    }

    public static void main(String []args) {
        System.out.println("try it");

        int n = 3;

        Matrix A = new Matrix(new double[][]{{1,2, 3}, {2,1,4}, {3,4,1}});//
        Matrix b = new Matrix(new double[][]{{1}, {2}, {3}});//



        Matrix help1 = new Matrix(n, 1);
        for(int i = 0; i < n; i++)
            help1.data[i][0] = 1;
        //Matrix b = A.times(help1);//
        b = Matrix.random(n,1);
        A = Matrix.randomsym(n, n);

        //Matrix approxim = new Matrix(new double[][]{{10}, {10}, {10}});//
        Matrix approxim =  Matrix.random(n, 1);


        SSORblock sb = new SSORblock(A, b, approxim, 1, 1);
        sb.printResultChangingW(null);

    }
}
