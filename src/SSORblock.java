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
            //get cur block i of b
            Matrix b_i = b.getBlock(sizeofBlock*i, sizeofBlock*(i+1), 0, 1);
            //get block ii of A
            Matrix a_ii = A.getBlock(sizeofBlock*i, sizeofBlock*(i+1), sizeofBlock*i, sizeofBlock*(i+1));

            //count sum on all j < i of a_ij * x_j^(k)
            Matrix sum1 = Matrix.Null(sizeofBlock, 1);
            for (int j = from; j != i; j += order) {
                Matrix a_ij = A.getBlock(i*sizeofBlock, (i+1)*sizeofBlock, sizeofBlock*j, sizeofBlock*(j+1));
                Matrix x_j = prevX.getBlock(sizeofBlock*j, sizeofBlock*(j+1), 0, 1);    //!!!in lectures - xk
                sum1.minus(a_ij.times(x_j));
            }

            //count sum on all j > i of a_ij * x_j^(k-1)
            Matrix sum2 = Matrix.Null(sizeofBlock, 1);
            for (int j = i + order; j != to; j += order) {
                Matrix a_ij = A.getBlock(i*sizeofBlock, (i+1)*sizeofBlock, sizeofBlock*j, sizeofBlock*(j+1));
                Matrix x_j = xk.getBlock(sizeofBlock*j, sizeofBlock*(j+1), 0, 1);   //!!!in lectures - x_kmin1
                sum2.minus(a_ij.times(x_j));
            }

            //count block i of xk
            Matrix firstElem = x_kmin1_i.mulOnConstant(1 - w);
            Matrix secondElem = sum1.plus(sum2).plus(b_i);

            Matrix xk_i = firstElem.plus(a_ii.degMin1().times(secondElem).mulOnConstant(w));

            //copy block xk_i to xk
            for (int z = 0; z < sizeofBlock; z++)
                xk.setElement(i*sizeofBlock+z, 0, xk_i.GetElement(z, 0));
        }
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

            //!!!must be first if, but it doesn't work
            //if (xk.minus(prevX).norm2() < accuracy.norm2()) {//!!!may be it's better to count |xk| and |precX| first
            //if (numOfSteps > 5){//
            if (b.minus(A.times(xk)).norm2() < accuracy) {//accuracy.norm2()) {
                goOn = false;
                System.out.println("check: b-Ax=");
                b.minus(A.times(xk)).show();
            }
            prevX = new Matrix(xk);

            //if algorithm doesn't work we finish manually - in case of infinity
            if (numOfSteps > 5000) {
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
         double accur = 0.01;

         for (double w = 0; w < 2; w += 0.01) {
             this.w = w;

             System.out.print("steps: (" + w + ","+ count(accur)+");\n ");

         }
    }

    public static void main(String []args) {
        System.out.println("try it");

        int n = 21;

        //Matrix A = new Matrix(new double[][]{{1,0, 0}, {0,1,0}, {0,0,1}});//
        //Matrix b = new Matrix(new double[][]{{1}, {2}, {3}});//

        Matrix A = Matrix.randomsym(n, n);
        Matrix b = Matrix.random(n,1);

        //Matrix approxim = new Matrix(new double[][]{{10}, {10}, {10}});//
        Matrix approxim =  Matrix.random(n, 1);


        SSORblock sb = new SSORblock(A, b, approxim, 1, 3);
        sb.printResultChangingW(null);

    }
}
