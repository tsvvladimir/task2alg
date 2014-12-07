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
        boolean flag  = false;//true;//

        if (flag) System.out.println("--------------------------------------------\nprevX=");
        if (flag) prevX.show();
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
            if (flag) {System.out.print("\nx^k-1_i:\n"); x_kmin1_i.show();}
            //get cur block i of b
            Matrix b_i = b.getBlock(sizeofBlock*i, sizeofBlock*(i+1), 0, 1);
            if (flag){ System.out.print("\nbi:\n"); b_i.show();}
            //get block ii of A
            Matrix a_ii = A.getBlock(sizeofBlock*i, sizeofBlock*(i+1), sizeofBlock*i, sizeofBlock*(i+1));
            if (flag) {System.out.print("\naii:\n"); a_ii.show();}

            //count sum on all j < i of a_ij * x_j^(k)
            Matrix sum1 = Matrix.Null(sizeofBlock, 1);
            for (int j = from; j != i; j += order) {
                Matrix a_ij = A.getBlock(i*sizeofBlock, (i+1)*sizeofBlock, sizeofBlock*j, sizeofBlock*(j+1));
                //System.out.print("\na ij:\n"); a_ij.show();
                Matrix x_j = xk.getBlock(sizeofBlock*j, sizeofBlock*(j+1), 0, 1);    //!!!in lectures - xk
                //System.out.print("\nx j:\n"); x_j.show();
                //a_ij.times(x_j).show();
                sum1 = sum1.minus(a_ij.times(x_j));
                //System.out.print("\nsum1:\n"); sum1.show();

            }
            if (flag){ System.out.print("\nsum1:\n"); sum1.show();}

            //count sum on all j > i of a_ij * x_j^(k-1)
            Matrix sum2 = Matrix.Null(sizeofBlock, 1);
            for (int j = i + order; j != to; j += order) {
                Matrix a_ij = A.getBlock(i*sizeofBlock, (i+1)*sizeofBlock, sizeofBlock*j, sizeofBlock*(j+1));
                //System.out.print("a["+i+"]["+j+"]  must be equal to " + A.GetElement(i, j) + ";\n");
                //System.out.print("\na ij:\n"); a_ij.show();
                Matrix x_j = prevX.getBlock(sizeofBlock*j, sizeofBlock*(j+1), 0, 1);   //!!!in lectures - x_kmin1
                //System.out.print("\nx j:\n"); x_j.show();
                sum2 = sum2.minus(a_ij.times(x_j));
                //sum2.show();
            }
            if (flag){ System.out.print("\nsum2:\n"); sum2.show();}

            //count block i of xk
            Matrix firstElem = x_kmin1_i.mulOnConstant(1 - w);
            if (flag) {System.out.print("I:\n"); firstElem.show();}
            Matrix secondElem = sum1.plus(sum2).plus(b_i);
            if (flag) {System.out.print("II:\n"); secondElem.show();}

            Matrix xk_i = firstElem.plus(a_ii.degMin1().times(secondElem).mulOnConstant(w));
            if (flag) { System.out.print("\nxk["+i+"]:\n"); xk_i.show(); }
            //System.out.print("\nx^k i:\n"); xk_i.show();

            //copy block xk_i to xk
            for (int z = 0; z < sizeofBlock; z++)
                xk.setElement(i*sizeofBlock+z, 0, xk_i.GetElement(z, 0));
            if (flag) { System.out.print("\n--iter finished--\n xk:\n"); xk.show();}
            if (flag) System.out.print("\n");
        }
        if (flag) System.out.println("xk=");
        if (flag) xk.show();
        return xk;
    }

    public int count(double accuracy) {
        boolean flag = false;
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

            prevX = new Matrix(xk);

            //!!!must be first if, but it doesn't work
            //if (xk.minus(prevX).norm2() < accuracy.norm2()) {//!!!may be it's better to count |xk| and |precX| first
            //if (numOfSteps > 5){//


            if (b.minus(A.times(xk)).norm2() < accuracy) {//accuracy.norm2()) {
                goOn = false;
                if (flag) {
                    System.out.println("check: b-Ax=");
                    System.out.println("check: b=");
                    b.show();
                    System.out.println();
                    Matrix must_be_B = A.times(xk);
                    System.out.println("check: Ax=");
                    must_be_B.show();
                }
                //b.minus(A.times(xk)).show();
                return numOfSteps;
            }

            //b.minus(A.times(xk)).show();
            //System.out.println();
            //if algorithm doesn't work we finish manually - in case of infinity
            if (numOfSteps > 1000) {

                /*System.out.println("ERR: check: b-Ax=");
                System.out.println("check: b=");
                b.show();
                System.out.println();
                Matrix must_be_B = A.times(xk);

                 System.out.println("A");
                 A.show();
                 System.out.println("xk");
                xk.show();

                System.out.println("check: Ax=");
                must_be_B.show();
                //b.minus(A.times(xk)).show();
                //xk.show();
                System.out.println("check: number of steps is more then "+numOfSteps+" , is it ok?;)");*/
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
        double accur = 0.0000001;

        for (double w = 0; w < 2; w += 0.01) {
            this.w = w;

            System.out.print(" (" + w + ";"+ count(accur)+")\n ");

        }
    }

    public static void main(String []args) {
        System.out.println("try it");

        int n = 21;

        Matrix A = new Matrix(new double[][]{{(0.9967460596589041), (0.024048068951804447), (0.5338653126932894)},
                {(0.024048068951804447), (0.8406518900095598), (0.17713536999280655)},
                {(0.5338653126932894), (0.17713536999280655), (0.38255269008212744)}}); //
        A = new Matrix    (new double[][]{{4,3, 0}, {3,4,-1}, {0,-1,4}});//
        Matrix b = new Matrix(new double[][]{{(0.961186269398937)},
                {(0.4675750079428419)},
                {(0.019207511703432045)}}); //
        b = new Matrix(new double[][]{{24}, {30}, {-24}});//



        Matrix help1 = new Matrix(n, 1);
        for(int i = 0; i < n; i++)
            help1.data[i][0] = 1;
        //Matrix b = A.times(help1);//
        b = Matrix.random(n,1);
        A = Matrix.randomsym(n, n);
        System.out.print("A, b generated:\n");
        A.show();
        System.out.print("b: \n");
        b.show();


        //Matrix approxim = new Matrix(new double[][]{{10}, {10}, {10}});//
        Matrix approxim =  Matrix.random(n, 1);
        approxim.show();


        SSORblock sb = new SSORblock(A, b, approxim, 1, 3);
        sb.printResultChangingW(null);

    }
}
