package Methods_functions;

public class SLAEIterationalMethods {
    private static int last_method_iterations = 0;
    public static double[][] SimpleIterationMethod(double[][] A, double[][] b, double accurateness) throws Exception {
        last_method_iterations = 0;
        double coefficient = 1/MathFunctions.norm1(A);
        double[][] B = MathFunctions.subtractMatrices(MathFunctions.getNeutralMatrix(A.length),
                                                      MathFunctions.multiplyMatrix(A, coefficient));
        double[][] c = MathFunctions.multiplyMatrix(b, coefficient);
        double normB = Math.min(MathFunctions.norm1(B),MathFunctions.normInfinity(B));

        double[][] x = new double[b.length][1];
        for (int i = 0; i < x.length; i++) {
            x[i][0] = 0;
        }

        if(normB >= 1) {
            b = MathFunctions.multiplyMatrices(MathFunctions.transpondMatrix(A), b);
            A = MathFunctions.multiplyMatrices(MathFunctions.transpondMatrix(A), A);

            coefficient = 1/MathFunctions.norm1(A);
            B = MathFunctions.subtractMatrices(MathFunctions.getNeutralMatrix(A.length),
                    MathFunctions.multiplyMatrix(A, coefficient));
            c = MathFunctions.multiplyMatrix(b, coefficient);
            normB = Math.min(MathFunctions.norm1(B),MathFunctions.normInfinity(B));
        }

        if(normB >=1) {
            while(MathFunctions.norm1(MathFunctions.subtractMatrices(MathFunctions.multiplyMatrices(A, x), b)) >= accurateness) {
                x = MathFunctions.sumMatrices(MathFunctions.multiplyMatrices(B, x), c);
                last_method_iterations += 1;
            }
        }
        else {
            double accuracyCoefficient = normB/(1-normB);
            double[][] newX = MathFunctions.sumMatrices(MathFunctions.multiplyMatrices(B, x), c);

            while(accuracyCoefficient*MathFunctions.norm1(MathFunctions.subtractMatrices(newX, x)) >= accurateness &&
                  accuracyCoefficient*MathFunctions.normInfinity(MathFunctions.subtractMatrices(newX, x)) >= accurateness) {
                x = newX;
                newX = MathFunctions.sumMatrices(MathFunctions.multiplyMatrices(B, x), c);
                last_method_iterations += 1;
            }
            x = newX;
        }
        return x;
    }
    public static double[][] ZeidelMethod(double[][] A, double[][] b, double accurateness) throws Exception{

        last_method_iterations = 0;
        if(MathFunctions.isDiagonalMajority(A) <= 0) {
            b = MathFunctions.multiplyMatrices(MathFunctions.transpondMatrix(A), b);
            A = MathFunctions.multiplyMatrices(MathFunctions.transpondMatrix(A), A);
        }

        double[][] x = new double[b.length][1];
        for (int i = 0; i < x.length; i++) {
            x[i][0] = 0;
        }
        double[][] c = new double[A.length][A[0].length];
        double[] d = new double[b.length];

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                if(i!=j) {
                    c[i][j] = -1 * A[i][j] / A[i][i];
                } else {
                    c[i][j] = 0;
                }
            }
        }
        for (int i = 0; i < b.length; i++) {
            d[i] = b[i][0]/A[i][i];
        }

        double[][] newX;
        while(MathFunctions.norm1(MathFunctions.subtractMatrices(MathFunctions.multiplyMatrices(A, x), b)) >= accurateness) {
            newX = new double[A.length][1];
            for (int i = 0; i < newX.length; i++) {
                newX[i][0] = d[i];
                for (int j = 0; j < i; j++) {
                    newX[i][0] += c[i][j]*newX[j][0];
                }
                for (int j = i; j < newX.length; j++) {
                    newX[i][0] += c[i][j]*x[j][0];
                }
            }
            x = newX;
            last_method_iterations += 1;
        }

        return x;
    }
    public static int getLastMethodIterations() {
        return last_method_iterations;
    }
}
