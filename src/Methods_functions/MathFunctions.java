package Methods_functions;


import java.util.Arrays;

public class MathFunctions {
    private static double sqrt(double val, double prevRes) {
        double newRes = 0.5 * (prevRes + val / prevRes);
        if (Math.abs(newRes - prevRes) < Math.pow(10, -9)) {
            return newRes;
        }
        else {
            return sqrt(val, newRes);
        }
    }
    public static double sqrt(double val) {
        if(val != 1) {
            return sqrt(val, val);
        }
        return 1;
    }


    public static double isDiagonalMajority(double[][] matrix) {
        int n = matrix.length;
        double minMajority = Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            double val = Math.abs(matrix[i][i]);

            for (int j = 0; j < n; j++) {
                if(j != i) {
                    val -= Math.abs(matrix[i][j]);
                }
            }

            if(val <= 0) return -1;
            if(val < minMajority) minMajority = val;
        }
        return minMajority;
    }
    public static double[][] transpondMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] transpondedMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                transpondedMatrix[i][j] = matrix[j][i];
            }
        }
        return transpondedMatrix;
    }
    public static double norm1(double[][] matrix) {
        double maxColumnValue = 0;
        for (int j = 0; j < matrix[0].length; j++) {
            double columnSum = 0;
            for (int i = 0; i < matrix.length; i++) {
                columnSum += Math.abs(matrix[i][j]);
            }
            if(columnSum > maxColumnValue) maxColumnValue = columnSum;
        }
        return maxColumnValue;
    }
    public static double normInfinity(double[][] matrix) {
        int n = matrix.length;
        if(matrix[0].length > 1) {
            double maxRowValue = 0;
            for (int i = 0; i < n; i++) {
                double rowSum = 0;
                for (int j = 0; j < n; j++) {
                    rowSum += Math.abs(matrix[i][j]);
                }
                if (rowSum > maxRowValue) maxRowValue = rowSum;
            }
            return maxRowValue;
        } else {
            double maxValue = 0;
            for (int i = 0; i < n; i++) {
                if(maxValue < Math.abs(matrix[i][0])) maxValue = Math.abs(matrix[i][0]);
            }
            return maxValue;
        }
    }
    public static double[][] sumMatrices(double[][] m1, double[][] m2) throws Exception {
        if(m1.length != m2.length || m1[0].length != m2[0].length) throw new Exception("matrices are not compatible for sum");
        double[][] resMatrix = new double[m1.length][m1[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m1[0].length; j++) {
                resMatrix[i][j] = m1[i][j] + m2[i][j];
            }
        }
        return resMatrix;
    }
    public static double[][] subtractMatrices(double[][] m1, double[][] m2) throws Exception{
        return sumMatrices(m1, multiplyMatrix(m2, -1));
    }
    public static double[][] multiplyMatrices(double[][] m1, double[][] m2) throws Exception {
        if(m1[0].length != m2.length) throw new Exception("matrices are not compatible for multiplication");
        double[][] resMatrix = new double[m1.length][m2[0].length];
        for (int i = 0; i < resMatrix.length ; i++) {
            for (int j = 0; j < resMatrix[0].length; j++) {
                double value = 0;
                for (int k = 0; k < m2.length; k++) {
                    value += m1[i][k]*m2[k][j];
                }
                resMatrix[i][j] = value;
            }
        }
        return resMatrix;
    }
    public static double[][] multiplyMatrix(double[][] m, double value) {
        double[][] resMatrix = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                resMatrix[i][j] = value * m[i][j];
            }
        }
        return resMatrix;
    }
    public static double[][] getNeutralMatrix(int n) {
       double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if(i == j) matrix[i][j] = 1;
                else matrix[i][j] = 0;
            }
        }
        return matrix;
    }

    public static double[][] getBadMatrix(int n, double epsilon) {
        double[][] badMatrix = getNeutralMatrix(n);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                badMatrix[i][j] = -1 - epsilon*(10);
            }
        }
        for (int i = 0; i < n; i++) {
            badMatrix[i][i] += epsilon*(10);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                badMatrix[i][j] = epsilon*( 10);
            }
        }
        return badMatrix;
    }

    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length ; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }
    
    public static double[][] straightGaussMove(double[][] matrix, double[][] b) {
        double[][] x = new double[b.length][1];
        int n = x.length;
        for (int i = n -1; i >= 0; i--) {
            x[i][0] = b[i][0];
            for (int j = i+1; j < n ; j++) {
                x[i][0] -= x[j][0]*matrix[i][j];
            }
            x[i][0] /= matrix[i][i];
        }
        
        return x;
    }
    public static double[][] reversedGaussMove(double[][] matrix, double[][] b) {
        double[][] x = new double[b.length][1];
        int n = x.length;
        for (int i = 0; i < n; i++) {
            x[i][0] = b[i][0];
            for (int j = i-1; j >=0 ; j--) {
                x[i][0] -= x[j][0]*matrix[i][j];
            }
            x[i][0] /= matrix[i][i];
        }

        return x;
    }

    public static double[][] replaceRows(double[][] matrix, int i, int j) {
        double[] help = Arrays.copyOf(matrix[i], matrix[i].length);
        matrix[i] = matrix[j];
        matrix[j] = help;
        return matrix;
    }

    public static double vectorNorm2(double[][] vector) throws Exception{
        if(vector[0].length != 1) throw new Exception("Matrix is not a vector-column");
        double res = 0;
        for (int i = 0; i < vector.length; i++) {
            res += Math.pow(vector[i][0], 2);
        }
        res = Math.sqrt(res);
        return res;
    }
}


