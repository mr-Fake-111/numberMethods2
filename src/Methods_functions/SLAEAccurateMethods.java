package Methods_functions;

import java.util.Arrays;

public class SLAEAccurateMethods {
    public static double[][] GaussMethod(double[][] A, double[][] b) throws Exception {
        int n = b.length;
        double[][] P = MathFunctions.getNeutralMatrix(n);
        double[][] M = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] = A[i][j];
            }
        } //M == A

        for (int i = 0; i < n; i++) {

            double max_i = Math.abs(A[i][i]);
            int max_i_index = i;
            for (int j = i; j < n; j++) {
                if(Math.abs(A[j][i]) > max_i){
                    max_i_index = j;
                    max_i = Math.abs(A[j][i]);
                }
            }
            M = MathFunctions.replaceRows(M, i, max_i_index);
            P = MathFunctions.replaceRows(P, i, max_i_index);

            for (int j = i+1; j < n; j++) {
                M[j][i] = M[j][i]/M[i][i];
                for (int k = i+1; k < n; k++) {
                    M[j][k] = M[j][k] - M[j][i]*M[i][k];
                }
            }
            
        }
        double[][] L = MathFunctions.getNeutralMatrix(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                L[i][j] = M[i][j];
            }
        }
        double[][] U = MathFunctions.subtractMatrices(MathFunctions.sumMatrices(M, MathFunctions.getNeutralMatrix(n)), L);

        b = MathFunctions.multiplyMatrices(P, b);
        double[][] y = MathFunctions.reversedGaussMove(L, b);
        double[][] x = MathFunctions.straightGaussMove(U, y);

        return x;
    }

    public static double[][] ReflectionMethod(double[][] A, double[][] b) throws Exception {
        int n = A.length;
        double[][] Q = MathFunctions.getNeutralMatrix(n);
        double[][] R = new double[n][n];
        for (int i = 0; i < n; i++) {
            R[i] = Arrays.copyOf(A[i], n);
        }

        for (int i = 0; i < n-1; i++) {

            double[][] z = new double[n][1];
            for (int j = 0; j < n; j++) {
                z[j][0] = 0;
            }
            z[i][0] = 1;

            double[][] y = new double[n][1];
            for (int j = 0; j < i; j++) {
                y[j][0] = 0;
            }
            for (int j = i; j < n; j++) {
                y[j][0] = R[j][i];
            }

            double[][] w = MathFunctions.subtractMatrices(y, MathFunctions.multiplyMatrix(z, MathFunctions.vectorNorm2(y)));
            w = MathFunctions.multiplyMatrix(w, 1.0/MathFunctions.vectorNorm2(w));
            double[][] transpondW = new double[1][n];
            for (int j = 0; j < n; j++) {
                transpondW[0][j] = w[j][0];
            }

            double[][] Q_i = MathFunctions.subtractMatrices(MathFunctions.getNeutralMatrix(n),
                                                            MathFunctions.multiplyMatrix(MathFunctions.multiplyMatrices(w, transpondW), 2));
            double[][] R_i = new double[n][n];
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    R_i[i][i] = 0;
                }
            }
            for (int j = i; j < n; j++) {
                for (int k = i; k < n; k++) {
                    R_i[j][k] = R[j][k];
                }
            }
            R_i = MathFunctions.multiplyMatrices(Q_i, R);
            for (int j = i; j < n; j++) {
                for (int k = i; k < n; k++) {
                    R[j][k] = R_i[j][k];
                }
            }

            Q = MathFunctions.multiplyMatrices(Q, Q_i);
        }

        b = MathFunctions.multiplyMatrices(MathFunctions.transpondMatrix(Q), b);
        double[][] x = MathFunctions.straightGaussMove(R, b);

        return x;
    }
}
