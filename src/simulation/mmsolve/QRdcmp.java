package simulation.mmsolve;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

/**
 * Performs QR matrix decomposition and solved linear SOEs
 * Used for Broyden's method.
 * Based on numerical recipes 3rd version and equations are references there
 */
public class QRdcmp {

    // Object for QR decomposition of a matrix A, and related functions.
    /** linear dimension of the matrix */
    private int n;
    /** Stores the decomposition, Q^T and R */
    private final double[][] qt, r;
    /** Indicates whether A is singular. */
    private boolean sing;

    // QRdcmp(MatDoub_I &a); Constructor from A.
    // void solve(VecDoub_I &b, VecDoub_O &x); Solve A  x D b for x.
    // void qtmult(VecDoub_I &b, VecDoub_O &x); Multiply QT  b D x.
    // void rsolve(VecDoub_I &b, VecDoub_O &x); Solve R  x D b for x.
    // void update(VecDoub_I &u, VecDoub_I &v); See next subsection.
    // void rotate(final int  i, final double a, final double b); Used by update.

    /** Returns the matrix Q^T*/
    public double[][] qt() {
        return qt;
    }

    /** Returns the matrix R */
    public double[][] r() {
        return r;
    }

    /** Returns whether the matrix A is singular */
    public boolean sing() {
        return sing;
    }

    /**
     * Copies an NxN matrix
     * @param m - the matrix to be copied
     * @return copy of the given matrix
     */
    public static final double[][] doub_mat(final double[][] m) {
        final double[][] r = new double[m.length][];
        for (int i = 0; i < m.length; i++) {
            r[i] = new double[m[i].length];
            System.arraycopy(m[i], 0, r[i], 0, m[i].length);
        }
        return r;
    }

    /**
     * Copies a vector (of doubles)
     * @param a - the vector to be copied
     * @return copy of tje given vector
     */
    public static final double[] doub_vec(final double[] a) {
        final double[] r = new double[a.length];
        System.arraycopy(a, 0, r, 0, a.length);
        return r;
    }

    /**
     * Returns a value with the same magnitude as a and the same sign as b.
     */
    public static double SIGN(final double a, final double b) {
        return (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
    }

    /**
     * Construct the QR decomposition of a[0..n-1][0..n-1]. The upper
     * triangular matrix R and the transpose of the orthogonal matrix
     * Q are stored. sing is set to true if a singularity is
     * encountered during the decomposition, but the decomposition is
     * still completed in this case; otherwise it is set to false.
     * @param a - the matrix to be decomposed
     */
    public QRdcmp(final double[][] a) {
        n = a.length;
        qt = new double[n][n];
        r = doub_mat(a);
        sing = (false);
        int i, j, k;
        final double[] c = new double[n], d = new double[n];
        double scale, sigma, sum, tau;
        for (k = 0; k < n - 1; k++) {
            scale = 0.0;
            for (i = k; i < n; i++)
                scale = Math.max(scale, abs(r[i][k]));
            if (scale == 0.0) { // Singular case.
                sing = true;
                c[k] = d[k] = 0.0;
            } else { // Form Q_k and Q_k*A.
                for (i = k; i < n; i++)
                    r[i][k] /= scale;
                for (sum = 0.0, i = k; i < n; i++)
                    sum += r[i][k]*r[i][k];
                sigma = SIGN(sqrt(sum), r[k][k]);
                r[k][k] += sigma;
                c[k] = sigma * r[k][k];
                d[k] = -scale * sigma;
                for (j = k + 1; j < n; j++) {
                    for (sum = 0.0, i = k; i < n; i++)
                        sum += r[i][k] * r[i][j];
                    tau = sum / c[k];
                    for (i = k; i < n; i++)
                        r[i][j] -= tau * r[i][k];
                }
            }
        }
        d[n - 1] = r[n - 1][n - 1];
        if (d[n - 1] == 0.0)
            sing = true;
        for (i = 0; i < n; i++) { // Form Q^T explicitly.
            for (j = 0; j < n; j++)
                qt[i][j] = 0.0;
            qt[i][i] = 1.0;
        }
        for (k = 0; k < n - 1; k++) {
            if (c[k] != 0.0) {
                for (j = 0; j < n; j++) {
                    sum = 0.0;
                    for (i = k; i < n; i++)
                        sum += r[i][k] * qt[i][j];
                    sum /= c[k];
                    for (i = k; i < n; i++)
                        qt[i][j] -= sum * r[i][k];
                }
            }
        }
        for (i = 0; i < n; i++) { // Form R explicitly.
            r[i][i] = d[i];
            for (j = 0; j < i; j++)
                r[i][j] = 0.0;
        }
    }

    // The next set of member functions is used to solve linear
    // systems. In many applications only the part (2.10.4) of
    // the algorithm is needed, so we put in separate routines
    // the multiplication QT b and the backsubstitution on R.

    /**
     * Solve the set of n linear equations A*x = b.
     * @param b - input as the right-hand side vector b[0..n-1]
     * @param x - returned as the solution vector x[0..n-1]
     * @throws Exception
     */
    public void solve(final double[] b, final double[] x) throws Exception {
        qtmult(b, x); // Form Q^T*b.
        rsolve(x, x); // Solve R*x = Q^T*b.
    }

    /**
     * Multiply Q^T*b and put the result in x. Since Q is
     * orthogonal, this is equivalent to solving
     * Q*x=b for x.
     * @param b - input as the right-hand side vector b[0..n-1]
     * @param x - returned as the solution vector x[0..n-1]
     */
    public void qtmult(final double[] b, final double[] x) {
        int i, j;
        double sum;
        for (i = 0; i < n; i++) {
            sum = 0.;
            for (j = 0; j < n; j++)
                sum += qt[i][j] * b[j];
            x[i] = sum;
        }
    }

    /**
     * Solve the triangular set of n linear equations R*x=b.
     * @param b b[0..n-1] is input as the right-hand side vector
     * @param x x[0..n-1] is returned as the solution vector.
     * @throws Exception
     */
    public void rsolve(final double[] b, final double[] x) throws Exception {
        int i, j;
        double sum;
        if (sing)
            throw new Exception("attempting solve in a singular QR");
        for (i = n - 1; i >= 0; i--) {
            sum = b[i];
            for (j = i + 1; j < n; j++)
                sum -= r[i][j] * x[j];
            x[i] = sum / r[i][i];
        }
    }

    // From Numerical Recipes (3rd edition) section 2.10.1:
    // Some numerical algorithms involve solving a succession of
    // linear systems each of which differs only slightly from its
    // predecessor. Instead of doing O(N^3) operations each time to
    // solve the equations from scratch, one can often update a
    // matrix factorization in O(N^2) operations and use the new
    // factorization to solve the next set of linear equations.
    // The LU decomposition is complicated to update because of
    // pivoting. However, QR turns out to be quite simple for a
    // very common kind of update,
    // A -> A + s \otimes t (2.10.7)
    // (compare equation 2.7.1). In practice it is more convenient
    // to work with the equivalent form
    // A = Q*R -> A' = Q'*R' = Q*(R + u \otimes v) (2.10.8)
    // One can go back and forth between equations (2.10.7) and
    // (2.10.8) using the fact that Q is orthogonal, giving
    // t = v and either s = Q*u or u = Q^T * s (2.10.9)
    // The algorithm has two phases. In the first we apply N-1
    // Jacobi rotations to reduce R + u \otimes v to upper
    // Hessenberg form. Another N - 1 Jacobi rotations transform
    // this upper Hessenberg matrix to the new upper triangular
    // matrix R'. The matrix Q' is simply the product of Q with
    // the 2(N-1) Jacobi rotations. In applications we usually
    // want Q^T, so the algorithm is arranged to work with this
    // matrix (which is stored in the QRdcmp object) instead of
    // with Q.

    /**
     * Starting from the stored QR decomposition A=Q*R, update
     * it to be the QR decomposition of the matrix Q*(R + u \otimes v).
     * @param u - left vector in outer product that is added to R
     * @param v - right vector in outer product that is added to R
     */
    public void update(final double[] u, final double[] v) {
        int i, k;
        final double[] w = doub_vec(u);
        for (k = n - 1; k >= 0; k--)
            // Find largest k such that u[k] ₪ 0.
            if (w[k] != 0.0)
                break;
        if (k < 0)
            k = 0;
        for (i = k - 1; i >= 0; i--) { // Transform RCu?v to upper Hessenberg.
            rotate(i, w[i], -w[i + 1]);
            if (w[i] == 0.0)
                w[i] = abs(w[i + 1]);
            else if (abs(w[i]) > abs(w[i + 1]))
                w[i] = abs(w[i]) * sqrt(1.0 + (w[i + 1] / w[i])*(w[i + 1] / w[i]));
            else
                w[i] = abs(w[i + 1]) * sqrt(1.0 + (w[i] / w[i + 1])*(w[i] / w[i + 1]));
        }
        for (i = 0; i < n; i++)
            r[0][i] += w[0] * v[i];
        for (i = 0; i < k; i++)
            // Transform upper Hessenberg matrix to upper tri
            rotate(i, r[i][i], -r[i + 1][i]); // angular.
        for (i = 0; i < n; i++)
            if (r[i][i] == 0.0)
                sing = true;
    }

    /**
     * Utility used by {@link QRdcmp#update(double[], double[])}. Given matrices r[0..n-1][0..n-1]
     * and qt[0..n-1][0..n-1], carry out a Jacobi rotation on rows
     * i and i+1 of each matrix.
     * the rotation: cos(\theta) = a/\sqrt{a^2 + b^2}, sin(\theta) = b/\sqrt{a^2 + b^2}
     * @param i - row to rotate
     * @param a - rotation parameter
     * @param b - rotation parameter
     */
    public void rotate(final int i, final double a, final double b) {
        int j;
        double c, fact, s, w, y;
        if (a == 0.0) { // Avoid unnecessary overflow or underflow.
            c = 0.0;
            s = (b >= 0.0 ? 1.0 : -1.0);
        } else if (abs(a) > abs(b)) {
            fact = b / a;
            c = SIGN(1.0 / sqrt(1.0 + (fact * fact)), a);
            s = fact * c;
        } else {
            fact = a / b;
            s = SIGN(1.0 / sqrt(1.0 + (fact * fact)), b);
            c = fact * s;
        }
        for (j = i; j < n; j++) { // Premultiply r by Jacobi rotation.
            y = r[i][j];
            w = r[i + 1][j];
            r[i][j] = c * y - s * w;
            r[i + 1][j] = s * y + c * w;
        }
        for (j = 0; j < n; j++) { // Premultiply qt by Jacobi rotation.
            y = qt[i][j];
            w = qt[i + 1][j];
            qt[i][j] = c * y - s * w;
            qt[i + 1][j] = s * y + c * w;
        }
    }

}


