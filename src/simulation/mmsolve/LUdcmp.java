
package simulation.mmsolve;

import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.math.BigDecimal;

import static java.lang.Math.abs;

/**
 * Performs LU matrix decomposition and solved linear SOEs
 * Used for Newton's method.
 * Based on numerical recipes 3rd version and equations are references there
 */
public class LUdcmp {

    /** linear dimension of the matrix */
    private int n;
    /** Stores the decomposition */
    private final double[][] lu;
    /** Stores the permutation */
    private final int[] indx;
    /** Used by det */
    private double d;
    // private void solve(VecDoub_I &b, VecDoub_O &x); // Solve for a single
    // right-hand side.
    // private void solve(MatDoub_I &b, MatDoub_O &x); // Solve for multiple
    // right-hand sides.
    // private void inverse(MatDoub_O &ainv); // Calculate matrix inverse A1.
    // private Doub det(); // Return determinant of A.
    // private void mprove(VecDoub_I &b, VecDoub_IO &x); // Discussed in 2.5.
    private final double[][] aref; // Used only by mprove.

    /**
     * Constructs a LUdcmp object.
     * Given a matrix a[0..n-1][0..n-1], this routine replaces it by the
     * LU decomposition of a rowwise permutation of itself.
     * On output, it is arranged as in equation (2.3.14);
     * indx[0..n-1] is an output vector that records the row permutation
     * effected by the partial pivoting; d is output as +/- 1 depending on
     * whether the number of row interchanges was even or odd, respectively.
     * This routine is used in combination with solve to solve linear
     * equations or invert a matrix.
     * @param a - the matrix A to decompose
     * @throws SingularMatrixException
     */
    public LUdcmp(final double[][] a) throws SingularMatrixException {
        n = a.length;
        lu = QRdcmp.doub_mat(a);
        aref = QRdcmp.doub_mat(a);
        indx = new int[n];

        final double TINY = 1.0e-40; // A small number.
        int i, imax = 0, j, k;
        double big, temp;
        final double[] vv = new double[n];  // vv stores the implicit scaling of each
                                            // row.
        d = 1.0; // No row interchanges yet.
        for (i = 0; i < n; i++) { // Loop over rows to get the implicit scaling
                                  // information
            big = 0.0;
            for (j = 0; j < n; j++)
                if ((temp = abs(lu[i][j])) > big)
                    big = temp;
            if (big == 0.0)
                throw new SingularMatrixException();
            // No nonzero largest element.
            vv[i] = 1.0 / big; // Save the scaling.
        }
        for (k = 0; k < n; k++) { // This is the outermost kij loop.
            big = 0.0; // Initialize for the search for largest pivot element.
            for (i = k; i < n; i++) {
                temp = vv[i] * abs(lu[i][k]);
                if (temp > big) { // Is the figure of merit for the pivot better
                                  // than the best so far?
                    big = temp;
                    imax = i;
                }
            }
            if (k != imax) {                // Do we need to interchange rows?
                for (j = 0; j < n; j++) {   // Yes, do so...
                    temp = lu[imax][j];
                    lu[imax][j] = lu[k][j];
                    lu[k][j] = temp;
                }
                d = -d;             // ...and change the parity of d.
                vv[imax] = vv[k];   // Also interchange the scale factor.
            }
            indx[k] = imax;
            if (lu[k][k] == 0.0)
                lu[k][k] = TINY;
            // If the pivot element is zero, the matrix is singular (at least
            // to the precision of the algorithm). For some applications on
            // singular matrices, it is desirable to substitute TINY for zero.
            for (i = k + 1; i < n; i++) {
                temp = lu[i][k] /= lu[k][k]; // Divide by the pivot element.
                for (j = k + 1; j < n; j++)
                    // Innermost loop: reduce remaining submatrix.
                    lu[i][j] -= temp * lu[k][j];
            }
        }
    }

    // Once the LUdcmp object is constructed, two functions implementing
    // equations (2.3.6) and (2.3.7) are available for solving linear
    // equations. The first solves a single right-hand side vector b for
    // a solution vector x. The second simultaneously solves multiple
    // right-hand vectors, arranged as the columns of a matrix B. In
    // otherwords, it calculates the matrix A^{-1}*B.

    /**
     * Solves the set of n linear equations A*x = b using the stored
     * LU decomposition of A. b and x may reference the
     * same vector, in which case the solution overwrites
     * the input. This routine takes into account the possibility that
     * b will begin with many zero elements, so it is efficient for use
     * in matrix inversion.
     * @param b - input as the right-hand side vector b[0..n-1]
     * @param x - returns the solution vector x
     * @throws MatrixDimensionMismatchException
     */
    public void solve(final double[] b, final double[] x) throws MatrixDimensionMismatchException {
        int i, ii = 0, ip, j;
        double sum;
        if (b.length != n || x.length != n)
            throw new MatrixDimensionMismatchException(b.length,x.length,n,n);
        for (i = 0; i < n; i++)
            x[i] = b[i];
        for (i = 0; i < n; i++) {   // When ii is set to a positive value, it will
                                    // become the
                                    // index of the first nonvanishing element of b. We now
                                    // do the forward substitution, equation (2.3.6). The
                                    // only new wrinkle is to unscramble the permutation
                                    // as we go.
            ip = indx[i];
            sum = x[ip];
            x[ip] = x[i];
            if (ii != 0)
                for (j = ii - 1; j < i; j++)
                    sum -= lu[i][j] * x[j];
            else if (sum != 0.0)    // A nonzero element was encountered, so from
                                    // now on we
                ii = i + 1;         // will have to do the sums in the loop above.
            x[i] = sum;
        }
        for (i = n - 1; i >= 0; i--) { // Now we do the backsubstitution,
                                       // equation (2.3.7).
            sum = x[i];
            for (j = i + 1; j < n; j++)
                sum -= lu[i][j] * x[j];
            x[i] = sum / lu[i][i]; // Store a component of the solution vector
                                   // X.
        } // All done!
    }

    /**
     * Solves m sets of n linear equations A*X = B using the stored
     * LU decomposition of A. The matrix b[0..n-1][0..m-1] inputs the
     * right-hand sides, while x[0..n-1][0..m-1] returns the solution
     * A^{-1}*B. b and x may reference the same matrix, in which case the
     * solution overwrites the input.
     * @param b - the matrix b[0..n-1][0..m-1] that inputs the right-hand sides
     * @param x - the matrix x[0..n-1][0..m-1] that stores the solution A^{-1}*B
     * @throws MatrixDimensionMismatchException
     */
    public void solve(final double[][] b, final double[][] x) throws MatrixDimensionMismatchException {
        int i, j, m = b[0].length;
        if (b.length != n || x.length != n ||
                b[0].length != x[0].length)
            throw new MatrixDimensionMismatchException(b.length, b[0].length, n, n);
        final double[] xx = new double[n];
        for (j = 0; j < m; j++) { // Copy and solve each column in turn.
            for (i = 0; i < n; i++)
                xx[i] = b[i][j];
            solve(xx, xx);
            for (i = 0; i < n; i++)
                x[i][j] = xx[i];
        }
    }

    /**
     * Using the stored LU decomposition, return in ainv the matrix inverse A^{-1}
     * @param ainv - stores (returns) the solution
     * @throws Exception
     */
    public void inverse(final double[][] ainv) throws Exception {
        int i, j;

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                ainv[i][j] = 0.;
            ainv[i][i] = 1.;
        }
        solve(ainv, ainv);
    }

    /**
     * Using the stored LU decomposition, calculates the determinant of the matrix A
     * @return the determinant of the matrix A which was LU decomposed
     */
    public double det() {
        double dd = d;
        for (int i = 0; i < n; i++)
            dd *= lu[i][i];
        return dd;
    }

    /**
     * Improves a solution vector x[0..n-1] of the linear set of equations
     * A*x = b. The vectors b[0..n-1] and x[0..n-1] are input. On output,
     * x[0..n-1] is modified, to an improved set of values.
     * @param b - the vector b[0..n-1] in the equation A*x = b
     * @param x - the vector x[0..n-1] in the equation A*x = b
     * @throws Exception
     */
    public void mprove(final double[] b, final double[] x) throws Exception {
        int i, j;
        final double[] r = new double[n];
        for (i = 0; i < n; i++) {                   // Calculate the right-hand side, accumulating
            BigDecimal sdp = new BigDecimal(-b[i]); // the residual in higher
                                                    // precision.
            for (j = 0; j < n; j++)
                sdp = sdp.add(new BigDecimal(aref[i][j]).multiply(new BigDecimal(x[j])));
            r[i] = sdp.doubleValue();
        }
        solve(r, r); // Solve for the error term,
        for (i = 0; i < n; i++)
            x[i] -= r[i]; // and subtract it from the old solution.
    }

}
