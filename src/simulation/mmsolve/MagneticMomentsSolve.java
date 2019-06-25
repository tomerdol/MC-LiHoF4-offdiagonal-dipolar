package simulation.mmsolve;


import simulation.montecarlo.*;
import org.apache.commons.math3.linear.MatrixDimensionMismatchException;
import org.apache.commons.math3.linear.SingularMatrixException;

import java.util.Random;

public class MagneticMomentsSolve {


    static void printMean(double[] arr){
        double mean=0;
        for (int i=0;i<arr.length;i++) mean+=arr[i];
        System.out.println(mean/arr.length);
    }

    public static double [] gauss_seidel(fi_xi f, double xx[], double alpha, int iter) throws ConvergenceException {
        int n=xx.length;
        //double xx[]=new double[n];
        double xx1[]=new double[n];
        double delta=0;
        //for(int i=0;i<n;i++) xx[i]=x[i];
        double[] g;
        iter=(int) (0.333*iter*xx.length);
        //int iter=50;
        int k;
        for(k=0;k<iter;k++) {
            for(int i=0;i<n;i++) xx1[i]=xx[i];
            delta=0;
            for(int j=0;j<n;j++) {
                g=f.func(xx);
                xx[j]=(1-alpha)*xx[j]+alpha*g[j];
                delta+=Math.abs(xx[j]-xx1[j]);
                //System.out.print("xx["+j+"] = "+xx[j]);
            }
            //System.out.println(" k = "+k+" delta = "+delta);
            //System.out.println(checkConverge(xx,f));

            if(delta< Constants.tol*1.0e-1) break;
        }
        //if (k==iter) System.out.println("reached max iteration");
        return xx;
    }

    public static double checkConverge(double[] x, fi_xi f) throws ConvergenceException{
        double converged = 0;
        double[] g = f.func(x);
        for (int i = 0; i < x.length; i++) {
            double shouldBe = g[i];
            double is = x[i];

            converged += Math.abs(is - shouldBe);
        }

        return Math.abs(converged/x.length);
    }

    public static double[] doub_vec(int n){
        return new double[n];
    }

    public static double[] broyden(fi_xi vecFunc, final double[] x) throws ConvergenceException {
        return broyden(vecFunc, x, false, false);
    }

    public static double[] broyden(fi_xi vecFunc, final double[] x, boolean initialJacobianIdentity, boolean changeJacobianIdentity) throws ConvergenceException {
        // Given an initial guess x[0..n-1] for a root in n dimensions, find the root by Broyden’s
        // method embedded in a globally convergent strategy. The vector of functions to be zeroed
        // called fvec[0..n-1] in the routine below, is returned by the user-supplied function or functor
        // vecfunc. On a successful run the x[0..n-1] which zeros fvec is returned.
        // ConvergenceException is thrown if the routine has converged to a local minimum of the function
        // fmin or if Broyden’s method can make no further progress. In this case try restarting from a
        // different initial guess.
        int MAXITS=250;
        double EPSILON = Math.ulp(1.0);
        double TOLF=Constants.tol, TOLX=EPSILON, STPMX=100.0, TOLMIN=Math.pow(Constants.tol,3/2);
        // Here MAXITS is the maximum number of iterations; EPS is the machine precision; TOLF
        // is the convergence criterion on function values; TOLX is the convergence criterion on delta_x;
        // STPMX is the scaled maximum step length allowed in line searches; and TOLMIN is used to
        // decide whether spurious convergence to a minimum of fmin has occurred.

        boolean restrt, skip, identityJacobian;
        int i, its, j, n = x.length;
        double den, fold, stpmax, sum, temp, test;
        double f = 0.0;

        double[] fvcold = doub_vec(n), g = doub_vec(n), p = doub_vec(n), s = doub_vec(n), t = doub_vec(n), w = doub_vec(n), xold = doub_vec(n);
        QRdcmp qr = null;
        fb fdjac = new fb(vecFunc);

        fmin fmin = new fmin(vecFunc);

        double[] fvec;
        try {
            f = fmin.calc(x);
        } catch (ConvergenceException e){
            e.setErrorLocation("Broyden (" + e.getErrorLocation() + ")");
            e.setIndex(0);
            throw e;
        }
        fvec=fmin.getFvec();

        test = 0.0;
        for (i = 0; i < n; i++) {
            // Test for initial guess being a root. Use more
            // stringent test than simply TOLF.
            if (Math.abs(fvec[i]) > test)
                test = Math.abs(fvec[i]);
        }
        if (test < 0.01 * TOLF) {
            // hurray! initial guess is a root!
            return x;
        }
        sum = 0.0;
        for (i = 0; i < n; i++)
            sum += x[i]*x[i]; // Calculate stpmax for line searches.
        stpmax = STPMX * Math.max(Math.sqrt(sum), (double)n);
        restrt = true; // Ensure initial Jacobian gets computed.
        identityJacobian = initialJacobianIdentity;    // first time use the identity matrix for the Jacobian

        for (its = 1; its <= MAXITS; its++) { // Start of iteration loop.
            //System.out.println(its +" " + restrt);
            if (restrt) { // Initialize or reinitialize Jacobian and QR decompose it.
                qr = new QRdcmp(fdjac.func(x, identityJacobian));
                if (changeJacobianIdentity) identityJacobian=!identityJacobian; // next time use actual Jacobian
                if (qr.sing()) {
                    throw new ConvergenceException.Builder("singular Jacobian in broydn. ", "Broyden")
                                                    .setIndex(its)
                                                    .setConvergenceDistance(test)
                                                    .build();
                }
            } else { // Carry out Broyden update.
                for (i = 0; i < n; i++)
                    s[i] = x[i] - xold[i]; // s = delta_x.
                for (i = 0; i < n; i++) { // t = R(dot)s
                    for (sum = 0.0, j = i; j < n; j++)
                        sum += qr.r()[i][j] * s[j];
                    t[i] = sum;
                }
                skip = true;
                for (i = 0; i < n; i++) { // w = delta_F - B (dot) s.
                    for (sum = 0.0, j = 0; j < n; j++)
                        sum += qr.qt()[j][i] * t[j];
                    w[i] = fvec[i] - fvcold[i] - sum;
                    if (Math.abs(w[i]) >= EPSILON * (Math.abs(fvec[i]) + Math.abs(fvcold[i])))
                        skip = false;
                        // Don’t update with noisy components of w.
                    else
                        w[i] = 0.0;
                }
                if (!skip) {
                    qr.qtmult(w, t); // t = QT (dot) w.
                    for (den = 0.0, i = 0; i < n; i++)
                        den += s[i]*s[i];
                    for (i = 0; i < n; i++)
                        s[i] /= den; // Store s/(s(dot)s) in s.
                    qr.update(t, s); // Update R and QT .
                    if (qr.sing())
                        throw new ConvergenceException.Builder("singular update in broydn. ", "Broyden")
                                .setIndex(its)
                                .setConvergenceDistance(test)
                                .build();
                }
            }
            qr.qtmult(fvec, p);
            for (i = 0; i < n; i++)
                // Right-hand side for linear equations is -QT (dot) F.
                p[i] = -p[i];
            for (i = n - 1; i >= 0; i--) { // Compute grad f ~= (Q(dot)R)^T (dot) F for the line search.
                for (sum = 0.0, j = 0; j <= i; j++)
                    sum -= qr.r()[j][i] * p[j];
                g[i] = sum;
            }
            for (i = 0; i < n; i++) { // Store x and F.
                xold[i] = x[i];
                fvcold[i] = fvec[i];
            }
            fold = f; // Store f.
            try {
                qr.rsolve(p, p); // Solve linear equations.
            } catch (Exception e){
                throw new ConvergenceException.Builder("Something went wrong with rsolve in QRdcmp. ", "Broyden")
                        .setIndex(its)
                        .setConvergenceDistance(test)
                        .setCause(e)
                        .build();
            }

            Double lnsrchRet;
            try {
                lnsrchRet = lnsrch(xold, fold, g, p, x, stpmax, EPSILON, fmin);
            } catch (ConvergenceException e){
                e.setErrorLocation("Broyden (" + e.getErrorLocation() + ")");
                e.setIndex(its);
                throw e;
            }
            fvec=fmin.getFvec();

            // lnsrch returns new x and f.

            test = 0.0; // Test for convergence on function values.
            for (i = 0; i < n; i++)
                if (Math.abs(fvec[i]) > test)
                    test = Math.abs(fvec[i]);
            if (test < TOLF) {
                qr = null;
                return x;
            }
            if (lnsrchRet==null) { // True if line search failed to find a new x.
                if (restrt) { // Failure; already tried reinitializing the Jacobian.
                    qr = null;
                    throw new ConvergenceException.Builder("Broyden failed to find solution even after reinitialization " +
                            "of the Jacobian. Perhaps try a different initial guess. ", "Broyden")
                            .setIndex(its)
                            .setConvergenceDistance(test)
                            .build();
                } else {
                    test = 0.0; // Check for gradient of f zero, i.e., spurious convergence.
                    den = Math.max(f, 0.5 * n);
                    for (i = 0; i < n; i++) {
                        temp = Math.abs(g[i]) * Math.max(Math.abs(x[i]), 1.0) / den;
                        if (temp > test)
                            test = temp;
                    }
                    if (test < TOLMIN) {
                        // converged.
                        qr = null;
                        return x;
                    } else
                        restrt = true; // Try reinitializing the Jacobian.
                }
            } else { // Successful step; will use Broyden update for next
                f=lnsrchRet.doubleValue();
                restrt = false; // step.
                test = 0.0; // Test for convergence on delta_x.
                for (i = 0; i < n; i++) {
                    temp = (Math.abs(x[i] - xold[i])) / Math.max(Math.abs(x[i]), 1.0);
                    if (temp > test)
                        test = temp;
                }
                if (test < TOLX) {
                    qr = null;
                    return x;
                }
            }
        }
        throw new ConvergenceException.Builder("MAXITS exceeded in broyden. ", "Broyden")
                                        .setIndex(its)
                                        .setConvergenceDistance(test)
                                        .build();

    }

    public static double[] newt(fi_xi vecfunc, final double[] x)
            throws ConvergenceException {
        // Given an initial guess x[0..n-1] for a root in n dimensions, find the
        // root by a globally convergent Newton’s method. The vector of
        // functions
        // to be zeroed, called fvec[0..n-1] in the routine below, is returned
        // by the user-supplied function or functor vecfunc (see text). The
        // output quantity check is false on a normal return and true if the
        // routine has converged to a local minimum of the function fmin defined
        // below. In this case try restarting from a different initial guess.
        final double EPSILON=Math.ulp(1.0);
        final int MAXITS = 200;
        final double TOLF = Constants.tol, TOLMIN = Math.pow(Constants.tol,3/2), STPMX = 100.0;
        final double TOLX = EPSILON; // numeric_limits<double>::epsilon();
        // Here MAXITS is the maximum number of iterations; TOLF sets the
        // convergence criterion on function values; TOLMIN sets the criterion
        // for deciding whether spurious convergence to a minimum of fmin has
        // occurred; STPMX is the scaled maximum step length allowed in line
        // searches; and TOLX is the convergence criterion on ix.
        int i, j, its, n = x.length;
        double den, fold, stpmax, sum, temp, test;
        double f = 0.0;
        double[] g = doub_vec(n), p = doub_vec(n), xold = doub_vec(n);
        double[][] fjac = new double[n][n];
        fb fdjac = new fb(vecfunc);
        fmin fmin = new fmin(vecfunc); // Set up fmin object.

        double[] fvec;
        try {
            f = fmin.calc(x);
        } catch (ConvergenceException e){
            e.setIndex(0);
            e.setErrorLocation("Newton (" + e.getErrorLocation() + ")");
            throw e;
        }
        fvec=fmin.getFvec();

        test = 0.0; // Test for initial guess being a root. Use
        for (i = 0; i < n; i++) {
            // more stringent test than simply TOLF.
            if (Math.abs(fvec[i]) > test) test = Math.abs(fvec[i]);
        }
        if (test < 0.01 * TOLF) {
            // hurray! initial guess is a root!
            return x;
        }

        sum = 0.0;
        for (i = 0; i < n; i++)
            sum += x[i]*x[i]; // Calculate stpmax for line searches.
        stpmax = STPMX * Math.max(Math.sqrt(sum), (double)n);

        for (its = 0; its < MAXITS; its++) { // Start of iteration loop.
            fjac = fdjac.func(x, false);

            for (i = 0; i < n; i++) { // Compute grad(f) for the line search.
                sum = 0.0;
                for (j = 0; j < n; j++)
                    sum += fjac[j][i] * fvec[j];
                g[i] = sum;
            }
            for (i = 0; i < n; i++)
                xold[i] = x[i]; // Store x,
            fold = f; // and f .

            for (i = 0; i < n; i++)
                p[i] = -fvec[i]; // Right-hand side for linear equations.
            try {
                LUdcmp alu = new LUdcmp(fjac); // Solve linear equations by LU
                // decomposition
                alu.solve(p, p);
            } catch (SingularMatrixException e){
                throw new ConvergenceException.Builder("Singular Jacobian matrix given for LU decomposition. ", "Newton")
                        .setIndex(its)
                        .setCause(e)
                        .build();
            } catch (MatrixDimensionMismatchException e){
                throw new ConvergenceException.Builder("Problem with matrix dimensions in LU decomposition. ", "Newton")
                        .setIndex(its)
                        .setCause(e)
                        .build();
            }
            Double lnsrchRet;
            try {
                lnsrchRet = lnsrch(xold, fold, g, p, x, stpmax, EPSILON, fmin);
            }catch (ConvergenceException e){
                e.setIndex(its);
                e.setErrorLocation("Newton (" + e.getErrorLocation() + ")");
                throw e;
            }
            fvec=fmin.getFvec();

            // lnsrch returns new x and f.

            test = 0.0; // Test for convergence on function values.
            for (i = 0; i < n; i++)
                if (Math.abs(fvec[i]) > test)
                    test = Math.abs(fvec[i]);
            if (test < TOLF) {
                return x;
            }
            if (lnsrchRet == null) { // Check for gradient of f zero, i.e., spurious convergence.
                test = 0.0;
                den = Math.max(f, 0.5 * n);
                for (i = 0; i < n; i++) {
                    temp = Math.abs(g[i]) * Math.max(Math.abs(x[i]), 1.0) / den;
                    if (temp > test)
                        test = temp;
                }
                if (test < TOLMIN){
                    // converged.
                    return x;
                }

            }
            test = 0.0; // Test for convergence on ix.
            for (i = 0; i < n; i++) {
                temp = (Math.abs(x[i] - xold[i])) / Math.max(Math.abs(x[i]), 1.0);
                if (temp > test)
                    test = temp;
            }
            if (test < TOLX)
                return x;
        }
        throw new ConvergenceException.Builder("MAXITS exceeded in newt. ", "Newton")
                                        .setIndex(its)
                                        .setConvergenceDistance(test)
                                        .build();
    }



    public static Double lnsrch(final double[] xold, final double fold,
                                final double[] g, final double[] p, final double[] x, final double stpmax,
                                final double EPSILON, final f_xi func)  throws ConvergenceException{
        // Given an n-dimensional point xold[0..n-1], the value of the function
        // and gradient there, fold and g[0..n-1], and a direction p[0..n-1],
        // finds a new point x[0..n-1] along the direction p from xold where the
        // function or functor func has decreased “sufficiently.” The new
        // function value is returned. stpmax is an input quantity that limits the
        // length of the steps so that you do not try to evaluate the function
        // in regions where it is undefined or subject to overflow. p is usually
        // the Newton direction.
        // If a proper step was taken then the value of f at the new x[0..n-1] is
        // returned. The other option is that the expected step is smaller than EPSILON,
        // and then we return null. In this case the calling function should check for
        // convergence (delta_x convergence).
        final double ALF = 1.0e-4; final double TOLX = EPSILON;
        // ALF ensures sufficient decrease in function value; TOLX is the
        // convergence criterion on delta_x.
        double a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
        double f; // value of the function func
        double rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
        int i, n = xold.length;

        for (i = 0; i < n; i++)
            sum += p[i] * p[i];
        sum = Math.sqrt(sum);
        if (sum > stpmax)
            for (i = 0; i < n; i++)
                p[i] *= stpmax / sum; // Scale if attempted step is too big.
        for (i = 0; i < n; i++)
            slope += g[i] * p[i];
        if (slope >= 0.0) {
            //printMean(x);
            throw new RuntimeException("Roundoff problem in lnsrch. slope="+slope);
        }
        test = 0.0; // Compute lambda_min.
        for (i = 0; i < n; i++) {
            temp = Math.abs(p[i]) / Math.max(Math.abs(xold[i]), 1.0);
            if (temp > test)
                test = temp;
        }
        alamin = TOLX / test;
        alam = 1.0; // Always try full Newton step first.
        for (int its=0;its<200;its++) { // Start of iteration loop.
            for (i = 0; i < n; i++)
                x[i] = xold[i] + alam * p[i];

            f = func.calc(x);
            if (alam < alamin) { // It appears there is convergence on delta_x. For zero finding,
                // the calling program should verify the convergence.
                for (i = 0; i < n; i++)
                    x[i] = xold[i];
                return null;
            } else if (f <= fold + ALF * alam * slope)
                return f; // Sufficient function decrease.
            else { // Backtrack.
                if (alam == 1.0)
                    tmplam = -slope / (2.0 * (f - fold - slope)); // First time.
                else { // Subsequent backtracks.
                    rhs1 = f - fold - alam * slope;
                    rhs2 = f2 - fold - alam2 * slope;
                    a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                    b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
                    if (a == 0.0)
                        tmplam = -slope / (2.0 * b);
                    else {
                        disc = b * b - 3.0 * a * slope;
                        if (disc < 0.0)
                            tmplam = 0.5 * alam;
                        else if (b <= 0.0)
                            tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
                        else
                            tmplam = -slope / (b + Math.sqrt(disc));
                    }
                    if (tmplam > 0.5 * alam)
                        tmplam = 0.5 * alam; // lambda <= 0.5*lambda_1.
                }
            }
            alam2 = alam;
            f2 = f;
            alam = Math.max(tmplam, 0.1 * alam);// lambda >= 0.1*lambda_1
        } // Try again.
        throw new ConvergenceException.Builder("lnssrch could not find a new x[] after 200 iterations. ", "lnsrch")
                .build();
    }

}
