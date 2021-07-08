package simulation.mmsolve;

import simulation.montecarlo.*;

/**
 * Jacobian function
 */
class fb extends fij_xi {
    /** The function for which a Jacobian is to be calculated. */
    final fi_xi f;

    /**
     * Constructs a Jacobian function (a function that receives a vector and returns a matrix).
     * @param fi - the function for which to calculate the Jacobian
     */
    public fb(fi_xi fi)
    {f=fi;}


    /**
     * Evaluates the Jacobian of f at the given point x[0..n-1].
     * @param x - vector specifying the point at which to evaluate the Jacobian
     * @param identity - if <code>true</code> the returned Jacobian is just the identity matrix
     * @return the Jacobian Matrix of the function at the point x[0..n-1]
     * @throws ConvergenceException only if {@code identity} is false
     * @throws IndexOutOfBoundsException if the given point is outside the bounds of the {@code FieldTable} (derivative via exact diagonalization is not supported).
     */
    public double[][] func(double x[], boolean identity) throws ConvergenceException, IndexOutOfBoundsException
    {
        if (!identity) return func(x);

        int n=x.length;
        double a[][]=new double[n][n];
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if (i==j) a[i][j]=1;
                else a[i][j]=0;
            }
        }
        return a;
    }

    /**
     * Evaluates the Jacobian of f at the given point x[0..n-1].
     * @param x - vector specifying the point at which to evaluate the Jacobian
     * @return the Jacobian Matrix of the function at the point x[0..n-1]
     * @throws ConvergenceException if the number of times an exact diagonalization was performed exceeds 20.
     * @throws IndexOutOfBoundsException if the given point is outside the bounds of the {@code FieldTable} (derivative via exact diagonalization is not supported).
     */
    public double[][] func(double x[]) throws ConvergenceException, IndexOutOfBoundsException
    {
        if (!(f instanceof func)) {
            System.err.println("Cannot calculate Jacobian since f is not func.");
            System.exit(1);
        }

        boolean homotopy=false;
        double frac=1.0;
        int flipSpin=-1;
        if (f instanceof funcForHomotopy){
            homotopy=true;
            frac=((funcForHomotopy) f).getPerc();
            flipSpin = ((funcForHomotopy) f).getFlippedSpin();
        }

        // get the properties of f
        FieldTable momentTable = ((func)f).momentTable;
        singleSpin[] arr = ((func)f).arr;
        double extBx = ((func)f).extBx;
        double extBy = ((func)f).extBy;
        double[][][] E = ((func)f).intTable;
        boolean suppressInternalTransFields = ((func)f).suppressInternalTransFields;

        int n=x.length;
        double[][] ret=new double[n][n];
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++){
                ret[i][j]=0;
                if (i==j) ret[i][j]=1;

                if (!homotopy || i!=flipSpin) {
                    // the derivative of the required magnetic moment at site i is a function of the 3 component field, which is a function of the other spins.
                    // therefore, we use the chain rule here
                    ret[i][j] = ret[i][j] - (E[0][i][j] * momentTable.getDerivative(2, getField(x, i, E, extBx, extBy, suppressInternalTransFields), arr[i].getSpin(), arr[i].getPrevBIndices())
                            + E[1][i][j] * momentTable.getDerivative(1, getField(x, i, E, extBx, extBy, suppressInternalTransFields), arr[i].getSpin(), arr[i].getPrevBIndices())
                            + E[2][i][j] * momentTable.getDerivative(0, getField(x, i, E, extBx, extBy, suppressInternalTransFields), arr[i].getSpin(), arr[i].getPrevBIndices()));
                } else{
                    // this should work, but was not thoroughly tested and is currently unused.
                    // it is the same as above, with the flipped spin having a mixed magnetic moment.
                    System.err.println("problem! func for homotopy seems to be used");
                    ret[i][j] = ret[i][j] - (1-frac)*(E[0][i][j] * momentTable.getDerivative(2, getField(x, i, E, extBx, extBy, suppressInternalTransFields), arr[i].getSpin(), arr[i].getPrevBIndices())
                            + E[1][i][j] * momentTable.getDerivative(1, getField(x, i, E, extBx, extBy, suppressInternalTransFields), arr[i].getSpin(), arr[i].getPrevBIndices())
                            + E[2][i][j] * momentTable.getDerivative(0, getField(x, i, E, extBx, extBy, suppressInternalTransFields), arr[i].getSpin(), arr[i].getPrevBIndices()))
                            - (frac)*(E[0][i][j] * momentTable.getDerivative(2, getField(x, i, E, extBx, extBy, suppressInternalTransFields), -1*arr[i].getSpin(), arr[i].getPrevBIndices())
                            + E[1][i][j] * momentTable.getDerivative(1, getField(x, i, E, extBx, extBy, suppressInternalTransFields), -1*arr[i].getSpin(), arr[i].getPrevBIndices())
                            + E[2][i][j] * momentTable.getDerivative(0, getField(x, i, E, extBx, extBy, suppressInternalTransFields), -1*arr[i].getSpin(), arr[i].getPrevBIndices()));
                }
                // there is really no reason for this to happen here, but just in case.
                if (((func)f).numManualCalc >=20){
                    ConvergenceException e = new ConvergenceException.Builder("the function surpassed the permitted threshold (20) for manual calculations. ",
                            "Function evaluation")
                            .setNumManualCalc(((func)f).numManualCalc)
                            .build();
                    throw e;
                }
            }
        }
        return ret;
    }

    /**
     * Calculates the local field at site i
     * @param x - vector of the magnetic moments of all of the spins
     * @param i - the spin at which to calculate the field
     * @param int_config_Matrix - the interactions table
     * @param extBx - external Bx
     * @param extBy - external By
     * @param suppressInternalTransFields - are internal transverse fields suppressed
     * @return the magnetic field at spin <i>i</i> in the format {Bz,By,Bx}
     */
    private double[] getField(final double[] x, int i, double[][][] int_config_Matrix, double extBx, double extBy, boolean suppressInternalTransFields){
        double[] B = new double[]{extBx, extBy, 0};

        for (int dim=0; dim<B.length; dim++){
            if (!suppressInternalTransFields || dim==2) {
                for (int j = 0; j < int_config_Matrix[dim][i].length; j++) {
                    B[dim] += int_config_Matrix[dim][i][j] * x[j];
                }
            }
        }
        // here we define B[z,y,x] as opposed to B[x,y,z]:
        double temp=B[0];
        B[0]=B[2];
        B[2]=temp;

        return B;
    }

    /**
     * Evaluates the Jacobian at point x[0..n-1] by finite difference
     * @param x - vector that signifies the point for the function f at which to evaluate the Jacobian
     * @return the Jacobian matrix at point x[0...n-1]
     */
    public double[][] funcNumerical(double x[]) throws ConvergenceException
    {
        int n=x.length;
        double a[][]=new double[n][n];
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                a[i][j]=derivative(x,i,j);
            }
        }
        return a;
    }

    /**
     * Numerically calculates a derivative by finite difference - alternative to {@link fb#funcNumerical}
     * @param x - The point at which the derivative is to be calculated
     * @return the Jacobian evaluated at point x[0..n-1]
     */
    public double[][] altFunc(double x[]) throws ConvergenceException {
        double[] fvec = f.func(x);
        final double EPS = Math.ulp(1.0);
        int n=x.length;
        double[][] df = new double[n][n];
        double[] xh = new double[n];
        System.arraycopy(x,0,xh,0,xh.length);
        for (int j=0;j<n;j++){
            double temp = xh[j];
            double h=EPS*Math.abs(temp);
            if (h==0.0) h=EPS;
            xh[j]=temp+h;
            h=xh[j]-temp;
            double[] f1=f.func(xh);
            xh[j]=temp;
            for (int i=0;i<n;i++){
                df[i][j]=(f1[i]-fvec[i])/h;
            }
        }
        return df;
    }

    /**
     * Numerically calculates a derivative by finite difference
     * @param x - The point at which the derivative is to be calculated
     * @param denklem_ref - which vector element of the function f should be differentiated
     * @param x_ref - the variable with respect to which the differentiation should be performed
     * @return
     */
    public double derivative(double x[],int denklem_ref,int x_ref) throws ConvergenceException {
        // df(x)/dx[x_ref]  x_ref=0...n
        double h0=0.0256;
        int i,m;
        int n=7;
        double f1[];
        f1=new double[x.length];
        double f2[];
        f2=new double[x.length];
        double x1[];
        x1=new double[x.length];
        double x2[];
        x2=new double[x.length];
        for(i=0;i<x.length;i++)
        {
            x1[i]=x[i];
            x2[i]=x[i];
        }
        //derivative of a simple function
        double T[][];
        T=new double[n][n];
        double h[];
        h=new double[n];

        for(i=0;i<n;i++)
        {
            h[i]=0;
            for(int j=0;j<n;j++)
                T[i][j]=0;
        }
        h[0]=h0;
        double r=0.5;
        for( i=1;i<n;i++)
        {
            h[i]=h0*Math.pow(r,i);
        }

        for(i=0;i<n;i++)
        {
            x1[x_ref]+=h[i];
            x2[x_ref]-=h[i];
            f1=f.func(x1);
            f2=f.func(x2);
            T[i][0]=( f1[denklem_ref] - f2[denklem_ref])/(2.0*h[i]);
            x1[x_ref]=x[x_ref];
            x2[x_ref]=x[x_ref];
        }
        for(m=1;m<n;m++)
        {
            for(i=0;i<n-m;i++)
            {
                T[i][m]=(h[i]*h[i]*T[i+1][m-1]-h[i+m]*h[i+m]*T[i][m-1])/(h[i]*h[i]-h[i+m]*h[i+m]);
            }
        }
        double xx=T[0][n-1];
        return xx;
    }
}

