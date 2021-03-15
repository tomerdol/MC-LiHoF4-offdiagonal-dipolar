package simulation.montecarlo;

import org.ojalgo.array.Array1D;
import org.ojalgo.matrix.decomposition.Eigenvalue;
import org.ojalgo.matrix.store.ElementsSupplier;
import org.ojalgo.matrix.store.GenericDenseStore;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.scalar.ComplexNumber;
import org.ojalgo.structure.Access1D;
import static org.ojalgo.function.constant.ComplexMath.*;

/**
 * Builds and holds the crystal field Hamiltonian and given some (Bx,By,Bz) combination,
 * adds a Zeeman term and diagonalizes to get eigenenergies and magnetic moments.
 */
public class CrystalField {
    // constants
    private static final double hbar, g_L, u_B=0.6717;
    private static final double J;
    private static final int deg_J;
    private static final PhysicalStore.Factory<ComplexNumber, GenericDenseStore<ComplexNumber>> storeFactory = GenericDenseStore.COMPLEX;
    private static final GenericDenseStore<ComplexNumber> H_cf, jx, jy, jz, jplus, jminus, I_J;

    static{
        if (System.getProperty("system").equals("LiHoF4")) {
            hbar = 1;
            g_L = 1.25;
            J = 8;
            deg_J = (int) (J * 2) + 1;
            jx = storeFactory.makeZero(deg_J, deg_J);
            jy = storeFactory.makeZero(deg_J, deg_J);
            I_J = storeFactory.makeEye(deg_J, deg_J);
            jplus = initJplus();
            jminus = initJminus();
            jx.fillMatching(jplus, ADD, jminus);
            jx.operateOnAll(MULTIPLY, ComplexNumber.of(0.5, 0)).supplyTo(jx);
            jy.fillMatching(jplus, SUBTRACT, jminus);
            jy.operateOnAll(MULTIPLY, ComplexNumber.of(0, -0.5)).supplyTo(jy);
            jz = initJz();
            H_cf = initH_cfLiHoF4();
        } else if (System.getProperty("system").equals("Fe8")){
            hbar = 1;
            g_L = 2.0;
            J = 10;
            deg_J = (int) (J * 2) + 1;
            jx = storeFactory.makeZero(deg_J, deg_J);
            jy = storeFactory.makeZero(deg_J, deg_J);
            I_J = storeFactory.makeEye(deg_J, deg_J);
            jplus = initJplus();
            jminus = initJminus();
            jx.fillMatching(jplus, ADD, jminus);
            jx.operateOnAll(MULTIPLY, ComplexNumber.of(0.5, 0)).supplyTo(jx);
            jy.fillMatching(jplus, SUBTRACT, jminus);
            jy.operateOnAll(MULTIPLY, ComplexNumber.of(0, -0.5)).supplyTo(jy);
            jz = initJz();
            H_cf = initH_cfLiHoF4();
        } else {
            throw new RuntimeException("Cannot initialize crystal field parameters. Illegal system name given.");
        }
    }

    // initialize the eigenvalue decomposition object given (Bx,By,Bz)
    private static Eigenvalue<ComplexNumber> initDecomp(double Bx, double By, double Bz){
        // initialize 3 zeeman matrices
        // *********************************************************
        GenericDenseStore<ComplexNumber> BxJx = storeFactory.makeZero(deg_J,deg_J), ByJy = storeFactory.makeZero(deg_J,deg_J), BzJz = storeFactory.makeZero(deg_J,deg_J);
        GenericDenseStore<ComplexNumber> tmpStore;

        tmpStore = storeFactory.copy(jx);
        tmpStore.operateOnAll(MULTIPLY, ComplexNumber.of(Bx*u_B*g_L,0)).supplyTo(BxJx);

        tmpStore = storeFactory.copy(jy);
        tmpStore.operateOnAll(MULTIPLY, ComplexNumber.of(By*u_B*g_L,0)).supplyTo(ByJy);

        tmpStore= storeFactory.copy(jz);
        tmpStore.operateOnAll(MULTIPLY, ComplexNumber.of(Bz*u_B*g_L,0)).supplyTo(BzJz);
        // *********************************************************

        ElementsSupplier<ComplexNumber> tmpZeeman = BxJx.operateOnMatching(ADD, ByJy);
        tmpZeeman = tmpZeeman.operateOnMatching(ADD, BzJz);
        ElementsSupplier<ComplexNumber> tmpFullH = tmpZeeman.operateOnMatching(SUBTRACT, H_cf);
        tmpFullH = tmpFullH.operateOnAll(MULTIPLY, ComplexNumber.NEG);
        Eigenvalue<ComplexNumber> decomp = Eigenvalue.COMPLEX.make(H_cf,true);
        decomp.decompose(tmpFullH);
        return decomp;
    }

    // returns a size 2 int array containing the indices of the lowest and second lowest eigenenergies (at places 0 and 1 respectively)
    private static int[] findTwoLowestEnergies(Array1D<ComplexNumber> eigenvalues ){
        // look for smallest eigenvalues (and then use corresponding eigenvectors)
        int smallestIndex = 0, secondSmallestIndex=1;

        if (eigenvalues.get(smallestIndex).getReal()>eigenvalues.get(secondSmallestIndex).getReal()){
            // reverse
            smallestIndex = 1;
            secondSmallestIndex = 0;
        }

        for (int i=2;i<eigenvalues.length;i++){
            if (eigenvalues.get(i).getReal() < eigenvalues.get(smallestIndex).getReal()){
                secondSmallestIndex=smallestIndex;
                smallestIndex=i;
            } else if (eigenvalues.get(i).getReal() < eigenvalues.get(secondSmallestIndex).getReal() && eigenvalues.get(i).getReal() > eigenvalues.get(smallestIndex).getReal()) {
                secondSmallestIndex = i;
            }
        }

        return new int[]{smallestIndex,secondSmallestIndex};
    }

    /**
     * Calculate the magnetic moment of a spin-up Holmium ion under crystal field potential and magnetic field (Bx,By,Bz).
     * @param Bx Bx value
     * @param By By value
     * @param Bz Bz value
     * @return Magnetic moment of the Ho ion under the given field and crystal field.
     */
    public static double getMagneticMoment(double Bx, double By, double Bz){
        // default call to getMagneticMoment should be for spin UP.
        return getMagneticMoment(Bx,By,Bz,1);
    }

    /**
     * Calculate the magnetic moment of a Holmium ion under crystal field potential and magnetic field (Bx,By,Bz).
     * @param Bx Bx value
     * @param By By value
     * @param Bz Bz value
     * @param spin spin direction (1 or -1)
     * @return Magnetic moment of the Ho ion under the given field and crystal field.
     */
    public static double getMagneticMoment(double Bx, double By, double Bz, int spin){
        Eigenvalue<ComplexNumber> decomp = initDecomp(Bx,By,Bz);
        MatrixStore<ComplexNumber> V = decomp.getV();
        MatrixStore<ComplexNumber> VT = V.conjugate();
        Access1D<ComplexNumber> magneticMoments = VT.multiply(jz).multiply(V).sliceDiagonal();
        Array1D<ComplexNumber> eigenvalues = decomp.getEigenvalues();

        int[] lowIndices = findTwoLowestEnergies(eigenvalues);

        if (spin==1){
            return Math.max(magneticMoments.get(lowIndices[0]).getReal(), magneticMoments.get(lowIndices[1]).getReal());
        } else if (spin==-1){
            return Math.min(magneticMoments.get(lowIndices[0]).getReal(), magneticMoments.get(lowIndices[1]).getReal());
        }
        throw new RuntimeException("received spin is neither 1 nor -1");
    }

    /**
     * Calculate the energy of a spin-up Holmium ion under crystal field potential and magnetic field (Bx,By,Bz).
     * @param Bx Bx value
     * @param By By value
     * @param Bz Bz value
     * @return Energy of the Ho ion under the given field and crystal field.
     */
    public static double getEnergy(double Bx, double By, double Bz){
        return getEnergy(Bx,By,Bz,1);
    }


    /**
     * Calculate the energy of a Holmium ion under crystal field potential and magnetic field (Bx,By,Bz).
     * @param Bx Bx value
     * @param By By value
     * @param Bz Bz value
     * @param spin spin direction (1 or -1)
     * @return Energy of the Ho ion under the given field and crystal field.
     */
    public static double getEnergy(double Bx, double By, double Bz, int spin){
        Eigenvalue<ComplexNumber> decomp = initDecomp(Bx,By,Bz);
        MatrixStore<ComplexNumber> V = decomp.getV();
        MatrixStore<ComplexNumber> VT = V.conjugate();
        Access1D<ComplexNumber> magneticMoments = VT.multiply(jz).multiply(V).sliceDiagonal();
        Array1D<ComplexNumber> eigenvalues = decomp.getEigenvalues();
        int[] lowIndices = findTwoLowestEnergies(eigenvalues);

        if (spin==1){
            return magneticMoments.get(lowIndices[0]).getReal() > magneticMoments.get(lowIndices[1]).getReal() ? eigenvalues.get(lowIndices[0]).getReal() : eigenvalues.get(lowIndices[1]).getReal();
        } else if (spin==-1){
            return magneticMoments.get(lowIndices[0]).getReal() < magneticMoments.get(lowIndices[1]).getReal() ? eigenvalues.get(lowIndices[0]).getReal() : eigenvalues.get(lowIndices[1]).getReal();
        }
        throw new RuntimeException("received spin is neither 1 nor -1");
    }


    public static void main(String[] args){
        // tests
        System.out.println(getMagneticMoment(0.03,1.6,0.05,1));
        System.out.println(getEnergy(0.03,1.6,0.05,1));
    }

    // matrix power (no special technique is used for this so it should not be used for very high powers)
    // returns a new instance of GenericDenseStore<ComplexNumber>.
    private static GenericDenseStore<ComplexNumber> power(GenericDenseStore<ComplexNumber> A, int n){
        GenericDenseStore<ComplexNumber> tmpA = A.copy();
        if (n==0){
            return storeFactory.makeEye(A.countRows(),A.countColumns());
        } else{
            for (int i=1;i<n;i++){
                tmpA.multiply(A).supplyTo(tmpA);
            }
        }
        return tmpA;
    }


    private static GenericDenseStore<ComplexNumber> initH_cfLiHoF4(){
        // crystal field parameters:
        double B02 = -0.696;
        double B04 = 4.06e-3;
        double B06 =  4.64e-6;
        double B44C = 0.0418;
        double B46C = 8.12e-4;
        double B46S = 1.137e-4;

        final GenericDenseStore<ComplexNumber> O02, O04, O44C, O06, O46C1, O46S1, O46C, O46S;

        //PhysicalStore.Factory<ComplexNumber, GenericDenseStore<ComplexNumber>> storeFactory = GenericDenseStore.COMPLEX;

        O02 = storeFactory.makeZero(deg_J,deg_J);
        O04 = storeFactory.makeZero(deg_J,deg_J);
        O44C = storeFactory.makeZero(deg_J,deg_J);
        O06 = storeFactory.makeZero(deg_J,deg_J);
        O46C1 = storeFactory.makeZero(deg_J,deg_J);
        O46S1 = storeFactory.makeZero(deg_J,deg_J);
        O46C = storeFactory.makeZero(deg_J,deg_J);
        O46S = storeFactory.makeZero(deg_J,deg_J);

        // O02
        O02.fillByMultiplying(jz,jz);   // jz^2
        O02.modifyAll(MULTIPLY.second(ComplexNumber.of(3,0)));  // jz^2 * 3
        O02.operateOnMatching(SUBTRACT,I_J.multiply(J*(J+1))).supplyTo(O02);  // jz^2 * 3 - I * J(J+1)

        // O04
        jz.supplyTo(O04);   //jz
        O04.operateOnAll(POWER,4).supplyTo(O04);  // jz^4 (only because jz is diagonal!)
        O04.operateOnAll(MULTIPLY,ComplexNumber.of(35,0)).supplyTo(O04);  // jz^4 * 35
        O04.operateOnMatching(SUBTRACT, jz.multiply(jz).multiply(ComplexNumber.of(30 * J * (J+1),0))).supplyTo(O04);      //jz^4 * 35 - 30J(J+1)*jz^2
        O04.operateOnMatching(ADD, jz.multiply(jz).multiply(ComplexNumber.of(25,0))).supplyTo(O04);   //jz^4 * 35 - 30J(J+1)*jz^2 + 25*jz^2
        O04.operateOnMatching(SUBTRACT, I_J.multiply(ComplexNumber.of(6 * J * (J+1),0))).supplyTo(O04);   //jz^4 * 35 - 30J(J+1)*jz^2 + 25*jz^2 - 6J(J+1) * I
        O04.operateOnMatching(ADD, I_J.multiply(ComplexNumber.of(3 * J*J * (J+1)*(J+1),0))).supplyTo(O04);


        // O44C
        O44C.fillByMultiplying(jplus.multiply(jplus),jplus.multiply(jplus));    // jplus^4
        O44C.operateOnMatching(ADD,power(jminus,4)).supplyTo(O44C);    //jplus^4 + jminus^4
        O44C.operateOnAll(MULTIPLY, ComplexNumber.of(0.5,0)).supplyTo(O44C);   //(jplus^4 + jminus^4) * 0.5


        // O06
        jz.supplyTo(O06);   //jz
        O06.operateOnAll(POWER,6).supplyTo(O06);  // jz^6 (only because jz is diagonal!)
        O06.operateOnAll(MULTIPLY,ComplexNumber.of(231,0)).supplyTo(O06);  // jz^6 * 231
        O06.operateOnMatching(SUBTRACT, power(jz,4).multiply(ComplexNumber.of(315*J*(J+1),0))).supplyTo(O06); // jz^6 * 231 - jz^4 * 315J(J+1)
        O06.operateOnMatching(ADD, power(jz,4).multiply(ComplexNumber.of(735,0))).supplyTo(O06);  // jz^6 * 231 - jz^4 * 315J(J+1) + 735 * jz^4
        O06.operateOnMatching(ADD, jz.multiply(jz).multiply(ComplexNumber.of(105 * J*J * (J+1)*(J+1),0))).supplyTo(O06);  // jz^6 * 231 - jz^4 * 315J(J+1) + 735 * jz^4 + 105 * J^2 * (J+1)^2 * jz^2
        O06.operateOnMatching(SUBTRACT, jz.multiply(jz).multiply(ComplexNumber.of(525*J*(J+1),0))).supplyTo(O06); // jz^6 * 231 - jz^4 * 315J(J+1) + 735 * jz^4 + 105 * J^2 * (J+1)^2 * jz^2 - 525J(J+1) * jz^2
        O06.operateOnMatching(ADD, jz.multiply(jz).multiply(ComplexNumber.of(294,0))).supplyTo(O06);  // jz^6 * 231 - jz^4 * 315J(J+1) + 735 * jz^4 + 105 * J^2 * (J+1)^2 * jz^2 - 525J(J+1) * jz^2 + jz^2 * 294
        O06.operateOnMatching(SUBTRACT, I_J.multiply(ComplexNumber.of(5 * J*J*J * (J+1)*(J+1)*(J+1),0))).supplyTo(O06);   // jz^6 * 231 - jz^4 * 315J(J+1) + 735 * jz^4 + 105 * J^2 * (J+1)^2 * jz^2 - 525J(J+1) * jz^2 + jz^2 * 294 - 5*J^3*(J+1)^3 * I
        O06.operateOnMatching(ADD, I_J.multiply(ComplexNumber.of(40 * J*J *(J+1)*(J+1),0))).supplyTo(O06);    // jz^6 * 231 - jz^4 * 315J(J+1) + 735 * jz^4 + 105 * J^2 * (J+1)^2 * jz^2 - 525J(J+1) * jz^2 + jz^2 * 294 - 5*J^3*(J+1)^3 * I + 40*J^2*(J+1)^2 * I
        O06.operateOnMatching(SUBTRACT, I_J.multiply(ComplexNumber.of(60*J*(J+1),0))).supplyTo(O06);


        // O46C1
        O46C1.fillByMultiplying((power(jplus,4).add(power(jminus,4))).multiply(ComplexNumber.of(0.25,0)),(jz.multiply(jz).multiply(ComplexNumber.of(11,0))).subtract(I_J.multiply(J*(J+1))).subtract(I_J.multiply(38)));

        // O46S1
        O46S1.fillByMultiplying((power(jplus,4).subtract(power(jminus,4))).multiply(ComplexNumber.of(0,-0.25)),jz.multiply(jz).multiply(ComplexNumber.of(11,0)).subtract(I_J.multiply(J*(J+1))).subtract(I_J.multiply(38)));

        // O46C
        O46C.fillMatching(O46C1, ADD, O46C1.conjugate());

        // O46S
        O46S.fillMatching(O46S1, ADD, O46S1.conjugate());

        GenericDenseStore<ComplexNumber> H_cf = storeFactory.makeZero(deg_J,deg_J);
        H_cf.fillMatching(O02.multiply(B02), ADD, O04.multiply(B04));
        H_cf.operateOnMatching(ADD, O06.multiply(B06)).supplyTo(H_cf);
        H_cf.operateOnMatching(ADD, O44C.multiply(B44C)).supplyTo(H_cf);
        H_cf.operateOnMatching(ADD, O46C.multiply(B46C)).supplyTo(H_cf);
        H_cf.operateOnMatching(ADD, O46S.multiply(B46S)).supplyTo(H_cf);

        return H_cf;
    }



    private static GenericDenseStore<ComplexNumber> initH_cfFe8(){
        // crystal field parameters:
        double D=0.294;
        double E=0.046;

        final GenericDenseStore<ComplexNumber> Sz2, Sx2, Sy2;

        Sz2 = storeFactory.makeZero(deg_J,deg_J);
        Sx2 = storeFactory.makeZero(deg_J,deg_J);
        Sy2 = storeFactory.makeZero(deg_J,deg_J);

        Sz2.fillByMultiplying(jz,jz);   // jz^2
        Sx2.fillByMultiplying(jx,jx);   // jx^2
        Sy2.fillByMultiplying(jy,jy);   // jy^2

        GenericDenseStore<ComplexNumber> H_cf = storeFactory.makeZero(deg_J,deg_J);
        H_cf.fillMatching(Sz2.multiply(-D), ADD, Sx2.multiply(E));
        H_cf.operateOnMatching(ADD, Sy2.multiply(-E)).supplyTo(H_cf);

        return H_cf;
    }


    private static GenericDenseStore<ComplexNumber> initJplus(){
        ComplexNumber[][] jplusArray = new ComplexNumber[deg_J][deg_J];
        double m=-J;
        for (int i=0;i<jplusArray.length;i++){
            for (int j=0;j<jplusArray[i].length;j++){
                if (i==j-1) {
                    jplusArray[i][j] = ComplexNumber.of(hbar*Math.sqrt(J*(J+1)-m*(m+1)),0);

                    m++;
                }else{
                    jplusArray[i][j] = ComplexNumber.ZERO;
                }
            }
        }
        PhysicalStore.Factory<ComplexNumber, GenericDenseStore<ComplexNumber>> storeFactory = GenericDenseStore.COMPLEX;
        GenericDenseStore<ComplexNumber> mat = storeFactory.rows(jplusArray);
        return mat;
    }

    private static GenericDenseStore<ComplexNumber> initJminus(){
        ComplexNumber[][] jminusArray = new ComplexNumber[deg_J][deg_J];
        double m=-J+1;
        for (int i=0;i<jminusArray.length;i++){
            for (int j=0;j<jminusArray[i].length;j++){
                if (i==j+1) {
                    jminusArray[i][j] = ComplexNumber.of(hbar*Math.sqrt(J*(J+1)-m*(m-1)),0);
                    m++;
                }else{
                    jminusArray[i][j] = ComplexNumber.ZERO;
                }
            }
        }
        PhysicalStore.Factory<ComplexNumber, GenericDenseStore<ComplexNumber>> storeFactory = GenericDenseStore.COMPLEX;
        GenericDenseStore<ComplexNumber> mat = storeFactory.rows(jminusArray);
        return mat;
    }

    private static GenericDenseStore<ComplexNumber> initJz(){
        ComplexNumber[][] jzArray = new ComplexNumber[deg_J][deg_J];
        double m=-J;
        for (int i=0;i<jzArray.length;i++){
            for (int j=0;j<jzArray[i].length;j++) {
                if (i==j) {
                    jzArray[i][j] = ComplexNumber.of(hbar * m, 0);
                    m++;
                }else{
                    jzArray[i][j] = ComplexNumber.ZERO;
                }
            }
        }

        PhysicalStore.Factory<ComplexNumber, GenericDenseStore<ComplexNumber>> storeFactory = GenericDenseStore.COMPLEX;
        GenericDenseStore<ComplexNumber> mat = storeFactory.rows(jzArray);
        return mat;

    }
}
