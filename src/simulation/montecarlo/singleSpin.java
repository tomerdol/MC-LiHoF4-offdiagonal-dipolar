package simulation.montecarlo;


import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.Serializable;

/**
 * A single spin object
 */
public class singleSpin implements Serializable{
	// for serialization, should not be changed
	private static final long serialVersionUID = 625318162648175439L;
	/** spin orientation: up=+1 and down=-1 */
	private int s;
	/** the magnetic moment of the ion (can change depending on the applied field) */
	private double spinSize;
	/** local applied magnetic field in the z direction (longitudinal) which is defined as the easy axis */
	private double localBz;
	/**
	 * the spin's ordinal number within its lattice.
	 * @see Lattice#generateIsingLattice(int, int, double, boolean[]) */
	private int n;
	/** local applied magnetic field in the x direction (transverse) */
	private double localBx;
	/** local applied magnetic field in the y direction (transverse) */
	private double localBy;
	/**
	 * first array contains the last indices of the local field (x,y,z) == (2,1,0) in the FieldTable
	 * second array is used to track whether consecutive calls (to find magneticMoment or energy) are correlated (separately for x,y,z)
	 * with 1 indicating correlation and 0 indicating otherwise.
	 */
	private int[][] prevBIndices;

	/**
	 * Constructs an empty {@code singleSpin} object with spin up and
	 * 0 in all other fields.
	 */
	public singleSpin(){
		this(1, 0, 0, 0, 0, 0);
	}

	/**
	 * Constructs a {@code singleSpin} object with the given parameters and without knowing the previous B indices
	 * @param s0 spin orientation (+1 for up, -1 for down)
	 * @param n0 spin index (ordinal number)
	 * @param localBz local magnetic field in the z direction
	 * @param localBy local magnetic field in the y direction
	 * @param localBx local magnetic field in the x direction
	 * @param spinSize0 magnetic moment of this spin (ion)
	 * @throws IllegalArgumentException if the given spin orientation is different than -1, 0 or +1.
	 */
	public singleSpin(int s0, int n0, double localBz, double localBy, double localBx, double spinSize0){
		n=n0;
		spinSize=spinSize0;
		if (s0==1 || s0==-1 || s0==0){
			s=s0;
		} else{
			throw new IllegalArgumentException("Invalid spin! spin should be 1 or -1 or 0");
		}
		
		// local fields will be determined only after the lattice is complete
		this.localBz = localBz;
		this.localBy = localBy;
		this.localBx = localBx;
		prevBIndices=new int[2][3];
	}

	/**
	 * Constructs a {@code singleSpin} object with the given parameters, including the prevBindices array
	 * @param s0 spin orientation (+1 for up, -1 for down)
	 * @param n0 spin index (ordinal number)
	 * @param localBz local magnetic field in the z direction
	 * @param localBy local magnetic field in the y direction
	 * @param localBx local magnetic field in the x direction
	 * @param spinSize0 magnetic moment of this spin (ion)
	 * @param prevBIndices0 array of the indices of the magnetic field (in the {@link FieldTable} object) from the last calculation
	 *                      they will be deep copies to this new object.
	 * @throws NullPointerException if the received {@code prevBIndices0} is {@code null}
	 * @throws IllegalArgumentException if the received {@code prevBIndices0} is not of the right size
	 */
	public singleSpin(int s0, int n0, double localBz, double localBy, double localBx, double spinSize0, int[][] prevBIndices0){
		// call standard constructor which also initializes int[][] prevBIndices0
		this(s0,n0,localBz,localBy,localBx,spinSize0);

		// verify the received arrays are legal
		if (prevBIndices0==null || prevBIndices0[0]==null || prevBIndices0[1]==null) {
			throw new NullPointerException("prevBIndices cannot be null. if it is not known it should not be passed and a new one will be initialized.");
		}
		if (prevBIndices0.length!=2 || prevBIndices0[0].length!=3 || prevBIndices0[1].length!=3){
			throw new IllegalArgumentException("received prevBIndices has wrong dimensions. should be 2x3");
		}

		//************************ copy the input array to the existing array in the new object ************************
		System.arraycopy(prevBIndices0[0], 0, prevBIndices[0], 0, prevBIndices0[0].length);
		System.arraycopy(prevBIndices0[1], 0, prevBIndices[1], 0, prevBIndices0[1].length);
		//**************************************************************************************************************
	}

	/**
	 * Constructs a {@code singleSpin} object with the given parameters w/o an applied field or spinSize
	 * @param s0 spin orientation (+1 for up, -1 for down)
	 * @param n0 spin index (ordinal number)
	 */
	public singleSpin(int s0, int n0){
		this(s0,n0,0,0,0,0);	// Bx,By,Bz and spinSize should be initialized later
	}

	/**
	 * Constructs a {@code singleSpin} object with the given parameters w/o an applied field
	 * @param s0 spin orientation (+1 for up, -1 for down)
	 * @param n0 spin index (ordinal number)
	 * @param spinSize the positive magnetic moment size (the sign will be determined according to the given {@code s0})
	 */
	public singleSpin(int s0, int n0, double spinSize){
		this(s0,n0,0,0,0,s0*(spinSize));	// Bx,By,Bz should be initialized later and spinSize is an initial guess
	}

	/**
	 * Copies a {@code singleSpin} object
	 * @param other Another {@code singleSpin} object to (deep) copy
	 */
	public singleSpin(singleSpin other){
		this(other.s, other.n, other.localBz, other.localBy, other.localBx, other.spinSize, other.prevBIndices);
	}
	
	// getters to return the required values
	public int getSpin(){
		return this.s;
	}
	public double getLocalBz(){
		return this.localBz;
	}

	public double getLocalBx(){
		return this.localBx;
	}

	public double getLocalBy(){
		return this.localBy;
	}
	public int getN(){
		return this.n;
	}
	public double getSpinSize() {
		return spinSize;
	}
	public int[][] getPrevBIndices() {
		return prevBIndices;
	}

	//setters

	/**
	 * Sets the spin orientation
	 * @param s0 spin orientation (+1 for up, -1 for down)
	 * @throws IllegalArgumentException if the given spin orientation is different than -1, 0 or +1.
	 */
	public void setSpin(int s0){
		this.s=s0;

	    if (s0==1 || s0==-1 || s0==0){
			this.s=s0;
		} else{
			throw new IllegalArgumentException("Invalid spin! spin should be 1 or -1");
		}

	}

	public void setLocalBz(double localBz){
		this.localBz=localBz;
	}
	public void setLocalBx(double localBx){
		this.localBx=localBx;
	}
	public void setLocalBy(double localBy){
		this.localBy=localBy;
	}
	public void setN(int n){
		this.n=n;
	}
	public void setSpinSize(double spinSize) {
		this.spinSize = spinSize;
	}

	//toString prints the spin & of the spin in the following form: 'Spin:-1 Position:(x,y,z)'

	/**
	 * Prints the {@code singleSpin} object
	 * @return a string containing the data in this object in the following form: "Spin:-1 Number: 25 Local Bz: -0.2 Local By: 1.2 Local Bx: 0.04 Spin Size: 5.51"
	 */
	public String toString(){
		return "Spin:"+s+" Number: "+this.n + " Local Bz:"+this.localBz+" Local By: "+this.localBy
				+ " Local Bx: " + this.localBx + " Spin Size: " + this.spinSize;
	}
	
	/**
	 * Calculates the distance between two spins with (naive) periodic boundary conditions
	 * @param another spin to which the distance will be calculated
	 * @param Lx number of unit cells in the x-direction
	 * @param Lz number of unit cells in the z-direction
	 * @return array of the form [ds,dx,dx,dz] where,
	 * 	ds = total Euclidean distance
	 * 	dx = Euclidean distance along the x direction
	 * 	dy = Euclidean distance along the y direction
	 * 	dz = Euclidean distance along the z direction
	 * @see <a href="https://en.wikipedia.org/wiki/Periodic_boundary_conditions#(A)_Restrict_particle_coordinates_to_the_simulation_box" target="_top">periodic boundary conditions wiki</a>
	 */
	public double[] distance(singleSpin another, int Lz, int Lx){

		double dx = this.getX(Lz,Lx) - another.getX(Lz,Lx);
		double dy = this.getY(Lz,Lx) - another.getY(Lz,Lx);
		double dz = this.getZ(Lz,Lx) - another.getZ(Lz,Lx);
		double xL = Lx*Constants.a;
		double zL = Lz*Constants.c;
		// implementation of the periodic boundary conditions for each dimension separately
		if (dx>xL*0.5)
			dx=dx-xL;
		if (dx<=-xL*0.5)
			dx=dx+xL;
		if (dy>xL*0.5)
			dy=dy-xL;
		if (dy<=-xL*0.5)
			dy=dy+xL;
		if (dz>zL*0.5)
			dz=dz-zL;
		if (dz<=-zL*0.5)
			dz=dz+zL;

		return new double[]{Math.sqrt(dx*dx+dy*dy+dz*dz),dx,dy,dz};
	}
	
	//flips the spin

	/**
	 * Flip this spin and indicate no correlation is expected in the next {@link FieldTable} calculation
	 */
	public void flipSpin(){
		s=s*(-1);	// flip the spin
		prevBIndices[1][0]=0;	// change the correlation of Bz to false
								// this is because a flipped spin acts like
								// the un-flipped spin with Bz -> -Bz
	}

	/**
	 * Calculates the (naive) dipolar interaction between two spins
	 * @param o spin with which the dipolar interaction will be calculated
	 * @param Lz number of unit cells in the z direction
	 * @param Lx number of unit cells in the x and y directions
	 * @return the dipolar interaction between this and the given {@code singleSpin}
	 * @deprecated because the diplor interaction is now calculated using the Ewlad method
	 * 				and read to an interactions table
	 */
	@Deprecated
	public double calcDipolarInteraction(singleSpin o, int Lz, int Lx){
		if (this==o)	return 0;
		double r=this.distance(o, Lz, Lx)[0];
		double rz=this.getZ(Lz,Lx)-o.getZ(Lz,Lx);
		double actual_height=Lz*Constants.c;
		double actual_length=Lx*Constants.a;
		// implementation of periodic boundary conditions
		if (rz>actual_height*0.5)
			rz=rz-actual_height;
		if (rz<=-actual_height*0.5)
			rz=rz+actual_height;
		
		return (this.s*o.s)*(r*r-3*rz*rz)/Math.pow(r, 5);
		
	}

	/**
	 * Gets the location of this spin as a {@link Vector3D} object.
	 * Uses the spins ordinal index number
	 * @param Lz number of unit cells in the z direction of the system
	 * @param Ly number of unit cells in the y direction of the system
	 * @return a 3D vector which indicates the location of this spin
	 */
	public Vector3D getLocation(int Lz, int Ly){
		int i = ((this.n/Constants.num_in_cell)/Lz)/Ly;
		int j = ((this.n/Constants.num_in_cell)/Lz)%Ly;
		int k = (this.n/Constants.num_in_cell)%Lz;
		int l = this.n%Constants.num_in_cell;

		int[] locCoordinates = new int[]{i, j, k};

		Vector3D loc = Vector3D.ZERO;
		for (int coor=0;coor<Constants.primitiveLatticeVectors.length;coor++){
			loc=loc.add(locCoordinates[coor],Constants.primitiveLatticeVectors[coor]);	// add primitive lattice vectors according to (i,j,k)
			loc=loc.add(Constants.basis[l][coor],Constants.primitiveLatticeVectors[coor]);	// add fraction of primitive lattice vector according to basis index l
		}
		return loc;
	}

	/**
	 * Gets the x component of the location of this spin
	 * Uses the spins ordinal index number
	 * @param Lz number of unit cells in the z direction of the system
	 * @param Ly number of unit cells in the y direction of the system
	 * @return the x component of the location of this spin
	 */
	public double getX(int Lz, int Ly){
		return getLocation(Lz, Ly).getX();
	}

	/**
	 * Gets the y component of the location of this spin
	 * Uses the spins ordinal index number
	 * @param Lz number of unit cells in the z direction of the system
	 * @param Ly number of unit cells in the y direction of the system
	 * @return the y component of the location of this spin
	 */
	public double getY(int Lz, int Ly){
		return getLocation(Lz, Ly).getY();
	}

	/**
	 * Gets the z component of the location of this spin
	 * Uses the spins ordinal index number
	 * @param Lz number of unit cells in the z direction of the system
	 * @param Ly number of unit cells in the y direction of the system
	 * @return the z component of the location of this spin
	 */
	public double getZ(int Lz, int Ly){
		return getLocation(Lz,Ly).getZ();
	}
	
	
	
}
