package simulation.montecarlo;


import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.Serializable;

public class singleSpin implements Serializable{
	private static final long serialVersionUID = 625318162648175439L;
	private int s;			// spin
	private double spinSize;	//spin size
	private double localBz;	// local field - comparable to EO fitness
	private int n;			// ordinal number
	private double localBx;	// local field in the x direction (only allowed [-4,4])
	private double localBy;	// local field in the y direction (only allowed [-4,4])
	private int[][] prevBIndices;	// first array contains the last indices of the local field (x,y,z)<->(2,1,0) in the FieldTable
									// second array is used to track whether consecutive calls (to find magneticMoment or energy) are correlated (separately for x,y,z)

	public singleSpin(){
		this(1, 0, 0, 0, 0, 0);
	}

	// usual constructor
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

	public singleSpin(int s0, int n0, double localBz, double localBy, double localBx, double spinSize0, int[][] prevBIndices0){
		this(s0,n0,localBz,localBy,localBx,spinSize0);

		if (prevBIndices0==null || prevBIndices0[0]==null || prevBIndices0[1]==null) {
			throw new NullPointerException("prevBIndices cannot be null. if it is not known it should not be passed and a new one will be initialized.");
		}
		if (prevBIndices0.length!=2 || prevBIndices0[0].length!=3 || prevBIndices0[1].length!=3){
			throw new IllegalArgumentException("received prevBIndices has wrong dimensions. should be 2x3");
		}

		//******************************* copy the input array *********************************************************
		System.arraycopy(prevBIndices0[0], 0, prevBIndices[0], 0, prevBIndices0[0].length);
		System.arraycopy(prevBIndices0[1], 0, prevBIndices[1], 0, prevBIndices0[1].length);
		//**************************************************************************************************************
	}


	public singleSpin(int s0, int n0){
		this(s0,n0,0,0,0,0);	// Bx,By,Bz and spinSize should be initialized later
	}

	public singleSpin(int s0, int n0, double spinSize){
		this(s0,n0,0,0,0,s0*(spinSize));	// Bx,By,Bz should be initialized later and spinSize is an initial guess
	}

	// copy constructor
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
	public double getLocalBx(){	return this.localBx; }
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

	//toString prints the spin & position of the spin in the following form: 'Spin:-1 Position:(x,y,z)'
	public String toString(){
		return "Spin:"+s+" Number: "+this.n + " Local Bz:"+this.localBz+" Local By: "+this.localBy
				+ " Local Bx: " + this.localBx + " Spin Size: " + this.spinSize;
	}
	
	/**
	 * distance function for periodic boundary conditions
	 * @param another - spin to which the distance will be calculated
	 * @param Lx - number of unit cells in the x-direction
	 * @param Lz - number of unit cells in the z-direction
	 * @return array of [ds,dx,dx,dz]
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
	public void flipSpin(){
		s=s*(-1);	// flip the spin
		prevBIndices[1][0]=0;	// change the correlation of Bz to false
	}

	@Deprecated
	public double calcDipolarInteraction(singleSpin o, int Lz, int Lx){
		if (this==o)	return 0;
		double r=this.distance(o, Lz, Lx)[0];
		double rz=this.getZ(Lz,Lx)-o.getZ(Lz,Lx);
		//double rx=i.getX()-j.getX();
		double actual_height=Lz*Constants.c;
		double actual_length=Lx*Constants.a;
		// implementation of periodic boundary conditions
		if (rz>actual_height*0.5)
			rz=rz-actual_height;
		if (rz<=-actual_height*0.5)
			rz=rz+actual_height;
		
		return (this.s*o.s)*(r*r-3*rz*rz)/Math.pow(r, 5);
		
	}

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

	// location getters (calculation)
	public double getX(int Lz, int Ly){
		return getLocation(Lz, Ly).getX();
	}

	public double getY(int Lz, int Ly){
		return getLocation(Lz, Ly).getY();
	}

	public double getZ(int Lz, int Ly){
		return getLocation(Lz,Ly).getZ();
	}
	
	
	
}
