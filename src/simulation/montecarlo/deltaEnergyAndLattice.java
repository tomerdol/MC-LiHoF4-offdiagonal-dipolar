package simulation.montecarlo;

/**
 * Object returned by the Monte Carlo step that holds the new lattice (spin configuration),
 * the energy change, and the method used for the self-consistent solution
 */
class deltaEnergyAndLattice {
	/** The new lattice object */
	private Lattice lattice;
	/** The energy difference due to the MC step (single spin flip) */
	private double deltaEnergy;
	/** The index of the method used for solving the self-consistent calculation
	 * following the spin flip */
	private int methodUsed;

	public deltaEnergyAndLattice(Lattice lattice, double deltaEnergy, int methodUsed) {
		this.lattice = lattice;
		this.deltaEnergy = deltaEnergy;
		this.methodUsed = methodUsed;
	}

	public Lattice getLattice() {
		return lattice;
	}

	public void setLattice(Lattice lattice) {
		this.lattice = lattice;
	}

	public double getDeltaEnergy() {
		return deltaEnergy;
	}

	public void setDeltaEnergy(double deltaEnergy) {
		this.deltaEnergy = deltaEnergy;
	}

	public int getMethodUsed() {
		return methodUsed;
	}
}
