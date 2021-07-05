package simulation.montecarlo;

class deltaEnergyAndLattice {
	private Lattice lattice;
	private double deltaEnergy;
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
