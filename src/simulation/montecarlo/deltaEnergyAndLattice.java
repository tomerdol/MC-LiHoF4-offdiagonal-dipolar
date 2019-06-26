package simulation.montecarlo;

class deltaEnergyAndLattice {
	private Lattice lattice;
	private double deltaEnergy;

	public deltaEnergyAndLattice(Lattice lattice, double deltaEnergy) {
		this.lattice = lattice;
		this.deltaEnergy = deltaEnergy;
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
}
