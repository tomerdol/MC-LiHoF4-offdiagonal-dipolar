# Change log

# v1.1.0b - 24/09/2020
* Added measurement of mk2 for k in the y and z directions.
* Added a module that reads old checkpoints and re-saves in new format (with extra column for additional observable).
* Added serialVersionUID to all serialized objects for better backwards compatibility.
* Added the fraction of transverse-field-maximizing nearest neighbor configurations to measured observables.
* Write full lattice state at the end of the simulation for enhanced spin-specific measurements to be performed after the simulation has finished.
    * Write transverse fields also for simulations where they are suppressed.
* Fixed possible decrease of maxSweeps when reading checkpoint while giving max_sweeps that is lower than the original value.
* Better observable printout (column alignment).
* Fixed path name for binned data.
5ec2b3c Fix subscript in README.md
ee13d67 * Add README.md * Add changelog.md
45e6a32 Remove tmp file after analysis is finished
2314d5a Fixed temperature schedule creation based on initial Tc in init_sim.sh
e410181 Added interpolation table name to submission scripts
9d3a4fb VERSION 1.0.1b
c93448d Fixed scripts for creating interpolation table
fe5dba0 * Name extension for field interpolation table file can be given as commandline argument. * Removed 'broyden' from file names
36aee54 (HEAD -> beta) Merge branch 'dev' into beta
a2da57c * sync with /beta/ * change delta in fitting to 0.04
a9b0f45 (dev) sync files with /dev/
ed3c52d Add measurement of mk2 for k in the y and z directions
08b8ce3 (origin/dev) Change submission file
006bd48 Add module that reads old checkpoints and re-saves in new format with extra column
1a6fe25 Add serialVersionUID to all serialized objects for backwards compatibility
f54c7a5 Add fraction of transverse field maximizing configurations to thermal averages
8715811 Count Bx-maximizing configurations
10e7042 * Add serialVersionUID for backwards compatibility * Experimental histogram plotting
38a110b Open lattice_output file only just before writing, once simulation is finished
8a491c7 Write full lattice state at the end of the simulation
81e24e5 Fix possible decrease of maxSweeps
6e9530f FIX - Add space in front of positive integer (sweep and spin) in printout for better column alignment FIX - When using SPIN verbocity, print transverse fields for simulations where they are suppressed as well
37c4fd0 Overwrite output file for SPIN output
8fc4169 Add space in front of positive integer (sweep and spin) in printout for better column alignment
c0c29df Remove @CreatesInconsistency warning for backwards compatibility
1e1f837 When using SPIN verbocity, print transverse fields for simulations where they are suppressed as well
c60e935 Add option to print entire lattice INSTEAD of running the simulation
1986967 cleanup
19789d1 Merge branch 'beta' into dev
6fc63a6 (origin/beta) Fix path name for binned data
c5b66d4 Fix path name for binned data
08a1a00 * Update submission scripts * Testing Ewald sum
5ec2b3c Fix subscript in README.md
ee13d67 * Add README.md * Add changelog.md
45e6a32 Remove tmp file after analysis is finished
2314d5a Fixed temperature schedule creation based on initial Tc in init_sim.sh
e410181 Added interpolation table name to submission scripts

# v1.0.1b - 15/09/2020
* Fixed scripts for creating interpolation table.
* Name extension for field interpolation table file can be given as commandline argument.
* Removed 'broyden' from file names.
* Added crystal+field+hamiltonian-transversal+field+broyden.py and initiating scripts that create the interpolation tables.

# Before 1.0.1b:
30689fa Merge branch 'beta' into dev
9be1aa6 Added version printing
39feb77 Added version printing
89b97a4 Attempt to exclude deployment.xml from mergers
2930acd Excluded sge output files from sync
9ee4f49 Added argument check to cont_sim.sh
67a3c26 * cont_sim.sh script updated to check if seed is already running before running new job * reverted parameters.properties to original obsPrintSweepNum values
4efe793 fixed line breaks in shell scripts
674ba9c bug fix
53fd8c4 Create separate IntelliJ module for Python scripts
2dfa986 update fields from abstract MonteCarloSimulation (maxSweeps, continueFromSave...) while reading from checkpoint
4db537a (origin/bug) moved check for parallelTemperingFalse before energies are calculated
7a2f646 moved check for parallelTemperingFalse before energies are calculated
e922a16 Merge branch 'master' into dev
3199d45 (origin/master, origin/beta, master) fixed deployment.xml
b6eb8e8 change default deployment server to /dev/
815782b change default deployment server to /stable/
044986c change default deployment server to /beta/
fd00354 change default deployment server to /beta/
ad2867f refactoring folder structure
60229b9 * Added example result files (checkpoints, data, temperature_schedules, interactions) * Added scripts folder (Python and shell scripts)
aa0c1bc refactored to IdeaProjects
572d4d7 refactored to IdeaProjects
2c3f78f refactored to IdeaProjects
6e3c1a8 seed included in file name for multiple parallel runs
bdc3bd9 seed included in file name for multiple parallel runs
8c8e2f4 add seed compatibility check when restarting from checkpoint
5d8c6ef back to regular magnetic moment calculation (interpolation)
230b2a4 calculating moments and energy manually for |Bz|<0.029 to check how much of an effect the interpolation has in that range
299bb48 table headers now written together with parameters and have preceding '#' when state is successfully read from file
6563bbf Fixed all TO-DO items
b1240da Switched to new equilibration criterion
f658a34 Added annotation to mark methods that create an inconsistency between moments and fields
bbf32fe moved simulation specific constants to the simulation objects
3c49dfb moved simulation specific constants to the simulation objects
a88f1bd before seeding change
f3e57dd initial parallelization working
9a25ecd minor fixes and now it is running
24fba8a finished transition