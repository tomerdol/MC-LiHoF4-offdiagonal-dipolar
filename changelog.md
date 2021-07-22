# Change log

## v1.5.1 - 22/07/2021
* BUGFIX in output of final simulation state:
  1. interactions table was read without dilution
  2. interactions table was only read at the end in the case of suppression of transverse fields, so the exchange interaction remained in the other case
* Added option to force binning
* Added exception when binned data is read that seeds that do not include a given temperature (perhaps they were created during the binning process) are just skipped.
* Added Jupyter notebooks
* Moved definition of Jex to the submitting scripts
* Added important directory info to the README
* Added MSc thesis to docs
* Extended documentation and removed unused scripts

## v1.5.0 - 13/07/2021
* Removed tMaxConfig column
* Added junit4.jar to lib
* Extended README.md
* Fixes to method statistics for self-consistent calculation
* Added summary of useful commands in /docs
* removed tMAxConfigs column from BIN output
* Extended documentation and Javadoc
* Added Ion_positions.pdf
* Switched to try-with-resources in receiveTemperatureSchedule
* Removed printing of percentage of longitudinal fields below some threshold (0.029) - no longer relevant with constant moment
* Added usage message in case no arguments are given
* Changed nnArray_test in ReadInteractionsTableLiHoF4.java to compact indices
* Moved handling of IndexOutOfBoundsException when calculating a derivative to the calling function so that it causes the method to fail and tries the next method instead of immediately aborting the MC step
* Fixed table formatting for BIN output type
* Ewald obtains its singleSpin array from lattice
* deprecated GenerateLattice

## v1.4.0b - 05/07/2021
* Add script to plot multiple histograms
* Changed crystal+field+hamiltonian-transversal+field+const.py plateau location back to 1.1
* Plot (using Mathematica) the lattice output of Fe8
* Added option to simulate diluted system
* Seed must be given to simulation (as commandline argument)
* BUGFIX call to find_inital_tc.py from init_sim.sh to include system name
* Print concentration (x) to list of parameters
* Add system name to get_avgs.py
* Average only uncommented lines in init_sim.sh
* Added serial version UID to singleSpin
* Added python exception for missing binned files in bin_data.py
* Simulations with no binned files are just dropped from the analysis
* Allow usage of last bin in "3 consecutive bins" as equilibrated bin
* Removed real time equilibration testing (unused feature)
* Removed counting of transverse-field-maximizing nearest neighbor configurations
* Added script archive_binned_data.sh
* Binning data using hdf5
* bin_data.py to accept commandline arguments like plot_bin.py
* analysis_tools.py to ignore non-existent simulations
* Print which simulations specifically are not equilibrated during plot_equilibration_pdf_bin.py
* resubmit_to_all.sh to check that the new job is running before grabbing the next queued job
* FieldTable to specify given field values when throwing IndexOutOfBoundsException
* Specify file name when unable to read FieldTable
* Number of temperatures (and therefore needed cores) is saved in a variable $nT in submission scripts
* Add statistics on methods used to solve self-consistent calculation
* Added verification of dilution compatibility when loading checkpoint

## v1.3.0b - 24/05/2021
* Restructure for multiple projects
* New module to parse arguments for all scripts and hold variables that would be accessible from other scripts
* Added config module with the name of the system and added it to all paths
* Removed old functions from phase_diagram_bin.py
* Removed MagneticMomentsSolveIter.java and MonteCarloMetropolis.java
* Add extBy
* New class for system specific nearest neighbors configurations
* Added option to plot histogram of specific run
* Fixed reading L=10 obsPrintSweepNum
* Added primitive lattice vectors printing to results file
* Added function to calculate simple dipolar interactions with PBC in ewaldSum
* Upgraded find_unfinished_simulations.sh to find max steps
* Added single submission to sub scripts
* Upgraded find_unfinished_simulations.sh
* Add delta x as command line argument to phase_diagram_bin.py
* Changed plot_bin.py to use latest bin
* Cosmetic fix to plot_equilibration_pdf_bin.py
* BUGFIX: plot_histogram.py can now process simulations with different numbers of independent runs
* Show warning when using unequilibrated data
* Add swap to binned data
* Enable acceptance rate plot with plot_bin.py
* Add script to create magnetic moment table according to the calculation in P. B. Chakraborty, P. Henelius, H. Kjønsberg, A. W. Sandvik, and S. M. Girvin, Phys. Rev. B 70, 144411 (2004).
* (31/5/2021) Fixed double multiplying by 0.5 for self-interactions

## v1.2.0b - 08/03/2021
* Write full lattice state at the end of the simulation
* Count Bx-maximizing configurations
* Add fraction of transverse field maximizing configurations to thermal averages
* Add script print_lattice_from_existing_checkpoints.sh to rerun finished simulations just to write the lattice state
* Add histogram plotting
* Tighter verification that no simulations are lagging behind.
* Add scripts to look for unfinished simulations: find_unfinished_simulations.sh
* Plot equilibrated bin vs. temperature when running plot_equilibration_pdf_bin.py
* Option to use correlation lengths along different axes for FSS
* Changed to more compact way of finding nearest neighbor numbers
* Add correlator plotting to plot_bin.py
* Add crystal+field+hamiltonian-transversal+field+const.py to create magnetic moment tables with constant Bz
* create_crystal_field_table.py: Create tables with hyperfine interactions (and rotation)
* phase_diagram_bin.py: delta around initial T_c is chosen symmetric even if initial T_c is close to the edge of the temperature range
* Add option to shift the T axis in plot_bin.py so that multiple Bx's could be compared
* Fix cluster name in run_jupyter.sh
* Add exchange interaction directly to dipolar interaction table (so that it participates in the self consistent calculation)
* CRUCIAL BUGFIX: call to manualCalcValue always called for energy and never for magnetic moment!
* Improved finding of initial fitting parameters in fit6.py
* get_avgs.py to get data from all files in /output and calculate runtimes by simulation type 
* Added new homotopic method for self-consistent calculation
* Added possibility to exclude system sizes in sub.sh and in init_sim.sh
* Get user input for system sizes in gen_temp_schedule.py

## v1.1.0b - 24/09/2020
* Added measurement of mk2 for k in the y and z directions.
* Added a module that reads old checkpoints and re-saves in new format (with extra column for additional observable).
* Added serialVersionUID to all serialized objects for better backwards compatibility.
* Added the fraction of transverse-field-maximizing nearest neighbor configurations to measured observables.
* Write full lattice state at the end of the simulation for enhanced spin-specific measurements to be performed after the simulation has finished.
    * Write transverse fields also for simulations where they are suppressed.
* Fixed possible decrease of maxSweeps when reading checkpoint while giving max_sweeps that is lower than the original value.
* Better observable printout (column alignment).
* Fixed path name for binned data.
* Remove tmp file after analysis has finished.
* Fixed temperature schedule creation based on initial Tc in init_sim.sh
* Added interpolation table name to submission scripts

## v1.0.1b - 15/09/2020
* Fixed scripts for creating interpolation table.
* Name extension for field interpolation table file can be given as commandline argument.
* Removed 'broyden' from file names.
* Added crystal+field+hamiltonian-transversal+field+broyden.py and initiating scripts that create the interpolation tables.

## Before 1.0.1b:
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