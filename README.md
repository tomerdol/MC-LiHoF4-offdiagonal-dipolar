# Monte Carlo study of LiHoF<sub>4</sub> under a weak transverse magnetic field 

**Version 1.5.1** - [Change log](changelog.md)

## Overview
Software used to study LiHoF<sub>4</sub> and other dipolar magnetic systems.\
Used for the paper *Effect of intrinsic quantum fluctuations on the phase diagram of anisotropic dipolar magnets* by T. Dollberg, J. C. Andresen, M. Schechter, [Phys. Rev. B 105, L180413 (2022)](https://link.aps.org/doi/10.1103/PhysRevB.105.L180413).\
This repository contains the main Java simulation software (in [src](./src)), Java dependencies (in [lib](./lib)) and various Python and shell scripts (in [scripts](./scripts)) used to monitor the simulations and analyze their results.\
Different studied systems use separate working directories (LiHoF4 and Fe8) where simulation data, figures and checkpoints are stored.\
This software is run under the Sun Grid Engine ([SGE](http://gridscheduler.sourceforge.net/htmlman/manuals.html)).

## How to compile?
The project was developed using IntelliJ IDEA (2021.*, Ultimate edition) and is most easily compiled inside the IDE.\
If unavailable, it is also possible to compile using the primary Java compiler, javac. All dependencies are kept as .jar files in [lib](./lib).\
* On Linux:
```
find src/ -name *.java > sources.txt  
javac -sourcepath src -d bin/production/LiHoF4_transverse_field_MC/ -classpath "./lib/*" @sources.txt
rm -f sources.txt
```  
* On Windows:
```
dir /s /B src\*.java > sources.txt
javac -sourcepath src -d bin/production/LiHoF4_transverse_field_MC/ -classpath "./lib/*" @sources.txt
del sources.txt
```

## How to run?
Typically, many independent runs of the simulation must be run to get meaningful results. This is done using the [sub.sh](./scripts/sub.sh) and [init_sim.sh](./scripts/init_sim.sh) srcipts.
The following is a minimal working example for running the Java simulation itself after compilation:
```
java -Dsystem=LiHoF4 -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L 4,4 -d 1.0 -extBx 0.0 -continue_from_save no -max_sweeps 4096 -name default -temp_schedule ./LiHoF4/temperature_schedules/temp_schedule_1_1_exact_test_const_0.0_false.txt -alpha 0.95 -mode s -seed 123456789 -Jex 1.16e-3 -interpolation_table_name _const -p -print_sweep_num 1 -tol 5.0e-4
```
For more information on command line arguments run
```
java -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main --help
```

### Important directories

### `LiHoF4` and `Fe8`
This software is capable of simulating more than a single system.
Currently, two types of systems are defined: LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub>
and the Fe<sub>8</sub> molecular magnet. Each such project's data and functional files 
are kept in a directory with the corresponding name. 

### `results`
Observables are measured following each MC sweep, and outputted to a corresponding 
text file located in `./$SYS_NAME/data/results/$project_name/table_$Lx_$Lz_$extBx_$T_$mech_$seed.txt`.
Within this file lines starting with '#' include simulation metadata, i.e. software version, 
parameter values, etc.
The first row not starting with '#' lists the names of the columns, separated by whitespace, 
and all subsequent rows contains the measured values.  
If the simulation is restarted following a termination, the metadata is printed again 
and the column headers are printed again but this time with a preceding '#'.  
This way, by removing all lines starting with '#' one can obtain a table of observables
measured along the simulation.

### `binned_data`
The analysis process includes storing the raw results in bins of logarithmically increasing
lengths and tracking their evolution over MC time to determine when the simulation has become
equilibrated. Then, only data from the last bin is used in any further analysis.  
Each directory in `results` also contains a subdirectory called `binned_data` which 
contains HDF5 files that hold the binned results. These are the files used for the analysis
by the scripts in the `scripts` directory.  
Once simulations are fully completed (no further sweeps are expected to be run), the results
and the checkpoint files are archived, and only the binned data files are easily accessible.

### `checkpoints`
The state of the simulation is periodically saved using the Java serialization process
so that it may be continued later in case of an unexpected crash or in case more MC sweeps
are required for equilibration.
These files are saved in subdirectories of `checkpoints` with their respective project names.

### `analysis`
While initiating a simulation, the standard process includes running a short simulation
with a wide range of temperatures, and automatically analysing its results to determine
the optimal temperature sequence to use in a longer simulation. This process, performed
by the `init_sim.sh` script outputs data to the `analysis` directory.

### `interactions`
Keeps data calculated once before the simulation and used during it. Includes: 
1. Interactions between pairs of spins in the presence of periodic boundary conditions using the Ewald summation
   method. file names: `intTable_$Lx_$Lz.txt`.
2. Interpolation table data for the single-spin magnetic moment and energy.
   file names: `magnetic_moment_up_arr_$extBx_const.txt`, `energy_up_arr_$extBx_const.txt`.

### `p_configs`
When the simulation encounters a configuration for which it cannot solve the 
self-consistent calculation, it saves the configuration preceding the last spin-flip so that
the problematic configuration can be replicated and studied.
This was used at the development stage but all files in this directory should be empty
in the production stage.

### `temperature_schedules`
Each simulation runs several temperatures in parallel (parallel tempering). The sequence of
temperatures used for each simulation is kept in text files in this directory,
separated by commas.
