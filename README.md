# Monte Carlo study of LiHoF<sub>4</sub> under a weak transverse magnetic field 

**Version 1.3.0b** - [Change log](changelog.md)

## Overview
Software used to study LiHoF<sub>4</sub> and other dipolar magnetic systems.\
This repository contains the main Java simulation software (in [src](./src)), Java dependencies (in [lib](./lib)) and various Python and shell scripts (in [scripts](./scripts)) used to monitor the simulations and analyze their results.\
Different studied systems use separate working directories (LiHoF4 and Fe8) where simulation data, figures and checkpoints are stored.\
This software is run under the Sun Grid Engine ([SGE](http://gridscheduler.sourceforge.net/htmlman/manuals.html)).

### How to compile?
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

### How to run?
Typically, many independent runs of the simulation must be run to get meaningful results. This is done using the [sub.sh](./scripts/sub.sh) and [init_sim.sh](./scripts/init_sim.sh) srcipts.
The following is a minimal working example for running the Java simulation itself after compilation:
```
java -Dsystem=LiHoF4 -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main -max_iter 80 -L 4,4 -d 1.0 -extBx 0.0 -continue_from_save no -max_sweeps 4096 -name default -temp_schedule ./LiHoF4/temperature_schedules/temp_schedule_1_1_exact_test_const_0.0_false.txt -alpha 0.95 -mode s -seed 123456789 -Jex 1.16e-3 -interpolation_table_name _const -p -print_sweep_num 1 -tol 5.0e-4
```
For more information on command line arguments run
```
java -classpath "./bin/production/LiHoF4_transverse_field_MC/:./lib/*" simulation.montecarlo.Main --help
```