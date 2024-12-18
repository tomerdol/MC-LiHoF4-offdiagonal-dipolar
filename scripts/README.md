##Scripts for controlling, monitoring and analyzing Monte Carlo simulations

### System name
Some of the scripts use an environmental variable named `SYS_NAME` which holds the name of the
system currently under consideration. Can be either 'LiHoF4' or 'Fe8'.  
This variable can be set to a default value by adding it to `~/.cshrc`:
```
vim ~/.cshrc
setenv SYS_NAME 'LiHoF4'
:wq
source ~/.cshrc
```

### Create interactions tables using the Ewald method
From the root path:
```
ewSum.sh $Lx $Lz $k_cutoff $real_cutoff
```
or (replace L with the linear system size)
```
qsub -cwd -q smoshe.q,lublin.q -V -S /bin/bash -N "ewSumL" ./scripts/ewSum.sh L L 12 10
```

### Create magnetic moment and energy tables
From the root path:
1. Edit `create_all_tables.sh` to have the desired external transverse fields
2. Run `create_all_tables.sh` (which runs `create_one_table.sh`, which in turn runs `crystal+field+hamiltonian-transversal+field+const.py`:
    ```
    qsub -cwd -q smoshe.q,lublin.q -V -S /bin/bash scripts/create_all_tables.sh
    ```
   
### python3 path
Some scripts use a specific python path:
```
/gpfs0/smoshe/projects/Python-3.8_old/bin/python3
```
This is done since python compiled on newer machines does not work on the older 
machines (sge190), so this way we ensure it will run wherever.  
If all machines have the same OS version, these references can just be changed to `python3`
that is in the path.