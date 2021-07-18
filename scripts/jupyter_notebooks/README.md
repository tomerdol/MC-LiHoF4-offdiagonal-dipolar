# Jupyter notebooks

This directory contains Jupyter notebooks that were/are used for data analysis,
plotting and testing new ideas.  
There are two ways to work with them:
1. Locally: start the jupyter notebook server in this directory on the PC
   and work on them from there -- easier but not very good for heavy computations.  
   In Windows, using WinPython, run the following in cmd:
   ```
   "C:\WPy64-3770\Jupyter Notebook.exe" --notebook-dir C:\Users\Tomer\IdeaProjects\LiHoF4_transverse_field_MC\scripts\jupyter_notebooks\
   ```
2. Remotely: start the jupyter notebook server on the SGE cluster and connect 
   to it from the PC. This is done by running the script `run_jupyter.sh` and following 
   the instructions (best to submit it w/ qsub and inspecting the output file for the instructions).
   This option is preferable for heavier computations.