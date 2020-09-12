#!/bin/bash
#$ -q smoshe.q@sge1081,smoshe.q@sge1082,lublin.q
#$ -S /bin/bash
#$ -V
#$ -cwd

# get tunneling info
node=$(hostname -s)
user=$(whoami)
cluster="hpcin"
port=8889

# print tunneling instructions jupyter-notebook
echo -e "
Command to create ssh tunnel:
ssh -N -f -L ${port}:${node}:${port} ${user}@${cluster}.bgu.ac.il

Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"

# Run Jupyter
jupyter-notebook --no-browser --port=8889 --ip=0.0.0.0
