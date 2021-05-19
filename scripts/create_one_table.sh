#python3 scripts/crystal+field+hamiltonian-transversal+field+broyden.py "$1"
#python3 scripts/crystal+field+hamiltonian-transversal+field+const.py "$1"
if [ $SYS_NAME == "LiHoF4" ]; then
#    python3 scripts/crystal+field+hamiltonian-transversal+field+broyden.py "$1"
    python3 scripts/crystal+field+hamiltonian-transversal+field+const.py "$1"
elif [ $SYS_NAME == "Fe8" ]; then
  python3 scripts/crystal+field+hamiltonian-transversal+field+broyden_Fe8.py "$1" "$2"
fi


